"""exp_envelope_peaks.py -- per-run envelope peaks of US-1 (Test 9) and US-2 (Test 12) from
EXP_DATA/processed_globalzero, with the 2.317348 mass factor removed.

INPUT   <proc_dir>/Test{9,12}Run<N>_processed_globalzero.xlsx
        columns: Time [s], U_avg [m, global zero], base_shear_N, base_shear_kN, mass_scale, m15_kg, m16_kg, m17_kg
FACTOR  F = -base_shear_kN / mass_scale, mass_scale read from each file (2.317348). The masses in the file
        (946.68 / 743.09 / 936.50 kg) already contain the factor; dividing restores 408.52 / 320.66 / 404.13 kg.
        Sign -1: restoring force positive with positive displacement (as in baseshear_from_exp.py).
        Checked per run: the correlation of d and F over the strong-motion window must be positive.
PER RUN (all channels low-passed 25 Hz, 4th-order zero-phase, identical on d and F so pairs stay simultaneous)
  d      = U_avg x 1000, GLOBAL zero as stored (accumulated residual included), mm
  F      = force, zeroed by its pre-event median (accelerometer offset), kN
  sync   = d vs F cross-correlation over +-4 s; a slip > 3 samples with poor zero-lag correlation is corrected
           by shifting d, only if the correlation at the best lag exceeds 0.5 (US-1 run 18: channels 1-7 lead
           channels 8-17 by ~1.36 s in the raw file; US-2 run 3 gives a false far lag on potentiometer noise)
  simultaneous pairs, per sign:  d_max and F at that instant  |  F_max and d at that instant
  non-simultaneous peaks:        max |d| and max |F| of the run (the Fig. 13 / Table 4 convention), also raw (unfiltered)
OUTPUT  <out>/envelope_peaks_US1.csv, envelope_peaks_US2.csv, fig_envelope_peaks.png
USAGE   python exp_envelope_peaks.py <proc_dir> <out_dir>
"""
import sys, os, re, glob, csv
import numpy as np
from scipy.signal import butter, filtfilt, correlate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LP_HZ = 25.0
SIGN = -1.0
TESTS = (("US1", "Test9"), ("US2", "Test12"))


def load(path):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True); hdr = list(next(it)); c = {h: j for j, h in enumerate(hdr)}
    rows = [r for r in it if r[c["Time"]] is not None]
    t = np.array([r[c["Time"]] for r in rows], float)
    U = np.array([r[c["U_avg"]] for r in rows], float)
    V = np.array([r[c["base_shear_kN"]] for r in rows], float)
    ms = float(rows[0][c["mass_scale"]])
    return t, U * 1e3, SIGN * V / ms, ms


def lp(x, fs):
    b, a = butter(4, LP_HZ / (0.5 * fs)); return filtfilt(b, a, x)


def bp(x, fs):
    b, a = butter(4, [0.5 / (0.5 * fs), 25.0 / (0.5 * fs)], "band"); return filtfilt(b, a, x - np.mean(x))


def shift(x, k):
    if k == 0: return x
    return np.concatenate([np.full(k, x[0]), x[:-k]]) if k > 0 else np.concatenate([x[-k:], np.full(-k, x[-1])])


def analyse(path):
    t, d_raw, F_raw, ms = load(path)
    fs = 1.0 / np.median(np.diff(t))
    uu, ff = bp(d_raw, fs), bp(F_raw, fs)
    cc = correlate(ff, uu, "full", method="fft"); L = np.arange(-len(uu) + 1, len(uu)); m = np.abs(L) <= int(4 * fs)
    kb = int(L[m][np.argmax(cc[m])]); c0 = cc[L == 0][0] / cc[m].max()
    rho_best = cc[m].max() / (np.linalg.norm(uu) * np.linalg.norm(ff))
    # a slip only if the shifted pair correlates well; on tiny runs (pot noise) a far lag is a false hit
    k = kb if (abs(kb) > 3 and c0 < 0.5 and rho_best > 0.5) else 0
    d_raw = shift(d_raw, k)
    E = np.cumsum(ff ** 2); E /= E[-1]
    i0, i1 = int(np.searchsorted(E, 0.01)), int(np.searchsorted(E, 0.99))
    pre = slice(max(i0 - int(2 * fs), 0), max(i0 - int(0.2 * fs), 1))
    F_raw = F_raw - np.median(F_raw[pre])
    d, F = lp(d_raw, fs), lp(F_raw, fs)
    w = slice(i0, i1)
    rho = float(np.corrcoef(bp(d, fs)[w], bp(F, fs)[w])[0, 1])
    r = dict(lag=k, rho=rho, mass_scale=ms, d_start=float(np.median(d[pre])), d_end=float(np.median(d[-int(2 * fs):])))
    for s, key in ((1, "pos"), (-1, "neg")):
        i = i0 + int(np.argmax(s * d[w])); j = i0 + int(np.argmax(s * F[w]))
        r[key] = (d[i], F[i], F[j], d[j])
    r["dabs_raw"] = float(np.max(np.abs(d_raw[w]))); r["Fabs_raw"] = float(np.max(np.abs(F_raw[w])))
    r["dabs"] = float(np.max(np.abs(d[w]))); r["Fabs"] = float(np.max(np.abs(F[w])))
    return r


def main():
    proc, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
    allres = {}
    for lab, test in TESTS:
        files = sorted(glob.glob(os.path.join(proc, "{}Run*_processed_globalzero.xlsx".format(test))),
                       key=lambda p: int(re.search(r"Run(\d+)_", p).group(1)))
        R = {}
        print("== {} ({}): {} runs".format(lab, test, len(files)))
        for p in files:
            n = int(re.search(r"Run(\d+)_", p).group(1)); R[n] = r = analyse(p)
            print("  run {:2d}  lag {:+4d}  rho {:+.2f}  +d {:6.2f} @F {:6.2f} | -d {:7.2f} @F {:6.2f} | Fmax +{:5.2f} -{:5.2f} | |d|max {:6.2f} |F|max {:5.2f} (raw {:5.2f})".format(
                n, r["lag"], r["rho"], r["pos"][0], r["pos"][1], r["neg"][0], r["neg"][1], r["pos"][2], -r["neg"][2],
                r["dabs"], r["Fabs"], r["Fabs_raw"]))
        allres[lab] = R
        with open(os.path.join(out, "envelope_peaks_{}.csv".format(lab)), "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["run", "sync_shift_samples", "rho_dF", "mass_scale_removed", "d_start_mm", "d_end_mm",
                         "d_max_pos_mm", "F_at_dmax_pos_kN", "d_max_neg_mm", "F_at_dmax_neg_kN",
                         "F_max_pos_kN", "d_at_Fmax_pos_mm", "F_max_neg_kN", "d_at_Fmax_neg_mm",
                         "abs_d_max_mm", "abs_F_max_kN", "abs_d_max_raw_mm", "abs_F_max_raw_kN"])
            for n, r in R.items():
                wr.writerow([n, r["lag"], "{:.3f}".format(r["rho"]), "{:.6f}".format(r["mass_scale"])] + ["{:.4f}".format(x) for x in (
                    r["d_start"], r["d_end"], r["pos"][0], r["pos"][1], r["neg"][0], r["neg"][1],
                    r["pos"][2], r["pos"][3], r["neg"][2], r["neg"][3], r["dabs"], r["Fabs"], r["dabs_raw"], r["Fabs_raw"])])
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
    for a, (lab, test) in zip(ax, TESTS):
        R = allres[lab]; runs = sorted(R)
        for s, k, c in ((1, "pos", "#1f6f8b"), (-1, "neg", "#1f6f8b")):
            dd = [R[n][k][0] for n in runs]; ff = [R[n][k][1] for n in runs]
            o = np.argsort(np.abs(dd))
            a.plot(np.r_[0, np.array(dd)[o]], np.r_[0, np.array(ff)[o]], "o-", color=c, ms=4, lw=1,
                   label="simultaneous: F at d_max" if s > 0 else None)
            dF = [R[n][k][3] for n in runs]; FF = [R[n][k][2] for n in runs]
            a.plot(dF, FF, "^", color="#4a8f5a", ms=4, mec="k", mew=0.3, label="simultaneous: F_max and d there" if s > 0 else None)
            nd = [R[n][k][0] for n in runs]; nF = [R[n][k][2] for n in runs]     # own-sign peaks, not simultaneous
            o = np.argsort(np.abs(nd))
            a.plot(np.r_[0, np.array(nd)[o]], np.r_[0, np.array(nF)[o]], "s--", color="#b4491a", ms=3.5, lw=0.9,
                   label="non-simultaneous: d_max, F_max per sign (Fig. 13 style)" if s > 0 else None)
        for n in runs:
            if n in (runs[-1], runs[-2]) or n % 4 == 0:
                a.annotate(str(n), (R[n]["pos"][0], R[n]["pos"][1]), fontsize=6, xytext=(3, 2), textcoords="offset points")
        a.axhline(0, color="k", lw=0.5); a.axvline(0, color="k", lw=0.5); a.grid(alpha=0.25)
        a.set_title("{} ({}), runs 1-{}, force / {:.4f}".format(lab, test, runs[-1], R[runs[0]]["mass_scale"]), loc="left", fontsize=10)
        a.set_xlabel("d at 2.06 m, global zero (mm)")
    ax[0].set_ylabel("F = -sum m_j a_j / mass_scale (kN)"); ax[0].legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_envelope_peaks.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
