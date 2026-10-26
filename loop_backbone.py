"""loop_backbone.py -- experimental force-displacement backbone of US-1 from SIMULTANEOUS (d, F) pairs.

Why: Fig. 13 of Moshfeghi et al. (2024) is a locus of (peak |d|, peak |F|) taken from different instants of
each run. On a 10 Hz wall the two peaks are not simultaneous, so its initial slope under-reads the stiffness
(3.2 kN/mm -> m_eff = 672 kg with T = 0.091 s). Here every point is a pair taken at ONE instant.

INPUT   <v2_dir>/Test9Run<N>_processed_globalzero.xlsx   (make_processed_globalzero_v2.py)
            Time [s], U_avg [m] (0.5(ch3+ch4), relative to table, global zero), base_shear_N (sum m_j a_j, absolute)
        --table5 us1_fig14_digitised.csv  (run, T_exp_s)            optional, for the mass check
        --fig13  us1_fig13_digitised.csv  (drift_pct, F_kN, d_mm)   optional, overlay, drawn / 2.317348
        --pushover <csv with d_mm, F_kN columns>                    optional, model backbone overlay
SIGNS   F = -base_shear (restoring force, same sign as d), as in baseshear_from_exp.py.

PER RUN
  1. sync check: U and F band-passed 0.5-25 Hz, cross-correlated over +-4 s. A lag beyond 3 samples with
     the zero-lag correlation below half the best is a DAQ slip and U is shifted onto the accelerometer clock.
     (US-1 run 18: channels 1-7 lead channels 8-17 by 272 samples = 1.36 s in the raw file.)
  2. strong-motion window: Husid D1-99 of F^2. Baselines: medians over [onset-2 s, onset-0.2 s].
  3. identical 4th-order zero-phase low-pass (25 Hz) on U and F, so the pairs stay simultaneous.
  4. per sign: (d_pk, F at d_pk)      zero-velocity point: no viscous term, F = restoring force if SDOF holds
                (d at F_pk, F_pk)      where the capacity is reached
                pp-secant = (F(i+) - F(i-)) / (d(i+) - d(i-)) at the two displacement peaks
                excursion branch: from the last zero crossing before the largest excursion to its peak
  5. K_pk = F(d_pk)/d_pk per sign; on the elastic runs m_eff = K (T_run/2 pi)^2 with T_run from Table 5.
BACKBONE  per sign, pooled branch samples binned on |d|; the per-bin 90th percentile of |F| is the envelope,
          the per-bin median is drawn as the central estimate. Written to loop_backbone_curve.csv.
OUTPUT  <out>/loop_backbone_points.csv, loop_backbone_curve.csv, loop_backbone_branches.npz,
        fig_loop_backbone.png, fig_loop_backbone_runs.png
USAGE   python loop_backbone.py <v2_dir> <out_dir> [--table5 f] [--fig13 f] [--pushover f]
"""
import sys, os, re, glob, csv
import numpy as np
from scipy.signal import butter, filtfilt, correlate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FS_GUESS = 200.0
LP_HZ = 25.0
MAXLAG_S = 4.0
D_EL_MM = 0.5                 # runs whose largest |d| stays below this are treated as elastic for the mass check
D_MIN_MM = 0.10               # ... and whose SMALLER peak exceeds this (below ~0.1 mm the potentiometers stick)
FIG13_FACTOR = 2.317348       # the digitised Fig. 13 is stored at the paper's scale (old x2.317 recipe)
MASSES = {"wall (13.16 kN)": 1340.0, "v2 sum m_j": 1731.6, "model m_eff": 1406.0}


def arg(flag):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else None


def load_v2(path):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    r = np.array([x[:3] for x in ws.iter_rows(min_row=2, values_only=True) if x[0] is not None], float)
    return r[:, 0], r[:, 1] * 1e3, -r[:, 2] / 1e3          # s, mm, kN (restoring sign)


def lowpass(x, fs, fc):
    b, a = butter(4, fc / (0.5 * fs)); return filtfilt(b, a, x)


def bandpass(x, fs):
    b, a = butter(4, [0.5 / (0.5 * fs), 25.0 / (0.5 * fs)], "band"); return filtfilt(b, a, x - np.mean(x))


def sync_lag(u, F, fs):
    """samples by which F lags U (positive: shift U later)."""
    uu, ff = bandpass(u, fs), bandpass(F, fs)
    c = correlate(ff, uu, "full", method="fft"); L = np.arange(-len(uu) + 1, len(uu))
    m = np.abs(L) <= int(MAXLAG_S * fs); k = int(L[m][np.argmax(c[m])])
    c0 = c[L == 0][0] / c[m].max()
    return (k if (abs(k) > 3 and c0 < 0.5) else 0), k, c0


def shift(x, k):
    """x delayed by k samples (edge-padded)."""
    if k == 0: return x
    return np.concatenate([np.full(k, x[0]), x[:-k]]) if k > 0 else np.concatenate([x[-k:], np.full(-k, x[-1])])


def branch(d, F, i_pk, i_lo):
    """from the last zero crossing of d before i_pk (not earlier than i_lo) to i_pk."""
    s = np.sign(d[i_pk]); j = i_pk
    while j > i_lo and s * d[j] > 0: j -= 1
    return d[j:i_pk + 1], F[j:i_pk + 1]


def analyse_run(path):
    t, u, F = load_v2(path)
    fs = 1.0 / np.median(np.diff(t))
    k_used, k_best, c0 = sync_lag(u, F, fs)
    u = shift(u, k_used)
    Fb = bandpass(F, fs); E = np.cumsum(Fb ** 2); E /= E[-1]
    i0, i1 = int(np.searchsorted(E, 0.01)), int(np.searchsorted(E, 0.99))
    pre = slice(max(i0 - int(2 * fs), 0), max(i0 - int(0.2 * fs), 1))
    u0, F0 = np.median(u[pre]), np.median(F[pre])
    d = lowpass(u - u0, fs, LP_HZ); f = lowpass(F - F0, fs, LP_HZ)
    res = {"t": t, "d": d, "F": f, "i0": i0, "i1": i1, "lag_used": k_used, "lag_best": k_best, "c0": c0,
           "fs": fs, "u0": u0, "resid": float(np.median(d[-int(5 * fs):]))}
    w = slice(i0, i1)
    for s, key in ((1, "pos"), (-1, "neg")):
        i_d = i0 + int(np.argmax(s * d[w])); i_F = i0 + int(np.argmax(s * f[w]))
        bd, bF = branch(d, f, i_d, i0)
        res[key] = dict(d_pk=d[i_d], F_at_dpk=f[i_d], d_at_Fpk=d[i_F], F_pk=f[i_F], t_dpk=t[i_d], t_Fpk=t[i_F],
                        K_pk=f[i_d] / d[i_d] if abs(d[i_d]) > 1e-6 else np.nan, bd=bd, bF=bF)
    p, n = res["pos"], res["neg"]
    res["pp_secant"] = (p["F_at_dpk"] - n["F_at_dpk"]) / (p["d_pk"] - n["d_pk"])
    res["dmax"] = max(p["d_pk"], -n["d_pk"])
    return res


def envelope(branches, sign, edges):
    dd = np.concatenate([b[0] for b in branches]) * sign; ff = np.concatenate([b[1] for b in branches]) * sign
    ok = dd > 0; dd, ff = dd[ok], ff[ok]
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (dd >= a) & (dd < b)
        if m.sum() >= 3:
            rows.append((np.sqrt(a * b), np.percentile(ff[m], 90), np.median(ff[m]), int(m.sum())))
    return np.array(rows)


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    v2_dir, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(v2_dir, "Test9Run*_processed_globalzero.xlsx")),
                   key=lambda p: int(re.search(r"Run(\d+)_", p).group(1)))
    T5 = {}
    if arg("--table5"):
        for r in csv.DictReader(open(arg("--table5"))): T5[int(r["run"])] = float(r["T_exp_s"])
    R = {}
    for p in files:
        n = int(re.search(r"Run(\d+)_", p).group(1)); R[n] = analyse_run(p)
        r = R[n]
        msg = "run {:2d}: lag best {:+4d} (rho0/best {:+.2f}) used {:+4d} | +({:6.2f} mm, {:6.2f} kN) -({:7.2f}, {:6.2f}) | Fpk +{:5.2f} -{:5.2f} | pp {:5.2f} kN/mm".format(
            n, r["lag_best"], r["c0"], r["lag_used"], r["pos"]["d_pk"], r["pos"]["F_at_dpk"], r["neg"]["d_pk"], r["neg"]["F_at_dpk"],
            r["pos"]["F_pk"], -r["neg"]["F_pk"], r["pp_secant"])
        print(msg)

    # ---- points table
    with open(os.path.join(out, "loop_backbone_points.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["run", "lag_used_samples", "dmax_mm", "resid_mm", "T_table5_s",
                     "d_pk_pos", "F_at_dpk_pos", "d_at_Fpk_pos", "F_pk_pos", "K_pk_pos",
                     "d_pk_neg", "F_at_dpk_neg", "d_at_Fpk_neg", "F_pk_neg", "K_pk_neg",
                     "pp_secant_kN_per_mm", "m_eff_pp_kg"])
        for n, r in R.items():
            T = T5.get(n, np.nan); p, q = r["pos"], r["neg"]
            wr.writerow([n, r["lag_used"]] + ["{:.4f}".format(x) for x in (
                r["dmax"], r["resid"], T, p["d_pk"], p["F_at_dpk"], p["d_at_Fpk"], p["F_pk"], p["K_pk"],
                q["d_pk"], q["F_at_dpk"], q["d_at_Fpk"], q["F_pk"], q["K_pk"],
                r["pp_secant"], r["pp_secant"] * 1e6 * (T / 2 / np.pi) ** 2)])

    # ---- mass check on the elastic runs
    el = [n for n, r in R.items() if r["dmax"] < D_EL_MM and min(r["pos"]["d_pk"], -r["neg"]["d_pk"]) > D_MIN_MM and n in T5]
    print("\nelastic runs ({} < |d| peaks < {} mm): {}".format(D_MIN_MM, D_EL_MM, el))
    mk = {}
    for lab, f in (("pp-secant", lambda r: r["pp_secant"]), ("K at d_pk, +", lambda r: r["pos"]["K_pk"]),
                   ("K at d_pk, -", lambda r: r["neg"]["K_pk"])):
        K = np.array([f(R[n]) for n in el]); m = np.array([f(R[n]) * 1e6 * (T5[n] / 2 / np.pi) ** 2 for n in el])
        mk[lab] = (K, m)
        print("  {:14s} K = {} kN/mm  -> m_eff = {} kg   (median {:.0f})".format(
            lab, np.round(K, 2), np.round(m).astype(int), np.median(m)))

    # ---- backbone
    edges = np.geomspace(0.05, 40, 40)
    brp = [(R[n]["pos"]["bd"], R[n]["pos"]["bF"]) for n in R]
    brn = [(R[n]["neg"]["bd"], R[n]["neg"]["bF"]) for n in R]
    Ep, En = envelope(brp, 1, edges), envelope(brn, -1, edges)
    with open(os.path.join(out, "loop_backbone_curve.csv"), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["sign", "d_mm", "F_p90_kN", "F_median_kN", "n_samples"])
        for lab, E in (("pos", Ep), ("neg", En)):
            for row in E: wr.writerow([lab] + ["{:.4f}".format(x) for x in row[:3]] + [int(row[3])])
    np.savez(os.path.join(out, "loop_backbone_branches.npz"),
             **{"run{}_{}_{}".format(n, s, k): R[n][s][k] for n in R for s in ("pos", "neg") for k in ("bd", "bF")})

    # ---- figure 1: backbone
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    cm = plt.cm.viridis
    for n, r in R.items():
        c = cm((n - 1) / 23)
        for s in ("pos", "neg"):
            ax[0].plot(r[s]["bd"], r[s]["bF"], lw=0.6, color=c, alpha=0.8)
        ax[0].plot([r["pos"]["d_pk"], r["neg"]["d_pk"]], [r["pos"]["F_at_dpk"], r["neg"]["F_at_dpk"]], "o", ms=3.5, color=c, mec="k", mew=0.3)
        ax[0].plot([r["pos"]["d_at_Fpk"], r["neg"]["d_at_Fpk"]], [r["pos"]["F_pk"], r["neg"]["F_pk"]], "^", ms=3.5, color=c, mec="k", mew=0.3)
    ax[0].plot(Ep[:, 0], Ep[:, 1], "k-", lw=2, label="branch envelope (p90)")
    ax[0].plot(-En[:, 0], -En[:, 1], "k-", lw=2)
    if arg("--fig13"):
        g = np.genfromtxt(arg("--fig13"), delimiter=",", names=True)
        ax[0].plot(g["d_mm"], g["F_kN"] / FIG13_FACTOR, "--", color="#c0392b", lw=1.4,
                   label="Fig. 13 envelope / {:.3f} (peaks, non-simultaneous)".format(FIG13_FACTOR))
    if arg("--pushover"):
        g = np.genfromtxt(arg("--pushover"), delimiter=",", names=True)
        ax[0].plot(g["d_mm"], g["F_kN"], "-", color="#2b5f8a", lw=1.6, label="model pushover")
    ax[0].axhline(0, color="k", lw=0.5); ax[0].axvline(0, color="k", lw=0.5); ax[0].grid(alpha=0.25)
    ax[0].set_xlabel("d at 2.06 m, relative to table (mm)"); ax[0].set_ylabel("F = -sum m_j a_j (kN)")
    ax[0].set_title("Excursion branches, zero-velocity points (o) and peak-force points (^)", loc="left", fontsize=10)
    ax[0].legend(fontsize=8, loc="lower right")
    # log-log, both signs folded
    for lab, E, c in (("+", Ep, "#1f6f8b"), ("-", En, "#b4491a")):
        ax[1].loglog(E[:, 0], E[:, 1], "-", color=c, lw=2, label="envelope p90, sign " + lab)
        ax[1].loglog(E[:, 0], E[:, 2], ":", color=c, lw=1.2, label="median, sign " + lab)
    for n, r in R.items():
        for s, mk_ in (("pos", "o"), ("neg", "s")):
            dd, ff = abs(r[s]["d_pk"]), abs(r[s]["F_at_dpk"])
            if dd > 0.02 and ff > 0.05: ax[1].loglog(dd, ff, mk_, ms=3.5, color=cm((n - 1) / 23), mec="k", mew=0.3)
    dd = np.geomspace(0.05, 40)
    for K, ls in ((3.2, "--"), (5.1, "-."), (7.6, "-")):
        ax[1].loglog(dd, K * dd, ls, color="grey", lw=0.8)
        ax[1].text(0.06, K * 0.06 * 1.15, "{} kN/mm".format(K), fontsize=7, color="grey")
    ax[1].set_xlim(0.05, 40); ax[1].set_ylim(0.1, 20); ax[1].grid(alpha=0.25, which="both")
    ax[1].set_xlabel("|d| (mm)"); ax[1].set_ylabel("|F| (kN)"); ax[1].legend(fontsize=7, loc="lower right")
    ax[1].set_title("Folded, log-log: 3.2 = Fig. 13 slope, 5.1 = model K1, 7.6 = loop pp-secant", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_loop_backbone.png"), dpi=150); plt.close(fig)

    # ---- figure 2: per-run loops
    fig, axs = plt.subplots(4, 6, figsize=(20, 12)); axs = axs.ravel()
    for a, (n, r) in zip(axs, R.items()):
        w = slice(r["i0"], r["i1"])
        a.plot(r["d"][w], r["F"][w], lw=0.4, color="#8a8f97")
        for s in ("pos", "neg"):
            a.plot(r[s]["bd"], r[s]["bF"], lw=1.2, color="#b4491a")
            a.plot(r[s]["d_pk"], r[s]["F_at_dpk"], "o", ms=4, color="k")
            a.plot(r[s]["d_at_Fpk"], r[s]["F_pk"], "^", ms=4, color="#1f6f8b")
        a.axhline(0, color="k", lw=0.4); a.axvline(0, color="k", lw=0.4)
        a.set_title("run {}{}".format(n, "  (U shifted {} smp)".format(r["lag_used"]) if r["lag_used"] else ""), fontsize=9, loc="left")
        a.tick_params(labelsize=7)
    fig.supxlabel("d (mm)"); fig.supylabel("F (kN)"); fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_loop_backbone_runs.png"), dpi=110); plt.close(fig)
    print("->", out)
    return R, Ep, En, mk


if __name__ == "__main__":
    main()
