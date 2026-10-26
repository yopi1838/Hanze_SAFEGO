"""ringdown_cycles.py -- amplitude-dependent period of US-1 / US-2 from the post-shaking ring-down, cycle by cycle,
from EXP_DATA/processed_globalzero (Time, U_avg, base_shear_N).

SIGNAL  Once the table has stopped, base_shear = sum m_j a_j is the wall's own acceleration (times a mass); its
        frequency content does not depend on the masses, the 2.317 factor or the table term. The potentiometers
        stick below ~0.1 mm, so the force is the ring-down signal; U_avg is read alongside where it moves.
TABLE STOP  the table term shows in the force below 4 Hz (low-pass 4 Hz, Hilbert envelope; on US-1 it tracks ch12 of
        the raw files within 0.03 s). Cycles are read from the END OF STRONG MOTION = last time that envelope
        exceeds 5 % of its peak, +0.1 s. The table then crawls at the ~1 % level for 10-15 s more; that true stop
        (1 % threshold) is recorded as t_table_true_stop, the paper's 45-60 s window starts there.
CYCLES  force band-passed 3-30 Hz (zero-phase), read from the main pulse onward. A cycle is FREE once the table term
        (low-pass 4 Hz) is below 5 % of its peak; earlier cycles are forced and are listed but not used in the bins.
        Period of a cycle = time between successive upward zero crossings,
        with linear interpolation of the crossing instant; amplitude = half the peak-to-peak force within the cycle.
        Cycles are kept while the amplitude exceeds 10x the pre-event noise of the band-passed force and the period
        lies in 0.05-0.25 s. Equivalent displacement amplitude d_eq = F_amp / (m w^2) with m = 1133 kg (the files'
        mass without the factor; a scale, not a measurement) -- and the measured U amplitude where > 0.1 mm.
OUTPUT  <out>/ringdown_cycles_<spec>.csv (every cycle), ringdown_summary.csv (per run: T of the first 3 cycles,
        T at the smallest kept amplitude, number of cycles, table-stop time), fig_ringdown_cycles.png
USAGE   python ringdown_cycles.py <proc_dir> <out_dir> [--fig14 us1.csv us2.csv]
"""
import sys, os, re, glob, csv
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_SCALE_KG = 1133.3
TESTS = (("US1", "Test9", 0.091), ("US2", "Test12", 0.087))


def load(path):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True); hdr = list(next(it)); c = {h: j for j, h in enumerate(hdr)}
    rows = [r for r in it if r[c["Time"]] is not None]
    t = np.array([r[c["Time"]] for r in rows], float); U = np.array([r[c["U_avg"]] for r in rows], float) * 1e3
    F = -np.array([r[c["base_shear_kN"]] for r in rows], float) / float(rows[0][c["mass_scale"]])
    return t, U, F


def filt(x, fs, lo=None, hi=None):
    if lo is None: b, a = butter(4, hi / (0.5 * fs))
    elif hi is None: b, a = butter(4, lo / (0.5 * fs), "high")
    else: b, a = butter(4, [lo / (0.5 * fs), hi / (0.5 * fs)], "band")
    return filtfilt(b, a, x)


def analyse(path):
    t, U, F = load(path); fs = 1.0 / np.median(np.diff(t)); F = F - np.median(F[: int(2 * fs)])
    fb = filt(F, fs, 3.0, 30.0); fl = filt(F, fs, hi=4.0)
    env_l = np.abs(hilbert(fl)); pre = slice(int(0.2 * fs), int(2.0 * fs))
    noise_l = np.std(fl[pre]); noise_b = np.std(fb[pre])
    # end of strong motion: last time the table term exceeds 5 % of its peak (+0.1 s). The table keeps crawling at
    # ~1 % level for 10-15 s after that (true stop, also recorded); the 10 Hz ring-down has decayed long before.
    on = np.where(env_l > 0.05 * env_l.max())[0]; i_stop = int(on[-1] + 0.1 * fs)
    t_true = t[np.where(env_l > 0.01 * env_l.max())[0][-1]]
    # cycles from the main pulse onward; a cycle is "free" once the table term is below 5 % of its peak
    i_start = int(np.argmax(np.abs(fb))); x = fb[i_start:]; tt = t[i_start:]
    k = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    tz = tt[k] - x[k] * (tt[k + 1] - tt[k]) / (x[k + 1] - x[k])
    cycles = []
    for a, b, ta, tb in zip(k[:-1], k[1:], tz[:-1], tz[1:]):
        T = tb - ta; seg = x[a:b + 1]; amp = 0.5 * (seg.max() - seg.min())
        if not (0.05 <= T <= 0.25): continue
        if amp < 10 * noise_b: break
        useg = U[i_start + a:i_start + b + 1]; uamp = 0.5 * (useg.max() - useg.min())
        w = 2 * np.pi / T; free = float(env_l[i_start + a:i_start + b + 1].max() < 0.05 * env_l.max())
        cycles.append((ta - t[i_stop], T, amp, amp * 1e3 / (M_SCALE_KG * w * w) * 1e3, uamp, free))
    return dict(t_stop=t[i_stop], t_true=t_true, cycles=np.array(cycles), noise_b=noise_b)


def main():
    proc, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
    fig14 = {}
    if "--fig14" in sys.argv:
        j = sys.argv.index("--fig14")
        for spec, p in (("US1", sys.argv[j + 1]), ("US2", sys.argv[j + 2])):
            g = np.genfromtxt(p, delimiter=",", names=True); fig14[spec] = dict(zip(g["run"].astype(int), g["T_exp_s"]))
    summary = []; allc = {}; allbin = {}; BINS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
    for spec, test, T1 in TESTS:
        files = sorted(glob.glob(os.path.join(proc, test + "Run*_processed_globalzero.xlsx")),
                       key=lambda p: int(re.search(r"Run(\d+)_", p).group(1)))
        print("== {} ({})".format(spec, test))
        with open(os.path.join(out, "ringdown_cycles_{}.csv".format(spec)), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["run", "t_after_strong_motion_end_s", "T_s", "F_amp_kN", "d_eq_mm_m1133", "U_amp_mm", "free"])
            for p in files:
                n = int(re.search(r"Run(\d+)_", p).group(1)); r = analyse(p); c = r["cycles"]; allc[(spec, n)] = c
                for row in c: w.writerow([n] + ["{:.5f}".format(v) for v in row])
                fr = c[c[:, 5] > 0] if len(c) else c
                BINS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
                binned = []
                for lo, hi in zip(BINS[:-1], BINS[1:]):
                    m = (fr[:, 3] >= lo) & (fr[:, 3] < hi) if len(fr) else np.zeros(0, bool)
                    binned.append(float(np.median(fr[m, 1])) if m.sum() >= 2 else np.nan)
                allbin[(spec, n)] = binned
                if len(c) == 0:
                    print("  run {:2d}: no ring-down cycles above noise"); continue
                T_first = np.mean(fr[: min(3, len(fr)), 1]) if len(fr) else np.nan; T_last = np.mean(fr[-min(3, len(fr)):, 1]) if len(fr) else np.nan
                summary.append([spec, n, r["t_stop"], r["t_true"], len(c), len(fr), c[0, 3], c[0, 4], fr[0, 3] if len(fr) else np.nan, T_first, T_last, fig14.get(spec, {}).get(n, np.nan)] + binned)
                print("  run {:2d}: {:3d} cycles ({:3d} free) | largest d_eq {:5.3f} mm (U {:5.3f}) | first free d_eq {:6.4f} T {:.4f} | bins {} | paper {}".format(
                    n, len(c), len(fr), c[0, 3], c[0, 4], fr[0, 3] if len(fr) else np.nan, T_first,
                    " ".join("{:.3f}".format(v) if np.isfinite(v) else "  -  " for v in binned),
                    "{:.4f}".format(fig14[spec][n]) if spec in fig14 and n in fig14[spec] else "-"))
    with open(os.path.join(out, "ringdown_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["specimen", "run", "t_strong_motion_end_s", "t_table_true_stop_s", "n_cycles", "n_free", "d_eq_largest_mm", "U_amp_largest_mm",
                    "d_eq_first_free_mm", "T_first3free_s", "T_last3free_s", "T_paper_s"] + ["T_deq_{}_{}".format(lo, hi) for lo, hi in zip(BINS[:-1], BINS[1:])])
        for r in summary: w.writerow(r[:2] + ["{:.5f}".format(v) for v in r[2:]])
    # figure
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    for col, (spec, test, T1) in enumerate(TESTS):
        runs = sorted(n for s_, n in allc if s_ == spec); cm = plt.cm.viridis
        a = ax[0, col]
        for n in runs:
            c = allc[(spec, n)]
            if not len(c): continue
            f_ = c[:, 5] > 0
            a.semilogx(c[f_, 3], c[f_, 1], ".", ms=4, color=cm((n - 1) / 24), alpha=0.8)
            a.semilogx(c[~f_, 3], c[~f_, 1], "o", ms=4, mfc="none", color=cm((n - 1) / 24), alpha=0.5)
        a.axhline(T1, color="grey", ls=":", lw=1); a.text(0.0012, T1 + 0.001, "Table 5 initial {:.3f} s".format(T1), fontsize=7, color="grey")
        a.set_xlim(0.001, 3); a.set_ylim(0.05, 0.25)
        a.set_xlabel("equivalent amplitude F/(m w^2), m = 1133 kg  (mm)"); a.set_ylabel("period of the cycle (s)")
        a.set_title("{}: cycle periods, filled = free, hollow = table still moving; colour = run 1 (dark) .. {} (light)".format(spec, runs[-1]), loc="left", fontsize=9)
        a.grid(alpha=0.25, which="both")
        a = ax[1, col]
        S = [r for r in summary if r[0] == spec]; rn = [r[1] for r in S]
        for j, (lo, hi) in enumerate(zip(BINS[:-1], BINS[1:])):
            a.plot(rn, [r[12 + j] for r in S], "o-", ms=3.5, lw=0.9, color=plt.cm.plasma(j / 5), label="free cycles, d_eq {}-{} mm".format(lo, hi))
        if spec in fig14: a.plot(rn, [r[11] for r in S], "^--", color="k", ms=4, label="paper Fig. 14 / Table 5 (PSD of the quiet 45-60 s)")
        a.set_xlabel("run"); a.set_ylabel("period (s)"); a.grid(alpha=0.25); a.legend(fontsize=7, ncol=2)
        a.set_title("{}: period per run, binned by amplitude".format(spec), loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_ringdown_cycles.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
