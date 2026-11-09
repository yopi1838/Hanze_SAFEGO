# -*- coding: utf-8 -*-
"""compare_pocket_pushover.py -- glued vs timber-mortar (pK11) pocket pushovers of the pk2 f_t 0.05 build, against the
test envelope (US-1, US-2: force at the displacement peaks of every run, 10 Hz low-pass, exp_processed_v3). Linear axes.
Run from the Hanze_SAFEGO root:  python compare_pocket_pushover.py
Output: pushover_pocket_compare/fig_pocket_pushover.png, pocket_pushover_metrics.csv, test_envelope.csv"""
import os, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASES = [("glued pocket (brick-mortar set)", "pushover_results_v5zup_pk2_ft0p05_TopMason", "#0F766E"),
         ("pK11 pocket (timber-mortar set)", "pushover_results_v5zup_pk2_ft0p05_pK11_TopMason", "#B91C1C")]
EXP = "exp_processed_v3"
OUT = "pushover_pocket_compare"
LP_HZ = 10.0


def read(path):
    lines = open(path).read().splitlines()
    mt = float(re.search(r"m_tot\s+([0-9.]+)", lines[0]).group(1))
    rows = list(csv.DictReader(l for l in lines[1:] if l.strip()))
    g = lambda k: np.array([float(r[k]) for r in rows])
    return dict(d=np.abs(g("d_ctrl_mm")), Fa=np.abs(g("a_mps2")) * mt / 1e3, Ft=np.abs(g("F_tot_kN")),
                pocket=g("pocket_open_mm"), jN=g("joist_N_kN"), jV=g("joist_Vz_kN"), phase=[r["phase"] for r in rows])


def metrics(c):
    d, F = c["d"], c["Fa"]
    i = int(np.argmax(F)); Fp, dp = F[i], d[i]
    k_in = np.polyfit(d[(d > 0) & (d <= 0.6)], F[(d > 0) & (d <= 0.6)], 1)[0] if np.sum((d > 0) & (d <= 0.6)) >= 2 else F[0] / d[0]
    post = np.where((d > dp) & (F <= 0.8 * Fp))[0]
    du = d[post[0]] if len(post) else np.nan
    at = lambda x: float(np.interp(x, d, F))
    po = lambda x: float(np.interp(x, d, c["pocket"]))
    return dict(K_init_kN_per_mm=round(k_in, 2), F_peak_kN=round(Fp, 2), d_at_peak_mm=round(dp, 2), d_u80_mm=round(du, 1),
                F_5mm=round(at(5), 2), F_10mm=round(at(10), 2), F_20mm=round(at(20), 2), F_30mm=round(at(30), 2),
                pocket_open_at_peak_mm=round(po(dp), 2), pocket_open_20mm=round(po(20), 2), d_end_mm=round(d[-1], 1))


def test_env():
    import pandas as pd
    from scipy.signal import butter, filtfilt
    pts = []
    for test, nmax in (("Test9", 24), ("Test12", 25)):
        for n in range(1, nmax + 1):
            f = os.path.join(EXP, "{}Run{}_processed_globalzero.xlsx".format(test, n))
            if not os.path.isfile(f):
                continue
            df = pd.read_excel(f); t = df.Time.values; fs = 1 / np.median(np.diff(t))
            d = df.U_avg.values * 1e3; d = d - d[0]
            b, a = butter(4, LP_HZ / (0.5 * fs)); F = filtfilt(b, a, df.base_shear_kN.values)
            F = F * np.sign(np.polyfit(d, F, 1)[0])
            w = int(0.05 * fs)
            for s, i in ((+1, int(np.argmax(d))), (-1, int(np.argmin(d)))):
                seg = F[max(0, i - w): i + w]
                pts.append((("US-1" if test == "Test9" else "US-2"), n, s, abs(d[i]), float(np.max(s * seg))))
    return pts


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {}
    rows = []
    for lab, fo, col in CASES:
        for side in ("pos", "neg"):
            p = os.path.join(fo, "pushover_{}.csv".format(side))
            if os.path.isfile(p):
                data[(lab, side)] = read(p)
                rows.append(dict(case=lab, side=side, **metrics(data[(lab, side)])))
    with open(os.path.join(OUT, "pocket_pushover_metrics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    pts = test_env()
    with open(os.path.join(OUT, "test_envelope.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["test", "run", "sign", "abs_d_peak_mm", "F_at_d_peak_kN_10Hz"]); w.writerows(pts)

    fig, ax = plt.subplots(2, 3, figsize=(18, 10.5))
    tst = {"US-1": dict(marker="o", mfc="white", mec="#111827", ls="none", ms=6),
           "US-2": dict(marker="^", mfc="white", mec="#1D4ED8", ls="none", ms=6)}
    for j, (side, sgn) in enumerate((("pos", +1), ("neg", -1))):
        for k, xmax in enumerate((40.0, 3.0)):
            a = ax[k, j]
            for lab, fo, col in CASES:
                c = data.get((lab, side))
                if c is None:
                    continue
                a.plot(c["d"], c["Fa"], color=col, lw=2.0, label=lab + " (a m_tot)")
                a.plot(c["d"], c["Ft"], color=col, lw=0.9, ls=":", label=lab + " (reactions)")
            for tname, st in tst.items():
                P = [(dd, ff) for (tn, n, s, dd, ff) in pts if tn == tname and s == sgn and dd <= xmax]
                if P:
                    a.plot(*zip(*P), label=tname + ": F at the run's peak d (10 Hz)", **st)
            a.set_xlim(0, xmax); a.set_ylim(bottom=0); a.grid(alpha=0.25)
            a.set_xlabel("|d| at 2.06 m (mm)"); a.set_ylabel("|F| (kN)")
            a.set_title("{} push{}".format("Positive" if sgn > 0 else "Negative", "" if k == 0 else ", first 3 mm"), loc="left")
            if j == 0 and k == 0:
                a.legend(fontsize=8)
    for lab, fo, col in CASES:
        for side, ls in (("pos", "-"), ("neg", "--")):
            c = data.get((lab, side))
            if c is not None:
                ax[0, 2].plot(c["d"], c["pocket"], color=col, ls=ls, lw=1.6, label="{}, {}".format(lab.split(" (")[0], side))
                ax[1, 2].plot(c["d"], c["jN"], color=col, ls=ls, lw=1.6, label="{}, {}".format(lab.split(" (")[0], side))
    ax[0, 2].set_title("Largest joist-masonry contact opening", loc="left"); ax[0, 2].set_ylabel("pocket_open_mm (mm)")
    ax[1, 2].set_title("Joist support normal force", loc="left"); ax[1, 2].set_ylabel("joist_N (kN)")
    for a in ax[:, 2]:
        a.set_xlabel("|d| at 2.06 m (mm)"); a.grid(alpha=0.25); a.legend(fontsize=8); a.set_xlim(0, 40)
    ax[0, 2].set_ylim(bottom=0)
    fig.suptitle("pk2 f_t 0.05 build: glued pocket vs timber-mortar pocket (pK11), pushovers against the test envelope (linear axes)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_pocket_pushover.png"), dpi=130)
    print("-> " + OUT)


if __name__ == "__main__":
    main()
