# -*- coding: utf-8 -*-
"""compare_stratF_old.py -- the older Strategy F folder stratF_full_results_US1 against the current z-up runs and the tests.

Run from the Hanze_SAFEGO root:  python compare_stratF_old.py
Reads the metrics written by postprocess_stratF_full25.py (same peak definitions and the same group force formula for all):
  stratF_postproc_US1folder/stratF_metrics.csv         stratF_full_results_US1 (y-up build, Part_I_MASON_v8.sav)
  stratF_postproc_pk2_ft0p05_5vs1p5/stratF_metrics.csv  z-up pk2 f_t 0.05 glued pocket, 1.5 % Maxwell
  stratF_postproc_pK11_1p5/stratF_metrics.csv           z-up pk2 f_t 0.05 pK11 pocket, 1.5 % Maxwell
Output: stratF_old_compare/fig_stratF_old.png, stratF_old_compare.csv. Linear axes.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "stratF_old_compare"
CASES = [("old", "stratF_postproc_US1folder/stratF_metrics.csv", "US1 folder, 31 Aug",
          dict(color="#9333EA", lw=2.2, alpha=0.6, marker="o", ms=4, label="stratF_full_results_US1 (y-up build)")),
         ("glued", "stratF_postproc_pk2_ft0p05_5vs1p5/stratF_metrics.csv", "pk2 f t 0.05, 1.5 %",
          dict(color="#0F766E", lw=2.2, alpha=0.6, marker="s", ms=4, label="z-up glued pocket, 1.5 %")),
         ("pK11", "stratF_postproc_pK11_1p5/stratF_metrics.csv", "pK11 pocket, 1.5 %",
          dict(color="#B91C1C", lw=2.2, alpha=0.6, marker="D", ms=4, label="z-up pK11 pocket, 1.5 %"))]
TESTS = dict(US1=dict(color="#111827", lw=1.0, marker="o", mfc="white", ms=5, label="US-1"),
             US2=dict(color="#1D4ED8", lw=1.0, ls="--", marker="^", mfc="white", ms=5, label="US-2"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def metrics(path, variant):
    with open(path) as f:
        return {int(r["run"]): r for r in csv.DictReader(f) if r["variant"] == variant}


def main():
    os.makedirs(OUT, exist_ok=True)
    M = {k: metrics(p, v) for k, p, v, _ in CASES}
    ref = M["pK11"]
    runs = list(range(1, 26))
    rows = []
    for n in runs:
        r = dict(run=n, d_US1=fnum(ref.get(n, {}).get("US1_peak_mm")), d_US2=fnum(ref.get(n, {}).get("US2_peak_mm")),
                 F_US1=fnum(ref.get(n, {}).get("US1_peak_F_kN")), F_US2=fnum(ref.get(n, {}).get("US2_peak_F_kN")))
        for k, _, _, _ in CASES:
            r["d_" + k] = fnum(M[k].get(n, {}).get("peak_incr_mm"))
            r["F_" + k] = fnum(M[k].get(n, {}).get("peak_F_kN"))
        rows.append(r)
    with open(os.path.join(OUT, "stratF_old_compare.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    x = np.array(runs)
    A = lambda k: np.array([r[k] for r in rows], float)
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.2))
    for a, q, yl in ((ax[0], "d", "peak |d| at 2.06 m, per run (mm)"), (ax[1], "F", "peak |F|, group formula (kN)")):
        for k, _, _, st in CASES:
            a.plot(x, A(q + "_" + k), **st)
        for t, st in TESTS.items():
            a.plot(x, A(q + "_" + t), **st)
        a.set_ylabel(yl); a.set_xlabel("run"); a.set_xticks(runs); a.set_xlim(0.5, 25.5); a.set_ylim(bottom=0)
        a.grid(color="#E8E8E8", lw=0.7); a.axvline(21.5, color="#9CA3AF", lw=0.8)
    ax[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Strategy F: stratF_full_results_US1 (y-up build) against the z-up pk2 f_t 0.05 runs (1.5 % Maxwell) "
                 "and US-1 / US-2", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_stratF_old.png"), dpi=150)

    def gm(a, b, sel):
        a, b = a[sel], b[sel]; k = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        return float(np.exp(np.mean(np.log(a[k] / b[k])))) if k.any() else np.nan
    for lab, sel in (("runs 1-14", x <= 14), ("runs 15-21", (x >= 15) & (x <= 21)), ("runs 22-24", (x >= 22) & (x <= 24))):
        print(lab, " | ".join("{} d/US1 {:.2f} F/US1 {:.2f} F/US2 {:.2f}".format(
            k, gm(A("d_" + k), A("d_US1"), sel), gm(A("F_" + k), A("F_US1"), sel), gm(A("F_" + k), A("F_US2"), sel))
            for k, _, _, _ in CASES))
    print("-> " + OUT)


if __name__ == "__main__":
    main()
