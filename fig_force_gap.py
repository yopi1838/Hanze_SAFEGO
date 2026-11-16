# -*- coding: utf-8 -*-
"""fig_force_gap.py -- is the model/test peak-force gap a constant offset or a factor? (16 Nov 2026)
Reads stratF_postproc_unclamp_vs_glued_mean/stratF_metrics.csv (postprocess_stratF_full25.py --tests mean).
Panels: as computed | test mean + c (c = mean gap over runs 1-21) | test mean x k (least-squares k over runs 1-21).
Linear axes. Output: fig_force_gap.png next to the CSV."""
import csv, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
SRC = sys.argv[1] if len(sys.argv) > 1 else "stratF_postproc_unclamp_vs_glued_mean"
VARS = [("unclamped pocket, 1.5 %", "#C2410C"), ("glued pocket, 1.5 %", "#0F766E")]
R = list(csv.DictReader(open(os.path.join(SRC, "stratF_metrics.csv"))))
def get(v):
    rr = [r for r in R if r["variant"] == v and r["ref_peak_F_kN"]]
    return (np.array([int(r["run"]) for r in rr]), np.array([float(r["peak_F_kN"]) for r in rr]),
            np.array([float(r["ref_peak_F_kN"]) for r in rr]))
n, m0, t = get(VARS[0][0]); s = n <= 21
c = np.mean(np.concatenate([get(v)[1][get(v)[0] <= 21] - t[s] for v, _ in VARS]))
k = np.mean([np.sum(get(v)[1][s] * t[s]) / np.sum(t[s] ** 2) for v, _ in VARS])
fig, ax = plt.subplots(3, 1, figsize=(12, 10.5), sharex=True, dpi=160)
for a, (tt, ttl) in zip(ax, ((t, "As computed: test mean of US-1 and US-2"),
                             (t + c, "Test mean shifted up by a constant {:.2f} kN (mean gap, runs 1-21)".format(c)),
                             (t * k, "Test mean scaled by {:.2f} (least squares, runs 1-21)".format(k)))):
    for v, col in VARS:
        nn, mm, _ = get(v)
        a.plot(nn, mm, "-s", color=col, lw=2.2, ms=5, mec="none", label=v)
    a.plot(n, tt, "-D", color="#111827", mfc="white", lw=1.3, ms=4.5, label="test (mean), adjusted")
    rms = np.sqrt(np.mean([(get(v)[1][s] - tt[s]) ** 2 for v, _ in VARS]))
    a.set_title(ttl + "   |   rms misfit runs 1-21: {:.2f} kN".format(rms), loc="left", fontweight="bold", fontsize=10)
    a.axvspan(21.5, 25.5, color="#F3E6D8", alpha=0.55, lw=0); a.grid(color="#E8E8E8", lw=0.7); a.set_axisbelow(True)
    a.set_ylabel("peak |F| (kN)"); a.set_xlim(0.5, 25.5); a.set_xticks(range(1, 26))
ax[0].legend(frameon=False, fontsize=8, loc="upper left"); ax[2].set_xlabel("run")
fig.tight_layout(); p = os.path.join(SRC, "fig_force_gap.png"); fig.savefig(p, bbox_inches="tight"); print("->", p, "c", c, "k", k)
