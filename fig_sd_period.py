# -*- coding: ascii -*-
"""
Paper Fig-4-style figure per model: adaptive Sd target + amplification
(top) and fundamental period evolution NUM vs US-1/US-2 PSD periods
(bottom, with linear trends). Periods from exp_TestN_period_psd.csv
(PSD/Welch method, monotonicity-corrected values).

Usage: python fig_sd_period.py
Outputs: figP2_nodamp.png, figP2_ratcheting.png in postproc/.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PP = HERE / "stratC_results_NODAMP_v6_NEW" / "postproc"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.grid": True, "grid.linestyle": ":", "grid.alpha": 0.6,
    "legend.frameon": True, "legend.framealpha": 1.0,
    "legend.edgecolor": "0.7"})

def fread(p):
    with open(str(p)) as f:
        return list(csv.DictReader(f))

def col(rows, k, cast=float):
    return np.array([cast(r[k]) for r in rows if r.get(k) not in ("", None)])

exp9 = fread(PP / "exp_Test9_period_psd.csv")
exp12 = fread(PP / "exp_Test12_period_psd.csv")
r9 = np.array([int(r["run"]) for r in exp9 if r["T1_best_s"]])
T9 = np.array([float(r["T1_best_s"]) for r in exp9 if r["T1_best_s"]])
r12 = np.array([int(r["run"]) for r in exp12 if r["T1_best_s"]])
T12 = np.array([float(r["T1_best_s"]) for r in exp12 if r["T1_best_s"]])

CASES = [("figP2_nodamp.png", "NUM: undamped",
          HERE / "stratC_results_NODAMP_v6_NEW"),
         ("figP2_ratcheting.png", "NUM: 3% Rayleigh",
          HERE / "stratC_results_RATCHETING")]

for fname, lbl, simdir in CASES:
    sim = fread(simdir / "strategy_C_summary.csv")
    runs = col(sim, "run")
    sd = col(sim, "Sd_target_mm")
    amp = col(sim, "amplification")
    t_end = col(sim, "T_end")
    ratio = col(sim, "T_end_over_Tinit")
    fallback = ratio == 1.0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.2, 8.2), sharex=True)
    axb = ax1.twinx()
    axb.bar(runs, amp, color="lightcoral", alpha=0.45, width=0.62, zorder=1)
    axb.set_ylabel("Amplification factor", color="red", fontsize=14)
    axb.tick_params(axis="y", colors="red")
    axb.set_ylim(0, max(3.0, amp.max() * 1.05))
    axb.axhline(1.0, color="red", ls=":", lw=0.8, alpha=0.5)
    axb.grid(False)
    ax1.plot(runs, sd, "s-", color="black", ms=6, lw=1.5,
             label="Adaptive $S_d$ target")
    ax1.bar([np.nan], [np.nan], color="lightcoral", alpha=0.45,
            label="Amplification factor")
    ax1.set_ylabel("$S_d$ target (mm)", fontsize=14)
    ax1.set_title(lbl, fontsize=13)
    ax1.legend(fontsize=11, loc="upper left")
    ax1.set_zorder(axb.get_zorder() + 1)
    ax1.patch.set_visible(False)

    # bottom: periods + linear trends (paper style)
    ax2.plot(runs[~fallback], t_end[~fallback], "s", color="red", ms=6,
             label="NUM")
    if fallback.any():
        ax2.plot(runs[fallback], t_end[fallback], "s", ms=7, mfc="white",
                 mec="red", mew=1.4, label="NUM (period-ID fallback)")
    ax2.plot(r9, T9, "o", color="royalblue", ms=6, label="US-1 (Test 9)")
    ax2.plot(r12, T12, "o", color="orange", ms=6, label="US-2 (Test 12)")
    for x, y, c in ((runs[~fallback], t_end[~fallback], "red"),
                    (r9, T9, "royalblue"), (r12, T12, "orange")):
        if len(x) > 2:
            cf = np.polyfit(x, y, 1)
            ax2.plot(x, np.polyval(cf, x), "-", color=c, lw=1.4)
    ax2.set_xlabel("Run number", fontsize=14)
    ax2.set_ylabel("Fundamental period (s)", fontsize=14)
    ax2.set_xlim(0, 26)
    ax2.set_xticks(range(0, 27, 2))
    ymax = max(0.2, float(t_end.max()) * 1.1)
    ax2.set_ylim(0.025, ymax)
    ax2.legend(fontsize=10, loc="upper left")
    fig.tight_layout()
    fig.savefig(str(PP / fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("->", fname)
