# -*- coding: ascii -*-
"""
Presentation figures: experiment (US-1, Test 9) vs numerical prediction
(stratC_results_NODAMP_v6_NEW), styled after the IMC2026 conference paper
(Figs. 4 and 5: serif fonts, red squares = simulation, blue circles = US-1).

  figP1_peak_disp.png   paper Fig. 5 style: peak rel OOP disp per run
  figP2_sd_period.png   paper Fig. 4 style: Sd target + amplification (top),
                        fundamental period evolution (bottom)
  figP3_residual.png    residual displacement and tilt (not in paper)

Usage: python make_presentation_figs.py [SIM_DIR] [OUT_DIR]
"""
import csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

SIM = Path(sys.argv[1] if len(sys.argv) > 1 else "stratC_results_NODAMP_v6_NEW")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "presentation_figs")
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix",
    "font.size": 13, "axes.labelsize": 15, "axes.titlesize": 15,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "axes.grid": True, "grid.linestyle": ":", "grid.alpha": 0.6,
    "legend.frameon": True, "legend.framealpha": 1.0,
    "legend.edgecolor": "0.7",
})
C_SIM, C_US1 = "red", "royalblue"

# ---------------------------------------------------------------- load data
def load_csv(path, key="run"):
    return {int(float(r[key])): r for r in csv.DictReader(open(str(path)))}

sim_sum = load_csv(SIM / "postproc" / "postproc_summary.csv")
drv     = load_csv(SIM / "strategy_C_summary.csv")
exp     = load_csv(SIM / "postproc" / "exp_Test9_metrics.csv")

# residuals from the all-channels table
sim_res_disp, sim_res_tilt = {}, {}
for r in csv.DictReader(open(str(SIM / "postproc" / "postproc_all_channels.csv"))):
    rn = int(r["run"])
    if r["channel"] == "rel_disp_top_mm":
        sim_res_disp[rn] = float(r["residual"])
    elif r["channel"] == "tilt_full_wall":
        sim_res_tilt[rn] = float(r["residual"])

runs_sim = sorted(sim_sum)
runs_exp = sorted(exp)
sim_peak = np.array([float(sim_sum[r]["peak_rel_disp_mm"]) for r in runs_sim])
exp_peak = np.array([float(exp[r]["peak_rel_mm"]) for r in runs_exp])

T_INIT = 0.092


def annotate_offscale(ax, x, y, ymax, color):
    """Clip a series at ymax; mark clipped points with arrows + labels."""
    for xi, yi in zip(x, y):
        if yi > ymax:
            ax.annotate("%.0f mm" % yi, xy=(xi, ymax * 0.97),
                        xytext=(xi, ymax * 0.80), ha="center", fontsize=11,
                        color=color,
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.2))


# ------------------------------------------------- figP1: peak displacement
YMAX1 = 35.0
fig, ax = plt.subplots(figsize=(11, 5.8))
ax.plot(runs_exp, exp_peak, "o--", color=C_US1, ms=7, lw=1.4,
        mec="darkblue", label="US-1 (Test 9)")
m = sim_peak <= YMAX1
ax.plot(np.array(runs_sim)[m], sim_peak[m], "s-", color=C_SIM, ms=7, lw=1.8,
        label="Simulation", zorder=3)
annotate_offscale(ax, runs_sim, sim_peak, YMAX1, C_SIM)
ax.set_xlim(0, 26); ax.set_ylim(-2, YMAX1)
ax.set_xticks(range(1, 26))
ax.set_xlabel("Run number")
ax.set_ylabel("Peak relative OOP displacement (mm)")
ax.legend(loc="upper left", fontsize=13)
fig.tight_layout()
fig.savefig(str(OUT / "figP1_peak_disp.png"), dpi=300, bbox_inches="tight")
plt.close(fig); print("-> figP1_peak_disp.png")

# ------------------------------------------------- figP2: Sd + period
fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(11, 9.5),
                               gridspec_kw={"height_ratios": [1, 1.15]})
rr   = sorted(drv)
sd   = np.array([float(drv[r]["Sd_target_mm"]) for r in rr])
amp  = np.array([float(drv[r]["amplification"]) for r in rr])
ax2 = ax1.twinx()
ax2.bar(rr, amp, color="lightcoral", alpha=0.55, width=0.62, zorder=1)
ax2.set_ylabel("Amplification factor", color="red")
ax2.tick_params(axis="y", colors="red")
ax2.grid(False)
ax1.plot(rr, sd, "s-", color="black", ms=7, lw=1.6,
         label=r"Adaptive $S_d$ target", zorder=3)
ax1.set_ylabel(r"$S_d$ target (mm)")
ax1.set_xlim(0, 26); ax1.set_xticks(range(1, 26))
ax1.set_zorder(ax2.get_zorder() + 1); ax1.patch.set_visible(False)
h1, l1 = ax1.get_legend_handles_labels()
import matplotlib.patches as mpatches
h1.append(mpatches.Patch(color="lightcoral", alpha=0.55,
                         label="Amplification factor"))
ax1.legend(handles=h1, loc="upper left", fontsize=12)

Tend = np.array([float(drv[r]["T_end"]) for r in rr])
fallback = np.array([float(drv[r]["T_end_over_Tinit"]) == 1.0 and r > 1
                     for r in rr])
ok = ~fallback
ax3.plot(np.array(rr)[ok], Tend[ok], "s", color=C_SIM, ms=8,
         label="Simulation", zorder=3)
if fallback.any():
    ax3.plot(np.array(rr)[fallback], Tend[fallback], "s", mfc="white",
             mec=C_SIM, ms=8, label="Simulation (period-ID fallback)",
             zorder=3)
ax3.plot(np.array(rr)[ok], Tend[ok], "-", color=C_SIM, lw=1.0, alpha=0.6,
         zorder=2)
ax3.axhline(T_INIT, color="gray", ls="--", lw=1.0)
ax3.text(25.6, T_INIT, r"$T_{init}$ = %.3f s" % T_INIT, fontsize=11,
         color="gray", ha="right", va="bottom")
ax3.set_xlim(0, 26); ax3.set_xticks(range(1, 26))
ax3.set_xlabel("Run number")
ax3.set_ylabel("Fundamental period (s)")
ax3.legend(loc="upper left", fontsize=12)
fig.tight_layout()
fig.savefig(str(OUT / "figP2_sd_period.png"), dpi=300, bbox_inches="tight")
plt.close(fig); print("-> figP2_sd_period.png")

# ------------------------------------------------- figP3: residuals
YMAX3 = 25.0
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.6))
exp_res  = np.array([float(exp[r]["resid_rel_mm"]) for r in runs_exp])
b0 = float(exp[runs_exp[0]]["raw_base4_mm"])
exp_cum  = np.array([float(exp[r]["raw_base4_mm"]) - b0 for r in runs_exp])
sr = sorted(sim_res_disp)
sim_res  = np.array([sim_res_disp[r] for r in sr])
axL.plot(runs_exp, exp_cum, "o--", color=C_US1, ms=6, lw=1.3,
         mec="darkblue", label="US-1: cumulative lean (baseline drift)")
axL.plot(runs_exp, exp_res, "^--", color="darkturquoise", ms=5, lw=1.0,
         alpha=0.8, label="US-1: per-run residual")
mm = sim_res <= YMAX3
axL.plot(np.array(sr)[mm], sim_res[mm], "s-", color=C_SIM, ms=6, lw=1.6,
         label="Simulation: cumulative residual", zorder=3)
annotate_offscale(axL, sr, sim_res, YMAX3, C_SIM)
axL.set_xlim(0, 26); axL.set_ylim(-3, YMAX3)
axL.set_xticks(range(1, 26, 2))
axL.set_xlabel("Run number")
axL.set_ylabel("Residual rel. OOP displacement (mm)")
axL.legend(loc="upper left", fontsize=10)

exp_rt = np.array([float(exp[r]["resid_tilt_deg"]) for r in runs_exp])
st = sorted(sim_res_tilt)
sim_rt = np.array([sim_res_tilt[r] for r in st])
YMAXT = 0.8
axR.plot(runs_exp, exp_rt, "^--", color="darkturquoise", ms=5, lw=1.0,
         alpha=0.8, label="US-1: per-run residual tilt")
mt = sim_rt <= YMAXT
axR.plot(np.array(st)[mt], sim_rt[mt], "s-", color=C_SIM, ms=6, lw=1.6,
         label="Simulation: cumulative residual tilt", zorder=3)
for xi, yi in zip(st, sim_rt):
    if yi > YMAXT:
        axR.annotate("%.1f deg" % yi, xy=(xi, YMAXT * 0.97),
                     xytext=(xi - 1.5, YMAXT * 0.80), ha="center",
                     fontsize=11, color=C_SIM,
                     arrowprops=dict(arrowstyle="->", color=C_SIM, lw=1.2))
axR.set_xlim(0, 26); axR.set_ylim(-0.1, YMAXT)
axR.set_xticks(range(1, 26, 2))
axR.set_xlabel("Run number")
axR.set_ylabel("Residual tilt (deg)")
axR.legend(loc="upper left", fontsize=10)
fig.tight_layout()
fig.savefig(str(OUT / "figP3_residual.png"), dpi=300, bbox_inches="tight")
plt.close(fig); print("-> figP3_residual.png")

print("Done. Output:", OUT)
