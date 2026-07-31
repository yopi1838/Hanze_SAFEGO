# -*- coding: ascii -*-
"""
Presentation figures: experiment (US-1 / Test 9, US-2 / Test 12) vs
simulation, styled after the IMC2026 paper figures.

Paths resolve through safego_paths.py, so this follows SAFEGO_SIM_DIR rather
than being pinned to one results folder:

    SAFEGO_SIM_DIR=stratC_results_NODAMP_v7 python exp_vs_sim_figs.py

Inputs:
  <SIM>/postproc/postproc_summary.csv       (from postprocess_stratC.py)
  <SIM>/postproc/postproc_all_channels.csv
  <SIM>/postproc/exp_Test9_metrics.csv      (experimental metrics, Test 9)

Outputs (into <SIM>/postproc/):
  figP1_peak_disp_exp_vs_sim.png   Fig.5 style: peak rel OOP disp per run
  figP2_sd_target_period.png       Fig.4 style: Sd target + amplification,
                                   period evolution (sim; exp period TBD)
  figP3_residual_tilt.png          residual disp + peak tilt, exp vs sim

The sim-derived inputs come from the active simulation's postproc/. The
derived experimental CSVs (exp_*) describe the specimens rather than any one
run, so they fall back to the canonical copies under
stratC_results_NODAMP_v6_NEW/postproc/ when the active sim has none -- see
safego_paths.exp_derived(). US-2 series are drawn when their CSV is present.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import safego_paths as sp

SIM = sp.sim_dir()
PP = sp.postproc_dir()
print("sim dir : {}".format(SIM))
print("postproc: {}".format(PP))

# ---- paper-like style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman"],
    "mathtext.fontset": "dejavuserif",
    "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.grid": True, "grid.linestyle": ":", "grid.alpha": 0.6,
    "legend.frameon": True, "legend.framealpha": 1.0,
    "legend.edgecolor": "0.7",
})
C_EXP1 = "royalblue"
C_EXP2 = "orange"
C_SIM = "red"

def fread(path):
    with open(str(path)) as f:
        return list(csv.DictReader(f))

def col(rows, k):
    return np.array([float(r[k]) for r in rows])

sim = fread(PP / "postproc_summary.csv")
exp = fread(sp.exp_derived("exp_Test9_metrics.csv"))
runs_s = col(sim, "run"); runs_e = col(exp, "run")
t12f = sp.exp_derived("exp_Test12_tilt.csv")
exp2 = fread(t12f) if t12f.exists() else None
if exp2 is not None:
    runs_e2 = col(exp2, "run")
    exp2_tilt_full = col(exp2, "peak_tilt_full_wall")
    # top rel disp recovered from the full-wall chord tilt (dy = 2.06 m)
    exp2_peak_mm = np.tan(np.radians(exp2_tilt_full)) * 2060.0

# ---------------------------------------------------------------- fig P1
sim_peak = col(sim, "peak_rel_disp_mm")
exp_peak = col(exp, "peak_rel_mm")
YCAP = 35.0

fig, ax = plt.subplots(figsize=(9.2, 5.0))
ax.plot(runs_e, exp_peak, "o--", color=C_EXP1, ms=6, lw=1.2,
        label="US-1 (Test 9)", zorder=3)
if exp2 is not None:
    ax.plot(runs_e2, exp2_peak_mm, "o--", color=C_EXP2, ms=6, lw=1.2,
            label="US-2 (Test 12)", zorder=3)
clipped = np.minimum(sim_peak, YCAP)
ax.plot(runs_s, clipped, "s-", color=C_SIM, ms=6, lw=1.6,
        label="NUM", zorder=4)
off = runs_s[sim_peak > YCAP]
for k, rn in enumerate(off):
    v = sim_peak[runs_s == rn][0]
    ax.plot(rn, YCAP, "s", ms=8, mfc="white", mec=C_SIM, mew=1.6, zorder=5)
    ax.annotate("Run %d: %.0f mm" % (rn, v), xy=(rn, YCAP),
                xytext=(rn - 9.5, YCAP - 3.5 - 3.2 * k), fontsize=10,
                color=C_SIM, ha="left",
                arrowprops=dict(arrowstyle="->", color=C_SIM, lw=1.0))
ax.set_xlim(0, 26); ax.set_ylim(-2, YCAP)
ax.set_xticks(range(0, 27, 2))
ax.set_xlabel("Run number", fontsize=14)
ax.set_ylabel("Peak relative OOP displacement (mm)", fontsize=14)
ax.legend(fontsize=11, loc="upper left")
fig.tight_layout()
fig.savefig(str(PP / "figP1_peak_disp_exp_vs_sim.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig); print("-> figP1_peak_disp_exp_vs_sim.png")

# ---------------------------------------------------------------- fig P2
sd = col(sim, "Sd_target_mm")
amp = col(sim, "amplification")
t_end = col(sim, "T_end")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.2, 8.2), sharex=True)
axb = ax1.twinx()
axb.bar(runs_s, amp, color="lightcoral", alpha=0.45, width=0.62, zorder=1)
axb.set_ylabel("Amplification factor", color="red", fontsize=14)
axb.tick_params(axis="y", colors="red")
axb.set_ylim(0, max(3.0, amp.max() * 1.05))
axb.axhline(1.0, color="red", ls=":", lw=0.8, alpha=0.5)
axb.grid(False)
ax1.plot(runs_s, sd, "s-", color="black", ms=6, lw=1.5,
         label="Adaptive $S_d$ target", zorder=3)
ax1.bar([np.nan], [np.nan], color="lightcoral", alpha=0.45,
        label="Amplification factor")
ax1.set_ylabel("$S_d$ target (mm)", fontsize=14)
ax1.legend(fontsize=11, loc="upper left")
ax1.set_zorder(axb.get_zorder() + 1)
ax1.patch.set_visible(False)

ratio = col(sim, "T_end_over_Tinit")
fallback = ratio == 1.0  # FFT period-ID fallback, not true recovery
ax2.plot(runs_s[~fallback], t_end[~fallback], "s", color=C_SIM, ms=6,
         label="NUM")
ax2.plot(runs_s[fallback], t_end[fallback], "s", ms=7, mfc="white",
         mec=C_SIM, mew=1.4, label="Period-ID fallback (excluded)")
ax2.plot(runs_s[~fallback], t_end[~fallback], "-", color=C_SIM,
         lw=1.0, alpha=0.6)
ax2.set_xlabel("Run number", fontsize=14)
ax2.set_ylabel("Fundamental period (s)", fontsize=14)
ax2.legend(fontsize=11, loc="upper left")
ax2.set_xlim(0, 26); ax2.set_xticks(range(0, 27, 2))
fig.tight_layout()
fig.savefig(str(PP / "figP2_sd_target_period.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig); print("-> figP2_sd_target_period.png")

# ---------------------------------------------------------------- fig P3
sim_res = col(sim, "residual_tilt_deg")
sim_peak_tilt = col(sim, "peak_tilt_deg")
# measured tilt from the raw Test9 dataset (exp_tilt_from_raw.py); falls
# back to the ch4-derived chord tilt if the raw extraction is absent
tiltf = sp.exp_derived("exp_Test9_tilt.csv")
if tiltf.exists():
    expt = fread(tiltf)
    runs_t = col(expt, "run")
    exp_peak_tilt = col(expt, "peak_tilt_full_wall")
else:
    expt = None
    runs_t = runs_e
    exp_peak_tilt = col(exp, "peak_tilt_deg")
exp_res_disp = col(exp, "resid_rel_mm")
exp_cum = col(exp, "raw_base4_mm") - col(exp, "raw_base4_mm")[0]

# sim residual displacement from the all-channels table
allch = fread(PP / "postproc_all_channels.csv")
rd = {int(float(r["run"])): float(r["residual"]) for r in allch
      if r["channel"] == "rel_disp_top_mm"}
sim_res_disp = np.array([rd[int(r)] for r in runs_s])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.0))
ax1.plot(runs_e, exp_cum, "o--", color=C_EXP1, ms=6, lw=1.2,
         label="US-1: cumulative lean")
ax1.plot(runs_e, exp_res_disp, "^--", color="darkturquoise", ms=5, lw=1.0,
         alpha=0.85, label="US-1: per-run residual")
CAP1 = 40.0
clip1 = np.minimum(sim_res_disp, CAP1)
ax1.plot(runs_s, clip1, "s-", color=C_SIM, ms=6, lw=1.6,
         label="NUM: cumulative residual")
for rn in runs_s[sim_res_disp > CAP1]:
    v = sim_res_disp[runs_s == rn][0]
    ax1.plot(rn, CAP1, "s", ms=8, mfc="white", mec=C_SIM, mew=1.6, zorder=5)
    ax1.annotate("Run %d: %.0f mm" % (rn, v), xy=(rn, CAP1),
                 xytext=(rn - 10.5, CAP1 - 6), fontsize=10, color=C_SIM,
                 arrowprops=dict(arrowstyle="->", color=C_SIM, lw=1.0))
ax1.set_xlim(0, 26); ax1.set_ylim(-6, 42)
ax1.set_xlabel("Run number", fontsize=14)
ax1.set_ylabel("Residual OOP displacement (mm)", fontsize=14)
ax1.legend(fontsize=10, loc="upper left")

ax2.plot(runs_t, exp_peak_tilt, "o--", color=C_EXP1, ms=6, lw=1.2,
         label="US-1 (Test 9)")
if exp2 is not None:
    ax2.plot(runs_e2, exp2_tilt_full, "o--", color=C_EXP2, ms=6, lw=1.2,
             label="US-2 (Test 12)")
CAP2 = 1.5
clip2 = np.minimum(sim_peak_tilt, CAP2)
ax2.plot(runs_s, clip2, "s-", color=C_SIM, ms=6, lw=1.6,
         label="NUM")
for rn in runs_s[sim_peak_tilt > CAP2]:
    v = sim_peak_tilt[runs_s == rn][0]
    ax2.plot(rn, CAP2, "s", ms=8, mfc="white", mec=C_SIM, mew=1.6, zorder=5)
    ax2.annotate("Run %d: %.1f deg" % (rn, v), xy=(rn, CAP2),
                 xytext=(rn - 10.5, CAP2 - 0.22), fontsize=10, color=C_SIM,
                 arrowprops=dict(arrowstyle="->", color=C_SIM, lw=1.0))
ax2.set_xlim(0, 26); ax2.set_ylim(-0.05, 1.58)
ax2.set_xlabel("Run number", fontsize=14)
ax2.set_ylabel("Peak wall tilt (deg)", fontsize=14)
ax2.legend(fontsize=10, loc="upper left")
fig.tight_layout()
fig.savefig(str(PP / "figP3_residual_tilt.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig); print("-> figP3_residual_tilt.png")
print("Done.")


# ---------------------------------------------------------------- fig P4
# Segment-by-segment tilt comparison (mechanism check). Sim channels from
# postproc_all_channels.csv; experiment from exp_Test9_tilt.csv.
if expt is not None:
    seg_names = ["tilt_bot_seg", "tilt_low_seg", "tilt_up_seg",
                 "tilt_full_wall"]
    seg_lbl = {"tilt_bot_seg": "Bottom segment (0 - 0.66 m)",
               "tilt_low_seg": "Lower segment (0.66 - 1.26 m)",
               "tilt_up_seg": "Upper segment (1.26 - 2.06 m)",
               "tilt_full_wall": "Full wall (0 - 2.06 m)"}
    sim_seg = {sn: {} for sn in seg_names}
    for r in allch:
        if r["channel"] in sim_seg:
            sim_seg[r["channel"]][int(float(r["run"]))] = float(r["peak"])
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), sharex=True)
    CAPS = 1.6
    for ax, sn in zip(axes.ravel(), seg_names):
        ye = col(expt, "peak_" + sn)
        ax.plot(runs_t, ye, "o--", color=C_EXP1, ms=5, lw=1.1,
                label="US-1 (Test 9)")
        if exp2 is not None:
            ax.plot(runs_e2, col(exp2, "peak_" + sn), "o--", color=C_EXP2,
                    ms=5, lw=1.1, label="US-2 (Test 12)")
        rs = sorted(sim_seg[sn])
        ys = np.array([sim_seg[sn][r] for r in rs])
        ax.plot(rs, np.minimum(ys, CAPS), "s-", color=C_SIM, ms=5, lw=1.4,
                label="NUM")
        for rr in [r for r, v in zip(rs, ys) if v > CAPS]:
            v = ys[rs.index(rr)]
            ax.plot(rr, CAPS, "s", ms=7, mfc="white", mec=C_SIM, mew=1.4)
            ax.annotate("%.1f" % v, xy=(rr, CAPS), xytext=(rr - 3.2, CAPS - 0.18),
                        fontsize=9, color=C_SIM)
        ax.set_title(seg_lbl[sn], fontsize=12)
        ax.set_ylim(-0.05, CAPS + 0.1)
        ax.set_xlim(0, 26)
        ax.legend(fontsize=9, loc="upper left")
    for ax in axes[1]:
        ax.set_xlabel("Run number", fontsize=13)
    for ax in axes[:, 0]:
        ax.set_ylabel("Peak segment tilt (deg)", fontsize=13)
    fig.tight_layout()
    fig.savefig(str(PP / "figP4_tilt_segments_exp_vs_sim.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("-> figP4_tilt_segments_exp_vs_sim.png")


# ---------------------------------------------------------------- fig P5
# Residual tilt angle vs run: US-1 tiltmeter (sensor 18, settled values,
# OOP = sensor Y axis) against the simulation cumulative residual tilt.
tmf = sp.exp_derived("exp_US1_tiltmeter.csv")
if tmf.exists():
    tm = fread(tmf)
    runs_m = col(tm, "run")
    exp_res_tilt = col(tm, "residual_deg")
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.plot(runs_m, exp_res_tilt, "o--", color=C_EXP1, ms=6, lw=1.2,
            label="US-1 tiltmeter (settled, cumulative)")
    tmf2 = sp.exp_derived("exp_US2_tiltmeter.csv")
    if tmf2.exists():
        tm2 = fread(tmf2)
        ax.plot(col(tm2, "run"), col(tm2, "residual_deg"), "o--",
                color=C_EXP2, ms=6, lw=1.2,
                label="US-2 tiltmeter (settled, cumulative)")
    ax.plot(runs_s, sim_res, "s-", color=C_SIM, ms=6, lw=1.6,
            label="NUM (full-wall residual, cumulative)")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlim(0, 26); ax.set_xticks(range(0, 27, 2))
    ax.set_xlabel("Run number", fontsize=14)
    ax.set_ylabel("Residual tilt angle (deg)", fontsize=14)
    ax.legend(fontsize=11, loc="upper left")
    fig.tight_layout()
    fig.savefig(str(PP / "figP5_tiltmeter_residual.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("-> figP5_tiltmeter_residual.png")
