import numpy as np, csv, math, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

SIM = Path("/sessions/zen-brave-heisenberg/mnt/Hanze_SAFEGO/stratC_results_NODAMP_v6_NEW")
OUT = SIM / "postproc"
REC = {1:"HU12",2:"HU12",3:"EC40",4:"HU12",5:"HU12",6:"EC40",7:"HU12",8:"EC40",
       9:"HU12",10:"HU12",11:"EC40",12:"HU12",13:"HU12",14:"HU12",15:"HU12",
       16:"HU12",17:"HU12",18:"HU12",19:"HU12",20:"HU12",21:"HU12",
       22:"FR76",23:"FR76",24:"FR76",25:"FR76"}

# --- experiment ---
exp = {int(r["run"]): r for r in csv.DictReader(open("/tmp/exp_metrics.csv"))}
xr = sorted(exp)
exp_peak  = np.array([float(exp[r]["peak_rel_mm"]) for r in xr])
exp_resid = np.array([float(exp[r]["resid_rel_mm"]) for r in xr])
# cumulative lean from raw sensor baseline drift (absolute reading)
b0 = float(exp[1]["raw_base4_mm"])
exp_cum = np.array([float(exp[r]["raw_base4_mm"]) - b0 for r in xr])

# --- simulation: peak and residual rel disp per run ---
sim_peak, sim_resid, run1_start = {}, {}, None
for rn in range(1, 26):
    fl = glob.glob(str(SIM / ("Run%02d_*" % rn) / "*rel_disp_top_mm*"))
    if not fl: continue
    d = np.genfromtxt(fl[0], skip_header=2)
    v = d[np.isfinite(d[:,1]), 1]
    if run1_start is None: run1_start = float(v[0])
    sim_peak[rn]  = float(np.max(np.abs(v - v[0])))
    sim_resid[rn] = float(np.mean(v[int(0.95*len(v)):])) - run1_start
sr = sorted(sim_peak)

# --- comparison csv ---
with open(str(OUT / "xval_exp_vs_sim.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run","record","exp_peak_mm","sim_peak_mm","sim_over_exp",
                "exp_resid_perrun_mm","exp_cum_lean_mm","sim_resid_cum_mm"])
    for rn in range(1, 26):
        e = exp.get(rn); i = rn - 1
        ep = float(e["peak_rel_mm"]) if e else None
        sp = sim_peak.get(rn)
        w.writerow([rn, REC[rn],
                    round(ep,2) if ep is not None else "",
                    round(sp,2) if sp is not None else "",
                    round(sp/ep,2) if (ep and sp and ep>0.05) else "",
                    round(float(e["resid_rel_mm"]),2) if e else "",
                    round(exp_cum[i],2) if e else "",
                    round(sim_resid[rn],2) if rn in sim_resid else ""])
print("wrote xval_exp_vs_sim.csv")

# --- figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
ax1.plot(xr, exp_peak, "bo-", ms=5, lw=0.9, label="Experiment Test9 (top-quarter rel. disp)")
ax1.plot(sr, [sim_peak[r] for r in sr], "rs-", ms=5, lw=0.9, label="Sim NODAMP v6 NEW (updated Maxwell)")
ax1.set_yscale("log")
ax1.set_xlabel("Run"); ax1.set_ylabel("Peak |rel. top OOP disp| (mm, log)")
ax1.set_title("Peak relative displacement")
ax1.grid(True, alpha=0.15, which="both"); ax1.legend(fontsize=9)
ax2.plot(xr, exp_cum, "bo-", ms=5, lw=0.9, label="Experiment: cumulative lean (sensor baseline drift)")
ax2.plot(xr, exp_resid, "c^--", ms=4, lw=0.8, alpha=0.7, label="Experiment: per-run residual")
ax2.plot(sr, [sim_resid[r] for r in sr], "rs-", ms=5, lw=0.9, label="Sim: cumulative residual")
ax2.set_ylim(-5, 30)
ax2.set_xlabel("Run"); ax2.set_ylabel("Residual rel. top OOP disp (mm)")
ax2.set_title("Residual (permanent) displacement -- sim Run 25 at %.0f mm (off scale)" % sim_resid[25])
ax2.grid(True, alpha=0.15); ax2.legend(fontsize=8)
fig.suptitle("Cross-validation vs Test9: sim over-ratchets -- onset Run 14 vs exp Run ~23, "
             "peaks ~2-3x exp from Run 15", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(str(OUT / "fig4_xval_exp_vs_sim.png"), dpi=300, bbox_inches="tight")
print("wrote fig4_xval_exp_vs_sim.png")

# console summary
print("\n run  rec   exp_peak  sim_peak  ratio   exp_cum  sim_resid")
for rn in range(1, 26):
    e = exp.get(rn)
    ep = float(e["peak_rel_mm"]) if e else float("nan")
    sp = sim_peak.get(rn, float("nan"))
    print("%4d %5s %9.2f %9.2f %6.2f %9.2f %9.2f" % (
        rn, REC[rn], ep, sp, sp/ep if ep>0.05 else float("nan"),
        exp_cum[rn-1] if e else float("nan"), sim_resid.get(rn, float("nan"))))
