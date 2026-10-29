# -*- coding: ascii -*-
"""fig_force_decomposition.py -- lateral force of the slip-work pushover, split into rigid rocking, measured slip
(equivalent force dW/dZ, Z = the rigid model's mass-weighted sway) and the remainder.
Usage: python fig_force_decomposition.py [DIR]   (default pushover_results_v5_ft_0p050_TopMason_slipwork)"""
import os, sys, io
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
DIR = sys.argv[1] if len(sys.argv) > 1 else "pushover_results_v5_ft_0p050_TopMason_slipwork"
sys.path.insert(0, os.path.abspath("mechanical_model_US1"))
import rigid_two_block_vw as vw
W = 1308.3 * 9.81 / 1e3; O = 10.35 + 212.986 * 9.81 / 1e3; H1, TB = 2.09, 2.629; MU = 0.75; PIVOT = 0.30
REG = ["pocket_low", "pocket_up", "base", "top", "mas_lowfloor", "mas_upfloor", "mas_other", "new"]   # supports excluded
INK, GRID, DEM, RIG, MOD, SLIP, REST = "#1f2937", "#e5e7eb", "#16A34A", "#0F766E", "#134E4A", "#CA8A04", "#9CA3AF"
fig, axs = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
rows = []
for ax, (br, side, title) in zip(axs, (("pos", "leeward", "Positive push (towards the joists)"), ("neg", "windward", "Negative push"))):
    L = [l for l in open(os.path.join(DIR, "pushover_{}.csv".format(br))) if l.strip() and not l.startswith("#")]
    p = pd.read_csv(io.StringIO("".join(L))); p = p[p.converged == 1]
    dD, FD = p.d_ctrl_mm.abs().values, p.a_mps2.abs().values * 1.3083
    th = np.linspace(1e-5, 0.12, 600)
    st = np.array([vw.state(a, H1, W, O, joists=vw.us1_joists(side, k_r=0.0), top=(vw.t / 2, TB), base_y=PIVOT) for a in th])
    U, Z, dr = st.T; dr = dr * 1e3; Fr = np.gradient(U, th) / np.gradient(Z, th)
    dm, Fm = vw.curve(H1, W, O, joists=vw.us1_joists(side, k_r=0.0, mu=MU, W=W, O=O), thmax=0.12, n=600,
                      top=(vw.t / 2, TB), base_y=PIVOT)
    s = pd.read_csv(os.path.join(DIR, "slipwork_{}.csv".format(br)), comment="#")
    s["dd"] = s.d_ctrl_mm.abs(); s = s.sort_values("dd").drop_duplicates("dd")
    Wt = sum(s[r + "_Wshear_J"] - (s[r + "_Eel_J"] - s[r + "_Eel_J"].iloc[0]) for r in REG).values
    xs = np.arange(10.0, min(dD.max(), s.dd.max()) - 2.5, 1.0)             # equivalent slip force, 5 mm window; from 10 mm (rocking)
    Zx = lambda x: np.interp(x, dr, Z)
    Fs = np.array([(np.interp(x + 2.5, s.dd, Wt) - np.interp(x - 2.5, s.dd, Wt)) / (Zx(x + 2.5) - Zx(x - 2.5)) / 1e3 for x in xs])
    Frx = np.interp(xs, dr, Fr); FDx = np.interp(xs, dD, FD)
    ax.fill_between(xs, Frx, Frx + Fs, color=SLIP, alpha=0.35, lw=0, label="measured slip work, as force dW/dZ")
    ax.fill_between(xs, Frx + Fs, np.maximum(FDx, Frx + Fs), color=REST, alpha=0.30, lw=0, label="remainder: not dissipated")
    ax.plot(dr, Fr, color=RIG, lw=1.6, label="rigid model, lower pivot {:.2f} m (no friction)".format(PIVOT))
    ax.plot(dm, Fm, color=MOD, lw=1.4, ls="--", label=r"rigid model + pocket friction ($\mu$ 0.75, N {:.1f} kN)".format(vw.pocket_clamp_N(W, O, 0.38, side)))
    ax.plot(dD, FD, color=DEM, lw=2.4, label="DEM, f_t = 0.05 MPa (slip-work run)")
    for x in (20, 40, 60, 80):
        if x <= xs.max() + 2.5:
            i = int(np.argmin(np.abs(xs - min(x, xs.max()))))
            rows.append(dict(branch=br, d_mm=x, DEM_kN=np.interp(x, dD, FD), rigid_kN=np.interp(x, dr, Fr),
                             slip_eq_kN=Fs[i], remainder_kN=np.interp(x, dD, FD) - np.interp(x, dr, Fr) - Fs[i],
                             model_friction_kN=np.interp(x, dm, Fm)))
    ax.set_title(title, loc="left", fontweight="bold", color=INK)
    ax.set_xlim(0, 85); ax.set_ylim(0, 14); ax.set_xlabel("displacement at the 2.06 m sensor (mm)")
    ax.grid(color=GRID); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
axs[0].set_ylabel("lateral force on the wall (kN)")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
fig.suptitle("Lateral force: rigid rocking + measured slip + remainder (DEM pushover with slip tracking)", x=0.01, ha="left", color=INK)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(os.path.join(DIR, "fig_force_decomposition.png"), dpi=150)
pd.DataFrame(rows).round(2).to_csv(os.path.join(DIR, "force_decomposition.csv"), index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False))
