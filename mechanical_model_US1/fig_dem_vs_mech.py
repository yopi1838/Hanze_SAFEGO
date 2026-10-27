"""DEM pushover (v5, f_t = 0, elastic top) against the rigid two-block model evaluated with the mechanism the DEM
actually forms (deformed_shape_3dec.py, 23 Oct): middle hinge at the upper-floor course, top pivot at the top-beam
centroid (the glued top joint rotates with the upper block), lower block kinked at the lower floor.

DEM force, two measures:
  F_tot as logged   base + top-joint SHEAR z-components. Once the top beam tilts by theta the top-joint normal force N
                    has a horizontal component N sin(theta) that this sum leaves out (pos 1.83 kN, neg 2.61 kN at 80 mm,
                    equal to the logged F_tot - |a| m_tot gap within 0.02 kN).
  |a| x m_tot       the applied push on the masonry at a settled servo point = the horizontal force the wall carries.
                    This is the quantity the rigid model's F = a M_wall is.
Rigid model: rigid_two_block_vw.curve with the DEM's loads (masonry 1308.3 kg, O = 10.35 spring + 1.485 top beam at
the top pivot, joists pinned with their bearing load), band = lower pivot at the base ... at 0.30 m."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rigid_two_block_vw as vw

HERE = os.path.dirname(os.path.abspath(__file__))
COMBINED = "--combined" in sys.argv    # panel letters (e)/(f) and no suptitle, for fig_mech_model_and_pushover.png
TB_Y = 2.629      # top-beam centroid height in the build (deformed_shape ref dump), mid-thickness
INK, MUTED, GRID, DEMC, BAND = "#1f2937", "#6b7280", "#e5e7eb", "#E07A2D", "#2B5F8A"


def read_po(path):
    m_tot = float(open(path).readline().split("m_tot")[1].split("kg")[0])
    d = np.genfromtxt(path, delimiter=",", names=True, skip_header=1, dtype=None, encoding=None)
    return np.abs(d["d_ctrl_mm"]), np.abs(d["F_tot_kN"]), np.abs(d["a_mps2"]) * m_tot / 1e3, m_tot


fig, axs = plt.subplots(1, 2, figsize=(13.5, 6.2), sharey=True)
rows = []
for ax, lab, fn, side, h1 in ((axs[0], "Positive push (towards the joists)", "pushover_v5_ft0_elastictop_pos.csv", "leeward", 2.09),
                              (axs[1], "Negative push", "pushover_v5_ft0_elastictop_neg.csv", "windward", 2.20)):
    d, Ftot, Fam, m_tot = read_po(os.path.join(HERE, fn))
    W, O = m_tot * 9.81 / 1e3, 10.35 + 212.986 * 9.81 / 1e3   # DEM: springs + top beam (213 kg in the build)
    J = vw.us1_joists(side, k_r=0.0)
    cur = {}
    for key, kw in (("tb0", dict(top=(vw.t / 2, TB_Y), base_y=0.0)), ("tb3", dict(top=(vw.t / 2, TB_Y), base_y=0.30)),
                    ("formula", dict())):
        hh = 1.51 if key == "formula" else h1
        dd, FF = vw.curve(hh, W, O, joists=J, thmax=0.2, n=500, **kw)
        k = np.argmax(FF < 0) if (FF < 0).any() else len(FF)
        cur[key] = (np.r_[dd[:k], np.interp(0, -FF[k - 1:k + 1], dd[k - 1:k + 1]) if k < len(FF) else dd[-1]],
                    np.r_[FF[:k], 0.0 if k < len(FF) else FF[-1]])
    xs = np.linspace(0.5, 150, 300)
    lo = np.interp(xs, *cur["tb0"], right=np.nan); hi = np.interp(xs, *cur["tb3"], right=np.nan)
    ax.fill_between(xs, lo, hi, color=BAND, alpha=0.18, lw=0)
    ax.plot(*cur["tb0"], color=BAND, lw=2.0, label="rigid, adjusted hinge (2.09 / 2.20 m)")
    ax.plot(*cur["tb3"], color=BAND, lw=1.2, ls="--", label="rigid, adjusted hinge, lower pivot at 0.30 m")
    ax.plot(*cur["formula"], color=INK, lw=2.2, ls=(0, (6, 3)), label="rigid, mid-height hinge (1.51 m)")
    zf, zt = cur["formula"][0][-1], cur["tb0"][0][-1]
    ax.annotate("mid-height hinge:\nzero at {:.0f} mm".format(zf), xy=(zf, 0), xytext=(zf - 52, 1.0), fontsize=8.5, color=INK,
                fontweight="bold", arrowprops=dict(arrowstyle="-", color=INK, lw=0.9))
    ax.annotate("adjusted hinge:\nzero at {:.0f} mm".format(zt), xy=(zt, 0), xytext=(zt + 4, 1.0), fontsize=8.5, color=BAND,
                ha="left", arrowprops=dict(arrowstyle="-", color=BAND, lw=0.9))
    ax.plot(d, Fam, color=DEMC, lw=2.4, label="DEM, |a|·m_tot")
    ax.plot(d, Ftot, color=DEMC, lw=1.2, ls=":", label="DEM, F_tot as logged")
    for x in (20, 40, 60, 80):
        rows.append(dict(direction=lab.split()[0].lower(), d_mm=x, DEM_applied=np.interp(x, d, Fam), DEM_Ftot=np.interp(x, d, Ftot),
                         rigid_TB_base=np.interp(x, *cur["tb0"]), rigid_TB_030=np.interp(x, *cur["tb3"]),
                         rigid_formula=np.interp(x, *cur["formula"])))
    ax.annotate("gap at 80 mm:\n{:.1f} kN ({:.1f}-{:.1f} vs {:.1f})".format(np.interp(80, d, Fam) - np.interp(80, *cur["tb3"]),
                np.interp(80, *cur["tb0"]), np.interp(80, *cur["tb3"]), np.interp(80, d, Fam)),
                xy=(80, np.interp(80, d, Fam)), xytext=(96, 8.6), fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_title(("(e) " if ax is axs[0] else "(f) ") + lab if COMBINED else lab, loc="left", fontsize=11.5 if COMBINED else 11,
                 fontweight="bold")
    ax.set_xlim(0, 170); ax.set_ylim(0, 14); ax.grid(color=GRID); ax.set_axisbelow(True)
    ax.set_xlabel("displacement at the 2.06 m sensor (mm)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
axs[0].set_ylabel("lateral force on the wall (kN)")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, frameon=False, fontsize=8.6, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
if not COMBINED:
  fig.suptitle("DEM pushover (f_t = 0, elastic top) vs the rigid two-block model",
             x=0.01, ha="left", fontsize=12)
fig.tight_layout(rect=(0, 0.11, 1, 1))
fig.savefig(os.path.join(HERE, "_pushover_panel.png" if COMBINED else "fig_dem_vs_mech.png"), dpi=float(os.environ.get("PO_DPI", 150)))
import csv
with open(os.path.join(HERE, "dem_vs_mech.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
    w.writerows([{k: (round(float(v), 2) if not isinstance(v, str) else v) for k, v in r.items()} for r in rows])
for r in rows:
    print(" ".join("{}={}".format(k, round(float(v), 2) if not isinstance(v, str) else v) for k, v in r.items()))
