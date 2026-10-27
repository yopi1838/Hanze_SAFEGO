"""Illustrations for MECHANICAL_MODEL_US1.md:
   fig_mech_schematic.png  (a) three-hinge kinematics and loads of the US-1 strip, (b) Onur's cavity split vs the solid strip
   fig_mech_curves.png     (c) how the force-displacement law is built, (d) effect of the hinge height
Run from mechanical_model_US1/ (needs rocking_model.py and rigid_two_block_vw.py next to it)."""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Rectangle, Circle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rocking_model import graziotti, doherty_trilinear
try:
    from rigid_two_block_vw import curve, config
except ImportError:
    from vw import curve, config

INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
BLUE, ORANGE, GREEN, PURPLE, RED = "#2a78d6", "#e07a2d", "#1baf7a", "#6250d6", "#c2410c"
t, h, zs, b = 0.21, 2.58, 2.06, 1.292
W, O = 13.15, 13.19
XS = 2.2   # horizontal exaggeration of the drawings (thickness and sway), stated on the figure


def arrow(ax, p0, p1, col=INK, lw=1.4, ms=10, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, **kw))


def X(z):
    return np.asarray(z) * XS


# ================= figure 1: schematic =================
fig, ax = plt.subplots(1, 2, figsize=(14, 8.2), gridspec_kw=dict(width_ratios=[1.15, 1]))

# ---- (a) kinematics of the three-hinge strip ----
a = ax[0]; h1 = 1.51; th1 = 0.055
lo, up, th2 = config(th1, h1)
# undeformed outline
a.add_patch(Rectangle((0, 0), X(t), h, fill=False, ls=(0, (4, 3)), lw=1, ec=MUTED))
# deformed blocks
corners_lo = [(0, 0), (t, 0), (t, h1), (0, h1)]
corners_up = [(0, h1), (t, h1), (t, h), (0, h)]
P_lo = np.array([lo(np.array(c)) for c in corners_lo]); P_up = np.array([up(np.array(c)) for c in corners_up])
for P, col in ((P_lo, "#dbeafe"), (P_up, "#fde7d4")):
    a.add_patch(Polygon(np.c_[X(P[:, 0]), P[:, 1]], closed=True, fc=col, ec=INK, lw=1.4))
# base slab, top beam, springs
a.add_patch(Rectangle((-0.25, -0.12), X(t) + 0.5, 0.12, fc="#d1d5db", ec=INK, lw=1))
top = up(np.array([t, h]))
a.add_patch(Rectangle((X(top[0]) - X(t) - 0.03, top[1]), X(t) + 0.06, 0.08, fc="#9ca3af", ec=INK, lw=1))
for xs_ in (X(top[0]) - X(t) * 0.8, X(top[0]) - X(t) * 0.2):
    yy = np.linspace(top[1] + 0.08, top[1] + 0.36, 9); xx = xs_ + 0.025 * np.where(np.arange(9) % 2, 1, -1); xx[0] = xx[-1] = xs_
    a.plot(xx, yy, color=INK, lw=1)
a.text(X(top[0]) + 0.08, top[1] + 0.26, "springs 10.35 kN\n+ beam (dead load O)", fontsize=8.5, color=INK, va="center")
# lateral restraint of the beam
a.plot([X(top[0]) + 0.02, X(top[0]) + 0.02], [top[1], top[1] + 0.08], color=INK, lw=2)
for k in range(4):
    a.plot([X(top[0]) + 0.02, X(top[0]) + 0.07], [top[1] + 0.02 * k, top[1] + 0.02 * k + 0.02], color=INK, lw=0.8)
a.text(X(top[0]) + 0.09, top[1] - 0.06, "beam held sideways,\nfree to rise", fontsize=8, color=MUTED, va="top")
# hinges
hinges = [(lo(np.array([t, 0.0])), "base pivot (leeward toe)"),
          (lo(np.array([0.0, h1])), "mid hinge at h1 (windward face)"),
          (up(np.array([t, h])), "top pivot (leeward corner)")]
for (p, lab), dx in zip(hinges, (0.12, 0.12, 0.12)):
    a.add_patch(Circle((X(p[0]), p[1]), 0.035, fc="white", ec=RED, lw=2, zorder=5))
    a.annotate(lab, (X(p[0]), p[1]), xytext=(X(p[0]) + (0.25 if dx > 0 else -0.25), p[1] + (0.09 if p[1] != h1 else -0.12)),
               fontsize=8.5, color=RED, ha="left" if dx > 0 else "right", va="center",
               arrowprops=dict(arrowstyle="-", color=RED, lw=0.8))
# distributed lateral load (mass proportional)
for y in np.linspace(0.12, h - 0.12, 12):
    pz = (lo(np.array([0.0, y])) if y <= h1 else up(np.array([0.0, y])))
    arrow(a, (X(pz[0]) - 0.32, y), (X(pz[0]) - 0.02, y), col=BLUE, lw=1.1, ms=8)
a.text(-1.3, 2.45, "F = a·m\nuniform over\nthe wall mass", color=BLUE, fontsize=9, ha="left", va="center")
# self weights
for P, lab in ((P_lo, "W·h1/h"), (P_up, "W·(h−h1)/h")):
    c = P.mean(axis=0); arrow(a, (X(c[0]), c[1] + 0.12), (X(c[0]), c[1] - 0.12), col=INK, lw=1.4)
    a.text(X(c[0]) - 0.04, c[1] + 0.17, lab, fontsize=8.5, color=INK, ha="center")
# sensor
ps = up(np.array([t, zs])) if zs > h1 else lo(np.array([t, zs]))
a.plot([X(t), X(ps[0])], [zs, zs], color=GREEN, lw=2.2)
a.annotate("u_s at 2.06 m\n(sensor, the EDP)", (X(ps[0]), zs), xytext=(X(t) + 0.55, zs - 0.28), fontsize=8.5, color=GREEN,
           arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.8))
# joists
for yj in (0.38, 2.14):
    pj = lo(np.array([t, yj])) if yj <= h1 else up(np.array([t, yj]))
    a.plot([X(pj[0]), X(pj[0]) + 0.75], [yj, yj], color="#a16207", lw=4, solid_capstyle="butt", alpha=0.8)
a.text(X(t) + 0.55, 0.42, "timber joist\n(in pocket)", fontsize=8, color="#a16207")
# dimensions
a.annotate("", (-0.62, 0), (-0.62, h1), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.8))
a.text(-0.66, h1 / 2, "h1", color=MUTED, ha="right", va="center", fontsize=9)
a.annotate("", (-0.62, h1), (-0.62, h), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.8))
a.text(-0.66, (h1 + h) / 2, "h − h1", color=MUTED, ha="right", va="center", fontsize=9)
arrow(a, (X(t) + 0.95, 1.2), (X(t) + 1.35, 1.2), col=INK, lw=1.6)
a.text(X(t) + 1.15, 1.26, "push", ha="center", fontsize=9)
a.set_xlim(-1.35, 2.25); a.set_ylim(-0.25, 3.1); a.set_aspect("equal"); a.axis("off")
a.set_title("(a) US-1 strip: three hinges, mass-proportional push, dead overburden", loc="left", fontsize=11, fontweight="bold")
a.text(-1.3, -0.22, "thickness and sway drawn x{:.1f}; rotation exaggerated".format(XS), fontsize=7.5, color=MUTED)

# ---- (b) Onur's cavity split vs the solid strip ----
c = ax[1]
def leaf(x0, w, y0, y1, fc, ec=INK, ls="-"):
    c.add_patch(Rectangle((x0, y0), w, y1 - y0, fc=fc, ec=ec, lw=1.2, ls=ls))
# cavity wall (left)
xo, xi, wl = 0.0, 0.42, 0.22
leaf(xo, wl, 0, 2.58, "#fef3c7")                    # outer leaf, full height
c.text(xo + wl / 2, 2.68, "outer leaf\nGraziotti strip\n(full height)", ha="center", fontsize=8, va="bottom")
for y0, y1, lab, fc in ((0, 0.454, "Casapulla\ncantilever\n+ joist slip F_Q", "#fde7d4"),
                        (0.454, 2.174, "Graziotti strip\nbetween floors", "#dbeafe"),
                        (2.174, 2.58, "Casapulla\ncantilever", "#fde7d4")):
    leaf(xi, wl, y0, y1, fc)
    c.text(xi + wl + 0.05, (y0 + y1) / 2, lab, fontsize=7.8, va="center")
for yj in (0.454, 2.174):
    c.plot([xi + wl, xi + wl + 0.0], [yj, yj], color="#a16207")
    c.add_patch(Rectangle((xi + wl * 0.4, yj - 0.04), 0.5, 0.08, fc="#a16207", ec="none", alpha=0.7))
for y in np.linspace(0.3, 2.3, 6):
    c.plot([xo + wl, xi], [y, y], color=PURPLE, lw=1, ls=(0, (2, 1.5)))
c.text((xo + wl + xi) / 2, 0.1, "ties\n(F_c)", color=PURPLE, fontsize=7.5, ha="center")
c.text(0.42, -0.35, "Onur: cavity wall C-US-1/2\nF = F_outer + F_inner,strip + F_cantilevers + F_c\n= 8.3 kN (e=0) / 10.4 (e=t/2) / 11.6 (+F_c)",
       ha="center", fontsize=8.3, va="top")
# solid wall (right)
xs0, ws = 2.0, 0.42
leaf(xs0, ws, 0, 2.58, "#dbeafe")
for yy, lab, col in ((1.51, "A1: optimum h1 = 1.51 m", RED), (2.0, "high hinge ~2.0 m\n(what the DEM forms)", ORANGE), (2.14, "A2: top joist 2.14 m", MUTED)):
    c.plot([xs0 - 0.05, xs0 + ws + 0.05], [yy, yy], color=col, lw=1.6, ls="--")
    c.plot([xs0 + ws + 0.05, xs0 + ws + 0.62], [yy, yy], color=col, lw=0.6, ls=":")
    c.text(xs0 + ws + 0.65, yy, lab, color=col, fontsize=7.8, va="center")
for yj in (0.38, 2.14):
    c.add_patch(Rectangle((xs0 + ws * 0.5, yj - 0.04), 0.5, 0.08, fc="#a16207", ec="none", alpha=0.45))
c.text(xs0 + ws / 2, 2.68, "US-1 solid leaf\none strip, full height", ha="center", fontsize=8, va="bottom")
c.text(xs0 + ws / 2 + 0.3, -0.35, "US-1: no ties, no joist slip,\nno split at the floors\nF = F_strip  (9.8 kN e=0 / 12.4 e=t/2 at h1 = 1.51 m)",
       ha="center", fontsize=8.3, va="top")
arrow(c, (1.25, 1.3), (1.85, 1.3), col=MUTED, lw=1.4)
c.text(1.55, 1.38, "solid wall", ha="center", fontsize=8.5, color=MUTED)
c.set_xlim(-0.2, 4.2); c.set_ylim(-1.0, 3.15); c.set_aspect("equal"); c.axis("off")
c.set_title("(b) Onur's cavity split vs the solid-wall strip", loc="left", fontsize=11, fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_mech_schematic.png"), dpi=160); plt.close(fig)

# ================= figure 2: building the F-u law =================
fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
a = ax[0]
r = graziotti(210, 2580, W * 1e3 / 1292, O * 1e3 / 1292, 0.0, 1, 1, 105.0)
F0 = r["F0"] * b; us = r["uins"] * (2580 - 2060) / (2580 - r["h1"])
K = 5.7; uy = F0 / (K + F0 / us); Fy = K * uy
rc = graziotti(210, 2580, W * 1e3 / 1292, O * 1e3 / 1292, 0.15, 1, 1, 105.0); Fcr = rc["Fcr"] * b
p = np.array(doherty_trilinear(F0, us, "new")); ucr = Fcr / K
def draw_c(a, labels):
    a.plot([0, uy, us], [0, Fy, 0], color=BLUE, lw=2.4, label="bilinear: elastic K = 5.7 kN/mm, then the rigid line" if labels else None)
    a.plot(p[:, 0], p[:, 1], color=GREEN, lw=2, label="Doherty trilinear (new wall)" if labels else None)
    a.plot([0, ucr, ucr * 1.6], [0, Fcr, F0 * (1 - ucr * 1.6 / us)], color=ORANGE, lw=1.6, ls="--",
           label="bonded joints (f_w 0.15 MPa): crack at F_cr, drop to the rigid line" if labels else None)
    a.plot([0, us], [F0, 0], color=PURPLE, lw=1.6, ls=(0, (5, 3)), zorder=5,
           label="rigid-body line F0 (1 - u/u_ins)" if labels else None)
draw_c(a, True)
for (x, y, s, dx, dy) in ((0, F0, "F0", 1.5, -0.9), (us, 0, "u_ins", -7, 0.5),
                          (p[1, 0], p[1, 1], "(Δ1, F_i)", 1.5, -1.0), (p[2, 0], p[2, 1], "(Δ2, F_i)", 1.5, 0.4)):
    a.plot(x, y, "o", ms=4, color=INK, zorder=6); a.text(x + dx, y + dy, s, fontsize=8.5, color=INK)
a.set_xlim(0, 110); a.set_ylim(0, 15.5); a.grid(color=GRID)
a.set_xlabel("displacement at the 2.06 m sensor (mm)"); a.set_ylabel("horizontal force on the wall mass (kN)")
a.set_title("(c) Building the force-displacement law (h1 = 1.51 m, e = t/2)", loc="left", fontsize=11, fontweight="bold")
ins = a.inset_axes([0.12, 0.06, 0.28, 0.34])
draw_c(ins, False)
for (x, y, s, dx, dy) in ((uy, Fy, "(u_y, F_y)", 0.25, -1.4), (ucr, Fcr, "F_cr", 0.25, 0.2), (0, F0, "F0", 0.15, 0.3)):
    ins.plot(x, y, "o", ms=3.5, color=INK, zorder=6); ins.text(x + dx, y + dy, s, fontsize=8, color=INK)
ins.set_xlim(0, 8); ins.set_ylim(0, 15); ins.grid(color=GRID); ins.tick_params(labelsize=7)
ins.set_title("zoom 0-8 mm: the elastic branch", fontsize=8, loc="left")
a.legend(frameon=False, fontsize=8, loc="upper right")
a = ax[1]
for h1_, col in ((1.51, PURPLE), (1.8, BLUE), (2.0, ORANGE), (2.14, MUTED)):
    d, F = curve(h1_, W, O, thmax=0.2); m = F >= 0
    a.plot(d[m], F[m], color=col, lw=2, label="hinge at {:.2f} m: F0 {:.1f} kN".format(h1_, F[0]))
a.set_xlim(0, 110); a.set_ylim(0, 22); a.grid(color=GRID)
a.set_xlabel("displacement at the 2.06 m sensor (mm)"); a.set_ylabel("horizontal force on the wall mass (kN)")
a.set_title("(d) Hinge height: higher hinge = stronger and flatter at the sensor", loc="left", fontsize=11, fontweight="bold")
a.legend(frameon=False, fontsize=8.5, loc="upper right")
a.text(2, 1.0, "large-displacement virtual work (rigid_two_block_vw.py), e = t/2, latest build loads", fontsize=8, color=MUTED)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_mech_curves.png"), dpi=160); plt.close(fig)
print("-> fig_mech_schematic.png, fig_mech_curves.png")
