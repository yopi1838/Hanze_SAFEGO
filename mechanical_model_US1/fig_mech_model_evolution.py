"""Figure for MECHANICAL_MODEL_US1.md: how the mechanical model evolved from Onur's cavity-wall model.
   (a) Onur's model for the cavity walls: inner leaf split at the floors (block + strip + block), outer leaf one strip,
       ties, joist slip friction; the floors are supports.
   (b) US-1, mid-height hinge: Onur's strip formulas on one leaf over the full height; floors move with the wall;
       hinge at the formula's optimum 1.51 m, top pivot on the wall's top corner.
   (c) US-1, adjusted hinge: the mechanism the DEM forms (and the test shows): hinge at the upper-floor course,
       top pivot at the top-beam centroid (elastic top: the beam turns with the upper block), lower pivot between
       the base and 0.30 m (kink at the lower floor). Orange dots: the DEM course displacements at 80 mm (pos),
       scaled to the drawing.
   (d) What changed.
Drawings: thickness and sway exaggerated (XS), rotations exaggerated; not to scale."""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Circle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rigid_two_block_vw as vw

INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
LOW, UPP, STRIP, OUTER = "#cfe0f1", "#fbe0c8", "#dbeafe", "#e5e7eb"
HINGE, DEMC, WOOD, BEAM, BLUEL = "#c2410c", "#E07A2D", "#a16207", "#9ca3af", "#2B5F8A"
t, h, zs = vw.t, vw.h, vw.zs
XS = 2.4                       # horizontal exaggeration
TB_Y, TB_H = 2.629, 0.10       # top-beam centroid and drawn height
Y_LO, Y_UP = 0.38, 2.14


def X(z, x0=0.0):
    return x0 + np.asarray(z) * XS


def arrow(ax, p0, p1, col=INK, lw=1.3, ms=9, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, **kw))


def poly(ax, pts, x0, **kw):
    P = np.array(pts)
    ax.add_patch(Polygon(np.c_[X(P[:, 0], x0), P[:, 1]], closed=True, **kw))


def ground(ax, x0, x1, y=0.0):
    ax.add_patch(Rectangle((x0, y - 0.10), x1 - x0, 0.10, fc="#d1d5db", ec=INK, lw=1))


def pivot(ax, p, x0, col=INK, r=0.028):
    ax.add_patch(Circle((X(p[0], x0), p[1]), r, fc="white", ec=col, lw=1.6, zorder=6))


def held(ax, x, y, side=1, n=4, col=INK):
    """support hatching: wall held sideways at (x, y)"""
    ax.plot([x, x], [y - 0.06, y + 0.06], color=col, lw=2)
    for k in range(n):
        yy = y - 0.06 + 0.04 * k
        ax.plot([x, x + side * 0.06], [yy, yy + 0.03], color=col, lw=0.8)


def push(ax, x0, y0, y1, n=7, L=0.16):
    for y in np.linspace(y0 + 0.12, y1 - 0.12, n):
        arrow(ax, (x0 - L - 0.05, y), (x0 - 0.05, y), col=MUTED, lw=1.0, ms=7)


def label(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 8.3); kw.setdefault("color", INK)
    ax.text(x, y, s, **kw)


fig = plt.figure(figsize=(20, 8.6))
gs = fig.add_gridspec(1, 4, width_ratios=[3.05, 2.2, 2.4, 4.3], wspace=0.04)
axs = [fig.add_subplot(gs[0, i]) for i in range(4)]
for a in axs:
    a.set_aspect("equal"); a.axis("off")

# ============================ (a) Onur ============================
a = axs[0]; x_in, x_out = 1.05, 0.0
ground(a, -0.35, X(t, x_in) + 0.75)
# outer leaf: one full-height strip, hinge at its optimum
a.add_patch(Rectangle((x_out, 0), X(t), h, fc=OUTER, ec=INK, lw=1.2))
a.plot([x_out - 0.04, X(t, x_out) + 0.04], [1.29, 1.29], color=HINGE, ls="--", lw=1.2)
label(a, x_out + X(t) / 2, h + 0.08, "outer\nleaf", ha="center", fontsize=8.5)
# inner leaf: split at the floors
for y0, y1, fc in ((0, Y_LO, "#d9f2e4"), (Y_LO, Y_UP, STRIP), (Y_UP, h, UPP)):
    a.add_patch(Rectangle((x_in, y0), X(t), y1 - y0, fc=fc, ec=INK, lw=1.2))
a.plot([x_in - 0.04, X(t, x_in) + 0.04], [1.36, 1.36], color=HINGE, ls="--", lw=1.2)
label(a, x_in + X(t) / 2, h + 0.08, "inner\nleaf", ha="center", fontsize=8.5)
# ties
for y in (0.7, 1.2, 1.7, 2.2):
    xx = np.linspace(X(t, x_out), x_in, 9); yy = y + 0.03 * np.where(np.arange(9) % 2, 1, -1); yy[0] = yy[-1] = y
    a.plot(xx, yy, color=MUTED, lw=0.9)
label(a, (X(t, x_out) + x_in) / 2, 2.38, "ties\nF_c", ha="center", fontsize=8, color=MUTED)
# floors held
for y in (Y_LO, Y_UP):
    a.add_patch(Rectangle((X(t, x_in), y - 0.035), 0.55, 0.07, fc=WOOD, ec="none", alpha=0.85))
    held(a, X(t, x_in) + 0.62, y, side=1)
label(a, X(t, x_in) + 0.08, Y_UP - 0.30, "floor =\nsupport", fontsize=8)
arrow(a, (X(t, x_in) + 0.35, Y_LO + 0.12), (X(t, x_in) + 0.05, Y_LO + 0.12), col=HINGE, lw=1.1, ms=7)
label(a, X(t, x_in) + 0.38, Y_LO + 0.08, "joist slip F_Q = μQ", fontsize=8, color=HINGE)
label(a, X(t, x_in) + 0.08, (Y_LO + Y_UP) / 2, "Graziotti\nstrip", fontsize=8.3, color=INK)
label(a, X(t, x_in) + 0.08, 2.42, "Casapulla block", fontsize=8, color=INK)
label(a, X(t, x_in) + 0.08, 0.10, "Casapulla block", fontsize=8, color=INK)
push(a, x_out, 0, h, n=8)
label(a, -0.35, -0.42, "F0 = strip + blocks (inner) + strip (outer) + F_c\nfloors held: the leaf only moves between them",
      fontsize=8.6, va="top")

a.set_xlim(-0.6, 2.45); a.set_ylim(-0.95, 3.15)

# ======================= (b) mid-height hinge ======================
a = axs[1]; x0 = 0.35; h1 = 1.51; th1 = 0.05
lo, up, th2 = vw.config(th1, h1)
ground(a, x0 - 0.35, X(t, x0) + 0.75)
a.add_patch(Rectangle((x0, 0), X(t), h, fill=False, ls=(0, (4, 3)), lw=0.9, ec=MUTED))
poly(a, [lo(np.array(c)) for c in ((0, 0), (t, 0), (t, h1), (0, h1))], x0, fc=LOW, ec=INK, lw=1.3)
poly(a, [up(np.array(c)) for c in ((0, h1), (t, h1), (t, h), (0, h))], x0, fc=UPP, ec=INK, lw=1.3)
ptop = up(np.array([t, h]))
a.add_patch(Rectangle((X(ptop[0] - t, x0) - 0.03, ptop[1]), X(t) + 0.06, TB_H, fc=BEAM, ec=INK, lw=1))   # beam: translates up
held(a, X(ptop[0], x0) + 0.06, ptop[1] + 0.05, side=1)
arrow(a, (X(t / 2, x0), ptop[1] + TB_H + 0.30), (X(t / 2, x0), ptop[1] + TB_H + 0.02), lw=1.3)
label(a, X(t / 2, x0) + 0.08, ptop[1] + TB_H + 0.22, "O", fontsize=9)
for p in (lo(np.array([t, 0.0])), lo(np.array([0.0, h1])), ptop):
    pivot(a, p, x0, col=HINGE)
a.plot([X(-0.05, x0), X(t, x0) + 0.62], [h1, h1], color=HINGE, ls="--", lw=0.9)
label(a, X(t, x0) + 0.30, h1 + 0.04, "hinge 1.51 m", color=HINGE)
label(a, X(t, x0) + 0.30, ptop[1] - 0.12, "top corner", color=HINGE)
label(a, X(t, x0) + 0.30, 0.04, "base toe", color=HINGE)
for y in (Y_LO, Y_UP):
    blk = lo if y <= h1 else up
    pl, pr = blk(np.array([t, y])), blk(np.array([t, y]))
    a.add_patch(Rectangle((X(pr[0], x0), pr[1] - 0.035), 0.42, 0.07, fc=WOOD, ec="none", alpha=0.85))
    arrow(a, (X(pr[0], x0) + 0.44, pr[1]), (X(pr[0], x0) + 0.60, pr[1]), col=WOOD, lw=1.0, ms=7)
label(a, X(t, x0) + 0.10, Y_UP + 0.17, "floors move", fontsize=8, color=WOOD)
push(a, x0, 0, h, n=8)
label(a, x0 - 0.35, -0.42, "Onur's strip formulas, one leaf, full height\n"
      "F0 = 9.8 (e = 0) to 12.4 kN (e = t/2)\nzero force at 100 mm (sensor)", fontsize=8.6, va="top")

a.set_xlim(-0.25, 1.95); a.set_ylim(-0.95, 3.15)

# ========================= (c) adjusted hinge ======================
a = axs[2]; x0 = 0.35; h1 = 2.09; th1 = 0.04
TOP = (t / 2, TB_Y)
lo, up, th2 = vw.config(th1, h1, top=TOP)
ground(a, x0 - 0.35, X(t, x0) + 0.85)
a.add_patch(Rectangle((x0, 0), X(t), h, fill=False, ls=(0, (4, 3)), lw=0.9, ec=MUTED))
a.add_patch(Rectangle((x0 - 0.03, h), X(t) + 0.06, TB_H, fill=False, ls=(0, (4, 3)), lw=0.9, ec=MUTED))
poly(a, [lo(np.array(c)) for c in ((0, 0), (t, 0), (t, h1), (0, h1))], x0, fc=LOW, ec=INK, lw=1.3)
poly(a, [up(np.array(c)) for c in ((0, h1), (t, h1), (t, h), (0, h))], x0, fc=UPP, ec=INK, lw=1.3)
poly(a, [up(np.array(c)) for c in ((-0.03 / XS, h), (t + 0.03 / XS, h), (t + 0.03 / XS, h + TB_H), (-0.03 / XS, h + TB_H))],
     x0, fc=BEAM, ec=INK, lw=1)                                                   # beam turns with the upper block
pt = up(np.array(TOP))
held(a, X(pt[0], x0) + X(t) / 2 + 0.07, pt[1], side=1)
arrow(a, (X(pt[0], x0), pt[1] + 0.38), (X(pt[0], x0), pt[1] + 0.10), lw=1.3)
label(a, X(pt[0], x0) + 0.08, pt[1] + 0.28, "O", fontsize=9)
for p in (lo(np.array([t, 0.0])), lo(np.array([0.0, h1])), pt):
    pivot(a, p, x0, col=HINGE)
a.plot([X(-0.05, x0), X(t, x0) + 0.70], [h1, h1], color=HINGE, ls="--", lw=0.9)
label(a, X(t, x0) + 0.22, h1 - 0.36, "hinge 2.09 m\n(upper-floor course)", color=HINGE)
label(a, X(t, x0) + 0.33, pt[1] + 0.10, "top-beam centroid\n(elastic top)", color=HINGE)
# lower pivot range
a.plot([X(t, x0) + 0.10] * 2, [0.0, 0.30], color=BLUEL, lw=2.4)
label(a, X(t, x0) + 0.33, 0.02, "lower pivot:\nbase … 0.30 m", color=BLUEL)
for y in (Y_LO, Y_UP):
    blk = lo if y <= h1 else up
    pr = blk(np.array([t, y]))
    a.add_patch(Rectangle((X(pr[0], x0), pr[1] - 0.035), 0.42, 0.07, fc=WOOD, ec="none", alpha=0.85))
    a.add_patch(Circle((X(pr[0], x0) + 0.02, pr[1]), 0.022, fc="white", ec=WOOD, lw=1.2, zorder=6))
label(a, X(t, x0) + 0.08, Y_LO + 0.08, "joists pinned,\nbearing load in pocket", fontsize=7.8, color=WOOD)
# DEM course displacements, scaled to the drawing at the sensor
try:
    rows = list(csv.DictReader(open(os.path.join(HERE, "dem_pos_end_courses.csv"))))
    y0 = np.array([float(r["y0"]) for r in rows]); dz = np.array([float(r["dz_mm"]) for r in rows])
    d_draw = (lo(np.array([t / 2, zs])) if zs <= h1 else up(np.array([t / 2, zs])))[0] - t / 2
    k = d_draw / (np.interp(zs, y0, dz) / 1e3)
    a.plot(X(t / 2 + dz / 1e3 * k, x0), y0, "o", ms=2.6, color=DEMC, zorder=7)
    label(a, x0 - 0.33, 1.05, "DEM\ncourses", fontsize=7.8, color=DEMC)
except Exception:
    pass
push(a, x0, 0, h, n=8)
label(a, x0 - 0.35, -0.42, "same strip, mechanism taken from the DEM / test\n"
      "F0 = 9.9 to 11.8 kN (pos), 10.7 to 12.6 kN (neg)\nzero force at 127-145 mm (sensor)", fontsize=8.6, va="top")

a.set_xlim(-0.25, 2.15); a.set_ylim(-0.95, 3.15)

# ========================= (d) what changed ========================
a = axs[3]; a.set_aspect("auto"); a.set_xlim(0, 1); a.set_ylim(0, 1)

rows = [("", "Onur (a)", "mid-height hinge (b)", "adjusted hinge (c)"),
        ("leaves", "inner + outer,\nties F_c", "one leaf, no ties", "same"),
        ("floors", "supports; leaf split\nat them; slip F_Q", "move with the wall;\nno split, no F_Q", "same; joists pinned,\nbearing load in pockets"),
        ("mid hinge", "strip optimum\n(between floors)", "formula optimum\n1.51 m", "upper-floor course\n2.09 / 2.20 m (DEM, test)"),
        ("top pivot", "corner of the leaf\n(e = 0 or t/2)", "wall's top corner\n(e = 0 or t/2)", "top-beam centroid,\n5 cm above the wall"),
        ("lower pivot", "base of each part", "base toe", "base … 0.30 m\n(lower-floor kink)"),
        ("overburden", "0.1 MPa + floor\nload (3.84 kN)", "springs + beam +\nfloor, at the top", "springs + 213 kg beam\nat the pivot"),
        ("F–u law", "linear line\nF0 (1 − u/u_ins)", "linear line, u_ins\nmoved to the sensor", "large-displacement\nvirtual work"),
        ("checked\nagainst", "cavity tests", "tests' peaks", "DEM pushover\n(|a|·m_tot)")]
xc = [0.0, 0.15, 0.43, 0.71]
y = 0.97; rh = [0.055] + [0.098] * (len(rows) - 1)
for i, r in enumerate(rows):
    if i % 2 == 1:
        a.add_patch(Rectangle((0, y - rh[i]), 1, rh[i], fc="#f3f4f6", ec="none"))
    for j, c in enumerate(r):
        bold = i == 0 or j == 0
        col = HINGE if (j == 3 and i in (3, 4, 5)) else INK
        a.text(xc[j] + 0.01, y - rh[i] / 2, c, fontsize=8.6, va="center", ha="left", color=col,
               fontweight="bold" if bold else "normal")
    y -= rh[i]
a.plot([0, 1], [0.97 - rh[0]] * 2, color=INK, lw=0.8)
a.text(0.01, y - 0.03, "Red: what the adjusted hinge takes from the DEM deformed shape (Appendix D).\n"
       "Drawings: thickness and sway exaggerated, not to scale.", fontsize=8, color=MUTED, va="top")
fig.canvas.draw()
ytop = max(ax_.get_position().y1 for ax_ in axs[:3]) + 0.01
ref = axs[2].get_position(); p3 = axs[3].get_position()
axs[3].set_position([p3.x0, ref.y0 - 0.02, p3.width, ytop - 0.03 - (ref.y0 - 0.02)])
for ax_, ttl in zip(axs, ("(a) Onur: cavity wall, split at the floors", "(b) US-1, mid-height hinge",
                         "(c) US-1, adjusted hinge", "(d) What changed")):
    fig.text(ax_.get_position().x0 + 0.005, ytop, ttl, fontsize=11.5, fontweight="bold", va="bottom")
fig.savefig(os.path.join(HERE, "fig_mech_model_evolution.png"), dpi=150, bbox_inches="tight")
print("-> fig_mech_model_evolution.png")
