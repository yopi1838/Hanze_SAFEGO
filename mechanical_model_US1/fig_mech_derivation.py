"""Free-body illustration of the Graziotti / Doherty rigid-body strip derivation -> fig_mech_derivation.png
(drawn schematically: thickness exaggerated, not to scale)."""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, Circle

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
Q, WC, NV, OC, PIV, LEV = "#2a78d6", "#1f2937", "#6250d6", "#c2410c", "#c2410c", "#6b7280"
t, h1, h2 = 0.62, 1.55, 1.05      # drawing units (thickness exaggerated)
h = h1 + h2
zO = 0.82 * t                    # top contact = t/2 + e


def arr(ax, p0, p1, col, lw=1.6, ms=11, z=4):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, zorder=z))


def dim(ax, p0, p1, txt, off=(0, 0), col=LEV, fs=9, ha="center", va="center"):
    ax.annotate("", p0, p1, arrowprops=dict(arrowstyle="<->", color=col, lw=0.9, shrinkA=0, shrinkB=0))
    ax.text((p0[0] + p1[0]) / 2 + off[0], (p0[1] + p1[1]) / 2 + off[1], txt, color=col, fontsize=fs, ha=ha, va=va)


def hinge(ax, p, lab=None, dx=0.08, dy=0.0, ha="left"):
    ax.add_patch(Circle(p, 0.045, fc="white", ec=PIV, lw=2, zorder=6))
    if lab:
        ax.text(p[0] + dx, p[1] + dy, lab, color=PIV, fontsize=10, fontweight="bold", ha=ha, va="center", zorder=7)


def q_arrows(ax, y0, y1, n, x=0.0, L=0.32):
    for y in np.linspace(y0 + 0.07, y1 - 0.07, n):
        arr(ax, (x - L, y), (x - 0.01, y), Q, lw=1.1, ms=8)


fig = plt.figure(figsize=(18, 9.2))
gs = fig.add_gridspec(2, 4, width_ratios=[1.05, 1.15, 1.15, 1.25], height_ratios=[1, 0.62], wspace=0.08, hspace=0.18)

# ---------------- (a) the whole strip ----------------
a = fig.add_subplot(gs[:, 0])
a.add_patch(Rectangle((0, 0), t, h1, fc="#dbeafe", ec=INK, lw=1.4))
a.add_patch(Rectangle((0, h1), t, h2, fc="#fde7d4", ec=INK, lw=1.4))
a.add_patch(Rectangle((-0.15, -0.12), t + 0.3, 0.12, fc="#d1d5db", ec=INK, lw=1))
a.add_patch(Rectangle((-0.05, h), t + 0.1, 0.12, fc="#9ca3af", ec=INK, lw=1))
q_arrows(a, 0, h, 12)
a.text(-0.36, h + 0.05, "q per unit height\n(F = q·h, uniform\nover the mass)", color=Q, fontsize=9, ha="center", va="bottom")
arr(a, (zO, h + 0.55), (zO, h + 0.13), OC, lw=2)
a.text(zO + 0.05, h + 0.42, "O", color=OC, fontsize=11, fontweight="bold")
a.plot([t / 2, t / 2], [h + 0.13, h + 0.3], color=MUTED, lw=0.8, ls=":")
dim(a, (t / 2, h + 0.24), (zO, h + 0.24), "e", off=(0, 0.07), col=OC)
for (yc, lab) in ((h1 / 2, "W1"), (h1 + h2 / 2, "W2")):
    arr(a, (t / 2, yc + 0.17), (t / 2, yc - 0.17), WC, lw=1.6)
    a.text(t / 2 + 0.05, yc + 0.05, lab, fontsize=10, color=WC)
hinge(a, (t, 0), "B", dx=0.08, dy=0.1)
hinge(a, (0, h1), "M", dx=-0.07, dy=0.12, ha="right")
hinge(a, (zO, h), "P", dx=0.1, dy=-0.12)
a.plot([-0.05, t + 0.05], [h1, h1], color=PIV, lw=1, ls="--")
dim(a, (t + 0.22, 0), (t + 0.22, h1), "h1", off=(0.09, 0))
dim(a, (t + 0.22, h1), (t + 0.22, h), "h2 = h − h1", off=(0.09, 0), ha="left")
dim(a, (0, -0.25), (t, -0.25), "t", off=(0, -0.08))
a.text(0.0, -0.55, "windward\nface (z = 0)", fontsize=8, color=MUTED, ha="center", va="top")
a.text(t, -0.55, "leeward\nface (z = t)", fontsize=8, color=MUTED, ha="center", va="top")
a.text(-0.7, -1.0, "B: base pivot (leeward toe)\nM: mid hinge (windward face)\nP: top contact, z_O = t/2 + e",
       fontsize=8.8, color=INK, va="top")
a.set_xlim(-0.75, 1.45); a.set_ylim(-1.45, h + 0.75); a.set_aspect("equal"); a.axis("off")
a.set_title("(a) The strip and its loads", loc="left", fontsize=11.5, fontweight="bold")

# ---------------- (b) lower block FBD ----------------
b = fig.add_subplot(gs[:, 1])
b.add_patch(Rectangle((0, 0), t, h1, fc="#dbeafe", ec=INK, lw=1.4))
q_arrows(b, 0, h1, 7)
arr(b, (-0.62, h1 / 2), (-0.36, h1 / 2), Q, lw=2.4, ms=13)
b.text(-0.66, h1 / 2 + 0.1, "q·h1\nat h1/2", color=Q, fontsize=9, ha="right")
arr(b, (t / 2, h1 / 2 + 0.2), (t / 2, h1 / 2 - 0.2), WC, lw=1.8)
b.text(t / 2 + 0.05, h1 / 2 + 0.08, "W1", fontsize=10, color=WC)
# forces from the upper block at M: N down, V back (towards -z)
arr(b, (0, h1 + 0.42), (0, h1 + 0.03), NV, lw=2)
b.text(0.05, h1 + 0.33, "N = O + W2\n(upper block\npressing down)", color=NV, fontsize=8.5, va="center")
arr(b, (0.38, h1 + 0.0), (0.04, h1 + 0.0), NV, lw=2)
b.text(0.40, h1 - 0.1, "V", color=NV, fontsize=10, fontweight="bold")
hinge(b, (t, 0), "B", dx=0.08, dy=0.1)
hinge(b, (0, h1), "M", dx=-0.07, dy=-0.13, ha="right")
b.add_patch(Polygon([(t, 0), (t - 0.08, -0.12), (t + 0.08, -0.12)], fc="#9ca3af", ec=INK, lw=1))
# lever arms about B
dim(b, (0, -0.25), (t, -0.25), "t  (N about B)", off=(0, -0.09))
dim(b, (t / 2, -0.48), (t, -0.48), "t/2  (W1)", off=(0, -0.09))
dim(b, (t + 0.2, 0), (t + 0.2, h1), "h1  (V)", off=(0.09, 0), ha="left")
dim(b, (t + 0.55, 0), (t + 0.55, h1 / 2), "h1/2\n(q·h1)", off=(0.09, 0), ha="left")
b.text(-0.85, -0.85, "Moments about B:\n"
       r"$\mathbf{q\,h_1^2/2 \;=\; W_1\,t/2 \;+\; N\,t \;+\; V\,h_1}$" "\n"
       "push  =  own weight + upper block\n            + hinge force", fontsize=9.5, va="top", color=INK)
b.set_xlim(-0.95, 1.6); b.set_ylim(-1.45, h + 0.75); b.set_aspect("equal"); b.axis("off")
b.set_title("(b) Lower block: moments about B", loc="left", fontsize=11.5, fontweight="bold")

# ---------------- (c) upper block FBD ----------------
c = fig.add_subplot(gs[:, 2])
y0 = h1
c.add_patch(Rectangle((0, y0), t, h2, fc="#fde7d4", ec=INK, lw=1.4))
q_arrows(c, y0, y0 + h2, 5)
arr(c, (-0.62, y0 + h2 / 2), (-0.36, y0 + h2 / 2), Q, lw=2.4, ms=13)
c.text(-0.66, y0 + h2 / 2 + 0.1, "q·h2\nat h2/2", color=Q, fontsize=9, ha="right")
arr(c, (t / 2, y0 + h2 / 2 + 0.17), (t / 2, y0 + h2 / 2 - 0.17), WC, lw=1.8)
c.text(t / 2 - 0.05, y0 + h2 / 2 + 0.08, "W2", fontsize=10, color=WC, ha="right")
# at M (bottom windward corner): N up, V forward (+z)
arr(c, (0, y0 - 0.42), (0, y0 - 0.03), NV, lw=2)
c.text(0.05, y0 - 0.3, "N", color=NV, fontsize=10, fontweight="bold")
arr(c, (-0.36, y0), (-0.04, y0), NV, lw=2)
c.text(-0.40, y0 - 0.1, "V", color=NV, fontsize=10, fontweight="bold", ha="right")
# at P: O down, H_T unknown
arr(c, (zO, y0 + h2 + 0.45), (zO, y0 + h2 + 0.03), OC, lw=2)
c.text(zO + 0.05, y0 + h2 + 0.33, "O", color=OC, fontsize=11, fontweight="bold")
arr(c, (zO + 0.42, y0 + h2), (zO + 0.05, y0 + h2), MUTED, lw=1.4)
c.text(zO + 0.45, y0 + h2 + 0.02, "H_T (beam reaction,\nunknown, but passes\nthrough P)", color=MUTED, fontsize=8, va="center")
hinge(c, (0, y0), "M", dx=-0.07, dy=0.15, ha="right")
hinge(c, (zO, y0 + h2), "P", dx=0.03, dy=-0.17)
# lever arms about P
dim(c, (0, y0 - 0.62), (zO, y0 - 0.62), "z_O  (N about P)", off=(0.55, 0), ha="left")
c.plot([zO, zO], [y0 - 0.62, y0 + h2], color=MUTED, lw=0.6, ls=":")
dim(c, (t / 2, y0 + h2 / 2 - 0.32), (zO, y0 + h2 / 2 - 0.32), "z_O − t/2", off=(0.0, -0.09), fs=8)
dim(c, (t + 0.3, y0), (t + 0.3, y0 + h2), "h2  (V)", off=(0.09, 0), ha="left")
c.text(-0.85, -0.85, "Moments about P (H_T, O drop out):\n"
       r"$\mathbf{N\,z_O \;=\; q\,h_2^2/2 \;+\; V\,h_2 \;+\; W_2\,(z_O - t/2)}$" "\n"
       "Vertical equilibrium:   N = O + W2", fontsize=9.5, va="top", color=INK)
c.set_xlim(-0.95, 1.7); c.set_ylim(-1.45, h + 0.75); c.set_aspect("equal"); c.axis("off")
c.set_title("(c) Upper block: moments about P", loc="left", fontsize=11.5, fontweight="bold")

# ---------------- (d) displaced: lever arms shrink ----------------
d = fig.add_subplot(gs[0, 3])
D = 0.22
th1 = D / h1
Bp = np.array([t, 0.0])
def rot(p, c0, a):
    v = np.array(p) - c0
    return c0 + np.array([v[0] * np.cos(a) + v[1] * np.sin(a), -v[0] * np.sin(a) + v[1] * np.cos(a)])
L = [rot(p, Bp, th1) for p in [(0, 0), (t, 0), (t, h1), (0, h1)]]
Mp = L[3]; Pp = np.array([zO, h])
U = [Mp, Mp + (L[2] - L[3]), Pp + np.array([t - zO, 0]), Pp - np.array([zO, 0])]
d.add_patch(Rectangle((0, 0), t, h, fill=False, ec=MUTED, ls=(0, (4, 3)), lw=0.9))
d.add_patch(Polygon(L, fc="#dbeafe", ec=INK, lw=1.3)); d.add_patch(Polygon(U, fc="#fde7d4", ec=INK, lw=1.3))
hinge(d, Bp, "B", dx=0.08, dy=0.08); hinge(d, Mp, "M", dx=-0.08, dy=0.1, ha="right"); hinge(d, Pp, "P", dx=0.08, dy=0.1)
dim(d, (Mp[0], -0.2), (t, -0.2), "t − Δ", off=(0, -0.1), col=PIV)
d.plot([Mp[0], Mp[0]], [-0.2, Mp[1]], color=PIV, lw=0.6, ls=":")
dim(d, (0, h1 + 0.03), (Mp[0], h1 + 0.03), "Δ", off=(0, 0.1), col=NV)
d.text(1.0, h * 0.55, "Push the mid hinge by Δ:\nevery horizontal lever arm\nshrinks linearly with Δ\n"
       "  t → t − Δ      (N about B)\n  t/2 → t/2 − Δ/2  (W1)\n  z_O → z_O − Δ   (N about P)\n"
       "The top contact P does not move\n(beam held sideways).", fontsize=8.6, va="center", color=INK)
d.set_xlim(-0.4, 2.6); d.set_ylim(-0.55, h + 0.3); d.set_aspect("equal"); d.axis("off")
d.set_title("(d) Displaced by Δ", loc="left", fontsize=11.5, fontweight="bold")

e = fig.add_subplot(gs[1, 3])
F0, k = 1.0, 1.0 / 1.0
uu = np.linspace(0, 1.0, 10)
e.plot(uu, F0 - uu, color="#6250d6", lw=2.2)
e.plot([0], [F0], "o", color=INK); e.plot([1.0], [0], "o", color=INK)
e.text(0.03, F0 + 0.04, "F0  (Δ = 0)", fontsize=9.5, va="bottom")
e.text(0.98, 0.07, "u_ins = F0 / k", fontsize=9.5, ha="right")
e.text(0.5, 0.58, "slope −k", fontsize=9.5, color="#6250d6")
e.set_xlim(0, 1.08); e.set_ylim(0, 1.25); e.set_xticks([]); e.set_yticks([])
e.set_xlabel("hinge displacement Δ"); e.set_ylabel("lateral force F")
for s in ("top", "right"):
    e.spines[s].set_visible(False)
e.set_title("(e) The rigid-body line", loc="left", fontsize=11.5, fontweight="bold")

fig.text(0.02, 0.025,
         "Eliminate V between (b) and (c), use N = O + W2, W1 = W·h1/h, W2 = W·h2/h, z_O = t/2 + e:      "
         r"$F_0 = \frac{2}{h_1}(W+O)\,t \;+\; \frac{O\,(t+2e)}{h-h_1}$" "        "
         r"$k = \frac{2}{h_1}(W+O) + \frac{2\,O}{h-h_1}$" "        "
         r"$u_{ins} = F_0/k$" "        (drawn schematically, thickness exaggerated)",
         fontsize=10.5, color=INK)
fig.savefig(os.path.join(HERE, "fig_mech_derivation.png"), dpi=150, bbox_inches="tight")
print("-> fig_mech_derivation.png")
