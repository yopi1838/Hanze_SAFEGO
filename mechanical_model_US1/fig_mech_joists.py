"""Joists in the rigid-body strip: illustration (fig_mech_joists_fbd.png) and curves (fig_mech_joists_curves.png)."""
import os, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Circle
import rigid_two_block_vw as vw

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
WOOD, CLAMP, LIFT, BEAR, INER = "#a16207", "#2a78d6", "#c2410c", "#6250d6", "#1baf7a"


def arr(ax, p0, p1, col, lw=1.8, ms=11):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, zorder=6))


def draw_wall(ax, side, title):
    """vw kinematics with an exaggerated thickness/rotation; side 'leeward' = push towards the joists (drawn to the
    right), 'windward' = push away (drawn mirrored, push to the left, joists still on the right)."""
    t_real = vw.t
    vw.t = 0.42; t = vw.t; h1 = 1.45; h = vw.h
    lo, up, th2 = vw.config(0.07, h1)
    s = 1.0 if side == "leeward" else -1.0                    # mirror for drawing
    X = lambda p: np.array([s * (p[0] if side == "leeward" else p[0]), p[1]]) if side == "leeward" else np.array([-p[0], p[1]])
    # undeformed outline
    xs = [0, t] if side == "leeward" else [0, -t]
    ax.add_patch(Rectangle((min(xs), 0), t, h, fill=False, ec=MUTED, ls=(0, (4, 3)), lw=0.9))
    for blk, rng, fc in ((lo, (0, h1), "#dbeafe"), (up, (h1, h), "#fde7d4")):
        P = np.array([X(blk(np.array(c))) for c in [(0, rng[0]), (t, rng[0]), (t, rng[1]), (0, rng[1])]])
        ax.add_patch(Polygon(P, fc=fc, ec=INK, lw=1.3, zorder=2))
    ax.add_patch(Rectangle((-0.6, -0.12), 1.2 + 0.0, 0.12, fc="#d1d5db", ec=INK, lw=1))
    # push
    for y in np.linspace(0.7, 1.9, 5):
        x0 = -0.75 if side == "leeward" else 0.75
        arr(ax, (x0, y), (x0 + (0.28 if side == "leeward" else -0.28), y), "#2a78d6", lw=1.0, ms=8)
    ax.text(-0.75 if side == "leeward" else 0.75, h + 0.05, "push", color="#2a78d6", ha="center", fontsize=9)
    # joists: bearing at the middle of the pocket on the joist side, far support L away
    L = 1.15
    for yj in (0.38, 2.14):
        blk = lo if yj <= h1 else up
        if side == "leeward":
            zb, zface, zf = t * 0.75, t, t + L
        else:
            zb, zface, zf = t * 0.25, 0.0, -L
        pb, pface, pfar = X(blk(np.array([zb, yj]))), X(blk(np.array([zface, yj]))), X(blk(np.array([zf, yj])))
        ax.plot([pb[0], pfar[0]], [pb[1], pfar[1]], color=WOOD, lw=6, solid_capstyle="butt", zorder=3, alpha=0.9)
        # far support (fixed, at the undeformed far end)
        sx = X(np.array([zf, yj]))[0]
        ax.add_patch(Rectangle((sx - 0.12, yj - 0.13), 0.24, 0.07, fc="#9ca3af", ec=INK, lw=0.8, zorder=1))
        rise = pfar[1] - yj
        if rise < 0:      # pressed: clamp reaction
            arr(ax, (sx, yj - 0.5), (sx, yj - 0.06), CLAMP, lw=2.6, ms=14)
            ax.text(sx + 0.1, yj - 0.33, "support pushes back:\nclamp moment k_r·θ\n(capped at M_cap)", color=CLAMP, fontsize=8.3, va="center")
        else:             # lifted off: the joist hangs off the wall
            arr(ax, (pfar[0] - 0.25 if side == "leeward" else pfar[0] + 0.25, pfar[1] + 0.02),
                (pfar[0] - 0.25 if side == "leeward" else pfar[0] + 0.25, pfar[1] - 0.38), LIFT, lw=2.2)
            ax.text(sx + 0.1, yj + 0.28, "lifted off: joist weight\nhangs off the wall\n(moment W_lever·L)", color=LIFT, fontsize=8.3, va="center")
        arr(ax, (pb[0], pb[1] + 0.32), (pb[0], pb[1] + 0.03), BEAR, lw=1.6, ms=9)
        ax.text(pb[0] + (0.05 if side == "leeward" else -0.05), pb[1] + 0.36, "R_wall", color=BEAR, fontsize=8,
                ha="left" if side == "leeward" else "right")
        # floor inertia
        xi = pface[0] + (0.35 if side == "leeward" else 0.35)
        arr(ax, (xi, yj + 0.12), (xi + (0.3 if side == "leeward" else -0.3), yj + 0.12), INER, lw=1.6, ms=9)
        ax.text(xi + (0.33 if side == "leeward" else -0.33), yj + 0.17, "m_f·a", color=INER, fontsize=8,
                ha="left" if side == "leeward" else "right")
    for p, lab in ((X(lo(np.array([t, 0.0]))), "B"), (X(lo(np.array([0.0, h1]))), "M"), (X(up(np.array([t, h]))), "P")):
        ax.add_patch(Circle(p, 0.04, fc="white", ec="#c2410c", lw=1.8, zorder=7))
        ax.text(p[0] + 0.07 * (1 if side == "leeward" else -1), p[1] + 0.06, lab, color="#c2410c", fontsize=9, fontweight="bold",
                ha="left" if side == "leeward" else "right")
    ax.set_xlim(-1.1, 2.2) if side == "leeward" else ax.set_xlim(-1.1, 2.2)
    ax.set_ylim(-0.75, h + 0.35); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    vw.t = t_real


fig = plt.figure(figsize=(17, 8))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.05)
a = fig.add_subplot(gs[0]); draw_wall(a, "leeward", "(a) Push towards the joists (+z in 3DEC)")
a.text(-1.05, -0.45, "lower floor PRESSES its support (clamp), upper floor LIFTS OFF", fontsize=9, color=INK)
b = fig.add_subplot(gs[1]); draw_wall(b, "windward", "(b) Push away from the joists (−z)")
b.text(-1.05, -0.45, "lower floor LIFTS OFF, upper floor PRESSES its support", fontsize=9, color=INK)
c = fig.add_subplot(gs[2])
v = np.linspace(-1, 1.4, 400); kv, Wl, fc = 1.0, 0.8, 0.5
for cap, col, lab in ((9, CLAMP, "clamped, no cap"), (fc, "#0f766e", "clamped, M_cap"), (0, MUTED, "pinned (k_r = 0)")):
    g = np.clip(kv * v, -cap, min(Wl, cap)) if cap > 0 else 0 * v
    c.plot(v, g, color=col, lw=2.2 if cap != 9 else 1.6, ls="-" if cap != 9 else "--", label=lab)
c.axhline(0, color=INK, lw=0.6); c.axvline(0, color=INK, lw=0.6)
c.annotate("lift-off: the joist weight\nW_lever hangs on the wall", (1.15, Wl), xytext=(0.2, 1.15), fontsize=8.5, color=LIFT,
           arrowprops=dict(arrowstyle="-", color=LIFT, lw=0.8))
c.annotate("slope k_v = k_r / L²", (0.35, 0.35), xytext=(-0.95, 0.75), fontsize=8.5, color=CLAMP,
           arrowprops=dict(arrowstyle="-", color=CLAMP, lw=0.8))
c.annotate("pocket turns: moment\ncapped at M_cap", (-0.8, -fc), xytext=(-0.95, -1.15), fontsize=8.5, color="#0f766e",
           arrowprops=dict(arrowstyle="-", color="#0f766e", lw=0.8))
c.set_xlabel("far-end rise v  (− pressed, + lifted)"); c.set_ylabel("restoring force at the far end")
c.set_xticks([]); c.set_yticks([]); c.set_ylim(-1.4, 1.4); c.grid(color=GRID)
c.legend(frameon=False, fontsize=8.5, loc="lower right")
c.set_title("(c) What the far support does", loc="left", fontsize=11, fontweight="bold")
fig.text(0.02, 0.02, "Virtual work, per floor j:  push  a·(M_wall/2 + Σ m_f,j ψ_j)·δΔ  =  strip terms  +  Σ R_wall,j·δv_bearing  "
         "+  Σ g_j(v)·δv_far          (ψ_j = sway of floor j per unit hinge displacement;  drawn with thickness and rotation exaggerated)",
         fontsize=9.5, color=INK)
fig.savefig(os.path.join(HERE, "fig_mech_joists_fbd.png"), dpi=150, bbox_inches="tight"); plt.close(fig)

# ======================= curves =======================
W, O = 13.15, 10.35 + 2.09           # overburden WITHOUT the joist share (now applied at the pockets)
O_lumped = O + 0.75
dem_pos = np.genfromtxt(os.path.join(HERE, "v5_pushover_digitised.csv"), delimiter=",", names=True) if os.path.isfile(
    os.path.join(HERE, "v5_pushover_digitised.csv")) else None
neg_csv = os.path.join(HERE, "pushover_neg_ft0.csv")
dem_neg = None
if os.path.isfile(neg_csv):
    rr = [r for r in csv.DictReader(open(neg_csv)) if r.get("converged", "1") != "0"]
    dem_neg = (np.array([abs(float(r["d_ctrl_mm"])) for r in rr]), np.array([abs(float(r["F_tot_kN"])) for r in rr]))
cases = [("bare strip (joist load lumped in O)", dict(O=O_lumped, joists=()), MUTED, ":", 1.6),
         ("pinned joists (test-like)", dict(k_r=0.0, M_cap=1e9), INK, "-", 2.0),
         ("clamped, M_cap 0.25 kN·m", dict(k_r=1.3e5, M_cap=0.25), "#0f766e", "-", 2.0),
         ("clamped, M_cap 1.0 kN·m", dict(k_r=1.3e5, M_cap=1.0), CLAMP, "-", 2.0),
         ("clamped, no cap", dict(k_r=1.3e5, M_cap=1e9), CLAMP, "--", 1.3),
         ("pinned + floor inertia", dict(k_r=0.0, M_cap=1e9, inertia=True), INER, "-.", 1.6)]
fig, ax = plt.subplots(2, 2, figsize=(14, 9.5), sharex=True)
for i, h1 in enumerate((1.51, 2.0)):
    for k, (side, lab) in enumerate((("leeward", "push towards the joists (+)"), ("windward", "push away from the joists (−)"))):
        a = ax[i, k]
        for name, kw, col, ls, lw in cases:
            if "joists" in kw:
                d, F = vw.curve(h1, W, kw["O"], thmax=0.12)
            else:
                d, F = vw.curve(h1, W, O, thmax=0.12, joists=vw.us1_joists(side, k_r=kw["k_r"], M_cap=kw["M_cap"]),
                                inertia=kw.get("inertia", False))
            m = F >= 0
            a.plot(d[m], F[m], color=col, ls=ls, lw=lw, label=name)
        dem = dem_pos if side == "leeward" else dem_neg
        if dem is not None:
            dx, dF = (dem["d_mm"], dem["F_kN"]) if side == "leeward" else dem
            a.plot(dx, dF, color="#e07a2d", lw=2.4, alpha=0.9, label="DEM pushover f_t = 0 ({})".format("+" if side == "leeward" else "−"))
        a.set_xlim(0, 100); a.set_ylim(0, 26); a.grid(color=GRID)
        a.set_title("hinge at {:.2f} m, {}".format(h1, lab), loc="left", fontsize=10.5, fontweight="bold")
        if k == 0: a.set_ylabel("force on the wall mass F = a·M_wall (kN)")
        if i == 1: a.set_xlabel("displacement at the 2.06 m sensor (mm)")
ax[0, 0].legend(frameon=False, fontsize=8.2, loc="upper right")
fig.suptitle("Rigid-body strip with the joists (US-1 DEM joists: 2 per floor, 406 kg, 76 kg on the wall; k_r 1.3e5 N·m/rad per floor)",
             fontweight="bold", fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(HERE, "fig_mech_joists_curves.png"), dpi=150); plt.close(fig)
print("-> fig_mech_joists_fbd.png, fig_mech_joists_curves.png")
