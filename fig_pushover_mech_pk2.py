# -*- coding: ascii -*-
"""fig_pushover_mech_pk2.py -- f_t = 0 pushovers (old frame, z-up, z-up pk2) against the rigid two-block mechanical
model (mechanical_model_US1/rigid_two_block_vw.py), with the joists at the OLD (0.38 / 2.14 m) and at the TEST pocket
heights (pk2: 0.3225 / 2.1989 m, pocket centres). Adjusted hinge 2.09 m, top pivot at the top-beam centroid, joists
pinned; band = lower pivot at the base toe (solid edge) to 0.30 m. Second band: + pocket friction (mu 0.75, lower pocket
clamped when pushed towards the joists) + elastic branch (K from T1 = 0.091 s, L 1200 kg).
Force: applied push |a| x m_masonry. Output: pushover_postproc_batch/fig_pushover_mech_pk2.png, pushover_mech_pk2.csv
Usage: python fig_pushover_mech_pk2.py"""
import os, sys, csv, io, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath("mechanical_model_US1")); import rigid_two_block_vw as vw
M_KG, H1, TB_Y, TB_KG, SPRING = 1308.3, 2.09, 2.629, 212.986, 10.35
MU, T1K, L = 0.75, 0.091, 1200.0
JOISTS = {"old": (0.38, 2.14), "pk2": (0.3225, 2.1989)}
OUT = "pushover_postproc_batch"; os.makedirs(OUT, exist_ok=True)
W = M_KG * 9.81 / 1e3; O = SPRING + TB_KG * 9.81 / 1e3; K = 4 * np.pi ** 2 * L / T1K ** 2 / 1e6


def joists(side, ys, mu):
    R, Wl = 0.186 * 2, (0.996 - 0.186) * 2
    Nl = vw.pocket_clamp_N(W, O, ys[0], side) if mu > 0 else 0.0
    return [vw.Joist(ys[0], R, Wl, k_r=0.0, side=side, mu=mu, N_clamp=Nl), vw.Joist(ys[1], R, Wl, k_r=0.0, side=side)]


def mech(side, ys, yb, mu, elastic):
    for thmax in (0.25, 0.2, 0.16, 0.12, 0.08):
        try:
            d, F = vw.curve(H1, W, O, joists=joists(side, ys, mu), thmax=thmax, n=600, top=(vw.t / 2, TB_Y), base_y=yb); break
        except ValueError:
            continue
    k = int(np.argmax(F < 0)) if (F < 0).any() else len(F)
    d, F = d[:k], F[:k]
    if elastic:
        d, F, _ = vw.with_elastic_branch(d, F, K)
    return d, F


def dem(p):
    Ls = open(p).read().splitlines()
    rows = list(csv.DictReader(io.StringIO("\n".join(l for l in Ls if l.strip() and not l.lstrip().lstrip('"').startswith("#")))))
    d = np.array([abs(float(r["d_ctrl_mm"])) for r in rows]); a = np.array([abs(float(r["a_mps2"])) for r in rows])
    ok = np.array([int(float(r["converged"])) == 1 for r in rows]); o = np.argsort(d[ok])
    return d[ok][o], (a[ok] * M_KG / 1e3)[o]


INK, MUT, GRID, SURF = "#0b0b0b", "#52514e", "#e7e6e2", "#fcfcfb"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
DEM = {"pos": [("DEM old frame, old pockets", "pushover_results_v5_ft0_masonrypush/pushover_pos.csv", C1),
               ("DEM z-up, old pockets", "pushover_results_v5zup_ft0_TopMason/pushover_pos.csv", C2),
               ("DEM z-up, pockets at test heights (pk2)", "pushover_results_v5zup_pk2_ft0_TopMason/pushover_pos.csv", C3)],
       "neg": [("DEM old frame, old pockets", "pushover_results_v5_ft0_masonrypush/pushover_neg.csv", C1),
               ("DEM z-up, pockets at test heights (pk2)", "pushover_results_v5zup_pk2_ft0_TopMason/pushover_neg.csv", C3)]}
rows = []
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True, facecolor=SURF)
for ax, br, side, title in ((axs[0], "pos", "leeward", "Positive push (towards the joists)"), (axs[1], "neg", "windward", "Negative push")):
    ax.set_facecolor(SURF)
    xs = np.linspace(0, 80, 161)
    for key, ls, alpha in (("pk2", "-", 0.22),):   # old heights: within 0.02 kN of pk2 (pinned joists act through weight only)
        for mu, el, fc, lab in ((0.0, False, "#8a8984", "rigid model"),
                                (MU, True, "#0f766e", "+ pocket friction + elastic branch" if br == "pos" else "+ elastic branch (no pocket clamp this way)")):
            if br == "neg" and mu > 0 and key == "old":
                continue
            c0 = mech(side, JOISTS[key], 0.0, mu, el); c3 = mech(side, JOISTS[key], 0.30, mu, el)
            f0 = np.interp(xs, *c0, right=np.nan); f3 = np.interp(xs, *c3, right=np.nan)
            name = "mech. {}".format(lab)
            if key == "pk2":
                ax.fill_between(xs, f0, f3, color=fc, alpha=alpha, lw=0, zorder=1)
                ax.plot(*c0, color=fc, lw=1.3, ls="-", zorder=2, label=name + ", lower pivot 0-0.30 m")
            else:
                ax.plot(*c0, color=fc, lw=1.1, ls=ls, zorder=2, label=name + ", lower pivot 0")
            for x in (20, 40, 80):
                rows.append(dict(branch=br, model=name, F_at_mm=x, pivot0_kN=round(float(np.interp(x, *c0)), 3),
                                 pivot030_kN=round(float(np.interp(x, *c3)), 3)))
    for lab, p, c in DEM[br]:
        d, F = dem(p); ax.plot(d, F, color=c, lw=2, zorder=3, label=lab)
        for x in (20, 40, 80):
            rows.append(dict(branch=br, model=lab, F_at_mm=x, pivot0_kN=round(float(np.interp(x, d, F)), 3) if x <= d.max() else "", pivot030_kN=""))
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_xlim(0, 82); ax.set_ylim(0, 14.5); ax.set_xlabel("displacement at the 2.06 m sensor (mm)", color=MUT)
    ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color("#b5b4ae")
    ax.tick_params(colors=MUT); ax.legend(fontsize=7.6, frameon=False, loc="lower left")
axs[0].set_ylabel("lateral force: applied push |a| x m_masonry (kN)", color=MUT)
fig.suptitle("f_t = 0 pushovers vs the rigid two-block model (hinge 2.09 m; shaded = lower pivot 0 to 0.30 m)",
             x=0.01, ha="left", color=INK, fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_pushover_mech_pk2.png"), dpi=150, facecolor=SURF)
with open(os.path.join(OUT, "pushover_mech_pk2.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["branch", "model", "F_at_mm", "pivot0_kN", "pivot030_kN"]); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
print("clamp N lower pocket, pk2 / old: {:.1f} / {:.1f} kN".format(vw.pocket_clamp_N(W, O, 0.3225, "leeward"), vw.pocket_clamp_N(W, O, 0.38, "leeward")))
