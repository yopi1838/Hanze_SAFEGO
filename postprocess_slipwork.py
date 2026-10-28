# -*- coding: ascii -*-
"""
postprocess_slipwork.py -- slip work per contact region (run_pushover_slipwork.py) against the mechanical model.

Reads <DIR>/slipwork_pos.csv and slipwork_neg.csv (slip_work_tracker.py) and draws, per direction:
    (a) slip work per region vs the sensor displacement (J)
    (b) the lower pocket's clamp: sum of its normal forces (kN) vs d, with the model's N
        (rigid_two_block_vw.pocket_clamp_N: about 23.4 kN pushed towards the joists, 0 the other way)
    (c) the lower pocket's slip work vs d, with the model's mu N (0.088 m) |dphi(d)| for the 0 and 0.30 m pivots
Writes slipwork_summary.csv (values at 20 / 40 / 80 mm) and fig_slipwork.png into DIR.
Usage: python postprocess_slipwork.py [DIR]     (default pushover_results_v5_ft_0p050_TopMason_slipwork)
"""
import os, sys, csv, io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR      = sys.argv[1] if len(sys.argv) > 1 else "pushover_results_v5_ft_0p050_TopMason_slipwork"
MECH_DIR = "mechanical_model_US1"
MU       = 0.75
H1, TB_Y = 2.09, 2.629
M_TOT, SPRING_KN, TB_KG = 1308.3, 10.35, 212.986
REGIONS  = [("pocket_low", "lower pocket", "#0F766E"), ("pocket_up", "upper pocket", "#14B8A6"),
            ("mas_lowfloor", "joints 0.25-0.55 m", "#2B5F8A"), ("mas_upfloor", "joints 2.03-2.25 m (hinge)", "#7C3AED"),
            ("base", "base joint", "#CA8A04"), ("top", "top joint", "#DC2626"),
            ("mas_other", "other joints", "#9CA3AF"), ("support", "joist supports (shear only)", "#6B7280"),
            ("new", "contacts formed later", "#DB2777")]
INK, GRID = "#1f2937", "#e5e7eb"


def read(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(io.StringIO("".join(l for l in f if not l.startswith("#")))))
    out = {k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows])
           for k in rows[0] if k != "tag"} if rows else None
    if out:
        out["tag"] = np.array([r["tag"] for r in rows])
        out["d"] = np.abs(out["d_ctrl_mm"])
    return out


def model_curves(side):
    """d (mm) and the lower-pocket friction work (J) of the mechanical model, for the two lower pivots"""
    sys.path.insert(0, os.path.abspath(MECH_DIR))
    import rigid_two_block_vw as vw
    W = M_TOT * 9.81 / 1e3; O = SPRING_KN + TB_KG * 9.81 / 1e3
    N = vw.pocket_clamp_N(W, O, 0.38, side)
    j = vw.us1_joists(side, k_r=0.0, mu=MU, W=W, O=O)[0]
    res = {}
    for yb in (0.0, 0.30):
        th = np.linspace(1e-5, 0.12, 300); d = []; wf = []
        for a in th:
            lo, up, th2 = vw.config(a, H1, (vw.t / 2, TB_Y), yb)
            pb = lo(np.array([j.zb(), j.y]))
            dphi = a - (pb[1] - j.y) / j.L * (1.0 if side == "leeward" else -1.0)
            d.append(vw.state(a, H1, W, O, joists=(), top=(vw.t / 2, TB_Y), base_y=yb)[2] * 1e3)
            wf.append(MU * N * j.slip_lever * abs(dphi) * 1e3)
        res[yb] = (np.array(d), np.array(wf))
    return N, res


def main():
    data = {br: read(os.path.join(DIR, "slipwork_{}.csv".format(br))) for br in ("pos", "neg")}
    if not any(v is not None for v in data.values()):
        raise SystemExit("no slipwork_*.csv in " + DIR)
    fig, axs = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    summ = []
    for j, (br, side, title) in enumerate((("pos", "leeward", "Positive push (towards the joists)"),
                                           ("neg", "windward", "Negative push"))):
        c = data[br]
        ax_a, ax_b, ax_c = axs[0, j], axs[1, j], axs[2, j]
        ax_a.set_title(title, loc="left", fontweight="bold", color=INK)
        if c is None:
            for ax in (ax_a, ax_b, ax_c):
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", color="#6b7280")
            continue
        o = np.argsort(c["d"]); d = c["d"][o]
        for key, lab, col in REGIONS:
            w = c[key + "_Wslip_J"][o]
            if np.nanmax(np.abs(w)) > 0.5:
                ax_a.plot(d, w, color=col, lw=1.8, label=lab)
        tot = sum(c[k + "_Wslip_J"][o] for k, _, _ in REGIONS)
        ax_a.plot(d, tot, color=INK, lw=1.2, ls="--", label="all regions")
        try:
            N, mod = model_curves(side)
        except Exception as e:
            print("  ! mechanical model not drawn: {}".format(e)); N, mod = None, {}
        ax_b.plot(d, c["pocket_low_N_kN"][o], color="#0F766E", lw=1.8, label="DEM: lower pocket, sum of normal forces")
        ax_b.plot(d, c["pocket_up_N_kN"][o], color="#14B8A6", lw=1.2, label="DEM: upper pocket")
        if N is not None:
            ax_b.axhline(N, color=INK, ls=":", lw=1.4, label="model clamp N = {:.1f} kN".format(N))
        ax_c.plot(d, c["pocket_low_Wslip_J"][o], color="#0F766E", lw=2.0, label="DEM: lower pocket slip work")
        for yb, ls in ((0.0, "--"), (0.30, "-")):
            if yb in mod:
                ax_c.plot(mod[yb][0], mod[yb][1], color=INK, ls=ls, lw=1.2,
                          label=r"model $\mu N$ (0.088 m) $|\Delta\phi|$, lower pivot " + "{:.2f} m".format(yb))
        for x in (20, 40, 80):
            if d.size and x <= d.max():
                r = dict(branch=br, d_mm=x, pocket_low_N_kN=float(np.interp(x, d, c["pocket_low_N_kN"][o])))
                for key, _, _ in REGIONS:
                    r[key + "_Wslip_J"] = float(np.interp(x, d, c[key + "_Wslip_J"][o]))
                for yb in mod:
                    r["model_pocket_low_W_J_pivot{:.2f}".format(yb)] = float(np.interp(x, *mod[yb]))
                summ.append(r)
        ax_c.set_xlim(0, max(85.0, 1.05 * float(np.nanmax(d))) if d.size else 85.0)
        ax_a.set_ylabel("slip work (J)"); ax_b.set_ylabel("normal force (kN)"); ax_c.set_ylabel("lower pocket slip work (J)")
        ax_c.set_xlabel("displacement at the 2.06 m sensor (mm)")
        for ax in (ax_a, ax_b, ax_c):
            ax.grid(color=GRID); ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Slip work per contact region (slip_work_tracker) vs the mechanical model's pocket friction", x=0.01,
                 ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(DIR, "fig_slipwork.png"), dpi=150)
    if summ:
        keys = []
        for r in summ:
            keys += [k for k in r if k not in keys]
        with open(os.path.join(DIR, "slipwork_summary.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(keys)
            for r in summ:
                w.writerow([round(r[k], 3) if isinstance(r.get(k), float) else r.get(k, "") for k in keys])
        for r in summ:
            print("{} {:3d} mm: lower pocket N {:6.2f} kN, W {:7.1f} J (model {}) | upper pocket {:6.1f} J | hinge course "
                  "{:6.1f} J | base {:6.1f} J | top {:6.1f} J".format(
                      r["branch"], r["d_mm"], r["pocket_low_N_kN"], r["pocket_low_Wslip_J"],
                      " / ".join("{:.1f}".format(r[k]) for k in r if k.startswith("model_")) or "-",
                      r["pocket_up_Wslip_J"], r["mas_upfloor_Wslip_J"], r["base_Wslip_J"], r["top_Wslip_J"]))
    print("-> {}".format(os.path.join(DIR, "fig_slipwork.png")))


if __name__ == "__main__":
    main()
