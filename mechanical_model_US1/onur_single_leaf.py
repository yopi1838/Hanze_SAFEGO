"""Onur Arslan's mechanical model (SafeGo_Mechanic_Model_v1.xlsx) for the US-1 solid wall: ONE leaf, no outer leaf,
no ties, two ways of treating the floors.

FLOORS MOVE (main model, `full_height`) -- the floors are not lateral supports, so the whole leaf spans base -> top
    beam as one Graziotti / Tomassetti three-hinge strip (Onur's treatment of his outer leaf, column C). The floor
    load on the wall is lumped into the top overburden; the joists slide with the wall, so there is no slip friction.
        F0, F_cr, h1, u_ins from graziotti();  u_ins mapped to the 2.06 m sensor by the strip kinematics.

FLOORS HELD (`split_at_floors`, Onur's inner leaf) -- the floors are lateral supports:
    upper block   Casapulla cantilever;  middle strip  Graziotti strip between the floors;
    lower block   Casapulla cantilever + joist slip friction F_Q = mu*Q (Onur's L7).

RB0: e = 0 (overburden at mid-thickness);  RBD: e = t/2 (overburden at the edge), u_ins,D = t at the hinge.
Units: kN, m (the strip formula per unit width in N/mm and mm, as in the sheet).
Run as a script: prints the US-1 cases, writes onur_single_leaf_US1.csv and draws fig_onur_single_leaf.png."""
import os, csv
import numpy as np
from rocking_model import graziotti, cantilever_a0

US1 = dict(t=0.21, b=1.292, h=2.58, W=13.16, O_top=10.35 + 1.485, Q_floor=1.18, y_lo=0.38, y_up=2.14, z_sensor=2.06)


def _strip(t, b, h, W, O, fw, e):
    return graziotti(t * 1e3, h * 1e3, W * 1e3 / (b * 1e3), O * 1e3 / (b * 1e3), fw, 1, 1, e * 1e3)


def full_height(t, b, h, W, O_top, Q_floor, z_sensor=2.06, fw=0.0, e_edge=False, **_):
    O = O_top + Q_floor                                   # upper-floor load lumped at the top (lower floor: < 0.1 kN, ignored)
    e = t / 2 if e_edge else 0.0
    r = _strip(t, b, h, W, O, fw, e)
    y_h = r["h1"] / 1e3
    u_h = t * 1e3 if e_edge else r["uins"]
    k_s = (h - z_sensor) / (h - y_h) if z_sensor > y_h else z_sensor / y_h
    return dict(model="floors move", F0=r["F0"] * b, F_cr=r["Fcr"] * b, y_hinge=y_h, u_ins_hinge_mm=u_h,
                u_ins_sensor_mm=u_h * k_s, F_strip=r["F0"] * b, F_upper=0.0, F_lower=0.0, F_Q=0.0)


def split_at_floors(t, b, h, W, O_top, Q_floor, y_lo, y_up, z_sensor=2.06, mu=0.6, fw=0.0, e_edge=False, joist_slip=True):
    hl, hs, hu = y_lo, y_up - y_lo, h - y_up
    Wl, Ws, Wu = W * hl / h, W * hs / h, W * hu / h
    O_strip = O_top + Wu + Q_floor                        # Onur B22
    O_lower = O_strip + Ws + Q_floor                      # Onur L6
    F_Q = mu * O_lower if joist_slip else 0.0             # Onur L7
    e = t / 2 if e_edge else 0.0
    r = _strip(t, b, hs, Ws, O_strip, fw, e)
    F_strip = r["F0"] * b
    F_up = cantilever_a0(Wu, O_top, t, hu, edge=e_edge) * Wu
    F_lo = cantilever_a0(Wl, O_lower, t, hl, FQ=F_Q, edge=e_edge) * Wl
    u_h = t * 1e3 if e_edge else r["uins"]
    y_h = y_lo + r["h1"] / 1e3
    k_s = (y_up - z_sensor) / (y_up - y_h) if z_sensor > y_h else (z_sensor - y_lo) / (y_h - y_lo)
    return dict(model="floors held", F0=F_strip + F_up + F_lo, F_cr=r["Fcr"] * b, y_hinge=y_h, u_ins_hinge_mm=u_h,
                u_ins_sensor_mm=u_h * k_s, F_strip=F_strip, F_upper=F_up, F_lower=F_lo, F_Q=F_Q)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    HERE = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for fn, nm in ((full_height, "floors move"), (split_at_floors, "floors held")):
        print("== US-1, one leaf, {} ==".format(nm))
        for fw in (0.0, 0.15):
            for edge in (False, True):
                slips = (True, False) if fn is split_at_floors else (False,)
                for slip in slips:
                    kw = dict(US1, fw=fw, e_edge=edge)
                    if fn is split_at_floors:
                        kw["joist_slip"] = slip
                    r = fn(**kw)
                    lab = "{}, {}{}, f_w {:.2f}".format(nm, "RBD" if edge else "RB0", ", F_Q on" if (fn is split_at_floors and slip) else "", fw)
                    print("  {:40s} F0 {:6.2f} kN (strip {:5.2f} + blocks {:4.2f}) | F_cr {:5.2f} | hinge {:.2f} m | u_ins {:4.0f} mm "
                          "at the hinge, {:5.1f} mm at 2.06 m".format(lab, r["F0"], r["F_strip"], r["F_upper"] + r["F_lower"], r["F_cr"],
                                                                       r["y_hinge"], r["u_ins_hinge_mm"], r["u_ins_sensor_mm"]))
                    rows.append(dict(case=lab, fw_MPa=fw, e="t/2" if edge else "0", **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
    with open(os.path.join(HERE, "onur_single_leaf_US1.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # ---------------- figure ----------------
    INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
    fig, ax = plt.subplots(1, 2, figsize=(16, 6.8), gridspec_kw=dict(width_ratios=[1.0, 1.45]))
    a = ax[0]; T = 0.42; g = US1
    rm = full_height(**US1); rs = split_at_floors(**US1)
    for x0, title, held in ((0.0, "floors MOVE (main)", False), (1.9, "floors HELD (Onur's split)", True)):
        if held:
            parts = ((0, g["y_lo"], "#d1fae5"), (g["y_lo"], g["y_up"], "#dbeafe"), (g["y_up"], g["h"], "#fde7d4"))
        else:
            parts = ((0, g["h"], "#dbeafe"),)
        for y0, y1, fc in parts:
            a.add_patch(Rectangle((x0, y0), T, y1 - y0, fc=fc, ec=INK, lw=1.3))
        for yj in (g["y_lo"], g["y_up"]):
            a.add_patch(Rectangle((x0 + T * 0.5, yj - 0.04), 0.75, 0.08, fc="#a16207", ec="none", alpha=0.75))
            if held:
                a.plot([x0 - 0.08, x0], [yj, yj], color=INK, lw=2)
                for k in range(4):
                    a.plot([x0 - 0.16, x0 - 0.08], [yj - 0.03 + 0.02 * k, yj - 0.01 + 0.02 * k], color=INK, lw=0.8)
        for y in (0.0, g["h"]):
            a.plot([x0 - 0.08, x0], [y, y], color=INK, lw=2)
        r = rs if held else rm
        a.plot([x0 - 0.05, x0 + T + 0.05], [r["y_hinge"]] * 2, color="#c2410c", ls="--", lw=1.2)
        a.text(x0 + T + 0.08, r["y_hinge"] - 0.12, "hinge {:.2f} m".format(r["y_hinge"]), color="#c2410c", fontsize=8.3)
        a.plot([x0 + T, x0 + T + 0.07], [g["z_sensor"]] * 2, color="#0f766e", lw=3)
        a.text(x0 + T / 2, g["h"] + 0.12, title, ha="center", fontsize=9.5, fontweight="bold")
        a.text(x0 + T / 2, -0.15, "F0 {:.1f} kN\nu_ins at 2.06 m: {:.0f} mm".format(r["F0"], r["u_ins_sensor_mm"]), ha="center", va="top", fontsize=8.8)
    a.text(0.0, -0.75, "Tests: the 2.06 m sensor reached about 29 mm (US-1 run 24)\nand the mode shape is straight to 2.06 m, so the floors move.",
           fontsize=8.6, color=INK, va="top")
    a.set_xlim(-0.4, 3.4); a.set_ylim(-1.15, 2.95); a.set_aspect("equal"); a.axis("off")
    a.set_title("(a) One leaf, two ways to treat the floors (RB0)", loc="left", fontsize=11, fontweight="bold")
    c = ax[1]
    for fw, edge, ls, col, lw in ((0.0, False, "-", INK, 2.2), (0.0, True, "--", INK, 2.0)):
        r = full_height(**US1, fw=fw, e_edge=edge)
        c.plot([0, r["u_ins_sensor_mm"]], [r["F0"], 0], color=col, ls=ls, lw=lw,
               label="floors move, {}: F0 {:.1f} kN, u_ins {:.0f} mm".format("RBΔ" if edge else "RB0", r["F0"], r["u_ins_sensor_mm"]))
    rb = full_height(**US1, fw=0.15)
    c.axhline(rb["F_cr"], color="#6250d6", ls=":", lw=1.6, label="cracking force with f_w = 0.15 MPa: {:.1f} kN".format(rb["F_cr"]))
    for edge, ls in ((False, "-"), (True, "--")):
        r = split_at_floors(**US1, e_edge=edge)
        c.plot([0, r["u_ins_sensor_mm"]], [r["F0"], 0], color="#9ca3af", ls=ls, lw=1.4,
               label="floors held, {}: F0 {:.1f} kN".format("RBΔ" if edge else "RB0", r["F0"]))
    dem = os.path.join(HERE, "pushover_v5_ft0_elastictop_pos.csv")
    if os.path.isfile(dem):                      # force the wall carries: |a| x m_tot (Appendix D.2), not F_tot as logged
        m_tot = float(open(dem).readline().split("m_tot")[1].split("kg")[0])
        d = np.genfromtxt(dem, delimiter=",", names=True, skip_header=1, dtype=None, encoding=None)
        c.plot(np.abs(d["d_ctrl_mm"]), np.abs(d["a_mps2"]) * m_tot / 1e3, color="#e07a2d", lw=2,
               label="DEM pushover f_t = 0 (+), |a| x m_tot")
    for fn, col, nm in (("envelope_peaks_US1.csv", "#2b5f8a", "US-1"), ("envelope_peaks_US2.csv", "#f0a020", "US-2")):
        pth = os.path.join(HERE, fn)
        if not os.path.isfile(pth):
            continue
        pts = []
        for rr in csv.DictReader(open(pth)):
            d0 = float(rr["d_start_mm"])
            for dk, fk in (("d_max_pos_mm", "F_at_dmax_pos_kN"), ("d_max_neg_mm", "F_at_dmax_neg_kN"),
                           ("d_at_Fmax_pos_mm", "F_max_pos_kN"), ("d_at_Fmax_neg_mm", "F_max_neg_kN")):
                pts.append((abs(float(rr[dk]) - d0), abs(float(rr[fk]))))
        pts = np.array(pts); c.plot(pts[:, 0], pts[:, 1], "o", ms=3.5, mfc="none", color=col, label="{} test, peak pairs per run".format(nm))
    c.set_xlim(0, 110); c.set_ylim(0, 24); c.grid(color=GRID)
    c.set_xlabel("displacement at the 2.06 m sensor (mm)"); c.set_ylabel("force (kN)")
    c.legend(frameon=False, fontsize=8, loc="upper right")
    c.set_title("(b) Rigid-body lines vs DEM and tests", loc="left", fontsize=11, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_onur_single_leaf.png"), dpi=150)
    print("-> onur_single_leaf_US1.csv, fig_onur_single_leaf.png")
