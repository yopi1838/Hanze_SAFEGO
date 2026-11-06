# -*- coding: utf-8 -*-
"""compare_pushover_test_mech.py -- the pk2 f_t 0.05 pushovers (glued pocket, pK11 pocket) against
  (a) the tests' WALL-ONLY base shear: the group's formula (Ihsan / Eleni, 26 Oct 2026) without its slab term,
      V_wall = 12.94 kN x (f1 a15 + f2 a16 + f3 a17), f from each test's sensor heights (US-1 0.66/1.26/2.06 m,
      US-2 0.60/1.26/2.06 m, read from TestX_Info.xlsx). That is the load a masonry-only uniform push applies.
      The full formula (with 165 kg on a13 and a14) is drawn faintly for reference.
      Points: force at each run's peak displacement (+ and -), 10 Hz zero-phase low-pass, displacement from run start.
      Sign: the test's + is 0.5(ch3+ch4) as recorded, force sign set so that F follows d (mapping of the test's + to the
      model's + = towards the joists is NOT verified; the tests are nearly symmetric).
  (b) the rigid three-hinge model (mechanical_model_US1/rigid_two_block_vw.py, Appendix D of MECHANICAL_MODEL_US1.md)
      with this build's loads: masonry 1308.3 kg, overburden O = springs 10.35 kN + top beam 213 kg, joists pinned
      (bearing load only), pk2 joist heights 0.32 / 2.20 m. Hinge: the weakest height at every d (what a rigid wall
      picks), or held at the upper pocket (2.11 m, where the DEM and US-1 crack). Top: rocking on the top-beam corner
      (mason top joint) or about the beam centroid (0.105, 2.63).
        "floors NOT pushed": the floors ride along without inertia -> like the masonry-only pushover
        "floors pushed"    : the floor masses (203 kg per floor) accelerate with the wall; F is still the force on the
                             wall mass -> like the tests' wall-only base shear (the floors' own inertia is the slab term)
      Each drawn with the elastic branch K = the pK11 pushover's initial stiffness.
Run from the Hanze_SAFEGO root:  python compare_pushover_test_mech.py   (needs EXP_DATA/ and mechanical_model_US1/)
Output: pushover_pocket_compare/fig_pushover_test_mech.png, test_envelope_wall_only.csv, mech_curves.csv. Linear axes."""
import os, sys, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "mechanical_model_US1")
import rigid_two_block_vw as rb

RAW = "EXP_DATA"
OUT = "pushover_pocket_compare"
LP_HZ = 10.0
G = 9.8104
H = 2.58
WALL_KN = 12.94
SLAB_KG = 165.0
CASES = [("DEM, glued pocket", "pushover_results_v5zup_pk2_ft0p05_TopMason", "#0F766E"),
         ("DEM, pK11 pocket", "pushover_results_v5zup_pk2_ft0p05_pK11_TopMason", "#B91C1C")]
M_MASONRY = 1308.3
O_KN = 10.35 + 213.0 * 9.81 / 1e3
FLOOR_KG = 203.0


def heights(test):
    import openpyxl
    ws = openpyxl.load_workbook(os.path.join(RAW, "{}_Info.xlsx".format(test)), read_only=True, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    return tuple(float(rows[r][5]) for r in (15, 16, 17))


def fractions(h):
    h1, h2, h3 = h
    return ((h1 + (h2 - h1) / 2) / H, ((h2 - h1) / 2 + (h3 - h2) / 2) / H, (H - h3 + (h3 - h2) / 2) / H)


def test_points():
    import pandas as pd
    from scipy.signal import butter, filtfilt
    pts = []
    for test, name, nmax in (("Test9", "US-1", 24), ("Test12", "US-2", 25)):
        f = fractions(heights(test))
        for n in range(1, nmax + 1):
            p = os.path.join(RAW, "{}Run{}.xlsx".format(test, n))
            if not os.path.isfile(p):
                continue
            df = pd.read_excel(p); df.columns = [str(c).strip() for c in df.columns]
            t = df["Time"].values.astype(float); fs = 1 / np.median(np.diff(t))
            b, a = butter(4, LP_HZ / (0.5 * fs))
            L = lambda c: filtfilt(b, a, df[c].values.astype(float) - df[c].values[:200].astype(float).mean())
            wall = WALL_KN * (f[0] * L("15") + f[1] * L("16") + f[2] * L("17"))
            slab = G / 1e3 * SLAB_KG * (L("13") + L("14"))
            d = 0.5 * (df["3"].values + df["4"].values).astype(float) * 1e3
            d = filtfilt(b, a, d - d[0])
            s = np.sign(np.polyfit(d, wall + slab, 1)[0]); wall *= s; slab *= s
            w = int(0.05 * fs)
            for sg, i in ((+1, int(np.argmax(d))), (-1, int(np.argmin(d)))):
                sl = slice(max(0, i - w), i + w)
                k = sl.start + int(np.argmax(sg * wall[sl]))
                pts.append(dict(test=name, run=n, sign=sg, d_mm=round(abs(d[i]), 3), F_wall_kN=round(sg * wall[k], 3),
                                F_total_kN=round(sg * (wall[k] + slab[k]), 3), slab_share=round(slab[k] / (wall[k] + slab[k]), 3)))
    return pts


def read_push(path):
    lines = open(path).read().splitlines()
    mt = float(re.search(r"m_tot\s+([0-9.]+)", lines[0]).group(1))
    rows = list(csv.DictReader(l for l in lines[1:] if l.strip()))
    return np.abs([float(r["d_ctrl_mm"]) for r in rows]), np.abs([float(r["a_mps2"]) for r in rows]) * mt / 1e3


def mech(side, inertia, K, top=None, h1=None):
    """rigid three-hinge curve (d mm, F kN on the wall mass). h1 None -> the weakest hinge height at every d (the
    mechanism a rigid wall chooses); h1 given -> the hinge held there. Elastic branch K from the origin."""
    W = M_MASONRY * 9.81 / 1e3
    js = [rb.Joist(0.3225, 0.186 * 2, (0.996 - 0.186) * 2, m_floor=FLOOR_KG, side=side),
          rb.Joist(2.1989, 0.186 * 2, (0.996 - 0.186) * 2, m_floor=FLOOR_KG, side=side)]
    dg = np.linspace(0.2, 45.0, 300)
    hs = [h1] if h1 is not None else list(np.arange(1.0, 2.451, 0.05))
    Fs = []
    for hh in hs:
        try:
            d, F = rb.curve(hh, W, O_KN, joists=js, inertia=inertia, top=top, thmax=0.08, n=260)
        except Exception:
            continue
        if d[-1] >= dg[-1]:
            Fs.append(np.interp(dg, d, F))
    F = np.min(np.array(Fs), axis=0)
    d, F, _ = rb.with_elastic_branch(dg, F, K)
    return d, F


MECH = [  # (label, inertia, top, h1, style)
    ("rigid, floors NOT pushed, top-beam corner, weakest hinge", False, None, None, dict(color="#374151", ls="--", lw=1.6)),
    ("rigid, floors NOT pushed, top-beam corner, hinge at upper pocket 2.11 m", False, None, 2.11, dict(color="#374151", ls="-.", lw=1.2)),
    ("rigid, floors pushed, top-beam corner, weakest hinge", True, None, None, dict(color="#7C3AED", ls="--", lw=1.6)),
    ("rigid, floors pushed, top at beam centroid, weakest hinge", True, (0.105, 2.63), None, dict(color="#7C3AED", ls=":", lw=1.8)),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    pts = test_points()
    with open(os.path.join(OUT, "test_envelope_wall_only.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pts[0])); w.writeheader(); w.writerows(pts)
    push = {(lab, s): read_push(os.path.join(fo, "pushover_{}.csv".format(s))) for lab, fo, _ in CASES for s in ("pos", "neg")}
    dk, Fk = push[("DEM, pK11 pocket", "pos")]
    K = float(np.polyfit(dk[dk <= 0.6], Fk[dk <= 0.6], 1)[0])
    mrows = []
    fig, ax = plt.subplots(2, 2, figsize=(16, 11))
    for j, (side_lab, s, sg, side) in enumerate((("Positive (towards the joists)", "pos", +1, "leeward"),
                                                ("Negative (away from the joists)", "neg", -1, "windward"))):
        for k, xmax in enumerate((40.0, 5.0)):
            a = ax[k, j]
            for lab, fo, col in CASES:
                d, F = push[(lab, s)]
                a.plot(d, F, color=col, lw=2.0, label=lab + " pushover, |a| m_tot")
            for ml, inertia, top, h1, sty in MECH:
                d, F = mech(side, inertia, K, top, h1)
                a.plot(d, F, label=ml, **sty)
                if k == 0:
                    for x in (5, 10, 20, 30, 40):
                        mrows.append(dict(side=s, model=ml, d_mm=x, F_kN=round(float(np.interp(x, d, F)), 2)))
            for tn, mk, mc in (("US-1", "o", "#111827"), ("US-2", "^", "#1D4ED8")):
                P = [p for p in pts if p["test"] == tn and p["sign"] == sg and p["d_mm"] <= xmax]
                a.plot([p["d_mm"] for p in P], [p["F_wall_kN"] for p in P], mk, mfc=mc, mec=mc, ms=6, ls="none",
                       label=tn + ": wall-only base shear at the run's peak d")
                a.plot([p["d_mm"] for p in P], [p["F_total_kN"] for p in P], mk, mfc="none", mec=mc, ms=5, ls="none",
                       alpha=0.35, label=tn + ": full formula (with slabs)")
            a.set_xlim(0, xmax); a.set_ylim(0, 16); a.grid(alpha=0.25)
            a.set_xlabel("|d| at 2.06 m (mm)"); a.set_ylabel("|F| (kN)")
            a.set_title(side_lab + ("" if k == 0 else ", first 5 mm"), loc="left")
            if j == 0 and k == 0:
                a.legend(fontsize=7.5, loc="lower right", ncol=1)
    fig.suptitle("Pushovers (pk2 f_t 0.05, glued vs pK11 pocket) against the tests' wall-only base shear and the rigid "
                 "three-hinge model; elastic branch of the rigid model K = {:.2f} kN/mm (pK11 pushover)".format(K), fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_pushover_test_mech.png"), dpi=130)
    with open(os.path.join(OUT, "mech_curves.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mrows[0])); w.writeheader(); w.writerows(mrows)
    print("K = {:.2f} kN/mm".format(K))
    for r in mrows:
        if r["d_mm"] in (10, 20, 30):
            print(r["side"], r["d_mm"], r["F_kN"], r["model"])
    import statistics as st
    for tn in ("US-1", "US-2"):
        for sg in (+1, -1):
            P = [p for p in pts if p["test"] == tn and p["sign"] == sg]
            big = [p for p in P if p["d_mm"] >= 5]
            print(tn, "+" if sg > 0 else "-", "slab share median {:.2f}".format(st.median(p["slab_share"] for p in P)),
                  "| d >= 5 mm: F_wall", [round(p["F_wall_kN"], 1) for p in big], "F_total", [round(p["F_total_kN"], 1) for p in big],
                  "at d", [round(p["d_mm"], 1) for p in big])
    print("-> " + OUT)


if __name__ == "__main__":
    main()
