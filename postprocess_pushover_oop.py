# -*- coding: ascii -*-
"""
postprocess_pushover_oop.py -- OOP pushover curves of the US-1 model after batch_ft0p15_pipeline.py.

WHAT IT READS (run from the 3DEC working folder; every path is relative to BASE_DIR)
    CASES        the pushover output folders of the batch (pushover_pos.csv / pushover_neg.csv in each):
                   f_t = 0,    mason top   pushover_results_v5_ft0_masonrypush_TopMason     (batch stages 0, 1)
                   f_t = 0.15, mason top   pushover_results_v5_ft0p15_masonrypush_TopMason  (batch stages 2, 3)
                   f_t = 0,    elastic top pushover_results_v5_ft0_masonrypush              (22 Oct reference)
                   forced hinges          pushover_results_v5_hingeft0_char_TopMason       (run_pushover_hinge_ft0.py)
                 A case or a branch that is missing is reported and skipped. pushover_*_prev_*.csv (kept by the
                 batch when it re-ran a branch) are ignored.
    BILINEAR_JSON  bilinear_idealisation_ft0p15/bilinear_periods.json (batch stage 4), drawn on the f_t = 0.15 curves
    BATCH_DIR      batch_ft0p15/batch_state.json: which stages finished (printed)
    ENVELOPE_CSV   us1_fig13_digitised.csv (test envelope, Fig 13), divided by FIG13_MASS_SCALE = 2.317, the
                   mass scale the processed test data remove (Fig 13 as printed peaks at ~30 kN, i.e. 2.3 x sum m a)
    TEST_PEAKS     exp_envelope_peaks/envelope_peaks_US1/US2.csv: peak force per run (sum m a) at its displacement

FORCE
    The v5 driver logs F_tot = base shear + top-joint SHEAR z-components. Once the top beam tilts, the top-joint
    normal force has a horizontal component |N| sin(tilt) that F_tot leaves out (deformed_shape_3dec.py, 23 Oct:
    +1.8 / +2.6 kN at 80 mm on the elastic-top runs). The force the wall carries is the applied push |a| x m_tot
    at a settled servo point (m_tot from the CSV's '#' line). FORCE = "applied" (default) plots and measures that;
    F_tot as logged can be drawn dotted for reference (SHOW_REACTION, off by default). Rows with converged = 0 (the
    runaway step that ends the force-controlled ascent) are not equilibria: left out of every metric and every curve
    (SHOW_UNCONVERGED = True draws them as hollow markers).

METRICS per branch (pushover_metrics.csv)
    peak F and d there; K_el = least-squares slope through the origin over K_EL_RANGE_MM -> T1 = 2 pi sqrt(L / K_el);
    T_eff = 2 pi sqrt(L d_u / F_y) with d_u = D_U_MM and F_y = peak (the Route 2 driver's law, L = M_EFF);
    secant at D_SEC_MM; force at the last point; largest joint / pocket opening; joist support load range;
    F_tot - |a| m_tot at the last point.
    Per case with both branches, the Route 2 "mean" combination (K_el pooled over both branches, F_y = mean of
    the two peaks): the T1 / T_eff route2_Tm_ihsan_full25_3dec.py computes from the same files (periods row).

MECHANICAL MODEL (MECH = True; needs mechanical_model_US1/rigid_two_block_vw.py, MECHANICAL_MODEL_US1.md Appendix D)
    One curve: the rigid two-block strip with the ADJUSTED HINGE at MECH_H1 = 2.09 m, the joint the DEM pushover
    (and the test) opens. Large-displacement virtual work with the DEM's loads: masonry m_tot of the CSV, springs
    10.35 kN + the 213 kg top beam, joists pinned with their bearing load in the pockets; top pivot at the top-beam
    centroid (the beam turns with the upper block, as in the DEM). Lower pivot: the line has it at the base toe,
    the shaded band runs to a lower pivot at 0.30 m (the kink at the lower floor in the DEM deformed shape).
    CRACKING FORCE (MECH_FCR): the bonded wall must first crack. Graziotti et al. (2016) / Tomassetti et al. (2018)
    strip: the lateral load q at which the bending moment at height y reaches the cracking moment
        M_c(y) = (f_w + N(y)/t) t^2/6,  N(y) = O + W (h - y)/h,
    with the base and top joints at their own cracking moments M_b, M_t:
        q(y) = 2 [M_c(y) + M_b (1 - y/h) + M_t y/h] / (y (h - y)),   F_cr(y) = q(y) h b.
    Onur's formula is the minimum over y (h1 = 1.36 m for f_w = 0.15); here y = MECH_H1, the crack forced at the
    adjusted hinge. Same loads as the rocking curve. Both values are printed and written to pushover_mech.csv.
    Values at 20 / 40 / 80 mm and the zero-force displacement -> pushover_mech.csv, next to the DEM at the same d.

OUTPUT (OUT_DIR)
    fig_pushover_curves.png       (a, b) both directions, all cases, test envelope, bilinear idealisation;
                                  (c, d) the first ZOOM_MM with the K_el fits
    fig_pushover_diagnostics.png  applied acceleration, joint / pocket opening, joist support load vs d
    pushover_metrics.csv, pushover_periods.csv
    and, printed, the snapshot saves found in each folder as a SAVES list for deformed_shape_3dec.py
Usage: python postprocess_pushover_oop.py [BASE_DIR]
"""
import os, sys, csv, io, glob, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
CASES = [   # (key, label, folder, colour, line style)
    ("ft0_mason",   "f_t = 0, mason top",          "pushover_results_v5_ft0_masonrypush_TopMason",    "#2B5F8A", "-"),
    ("ft015_mason", "f_t = 0.15 MPa, mason top",   "pushover_results_v5_ft0p15_masonrypush_TopMason", "#E07A2D", "-"),
    ("ft0_elastic", "f_t = 0, elastic top (ref.)", "pushover_results_v5_ft0_masonrypush",             "#9CA3AF", "--"),
    # forced mechanism (run_pushover_hinge_ft0.py / ANALYSIS_PART_I_MASON_hingeft0.dat): f_t = 0 only at the base,
    # the 2.110 m bed joint, the top joint and the joist pockets; 0.344 MPa / c 0.75 MPa everywhere else.
    # Compare with the solid edge of the mechanical-model band (lower pivot at the base; 2.09 m is within 0.2 kN of 2.11 m)
    ("hinge_char",  "forced hinges (f_t = 0 at base, 2.11 m, top)", "pushover_results_v5_hingeft0_char_TopMason", "#7C3AED", "-."),
]
BILINEAR_JSON = "bilinear_idealisation_ft0p15/bilinear_periods.json"
BILINEAR_CASE = "ft015_mason"            # the case the bilinear idealisation was made from
BATCH_DIR     = "batch_ft0p15"
ENVELOPE_CSV  = None                     # was "us1_fig13_digitised.csv": the paper's Fig 13 is on the OLD basis (2.317
                                         # factor, table acceleration added to absolute channels); off since 27 Oct
FIG13_MASS_SCALE = 2.317                 # Fig 13's base shear carries the paper's mass scale; the processed test data
                                         # divide it out (envelope_peaks_*.csv: mass_scale_removed 2.317). Same here,
                                         # so the envelope is on the model's basis (sum m a on the wall). 1.0 = as printed.
TEST_PEAKS    = {"US-1": "exp_envelope_peaks_v3/envelope_peaks_US1.csv",  # peak force per run and the displacement
                 "US-2": "exp_envelope_peaks_v3/envelope_peaks_US2.csv"}  # there; {} = off. v3 = the group's base shear:
                                         # python exp_envelope_peaks.py exp_processed_v3 exp_envelope_peaks_v3
                                         # (the old exp_envelope_peaks/ files are about 35 % high)
FORCE         = "applied"                # "applied": |a| x m_tot | "reaction": F_tot_kN as logged
SHOW_REACTION = False                    # True: also draw F_tot as logged, dotted
SHOW_UNCONVERGED = False                 # True: draw the unconverged (runaway) steps as hollow markers
M_TOT_DEFAULT = 1308.3                   # kg, if the CSV has no '# ... m_tot ...' line
M_EFF         = 1200.0                   # kg, L = sum m_j phi_j (as in Route 2 and the batch's bilinear stage)
K_EL_RANGE_MM = (0.1, 1.5)               # K_el fit range (Route 2)
D_U_MM        = 27.4                     # d_u of T_eff (Route 2: end of the US-1 envelope)
D_SEC_MM      = 29.5                     # secant evaluated here (US-1 run 24)
ZOOM_MM       = 5.0
POCKET_EXIT_MM = 50.0                    # the v5 driver aborts beyond this
OUT_DIR       = "pushover_postproc_batch"

# ---- mechanical model (rigid two-block strip) ---------------------------------------------------------------
MECH            = True
MECH_DIR        = "mechanical_model_US1"     # relative to BASE_DIR, or a full path
MECH_H1         = 2.09                       # m, adjusted hinge: the joint the DEM pushover opens (upper-floor course)
MECH_TB_Y       = 2.629                      # m, top-beam centroid in the build (top pivot)
MECH_TB_KG      = 212.986                    # kg, top beam in the build
MECH_SPRING_KN  = 10.35
MECH_FCR        = True                       # cracking force of the bonded wall with its crack at MECH_H1 (drawn on the
MECH_FW         = 0.15                       #   f_t = 0.15 case): bond strength f_w in MPa
MECH_B_M        = 1.292                      # m, wall width (per-unit-width strip formula -> kN)
MECH_BASE_Y     = (0.0, 0.30)                # m, lower pivot: the line (base toe) and the far edge of the band
                                             #   (kink at the lower floor, DEM deformed shape)
MECH_COL        = "#0F766E"

INK, MUTED, GRID, TEST = "#1f2937", "#6b7280", "#e5e7eb", "#374151"
P = lambda *p: os.path.join(BASE_DIR, *p)
os.makedirs(P(OUT_DIR), exist_ok=True)


# ============================ READERS ================================
def read_v5(path):
    """v5 pushover CSV -> dict of arrays (|d|, applied and reaction force, diagnostics) or None"""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        lines = f.readlines()
    head = " ".join(l for l in lines if l.lstrip().startswith("#"))
    rows = list(csv.DictReader(io.StringIO("".join(l for l in lines if l.strip() and not l.lstrip().startswith("#")))))
    if not rows:
        return None
    m_tot = M_TOT_DEFAULT
    if "m_tot" in head:
        try:
            m_tot = float(head.split("m_tot")[1].split("kg")[0])
        except ValueError:
            pass
    meta = {}
    for key in ("save", "top joint"):
        if key in head:
            meta[key] = head.split(key, 1)[1].split(";")[0].strip()

    def col(name, default=np.nan):
        out = []
        for r in rows:
            try:
                out.append(float(r.get(name, "")))
            except (TypeError, ValueError):
                out.append(default)
        return np.array(out)
    a = np.abs(col("a_mps2"))
    c = dict(n=len(rows), m_tot=m_tot, meta=meta, path=path,
             phase=np.array([r.get("phase", "") for r in rows]),
             d=np.abs(col("d_ctrl_mm")), a=a,
             F_app=a * m_tot / 1e3, F_rea=np.abs(col("F_tot_kN")),
             ok=col("converged", 1.0) >= 0.5,
             open=col("max_open_mm"), pocket=col("pocket_open_mm"), joistN=col("joist_N_kN"))
    c["F"] = c["F_app"] if FORCE == "applied" else c["F_rea"]
    return c


def read_envelope(path):
    if not os.path.isfile(path):
        print("  ! {} not found -- no test envelope".format(path)); return {}
    e = np.genfromtxt(path, delimiter=",", names=True)
    out = {}
    for br, m in (("pos", e["d_mm"] > 0), ("neg", e["d_mm"] < 0)):
        d = np.abs(e["d_mm"][m]); F = np.abs(e["F_kN"][m]) / FIG13_MASS_SCALE; o = np.argsort(d)
        out[br] = (d[o], F[o])
    return out


def read_peaks(path):
    """peak force per run and the displacement there (relative to the run's start), per sign"""
    if not os.path.isfile(path):
        return None
    out = {"pos": [], "neg": []}
    for r in csv.DictReader(open(path)):
        try:
            d0 = float(r["d_start_mm"])
            out["pos"].append((abs(float(r["d_at_Fmax_pos_mm"]) - d0), abs(float(r["F_max_pos_kN"]))))
            out["neg"].append((abs(float(r["d_at_Fmax_neg_mm"]) - d0), abs(float(r["F_max_neg_kN"]))))
        except (KeyError, ValueError):
            continue
    return {k: np.array(v) for k, v in out.items()}


def read_bilinear(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        j = json.load(f)
    return j.get("model")


# ============================ METRICS ================================
def k_el(d, F):
    m = (d >= K_EL_RANGE_MM[0]) & (d <= K_EL_RANGE_MM[1])
    if m.sum() < 2:
        m = np.argsort(d)[:3]
    return float(np.sum(d[m] * F[m]) / np.sum(d[m] ** 2)), d[m], F[m]


def periods(K, Fy):
    T1 = 2 * np.pi * np.sqrt(M_EFF / (K * 1e6)) if K > 0 else np.nan
    Te = 2 * np.pi * np.sqrt(M_EFF * D_U_MM * 1e-3 / (Fy * 1e3)) if Fy > 0 else np.nan
    return T1, Te


def branch_metrics(c):
    ok = c["ok"]; d, F = c["d"][ok], c["F"][ok]
    o = np.argsort(d, kind="stable"); d, F = d[o], F[o]
    i = int(np.argmax(F)); K, _, _ = k_el(d, F); T1, Te = periods(K, F[i])
    Fs = float(np.interp(D_SEC_MM, d, F)) if d.max() >= D_SEC_MM else np.nan
    end = int(np.argmax(c["d"]))
    return dict(peak_F_kN=float(F[i]), d_at_peak_mm=float(d[i]), K_el_kNmm=K, T1_s=T1, Teff_s=Te,
                F_sec_kN=Fs, K_sec_kNmm=Fs / D_SEC_MM if np.isfinite(Fs) else np.nan,
                d_end_mm=float(c["d"][end]), F_end_kN=float(c["F"][end]),
                Ftot_minus_applied_end_kN=float(c["F_rea"][end] - c["F_app"][end]),
                max_joint_open_mm=float(np.nanmax(c["open"])) if np.isfinite(c["open"]).any() else np.nan,
                max_pocket_open_mm=float(np.nanmax(c["pocket"])) if np.isfinite(c["pocket"]).any() else np.nan,
                joistN_min_kN=float(np.nanmin(c["joistN"])) if np.isfinite(c["joistN"]).any() else np.nan,
                joistN_max_kN=float(np.nanmax(c["joistN"])) if np.isfinite(c["joistN"]).any() else np.nan,
                n_rows=c["n"], n_unconverged=int((~c["ok"]).sum()), m_tot_kg=c["m_tot"],
                top_joint=c["meta"].get("top joint", ""), save=c["meta"].get("save", ""))


def combined_periods(cp, cn):
    """Route 2 'mean': K_el pooled over both branches, F_y = mean of the two peaks"""
    dd, ff = [], []
    for c in (cp, cn):
        ok = c["ok"]; _, d, F = k_el(c["d"][ok], c["F"][ok]); dd.append(d); ff.append(F)
    dd, ff = np.concatenate(dd), np.concatenate(ff)
    K = float(np.sum(dd * ff) / np.sum(dd * dd))
    Fy = 0.5 * (cp["F"][cp["ok"]].max() + cn["F"][cn["ok"]].max())
    T1, Te = periods(K, Fy)
    return dict(K_el_kNmm=K, F_y_kN=Fy, T1_s=T1, Teff_s=Te)


# ============================ COLLECT ================================
print("=" * 78)
print("pushover post-processing  [{}]  force: {}".format(os.path.abspath(BASE_DIR), "|a| x m_tot" if FORCE == "applied" else "F_tot as logged"))
st = P(BATCH_DIR, "batch_state.json")
if os.path.isfile(st):
    with open(st) as f:
        done = json.load(f).get("done", [])
    print("batch stages finished: {}".format(sorted(done) or "none"))
else:
    print("batch state: {} not found".format(st))

data = {}            # (case key, branch) -> curve
snaps = {}
for key, label, folder, colr, ls in CASES:
    for br in ("pos", "neg"):
        c = read_v5(P(folder, "pushover_{}.csv".format(br)))
        tag = "{:28s} {}".format(label, br)
        if c is None:
            print("  {}  NOT FOUND  [{}]".format(tag, os.path.join(folder, "pushover_{}.csv".format(br)))); continue
        data[(key, br)] = c
        top = c["meta"].get("top joint", "?")
        print("  {}  {} rows ({} unconverged), d to {:.1f} mm, peak {:.2f} kN, m_tot {:.1f} kg, top joint {}".format(
            tag, c["n"], int((~c["ok"]).sum()), c["d"].max(), c["F"][c["ok"]].max(), c["m_tot"], top))
        if ("elastic" in key) != (top == "elastic") and top != "?":
            print("    ! the CSV says top joint '{}' -- check that the folder holds the case its label says".format(top))
        if np.nanmax(c["pocket"]) > POCKET_EXIT_MM:
            print("    ! pocket opening {:.0f} mm > {:.0f} mm: a joist left its pocket".format(np.nanmax(c["pocket"]), POCKET_EXIT_MM))
    found = sorted(glob.glob(P(folder, "pushover_*_d[0-9]*.sav")) + glob.glob(P(folder, "pushover_*_end.sav")))
    if found:
        snaps[key] = (folder, [os.path.splitext(os.path.basename(s))[0] for s in found])
env = read_envelope(P(ENVELOPE_CSV)) if ENVELOPE_CSV else {}
peaks = {k: read_peaks(P(v)) for k, v in TEST_PEAKS.items()}
for k, v in peaks.items():
    print("test peaks {}: {}".format(k, "{} runs".format(len(v["pos"])) if v is not None else "not found ({})".format(TEST_PEAKS[k])))
bil = read_bilinear(P(BILINEAR_JSON))
print("bilinear idealisation: {}".format(
    "{} (K1 {:.2f} kN/mm, mu {:.2f}, T1 {:.4f} s, T_eff {:.4f} s, m_eff {:.0f})".format(
        BILINEAR_JSON, bil["K1_kNmm"], bil["mu"], bil["T1_s"], bil["Teff_s"],
        list(bil["per_branch"].values())[0].get("m_eff", np.nan)) if bil else "not found"))
if not data:
    raise SystemExit("no pushover CSV found in any case folder under " + os.path.abspath(BASE_DIR))

# ---- metrics
rows, prow = [], []
for key, label, folder, colr, ls in CASES:
    for br in ("pos", "neg"):
        if (key, br) in data:
            rows.append(dict(case=label, branch=br, folder=folder, **branch_metrics(data[(key, br)])))
    if (key, "pos") in data and (key, "neg") in data:
        prow.append(dict(case=label, combine="mean (Route 2)", **combined_periods(data[(key, "pos")], data[(key, "neg")])))
for br in ("pos", "neg"):
    if br in env:
        d, F = env[br]; i = int(np.argmax(F)); K, _, _ = k_el(d, F); T1, Te = periods(K, F[i])
        rows.append(dict(case="US-1 test (Fig 13)", branch=br, folder=ENVELOPE_CSV, peak_F_kN=float(F[i]), d_at_peak_mm=float(d[i]),
                         K_el_kNmm=K, T1_s=T1, Teff_s=Te))


def write(path, rr):
    keys = []
    for r in rr:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rr:
            w.writerow(["" if r.get(k) is None else (round(r[k], 4) if isinstance(r.get(k), float) and np.isfinite(r[k]) else
                        ("" if isinstance(r.get(k), float) else r.get(k, ""))) for k in keys])


write(P(OUT_DIR, "pushover_metrics.csv"), rows)

# ---- mechanical model: adjusted hinge
mech = {}
MECH_LABEL = "mechanical model, adjusted hinge {:.2f} m".format(MECH_H1)
if MECH:
    try:
        sys.path.insert(0, os.path.abspath(P(MECH_DIR)))
        import rigid_two_block_vw as vw
        m_ref = next((c["m_tot"] for c in data.values()), M_TOT_DEFAULT)
        W, O = m_ref * 9.81 / 1e3, MECH_SPRING_KN + MECH_TB_KG * 9.81 / 1e3
        for br, side in (("pos", "leeward"), ("neg", "windward")):
            for j, yb in enumerate(MECH_BASE_Y):
                for thmax in (0.25, 0.2, 0.16, 0.12, 0.08):     # largest rotation the geometry allows
                    try:
                        dd, FF = vw.curve(MECH_H1, W, O, joists=vw.us1_joists(side, k_r=0.0), thmax=thmax, n=600,
                                          top=(vw.t / 2, MECH_TB_Y), base_y=yb)
                        break
                    except ValueError:
                        continue
                k = int(np.argmax(FF < 0)) if (FF < 0).any() else len(FF)
                dz = float(np.interp(0.0, [-FF[k - 1], -FF[k]], [dd[k - 1], dd[k]])) if k < len(FF) else np.nan
                mech[(br, j)] = (np.r_[dd[:k], dz] if k < len(FF) else dd, np.r_[FF[:k], 0.0] if k < len(FF) else FF, dz)
        print("mechanical model: adjusted hinge {:.2f} m, W {:.2f} kN, O {:.2f} kN at the top-beam centroid, joists pinned".format(MECH_H1, W, O))
    except Exception as e:
        print("  ! mechanical model skipped: {} ({})".format(e, P(MECH_DIR, "rigid_two_block_vw.py")))
        mech = {}
def cracking_force(fw, y, W_kN, O_kN, t=0.21, h=2.58, b=MECH_B_M):
    """F_cr (kN) with the crack at height y (m); fw in MPa; returns (F_cr at y, min F_cr, its height)"""
    tm, hm = t * 1e3, h * 1e3; Wl, Ol = W_kN * 1e3 / (b * 1e3), O_kN * 1e3 / (b * 1e3)   # N/mm
    Mb = (fw + (Wl + Ol) / tm) * tm ** 2 / 6; Mt = (fw + Ol / tm) * tm ** 2 / 6

    def F(yy):
        yy = np.asarray(yy, float) * 1e3
        Mc = (fw + (Ol + Wl * (hm - yy) / hm) / tm) * tm ** 2 / 6
        q = 2 * (Mc + Mb * (1 - yy / hm) + Mt * yy / hm) / (yy * (hm - yy))
        return q * hm * (b * 1e3) / 1e3              # (N/mm per mm of width) x h (mm) x b (mm) -> N -> kN
    ys = np.linspace(0.05, h - 0.05, 2000); Fs = F(ys); i = int(np.argmin(Fs))
    return float(F(y)), float(Fs[i]), float(ys[i])


FCR = None
if mech and MECH_FCR:
    FCR = cracking_force(MECH_FW, MECH_H1, W, O)
    print("cracking force, f_w {:.2f} MPa: {:.2f} kN with the crack at {:.2f} m (adjusted hinge); "
          "minimum {:.2f} kN at {:.2f} m (intact-wall optimum)".format(MECH_FW, FCR[0], MECH_H1, FCR[1], FCR[2]))
if mech:
    mrows = []
    if FCR:
        mrows.append(dict(model="cracking force f_w {:.2f} MPa, crack at {:.2f} m".format(MECH_FW, MECH_H1), branch="both", F0_kN=FCR[0]))
        mrows.append(dict(model="cracking force f_w {:.2f} MPa, minimum (crack at {:.2f} m)".format(MECH_FW, FCR[2]), branch="both", F0_kN=FCR[1]))
    for br in ("pos", "neg"):
        for j, yb in enumerate(MECH_BASE_Y):
            dd, FF, dz = mech[(br, j)]
            r = dict(model="{}, lower pivot {:.2f} m".format(MECH_LABEL, yb), branch=br, F0_kN=float(FF[0]), zero_force_at_mm=dz)
            for x in (20, 40, 80):
                r["F_{}mm_kN".format(x)] = float(np.interp(x, dd, FF)) if x <= dd[-1] else 0.0
            mrows.append(r)
        for ckey, label, folder, colr, ls in CASES:      # the DEM at the same displacements
            c = data.get((ckey, br))
            if c is None:
                continue
            ok = c["ok"]; o = np.argsort(c["d"][ok]); d_, F_ = c["d"][ok][o], c["F"][ok][o]
            r = dict(model="DEM: " + label, branch=br, F0_kN=float(F_.max()))
            for x in (20, 40, 80):
                r["F_{}mm_kN".format(x)] = float(np.interp(x, d_, F_)) if x <= d_.max() else np.nan
            mrows.append(r)
    write(P(OUT_DIR, "pushover_mech.csv"), mrows)
    print("\n{:60s} {:4s} {:>6s} {:>6s} {:>6s} {:>6s} {:>8s}".format("mechanical model vs DEM (F0 / peak, F at d)", "br", "F0", "20mm", "40mm", "80mm", "zero at"))
    for r in mrows:
        print("{:60s} {:4s} {:6.2f} {:>6s} {:>6s} {:>6s} {:>8s}".format(
            r["model"][:60], r["branch"], r["F0_kN"],
            *["{:.2f}".format(r["F_{}mm_kN".format(x)]) if np.isfinite(r.get("F_{}mm_kN".format(x), np.nan)) else "-" for x in (20, 40, 80)],
            "{:.0f} mm".format(r["zero_force_at_mm"]) if np.isfinite(r.get("zero_force_at_mm", np.nan)) else ""))
    print("-> " + P(OUT_DIR, "pushover_mech.csv"))
if prow:
    write(P(OUT_DIR, "pushover_periods.csv"), prow)
print("\n{:28s} {:4s} {:>7s} {:>6s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s}".format(
    "case", "br", "peak kN", "at mm", "K_el", "T1 s", "T_eff s", "d_end", "F_end", "pocket"))
for r in rows:
    print("{:28s} {:4s} {:7.2f} {:6.1f} {:7.2f} {:7.4f} {:7.4f} {:>7s} {:>7s} {:>7s}".format(
        r["case"][:28], r["branch"], r["peak_F_kN"], r["d_at_peak_mm"], r["K_el_kNmm"], r["T1_s"], r["Teff_s"],
        "{:.1f}".format(r["d_end_mm"]) if "d_end_mm" in r else "", "{:.2f}".format(r["F_end_kN"]) if "F_end_kN" in r else "",
        "{:.1f}".format(r["max_pocket_open_mm"]) if "max_pocket_open_mm" in r else ""))
for r in prow:
    print("{:28s} mean: K_el {:.2f} kN/mm, F_y {:.2f} kN -> T1 {:.4f} s, T_eff {:.4f} s  (L {:.0f} kg, d_u {:.1f} mm)".format(
        r["case"][:28], r["K_el_kNmm"], r["F_y_kN"], r["T1_s"], r["Teff_s"], M_EFF, D_U_MM))
print("-> " + P(OUT_DIR, "pushover_metrics.csv") + ("  " + P(OUT_DIR, "pushover_periods.csv") if prow else ""))

if snaps:
    print("\nsnapshot saves (paste into deformed_shape_3dec.py, with REF_SAV = the save of that case):")
    for key, (folder, names) in snaps.items():
        print("  # {}".format(key))
        for n in names:
            print('    ("{}", "{}/{}", {}),'.format(n.replace("pushover_", ""), folder, n, "+1" if "_pos_" in n else "-1"))


# ============================ FIG 1: curves =========================
def style_axes(ax):
    ax.grid(color=GRID); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


fig, axs = plt.subplots(2, 2, figsize=(14, 9.8), gridspec_kw=dict(height_ratios=[1.6, 1]))
for j, (br, title) in enumerate((("pos", "Positive push (towards the joists)"), ("neg", "Negative push"))):
    for row, zoom in ((0, False), (1, True)):
        ax = axs[row, j]
        if br in env:
            d, F = env[br]
            ax.plot(d, F, "-", color=TEST, lw=1.6, label="US-1 envelope (Fig 13 / {:g}), same sign".format(FIG13_MASS_SCALE))
            if not zoom:
                i = int(np.argmax(F)); ax.plot(d[i], F[i], "o", color=TEST, ms=5)
                ax.annotate("{:.1f} kN".format(F[i]), (d[i], F[i]), xytext=(4, 5), textcoords="offset points", fontsize=8, color=TEST)
        if mech and not zoom:
            dm, fm, zm = mech[(br, 0)]
            ax.plot(dm, fm, color=MECH_COL, lw=2.0, zorder=1.6, label=MECH_LABEL + " (lower pivot at the base)")
            if len(MECH_BASE_Y) > 1:
                d2, f2, z2 = mech[(br, len(MECH_BASE_Y) - 1)]
                xs = np.linspace(0.3, max(zm, z2) if np.isfinite(z2) else dm[-1], 300)
                ax.fill_between(xs, np.interp(xs, dm, fm, right=np.nan), np.interp(xs, d2, f2, right=np.nan),
                                color=MECH_COL, alpha=0.14, lw=0, zorder=1,
                                label="  band: lower pivot from the base to {:.2f} m".format(MECH_BASE_Y[-1]))
                ax.plot(d2, f2, color=MECH_COL, lw=0.9, ls="--", zorder=1.5)
            if FCR and ("ft015_mason", br) in data:
                ax.plot([0, 18], [FCR[0], FCR[0]], ":", color=MECH_COL, lw=1.8, zorder=1.6,
                        label="mechanical model: cracking force, crack at {:.2f} m, f_w = {:.2f} MPa".format(MECH_H1, MECH_FW))
                ax.text(18.5, FCR[0], "{:.1f} kN".format(FCR[0]), color=MECH_COL, fontsize=8.3, va="center")
            if np.isfinite(zm):
                ax.annotate("zero at {:.0f} mm".format(zm), xy=(zm, 0), xytext=(zm - 22, 3.2), fontsize=8.3, color=MECH_COL, ha="right",
                            arrowprops=dict(arrowstyle="-", color=MECH_COL, lw=0.8))
        for (tk, pk), mk in zip(peaks.items(), ("o", "^")):
            if pk is not None and len(pk[br]):
                ax.plot(pk[br][:, 0], pk[br][:, 1], mk, ms=3.6, mfc="none", mec=TEST, mew=0.8, alpha=0.7,
                        label="{} test, peak force per run".format(tk))
        for key, label, folder, colr, ls in CASES:
            c = data.get((key, br))
            if c is None:
                continue
            ok = c["ok"]
            ax.plot(c["d"][ok], c["F"][ok], ls, color=colr, lw=2.0 if ls == "-" else 1.4, label=label, zorder=3 if ls == "-" else 2)
            if SHOW_UNCONVERGED and (~ok).any():
                ax.plot(c["d"][~ok], c["F"][~ok], "o", mfc="white", mec=colr, ms=5, label=None)
            if SHOW_REACTION and FORCE == "applied" and not zoom:
                ax.plot(c["d"][ok], c["F_rea"][ok], ":", color=colr, lw=1.0, label=None)
            if zoom:
                K, de, Fe = k_el(c["d"][ok], c["F"][ok])
                xs = np.linspace(0, ZOOM_MM, 2)
                ax.plot(xs, K * xs, ls=(0, (2, 2)), color=colr, lw=1.0)
        if bil and BILINEAR_CASE and (BILINEAR_CASE, br) in data and br in bil.get("per_branch", {}):
            b = bil["per_branch"][br]; colr = dict((k, cc) for k, _, _, cc, _ in CASES)[BILINEAR_CASE]
            ax.plot([0, b["dy_mm"], b["du_mm"]], [0, b["Vu_kN"], b["Vu_kN"]], "-", color=colr, lw=1.0, alpha=0.6,
                    label="bilinear idealisation (f_t = 0.15)")
            if not zoom:
                ax.plot(b["du_mm"], 0.8 * b["Vu_kN"], "x", color=colr, ms=7)
        if zoom:
            ax.set_xlim(0, ZOOM_MM); ax.set_ylim(0, 16)
            ax.set_title("{}: first {:g} mm (dashed: K_el fit over {:g}-{:g} mm)".format(
                "(c)" if j == 0 else "(d)", ZOOM_MM, *K_EL_RANGE_MM), loc="left", fontsize=10)
        else:
            ax.axvline(D_U_MM, color="#cbd5e1", lw=0.9, ls=":")
            ax.text(D_U_MM + 0.6, 0.6, "d_u {:g} mm".format(D_U_MM), fontsize=8, color=MUTED)
            dmax = max([c["d"].max() for (k, b), c in data.items() if b == br] + [10.0])
            ax.set_xlim(0, max(dmax * 1.1, 160.0) if mech else dmax * 1.05); ax.set_ylim(0, 20)
            ax.set_title("{} {}".format("(a)" if j == 0 else "(b)", title), loc="left", fontsize=11, fontweight="bold")
        ax.set_xlabel("displacement at the 2.06 m sensor (mm)")
        style_axes(ax)
axs[0, 0].set_ylabel("lateral force on the wall (kN)"); axs[1, 0].set_ylabel("force (kN)")
h, l = [], []
for ax in axs[0]:
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        if ll not in l:
            h.append(hh); l.append(ll)
if SHOW_REACTION and FORCE == "applied":
    h.append(plt.Line2D([], [], color=MUTED, ls=":", lw=1.0)); l.append("dotted: F_tot as logged (misses |N| sin(tilt) at the top)")
if SHOW_UNCONVERGED:
    h.append(plt.Line2D([], [], color=MUTED, ls="none", marker="o", mfc="white")); l.append("hollow: unconverged step (not an equilibrium)")
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8.4, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("OOP pushover, US-1 model (masonry-only push, force = {}) vs the test envelope".format(
    "|a| x m_tot" if FORCE == "applied" else "F_tot"), x=0.01, ha="left", fontsize=12)
fig.tight_layout(rect=(0, 0.095, 1, 1))
fig.savefig(P(OUT_DIR, "fig_pushover_curves.png"), dpi=150); plt.close(fig)
print("-> " + P(OUT_DIR, "fig_pushover_curves.png"))

# ============================ FIG 2: diagnostics =====================
fig, axs = plt.subplots(3, 2, figsize=(14, 10), sharex="col")
for j, br in enumerate(("pos", "neg")):
    for key, label, folder, colr, ls in CASES:
        c = data.get((key, br))
        if c is None:
            continue
        ok = c["ok"]; lw = 1.6 if ls == "-" else 1.2
        axs[0, j].plot(c["d"][ok], c["a"][ok] / 9.81, ls, color=colr, lw=lw, label=label)
        if SHOW_UNCONVERGED and (~ok).any():
            axs[0, j].plot(c["d"][~ok], c["a"][~ok] / 9.81, "o", mfc="white", mec=colr, ms=5)
        axs[1, j].plot(c["d"], c["open"], ls, color=colr, lw=lw)
        axs[1, j].plot(c["d"], c["pocket"], ls=(0, (1, 1.5)), color=colr, lw=lw)
        axs[2, j].plot(c["d"], c["joistN"], ls, color=colr, lw=lw, label=label)
    axs[0, j].set_title("{} {}".format("(a)" if j == 0 else "(b)", "Positive push" if br == "pos" else "Negative push"),
                        loc="left", fontsize=11, fontweight="bold")
    axs[1, j].axhline(POCKET_EXIT_MM, color="#ef4444", lw=0.8, ls=":")
    axs[1, j].text(0.99, POCKET_EXIT_MM + 1, "pocket exit ({:g} mm)".format(POCKET_EXIT_MM), fontsize=8, color="#ef4444",
                   ha="right", transform=axs[1, j].get_yaxis_transform())
    axs[2, j].set_xlabel("displacement at the 2.06 m sensor (mm)")
    for i in range(3):
        style_axes(axs[i, j])
axs[0, 0].set_ylabel("applied acceleration on the masonry (g)")
axs[1, 0].set_ylabel("opening (mm)")
axs[2, 0].set_ylabel("joist support load (kN)")
axs[0, 0].legend(frameon=False, fontsize=8.5)
axs[1, 0].legend([plt.Line2D([], [], color=MUTED, lw=1.4), plt.Line2D([], [], color=MUTED, lw=1.4, ls=(0, (1, 1.5)))],
                 ["largest masonry joint opening", "pocket opening (joist-masonry)"], frameon=False, fontsize=8.5, loc="upper left")
fig.tight_layout()
fig.savefig(P(OUT_DIR, "fig_pushover_diagnostics.png"), dpi=150); plt.close(fig)
print("-> " + P(OUT_DIR, "fig_pushover_diagnostics.png"))
