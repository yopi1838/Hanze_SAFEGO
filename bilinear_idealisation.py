# -*- coding: ascii -*-
"""
bilinear_idealisation.py -- T1 and T_eff from a pushover / envelope curve by
the sketch's rule (Tomazevic-type bilinear idealisation):

    Vu      : peak force of the branch
    d50     : displacement where the RISING branch first reaches 0.5 Vu
    K1      : 0.5 Vu / d50                      (initial secant stiffness)
    dy      : Vu / K1                           (yield of the elasto-plastic idealisation)
    du      : displacement on the DESCENDING branch where F first drops to 0.8 Vu
    mu      : du / dy
    T1      : 2 pi sqrt(m_eff / K1)
    T_eff   : T1 sqrt(mu / 0.8)                 (secant to (du, 0.8 Vu), as agreed)
Two branches (pos / neg) are idealised separately and then K1 and mu are
AVERAGED; T1 and T_eff come from the averages. With one branch only, that
branch is used and the output says so.

SOURCES (each produces its own block in the outputs)
    model : pushover_pos*.csv (+ pushover_neg*.csv if present) in PUSHOVER_DIR,
            F = |F_tot|, converged rows only, curve cut at MODEL_D_CUT_MM.
            Lines starting with '#' (the v5 driver writes one) are skipped.
    fig13 : US1_fig13_digitised.csv, both branches (skipped if absent).

OUTPUTS (OUT_DIR)
    bilinear_periods.json : per source: per-branch quantities + averaged K1, mu, T1, T_eff
                            (read by the Method 2 / Route 2 drivers as T1_JSON)
    bilinear_periods.csv  : the same as a table
    fig_bilinear_<source>.png : the construction, one panel per branch
Usage: python bilinear_idealisation.py        (edit the CONFIG block)
"""
import os, csv, json, glob, io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
PUSHOVER_DIR   = "pushover_results_v5_ft0_masonrypush"   # the pushover driver's OUT_DIR
OUT_DIR        = "bilinear_idealisation_v12"             # where the json / csv / figures go
MODEL_CSVS     = {"pos": None, "neg": None}   # explicit paths; None = pushover_pos*.csv / pushover_neg*.csv in PUSHOVER_DIR
M_EFF          = 1635.0      # kg, as agreed (Table 4 T_sec reproduced); sensitivity printed for M_ALT
M_ALT          = 1318.0      # wall mass, for the sensitivity line only
K_FRAC         = 0.5         # secant point on the rising branch
DROP_FRAC      = 0.8         # ultimate: force fallen to this fraction of Vu
SECANT_TO_DROP = True        # T_eff = T1 sqrt(mu / DROP_FRAC)  (False: T1 sqrt(mu))
MODEL_D_CUT_MM = 90.0        # ignore the model curve beyond this (abrupt tail)
DU_FALLBACK    = "end"     # "error": stop if the descending branch never reaches 0.8 Vu
                             # "end"  : take the last point as du and FLAG it everywhere (json, csv, figure)
FIG13_CSV      = "US1_fig13_digitised.csv"
COL = {"model": "#5b7f74", "fig13": "#3b5fa8", "k1": "#b85042", "keff": "#d98c21"}

if not os.path.isdir(PUSHOVER_DIR):
    raise SystemExit("pushover directory not found: " + PUSHOVER_DIR)
os.makedirs(OUT_DIR, exist_ok=True)

# ============================ THE RULE ===============================
def idealise(x, F, m_eff=M_EFF):
    """x (mm, >=0, increasing), F (kN, >=0). Returns dict or raises."""
    x = np.asarray(x, float); F = np.asarray(F, float)
    o = np.argsort(x); x, F = x[o], F[o]
    ip = int(np.argmax(F)); Vu, d_pk = float(F[ip]), float(x[ip])
    rise_x, rise_F = x[:ip + 1], F[:ip + 1]
    j = int(np.argmax(rise_F >= K_FRAC * Vu))
    if j == 0:
        raise ValueError("first point already above K_FRAC*Vu; curve starts too late")
    d50 = float(rise_x[j - 1] + (K_FRAC * Vu - rise_F[j - 1]) / (rise_F[j] - rise_F[j - 1]) * (rise_x[j] - rise_x[j - 1]))
    K1 = K_FRAC * Vu / d50; dy = Vu / K1
    post = np.where((np.arange(len(x)) > ip) & (F <= DROP_FRAC * Vu))[0]
    du_flag = ""
    if len(post) == 0:
        msg = "descending branch never reaches {:.0f}% of Vu (curve ends at {:.1f} mm, {:.1f} kN = {:.0f}% of Vu)".format(
            100 * DROP_FRAC, x[-1], F[-1], 100 * F[-1] / Vu)
        if DU_FALLBACK == "end":
            du = float(x[-1]); du_flag = "du = curve end (0.8Vu not reached; F_end/Vu = {:.2f})".format(F[-1] / Vu)
            br_x, br_F = (float(x[-2]), float(x[-1])), (float(F[-2]), float(F[-1]))
            print("  !! " + msg + " -> " + du_flag)
        else:
            raise ValueError(msg + ". Run the pushover further, or set DU_FALLBACK = 'end'.")
    else:
        k = int(post[0])
        du = float(x[k - 1] + (F[k - 1] - DROP_FRAC * Vu) / (F[k - 1] - F[k]) * (x[k] - x[k - 1]))
        br_x, br_F = (float(x[k - 1]), float(x[k])), (float(F[k - 1]), float(F[k]))
    mu = du / dy
    T1 = 2 * np.pi * np.sqrt(m_eff / (K1 * 1e6))
    Te = T1 * np.sqrt(mu / DROP_FRAC if SECANT_TO_DROP else mu)
    return dict(Vu_kN=Vu, d_peak_mm=d_pk, d50_mm=d50, K1_kNmm=K1, dy_mm=dy, du_mm=du, mu=mu,
                du_bracket_mm=br_x, F_bracket_kN=br_F, du_flag=du_flag,
                T1_s=T1, Teff_s=Te, m_eff=m_eff, x=x, F=F)

def periods_from(K1, mu, m_eff):
    T1 = 2 * np.pi * np.sqrt(m_eff / (K1 * 1e6))
    return T1, T1 * np.sqrt(mu / DROP_FRAC if SECANT_TO_DROP else mu)

# ============================ SOURCES ================================
def read_csv_skip_comments(path):
    with open(path) as f:
        text = "".join(line for line in f if not line.lstrip().startswith("#"))
    return list(csv.DictReader(io.StringIO(text)))

def load_model_branch(path):
    rows = [r for r in read_csv_skip_comments(path) if r.get("converged", "1") == "1"]
    x = np.array([float(r["d_ctrl_mm"]) for r in rows]); F = np.abs(np.array([float(r["F_tot_kN"]) for r in rows]))
    keep = np.abs(x) <= MODEL_D_CUT_MM
    return np.abs(x[keep]), F[keep]

def load_fig13(path):
    d = np.genfromtxt(path, delimiter=",", skip_header=1); F, x = d[:, 1], d[:, 2]
    out = {}
    for name, m, s in (("pos", x > 0, 1.0), ("neg", x < 0, -1.0)):
        out[name] = (s * x[m], s * F[m])
    return out

def find_branch(br):
    if MODEL_CSVS.get(br):
        return MODEL_CSVS[br] if os.path.isfile(MODEL_CSVS[br]) else None
    cands = sorted(glob.glob(os.path.join(PUSHOVER_DIR, "pushover_{}*.csv".format(br))))
    return cands[0] if cands else None

branches = {"model": {}, "fig13": {}}
for br in ("pos", "neg"):
    hit = find_branch(br)
    if hit:
        x, F = load_model_branch(hit)
        print("model {} branch: {}  ({} converged rows, {:.1f}..{:.1f} mm, peak {:.1f} kN)".format(br, hit, len(x), x.min(), x.max(), F.max()))
        if len(x) < 5:
            print("  !! {} has only {} usable rows".format(hit, len(x)))
        branches["model"][br] = (x, F)
    elif br == "pos":
        raise SystemExit("no pushover_pos*.csv in {} (or MODEL_CSVS['pos'] path invalid)".format(PUSHOVER_DIR))
    else:
        print("model neg branch: none in {} -> positive branch only".format(PUSHOVER_DIR))
if os.path.isfile(FIG13_CSV):
    branches["fig13"] = load_fig13(FIG13_CSV)
else:
    print("fig13: {} not found, envelope block skipped".format(FIG13_CSV))

# ============================ RUN ===================================
result = {"rule": dict(K_FRAC=K_FRAC, DROP_FRAC=DROP_FRAC, secant_to_drop=SECANT_TO_DROP, m_eff=M_EFF,
                       model_d_cut_mm=MODEL_D_CUT_MM, du_fallback=DU_FALLBACK, pushover_dir=PUSHOVER_DIR)}
table = []
for src, brs in branches.items():
    if not brs:
        print("{}: no data".format(src)); continue
    per = {}
    for br, (x, F) in brs.items():
        try:
            per[br] = idealise(x, F)
        except ValueError as e:
            print("{} {}: {}".format(src, br, e)); continue
    if not per:
        if src == "model":
            raise SystemExit("model pushover could not be idealised (see message above)")
        continue
    K1m = float(np.mean([r["K1_kNmm"] for r in per.values()])); mum = float(np.mean([r["mu"] for r in per.values()]))
    T1m, Tem = periods_from(K1m, mum, M_EFF); T1a, Tea = periods_from(K1m, mum, M_ALT)
    flags = {br: r["du_flag"] for br, r in per.items() if r["du_flag"]}
    result[src] = {"branches_used": sorted(per), "K1_kNmm": K1m, "mu": mum, "T1_s": T1m, "Teff_s": Tem,
                   "T1_s_m_alt": T1a, "Teff_s_m_alt": Tea, "m_alt": M_ALT, "du_flags": flags,
                   "per_branch": {br: {k: v for k, v in r.items() if k not in ("x", "F")} for br, r in per.items()}}
    print("== {} ({} branch{}) ==".format(src, "+".join(sorted(per)), "es averaged" if len(per) > 1 else " only"))
    for br, r in per.items():
        print("  {}: Vu {:.2f} kN at {:.1f} mm | 0.5Vu at {:.2f} mm -> K1 {:.2f} kN/mm, dy {:.2f} mm | "
              "0.8Vu between {:.1f} and {:.1f} mm ({:.1f}->{:.1f} kN) -> du {:.2f} mm{} | mu {:.2f} | T1 {:.4f} s  Teff {:.4f} s".format(
            br, r["Vu_kN"], r["d_peak_mm"], r["d50_mm"], r["K1_kNmm"], r["dy_mm"], r["du_bracket_mm"][0], r["du_bracket_mm"][1],
            r["F_bracket_kN"][0], r["F_bracket_kN"][1], r["du_mm"], " [FLAGGED]" if r["du_flag"] else "", r["mu"], r["T1_s"], r["Teff_s"]))
        table.append([src, br, r["Vu_kN"], r["d_peak_mm"], r["d50_mm"], r["K1_kNmm"], r["dy_mm"], r["du_mm"], r["mu"], r["T1_s"], r["Teff_s"], r["du_flag"]])
    print("  AVERAGE: K1 {:.2f} kN/mm, mu {:.2f} -> T1 {:.4f} s, Teff {:.4f} s (m_eff {:.0f});  with m {:.0f}: T1 {:.4f}, Teff {:.4f}".format(
        K1m, mum, T1m, Tem, M_EFF, M_ALT, T1a, Tea))
    if flags:
        print("  !! T_eff for this source rests on a FLAGGED du: " + "; ".join("{}: {}".format(b, f) for b, f in flags.items()))
    table.append([src, "average", "", "", "", K1m, "", "", mum, T1m, Tem, "; ".join(flags.values())])

    # ----