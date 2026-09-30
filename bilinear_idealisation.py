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
    model : pushover_pos*.csv found here or in pushover_results/ (+ pushover_neg*.csv
            if present; the search list is MODEL_SEARCH), F = |F_tot|,
            converged rows only, curve cut at MODEL_D_CUT_MM (the user asked to
            ignore the abrupt tail).
    fig13 : US1_fig13_digitised.csv, both branches.

OUTPUTS (OUT_DIR = bilinear_idealisation/)
    bilinear_periods.json : per source: per-branch quantities + averaged
                            K1, mu, T1, T_eff (read by route2_revised_3dec.py
                            when PERIOD_SOURCE = "bilinear")
    bilinear_periods.csv  : the same as a table
    fig_bilinear_<source>.png : the construction, one panel per branch
Usage: python bilinear_idealisation.py
"""
import os, csv, json, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
M_EFF          = 1635.0      # kg, as agreed (Table 4 T_sec reproduced); sensitivity printed for M_ALT
M_ALT          = 1318.0      # wall mass, for the sensitivity line only
K_FRAC         = 0.5         # secant point on the rising branch
DROP_FRAC      = 0.8         # ultimate: force fallen to this fraction of Vu
SECANT_TO_DROP = True        # T_eff = T1 sqrt(mu / DROP_FRAC)  (False: T1 sqrt(mu))
MODEL_CSVS     = {"pos": None, "neg": None}   # explicit paths; None = search MODEL_SEARCH below
MODEL_SEARCH   = {"pos": ["pushover_pos_v3.csv", "pushover_pos.csv", "pushover_results_/pushover_pos.csv", "pushover_results/pushover_pos_v3.csv"],
                  "neg": ["pushover_neg_v3.csv", "pushover_neg.csv", "pushover_results/pushover_neg.csv", "pushover_results/pushover_neg_v3.csv"]}
MODEL_D_CUT_MM = 70.0        # ignore the model curve beyond this (abrupt tail)
FIG13_CSV      = "US1_fig13_digitised.csv"
OUT_DIR        = "bilinear_idealisation"
os.makedirs(OUT_DIR, exist_ok=True)
COL = {"model": "#5b7f74", "fig13": "#3b5fa8", "k1": "#b85042", "keff": "#d98c21"}

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
    if len(post) == 0:
        raise ValueError("descending branch never reaches {:.0f}% of Vu (curve ends at {:.1f} mm, {:.1f} kN)".format(
            100 * DROP_FRAC, x[-1], F[-1]))
    k = int(post[0])
    du = float(x[k - 1] + (F[k - 1] - DROP_FRAC * Vu) / (F[k - 1] - F[k]) * (x[k] - x[k - 1]))
    mu = du / dy
    T1 = 2 * np.pi * np.sqrt(m_eff / (K1 * 1e6))
    Te = T1 * np.sqrt(mu / DROP_FRAC if SECANT_TO_DROP else mu)
    return dict(Vu_kN=Vu, d_peak_mm=d_pk, d50_mm=d50, K1_kNmm=K1, dy_mm=dy, du_mm=du, mu=mu,
                du_bracket_mm=(float(x[k - 1]), float(x[k])), F_bracket_kN=(float(F[k - 1]), float(F[k])),
                T1_s=T1, Teff_s=Te, m_eff=m_eff, x=x, F=F)

def periods_from(K1, mu, m_eff):
    T1 = 2 * np.pi * np.sqrt(m_eff / (K1 * 1e6))
    return T1, T1 * np.sqrt(mu / DROP_FRAC if SECANT_TO_DROP else mu)

# ============================ SOURCES ================================
def load_model_branch(path):
    rows = [r for r in csv.DictReader(open(path)) if r.get("converged", "1") == "1"]
    x = np.array([float(r["d_ctrl_mm"]) for r in rows]); F = np.abs(np.array([float(r["F_tot_kN"]) for r in rows]))
    keep = np.abs(x) <= MODEL_D_CUT_MM
    return np.abs(x[keep]), F[keep]

def load_fig13(path):
    d = np.genfromtxt(path, delimiter=",", skip_header=1); F, x = d[:, 1], d[:, 2]
    out = {}
    for name, m, s in (("pos", x > 0, 1.0), ("neg", x < 0, -1.0)):
        out[name] = (s * x[m], s * F[m])
    return out

branches = {"model": {}, "fig13": {}}
for br in ("pos", "neg"):
    cands = [MODEL_CSVS[br]] if MODEL_CSVS.get(br) else MODEL_SEARCH[br] + sorted(glob.glob("**/pushover_{}*.csv".format(br), recursive=True))
    hit = next((c for c in cands if c and os.path.isfile(c)), None)
    if hit:
        x, F = load_model_branch(hit)
        print("model {} branch: {}  ({} converged rows, {:.1f}..{:.1f} mm, peak {:.1f} kN)".format(br, hit, len(x), x.min(), x.max(), F.max()))
        if len(x) < 5:
            print("  !! {} has only {} usable rows -- is this the old 2-row v2.0 file? Use the v3 (Nicolo + servo) pushover_pos.csv.".format(hit, len(x)))
        branches["model"][br] = (x, F)
    elif br == "pos":
        raise SystemExit("no model pushover found. Put the v3 pushover_pos.csv next to this script (or in pushover_results/), "
                         "or set MODEL_CSVS['pos'] to its path. Searched: " + ", ".join(MODEL_SEARCH["pos"]))
if os.path.isfile(FIG13_CSV):
    branches["fig13"] = load_fig13(FIG13_CSV)
else:
    print("fig13: {} not found, envelope block skipped".format(FIG13_CSV))

# ============================ RUN ===================================
result = {"rule": dict(K_FRAC=K_FRAC, DROP_FRAC=DROP_FRAC, secant_to_drop=SECANT_TO_DROP, m_eff=M_EFF, model_d_cut_mm=MODEL_D_CUT_MM)}
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
            raise SystemExit("model pushover could not be idealised (see message above): no descending branch to 0.8Vu?")
        continue
    K1m = float(np.mean([r["K1_kNmm"] for r in per.values()])); mum = float(np.mean([r["mu"] for r in per.values()]))
    T1m, Tem = periods_from(K1m, mum, M_EFF); T1a, Tea = periods_from(K1m, mum, M_ALT)
    result[src] = {"branches_used": sorted(per), "K1_kNmm": K1m, "mu": mum, "T1_s": T1m, "Teff_s": Tem,
                   "T1_s_m_alt": T1a, "Teff_s_m_alt": Tea, "m_alt": M_ALT,
                   "per_branch": {br: {k: v for k, v in r.items() if k not in ("x", "F")} for br, r in per.items()}}
    print("== {} ({} branch{}) ==".format(src, "+".join(sorted(per)), "es averaged" if len(per) > 1 else " only"))
    for br, r in per.items():
        print("  {}: Vu {:.2f} kN at {:.1f} mm | 0.5Vu at {:.2f} mm -> K1 {:.2f} kN/mm, dy {:.2f} mm | 0.8Vu between {:.1f} and {:.1f} mm ({:.1f}->{:.1f} kN) -> du {:.2f} mm | mu {:.2f} | T1 {:.4f} s  Teff {:.4f} s".format(
            br, r["Vu_kN"], r["d_peak_mm"], r["d50_mm"], r["K1_kNmm"], r["dy_mm"], r["du_bracket_mm"][0], r["du_bracket_mm"][1],
            r["F_bracket_kN"][0], r["F_bracket_kN"][1], r["du_mm"], r["mu"], r["T1_s"], r["Teff_s"]))
        table.append([src, br, r["Vu_kN"], r["d_peak_mm"], r["d50_mm"], r["K1_kNmm"], r["dy_mm"], r["du_mm"], r["mu"], r["T1_s"], r["Teff_s"]])
    print("  AVERAGE: K1 {:.2f} kN/mm, mu {:.2f} -> T1 {:.4f} s, Teff {:.4f} s (m_eff {:.0f});  with m {:.0f}: T1 {:.4f}, Teff {:.4f}".format(
        K1m, mum, T1m, Tem, M_EFF, M_ALT, T1a, Tea))
    table.append([src, "average", "", "", "", K1m, "", "", mum, T1m, Tem])

    # ---- figure: the sketch, one panel per branch
    fig, axs = plt.subplots(1, len(per), figsize=(max(9.0, 7.2 * len(per)), 5.4), dpi=140, squeeze=False)
    for ax, (br, r) in zip(axs[0], per.items()):
        x, F, Vu = r["x"], r["F"], r["Vu_kN"]
        ax.plot(x, F, color=COL[src], lw=1.6, label="{} backbone ({})".format("model pushover" if src == "model" else "Fig 13 envelope", br))
        xe = max(r["du_mm"] * 1.6, r["d_peak_mm"] * 1.3)
        ax.plot([0, r["dy_mm"], xe], [0, Vu, Vu], color=COL["k1"], lw=1.4, ls="--", label="K1 = 0.5Vu/d50, plateau at Vu")
        ax.plot([0, r["du_mm"]], [0, DROP_FRAC * Vu], color=COL["keff"], lw=1.4, ls="-.", label="K_eff = 0.8Vu/du")
        ax.plot([r["d50_mm"]], [K_FRAC * Vu], "o", ms=7, color=COL["k1"], mec="white")
        ax.plot([r["du_mm"]], [DROP_FRAC * Vu], "s", ms=7, color=COL["keff"], mec="white")
        ax.plot([r["d_peak_mm"]], [Vu], "^", ms=7, color=COL[src], mec="white")
        for xx, lab in ((r["dy_mm"], "dy {:.1f}".format(r["dy_mm"])), (r["du_mm"], "du {:.1f}".format(r["du_mm"]))):
            ax.axvline(xx, color="#bbbbbb", lw=0.7, ls=":"); ax.text(xx, -0.02 * Vu, lab, ha="center", va="top", fontsize=8.5, color="#444")
        ax.text(r["d50_mm"] + 0.02 * xe, K_FRAC * Vu, "0.5 Vu", fontsize=8.5, va="center", color=COL["k1"])
        ax.text(r["du_mm"] + 0.02 * xe, DROP_FRAC * Vu, "0.8 Vu", fontsize=8.5, va="center", color=COL["keff"])
        ax.text(0.98, 0.30, "Vu {:.1f} kN\nK1 {:.2f} kN/mm\nmu = du/dy = {:.2f}\nT1 {:.3f} s   T_eff {:.3f} s\n(m_eff {:.0f} kg)".format(
            Vu, r["K1_kNmm"], r["mu"], r["T1_s"], r["Teff_s"], M_EFF), transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(fc="white", ec="#dddddd"))
        ax.set_xlim(0, xe); ax.set_ylim(-0.08 * Vu, 1.12 * Vu)
        ax.set_xlabel("OOP displacement at 2.06 m (mm)"); ax.set_ylabel("force (kN)")
        ax.set_title("{} branch".format(br), loc="left", fontsize=10.5, fontweight="bold")
        ax.grid(color="#eeeeee"); ax.legend(fontsize=8, frameon=False, loc="upper right")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.suptitle("Bilinear idealisation, {}\naverage K1 {:.2f} kN/mm, mu {:.2f}  ->  T1 {:.3f} s, T_eff {:.3f} s".format(
        "model pushover" if src == "model" else "Fig 13 envelope", K1m, mum, T1m, Tem), fontsize=11, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.92)); p = os.path.join(OUT_DIR, "fig_bilinear_{}.png".format(src)); plt.savefig(p); plt.close(fig); print("  -> " + p)

json.dump(result, open(os.path.join(OUT_DIR, "bilinear_periods.json"), "w"), indent=1)
with open(os.path.join(OUT_DIR, "bilinear_periods.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["source", "branch", "Vu_kN", "d_peak_mm", "d50_mm", "K1_kNmm", "dy_mm", "du_mm", "mu", "T1_s", "Teff_s"])
    for row in table:
        w.writerow([("{:.4f}".format(v) if isinstance(v, float) else v) for v in row])
print("-> {}/bilinear_periods.json, .csv".format(OUT_DIR))