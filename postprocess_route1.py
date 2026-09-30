# -*- coding: ascii -*-
"""
postprocess_route1_extract.py -- Route 1 BY EXTRACTION (the 1.64 s window of
the FR76 record) vs the full record on the same model state vs US-1, runs
23 and 24.

Reads what route1_extract_3dec.py wrote:
    ROUTE_DIR/Run23_FR76_s1p50/, Run24_FR76_s1p75/   channel CSVs
    ROUTE_DIR/route1_extract_log.csv                  per-run log
    ROUTE_DIR/case.txt                                window, baseline, damping
and the record control runs 23 / 24 from CONTROL_DIR (first that exists of
CONTROL_CANDIDATES), and US-1 runs 23 / 24 from EXP_DATA/processed_globalzero
(same reader as postprocess_routes.py: U_avg [m], base_shear_kN).

Per run, WITHIN-RUN metrics as in postprocess_stratC.py:
    peak OOP rel. disp  = max |u - u(t0)|, u = rel_disp_top_exp_mm
                          (else rel_disp_top_mm, else 0.5(Ch3+Ch4)-Ch5)   [mm]
    peak base shear     = max |F - F(t0)|, F = cstav (+ topj_shear if exported) [kN]
    residual disp       = mean(last 5%) - u(t0)                           [mm]
    peak / residual tilt from tilt_full_wall if exported                  [deg]

ONE figure, 2 rows (run 23, run 24) x 3:
    (a) displacement histories, the extraction SHIFTED onto the record's
        time axis (t + window start) so the peaks line up with the lobe
        it was cut from
    (b) base shear vs displacement loops: record, Route 1, US-1
    (c) peak displacement and peak base shear bars: record | Route 1 | US-1
Base shear on the model is cstav (+topj_shear if exported) -- support-only,
NOT the same quantity as US-1's sum(m a). The loop and bar panels say so.

Usage: python postprocess_route1_extract.py            (case from CASE below)
       ROUTE1_CASE=LS_MAXWELL_1p5 python postprocess_route1_extract.py
Output: ROUTE_DIR/postproc/fig_route1_extract.png, route1_extract_metrics.csv
"""
import os, re, csv, glob, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
CASE = os.environ.get("ROUTE1_CASE", "SS_NODAMP").strip()
ROUTE_DIR = "route1_extract_" + CASE
CONTROL_CANDIDATES = {          # first folder that exists wins; the list is per case
    "SS_NODAMP":      ["stratF_full_SS_NODAMP", "stratF_full_results_US1"],
    "SS_MAXWELL_1p5": ["stratF_full_SS_MAXWELL_1p5"],
    "LS_NODAMP":      ["stratF_full_LS_NODAMP", "stratF_full_results_US1"],
    "LS_MAXWELL_1p5": ["stratF_full_LS_MAXWELL_1p5"],
}.get(CASE, ["stratF_full_results_US1"])
RUNS = [(23, 1.50), (24, 1.75)]
EXP_ROOT, EXP_GZ_DIR, EXP_TEST = "EXP_DATA", "processed_globalzero", "Test9"
EXP_LABEL = "US-1"
EXP_SHEAR_SIGN = -1.0        # sum(m a) is the inertia force; flipped to the Fig 13 sense
MODEL_SHEAR_SIGN = -1.0      # cstav (+topj_shear) is the support reaction; flipped so +F goes with +u
TAIL_FRAC = 0.05
KEY_DISP_EXP, KEY_DISP, KEY_TILT, KEY_SHEAR, KEY_SHEAR_TOP = \
    "rel_disp_top_exp_mm", "rel_disp_top_mm", "tilt_full_wall", "cstav", "topj_shear"
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
OUT_DIR = os.path.join(ROUTE_DIR, "postproc")
OUT_PNG = os.path.join(OUT_DIR, "fig_route1_extract.png")
OUT_CSV = os.path.join(OUT_DIR, "route1_extract_metrics.csv")
# validated categorical set (dataviz validate_palette.js, light surface): all checks PASS
COL_CTRL, COL_R1, COL_EXP, C_INK, C_MUTE = "#15803d", "#c2410c", "#1d4ed8", "#3f3f46", "#a1a1aa"

# ============================ READERS ================================
def read_hist_csv(p):
    try:
        d = np.genfromtxt(str(p), skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2:
        return None, None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return (d[:, 0], d[:, 1]) if len(d) >= 5 else (None, None)

def find_channel(folder, key):
    hits = sorted(Path(folder).glob("*" + key + "*.csv"))
    return read_hist_csv(hits[0]) if hits else (None, None)

def edp_series(folder):
    for key in (KEY_DISP_EXP, KEY_DISP):
        t, u = find_channel(folder, key)
        if u is not None:
            return t, u, key
    t3, c3 = find_channel(folder, CH3); _, c4 = find_channel(folder, CH4); _, c5 = find_channel(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, None, None
    n = min(len(c3), len(c4), len(c5))
    return t3[:n], (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0, "0.5(Ch3+Ch4)-Ch5"

def shear_series(folder):
    t, F = find_channel(folder, KEY_SHEAR)
    if F is None:
        return None, None, None
    _, Ft = find_channel(folder, KEY_SHEAR_TOP)
    if Ft is not None:
        n = min(len(F), len(Ft)); return t[:n], MODEL_SHEAR_SIGN * (F[:n] + Ft[:n]), "cstav+topj_shear"
    return t, MODEL_SHEAR_SIGN * F, "cstav (base+joist only)"

def run_metrics(folder):
    out = {"folder": str(folder)}
    t, u, src = edp_series(folder)
    if u is not None:
        u0 = u - u[0]; n = max(1, int(TAIL_FRAC * len(u0)))
        out.update(peak_disp_mm=float(np.max(np.abs(u0))), resid_disp_mm=float(np.mean(u0[-n:])),
                   edp_src=src, t_u=t - t[0], u=u0)
    t, F, ssrc = shear_series(folder)
    if F is not None:
        out.update(peak_shear_kN=float(np.max(np.abs(F - F[0]))), shear_src=ssrc, t_F=t - t[0], F=F - F[0])
    t, th = find_channel(folder, KEY_TILT)
    if th is not None:
        out.update(max_tilt_deg=float(np.max(np.abs(th - th[0]))),
                   resid_tilt_deg=float(np.mean(th[int((1 - TAIL_FRAC) * len(th)):]) - th[0]))
    return out

def find_run_folder(root, run_no):
    hits = [d for d in glob.glob(os.path.join(root, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    return hits[0] if hits else None

def exp_run(run_no):
    """(u_mm rel. to run start, F_kN in the Fig 13 sense, peak_mm, peak_kN) or None"""
    for root in (Path(EXP_ROOT), Path.cwd() / EXP_ROOT, Path(__file__).resolve().parent / EXP_ROOT):
        for nm in ("{}Run{:02d}_processed_globalzero.xlsx", "{}Run{}_processed_globalzero.xlsx"):
            p = root / EXP_GZ_DIR / nm.format(EXP_TEST, run_no)
            if p.is_file():
                try:
                    import openpyxl
                except ImportError:
                    print("  ! openpyxl missing -- experiment skipped"); return None
                ws = openpyxl.load_workbook(str(p), read_only=True, data_only=True).worksheets[0]
                rows = ws.iter_rows(values_only=True); hdr = [str(h) for h in next(rows)]
                d = np.array([r for r in rows if r and r[0] is not None], dtype=float)
                col = {h: d[:, j] for j, h in enumerate(hdr)}
                if "U_avg" not in col or "base_shear_kN" not in col:
                    print("  ! {} lacks U_avg / base_shear_kN".format(p.name)); return None
                u = col["U_avg"] * 1000.0; F = col["base_shear_kN"]
                print("  experiment run {}: {}".format(run_no, p))
                return u - u[0], EXP_SHEAR_SIGN * (F - F[0]), float(np.max(np.abs(u - u[0]))), float(np.max(np.abs(F - F[0])))
    print("  ! {}Run{}_processed_globalzero.xlsx not found under {}/{}".format(EXP_TEST, run_no, EXP_ROOT, EXP_GZ_DIR))
    return None

# ============================ COLLECT ================================
if not os.path.isdir(ROUTE_DIR):
    raise SystemExit("route folder not found: {}  (set CASE / ROUTE1_CASE)".format(ROUTE_DIR))
os.makedirs(OUT_DIR, exist_ok=True)
CONTROL_DIR = next((d for d in CONTROL_CANDIDATES if os.path.isdir(d)), None)
case = {}
if os.path.isfile(os.path.join(ROUTE_DIR, "case.txt")):
    try:
        case = json.load(open(os.path.join(ROUTE_DIR, "case.txt")))
    except Exception:
        pass
WIN_A = float(case.get("window", [12.461, 14.106])[0])
log = {}
if os.path.isfile(os.path.join(ROUTE_DIR, "route1_extract_log.csv")):
    for row in csv.DictReader(open(os.path.join(ROUTE_DIR, "route1_extract_log.csv"))):
        try:
            log[int(float(row["run"]))] = row
        except Exception:
            pass

print("Route 1 by extraction, case {}".format(CASE))
print("  route   : {}   window {:.3f}-{:.3f} s ({:.3f} s), baseline {}, damping: {}".format(
    ROUTE_DIR, *case.get("window", [np.nan, np.nan]), case.get("dur_s", np.nan),
    "zeroed" if case.get("baseline_zero", True) else "RAW", case.get("damping") or "none"))
print("  control : {}".format(CONTROL_DIR or "NONE of " + ", ".join(CONTROL_CANDIDATES)))
if case.get("start_how") == "FALLBACK":
    print("  ! the route started from the UNDAMAGED save (fallback) -- runs 23-24 on a virgin wall, not comparable to the control")
rows = []
for run_no, sc in RUNS:
    rf = find_run_folder(ROUTE_DIR, run_no); cf = find_run_folder(CONTROL_DIR, run_no) if CONTROL_DIR else None
    r1 = run_metrics(rf) if rf else None
    ct = run_metrics(cf) if cf else None
    ex = exp_run(run_no)
    print("  run {:02d} x{:.2f}: Route 1 {} | record {} | US-1 {}".format(
        run_no, sc,
        "{:.2f} mm ({})".format(r1["peak_disp_mm"], r1["edp_src"]) if r1 and "u" in r1 else "NOT FOUND",
        "{:.2f} mm".format(ct["peak_disp_mm"]) if ct and "u" in ct else "NOT FOUND",
        "{:.2f} mm".format(ex[2]) if ex else "n/a"))
    rows.append((run_no, sc, r1, ct, ex))
if not any(r1 for _, _, r1, _, _ in rows):
    raise SystemExit("no Route 1 run folders with an EDP under {}".format(ROUTE_DIR))

METRICS = ["peak_disp_mm", "resid_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg"]
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "scale", "series"] + METRICS + ["edp_src", "shear_src", "ratio_route1_over_record_disp",
                "ratio_route1_over_US1_disp", "vpeak_applied_mps", "net_base_disp_mm", "window_s", "case"])
    for run_no, sc, r1, ct, ex in rows:
        lg = log.get(run_no, {})
        for nm, m in (("record", ct), ("Route 1 extraction", r1)):
            if not m:
                continue
            rr = (r1["peak_disp_mm"] / ct["peak_disp_mm"]) if (nm != "record" and ct and "peak_disp_mm" in ct and "peak_disp_mm" in r1 and ct["peak_disp_mm"]) else ""
            re_ = (r1["peak_disp_mm"] / ex[2]) if (nm != "record" and ex and "peak_disp_mm" in r1 and ex[2]) else ""
            w.writerow([run_no, sc, nm] + [round(m[k], 4) if k in m else "" for k in METRICS] +
                       [m.get("edp_src", ""), m.get("shear_src", ""),
                        round(rr, 3) if rr != "" else "", round(re_, 3) if re_ != "" else "",
                        lg.get("vpeak_applied_mps", "") if nm != "record" else "",
                        lg.get("net_base_disp_mm", "") if nm != "record" else "",
                        case.get("dur_s", "") if nm != "record" else "", CASE])
        if ex:
            w.writerow([run_no, sc, EXP_LABEL, round(ex[2], 3), "", round(ex[3], 3), "", "", "", "", "", "", "", "", "", ""])
print("-> " + OUT_CSV)

# ============================ FIGURE =================================
fig, axs = plt.subplots(len(rows), 3, figsize=(18.5, 5.4 * len(rows)), dpi=140, squeeze=False,
                        gridspec_kw={"width_ratios": [1.5, 1.0, 0.9]})
for i, (run_no, sc, r1, ct, ex) in enumerate(rows):
    a0, a1, a2 = axs[i]
    # (a) histories, extraction shifted onto the record's clock
    if ct and "u" in ct:
        a0.plot(ct["t_u"], ct["u"], lw=0.8, color=COL_CTRL, label="record on the model, run {} ({:.2f} mm)".format(run_no, ct["peak_disp_mm"]))
    if r1 and "u" in r1:
        a0.plot(r1["t_u"] + WIN_A, r1["u"], lw=1.4, color=COL_R1, label="Route 1 extraction ({:.2f} mm)".format(r1["peak_disp_mm"]))
        a0.axvspan(WIN_A, WIN_A + float(case.get("dur_s", 1.64)), color="#f3e6d8", alpha=0.6, lw=0)
    if ex:
        a0.axhline(ex[2], color=COL_EXP, lw=1.0, ls=":", label="{} run {} peak {:.2f} mm".format(EXP_LABEL, run_no, ex[2]))
        a0.axhline(-ex[2], color=COL_EXP, lw=1.0, ls=":")
    a0.axhline(0, color="k", lw=0.5); a0.grid(color="#eee")
    a0.set_xlabel("time on the record's clock (s)"); a0.set_ylabel("OOP rel. displacement at 2.06 m (mm)")
    a0.set_title("(a) run {} FR76 x{:.2f}: the extraction placed where it was cut from".format(run_no, sc), fontsize=10.5, loc="left", fontweight="bold")
    a0.legend(fontsize=8.5, frameon=False, loc="upper left")
    # (b) loops
    if ex:
        a1.plot(ex[0], ex[1], lw=0.5, color=COL_EXP, alpha=0.7, label="{} run {} (sum(m a) x{:+.0f})".format(EXP_LABEL, run_no, EXP_SHEAR_SIGN))
    for nm, m, c, ls, lw in (("record run {}".format(run_no), ct, COL_CTRL, ":", 0.8), ("Route 1 extraction", r1, COL_R1, "-", 1.1)):
        if m and "u" in m and "F" in m:
            n = min(len(m["u"]), len(m["F"])); a1.plot(m["u"][:n], m["F"][:n], lw=lw, ls=ls, color=c, label=nm, alpha=0.9)
    a1.axhline(0, color="k", lw=0.5); a1.axvline(0, color="k", lw=0.5); a1.grid(color="#eee")
    a1.set_xlabel("OOP rel. displacement (mm)"); a1.set_ylabel("base shear (kN)")
    a1.set_title("(b) base shear vs displacement", fontsize=10.5, loc="left", fontweight="bold"); a1.legend(fontsize=8, frameon=False)
    # (c) bars
    names, vals, cols = [], [], []
    for nm, v, c in (("record", ct.get("peak_disp_mm") if ct else None, COL_CTRL),
                     ("Route 1", r1.get("peak_disp_mm") if r1 else None, COL_R1),
                     (EXP_LABEL, ex[2] if ex else None, COL_EXP)):
        if v is not None:
            names.append(nm); vals.append(v); cols.append(c)
    x = np.arange(len(vals)); b = a2.bar(x, vals, color=cols, width=0.65)
    for bb, v in zip(b, vals):
        a2.text(bb.get_x() + bb.get_width() / 2, v, "{:.1f}".format(v), ha="center", va="bottom", fontsize=9, color=C_INK)
    a2.set_xticks(x); a2.set_xticklabels(names, fontsize=9); a2.set_ylabel("peak OOP rel. displacement (mm)")
    a2.set_title("(c) peak displacement, run {}".format(run_no), fontsize=10.5, loc="left", fontweight="bold")
    a2.grid(axis="y", color="#eee")
    if r1 and ct and "peak_disp_mm" in r1 and "peak_disp_mm" in ct and ct["peak_disp_mm"]:
        a2.set_xlabel("Route 1 / record = {:.2f}".format(r1["peak_disp_mm"] / ct["peak_disp_mm"]) +
                      ("   Route 1 / US-1 = {:.2f}".format(r1["peak_disp_mm"] / ex[2]) if ex else ""),
                      fontsize=9.5, color=C_INK)
        a2.set_ylim(0, max(vals) * 1.15)
    for a in (a0, a1, a2):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
shear_src = next((m["shear_src"] for _, _, r1, ct, _ in rows for m in (r1, ct) if m and "shear_src" in m), "cstav")
fig.suptitle("Route 1 by extraction ({:.2f} s window of FR76) vs the record on the same state vs US-1  --  case {}".format(
    float(case.get("dur_s", 1.64)), CASE), fontsize=12, fontweight="bold")
fig.text(0.01, 0.01, "Within-run metrics, all relative to the run's own start. Model base shear = {} (support-only): not the same quantity as US-1's sum(m a). "
                     "Baseline {}; control {}.".format(shear_src, "zeroed" if case.get("baseline_zero", True) else "RAW (base drifts)", CONTROL_DIR or "none"),
         fontsize=8.5, color="#555")
plt.tight_layout(rect=(0, 0.03, 1, 0.96)); plt.savefig(OUT_PNG); print("-> " + OUT_PNG)