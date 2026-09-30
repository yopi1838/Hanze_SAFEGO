# -*- coding: ascii -*-
"""
postprocess_route2_from_save.py -- Route 2 on restored Strategy F states vs
the record on the same states.

Pairs (edit CASES):
    route2_from24/route2_FR76_s2p00   vs   stratF_full_results_US1/Run25_*   (scale 2.0)
    route2_from23/route2_FR76_s1p75   vs   stratF_full_results_US1/Run24_*   (scale 1.75)
plus US-1 run 24 from EXP_DATA/processed_globalzero as the experimental
reference on the displacement / shear panels (same reader conventions as
postprocess_routes.py).

Per run, definitions as in postprocess_stratC.py, all WITHIN-RUN (the
Route 2 runs are single continuations, so the control run is treated the
same way for a like-for-like pair):
    peak OOP rel. disp  = max |u - u(t0)|, u = rel_disp_top_exp_mm
                          (else rel_disp_top_mm, else 0.5(Ch3+Ch4)-Ch5)   [mm]
    peak base shear     = max |F - F(t0)|, F = cstav (+ topj_shear if exported) [kN]
    peak tilt           = max |tilt_full_wall - tilt(t0)|                 [deg]
    residual tilt (run) = mean(last 5%) - tilt(t0)                        [deg]

Figures
    (all written to OUT_DIR = route2_from_save_postproc/)
    fig_r2save_metrics.png    : 4 metric panels, bars grouped per pair
                                (control run | Route 2 | US-1 run 24 where applicable)
    fig_r2save_histories.png  : per pair, rel. displacement time histories
                                (control run and Route 2 on their own time axes)
                                and base shear vs displacement loops overlaid
    route2_from_save_metrics.csv
Usage: python postprocess_route2_from_save.py
"""
import os, re, csv, glob
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
CONTROL_DIR = "stratF_full_results_US1"
CASES = [   # (label, Route 2 run folder, control run number, scale)
    ("from run 24, x2.0",  "route2_from24/route2_FR76_s2p00", 25, 2.0),
    ("from run 23, x1.75", "route2_from23/route2_FR76_s1p75", 24, 1.75),
]
SUMMARY_CSVS = ["route2_from24/route2_from_save_summary.csv",
                "route2_from23/route2_from_save_summary.csv"]
EXP_ROOT, EXP_GZ_DIR, EXP_TEST, EXP_RUN = "EXP_DATA", "processed_globalzero", "Test9", 24
EXP_LABEL = "US-1 run 24"
EXP_SHEAR_SIGN = -1.0        # sum(m a) is the inertia force; flipped to the Fig 13 sense
OUT_DIR = "route2_from_save_postproc"     # all outputs go here (created if missing)
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG_M = os.path.join(OUT_DIR, "fig_r2save_metrics.png")
OUT_PNG_H = os.path.join(OUT_DIR, "fig_r2save_histories.png")
OUT_CSV   = os.path.join(OUT_DIR, "route2_from_save_metrics.csv")
TAIL_FRAC = 0.05
KEY_DISP_EXP, KEY_DISP, KEY_TILT, KEY_SHEAR, KEY_SHEAR_TOP = \
    "rel_disp_top_exp_mm", "rel_disp_top_mm", "tilt_full_wall", "cstav", "topj_shear"
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
COL_CTRL, COL_R2, COL_EXP = "#5b7f74", "#b85042", "royalblue"

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
        n = min(len(F), len(Ft)); return t[:n], F[:n] + Ft[:n], "cstav+topj_shear"
    return t, F, "cstav (base+joist only)"

def run_metrics(folder):
    """within-run metrics + the series needed for the figures"""
    out = {"folder": str(folder)}
    t, u, src = edp_series(folder)
    if u is not None:
        out.update(peak_disp_mm=float(np.max(np.abs(u - u[0]))), edp_src=src, t_u=t, u=u - u[0])
    t, F, ssrc = shear_series(folder)
    if F is not None:
        out.update(peak_shear_kN=float(np.max(np.abs(F - F[0]))), shear_src=ssrc, t_F=t, F=F - F[0])
    t, th = find_channel(folder, KEY_TILT)
    if th is not None:
        out.update(max_tilt_deg=float(np.max(np.abs(th - th[0]))),
                   resid_tilt_deg=float(np.mean(th[int((1 - TAIL_FRAC) * len(th)):]) - th[0]))
    return out

def control_folder(run_no):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    return hits[0] if hits else None

def exp_run24():
    """(u_mm rel. to run start, F_kN signed to the Fig 13 sense, peak_mm, peak_kN) or None"""
    for root in (Path(EXP_ROOT), Path.cwd() / EXP_ROOT, Path(__file__).resolve().parent / EXP_ROOT):
        for nm in ("{}Run{:02d}_processed_globalzero.xlsx", "{}Run{}_processed_globalzero.xlsx"):
            p = root / EXP_GZ_DIR / nm.format(EXP_TEST, EXP_RUN)
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
                print("  experiment: {}".format(p))
                return u - u[0], EXP_SHEAR_SIGN * (F - F[0]), float(np.max(np.abs(u - u[0]))), float(np.max(np.abs(F - F[0])))
    print("  ! {}Run{}_processed_globalzero.xlsx not found under {}/{} -- no experimental reference".format(
        EXP_TEST, EXP_RUN, EXP_ROOT, EXP_GZ_DIR))
    return None

# ============================ COLLECT ================================
print("discovery:")
pairs = []
for label, r2_folder, ctrl_run, scale in CASES:
    cf = control_folder(ctrl_run)
    r2 = run_metrics(r2_folder) if os.path.isdir(r2_folder) else None
    ct = run_metrics(cf) if cf else None
    print("  {:20s} Route 2: {}   control run {}: {}".format(
        label, "ok ({}, {})".format(r2.get("edp_src", "-"), r2.get("shear_src", "-")) if r2 and "u" in r2 else "NOT FOUND / no EDP  [{}]".format(r2_folder),
        ctrl_run, "ok" if ct and "u" in ct else "NOT FOUND [{}/Run{:02d}_*]".format(CONTROL_DIR, ctrl_run)))
    pairs.append((label, ctrl_run, scale, r2, ct))
exp = exp_run24()

stage = {}
for sp in SUMMARY_CSVS:
    if os.path.isfile(sp):
        for row in csv.DictReader(open(sp)):
            try:
                stage[float(row["scale"])] = row
            except Exception:
                pass

METRICS = [("peak_disp_mm", "Peak OOP rel. displacement (mm)"),
           ("peak_shear_kN", "Peak base shear (kN)"),
           ("max_tilt_deg", "Peak tilt, full wall (deg)"),
           ("resid_tilt_deg", "Residual tilt within run (deg)")]

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pair", "series", "run_or_scale", "peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg",
                "ratio_route2_over_record_disp", "stage_A_T_s", "stage_A_PGA_g", "stage_B_T_s", "stage_B_PGA_g"])
    for label, ctrl_run, scale, r2, ct in pairs:
        st = stage.get(scale, {})
        for nm, m in (("control record", ct), ("Route 2", r2)):
            if not m:
                continue
            ratio = (r2["peak_disp_mm"] / ct["peak_disp_mm"]) if (nm == "Route 2" and ct and "peak_disp_mm" in ct and "peak_disp_mm" in r2 and ct["peak_disp_mm"]) else ""
            w.writerow([label, nm, ctrl_run if nm == "control record" else scale] +
                       [round(m[k], 4) if k in m else "" for k, _ in METRICS] +
                       [round(ratio, 3) if ratio != "" else "",
                        st.get("T_A_s", ""), st.get("PGA_A_g", ""), st.get("T_B_s", ""), st.get("PGA_B_g", "")])
    if exp:
        w.writerow(["experiment", EXP_LABEL, EXP_RUN, round(exp[2], 3), round(exp[3], 3), "", "", "", "", "", "", ""])
print("-> " + OUT_CSV)

# ============================ FIG 1: metric bars per pair =============
fig, axs = plt.subplots(1, 4, figsize=(19, 5), dpi=140)
for ax, (key, ttl) in zip(axs, METRICS):
    names, vals, cols, x = [], [], [], []
    pos = 0
    for label, ctrl_run, scale, r2, ct in pairs:
        group = [("record\nrun {}".format(ctrl_run), ct, COL_CTRL), ("Route 2\nx{:g}".format(scale), r2, COL_R2)]
        if exp and key in ("peak_disp_mm", "peak_shear_kN") and ctrl_run == EXP_RUN:
            group.append((EXP_LABEL.replace(" run", "\nrun"), {key: exp[2] if key == "peak_disp_mm" else exp[3]}, COL_EXP))
        for nm, m, c in group:
            if m and key in m:
                names.append(nm); vals.append(m[key]); cols.append(c); x.append(pos); pos += 1
        pos += 0.8
    if not vals:
        ax.set_title(ttl + "  (no data)", fontsize=10, loc="left"); continue
    bars = ax.bar(x, vals, color=cols, width=0.7)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, "{:.2f}".format(v), ha="center",
                va="bottom" if v >= 0 else "top", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8)
    ax.set_title(ttl, fontsize=10.5, fontweight="bold", loc="left"); ax.grid(axis="y", color="#e5e5e5")
    ax.axhline(0, color="k", lw=0.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
# pair labels under the groups
fig.text(0.5, 0.01, "   |   ".join("{}: record run {} vs Route 2 x{:g}".format(l, r, s) for l, r, s, _, _ in pairs),
         ha="center", fontsize=9, color="#444")
fig.suptitle("Route 2 on restored Strategy F states vs the record on the same states  (within-run metrics; shear = {})".format(
    next((m.get("shear_src") for _, _, _, r2, ct in pairs for m in (r2, ct) if m and "shear_src" in m), "cstav")),
    fontsize=11, fontweight="bold")
plt.tight_layout(rect=(0, 0.04, 1, 0.95)); plt.savefig(OUT_PNG_M, dpi=140); print("-> " + OUT_PNG_M)

# ============================ FIG 2: histories + loops per pair =======
fig, axs = plt.subplots(len(pairs), 3, figsize=(17, 4.6 * len(pairs)), dpi=140, squeeze=False)
for i, (label, ctrl_run, scale, r2, ct) in enumerate(pairs):
    a0, a1, a2 = axs[i]
    if ct and "u" in ct:
        a0.plot(ct["t_u"], ct["u"], lw=0.7, color=COL_CTRL)
        a0.set_title("record run {} (FR76 x{:g}), peak {:.2f} mm".format(ctrl_run, scale, ct["peak_disp_mm"]), fontsize=10, loc="left")
    else:
        a0.text(0.5, 0.5, "control run {} not found".format(ctrl_run), transform=a0.transAxes, ha="center", color="#888")
    if r2 and "u" in r2:
        a1.plot(r2["t_u"], r2["u"], lw=0.8, color=COL_R2)
        st = stage.get(scale, {})
        a1.set_title("Route 2 {}, peak {:.2f} mm".format(label, r2["peak_disp_mm"]), fontsize=10, loc="left")
        a1.text(0.98, 0.04, "stage A: {} s, {} g\nstage B: {} s, {} g".format(
            st.get("T_A_s", "?"), st.get("PGA_A_g", "?"), st.get("T_B_s", "?"), st.get("PGA_B_g", "?")),
            transform=a1.transAxes, ha="right", va="bottom", fontsize=8, color="#444")
    else:
        a1.text(0.5, 0.5, "Route 2 folder not found", transform=a1.transAxes, ha="center", color="#888")
    for a in (a0, a1):
        a.set_xlabel("t (s)"); a.set_ylabel("OOP rel. disp. at 2.06 m (mm)"); a.axhline(0, color="k", lw=0.5); a.grid(color="#eee")
    ylim = max([abs(v) for m in (r2, ct) if m and "u" in m for v in (m["u"].min(), m["u"].max())] + [1.0]) * 1.1
    a0.set_ylim(-ylim, ylim); a1.set_ylim(-ylim, ylim)
    if exp and ctrl_run == EXP_RUN:
        a2.plot(exp[0], exp[1], lw=0.6, color=COL_EXP, alpha=0.8, label=EXP_LABEL + " (sum(m a) x{:+.0f})".format(EXP_SHEAR_SIGN))
    for nm, m, c, ls in (("record run {}".format(ctrl_run), ct, COL_CTRL, ":"), ("Route 2 " + label, r2, COL_R2, "-")):
        if m and "u" in m and "F" in m:
            n = min(len(m["u"]), len(m["F"])); a2.plot(m["u"][:n], m["F"][:n], lw=0.8, ls=ls, color=c, label=nm, alpha=0.9)
    a2.axhline(0, color="k", lw=0.5); a2.axvline(0, color="k", lw=0.5); a2.grid(color="#eee")
    a2.set_xlabel("OOP rel. displacement (mm)"); a2.set_ylabel("base shear (kN)"); a2.legend(fontsize=8, frameon=False)
    a2.set_title("base shear vs displacement", fontsize=10, loc="left")
    for a in (a0, a1, a2):
        for s in ("top", "right"): a.spines[s].set_visible(False)
fig.suptitle("Route 2 vs record, same restored state: displacement histories and loops", fontsize=11.5, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.96)); plt.savefig(OUT_PNG_H, dpi=140); print("-> " + OUT_PNG_H)
