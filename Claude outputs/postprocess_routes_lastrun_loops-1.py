# -*- coding: ascii -*-
"""
postprocess_routes_lastrun_loops.py -- base shear vs OOP displacement loops,
LAST RUN ONLY: the record on the model (Strategy F run 25) against Route 1
and Route 2 applied to the same restored state (end of run 24) at the same
level (FR76 x 2.0).

Series (edit ROUTES; a missing folder is skipped with a note):
    control      stratF_full_results_US1/Run25_*            (the record)
    Route 1      route1_from24/route1_FR76_s2p00_N1p0
    Route 2      route2_from24/route2_FR76_s2p00             (N_B = 1)
    Route 2 N=3  route2_from24_NB3/route2_FR76_s2p00         (optional)
    Route 2 N=1.5 route2_from24_NB1p5/route2_FR76_s2p00      (optional)
Experiment: US-1 run 24 loop from EXP_DATA/processed_globalzero (thin, grey)
as a reference only -- run 25 was not applied to US-1.

Displacement = rel_disp_top_exp_mm (else rel_disp_top_mm, else
0.5(Ch3+Ch4)-Ch5), zeroed at the start of the run. Shear = cstav + topj_shear
when both exist, else cstav (base + joists only; stated in the title). Both
zeroed at the start of the run.

Figures (OUT_DIR = routes_lastrun_postproc/):
    fig_lastrun_loops.png   : panel per route, record loop underneath in every panel
    fig_lastrun_overlay.png : all on one axis
    lastrun_loop_metrics.csv: peak |u|, peak |F|, F at peak |u|, u at peak |F|
Usage: python postprocess_routes_lastrun_loops.py
"""
import os, csv, glob
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
LAST_RUN    = 25
CONTROL_DIR = "stratF_full_results_US1"
ROUTES = [   # (label, folder, colour)
    ("Route 1 (PGV sine, N=1)",     "route1_from24/route1_FR76_s2p00_N1p0", "#d98c21"),
    ("Route 2 (two-stage, N_B=1)",  "route2_from24/route2_FR76_s2p00",      "#b85042"),
    ("Route 2 (two-stage, N_B=3)",  "route2_from24_NB3/route2_FR76_s2p00",  "#7a2e22"),
    ("Route 2 (two-stage, N_B=1.5)","route2_from24_NB1p5/route2_FR76_s2p00","#e0a48a"),
]
COL_CTRL = "#5b7f74"
EXP_ROOT, EXP_GZ_DIR, EXP_TEST, EXP_RUN = "EXP_DATA", "processed_globalzero", "Test9", 24
EXP_SHEAR_SIGN = -1.0
MODEL_SHEAR_SIGN = -1.0      # cstav (+topj_shear) is the support reaction; flipped so +F goes with +u like the US-1 loops (Fig 13 sense)
SHOW_EXP = True
OUT_DIR = "routes_lastrun_postproc"
os.makedirs(OUT_DIR, exist_ok=True)
KEY_DISP_EXP, KEY_DISP, KEY_SHEAR, KEY_SHEAR_TOP = "rel_disp_top_exp_mm", "rel_disp_top_mm", "cstav", "topj_shear"
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"

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

def loop(folder):
    """(u_mm, F_kN, shear_src) zeroed at the run start, or None"""
    u = None
    for key in (KEY_DISP_EXP, KEY_DISP):
        _, u = find_channel(folder, key)
        if u is not None:
            break
    if u is None:
        _, c3 = find_channel(folder, CH3); _, c4 = find_channel(folder, CH4); _, c5 = find_channel(folder, CH5)
        if c3 is None or c4 is None or c5 is None:
            return None
        n = min(len(c3), len(c4), len(c5)); u = (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0
    _, F = find_channel(folder, KEY_SHEAR)
    if F is None:
        return None
    _, Ft = find_channel(folder, KEY_SHEAR_TOP)
    src = "cstav (base+joists only)"
    if Ft is not None:
        n = min(len(F), len(Ft)); F = F[:n] + Ft[:n]; src = "cstav + topj_shear (total)"
    n = min(len(u), len(F))
    return u[:n] - u[0], MODEL_SHEAR_SIGN * (F[:n] - F[0]), src

def control_folder(run_no):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    return hits[0] if hits else None

def exp_loop():
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
                    return None
                u = col["U_avg"] * 1000.0; F = col["base_shear_kN"]
                print("  experiment: {}".format(p))
                return u - u[0], EXP_SHEAR_SIGN * (F - F[0])
    print("  ! no {}Run{} workbook under {}/{} -- experiment omitted".format(EXP_TEST, EXP_RUN, EXP_ROOT, EXP_GZ_DIR))
    return None

def metrics(u, F):
    i = int(np.argmax(np.abs(u))); j = int(np.argmax(np.abs(F)))
    return dict(peak_u_mm=float(abs(u[i])), F_at_peak_u_kN=float(F[i]), peak_F_kN=float(abs(F[j])), u_at_peak_F_mm=float(u[j]))

# ============================ COLLECT ================================
print("discovery:")
cf = control_folder(LAST_RUN)
ctrl = loop(cf) if cf else None
print("  control run {}: {}".format(LAST_RUN, "ok [{}]".format(cf) if ctrl else "NOT FOUND [{}/Run{:02d}_*]".format(CONTROL_DIR, LAST_RUN)))
series = []
for label, folder, colr in ROUTES:
    L = loop(folder) if os.path.isdir(folder) else None
    print("  {:30s} {}".format(label, "ok ({})".format(L[2]) if L else "not found / no channels [{}] -- skipped".format(folder)))
    if L:
        series.append((label, colr, L))
exp = exp_loop() if SHOW_EXP else None
shear_lbl = ctrl[2] if ctrl else (series[0][2][2] if series else KEY_SHEAR)

with open(os.path.join(OUT_DIR, "lastrun_loop_metrics.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["series", "peak_u_mm", "F_at_peak_u_kN", "peak_F_kN", "u_at_peak_F_mm", "shear"])
    if ctrl:
        m = metrics(ctrl[0], ctrl[1]); w.writerow(["record run {}".format(LAST_RUN)] + [round(m[k], 3) for k in ("peak_u_mm", "F_at_peak_u_kN", "peak_F_kN", "u_at_peak_F_mm")] + [ctrl[2]])
    for label, colr, L in series:
        m = metrics(L[0], L[1]); w.writerow([label] + [round(m[k], 3) for k in ("peak_u_mm", "F_at_peak_u_kN", "peak_F_kN", "u_at_peak_F_mm")] + [L[2]])
    if exp is not None:
        m = metrics(exp[0], exp[1]); w.writerow(["US-1 run {} (test)".format(EXP_RUN)] + [round(m[k], 3) for k in ("peak_u_mm", "F_at_peak_u_kN", "peak_F_kN", "u_at_peak_F_mm")] + ["sum(m a) x{:+.0f}".format(EXP_SHEAR_SIGN)])
print("-> lastrun_loop_metrics.csv")

# shared limits
allu = [np.max(np.abs(s[2][0])) for s in series] + ([np.max(np.abs(ctrl[0]))] if ctrl else []) + ([np.max(np.abs(exp[0]))] if exp is not None else [])
allF = [np.max(np.abs(s[2][1])) for s in series] + ([np.max(np.abs(ctrl[1]))] if ctrl else []) + ([np.max(np.abs(exp[1]))] if exp is not None else [])
ulim = 1.08 * max(allu + [1.0]); Flim = 1.08 * max(allF + [1.0])

# ============================ FIG 1: panel per route ==================
n = max(1, len(series))
fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 5.0), dpi=140, sharex=True, sharey=True, squeeze=False)
for ax, (label, colr, L) in zip(axs[0], series if series else [("no route folders found", "k", None)]):
    if exp is not None:
        ax.plot(exp[0], exp[1], lw=0.5, color="#999999", alpha=0.7, label="US-1 run {} (test, reference)".format(EXP_RUN))
    if ctrl:
        ax.plot(ctrl[0], ctrl[1], lw=0.8, color=COL_CTRL, alpha=0.9, label="record, run {}".format(LAST_RUN))
    if L:
        ax.plot(L[0], L[1], lw=1.0, color=colr, label=label)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.grid(color="#eeeeee")
    ax.set_xlim(-ulim, ulim); ax.set_ylim(-Flim, Flim)
    ax.set_xlabel("OOP rel. displacement at 2.06 m (mm)"); ax.set_title(label, fontsize=10.5, fontweight="bold", loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0][0].set_ylabel("base shear (kN)")
fig.suptitle("Last run (FR76 x2.0, from the end-of-run-24 state): record vs routes   [shear = {}]".format(shear_lbl),
             fontsize=11, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.95)); plt.savefig(os.path.join(OUT_DIR, "fig_lastrun_loops.png")); print("-> fig_lastrun_loops.png")

# ============================ FIG 2: overlay ==========================
fig, ax = plt.subplots(figsize=(7.5, 6), dpi=140)
if exp is not None:
    ax.plot(exp[0], exp[1], lw=0.5, color="#999999", alpha=0.7, label="US-1 run {} (test, reference)".format(EXP_RUN))
if ctrl:
    ax.plot(ctrl[0], ctrl[1], lw=0.9, color=COL_CTRL, label="record, run {}".format(LAST_RUN))
for label, colr, L in series:
    ax.plot(L[0], L[1], lw=0.9, color=colr, label=label, alpha=0.9)
ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.grid(color="#eeeeee")
ax.set_xlim(-ulim, ulim); ax.set_ylim(-Flim, Flim)
ax.set_xlabel("OOP rel. displacement at 2.06 m (mm)"); ax.set_ylabel("base shear (kN)")
ax.set_title("Last run: base shear vs displacement, all series   [shear = {}]".format(shear_lbl), fontsize=10.5, fontweight="bold", loc="left")
ax.legend(fontsize=8.5, frameon=False, loc="upper left")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "fig_lastrun_overlay.png")); print("-> fig_lastrun_overlay.png")
