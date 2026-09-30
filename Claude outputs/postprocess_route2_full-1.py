# -*- coding: ascii -*-
"""
postprocess_route2_full.py -- Route 2 over the FULL 25-run protocol
(route2_full25_3dec.py, or any driver that writes RunNN_REC_sXpYY folders
plus a per-run log) against the Strategy F control (record) and US-1.

Three series, run by run
    US-1        : EXP_DATA/processed_globalzero/Test9RunNN_processed_globalzero.xlsx
                  (U_avg [m] -> mm, base_shear_kN = sum(m a)); runs 1-24 only,
                  run 25 was not applied to US-1 (Table 2 asterisk).
    control     : CONTROL_DIR/RunNN_*   (the record on the model, Strategy F)
    Route 2 full: ROUTE_DIR/RunNN_*     (the two-stage pulses on the model)

Definitions, identical to postprocess_stratC.py (all WITHIN-RUN):
    peak OOP rel. disp  = max |u - u(t0)|, u = rel_disp_top_exp_mm
                          (else rel_disp_top_mm, else 0.5(Ch3+Ch4)-Ch5)   [mm]
    peak base shear     = max |F - F(t0)|, F = cstav (+ topj_shear if exported) [kN]
    peak tilt           = max |tilt - tilt(t0)|                            [deg]
    residual tilt       = mean(last 5%) - tilt(t0 of Run 1)                [deg]
NOTE the driver's own log column peak_edp_mm is computed from Record_Disp
(the integrated input table), not from Ch5; it is carried into the table as
log_peak_edp_mm for reference but every figure uses the recomputed EDP.

The per-run log (ROUTE_LOG) supplies the stage data: T_A_s, T_B_s, PGA_A_g,
PGA_B_g, two_stage, d_hist_mm, SdA_mm, SdB_mm. Missing columns are skipped.

Outputs (OUT_DIR = route2_full_postproc/)
    route2_full_metrics.csv      one row per run, three series side by side + stage data
    fig_r2full_sequence.png      1x3: peak disp / peak shear / tilt vs run, all three series.
                                 Tilt: if exp_Test9_tilt.csv exists (made by exp_tilt_from_raw.py from
                                 the raw Test9RunNN.xlsx channels 1-4), the panel compares the model's
                                 tilt_full_wall channel with US-1's reconstructed full-wall tilt, peak
                                 and residual. Otherwise it falls back to a proxy computed the same way
                                 for all three: atan(u / H_EDP_M) from the 2.06 m displacement.
                                 The stage data (periods, PGA, Sd targets) are in the CSV only.
    fig_r2full_loops.png         base shear vs displacement for HYST_RUNS (US-1 | control | Route 2)
    fig_r2full_histories.png     displacement histories for HIST_RUNS (control and Route 2)
Usage: python postprocess_route2_full.py
"""
import os, re, csv, glob
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
ROUTE_DIR   = "route2_full25_revised"     # or "route2_full25" for the old driver
ROUTE_LOG   = os.path.join(ROUTE_DIR, "route2_full_log.csv")
ROUTE_LABEL = "Route 2 (full sequence)"
CONTROL_DIR = "stratF_full_results_US1"
CONTROL_LABEL = "record (Strategy F)"
EXP_ROOT, EXP_GZ_DIR, EXP_TEST = "EXP_DATA", "processed_globalzero", "Test9"
EXP_LABEL   = "US-1 (Test 9)"
EXP_RUNS    = range(1, 25)          # run 25 not applied to US-1
EXP_SHEAR_SIGN = -1.0               # sum(m a) is the inertia force; flipped to the Fig 13 sense
MODEL_SHEAR_SIGN = -1.0      # cstav (+topj_shear) is the support reaction; flipped so +F goes with +u like the US-1 loops (Fig 13 sense)
H_EDP_M     = 2.06                  # height of the displacement channel: tilt proxy = atan(u / H_EDP_M)
EXP_TILT_CSV = "exp_Test9_tilt.csv" # OPTIONAL measured US-1 tilt: columns run, peak_tilt_full_wall, resid_tilt_full_wall (deg)
HYST_RUNS   = (21, 22, 23, 24, 25)  # loop panels
HIST_RUNS   = (22, 23, 24, 25)      # displacement-history panels
OUT_DIR     = "route2_full_postproc"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_CSV     = os.path.join(OUT_DIR, "route2_full_metrics.csv")
OUT_SEQ     = os.path.join(OUT_DIR, "fig_r2full_sequence.png")
OUT_LOOPS   = os.path.join(OUT_DIR, "fig_r2full_loops.png")
OUT_HIST    = os.path.join(OUT_DIR, "fig_r2full_histories.png")
TAIL_FRAC   = 0.05
RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")
KEY_DISP_EXP, KEY_DISP, KEY_TILT, KEY_SHEAR, KEY_SHEAR_TOP = \
    "rel_disp_top_exp_mm", "rel_disp_top_mm", "tilt_full_wall", "cstav", "topj_shear"
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
COL = {"exp": "royalblue", "ctrl": "#5b7f74", "r2": "#b85042", "A": "#5b7f74", "B": "#b85042"}
PROTOCOL = {1: ("HU12", 0.50), 2: ("HU12", 0.75), 3: ("EC40", 0.20), 4: ("HU12", 1.00), 5: ("HU12", 1.25),
            6: ("EC40", 0.30), 7: ("HU12", 1.50), 8: ("EC40", 0.40), 9: ("HU12", 1.75), 10: ("HU12", 2.00),
            11: ("EC40", 0.50), 12: ("HU12", 2.25), 13: ("HU12", 2.50), 14: ("HU12", 2.75), 15: ("HU12", 3.00),
            16: ("HU12", 3.50), 17: ("HU12", 4.00), 18: ("HU12", 4.50), 19: ("HU12", 5.00), 20: ("HU12", 5.50),
            21: ("HU12", 6.00), 22: ("FR76", 1.00), 23: ("FR76", 1.50), 24: ("FR76", 1.75), 25: ("FR76", 2.00)}

# ============================ READERS (postprocess_stratC conventions) ==
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

def discover(results_dir):
    out = {}
    if not Path(results_dir).is_dir():
        return out
    for d in sorted(Path(results_dir).iterdir()):
        m = RUN_RE.match(d.name)
        if d.is_dir() and m:
            out[int(m.group(1))] = d
    return out

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

def case_metrics(run_folders):
    """{run: metrics + series}; residual tilt relative to the start of the first run of THIS case."""
    out, tilt0, u0, edp_src, shear_src = {}, None, None, set(), set()
    for rn in sorted(run_folders):
        f = run_folders[rn]; m = {}
        t, u, src = edp_series(f)
        if u is not None:
            if u0 is None:
                u0 = float(u[0])                      # the channel is absolute: run 1's start is the global zero
            pk = float(np.max(np.abs(u - u[0]))); rs = float(np.mean(u[int((1 - TAIL_FRAC) * len(u)):]) - u0)
            m.update(peak_disp_mm=pk, resid_disp_mm=rs, t_u=t, u=u - u[0],
                     peak_tiltp_deg=np.degrees(np.arctan(pk / 1000.0 / H_EDP_M)), resid_tiltp_deg=np.degrees(np.arctan(rs / 1000.0 / H_EDP_M)))
            edp_src.add(src)
        t, th = find_channel(f, KEY_TILT)
        if th is not None:
            if tilt0 is None:
                tilt0 = float(th[0])
            m.update(max_tilt_deg=float(np.max(np.abs(th - th[0]))),
                     resid_tilt_deg=float(np.mean(th[int((1 - TAIL_FRAC) * len(th)):]) - tilt0))
        t, F, src = shear_series(f)
        if F is not None:
            m.update(peak_shear_kN=float(np.max(np.abs(F - F[0]))), t_F=t, F=F - F[0]); shear_src.add(src)
        out[rn] = m
    return out, edp_src, shear_src

# ============================ EXPERIMENT ============================
def exp_gz_path(run):
    for root in (Path(EXP_ROOT), Path.cwd() / EXP_ROOT, Path(__file__).resolve().parent / EXP_ROOT):
        for nm in ("{}Run{:02d}_processed_globalzero.xlsx", "{}Run{}_processed_globalzero.xlsx"):
            for base in (root / EXP_GZ_DIR, root):
                p = base / nm.format(EXP_TEST, run)
                if p.is_file():
                    return p
    return None

def exp_run(run):
    p = exp_gz_path(run)
    if p is None:
        return None
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
    pk = float(np.max(np.abs(u - u[0]))); rs = float(np.mean(u[int((1 - TAIL_FRAC) * len(u)):]))   # globalzero: tail mean = residual from run 1
    return dict(peak_disp_mm=pk, resid_disp_mm=rs, peak_shear_kN=float(np.max(np.abs(F - F[0]))),
                peak_tiltp_deg=np.degrees(np.arctan(pk / 1000.0 / H_EDP_M)), resid_tiltp_deg=np.degrees(np.arctan(rs / 1000.0 / H_EDP_M)),
                u=u - u[0], F=EXP_SHEAR_SIGN * (F - F[0]), path=p)

def exp_tilt_table():
    """optional measured US-1 tilt per run -> {run: {max_tilt_deg, resid_tilt_deg}}"""
    for cand in (Path(EXP_TILT_CSV), Path(EXP_ROOT) / EXP_TILT_CSV, Path(__file__).resolve().parent / EXP_TILT_CSV):
        if cand.is_file():
            out = {}
            for row in csv.DictReader(open(cand)):
                try:
                    rn = int(float(row["run"]))
                except Exception:
                    continue
                for key, colname in (("max_tilt_deg", "peak_tilt_full_wall"), ("resid_tilt_deg", "resid_tilt_full_wall")):
                    try:
                        out.setdefault(rn, {})[key] = float(row[colname])
                    except Exception:
                        pass
            print("  US-1 tilt: {} ({} runs)".format(cand, len(out))); return out
    print("  US-1 tilt: no {} -- run exp_tilt_from_raw.py on the raw Test9RunNN.xlsx to make it; the tilt panel falls back to the displacement proxy".format(EXP_TILT_CSV)); return {}

# ============================ COLLECT ================================
print("discovery:")
ctrl_f, r2_f = discover(CONTROL_DIR), discover(ROUTE_DIR)
ctrl, ctrl_edp, ctrl_sh = case_metrics(ctrl_f)
r2, r2_edp, r2_sh = case_metrics(r2_f)
print("  control  : {:2d} run folders in '{}'  EDP={}  shear={}".format(len(ctrl), CONTROL_DIR, "/".join(sorted(ctrl_edp)) or "-", "/".join(sorted(ctrl_sh)) or "-"))
print("  Route 2  : {:2d} run folders in '{}'  EDP={}  shear={}".format(len(r2), ROUTE_DIR, "/".join(sorted(r2_edp)) or "-", "/".join(sorted(r2_sh)) or "-"))
if not r2:
    raise SystemExit("no RunNN_* folders under '{}' -- fix ROUTE_DIR".format(ROUTE_DIR))
exp = {}
for rn in EXP_RUNS:
    e = exp_run(rn)
    if e:
        exp[rn] = e
print("  US-1     : {:2d} runs from {}".format(len(exp), exp[min(exp)]["path"].parent if exp else "(not found)"))
exp_tilt = exp_tilt_table()
for rn, d in exp_tilt.items():
    if rn in exp:
        exp[rn].update(d)

log = {}
if os.path.isfile(ROUTE_LOG):
    for row in csv.DictReader(open(ROUTE_LOG)):
        try:
            log[int(float(row["run"]))] = row
        except Exception:
            pass
    print("  log      : {} rows in {}  columns {}".format(len(log), ROUTE_LOG, list(next(iter(log.values())).keys()) if log else "-"))
else:
    print("  log      : {} not found -- stage panel will be empty".format(ROUTE_LOG))

def lf(rn, key):
    try:
        return float(log[rn][key])
    except Exception:
        return np.nan

# ============================ TABLE =================================
runs = sorted(set(r2) | set(ctrl) | set(exp))
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "record", "scale",
                "exp_peak_disp_mm", "exp_resid_disp_mm", "exp_peak_shear_kN", "exp_peak_tiltproxy_deg", "exp_resid_tiltproxy_deg", "exp_max_tilt_deg", "exp_resid_tilt_deg",
                "ctrl_peak_disp_mm", "ctrl_resid_disp_mm", "ctrl_peak_shear_kN", "ctrl_peak_tiltproxy_deg", "ctrl_resid_tiltproxy_deg", "ctrl_max_tilt_deg", "ctrl_resid_tilt_deg",
                "r2_peak_disp_mm", "r2_resid_disp_mm", "r2_peak_shear_kN", "r2_peak_tiltproxy_deg", "r2_resid_tiltproxy_deg", "r2_max_tilt_deg", "r2_resid_tilt_deg",
                "r2_over_ctrl_disp", "r2_over_exp_disp", "ctrl_over_exp_disp",
                "two_stage", "d_hist_mm", "T_A_s", "SdA_mm", "PGA_A_g", "T_B_s", "SdB_mm", "PGA_B_g", "log_peak_edp_mm", "cap_note"])
    for rn in runs:
        rec, sc = PROTOCOL.get(rn, ("?", ""))
        e, c, r = exp.get(rn, {}), ctrl.get(rn, {}), r2.get(rn, {})
        def g(d, k): return round(d[k], 4) if k in d else ""
        def ratio(a, b): return round(a[ "peak_disp_mm"] / b["peak_disp_mm"], 3) if "peak_disp_mm" in a and b.get("peak_disp_mm") else ""
        L = log.get(rn, {})
        w.writerow([rn, rec, sc] +
                   [g(d, k) for d in (e, c, r) for k in ("peak_disp_mm", "resid_disp_mm", "peak_shear_kN", "peak_tiltp_deg", "resid_tiltp_deg", "max_tilt_deg", "resid_tilt_deg")] +
                   [
                    ratio(r, c), ratio(r, e), ratio(c, e),
                    L.get("two_stage", ""), L.get("d_hist_mm", ""), L.get("T_A_s", ""), L.get("SdA_mm", ""), L.get("PGA_A_g", ""),
                    L.get("T_B_s", ""), L.get("SdB_mm", ""), L.get("PGA_B_g", ""), L.get("peak_edp_mm", ""), L.get("cap_note", "")])
print("-> " + OUT_CSV)

# ============================ FIG 1: sequence (validation row) ======
fig, axs = plt.subplots(1, 3, figsize=(19, 5.4), dpi=140)
def series(d, key):
    ks = [rn for rn in sorted(d) if key in d[rn]]; return ks, [d[rn][key] for rn in ks]
THREE = ((exp, EXP_LABEL, COL["exp"], "o"), (ctrl, CONTROL_LABEL, COL["ctrl"], "s"), (r2, ROUTE_LABEL, COL["r2"], "^"))
ax = axs[0]
for d, lab, c, mk in THREE:
    ks, vs = series(d, "peak_disp_mm")
    if ks: ax.plot(ks, vs, marker=mk, ms=5, lw=1.4, color=c, label=lab)
ax.set_ylabel("peak OOP rel. displacement at 2.06 m (mm)"); ax.set_title("Peak displacement per run", loc="left", fontweight="bold", fontsize=10.5)
ax.axhline(29.5, color=COL["exp"], lw=0.7, ls=":"); ax.text(1, 29.5, " 29.5 mm (US-1 run 24)", fontsize=8, va="bottom", color=COL["exp"])
ax = axs[1]
for d, lab, c, mk in THREE:
    ks, vs = series(d, "peak_shear_kN")
    if ks: ax.plot(ks, vs, marker=mk, ms=5, lw=1.4, color=c, label=lab)
ax.set_ylabel("peak base shear (kN)"); ax.set_title("Peak base shear per run  (model: {})".format("/".join(sorted(r2_sh | ctrl_sh)) or "-"), loc="left", fontweight="bold", fontsize=10.5)
ax = axs[2]
if exp_tilt:      # like for like: model tilt_full_wall channel vs US-1 tilt_full_wall reconstructed from channels 1-4 (exp_tilt_from_raw.py)
    for d, lab, c, mk in THREE:
        ks, vs = series(d, "max_tilt_deg")
        if ks: ax.plot(ks, vs, marker=mk, ms=5, lw=1.4, color=c, label=lab + ", peak")
        ks, vs = series(d, "resid_tilt_deg")
        if ks: ax.plot(ks, vs, marker=mk, ms=4, lw=1.0, ls="--", color=c, alpha=0.75, label=lab + ", residual (from run 1)")
    ax.set_ylabel("full-wall tilt (deg)")
    ax.set_title("Tilt per run: model tilt_full_wall vs US-1 (from channels 1-4)", loc="left", fontweight="bold", fontsize=10.5)
else:             # fallback: the same rigid-rocking proxy for all three
    for d, lab, c, mk in THREE:
        ks, vs = series(d, "peak_tiltp_deg")
        if ks: ax.plot(ks, vs, marker=mk, ms=5, lw=1.4, color=c, label=lab + ", peak")
        ks, vs = series(d, "resid_tiltp_deg")
        if ks: ax.plot(ks, vs, marker=mk, ms=4, lw=1.0, ls="--", color=c, alpha=0.75, label=lab + ", residual (from run 1)")
    ax.set_ylabel("tilt from the 2.06 m displacement, atan(u/{:.2f} m) (deg)".format(H_EDP_M))
    ax.set_title("Tilt per run: same proxy for all three (no exp_Test9_tilt.csv found)", loc="left", fontweight="bold", fontsize=10.5)
ax.axhline(0, color="k", lw=0.5)
for ax in axs:
    ax.set_xlabel("run"); ax.set_xticks(range(1, 26)); ax.tick_params(axis="x", labelsize=7.5); ax.grid(color="#eeeeee")
    ax.axvspan(21.5, 25.5, color="#f3e6d8", alpha=0.5, lw=0); ax.legend(fontsize=8, frameon=False)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0].text(23.5, 0.02, "FR76", ha="center", fontsize=8, color="#8a6a4a", transform=axs[0].get_xaxis_transform())
fig.suptitle("Route 2 full sequence vs the record on the model vs US-1  (within-run metrics; residuals from the run-1 zero)", fontsize=12, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.94)); plt.savefig(OUT_SEQ); plt.close(fig); print("-> " + OUT_SEQ)

# ============================ FIG 2: loops ==========================
hr = [rn for rn in HYST_RUNS if rn in r2 or rn in ctrl]
if hr:
    fig, axs = plt.subplots(1, len(hr), figsize=(5.2 * len(hr), 4.8), dpi=140, squeeze=False)
    for ax, rn in zip(axs[0], hr):
        rec, sc = PROTOCOL.get(rn, ("?", ""))
        if rn in exp:
            ax.plot(exp[rn]["u"], exp[rn]["F"], lw=0.6, color=COL["exp"], alpha=0.8, label=EXP_LABEL)
        for d, lab, c, ls in ((ctrl, CONTROL_LABEL, COL["ctrl"], ":"), (r2, ROUTE_LABEL, COL["r2"], "-")):
            m = d.get(rn, {})
            if "u" in m and "F" in m:
                n = min(len(m["u"]), len(m["F"])); ax.plot(m["u"][:n], m["F"][:n], lw=0.8, ls=ls, color=c, label="{} ({:.1f} mm)".format(lab, m["peak_disp_mm"]))
        ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.grid(color="#eee")
        ax.set_title("run {} ({} x{:g}){}".format(rn, rec, sc, "" if rn in exp or rn > 24 else "  [US-1 file missing]"), loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("OOP rel. displacement (mm)"); ax.legend(fontsize=7.5, frameon=False)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    axs[0][0].set_ylabel("base shear (kN)")
    fig.suptitle("Base shear vs displacement, last runs: US-1 | record on the model | Route 2 full sequence", fontsize=11.5, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.94)); plt.savefig(OUT_LOOPS); plt.close(fig); print("-> " + OUT_LOOPS)

# ============================ FIG 3: histories ======================
hh = [rn for rn in HIST_RUNS if rn in r2 or rn in ctrl]
if hh:
    fig, axs = plt.subplots(len(hh), 2, figsize=(14, 3.4 * len(hh)), dpi=140, squeeze=False)
    for i, rn in enumerate(hh):
        rec, sc = PROTOCOL.get(rn, ("?", ""))
        ylim = max([abs(v) for d in (ctrl, r2) if "u" in d.get(rn, {}) for v in (d[rn]["u"].min(), d[rn]["u"].max())] + [1.0]) * 1.1
        for ax, d, lab, c in ((axs[i, 0], ctrl, CONTROL_LABEL, COL["ctrl"]), (axs[i, 1], r2, ROUTE_LABEL, COL["r2"])):
            m = d.get(rn, {})
            if "u" in m:
                ax.plot(m["t_u"], m["u"], lw=0.7, color=c)
                ax.set_title("run {} ({} x{:g}) {}: peak {:.2f} mm".format(rn, rec, sc, lab, m["peak_disp_mm"]), loc="left", fontsize=9.5)
            else:
                ax.text(0.5, 0.5, "run {} not found".format(rn), transform=ax.transAxes, ha="center", color="#888")
            if d is r2 and rn in log:
                ax.text(0.98, 0.04, "A: {} s, {} g   B: {} s, {} g".format(log[rn].get("T_A_s", "?"), log[rn].get("PGA_A_g", "?"),
                        log[rn].get("T_B_s", "-") or "-", log[rn].get("PGA_B_g", "-") or "-"), transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#444")
            if rn in exp:
                ax.axhline(exp[rn]["peak_disp_mm"], color=COL["exp"], lw=0.7, ls=":"); ax.axhline(-exp[rn]["peak_disp_mm"], color=COL["exp"], lw=0.7, ls=":")
            ax.set_ylim(-ylim, ylim); ax.axhline(0, color="k", lw=0.5); ax.grid(color="#eee")
            ax.set_ylabel("rel. disp. (mm)"); ax.set_xlabel("t (s)")
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.suptitle("Displacement histories, FR76 runs: record (left) vs Route 2 (right); dotted = US-1 peak", fontsize=11.5, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.96)); plt.savefig(OUT_HIST); plt.close(fig); print("-> " + OUT_HIST)

# ============================ SUMMARY LINE ==========================
for rn in (21, 24, 25):
    e, c, r = exp.get(rn, {}).get("peak_disp_mm"), ctrl.get(rn, {}).get("peak_disp_mm"), r2.get(rn, {}).get("peak_disp_mm")
    print("run {:2d}: US-1 {}  record-on-model {}  Route 2 {}".format(
        rn, "{:.1f} mm".format(e) if e else "n/a", "{:.1f} mm".format(c) if c else "n/a", "{:.1f} mm".format(r) if r else "n/a"))
