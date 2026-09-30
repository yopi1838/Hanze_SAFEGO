# -*- coding: ascii -*-
"""
postprocess_routes.py -- Route 1 / Route 2 / Strategy-F control vs US-1.

Simulation metrics per run (definitions copied from postprocess_stratC.py):
    peak OOP rel. disp  : max |u(t) - u(t0)|                      [mm]
                          EDP = rel_disp_top_exp_mm (Eq. 1 of Moshfeghi et al.)
                          else rel_disp_top_mm, else 0.5(Ch3+Ch4)-Ch5 rebuilt
                          from the exported CSVs (extract_topq_avg_to_table).
    peak tilt           : max |tilt(t) - tilt(t0)|               [deg]
    residual tilt       : mean(last 5%) - tilt(t0 of Run 1)      [deg]
                          (compute_metrics: run1_start per channel; a one-shot
                          route folder has only itself as "Run 1")
    peak base shear     : max |F(t) - F(t0)|                      [kN]
                          F = cstav + topj_shear when the run exported both
                          (total base shear, Fig 13 quantity); else cstav alone
                          (base + joist reaction only, ~55% of the total).

Experimental data (US-1 = Test 9), read from the processed workbooks:
    EXP_DATA/processed_globalzero/Test9RunNN_processed_globalzero.xlsx
        Time | U_avg [m] (mean(Ch3,Ch4) - table, zeroed at run-1 start) |
        base_shear_kN = sum(m_i a_i), i = 15/16/17 | mass_scale | m15..m17_kg
      -> peak OOP rel. disp = max|U_avg - U_avg(t0)| * 1000 (= 29.48 mm at run 24)
      -> peak base shear    = max|F - F(t0)|            (= 29.50 kN at run 24)
      -> run-24 loop (u - u(t0), EXP_SHEAR_SIGN * F) for the shear-displacement panels
    exp_Test9_tilt.csv (optional, postprocess_stratC convention) -> peak / residual
      tilt; the workbooks carry no tilt channel, so without it the tilt panels
      show simulation only.
    US1_fig13_digitised.csv (optional) -> thin envelope overlay, sign flipped
      (EXP_FIG13_SIGN) into the raw-data / model convention.
Figure 2 (1x3): US-1 run 24 loop vs control runs 24/25 | vs Route 1 | vs Route 2.

Usage: python postprocess_routes.py
"""
import os, re, csv, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import safego_paths as sp
except ImportError:
    sp = None

# ============================ CONFIG =================================
CASES = [   # (label, results dir, style)
    ("Control (Strategy F)", "stratF_full_results_US1", dict(color="#5b7f74", marker="s")),
    ("Route 2 (two-stage)",  "route2_results",          dict(color="#b85042", marker="o")),
    ("Route 1 (velocity)",   "route1_results",          dict(color="#d98c21", marker="^")),
]
MODE = "single"            # "single": bars at SINGLE_SHOT_RUN; "sequence": metrics vs run
SINGLE_SHOT_RUN = 24       # route1_* / route2_* one-shot folders are compared here
CONTROL_RUNS = (24, 25)    # Strategy F runs shown beside the one-shot routes
HYST_RUNS = (24,)
OUT_PNG_M = "fig_routes_metrics.png"
OUT_PNG_H = "fig_routes_hysteresis.png"
OUT_CSV   = "routes_metrics.csv"

# --- conventions from postprocess_stratC.py (verbatim names) -----------
RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")
ROUTE_RE = re.compile(r"^route\d_[A-Za-z0-9]+_s\d+p\d+")
KEY_TILT = "tilt_full_wall"
KEY_DISP = "rel_disp_top_mm"
KEY_DISP_EXP = "rel_disp_top_exp_mm"
KEY_SHEAR = "cstav"             # base + joist reaction (dynamic exports)
KEY_SHEAR_TOP = "topj_shear"    # top-joint (T_B) reaction; added to cstav when the run exported it
EXP_FIG13 = "US1_fig13_digitised.csv"   # digitised Fig 13 envelope (drift_pct, F_kN, d_mm); same search order as exp_csv
EXP_FIG13_SIGN = -1.0   # Fig 13 uses the paper's sign; raw data / simulation EDP are opposite
                        # (postprocess_stratC: raw +29.48/-25.52 mm vs Table 4 25/-30 mm)
CH3_NAME = "Channel_3_DispTopQLeft"
CH4_NAME = "Channel_4_DispTopQRight"
CH_TABLE = "Channel_5_DispTable"
TAIL_FRAC = 0.05
# --- experimental data: EXP_DATA/processed_globalzero/Test9RunNN_processed_globalzero.xlsx
# columns: Time [s] | U_avg [m] = mean(Ch3, Ch4) - table, zeroed at the start of RUN 1
#          ("globalzero") | base_shear_N | base_shear_kN = sum(m_i * a_i) over the
# accelerometers 15/16/17 (bot quarter 0.66 m, mid 1.26 m, top quarter 2.06 m)
# | mass_scale | m15_kg | m16_kg | m17_kg  (constant per run; m_i already scaled)
EXP_ROOT   = "EXP_DATA"                 # also tried relative to cwd and to this script
EXP_GZ_DIR = "processed_globalzero"
EXP_TEST   = "Test9"                    # US-1
EXP_LABEL  = "US-1 (Test 9)"
EXP_COLOR  = "royalblue"
EXP_COL_T, EXP_COL_U, EXP_COL_F = "Time", "U_avg", "base_shear_kN"
EXP_SHEAR_SIGN = -1.0   # base_shear = sum(m a) is the inertia force: +F at -u (run 24: F max 29.5 kN
                        # at u = -11 mm, corr(u,F) = -0.59). Fig 13 draws base shear +F at +d, so the
                        # loop is negated to plot in the Fig 13 sense. Set +1 to plot as stored.
EXP_TILT_CSV = "exp_Test9_tilt.csv"     # OPTIONAL: tilt is not in the workbooks; read via exp_csv()
                                        # as postprocess_stratC.py does, columns TILT_SEG_MAP
TILT_SEG_MAP = {"tilt_full_wall": ("peak_tilt_full_wall", "resid_tilt_full_wall")}
SHOW_FIG13 = True                       # thin Fig 13 envelope overlay on the loop panels

# ============================ READERS (postprocess_stratC.py) ==========
def read_hist_csv(fpath):
    """3DEC history export: 2 header lines, whitespace-separated."""
    try:
        d = np.genfromtxt(str(fpath), skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2:
        return None, None
    m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    d = d[m]
    if len(d) < 5:
        return None, None
    return d[:, 0], d[:, 1]

def find_channel(folder, key):
    hits = sorted(Path(folder).glob("*" + key + "*.csv"))
    return read_hist_csv(hits[0]) if hits else (None, None)

def discover(results_dir):
    """{run_no: Path}. RunNN_REC_sXpYY -> NN; a routeN_REC_sXpYY one-shot folder -> SINGLE_SHOT_RUN."""
    out = {}
    if not Path(results_dir).is_dir():
        return out
    for d in sorted(Path(results_dir).iterdir()):
        if not d.is_dir():
            continue
        m = RUN_RE.match(d.name)
        if m:
            out[int(m.group(1))] = d
        elif ROUTE_RE.match(d.name):
            out[SINGLE_SHOT_RUN] = d
    return out

def edp_series(folder):
    """(t, u[mm], source): rel_disp_top_exp_mm > rel_disp_top_mm > 0.5(Ch3+Ch4)-Ch5."""
    for key in (KEY_DISP_EXP, KEY_DISP):
        t, u = find_channel(folder, key)
        if u is not None:
            return t, u, key
    t3, c3 = find_channel(folder, CH3_NAME); _, c4 = find_channel(folder, CH4_NAME)
    _, c5 = find_channel(folder, CH_TABLE)
    if c3 is None or c4 is None or c5 is None:
        return None, None, None
    n = min(len(c3), len(c4), len(c5))
    return t3[:n], (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0, "0.5(Ch3+Ch4)-Ch5"

def shear_series(folder):
    """(t, F[kN], label). cstav + topj_shear when both exist (total base shear,
    comparable to Fig 13); otherwise cstav alone (base + joist only, ~55% of total)."""
    t, F = find_channel(folder, KEY_SHEAR)
    if F is None:
        return None, None, None
    _, Ft = find_channel(folder, KEY_SHEAR_TOP)
    if Ft is not None:
        n = min(len(F), len(Ft))
        return t[:n], F[:n] + Ft[:n], "cstav+topj_shear (total)"
    return t, F, "cstav (base+joist only)"

def compute_metrics(t, x, ch, run1_start):
    """Verbatim from postprocess_stratC.py."""
    x0 = float(x[0])
    if ch not in run1_start:
        run1_start[ch] = x0
    tail = x[int((1.0 - TAIL_FRAC) * len(x)):]
    end_mean = float(np.mean(tail))
    return {"start": x0, "peak": float(np.max(np.abs(x - x0))),
            "end_mean": end_mean, "residual": end_mean - run1_start[ch]}

def case_metrics(run_folders):
    """{run: {peak_disp_mm, max_tilt_deg, resid_tilt_deg, peak_shear_kN}} in run order,
    so run1_start is the first run of THIS case (control: Run 1; one-shot route: itself)."""
    out, run1_start, edp_src, shear_src = {}, {}, set(), set()
    for rn in sorted(run_folders):
        f = run_folders[rn]; m = {}
        t, u, src = edp_series(f)
        if u is not None:
            m["peak_disp_mm"] = compute_metrics(t, u, "edp", run1_start)["peak"]; edp_src.add(src)
        t, th = find_channel(f, KEY_TILT)
        if th is not None:
            c = compute_metrics(t, th, KEY_TILT, run1_start)
            m["max_tilt_deg"], m["resid_tilt_deg"] = c["peak"], c["residual"]
        t, F, src = shear_series(f)
        if F is not None:
            m["peak_shear_kN"] = compute_metrics(t, F, KEY_SHEAR, run1_start)["peak"]; shear_src.add(src)
        out[rn] = m
    return out, edp_src, shear_src

def hysteresis(folder):
    """(u - u(t0) [mm], F - F(t0) [kN]) -- both zeroed at the start of the run."""
    _, u, _ = edp_series(folder); _, F, _ = shear_series(folder)
    if u is None or F is None:
        return None
    n = min(len(u), len(F))
    return u[:n] - u[0], F[:n] - F[0]

def load_fig13(out_dir):
    """Digitised Fig 13 envelope of US-1 in the raw-data / simulation sign convention."""
    p = exp_csv(EXP_FIG13, out_dir)
    if p is None:
        print("  ! {} not found -- loop panels will show simulation only".format(EXP_FIG13)); return None
    d = np.genfromtxt(str(p), delimiter=",", names=True)
    if "d_mm" not in d.dtype.names or "F_kN" not in d.dtype.names:
        print("  ! {} lacks d_mm / F_kN columns -- skipped".format(EXP_FIG13)); return None
    print("  experiment 'US-1 Fig 13 envelope': {} ({} points, sign x{:+.0f})".format(p, len(d), EXP_FIG13_SIGN))
    return EXP_FIG13_SIGN * d["d_mm"], EXP_FIG13_SIGN * d["F_kN"]

# ============================ EXPERIMENT (postprocess_stratC.py) =======
def exp_csv(name, out_dir):
    local = Path(out_dir) / name
    if local.is_file():
        return local
    for cand in (Path.cwd() / name, Path(__file__).resolve().parent / name):
        if cand.is_file():
            return cand
    if sp is not None:
        try:
            p = Path(sp.exp_derived(name))
            if p.is_file():
                return p
        except Exception:
            pass
    return None

def read_exp_csv(path):
    """-> {column_name: {run_no: value}}, skipping blanks and non-numerics."""
    out = {}
    with open(str(path)) as f:
        for row in csv.DictReader(f):
            try:
                rn = int(float(row["run"]))
            except (TypeError, ValueError, KeyError):
                continue
            for k, v in row.items():
                if k == "run" or v in ("", None):
                    continue
                try:
                    out.setdefault(k, {})[rn] = float(v)
                except (TypeError, ValueError):
                    pass
    return out

def exp_gz_path(run):
    """Locate Test9RunNN_processed_globalzero.xlsx (NN zero-padded or not)."""
    roots = [Path(EXP_ROOT), Path.cwd() / EXP_ROOT, Path(__file__).resolve().parent / EXP_ROOT]
    if sp is not None:
        for attr in ("exp_data", "EXP_DATA"):
            try:
                roots.append(Path(getattr(sp, attr)() if callable(getattr(sp, attr)) else getattr(sp, attr)))
            except Exception:
                pass
    for root in roots:
        for nm in ("{}Run{:02d}_processed_globalzero.xlsx", "{}Run{}_processed_globalzero.xlsx"):
            for base in (root / EXP_GZ_DIR, root):
                p = base / nm.format(EXP_TEST, run)
                if p.is_file():
                    return p
    return None

def read_gz(path):
    """-> {column: np.array} from the single-sheet workbook (header row = column names)."""
    try:
        import openpyxl
    except ImportError:
        print("  ! openpyxl not installed -- pip install openpyxl to read the experimental workbooks"); return None
    ws = openpyxl.load_workbook(str(path), read_only=True, data_only=True).worksheets[0]
    rows = ws.iter_rows(values_only=True)
    hdr = [str(h) for h in next(rows)]
    data = np.array([r for r in rows if r and r[0] is not None], dtype=float)
    return {h: data[:, j] for j, h in enumerate(hdr)}

_GZ_CACHE = {}
def exp_run(run):
    """(t, u[mm], F[kN], info) for one experimental run, or None. u is globalzero (run-1 start)."""
    if run in _GZ_CACHE:
        return _GZ_CACHE[run]
    p = exp_gz_path(run)
    if p is None:
        _GZ_CACHE[run] = None; return None
    d = read_gz(p)
    if d is None or any(c not in d for c in (EXP_COL_T, EXP_COL_U, EXP_COL_F)):
        print("  ! {}: missing one of {} -- columns are {}".format(p.name, (EXP_COL_T, EXP_COL_U, EXP_COL_F), list(d or {})))
        _GZ_CACHE[run] = None; return None
    info = {k: float(d[k][0]) for k in ("mass_scale", "m15_kg", "m16_kg", "m17_kg") if k in d}
    _GZ_CACHE[run] = (d[EXP_COL_T], d[EXP_COL_U] * 1000.0, d[EXP_COL_F], info, p)
    return _GZ_CACHE[run]

def exp_metrics(run):
    """Same definitions as compute_metrics: peak = max|x - x(t0)| within the run,
    residual = mean(last 5%) - x(t0 of Run 1) = tail mean, since the workbooks are globalzero."""
    r = exp_run(run)
    if r is None:
        return {}
    t, u, F, info, p = r
    m = {"peak_disp_mm": float(np.max(np.abs(u - u[0]))),
         "resid_disp_mm": float(np.mean(u[int((1.0 - TAIL_FRAC) * len(u)):])),
         "peak_shear_kN": float(np.max(np.abs(F - F[0])))}
    m.update(info)
    return m

def exp_loop(run):
    """(u - u(t0) [mm], sign * (F - F(t0)) [kN]) -- zeroed at the start of the run like hysteresis()."""
    r = exp_run(run)
    if r is None:
        return None
    t, u, F, info, p = r
    return u - u[0], EXP_SHEAR_SIGN * (F - F[0])

def exp_tilt_table(out_dir):
    """Optional {run: {max_tilt_deg, resid_tilt_deg}} from exp_Test9_tilt.csv (postprocess_stratC convention)."""
    p = exp_csv(EXP_TILT_CSV, out_dir)
    if p is None:
        print("  ! {} not found -- tilt panels show simulation only (tilt is not in the workbooks)".format(EXP_TILT_CSV))
        return {}
    tab = read_exp_csv(p); pk_col, rs_col = TILT_SEG_MAP[KEY_TILT]; out = {}
    for key, col in (("max_tilt_deg", pk_col), ("resid_tilt_deg", rs_col)):
        if col not in tab:
            print("  ! {} has no column '{}'".format(EXP_TILT_CSV, col)); continue
        for rn, v in tab[col].items():
            out.setdefault(rn, {})[key] = v
    print("  experiment tilt: {}".format(p))
    return out

# ============================ COLLECT ================================
table, folders, SHEAR_SRC = {}, {}, set()
print("discovery:")
for label, d, _ in CASES:
    if not Path(d).is_dir():
        print("  {:22s} dir '{}' NOT FOUND -- fix CASES".format(label, d)); folders[label] = {}; continue
    folders[label] = discover(d)
    mets, edp_src, shear_src = case_metrics(folders[label])
    SHEAR_SRC |= shear_src
    missing = {"disp": 0, "tilt": 0, "shear": 0}
    for rn, m in mets.items():
        table[(label, rn)] = m
        for k, key in (("disp", "peak_disp_mm"), ("tilt", "max_tilt_deg"), ("shear", "peak_shear_kN")):
            if key not in m: missing[k] += 1
    print("  {:22s} {:2d} run folders in '{}'  EDP={}  shear={}  (missing: disp {} tilt {} shear {})".format(
        label, len(mets), d, "/".join(sorted(edp_src)) or "-", "/".join(sorted(shear_src)) or "-",
        missing["disp"], missing["tilt"], missing["shear"]))
SHEAR_LBL = " | ".join(sorted(SHEAR_SRC)) or KEY_SHEAR

ctrl_dir = next((d for l, d, _ in CASES if l.lower().startswith("control")), ".")
EXP_RUNS = range(1, 26) if MODE == "sequence" else (SINGLE_SHOT_RUN,)
exp = {}
for rn in EXP_RUNS:
    m = exp_metrics(rn)
    if m:
        exp[rn] = m
if exp:
    ex = exp_run(SINGLE_SHOT_RUN if SINGLE_SHOT_RUN in exp else sorted(exp)[0])
    print("  experiment '{}': {} run workbook(s) under '{}'  e.g. {}  (mass_scale {:.3f}, m15+m16+m17 = {:.0f} kg)".format(
        EXP_LABEL, len(exp), ex[4].parent, ex[4].name, ex[3].get("mass_scale", float("nan")),
        sum(ex[3].get(k, 0.0) for k in ("m15_kg", "m16_kg", "m17_kg"))))
else:
    print("  ! no {}RunNN_processed_globalzero.xlsx found under '{}/{}' -- figures show simulation only".format(
        EXP_TEST, EXP_ROOT, EXP_GZ_DIR))
for rn, m in exp_tilt_table(Path(ctrl_dir) / "postproc").items():
    if rn in exp or MODE == "sequence":
        exp.setdefault(rn, {}).update(m)
fig13 = load_fig13(Path(ctrl_dir) / "postproc") if SHOW_FIG13 else None
EXP_LOOP = exp_loop(SINGLE_SHOT_RUN)

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "run", "peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg"])
    for rn in sorted(exp):
        w.writerow([EXP_LABEL, rn] + [round(exp[rn][k], 4) if k in exp[rn] else "" for k in
                                             ("peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg")])
    for (label, rn), m in sorted(table.items(), key=lambda x: (x[0][0], x[0][1])):
        w.writerow([label, rn] + [round(m[k], 4) if k in m else "" for k in
                                  ("peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg")])
print("-> " + OUT_CSV)

PANELS = [("peak_disp_mm", "Absolute max OOP rel. displacement (mm)"),
          ("peak_shear_kN", "Peak base shear (kN)"),
          ("max_tilt_deg", "Peak tilt, full wall (deg)"),
          ("resid_tilt_deg", "Residual tilt, full wall (deg)")]

# ============================ FIG 1 (single): bars at SINGLE_SHOT_RUN ==
if MODE == "single":
    ref = exp.get(SINGLE_SHOT_RUN, {})
    fig, axs = plt.subplots(1, 4, figsize=(18.5, 4.8), dpi=140)
    for ax, (key, ttl) in zip(axs, PANELS):
        names, vals, cols = [], [], []
        if key in ref:
            names.append("US-1\nrun {}".format(SINGLE_SHOT_RUN)); vals.append(ref[key]); cols.append(EXP_COLOR)
        for label, _, st in CASES:
            is_ctrl = label.lower().startswith("control")
            for rn in (CONTROL_RUNS if is_ctrl else (SINGLE_SHOT_RUN,)):
                m = table.get((label, rn))
                if m and key in m:
                    nm = ("Control\nrun {}".format(rn) if is_ctrl else label.replace(" (", "\n("))
                    names.append(nm); vals.append(m[key])
                    cols.append(st["color"] if not (is_ctrl and rn != CONTROL_RUNS[-1]) else "#9fb3ad")
        if not vals:
            ax.set_title(ttl + "  (no data)", fontsize=10, loc="left"); continue
        bars = ax.bar(range(len(vals)), vals, color=cols, width=0.62)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, "{:.2f}".format(v), ha="center", va="bottom", fontsize=9)
        if key in ref:
            ax.axhline(ref[key], color=EXP_COLOR, lw=0.9, ls=":")
        ax.set_xticks(range(len(vals))); ax.set_xticklabels(names, fontsize=8)
        ax.set_title(ttl, fontsize=10.5, fontweight="bold", loc="left")
        ax.grid(axis="y", color="#e5e5e5")
        for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
    fig.suptitle("One-shot routes (FR76 x2.0) vs Strategy F control (runs {}) vs US-1 run {}   "
                 "(base shear = {})".format("/".join(map(str, CONTROL_RUNS)), SINGLE_SHOT_RUN, SHEAR_LBL),
                 fontsize=11, fontweight="bold")
    plt.tight_layout(); plt.savefig(OUT_PNG_M, dpi=140); print("-> " + OUT_PNG_M)

# ============================ FIG 1 (sequence): metrics vs run =========
if MODE == "sequence":
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5), dpi=140)
    for ax, (key, ttl) in zip(axs.ravel(), PANELS):
        xs = sorted(r for r in exp if key in exp[r])
        if xs:
            ax.plot(xs, [exp[r][key] for r in xs], marker="o", ls="--",
                    color=EXP_COLOR, mfc="none", ms=5, lw=1.0, label=EXP_LABEL)
        for label, _, st in CASES:
            xs = sorted(r for (l, r) in table if l == label and key in table[(l, r)])
            if not xs:
                continue
            ys = [table[(label, r)][key] for r in xs]
            if len(xs) == 1:
                ax.plot(xs, ys, linestyle="none", ms=11, mew=2, label=label, **st)
            else:
                ax.plot(xs, ys, "-", lw=1.6, ms=5, label=label, **st)
        ax.axvline(21.5, color="#bbb", lw=0.8); ax.set_xlim(0, 26)
        ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left")
        ax.set_xlabel("run"); ax.grid(axis="y", color="#e5e5e5")
        for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
    axs[0, 0].legend(fontsize=9, frameon=False)
    fig.suptitle("Route 1 / Route 2 / Control vs US-1, run by run  (base shear = {})".format(SHEAR_LBL),
                 fontsize=12, fontweight="bold")
    plt.tight_layout(); plt.savefig(OUT_PNG_M, dpi=140); print("-> " + OUT_PNG_M)

# ============================ FIG 2: base shear vs OOP displacement, 1x3 =
# col 1: US-1 Fig 13 envelope + control runs CONTROL_RUNS (run 24 dotted, run 25 solid)
# col 2: US-1 + Route 1 one-shot;  col 3: US-1 + Route 2 one-shot
ctrl = next((c for c in CASES if c[0].lower().startswith("control")), None)
routes = sorted((c for c in CASES if not c[0].lower().startswith("control")), key=lambda c: c[0])  # Route 1, then Route 2
cols = ([("Control", [(ctrl, r_, ls_) for r_, ls_ in zip(CONTROL_RUNS, (":", "-"))])] if ctrl else []) + \
       [(c[0], [(c, SINGLE_SHOT_RUN, "-")]) for c in routes]
fig, axs = plt.subplots(1, len(cols), figsize=(5.4 * len(cols), 5.0), dpi=140, sharex=True, sharey=True)
axs = np.atleast_1d(axs)
for ax, (ttl, series) in zip(axs, cols):
    if fig13 is not None:
        ax.plot(fig13[0], fig13[1], "--", color="#8aa4d6", lw=1.2, alpha=0.9,
                label="US-1 Fig 13 envelope (digitised)", zorder=3)
    if EXP_LOOP is not None:
        ax.plot(EXP_LOOP[0], EXP_LOOP[1], "-", color=EXP_COLOR, lw=0.9, alpha=0.95,
                label="{} run {} (processed_globalzero)".format(EXP_LABEL, SINGLE_SHOT_RUN), zorder=4)
    for (label, _, st), r_, ls_ in series:
        f = folders.get(label, {}).get(r_)
        h = hysteresis(f) if f else None
        if h is None:
            ax.text(0.5, 0.5, "{} run {}: no data".format(label, r_), transform=ax.transAxes,
                    ha="center", fontsize=9, color="#888"); continue
        lab = "{} run {}".format(label, r_) if series is cols[0][1] and ctrl else "{} (one shot, FR76 x2.0)".format(label)
        ax.plot(h[0], h[1], lw=0.8, ls=ls_, color=st["color"], label=lab, alpha=0.9, zorder=2)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_title("US-1 vs {}".format(ttl), fontsize=11, fontweight="bold", loc="left")
    ax.set_xlabel("OOP rel. displacement at Z = 2.06 m (mm)")
    ax.grid(color="#eeeeee"); ax.legend(fontsize=8, frameon=False, loc="upper left")
    for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
axs[0].set_ylabel("base shear (kN)")
fig.suptitle("Base shear vs OOP displacement at Z = 2.06 m   (simulation shear = {}; US-1 shear = sum(m a) x{:+.0f})".format(
    SHEAR_LBL, EXP_SHEAR_SIGN), fontsize=11, fontweight="bold")
plt.tight_layout(); plt.savefig(OUT_PNG_H, dpi=140); print("-> " + OUT_PNG_H)