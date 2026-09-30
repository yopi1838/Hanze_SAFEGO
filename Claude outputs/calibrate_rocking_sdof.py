# -*- coding: ascii -*-
"""
calibrate_rocking_sdof.py -- fit the rocking oscillator to the MODEL's own
record responses (Strategy F), so the amplitude rule is built on a reduced
model identified from the full model, not from the experiment.

What is fitted: m_eff, restitution R, viscous ratio xi (grid search), with
the backbone fixed (BACKBONE_CSV: the model pushover, else Fig 13). Objective:
mean |log(d_sdof / d_model)| over the runs in FIT_RUNS, using the record
tables the control was driven with (vel_HU / vel_EC / vel_FR x Table-2 scale)
and the control's exported peaks (rel_disp_top_exp_mm, else rel_disp_top_mm,
else 0.5(Ch3+Ch4)-Ch5).

Outputs (OUT_DIR): rocking_sdof_calibration.csv (all grid points),
rocking_sdof_best.json (the chosen set, read by strategy_R_3dec.py),
fig_rocking_calibration.png (per-run model vs oscillator for the best set).

Note on what to expect: the 3DEC model is undamped, so R ~ 1 and xi ~ 0 are
the physical end of the grid; with damping the oscillator falls to 2-4 mm on
every FR76 run. The oscillator is chaotic near the cliff: judge the fit on
the whole sequence, not on one run.
Usage: python calibrate_rocking_sdof.py
"""
import os, csv, json, glob, itertools
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rocking_sdof as rs

# ============================ CONFIG =================================
CONTROL_DIR  = "stratF_full_results_US1"
BACKBONE_CSV = "pushover_results/pushover_pos.csv"    # model pushover (v2.2/v3); falls back to Fig 13
FALLBACK_CSV = "US1_fig13_digitised.csv"
BACKBONE_BRANCH = "neg"                               # only used for the Fig 13 fallback
TABLES = {"HU12": "vel_HU.txt", "EC40": "vel_EC.txt", "FR76": "vel_FR.txt"}
PROTOCOL = [( 1,"HU12",0.50),( 2,"HU12",0.75),( 3,"EC40",0.20),( 4,"HU12",1.00),( 5,"HU12",1.25),
            ( 6,"EC40",0.30),( 7,"HU12",1.50),( 8,"EC40",0.40),( 9,"HU12",1.75),(10,"HU12",2.00),
            (11,"EC40",0.50),(12,"HU12",2.25),(13,"HU12",2.50),(14,"HU12",2.75),(15,"HU12",3.00),
            (16,"HU12",3.50),(17,"HU12",4.00),(18,"HU12",4.50),(19,"HU12",5.00),(20,"HU12",5.50),
            (21,"HU12",6.00),(22,"FR76",1.00),(23,"FR76",1.50),(24,"FR76",1.75),(25,"FR76",2.00)]
FIT_RUNS = list(range(7, 26))          # skip the tiny early runs (sub-mm, dominated by noise)
GRID = dict(m_eff=[1318, 1635, 2000, 2600, 3200],
            R=[1.0, 0.98, 0.95, 0.90],
            xi=[0.0, 0.01, 0.02])
SUBSTEPS = 6
OUT_DIR = "rocking_calibration"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================ READ MODEL PEAKS ========================
def read_hist(p):
    try:
        d = np.genfromtxt(p, skip_header=2); d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
        return d[:, 1] if d.ndim == 2 and len(d) > 5 else None
    except Exception:
        return None
def find(folder, key):
    h = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv"))); return read_hist(h[0]) if h else None
def model_peak(run):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run))) if os.path.isdir(d)]
    if not hits:
        return None
    f = hits[0]
    for key in ("rel_disp_top_exp_mm", "rel_disp_top_mm"):
        x = find(f, key)
        if x is not None:
            return float(np.max(np.abs(x - x[0])))
    c3, c4, c5 = find(f, "Channel_3_DispTopQLeft"), find(f, "Channel_4_DispTopQRight"), find(f, "Channel_5_DispTable")
    if c3 is None or c4 is None or c5 is None:
        return None
    n = min(len(c3), len(c4), len(c5)); r = 0.5 * (c3[:n] + c4[:n]) - c5[:n]
    return float(np.max(np.abs(r - r[0]))) * 1000.0

# ============================ RECORDS ================================
REC = {}
for k, fn in TABLES.items():
    if not os.path.isfile(fn):
        raise RuntimeError("missing table " + fn)
    d = np.loadtxt(fn, skiprows=2); t, v = d[:, 0], d[:, 1]; dt = float(np.median(np.diff(t)))
    REC[k] = (np.gradient(v, dt), dt)

bb_path = BACKBONE_CSV if os.path.isfile(BACKBONE_CSV) else FALLBACK_CSV
bb = rs.Backbone(bb_path, branch=BACKBONE_BRANCH)
print("backbone: {}  K1 {:.2f} kN/mm  peak {:.1f} kN at {:.1f} mm  instability {} mm".format(
    bb_path, bb.K1 / 1e6, bb.F_peak / 1e3, bb.d_peak * 1e3, "inf" if bb.d_inst == float("inf") else "{:.0f}".format(bb.d_inst * 1e3)))

targets = {}
for run, rec, sc in PROTOCOL:
    if run in FIT_RUNS:
        p = model_peak(run)
        if p is not None:
            targets[run] = p
print("model peaks available for {} of {} fit runs".format(len(targets), len(FIT_RUNS)))
if not targets:
    raise SystemExit("no control peaks found under " + CONTROL_DIR)

# ============================ GRID SEARCH ============================
rows = []; best = None
for m, R, xi in itertools.product(GRID["m_eff"], GRID["R"], GRID["xi"]):
    errs = []; pred = {}
    for run, rec, sc in PROTOCOL:
        if run not in targets:
            continue
        ag, dt = REC[rec]
        pk, _, col = rs.rocking_response(ag * sc, dt, bb, m, restitution=R, xi=xi, substeps=SUBSTEPS)
        d_s = (bb.d_inst * 1e3 if col and bb.d_inst < float("inf") else pk * 1e3) if col else pk * 1e3
        pred[run] = d_s
        errs.append(abs(np.log(max(d_s, 0.05) / max(targets[run], 0.05))))
    e = float(np.mean(errs))
    rows.append((m, R, xi, e, pred))
    print("  m {:5.0f}  R {:.2f}  xi {:.2f}  mean|log ratio| {:.3f}".format(m, R, xi, e))
    if best is None or e < best[3]:
        best = (m, R, xi, e, pred)

with open(os.path.join(OUT_DIR, "rocking_sdof_calibration.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["m_eff", "R", "xi", "mean_abs_log_ratio"] + ["run{:02d}_mm".format(r) for r in sorted(targets)])
    for m, R, xi, e, pred in rows:
        w.writerow([m, R, xi, round(e, 4)] + [round(pred[r], 2) for r in sorted(targets)])
m, R, xi, e, pred = best
json.dump({"backbone": bb_path, "branch": BACKBONE_BRANCH, "m_eff": m, "R": R, "xi": xi,
           "mean_abs_log_ratio": e, "fit_runs": sorted(targets)}, open(os.path.join(OUT_DIR, "rocking_sdof_best.json"), "w"), indent=1)
print("\nBEST: m_eff {:.0f} kg, R {:.2f}, xi {:.2f}  (mean |log ratio| {:.3f}, i.e. typical factor {:.2f})".format(m, R, xi, e, np.exp(e)))
print("-> " + os.path.join(OUT_DIR, "rocking_sdof_best.json"))

# ============================ FIGURE ================================
runs = sorted(targets)
fig, ax = plt.subplots(figsize=(9, 4.6), dpi=140)
ax.plot(runs, [targets[r] for r in runs], "s-", color="#5b7f74", label="3DEC model, record (Strategy F)")
ax.plot(runs, [pred[r] for r in runs], "o--", color="#b85042", label="rocking SDOF, best fit (m {:.0f} kg, R {:.2f}, xi {:.2f})".format(m, R, xi))
ax.set_xlabel("run"); ax.set_ylabel("peak displacement at 2.06 m (mm)"); ax.set_yscale("log")
ax.set_title("Rocking oscillator identified against the model's own record responses", loc="left", fontweight="bold")
ax.grid(color="#eeeeee", which="both"); ax.legend(fontsize=9, frameon=False)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "fig_rocking_calibration.png")); print("-> fig_rocking_calibration.png")
