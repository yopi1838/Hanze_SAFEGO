# -*- coding: ascii -*-
"""
strategy_R_3dec.py -- STRATEGY R: adaptive single-pulse sequence, period from
the state of the model, amplitude from a ROCKING oscillator match.

THE LAW (fixed once, never retuned per run)
    state   : d_hist = largest displacement the model has reached so far
              (from its own exported channels, 0.5(Ch3+Ch4)-Ch5).
    period  : T_n = secant period of the BACKBONE at d_hist,
              T = 2 pi sqrt(M_EFF / (F(d_hist)/d_hist)); T_1 = elastic period
              of the backbone (= 2 pi sqrt(M_EFF/K1)) before any displacement.
    pulse   : v(t) = Vamp sin(2 pi t / T_n), N_CYCLES cycles, tail.
    amplitude: Vamp such that the rocking oscillator (rocking_sdof.py:
              backbone, M_EFF, R, XI) reaches under the pulse the same peak
              displacement it reaches under the run's record x scale.
    Nothing is read from the experiment when BACKBONE_CSV is the model's own
    pushover; with the Fig 13 fallback the run is VALIDATION mode, and the
    log says which.

WHY A ROCKING OSCILLATOR
    Sd-matching at one period gave 85 mm (N=1) and 5 mm (N=3) for the same
    target on the cracked wall. The oscillator sees the pulse the way the
    wall does, so N sits inside the match on both sides. Whether that holds
    on the 3DEC model is what this driver tests: every run logs the
    oscillator's prediction next to the model's actual peak
    (d_rec_sdof, d_pulse_sdof, peak_edp_mm).

PARAMETERS OF THE OSCILLATOR come from rocking_calibration/rocking_sdof_best.json
(calibrate_rocking_sdof.py, fitted to the model's Strategy F responses) when
present, else from the CONFIG defaults below.

SETUP mirrors the Strategy F control (small strain, no damping, joist
contact groups, histories, free S / T_B); checkpoint/resume as route2_full25.

OUTPUT  <OUT_DIR>/RunNN_REC_sXpYY/ channel CSVs, strategyR_run_NN.sav,
        strategyR_log.csv: run, record, scale, d_hist_mm, T_pulse_s, N,
        Vamp_mps, PGA_g, stroke_mm, d_rec_sdof_mm, d_pulse_sdof_mm,
        rule_note, peak_edp_mm, ratio_model_over_sdof
"""

import itasca as it
import os, csv, math, sys, json, glob
import numpy as np

for _p in (os.getcwd(), os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)
import rocking_sdof as rs        # self-tests on import

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ CONFIG =================================
BASE_SAVE    = "Part_I_MASON_v8.sav"          # the control's undamaged start
OUT_DIR      = "strategyR_results"
BACKBONE_CSV = "pushover_results/pushover_pos.csv"   # model pushover (prediction mode)
FALLBACK_CSV = "US1_fig13_digitised.csv"             # Fig 13 (validation mode)
BACKBONE_BRANCH = "neg"
CALIB_JSON   = "rocking_calibration/rocking_sdof_best.json"
M_EFF, R_REST, XI = 1635.0, 1.0, 0.0                 # defaults if no calibration file
N_CYCLES     = 1.0
DELTA_T      = 0.005
TAIL_SEC     = 2.5
DRIVE_GROUPS = ["S", "T_B"]
T_MAX        = 1.0                                   # cap on the pulse period
PGA_WARN_G   = 0.90
VMAX         = 3.0                                   # bisection ceiling, m/s
TABLES = {"HU12": "vel_HU.txt", "EC40": "vel_EC.txt", "FR76": "vel_FR.txt"}
PROTOCOL = [( 1,"HU12",0.50),( 2,"HU12",0.75),( 3,"EC40",0.20),( 4,"HU12",1.00),( 5,"HU12",1.25),
            ( 6,"EC40",0.30),( 7,"HU12",1.50),( 8,"EC40",0.40),( 9,"HU12",1.75),(10,"HU12",2.00),
            (11,"EC40",0.50),(12,"HU12",2.25),(13,"HU12",2.50),(14,"HU12",2.75),(15,"HU12",3.00),
            (16,"HU12",3.50),(17,"HU12",4.00),(18,"HU12",4.50),(19,"HU12",5.00),(20,"HU12",5.50),
            (21,"HU12",6.00),(22,"FR76",1.00),(23,"FR76",1.50),(24,"FR76",1.75),(25,"FR76",2.00)]
os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "strategyR_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "strategyR_log.csv")

# ============================ BACKBONE + OSCILLATOR ==================
if os.path.isfile(CALIB_JSON):
    cal = json.load(open(CALIB_JSON))
    M_EFF, R_REST, XI = float(cal["m_eff"]), float(cal["R"]), float(cal["xi"])
    print("oscillator parameters from {}: m_eff {:.0f} kg, R {:.2f}, xi {:.2f}".format(CALIB_JSON, M_EFF, R_REST, XI))
else:
    print("no calibration file -- defaults m_eff {:.0f}, R {:.2f}, xi {:.2f}".format(M_EFF, R_REST, XI))
bb_path = BACKBONE_CSV if os.path.isfile(BACKBONE_CSV) else FALLBACK_CSV
MODE = "PREDICTION (model pushover)" if bb_path == BACKBONE_CSV else "VALIDATION (Fig 13 envelope)"
bb = rs.Backbone(bb_path, branch=BACKBONE_BRANCH)
T1 = 2 * math.pi * math.sqrt(M_EFF / bb.K1)
print("backbone {} [{}]: K1 {:.2f} kN/mm -> T1 {:.4f} s ; peak {:.1f} kN at {:.1f} mm ; instability {}".format(
    bb_path, MODE, bb.K1 / 1e6, T1, bb.F_peak / 1e3, bb.d_peak * 1e3,
    "none (flat)" if bb.d_inst == float("inf") else "{:.0f} mm".format(bb.d_inst * 1e3)))
if bb.d_inst == float("inf"):
    print("  ** the backbone has no descending branch: the oscillator cannot collapse and the")
    print("     large-displacement side of the rule is unconstrained. Use the v3 pushover output.")

REC = {}
for k, fn in TABLES.items():
    if not os.path.isfile(fn):
        raise RuntimeError("missing table " + fn)
    d = np.loadtxt(fn, skiprows=2); t, v = d[:, 0], d[:, 1]; dt = float(np.median(np.diff(t)))
    REC[k] = (np.gradient(v, dt), dt, float(np.max(np.abs(v))))     # accel, dt, PGV at x1.0

# ============================ HELPERS ================================
def cmd_path(p): return os.path.normpath(p).replace("\\", "/")
def run_label(n, rec, sc): return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def save_path(n): return os.path.join(OUT_DIR, "strategyR_run_{:02d}.sav".format(n))

def read_hist(p):
    try:
        d = np.genfromtxt(p, skip_header=2); d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
        return d[:, 1] if d.ndim == 2 and len(d) > 5 else None
    except Exception:
        return None
def find(folder, key):
    h = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv"))); return read_hist(h[0]) if h else None
def peak_edp_mm(folder):
    for key in ("rel_disp_top_exp_mm", "rel_disp_top_mm"):
        x = find(folder, key)
        if x is not None:
            return float(np.max(np.abs(x - x[0])))
    c3, c4, c5 = find(folder, "Channel_3_DispTopQLeft"), find(folder, "Channel_4_DispTopQRight"), find(folder, "Channel_5_DispTable")
    if c3 is None or c4 is None or c5 is None:
        return None
    n = min(len(c3), len(c4), len(c5)); r = 0.5 * (c3[:n] + c4[:n]) - c5[:n]
    return float(np.max(np.abs(r - r[0]))) * 1000.0

def setup_dynamic(save_file):
    it.command("model restore '{}'".format(cmd_path(save_file).replace(".sav", "")))
    it.command("python-reset-state false")
    it.command("model large-strain off")
    it.command("model dynamic active on")
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
    it.command("call 'instrument_history_new.dat'")
    it.command("block free velocity-z range group 'S'")
    it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")

def load_ckpt():
    if os.path.isfile(CKPT):
        s = json.load(open(CKPT)); return s["last_run"], s["summary"]
    return 0, []
def save_ckpt(last, summary):
    json.dump({"last_run": last, "summary": summary}, open(CKPT, "w"), indent=1)

# ============================ START / RESUME =========================
last_done, summary = load_ckpt()
if last_done == 0:
    if not os.path.isfile(BASE_SAVE):
        raise RuntimeError("missing " + BASE_SAVE)
    setup_dynamic(BASE_SAVE); print("fresh start from", BASE_SAVE)
else:
    setup_dynamic(save_path(last_done)); print("resuming after run", last_done)
log_new = (last_done == 0) or not os.path.isfile(LOG)
logf = open(LOG, "w" if log_new else "a", newline=""); logw = csv.writer(logf)
COLS = ["run", "record", "scale", "mode", "d_hist_mm", "T_pulse_s", "N", "Vamp_mps", "PGA_g", "stroke_mm",
        "d_rec_sdof_mm", "d_pulse_sdof_mm", "rule_note", "peak_edp_mm", "ratio_model_over_sdof"]
if log_new:
    logw.writerow(COLS)

# ============================ MAIN LOOP ==============================
for run_no, rec, sc in PROTOCOL[last_done:]:
    d_hist = max([r["peak_edp_mm"] for r in summary if r.get("peak_edp_mm") is not None] + [0.0])
    T_n = min(bb.secant_period(d_hist / 1000.0, M_EFF), T_MAX) if d_hist > 0 else T1
    ag, dt, pgv1 = REC[rec]
    rule = rs.amplitude_for_state(ag * sc, dt, T_n, N_CYCLES, bb, M_EFF, restitution=R_REST, xi=XI,
                                  tail_s=TAIL_SEC, vmax=VMAX)
    note = rule["note"]
    if not rule["converged"] or not np.isfinite(rule["Vamp"]):
        # fall back: keep the wall moving with the record's PGV x scale so the sequence continues,
        # but flag it -- this run is NOT a test of the rule
        Vamp = pgv1 * sc
        note = (note + " | " if note else "") + "FALLBACK Vamp = PGV x scale"
    else:
        Vamp = rule["Vamp"]
    Tg = rule["T"]
    spc = int(round(Tg / DELTA_T)); n = int(round(N_CYCLES * spc)); w = 2 * math.pi / Tg
    tt = np.arange(n + 1) * DELTA_T; v_p = Vamp * np.sin(w * tt)
    v_all = np.r_[v_p, np.zeros(int(round(TAIL_SEC / DELTA_T)))]; t_all = np.arange(len(v_all)) * DELTA_T
    pga = Vamp * w / 9.81; stroke = float(np.max(np.abs(np.cumsum(v_p) * DELTA_T))) * 1e3

    lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(t_all)))
        for a_, b_ in zip(t_all, v_all):
            f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

    print("\n" + "=" * 70)
    print("Run {:02d}  {} x {:.2f}   d_hist {:.2f} mm -> T {:.4f} s   [{}]".format(run_no, rec, sc, d_hist, Tg, MODE))
    print("  oscillator: d_rec {:.2f} mm  ->  pulse N {:g}: Vamp {:.4f} m/s  PGA {:.3f} g  stroke {:.1f} mm  d_pulse {:.2f} mm  {}".format(
        rule["d_rec"] * 1e3, N_CYCLES, Vamp, pga, stroke, rule["d_pulse"] * 1e3 if np.isfinite(rule["d_pulse"]) else float("nan"), note))
    if pga > PGA_WARN_G:
        print("  ** PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(pga, PGA_WARN_G))
    print("=" * 70)

    tbl = "rR_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(float(tt[-1])))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("model save '{}'".format(cmd_path(save_path(run_no))))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(folder)))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call 'instrument_history_export_v2.dat'")
    peak = peak_edp_mm(folder)
    ratio = (peak / (rule["d_rec"] * 1e3)) if (peak is not None and rule["d_rec"] > 0) else None
    print("  MODEL peak = {}   (oscillator said {:.2f} mm; model/oscillator = {})".format(
        "{:.2f} mm".format(peak) if peak is not None else "n/a", rule["d_rec"] * 1e3,
        "{:.2f}".format(ratio) if ratio else "n/a"))

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")

    row = dict(run=run_no, record=rec, scale=sc, mode=MODE.split(" ")[0], d_hist_mm=round(d_hist, 3),
               T_pulse_s=round(Tg, 5), N=N_CYCLES, Vamp_mps=round(Vamp, 5), PGA_g=round(pga, 4),
               stroke_mm=round(stroke, 2), d_rec_sdof_mm=round(rule["d_rec"] * 1e3, 3),
               d_pulse_sdof_mm=round(rule["d_pulse"] * 1e3, 3) if np.isfinite(rule["d_pulse"]) else "",
               rule_note=note.replace(",", ";"), peak_edp_mm=round(peak, 3) if peak is not None else None,
               ratio_model_over_sdof=round(ratio, 3) if ratio else "")
    summary.append(row)
    logw.writerow([row[k] if row[k] is not None else "" for k in COLS]); logf.flush()
    save_ckpt(run_no, summary)

logf.close()
print("\nStrategy R complete. Log -> {}".format(LOG))
print("Read ratio_model_over_sdof per run: near 1 = the oscillator is a faithful reduced model of the wall;")
print("then compare peak_edp_mm run by run with stratF_full_results_US1 (postprocess_routes.py, MODE='sequence').")
