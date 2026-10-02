# -*- coding: ascii -*-
"""
route2_full25_revised_3dec.py -- REVISED Route 2 over the FULL 25-run
protocol, damage state carried between runs (as Strategy F did with the
records). SMALL STRAIN, no damping, same base save / drive groups /
instrumentation as the control.

THE SIGNAL, every run (route2_revised_3dec.py, unchanged):
    one velocity cycle at T1, then one at T_eff, no gap, same phase.
    T1, T_eff  : FIXED properties of the wall from bilinear_idealisation/
                 bilinear_periods.json (block BILINEAR_KEY; "model" = the
                 v3 pushover: 0.112 / 0.241 s with m_eff 1635 kg).
    amplitudes : each cycle's 5%-damped Sd at its own period = run scale x
                 the record's Sd at that period (HU12 / EC40 / FR76 spectra).
Nothing measured on the wall enters the loop: the periods are the same on
run 1 and run 25, the amplitude follows the protocol scale. There is no
gating (the old driver held stage B back until the wall had moved 2 mm).

PARTIAL RE-RUN (added): to re-run a tail of the protocol with extra
instrumentation without touching the finished 25-run results, set
    RESUME_AFTER = 20                                   (run 21 is the first re-run)
    RESUME_SAVE  = "route2_full25_revised/route2_run_20.sav"
    OUT_DIR      = "route2_full25_revised_bs"           (a NEW folder)
The new folder gets its own checkpoint/log holding runs 21..25 only; the
original folder is never written to. Set RESUME_AFTER = 0 for the normal
full run from BASE_SAVE.

PEAK EDP per run = rel_disp_top_exp_mm channel (0.5(Ch3+Ch4) - Ch5), else
rebuilt from Channel_3/4/5 -- NOT Record_Disp.

RESUME: OUT_DIR/route2_checkpoint.json holds the last completed run; a
restart restores that run's save and continues. Delete it to start over.

OUTPUT (same layout as route2_full25_3dec.py):
    OUT_DIR/RunNN_REC_sXpYY/        channel CSVs
    OUT_DIR/route2_run_NN.sav
    OUT_DIR/RunNN_REC_sXpYY_vel.txt the applied table
    OUT_DIR/route2_full_log.csv
REQUIRES next to this script: bilinear_idealisation/bilinear_periods.json,
spectrum_HU12/EC40/FR76.csv, Part_I_MASON_v8_SmallStrain.sav,
instrument_history_new.dat, instrument_baseshear_exp.dat,
instrument_history_export_v2.dat (with the vz_acc15/16/17 export lines added).
"""

import itasca as it
import os, csv, math, sys, json, glob
import numpy as np

# ============================ CONFIG =================================
OUT_DIR       = "route2_full25_revised_bs"       # NEW folder for the re-run (originals untouched)
BASE_SAVE     = "Part_I_MASON_v8_SmallStrain.sav"
RESUME_AFTER  = 20                               # 0 = full run from BASE_SAVE; N = start at run N+1 from RESUME_SAVE
RESUME_SAVE   = "route2_full25_revised/route2_run_20.sav"
BILINEAR_JSON = "bilinear_idealisation/bilinear_periods.json"
BILINEAR_KEY  = "model"        # "model" (v3 pushover, prediction mode) or "fig13" (envelope)
STAGE_B_PHASE = "same"         # "same" or "inverted"
XI            = 0.05
DELTA_T       = 0.005
TAIL_SEC      = 2.5            # ring-down after the table stops
DRIVE_GROUPS  = ["S", "T_B"]
PGA_WARN_G    = 0.90           # the table reached 0.78 g; above this is extrapolation
INSTRUMENT_FILES = ["instrument_history_new.dat", "instrument_baseshear_exp.dat"]   # called together, every time
SPECTRA = {"HU12": "spectrum_HU12.csv", "EC40": "spectrum_EC40.csv", "FR76": "spectrum_FR76.csv"}
PROTOCOL = [
    ( 1, "HU12", 0.50), ( 2, "HU12", 0.75), ( 3, "EC40", 0.20),
    ( 4, "HU12", 1.00), ( 5, "HU12", 1.25), ( 6, "EC40", 0.30),
    ( 7, "HU12", 1.50), ( 8, "EC40", 0.40), ( 9, "HU12", 1.75),
    (10, "HU12", 2.00), (11, "EC40", 0.50), (12, "HU12", 2.25),
    (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00),
    (16, "HU12", 3.50), (17, "HU12", 4.00), (18, "HU12", 4.50),
    (19, "HU12", 5.00), (20, "HU12", 5.50), (21, "HU12", 6.00),
    (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75),
    (25, "FR76", 2.00),
]
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
if STAGE_B_PHASE not in ("same", "inverted"):
    raise RuntimeError("STAGE_B_PHASE must be 'same' or 'inverted'")
PHASE_SIGN = 1.0 if STAGE_B_PHASE == "same" else -1.0
os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "route2_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "route2_full_log.csv")
LOG_COLS = ["run", "record", "scale", "d_hist_mm", "two_stage", "T_A_s", "SdA_mm", "VampA_mps", "PGA_A_g",
            "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g", "gap_s", "cap_note", "peak_edp_mm",
            "join_step_g", "SdA_total_mm", "SdB_total_mm", "period_source", "edp_source"]

# ============================ Sd INTEGRATOR (inlined, self-tested) ====
def sd_spectrum(ag, dt, Ts, xi=0.05):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2.0 * np.pi / Ts; c = 2.0 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0] - c * v - w2 * u; umax = np.abs(u)
    for i in range(1, len(ag)):
        u_new = u + dt * v + 0.25 * dt * dt * a
        v_new = v + 0.5 * dt * a
        a_new = (-ag[i] - c * v_new - w2 * u_new) / (1.0 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = u_new + 0.25 * dt * dt * a_new; v = v_new + 0.5 * dt * a_new; a = a_new
        umax = np.maximum(umax, np.abs(u))
    return umax
def _sd_selftest():
    T, xi, A, dt = 0.4, 0.05, 1.0, 0.001
    t = np.arange(int(60 * T / dt)) * dt
    got = sd_spectrum(A * np.sin(2 * np.pi / T * t), dt, [T], xi)[0]
    want = A / ((2 * np.pi / T) ** 2 * 2 * xi)
    if abs(got - want) / want > 0.01:
        raise RuntimeError("Sd integrator self-test failed: {} vs {}".format(got, want))
_sd_selftest()

# ============================ PERIODS (fixed) + SPECTRA ==============
if not os.path.isfile(BILINEAR_JSON):
    raise RuntimeError("run bilinear_idealisation.py first: {} not found".format(BILINEAR_JSON))
_b = json.load(open(BILINEAR_JSON))
if BILINEAR_KEY not in _b:
    raise RuntimeError("{} has no '{}' block (present: {})".format(BILINEAR_JSON, BILINEAR_KEY, ", ".join(k for k in _b if k != "rule")))
_b = _b[BILINEAR_KEY]
T1, T_EFF = float(_b["T1_s"]), float(_b["Teff_s"])
PERIOD_SOURCE = "bilinear/{} ({}; K1 {:.2f} kN/mm, mu {:.2f}, m_eff {:.0f})".format(
    BILINEAR_KEY, "+".join(_b["branches_used"]), _b["K1_kNmm"], _b["mu"], json.load(open(BILINEAR_JSON))["rule"]["m_eff"])
SPEC = {}
for k, f in SPECTRA.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing spectrum file: {}".format(f))
    d = np.genfromtxt(f, delimiter=",", skip_header=1); SPEC[k] = (d[:, 0], d[:, 1])
def sd_record(rec, T):
    Ts, Sd = SPEC[rec]; return float(np.interp(T, Ts, Sd))
print("REVISED ROUTE 2, full protocol: T1 {:.4f} s, T_eff {:.4f} s on every run  [{}]; stage B phase {}".format(
    T1, T_EFF, PERIOD_SOURCE, STAGE_B_PHASE))

# ============================ SIGNAL BUILDERS =======================
def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")
def grid_period(T):
    spc = max(1, int(round(T / DELTA_T))); return spc * DELTA_T, spc
def one_cycle(Vamp, T):
    Tg, spc = grid_period(T); tt = np.arange(spc + 1) * DELTA_T
    return tt, Vamp * np.sin(2.0 * math.pi * tt / Tg), Tg
def calibrate_vel(T, Sd_target):
    _, v, Tg = one_cycle(1.0, T)
    sd0 = float(sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [Tg], XI)[0])
    if sd0 <= 0:
        raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg
def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81)
def build_signal(rec, scale):
    SdA, SdB = scale * sd_record(rec, T1), scale * sd_record(rec, T_EFF)
    VA, TgA = calibrate_vel(T1, SdA); VB, TgB = calibrate_vel(T_EFF, SdB)
    _, vA, _ = one_cycle(VA, TgA); _, vB, _ = one_cycle(PHASE_SIGN * VB, TgB)
    v = np.r_[vA, vB[1:]]; t = np.arange(len(v)) * DELTA_T
    step = abs((vB[1] - vB[0]) - (vA[-1] - vA[-2])) / DELTA_T / 9.81
    tot = sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [TgA, TgB], XI)
    return dict(t=t, v=v, TgA=TgA, TgB=TgB, SdA=SdA, SdB=SdB, VA=VA, VB=VB,
                pgaA=pga_g(vA), pgaB=pga_g(vB), step=step, SdA_tot=float(tot[0]), SdB_tot=float(tot[1]))

# ============================ MODEL SETUP ============================
def call_instruments():
    """every instrumentation file, in order -- used at setup AND after each run's `history delete`."""
    for f in INSTRUMENT_FILES:
        if not os.path.isfile(f):
            raise RuntimeError("missing instrumentation file: " + f)
        it.command("call '{}'".format(f))

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
    call_instruments()
    it.command("block free velocity-z range group 'S'")
    it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")

def run_label(n, rec, sc):
    return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def save_path(n):
    return os.path.join(OUT_DIR, "route2_run_{:02d}.sav".format(n))

def read_hist(path):
    try:
        d = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None
def find_ch(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else None
def peak_edp_mm(folder):
    x = find_ch(folder, KEY_DISP_EXP)
    if x is not None:
        return float(np.max(np.abs(x - x[0]))), KEY_DISP_EXP
    c3, c4, c5 = find_ch(folder, CH3), find_ch(folder, CH4), find_ch(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, "n/a"
    n = min(len(c3), len(c4), len(c5)); rel = 0.5 * (c3[:n] + c4[:n]) - c5[:n]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0, "0.5(Ch3+Ch4)-Ch5"

def load_ckpt():
    if os.path.isfile(CKPT):
        with open(CKPT) as f:
            s = json.load(f)
        return s["last_run"], s["summary"]
    return 0, []
def save_ckpt(last, summary):
    with open(CKPT, "w") as f:
        json.dump({"last_run": last, "summary": summary, "resume_after": RESUME_AFTER, "resume_save": RESUME_SAVE}, f, indent=1)

# ============================ CHECKPOINT =============================
it.command("python-reset-state false")
it.command("program automatic-model-save active off")
last_done, summary = load_ckpt()
if last_done == 0:
    if RESUME_AFTER > 0:
        # partial re-run: state after run RESUME_AFTER, taken from another folder's save
        if not os.path.isfile(RESUME_SAVE):
            raise RuntimeError("RESUME_SAVE not found: " + RESUME_SAVE)
        if os.path.abspath(os.path.dirname(RESUME_SAVE)) == os.path.abspath(OUT_DIR):
            raise RuntimeError("OUT_DIR must be a NEW folder, not the one holding RESUME_SAVE")
        setup_dynamic(RESUME_SAVE); last_done = RESUME_AFTER
        print("partial re-run: state after run {} from {}; runs {}..{} go to {}".format(
            RESUME_AFTER, RESUME_SAVE, RESUME_AFTER + 1, PROTOCOL[-1][0], OUT_DIR))
    else:
        if not os.path.isfile(BASE_SAVE):
            raise RuntimeError("missing " + BASE_SAVE)
        setup_dynamic(BASE_SAVE); print("fresh start from", BASE_SAVE)
else:
    if last_done < RESUME_AFTER:
        raise RuntimeError("checkpoint says last_run {} but RESUME_AFTER is {}: delete {} or fix the config".format(last_done, RESUME_AFTER, CKPT))
    setup_dynamic(save_path(last_done)); print("resuming after run", last_done)

log_new = (len(summary) == 0) or not os.path.isfile(LOG)
logf = open(LOG, "w" if log_new else "a", newline="")
logw = csv.writer(logf)
if log_new:
    logw.writerow(LOG_COLS)

# ============================ MAIN LOOP ==============================
for run_no, rec, sc in PROTOCOL[last_done:]:
    d_hist = max([r["peak_edp_mm"] for r in summary if r.get("peak_edp_mm") is not None] + [0.0])   # logged only, not used
    S = build_signal(rec, sc)
    lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(S["t"])))
        for a_, b_ in zip(S["t"], S["v"]):
            f.write("{:.6f}\t{:.9e}\n".format(a_, b_))
    pulse_dur = float(S["t"][-1])

    print("\n" + "=" * 66)
    print("Run {:02d}  {} x {:.2f}   (largest peak so far {:.2f} mm, informative only)".format(run_no, rec, sc, d_hist))
    print("  stage A: T {:.3f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.2f}g".format(S["TgA"], S["SdA"] * 1e3, S["VA"], S["pgaA"]))
    print("  stage B: T {:.3f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.2f}g   join step {:.2f} g".format(S["TgB"], S["SdB"] * 1e3, S["VB"], S["pgaB"], S["step"]))
    print("  whole signal: Sd(T1) {:.2f} mm, Sd(T_eff) {:.2f} mm ; table {:.3f} s + tail {:.1f} s".format(S["SdA_tot"] * 1e3, S["SdB_tot"] * 1e3, pulse_dur, TAIL_SEC))
    for tag, p in (("A", S["pgaA"]), ("B", S["pgaB"])):
        if p > PGA_WARN_G:
            print("  ** stage {} PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(tag, p, PGA_WARN_G))
    print("=" * 66)

    tbl = "r2_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(pulse_dur))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("model save '{}'".format(cmd_path(save_path(run_no))))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(folder)))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call 'instrument_history_export_v2.dat'")
    peak, src = peak_edp_mm(folder)
    print("  PEAK EDP ({}) = {}".format(src, "{:.2f} mm".format(peak) if peak is not None else "n/a"))
    if find_ch(folder, "vz_acc15") is None:
        print("  ! vz_acc15 not exported -- add the export lines from instrument_baseshear_exp.dat to instrument_history_export_v2.dat")

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    call_instruments()

    row = dict(run=run_no, record=rec, scale=sc, d_hist_mm=round(d_hist, 3), two_stage=1,
               T_A_s=round(S["TgA"], 5), SdA_mm=round(S["SdA"] * 1e3, 3), VampA_mps=round(S["VA"], 5), PGA_A_g=round(S["pgaA"], 4),
               T_B_s=round(S["TgB"], 5), SdB_mm=round(S["SdB"] * 1e3, 3), VampB_mps=round(S["VB"], 5), PGA_B_g=round(S["pgaB"], 4),
               gap_s=0.0, cap_note="fixed_periods", peak_edp_mm=round(peak, 3) if peak is not None else None,
               join_step_g=round(S["step"], 3), SdA_total_mm=round(S["SdA_tot"] * 1e3, 3), SdB_total_mm=round(S["SdB_tot"] * 1e3, 3),
               period_source=PERIOD_SOURCE.replace(",", ";"), edp_source=src)
    summary.append(row)
    logw.writerow([row[k] if row[k] is not None else "" for k in LOG_COLS])
    logf.flush()
    save_ckpt(run_no, summary)

logf.close()
print("\ndone: {} runs in {}. Loops: python baseshear_exp_from_model.py {} 21 22 23 24 25".format(len(summary), OUT_DIR, OUT_DIR))
