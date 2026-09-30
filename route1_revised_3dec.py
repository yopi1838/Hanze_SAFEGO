# -*- coding: ascii -*-
"""
route1_extract_3dec.py -- ROUTE 1 BY EXTRACTION, runs 23 and 24.

THE SIGNAL: the FR76 record itself over its significant velocity window,
zero-padded after -- NOT a sinusoid fitted to it. No period definition, no
amplitude convention, no symmetry assumption.

    window  12.461 -> 14.106 s of vel_FR.txt  (1.640 s)
            = the two big lobes (the boxed window) + the negative lobe after
    delivers 101% of the full record's Sd at T1 and 98% at T_eff, from 1.64 s
    of table instead of 17.35 s -- a 10x saving.
    For comparison, the best SINUSOID over the same window gave 31% and 20%.

SOURCE FILE: vel_FR.txt, the MEASURED table velocity -- the same signal
Strategy F drives the model with and the one spectrum_FR76.csv was built from.
route1_velocity_pulse_3dec.py read Friuli_x1.0.csv (the COMMAND record, PGV
0.2202 vs 0.2524 m/s), so every earlier Route 1 pulse was 13% under-scaled.
That is fixed here by construction.

BASELINE CORRECTION (BASELINE_ZERO, default ON). A truncated window does not
integrate to zero: raw, it leaves the base 27 mm displaced at x1.00 and 47 mm
at run 24's x1.75, which would contaminate every table-referenced EDP. The
correction subtracts k*sin(pi*(t-t0)/Tw), which is zero at both ends so it
adds no velocity step, with k set to zero the net displacement. It costs
almost nothing spectrally: Sd at T_eff 98% -> 97%, at T1 101% -> 102%, and
peak velocity 0.2524 -> 0.2276 m/s.
    NB the full vel_FR.txt ALSO fails to return to zero (+18 mm), so the
    Strategy F control carries an offset of the same order. Worth checking
    against the experiment's own table displacement channel before trusting
    any residual/permanent-displacement result from either.

SCALING: the table is written at 100% and each run applies its Table-2 factor
through the velocity-z multiplier, exactly as Strategy F does.
    run 23  FR76 x 1.50        run 24  FR76 x 1.75
Damage is carried from 23 into 24.

START STATE: runs 23 and 24 are late in the protocol, so starting from the
undamaged save is NOT the same experiment. START_SAVE should be the state
after run 22 from the matching Strategy F case. If it is not found the script
falls back to BASE_SAVE and says so loudly -- results from that fallback are
not comparable to the control.

Usage (3DEC):  model new ; python-reset-state false ; call 'route1_extract_3dec.py'
Output: OUT_DIR/RunNN_FR76_sXpYY/  channel CSVs
        OUT_DIR/route1_run_NN.sav
        OUT_DIR/route1_extract_vel.txt   the applied table (100%)
        OUT_DIR/route1_extract_log.csv
        OUT_DIR/case.txt
"""

import itasca as it
import os, csv, math, glob, json
import numpy as np

# ============================ THE CASE ===============================
# Must match whichever control this is being compared against.
CASE = "SS_NODAMP"      # SS_NODAMP | SS_MAXWELL_1p5 | LS_NODAMP | LS_MAXWELL_1p5
CASE = os.environ.get("ROUTE1_CASE", CASE).strip()
MAXWELL_CMD = "block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 0.0104 40.0000"
SAVE_SS = "Part_I_MASON_v8_SmallStrain.sav"
SAVE_LS = "part_I_mason_LS_strong.sav"
CASES = {
    "SS_NODAMP":      dict(ls=False, damp=False, save=SAVE_SS, ctrl="stratF_full_SS_NODAMP"),
    "SS_MAXWELL_1p5": dict(ls=False, damp=True,  save=SAVE_SS, ctrl="stratF_full_SS_MAXWELL_1p5"),
    "LS_NODAMP":      dict(ls=True,  damp=False, save=SAVE_LS, ctrl="stratF_full_LS_NODAMP"),
    "LS_MAXWELL_1p5": dict(ls=True,  damp=True,  save=SAVE_LS, ctrl="stratF_full_LS_MAXWELL_1p5"),
}
if CASE not in CASES:
    raise RuntimeError("CASE must be one of {} (got '{}')".format(", ".join(sorted(CASES)), CASE))
_C = CASES[CASE]
LS_ON, DAMP_ON, BASE_SAVE, CONTROL_DIR = _C["ls"], _C["damp"], _C["save"], _C["ctrl"]
DAMP_CMD = MAXWELL_CMD if DAMP_ON else ""

# ============================ CONFIG =================================
SRC_TABLE     = "vel_FR.txt"        # measured table velocity, 100% FR76
WIN_A, WIN_B  = 12.461, 14.106      # the extraction window, s
BASELINE_ZERO = True                # zero the net base displacement (see docstring)
RUNS          = [(23, 1.50), (24, 1.75)]
PRIOR_RUN     = 22                  # the save this sequence must start from
DELTA_T       = 0.005
TAIL_SEC      = 2.5
DRIVE_GROUPS  = ["S", "T_B"]
OUT_DIR       = "route1_extract_" + CASE
START_SAVE    = None                # None = auto-discover run-22 save, else give a path
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
os.makedirs(OUT_DIR, exist_ok=True)
LOG = os.path.join(OUT_DIR, "route1_extract_log.csv")
LOG_COLS = ["run", "record", "scale", "window_a_s", "window_b_s", "dur_s", "baseline_zero",
            "vpeak_applied_mps", "net_base_disp_mm", "peak_edp_mm", "resid_edp_mm",
            "edp_source", "case", "large_strain", "damping_cmd", "start_save"]

def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")

# ============================ BUILD THE TABLE ========================
if not os.path.isfile(SRC_TABLE):
    raise RuntimeError("missing source table: " + SRC_TABLE)
_d = np.genfromtxt(SRC_TABLE, skip_header=2)
_t, _v = _d[:, 0], _d[:, 1]
_m = (_t >= WIN_A) & (_t <= WIN_B)
if _m.sum() < 10:
    raise RuntimeError("window {}-{} s selects {} samples of {}".format(WIN_A, WIN_B, _m.sum(), SRC_TABLE))
tw, vw = _t[_m] - _t[_m][0], _v[_m]
TW = float(tw[-1])
trapz = getattr(np, "trapezoid", None) or np.trapz      # numpy 2.x removed np.trapz
NET_RAW = float(trapz(vw, tw))
if BASELINE_ZERO:
    # zero at both ends -> no velocity step is introduced at either edge
    k = NET_RAW / (2.0 * TW / math.pi)
    vw = vw - k * np.sin(math.pi * tw / TW)
else:
    k = 0.0
# force exact zero endpoints so the table starts and ends at rest
vw[0] = 0.0; vw[-1] = 0.0
NET_FIN = float(trapz(vw, tw))
VPEAK = float(np.max(np.abs(vw)))

VEL_PATH = os.path.join(OUT_DIR, "route1_extract_vel.txt")
with open(VEL_PATH, "w", newline="\n") as f:
    f.write("route1_extract\n{}\t0\n".format(len(tw)))
    for a_, b_ in zip(tw, vw):
        f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

# ============================ START STATE ============================
def find_start():
    if START_SAVE:
        return START_SAVE, "explicit"
    pats = [os.path.join(CONTROL_DIR, "stratC_run_{:02d}.sav".format(PRIOR_RUN)),
            os.path.join(CONTROL_DIR, "*run_{:02d}.sav".format(PRIOR_RUN)),
            os.path.join("stratF_full_*", "*run_{:02d}.sav".format(PRIOR_RUN))]
    for p in pats:
        hits = sorted(glob.glob(p))
        if hits:
            return hits[0], "auto"
    return BASE_SAVE, "FALLBACK"
START, START_HOW = find_start()

print("=" * 72)
print("ROUTE 1 BY EXTRACTION -- runs {} and {}   [case {}]".format(RUNS[0][0], RUNS[1][0], CASE))
print("  source      : {}  (measured table velocity, 100% FR76)".format(SRC_TABLE))
print("  window      : {:.3f} -> {:.3f} s  ({:.3f} s, {} samples)".format(WIN_A, WIN_B, TW, len(tw)))
print("  baseline    : {}  net base displacement {:+.1f} mm -> {:+.3f} mm".format(
    "ZEROED" if BASELINE_ZERO else "raw (NOT zeroed)", NET_RAW * 1000, NET_FIN * 1000))
if not BASELINE_ZERO:
    print("  ** the base will end {:+.0f} mm off at run 24's x1.75 -- every table-referenced".format(
        NET_RAW * 1000 * RUNS[-1][1]))
    print("     EDP will carry that offset. Set BASELINE_ZERO = True unless you want this.")
print("  peak vel    : {:.4f} m/s at 100%  ->  {:.4f} m/s at x{:.2f}".format(
    VPEAK, VPEAK * RUNS[-1][1], RUNS[-1][1]))
print("  geometry    : large-strain {}".format("on" if LS_ON else "off"))
print("  damping     : {}".format(DAMP_CMD if DAMP_CMD else "none (contact dissipation only)"))
print("  start save  : {}   [{}]".format(START, START_HOW))
if START_HOW == "FALLBACK":
    print("  " + "!" * 68)
    print("  ! No run-{:02d} save found under {}/.".format(PRIOR_RUN, CONTROL_DIR))
    print("  ! Starting from the UNDAMAGED base save instead. Runs 23-24 applied to a")
    print("  ! virgin wall are NOT the protocol's runs 23-24 and are NOT comparable to")
    print("  ! the control. Set START_SAVE explicitly, or run the control first.")
    print("  " + "!" * 68)
print("  output      : {}".format(OUT_DIR))
print("  control     : {} (same window, same case -> loading is the only difference)".format(CONTROL_DIR))
print("=" * 72)
with open(os.path.join(OUT_DIR, "case.txt"), "w") as f:
    json.dump(dict(case=CASE, large_strain="on" if LS_ON else "off", damping=DAMP_CMD,
                   base_save=BASE_SAVE, start_save=START, start_how=START_HOW,
                   source=SRC_TABLE, window=[WIN_A, WIN_B], dur_s=TW,
                   baseline_zero=BASELINE_ZERO, net_raw_mm=NET_RAW * 1000,
                   net_final_mm=NET_FIN * 1000, vpeak_mps=VPEAK, control=CONTROL_DIR), f, indent=1)

# ============================ MODEL SETUP ============================
def setup_dynamic(save_file):
    it.command("model restore '{}'".format(cmd_path(save_file).replace(".sav", "")))
    it.command("python-reset-state false")
    it.command("model large-strain {}".format("on" if LS_ON else "off"))
    it.command("model dynamic active on")
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    if DAMP_CMD:
        it.command(DAMP_CMD)
        print("  damping applied: {}".format(DAMP_CMD))
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

def read_hist(p):
    try:
        d = np.genfromtxt(p, skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None
def find_ch(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else None
def edp(folder):
    x = find_ch(folder, KEY_DISP_EXP)
    src = KEY_DISP_EXP
    if x is None:
        c3, c4, c5 = find_ch(folder, CH3), find_ch(folder, CH4), find_ch(folder, CH5)
        if c3 is None or c4 is None or c5 is None:
            return None, None, "n/a"
        n = min(len(c3), len(c4), len(c5))
        x = (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0
        src = "0.5(Ch3+Ch4)-Ch5"
    n = max(1, int(0.05 * len(x)))
    return float(np.max(np.abs(x - x[0]))), float(np.mean(x[-n:]) - x[0]), src

# ============================ RUN ====================================
it.command("python-reset-state false")
it.command("program automatic-model-save active off")
if not os.path.isfile(START):
    raise RuntimeError("start save not found: " + START)
setup_dynamic(START)
it.command("table 'r1_ext' import '{}'".format(cmd_path(VEL_PATH)))

logf = open(LOG, "w", newline="")
logw = csv.writer(logf); logw.writerow(LOG_COLS)
for run_no, sc in RUNS:
    lbl = "Run{:02d}_FR76_s{}".format(run_no, "{:.2f}".format(sc).replace(".", "p"))
    print("\n" + "=" * 66)
    print("Run {:02d}  FR76 x {:.2f}   extracted window {:.3f} s + tail {:.1f} s".format(
        run_no, sc, TW, TAIL_SEC))
    print("  applied peak velocity {:.4f} m/s".format(VPEAK * sc))
    print("=" * 66)
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z {:g} table 'r1_ext' range group '{}'".format(sc, g))
    it.command("model solve dynamic time {:.6f}".format(TW))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    sav = os.path.join(OUT_DIR, "route1_run_{:02d}.sav".format(run_no))
    it.command("model save '{}'".format(cmd_path(sav)))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(folder)))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call 'instrument_history_export_v2.dat'")
    pk, rs, src = edp(folder)
    print("  PEAK EDP ({}) = {}{}".format(
        src, "{:.2f} mm".format(pk) if pk is not None else "n/a",
        "   settled {:+.2f} mm".format(rs) if rs is not None else ""))

    logw.writerow([run_no, "FR76", sc, WIN_A, WIN_B, round(TW, 4), int(BASELINE_ZERO),
                   round(VPEAK * sc, 5), round(NET_FIN * 1000 * sc, 3),
                   round(pk, 3) if pk is not None else "", round(rs, 3) if rs is not None else "",
                   src, CASE, "on" if LS_ON else "off", DAMP_CMD.replace(",", ";"), START])
    logf.flush()
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")
logf.close()

print("\ndone. {} and the run folders are in {}.".format(os.path.basename(LOG), OUT_DIR))
print("US-1 measured: run 23 -> 17.0 mm (approx), run 24 -> 29.5 mm at Z = 2.06 m.")
print("Compare against {} runs 23-24 for the loading-only difference.".format(CONTROL_DIR))