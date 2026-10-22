# -*- coding: ascii -*-
"""
Strategy F (FULL) -- recorded table motions, runs 1-25: the CONTROL on the corrected wall.

The fidelity counterpart to the pulse strategies, run as the FULL cumulative
sequence from run 1. There are THREE base velocity tables -- the HU12 / EC40 /
FR76 records at 100% intensity (vel_HU.txt / vel_EC.txt / vel_FR.txt, written by
integrate_table_accel.py) -- and each run applies its record's table SCALED by
the Table-2 factor through the 3DEC velocity-z multiplier. Nothing is
synthesised, no spectrum is read, no intensity is prescribed; the only per-run
quantity is the linear scale factor from Table 2.

SCALING (Table 2)
    Run 1  -> velocity-z 0.50 table 'vel_HU'      (HU12 x 0.50)
    Run 3  -> velocity-z 0.20 table 'vel_EC'      (EC40 x 0.20)
    Run 22 -> velocity-z 1.00 table 'vel_FR'      (FR76 x 1.00)
    ... 25 runs, HU12 0.50-6.00, EC40 0.20-0.50, FR76 1.00-2.00. See
    PROTOCOL_TABLE2.

TRUNCATION (done upstream in integrate_table_accel.py)
    Each base table is the channel-12 acceleration integrated to velocity, kept
    from t = 0 and cut at the velocity zero-crossing nearest its 95% Arias time,
    so the applied motion starts and ends at rest with no velocity step.

REPORTED
    T_end from the ring-down after each record, as a damage measure. It drives
    nothing. Run folders are Run{NN}_{record}_s{scale}, so postprocess_stratC.py
    parses them unchanged and the comparison figures work as-is.

GENERATED FROM strategy_F_3dec.py: only the config block differs. The
pulse-synthesis helpers remain in the file but are unused here.
run_strategy_C() is kept as the entry-point name only, for diff-ability.
"""

# =====================================================================
# STANDALONE CASE: F_LS_MAXWELL15_v12  -- the CONTROL on the corrected wall
# =====================================================================
# This file is a COMPLETE copy of the driver with the case constants baked in.
# It reads no overrides file. Run it directly:
#
#     python-reset-state false
#     call '<this file>'
#
#   build      : 10.35 kN spring load as an APPLIED FORCE on T_B (no equivalent
#                density, no follower traction), f_t ~ 0 / c = 0.2 MPa joints,
#                settle Part_I_MASON_v8_LargeStrain_springforce_ft0
#   geometry   : large-strain on
#   damping    : Maxwell, three units at 0.0085 (1.5 / 6 / 25 Hz) = 1.5 %
#                MEASURED on this wall by ring-down (ringdown_check v2).
#                NOT the 1-40 Hz triplet of strategy_C_3dec_maxwell.py, which
#                was never verified on this build.
#   period ID  : experiment   (T1_init = 0.103 s; undamped ring-down 9.7 Hz)
#   output     : stratF_full_LS_MAXWELL15_v12
#
# NOT comparable with stratF_full_SS_NODAMP_LATEST (small strain, undamped,
# 2x springs as density, f_t = 0.344). This is a NEW BASELINE: every surrogate
# (Method 2, Route 2, ...) is scored against this sequence.
#
# A fresh start aborts if the save's joint tension is above TENSION_MAX_PA, so
# pointing this at the bonded build fails before any solving.
# =====================================================================

import itasca as it
import numpy as np
import os, csv, math, json, sys

# Experiment-matched period identification (Moshfeghi et al. 2024 sec. 3.2).
try:
    import period_id_exp
    _HAVE_PERIOD_ID_EXP = True
except Exception as _e:                                  # pragma: no cover
    period_id_exp = None
    _HAVE_PERIOD_ID_EXP = False
    print("WARNING: period_id_exp not importable ({}). Falling back to the "
          "legacy single-channel displacement estimator.".format(_e))

PERIOD_ID_METHOD = "experiment"
_T_TRACK = None

# >>> RATCHETING : make ratcheting_pulse.py importable and import it
_here = os.getcwd()
for _p in (_here, os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else _here):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)
try:
    from ratcheting_pulse import (RECORD_ASYM, beta_eff,
                                  calibrate_amplitude_asym, build_velocity_asym)
except Exception as _e:
    raise RuntimeError("Cannot import ratcheting_pulse.py (put it next to this "
                       "script). Original error: {}".format(_e))
# <<< RATCHETING

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# =====================================================================
# 1.  PARAMETERS
# =====================================================================
T1_init     = 0.103      # s. Corrected wall rings at 9.7 Hz undamped (was 0.0948).
                         # Feeds the period-ID windows and the T_end/T_init
                         # column only; never the loading.
xi          = 0.05
delta_t     = 0.005
n_cycles    = 1.5
tail_sec    = 2.5
phase_accel = math.pi / 2
inter_run_gap = 0.5
bLength     = 1.292

FFT_F_MIN  = 2.0
FFT_F_MAX  = 50.0
T_MIN_PHYS = 0.02
T_MAX_PHYS = 1.0

OUT_DIR = None   # set just below, once RECORD_TEST exists
STATE_FILE_NAME = "stratC_checkpoint.json"

# >>> RATCHETING : controls
USE_RATCHETING = False   # OFF: symmetric pulse helpers, unused in Strategy F
ASYM_K         = 1.0
# Rayleigh damping OFF. Damping is Maxwell (below).
DAMP_RATIO = 0.0
DAMP_TYPE  = ""

# Maxwell (frequency-band) damping, issued verbatim after local/global damping
# are zeroed. Empty string = no Maxwell damping.
#
# 1.5 % over 1.5-25 Hz on THIS build, verified by ring-down (ringdown_check v2):
# three units at 0.0085 each superpose to xi = 0.015 +- 0.001 at ~10 Hz
# (0.015 each gave 0.0266). Adds ~1.6 % to the elastic frequency, < 1 % in the
# rocking range; static stiffness unaffected.
#
# CAUTION (from the JModelMason::run() note): getMaxShearStiffness() returns a
# constant ks_ regardless of contact state, so if 3DEC applies the Maxwell shear
# dashpot to open subcontacts a separated joint transmits viscous shear across
# the gap. A full record keeps cracked joints open far longer than a pulse.
# Check: after run 19 exists, repeat run 19 from stratC_run_18.sav with
# MAXWELL_CMD = '' and compare peaks; more than ~5 % difference means the
# dashpot is acting across open joints and the damping model must change.
MAXWELL_CMD = 'block mech damp maxwell 0.0112 1.0005 0.0081 7.4315 0.0096 25.0004'   # 1.5 % flat, 1-25 Hz (MATLAB-optimised)

# Geometry mode. Applied after the restore, so the static build is untouched.
LARGE_STRAIN = "on"

# >>> STRATEGY F (FULL) : recorded table motions, runs 1-25.
RECORD_TEST   = 9          # only used for the output-folder suffix
TABLE_DIR     = "."        # folder holding vel_HU.txt / vel_EC.txt / vel_FR.txt
BASE_TABLES = {
    "HU12": ("vel_HU", "vel_HU.txt"),
    "EC40": ("vel_EC", "vel_EC.txt"),
    "FR76": ("vel_FR", "vel_FR.txt"),
}
RUN_FROM      = 1          # first run to apply (inclusive)
RUN_TO        = 25         # last run to apply (inclusive)
START_SAVE    = "Part_I_MASON_v8_LargeStrain_springforce.sav"   # f_t ~ 0 / c = 0.2 build
TAIL_SEC_RECORD = 2.5      # ring-down appended after each record (period ID)
TENSION_MAX_PA  = 0.05e6   # fresh start aborts if the save's joint tension exceeds this
# <<< STRATEGY F (FULL)

OUT_DIR = 'stratF_full_LS_MAXWELL15_v12'
CASE_ID = 'F_LS_MAXWELL15_v12'

# Overrides DISABLED: standalone case file.
_OVERRIDABLE_F = ("OUT_DIR", "START_SAVE", "RUN_FROM", "RUN_TO", "CASE_ID",
                  "DAMP_RATIO", "DAMP_TYPE", "MAXWELL_CMD", "LARGE_STRAIN",
                  "PERIOD_ID_METHOD", "T1_init", "TAIL_SEC_RECORD",
                  "RECORD_TEST")
_OVR_FILE_F = "stratF_overrides.json"
print("\n(standalone case {}: overrides file ignored)".format(CASE_ID))

# >>> RETIRED : the Strategy C feedback controls. None is read here.
USE_SECANT_PERIOD = False
M_EFF_T           = 1.47
K_INIT_KNPM       = 6875.0
T_SEC_MAX         = 0.60
SECANT_MONOTONIC  = True
SD_GROWTH_CAP     = 3.0
SCALE_MULT        = 1.0
# <<< SECANT PERIOD

# >>> GI VARIANT : joint property overrides applied after every model restore.
G_I_NEW      = None      # leave the base save's properties alone
G_II_RATIO   = 10.0
COH_RESIDUAL = 0.0
# <<< GI VARIANT

# =====================================================================
# 2.  PROTOCOL (Table 2)
# =====================================================================
PROTOCOL_TABLE2 = [
    ( 1, "HU12", 0.50), ( 2, "HU12", 0.75), ( 3, "EC40", 0.20), ( 4, "HU12", 1.00),
    ( 5, "HU12", 1.25), ( 6, "EC40", 0.30), ( 7, "HU12", 1.50), ( 8, "EC40", 0.40),
    ( 9, "HU12", 1.75), (10, "HU12", 2.00), (11, "EC40", 0.50), (12, "HU12", 2.25),
    (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00), (16, "HU12", 3.50),
    (17, "HU12", 4.00), (18, "HU12", 4.50), (19, "HU12", 5.00), (20, "HU12", 5.50),
    (21, "HU12", 6.00), (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75),
    (25, "FR76", 2.00),
]


def read_table_duration(path):
    last_t, n = 0.0, 0
    with open(path) as f:
        f.readline(); f.readline()
        for line in f:
            q = line.split()
            if len(q) >= 2:
                try:
                    last_t = float(q[0]); n += 1
                except ValueError:
                    pass
    return last_t, n


def build_protocol():
    dur = {}
    for rec, (name, fname) in BASE_TABLES.items():
        pth = os.path.join(TABLE_DIR, fname)
        if os.path.isfile(pth):
            dur[rec] = read_table_duration(pth)
    rows = []
    for rn, rec, scale in PROTOCOL_TABLE2:
        if rn < RUN_FROM or rn > RUN_TO:
            continue
        name, fname = BASE_TABLES[rec]
        d, n = dur.get(rec, (0.0, 0))
        rows.append({"run": rn, "record": rec, "scale": scale,
                     "table_name": name, "table_file": fname,
                     "dur_s": d, "n_samples": n})
    if not rows:
        raise RuntimeError("no runs in [{}, {}]".format(RUN_FROM, RUN_TO))
    return rows


def import_base_tables():
    for rec, (name, fname) in BASE_TABLES.items():
        pth = os.path.join(TABLE_DIR, fname)
        if os.path.isfile(pth):
            it.command("table '{}' import '{}'".format(name, cmd_path(pth)))
            print("  imported base table '{}' <- {}".format(name, fname))
        else:
            print("  ! base table missing: {}".format(pth))


PROTOCOL = build_protocol()

FISH_HISTORIES = {
    1: "Record_Disp", 2: "Bot_Quarter_Disp", 3: "Mid_Disp",
    4: "Top_Quarter_A_Disp", 5: "Top_Quarter_B_Disp", 6: "cstav",
    7: "joist_s1_shear", 8: "joist_s2_shear", 9: "total_shear", 10: "tilt_angles",
    11: "tilt_bot_seg", 12: "tilt_low_seg", 13: "tilt_up_seg",
    14: "tilt_beam_seg", 15: "tilt_full_wall", 16: "rel_disp_top_mm",
}

# =====================================================================
# 3.  PATH HELPERS
# =====================================================================
def cmd_path(path):
    return os.path.normpath(path).replace("\\", "/")

def save_file_path(run_no):
    return os.path.join(OUT_DIR, "stratC_run_{:02d}.sav".format(run_no))

def _find_channel_csv(folder, needle="channel_19"):
    if not os.path.isdir(folder):
        return None
    hits = [os.path.join(folder, fn) for fn in os.listdir(folder)
            if fn.lower().endswith(".csv") and needle in fn.lower()]
    return sorted(hits)[0] if hits else None

def run_label(run_no, record, scale):
    return "Run{:02d}_{}_s{}".format(
        run_no, record, "{:.2f}".format(scale).replace(".", "p"))

def run_dir(run_no, record=None, scale=None):
    if record is not None and scale is not None:
        d = os.path.join(OUT_DIR, run_label(run_no, record, scale))
        if os.path.isdir(d):
            return d
    prefix = "Run{:02d}_".format(run_no)
    if os.path.isdir(OUT_DIR):
        for d in sorted(os.listdir(OUT_DIR)):
            if d.startswith(prefix) and os.path.isdir(os.path.join(OUT_DIR, d)):
                return os.path.join(OUT_DIR, d)
    return OUT_DIR

def ch19_csv_path(run_no, record=None, scale=None):
    if record is not None and scale is not None:
        scale_str = "{:.2f}".format(scale).replace(".", "p")
        run_folder = os.path.join(OUT_DIR, "Run{:02d}_{}_s{}".format(run_no, record, scale_str))
        hit = _find_channel_csv(run_folder)
        if hit:
            return hit
    prefix = "Run{:02d}_".format(run_no)
    if os.path.isdir(OUT_DIR):
        for d in sorted(os.listdir(OUT_DIR)):
            if d.startswith(prefix) and os.path.isdir(os.path.join(OUT_DIR, d)):
                hit = _find_channel_csv(os.path.join(OUT_DIR, d))
                if hit:
                    return hit
    return os.path.join(OUT_DIR,
        "Channel_19_DispTopQRight_run{:02d}.csv".format(run_no))

def run_folder_path(run_no, record, scale):
    scale_str = "{:.2f}".format(scale).replace(".", "p")
    return os.path.join(OUT_DIR,
                        "Run{:02d}_{}_s{}".format(run_no, record, scale_str))


def _read_hist(path):
    try:
        a = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None
    if a.ndim < 2 or a.shape[1] < 2:
        return None
    m = np.isfinite(a[:, 0]) & np.isfinite(a[:, 1])
    a = a[m]
    return a if len(a) >= 5 else None


def secant_period(run_no, record, scale):
    """Retired in Strategy F; kept for diff-ability with the pulse drivers."""
    info = {"V_peak_kN": None, "disp_peak_mm": None, "k_sec_kNpm": None}
    folder = run_folder_path(run_no, record, scale)
    fv = _find_channel_csv(folder, "total_shear")
    fd = _find_channel_csv(folder, "rel_disp_top_mm")
    if not (fv and fd):
        return None, info
    av, ad = _read_hist(fv), _read_hist(fd)
    if av is None or ad is None:
        return None, info
    v = av[:, 1]; u = ad[:, 1]
    V_peak = float(np.max(np.abs(v - v[0])))
    d_peak = float(np.max(np.abs(u - u[0]))) / 1000.0
    info["V_peak_kN"] = round(V_peak, 4)
    info["disp_peak_mm"] = round(d_peak * 1000.0, 4)
    if d_peak <= 0 or V_peak <= 0:
        return None, info
    k_sec = V_peak / d_peak
    info["k_sec_kNpm"] = round(k_sec, 2)
    return 2.0 * math.pi * math.sqrt(M_EFF_T / k_sec), info


def state_file_path():
    return os.path.join(OUT_DIR, STATE_FILE_NAME)

# =====================================================================
# 4.  CHECKPOINT / RESUME
# =====================================================================
def save_checkpoint(run_no, T_current, summary):
    state = {"last_completed_run": run_no, "T_current": T_current, "summary": summary,
             "case": CASE_ID, "start_save": START_SAVE,
             "large_strain": LARGE_STRAIN, "maxwell": MAXWELL_CMD}
    with open(state_file_path(), "w", newline="\n") as f:
        json.dump(state, f, indent=2)

def load_checkpoint():
    sf = state_file_path()
    if os.path.exists(sf):
        try:
            with open(sf, "r") as f:
                state = json.load(f)
            if state.get("case", CASE_ID) != CASE_ID or state.get("maxwell", MAXWELL_CMD) != MAXWELL_CMD \
                    or state.get("large_strain", LARGE_STRAIN) != LARGE_STRAIN:
                raise RuntimeError("checkpoint in {} belongs to a different case / damping / geometry; "
                                   "use a new OUT_DIR".format(OUT_DIR))
            last = state["last_completed_run"]
            T_cur = state["T_current"]
            summary = state.get("summary", [])
            print("  Checkpoint JSON found: last completed = run {:02d}, "
                  "T_current = {:.4f} s".format(last, T_cur))
            return last + 1, T_cur, summary
        except RuntimeError:
            raise
        except Exception as e:
            print("  WARNING: checkpoint corrupt ({}), scanning save files...".format(e))
    last_found = 0
    for run_no in range(len(PROTOCOL), 0, -1):
        if os.path.exists(save_file_path(run_no)):
            last_found = run_no
            break
    if last_found == 0:
        print("  No prior runs found. Starting fresh.")
        return 1, T1_init, []
    ch19_path = ch19_csv_path(last_found)
    if os.path.exists(ch19_path):
        sine_dur_approx = n_cycles * T1_init
        T_recovered = identify_Tend_from_csv(ch19_path, sine_dur_approx)
        print("  Recovered: last run = {:02d}, T_end = {:.4f} s".format(last_found, T_recovered))
        return last_found + 1, T_recovered, []
    if last_found > 1:
        prev_ch19 = ch19_csv_path(last_found - 1)
        if os.path.exists(prev_ch19):
            sine_dur_approx = n_cycles * T1_init
            T_recovered = identify_Tend_from_csv(prev_ch19, sine_dur_approx)
            print("  Run {:02d} CSV missing; recovered T from run {:02d}: {:.4f} s".format(
                  last_found, last_found-1, T_recovered))
            return last_found, T_recovered, []
    print("  WARNING: cannot recover T_current. Using T1_init.")
    return last_found, T1_init, []

def check_run_complete(run_no):
    return (os.path.exists(save_file_path(run_no)) and os.path.exists(ch19_csv_path(run_no)))

# =====================================================================
# 5.  LOAD PRECOMPUTED Sd SPECTRA  (unused in Strategy F; kept for diff-ability)
# =====================================================================
SPECTRA = {}
def load_spectra():
    for name in ["HU12", "EC40", "FR76"]:
        fpath = "spectrum_{}.csv".format(name)
        data = np.genfromtxt(fpath, delimiter=',', skip_header=1)
        SPECTRA[name] = (data[:, 0], data[:, 1])
        print("  Loaded {}: {} periods".format(fpath, len(data)))

def interpolate_sd(record, T):
    T_arr, Sd_arr = SPECTRA[record]
    return float(np.interp(T, T_arr, Sd_arr))

# =====================================================================
# 6.  SINUSOID GENERATION  (unused in Strategy F; kept for diff-ability)
# =====================================================================
def newmark_sd(a_g, dt, T, xi):
    if T <= 0:
        return 0.0
    w = 2*math.pi/T; c = 2*xi*w; w2 = w*w
    u = 0.0; v = 0.0; acc = -a_g[0]; u_max = 0.0
    for i in range(1, len(a_g)):
        u_new = u + dt*v + dt**2 * 0.25 * acc
        v_new = v + dt * 0.5 * acc
        denom = 1.0 + 0.5*dt*c + 0.25*dt**2*w2
        acc_new = (-a_g[i] - c*v_new - w2*u_new) / denom
        u_new += 0.25*dt**2*acc_new
        v_new += 0.5*dt*acc_new
        u = u_new; v = v_new; acc = acc_new
        if abs(u) > u_max:
            u_max = abs(u)
    return u_max

def calibrate_amplitude(T, Sd_target):
    A_trial = 1.0
    w = 2.0 * math.pi / T
    n = int(round(n_cycles * T / delta_t)) + 1
    t_arr = np.arange(n) * delta_t
    a0 = A_trial * np.cos(w * t_arr)
    sd0 = newmark_sd(a0, delta_t, T, xi)
    if sd0 == 0.0:
        raise RuntimeError("Trial Sd=0 at T={:.4f}".format(T))
    return (Sd_target / sd0) * A_trial

def pulse_arrays(A, T):
    w = 2.0 * math.pi / T
    n = int(round(n_cycles * T / delta_t)) + 1
    t = np.arange(n) * delta_t
    v = (A / w) * np.sin(w * t)
    ramp_sec = max(0.5 * T, delta_t)
    n_tail = int(round(tail_sec / delta_t))
    n_ramp = min(int(round(ramp_sec / delta_t)), n_tail)
    t_tail = t[-1] + delta_t + np.arange(n_tail) * delta_t
    t = np.r_[t, t_tail]
    v_end = float(v[-1])
    if n_ramp > 1:
        s = np.linspace(0.0, 1.0, n_ramp)
        v_ramp = v_end * 0.5 * (1.0 + np.cos(math.pi * s))
    else:
        v_ramp = np.array([0.0])
    v_tail = np.zeros(n_tail)
    v_tail[:len(v_ramp)] = v_ramp
    return t, np.r_[v, v_tail]

def table_displacement(t, v):
    return np.concatenate(([0.0],
                           np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(t))))

def calibrate_amplitude_pgd(T, PGD_target):
    t, v = pulse_arrays(1.0, T)
    d = table_displacement(t, v)
    amp = 0.5 * (float(np.max(d)) - float(np.min(d)))
    if amp <= 0:
        raise RuntimeError("zero trial table displacement at T={:.4f}".format(T))
    return PGD_target / amp

def build_velocity_file(A, T, run_no, out_dir):
    t, v = pulse_arrays(A, T)
    fname = "vel_run_{:02d}.txt".format(run_no)
    fpath = os.path.join(out_dir, fname)
    N = len(t)
    with open(fpath, "w", newline="\n") as f:
        f.write("StratD_run{:02d}_T{:.4f}\n".format(run_no, T))
        f.write("{}\t0\n".format(N))
        for ti, vi in zip(t, v):
            f.write("{:.6f}\t{:.9e}\n".format(ti, vi))
    duration = float(t[-1])
    v_peak = float(np.max(np.abs(v)))
    return fpath, duration, v_peak

def generate_pulse(record, T, Sd_target, run_no, out_dir):
    if not USE_RATCHETING:
        A = calibrate_amplitude(T, Sd_target)
        vp, dur, vpk = build_velocity_file(A, T, run_no, out_dir)
        return vp, dur, vpk, {"A": A, "beta": 1.0, "s": 0}
    rec = RECORD_ASYM[record]
    s = rec["s"]
    beta_target = rec["beta"] ** ASYM_K
    beta_now = beta_eff(T, T1_init, beta_target)
    A, _ = calibrate_amplitude_asym(T, Sd_target, beta_now, s,
                                    n_cycles, delta_t, xi, k=1.0)
    vp, dur, vpk = build_velocity_asym(A, T, beta_now, s, run_no, out_dir,
                                       n_cycles, delta_t, tail_sec)
    return vp, dur, vpk, {"A": A, "beta": beta_now, "s": s}

# =====================================================================
# 7.  PERIOD IDENTIFICATION (segment + window + autocorrelation)
# =====================================================================
RMS_ALIVE_MM  = 0.05
RMS_CHUNK_S   = 0.15
MAX_RD_WINDOW = 0.50

def extract_last_segment(t_all, u_all):
    dt = np.diff(t_all)
    neg_idx = np.where(dt < -0.5)[0]
    n_segments = len(neg_idx) + 1
    last_start = neg_idx[-1] + 1 if len(neg_idx) > 0 else 0
    return t_all[last_start:], u_all[last_start:], n_segments

def find_live_window(t, u, sine_end_time):
    mask = t >= sine_end_time
    t_rd = t[mask]; u_rd = u[mask]
    if len(u_rd) < 16:
        return t_rd, u_rd, 0.0
    dt_hist = float(np.median(np.diff(t_rd)))
    if dt_hist <= 0:
        dt_hist = delta_t
    u_dm = u_rd - np.mean(u_rd)
    chunk_n = max(1, int(RMS_CHUNK_S / dt_hist))
    cutoff_idx = len(u_dm)
    for start in range(0, len(u_dm) - chunk_n, chunk_n):
        chunk = u_dm[start:start + chunk_n]
        rms_mm = math.sqrt(float(np.mean(chunk**2))) * 1000
        if rms_mm < RMS_ALIVE_MM:
            cutoff_idx = start
            break
    max_idx = min(cutoff_idx, int(MAX_RD_WINDOW / dt_hist))
    max_idx = max(max_idx, int(0.2 / dt_hist))
    t_live = t_rd[:max_idx]; u_live = u_rd[:max_idx]
    window_dur = float(t_live[-1] - t_live[0]) if len(t_live) > 1 else 0
    return t_live, u_live, window_dur

def identify_Tend_from_csv(ch19_csv, sine_end_time):
    try:
        try:
            data = np.genfromtxt(ch19_csv, delimiter=',', skip_header=1)
            if data.ndim == 1 or data.shape[1] < 2:
                raise ValueError
        except (ValueError, IndexError):
            data = np.genfromtxt(ch19_csv, skip_header=2)
    except Exception as e:
        print("  WARNING: cannot read {}: {}".format(ch19_csv, e))
        return T1_init
    t_all = data[:, 0]; u_all = data[:, 1]
    t_seg, u_seg, n_segments = extract_last_segment(t_all, u_all)
    if n_segments > 1:
        print("    CSV has {} segments, using last ({} pts)".format(n_segments, len(t_seg)))
    t_live, u_live, win_dur = find_live_window(t_seg, u_seg, sine_end_time)
    if len(u_live) < 16:
        print("  WARNING: live ring-down too short ({} pts)".format(len(u_live)))
        return T1_init
    dt_hist = float(np.median(np.diff(t_live)))
    if dt_hist <= 0:
        dt_hist = delta_t
    u_dm = u_live - np.mean(u_live)
    u_win = u_dm * np.hanning(len(u_dm))
    N_orig = len(u_win)
    N_padded = N_orig * 16
    freqs = np.fft.rfftfreq(N_padded, d=dt_hist)
    amps = np.abs(np.fft.rfft(u_win, n=N_padded))
    band_mask = (freqs >= FFT_F_MIN) & (freqs <= FFT_F_MAX)
    T_fft = T1_init
    if np.any(band_mask):
        band_idx = np.where(band_mask)[0]
        i_peak = band_idx[np.argmax(amps[band_idx])]
        f_peak = float(freqs[i_peak])
        if i_peak > 0 and i_peak < len(amps) - 1:
            alpha = float(np.log(amps[i_peak - 1] + 1e-30))
            beta  = float(np.log(amps[i_peak]     + 1e-30))
            gamma = float(np.log(amps[i_peak + 1] + 1e-30))
            denom = alpha - 2*beta + gamma
            if abs(denom) > 1e-12:
                df = float(freqs[1] - freqs[0]) if len(freqs) > 1 else 1.0
                f_peak = f_peak + 0.5 * (alpha - gamma) / denom * df
        T_fft = 1.0 / f_peak if f_peak > 0 else T1_init
    u_ac = u_dm / (np.max(np.abs(u_dm)) + 1e-30)
    acorr = np.correlate(u_ac, u_ac, mode='full')
    acorr = acorr[len(u_ac)-1:]
    acorr = acorr / (acorr[0] + 1e-30)
    lag_min = max(1, int(T_MIN_PHYS / dt_hist))
    lag_max = min(len(acorr) - 2, int(T_MAX_PHYS / dt_hist))
    T_acorr = T1_init
    if lag_max > lag_min + 2:
        peaks = []
        for j in range(lag_min + 1, lag_max):
            if acorr[j] > acorr[j-1] and acorr[j] > acorr[j+1]:
                peaks.append((j, acorr[j]))
        if peaks:
            best_lag, _ = max(peaks, key=lambda x: x[1])
            if best_lag > 0 and best_lag < len(acorr) - 1:
                y0 = float(acorr[best_lag - 1]); y1 = float(acorr[best_lag]); y2 = float(acorr[best_lag + 1])
                d2 = y0 - 2*y1 + y2
                lag_ref = best_lag + 0.5*(y0 - y2)/d2 if abs(d2) > 1e-12 else float(best_lag)
            else:
                lag_ref = float(best_lag)
            T_acorr = lag_ref * dt_hist
    if T_acorr > T_MIN_PHYS and T_acorr < T_MAX_PHYS:
        T_end = T_acorr; method = "autocorr"
    else:
        T_end = T_fft; method = "FFT"
    T_end = max(T_MIN_PHYS, min(T_MAX_PHYS, T_end))
    print("  Period ID: T_acorr={:.5f}s, T_fft={:.5f}s -> T_end={:.5f}s ({})".format(
        T_acorr, T_fft, T_end, method))
    print("    (live window={:.3f}s, N_live={}, segments={})".format(win_dur, len(u_live), n_segments))
    return T_end

# =====================================================================
# 8.  EXPORT HISTORIES
# =====================================================================
def export_all_histories(run_no, record, scale, out_dir):
    scale_str = "{:.2f}".format(scale).replace(".", "p")
    run_label = "Run{:02d}_{}_s{}".format(run_no, record, scale_str)
    full_folder = os.path.join(out_dir, run_label)
    os.makedirs(full_folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(full_folder)))
    it.command("[runlabel='{}']".format(run_label))
    it.command("call 'instrument_history_export_v2.dat'")

# =====================================================================
# 9.  MODEL SETUP (fresh or resumed)
# =====================================================================
def setup_model_for_dynamic(save_file, verify_build=False):
    it.command("model restore '{}'".format(cmd_path(save_file)))
    if str(LARGE_STRAIN).lower() in ("on", "true", "1"):
        it.command("model large-strain on")
        print("  GEOMETRY: large-strain ON (positions updated, contacts "
              "re-detected)")
        print("           NB contacts CREATED during the run inherit the "
              "material-table default (tension 0, cohesion 0, friction), which "
              "is the intended cracked-joint behaviour on this build.")
    else:
        it.command("model large-strain off")
        print("  GEOMETRY: small-strain (default)")
    it.command("model dynamic active on")
    it.command("""
    block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
    block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
    """)
    it.command("call 'instrument_history_new.dat'")
    it.command("call 'instrument_baseshear_exp.dat'")
    it.command("""
    block free velocity-z range group 'S'
    """)
    it.command("""
    block free rotation-y range group 'T_B'
    block free rotation-z range group 'T_B'
    block mech damp local 0.0
    block mech damp global 0.0
    """)
    if MAXWELL_CMD:
        it.command(MAXWELL_CMD)
        print("  Maxwell damping applied: {}".format(MAXWELL_CMD))
    if DAMP_RATIO > 0:
        cmd = "block mechanical damping rayleigh {ratio} {freq} {dtype}".format(
            ratio=DAMP_RATIO, freq=1.0 / T1_init, dtype=DAMP_TYPE).strip()
        it.command(cmd)
        print("  Rayleigh damping applied: {}".format(cmd))
    elif not MAXWELL_CMD:
        print("  (no damping applied; contact dissipation only)")

    if G_I_NEW is not None:
        it.command("block contact property G_I {:g} G_II {:g}".format(
            G_I_NEW, G_II_RATIO * G_I_NEW))
        print("  G_I -> {:g} J/m2 , G_II -> {:g} J/m2".format(G_I_NEW, G_II_RATIO * G_I_NEW))
    if COH_RESIDUAL:
        it.command("block contact property cohesion-residual {:g}".format(COH_RESIDUAL))
        print("  cohesion-residual -> {:g} Pa".format(COH_RESIDUAL))

    # Build guard: the control must not run on the bonded (f_t = 0.344) build.
    if verify_build:
        try:
            _ft = float(it.fish.get("tension"))
        except Exception as _e:
            _ft = None
            print("  WARNING: cannot read FISH 'tension' from the save ({}); "
                  "build NOT verified.".format(_e))
        if _ft is not None:
            print("  joint tension in this save: {:.3f} MPa".format(_ft / 1e6))
            if _ft > TENSION_MAX_PA:
                raise RuntimeError(
                    "save has joint tension {:.3f} MPa > {:.3f} MPa: this is the bonded "
                    "build, not the f_t ~ 0 control build. Fix START_SAVE.".format(
                        _ft / 1e6, TENSION_MAX_PA / 1e6))
        try:
            _fapp = 0.0
            it.command("[global _fapp_tb = 0.0]")
            it.command("""
fish define _sum_fapp
    _fapp_tb = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'T_B')
            _fapp_tb = _fapp_tb + comp.y(block.force.app(b))
        endif
    endloop
end
[_sum_fapp]
""")
            _fapp = float(it.fish.get("_fapp_tb"))
            print("  applied spring force on T_B: {:.0f} N (expected -10350)".format(_fapp))
            if abs(_fapp + 10350.0) > 100.0:
                raise RuntimeError("applied spring force on T_B is {:.0f} N, not -10350: "
                                   "wrong build or force lost in the save.".format(_fapp))
        except RuntimeError:
            raise
        except Exception as _e:
            print("  WARNING: cannot read the applied force on T_B ({}).".format(_e))
    print("  Model setup complete (BCs + damping + joist contact groups applied).")

# =====================================================================
# 10.  EXECUTE ONE RUN
# =====================================================================
def execute_run(entry):
    run_no = entry["run"]
    record = entry["record"]
    scale  = entry["scale"]
    name   = entry["table_name"]
    dur    = entry["dur_s"]

    print("\n" + "=" * 70)
    print("  Run {:02d}: {} x {:.2f}   table '{}'   duration = {:.2f} s "
          "(+{:.1f} s ring-down)".format(
              run_no, record, scale, name, dur, TAIL_SEC_RECORD))
    print("=" * 70)

    it.command("model dynamic time-total 0")
    for grp in ["S", "T_B"]:
        it.command("block apply velocity-z {:g} table '{}' range group '{}'"
                   .format(scale, name, grp))
    it.command("model solve dynamic time {:.6f}".format(dur))
    for grp in ["S", "T_B"]:
        it.command("block gridpoint apply-remove velocity-z range group '{}'"
                   .format(grp))
    if TAIL_SEC_RECORD > 0:
        it.command("model solve dynamic time {:.6f}".format(TAIL_SEC_RECORD))

    it.command("model save '{}'".format(cmd_path(save_file_path(run_no))))
    export_all_histories(run_no, record, scale, OUT_DIR)

    # ---- peak / residual displacement and like-for-like base shear ----
    _folder = run_folder_path(run_no, record, scale)
    _peak_mm = _res_mm = _vexp_kN = float("nan")
    _fd = _find_channel_csv(_folder, "rel_disp_top_exp_mm")
    _ad = _read_hist(_fd) if _fd else None
    if _ad is not None:
        _u = _ad[:, 1]
        _peak_mm = float(np.max(np.abs(_u - _u[0])))
        _n = max(1, int(0.02 * len(_u)))
        _res_mm = float(np.mean(_u[-_n:]) - _u[0])
    _fb = _find_channel_csv(_folder, "base_shear_exp_kn")
    _ab = _read_hist(_fb) if _fb else None
    if _ab is not None:
        _vexp_kN = float(np.max(np.abs(_ab[:, 1] - _ab[0, 1])))

    # ---- period identification (damage measure only) ----
    global _T_TRACK
    T_legacy = identify_Tend_from_csv(
        ch19_csv_path(run_no, record, scale), dur)
    _exp = {"T": float("nan"), "ok": False, "spread": float("nan"),
            "window_s": 0.0, "T_welch": float("nan"), "note": "module missing"}
    if _HAVE_PERIOD_ID_EXP:
        _exp = period_id_exp.identify_period_experiment(
            run_dir(run_no, record, scale), run_label(run_no, record, scale),
            pulse_end_s=dur)
        _cl = _exp.get("candidates") or []
        if _cl:
            if _T_TRACK is None or not np.isfinite(_T_TRACK):
                _T_sel = period_id_exp.seed_period(_cl); _sel_rule = "seed"
            else:
                _T_sel, _rec = period_id_exp.track_step(_cl, _T_TRACK)
                _sel_rule = "track"
            if np.isfinite(_T_sel):
                _T_TRACK = _T_sel
                _exp["T"] = _T_sel
                _exp["ok"] = True
                _exp["note"] = (_exp["note"] + " | " + _sel_rule).strip(" |")

    if PERIOD_ID_METHOD == "experiment" and _exp["ok"]:
        T_end = _exp["T"]; _tid_src = "exp_psd_svd"
    elif PERIOD_ID_METHOD == "experiment":
        T_end = T_legacy
        _tid_src = "legacy_fallback:" + str(_exp["note"])
        print("  ** experiment-method identification rejected ({}). "
              "Using legacy {:.5f} s.".format(_exp["note"], T_legacy))
    else:
        T_end = T_legacy; _tid_src = "legacy"

    print("  peak rel. OOP {:.3f} mm   residual {:+.3f} mm   peak sum(m a) {:.2f} kN".format(
        _peak_mm, _res_mm, _vexp_kN))
    print("  T_end = {:.4f} s ({:.2f}x T1_init)   [damage measure, {}]".format(
        T_end, T_end / T1_init, _tid_src))
    if _HAVE_PERIOD_ID_EXP and np.isfinite(_exp["T"]):
        print("     legacy(disp,1ch) = {:.5f} s   experiment(acc,3ch) = {:.5f} s"
              .format(T_legacy, _exp["T"]))

    # NB: do NOT delete the base table -- it is reused every run.
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")
    it.command("call 'instrument_baseshear_exp.dat'")

    return {
        "run": run_no, "record": record, "scale": round(scale, 4),
        "table_file": entry["table_file"],
        "record_dur_s": round(dur, 4),
        "peak_rel_OOP_mm": round(_peak_mm, 4) if np.isfinite(_peak_mm) else "",
        "residual_OOP_mm": round(_res_mm, 4) if np.isfinite(_res_mm) else "",
        "peak_baseshear_exp_kN": round(_vexp_kN, 3) if np.isfinite(_vexp_kN) else "",
        "T_end": round(T_end, 6),
        "T_end_over_Tinit": round(T_end / T1_init, 4),
        "T_legacy_disp": round(T_legacy, 6),
        "T_exp_psd": round(_exp["T"], 6) if np.isfinite(_exp["T"]) else "",
        "T_exp_welch": round(_exp["T_welch"], 6) if np.isfinite(_exp["T_welch"]) else "",
        "T_exp_ch_spread": round(_exp["spread"], 4) if np.isfinite(_exp["spread"]) else "",
        "T_exp_window_s": round(_exp["window_s"], 4),
        "tid_source": _tid_src,
    }

# =====================================================================
# 11.  MAIN DRIVER
# =====================================================================
SPECTRA_FILES = [os.path.join(TABLE_DIR, f) for (_n, f) in BASE_TABLES.values()]
DAT_FILES     = ["instrument_history_new.dat",
                 "instrument_baseshear_exp.dat",
                 "instrument_history_export_v2.dat"]
BASE_SAVE     = START_SAVE

def preflight_checks():
    print("\n" + "=" * 70)
    print("  PRE-FLIGHT CHECKS  (working dir: {})".format(os.getcwd()))
    print("=" * 70)
    fresh = (not os.path.exists(state_file_path()) and
             not any(os.path.exists(save_file_path(r)) for r in range(1, len(PROTOCOL) + 1)))
    required = list(SPECTRA_FILES) + list(DAT_FILES) + ([BASE_SAVE] if fresh else [])
    missing = []
    for f in required:
        ok = os.path.isfile(f)
        if not ok:
            missing.append(f)
        print("   [{:^7}] {}".format("OK" if ok else "MISSING", f))
    print("   start mode  : {}".format("FRESH" if fresh else "RESUME"))
    if MAXWELL_CMD:
        print("   damping     : Maxwell  [{}]  (three units at 0.0085 = 1.5 % measured)".format(MAXWELL_CMD))
    if DAMP_RATIO > 0:
        print("   damping     : Rayleigh {:.1f}% at {:.2f} Hz  [{}]".format(
            DAMP_RATIO * 100.0, 1.0 / T1_init, DAMP_TYPE))
    if not MAXWELL_CMD and DAMP_RATIO <= 0:
        print("   damping     : NONE (contact dissipation only)")
    print("   geometry    : large-strain {}".format(LARGE_STRAIN))
    print("   case        : {}".format(CASE_ID))
    print("   T1_init     : {:.4f} s (period-ID windows only)".format(T1_init))
    print("   output dir  : {}".format(OUT_DIR))
    if missing:
        raise RuntimeError("Missing required input file(s) in the working "
                           "directory: {}.".format(", ".join(missing)))
    print("   all required inputs found.")
    print("=" * 70)

def run_strategy_C():
    os.makedirs(OUT_DIR, exist_ok=True)
    preflight_checks()
    print("\n--- Checking for existing runs ---")
    resume_from, T_current, summary = load_checkpoint()
    if resume_from > len(PROTOCOL):
        print("\nAll {} runs already completed.".format(len(PROTOCOL)))
        print("Delete '{}' to restart.".format(state_file_path()))
        return summary
    _tot = sum(e["dur_s"] for e in PROTOCOL) + len(PROTOCOL) * TAIL_SEC_RECORD
    print("\n  Strategy F: {} records, {:.0f} s of model time total.".format(
        len(PROTOCOL), _tot))
    print("  Large strain re-detects contacts every step; read the wall-clock of "
          "run 1 before extrapolating (the old ~230 s wall per s model was small strain).")
    if resume_from == 1:
        print("\n--- Starting fresh from {} ---".format(START_SAVE))
        setup_model_for_dynamic(START_SAVE, verify_build=True)
    else:
        last_done = resume_from - 1
        print("\n--- Resuming: restoring run {:02d}, will execute run {:02d} next ---".format(
            last_done, resume_from))
        print("  T_current = {:.4f} s ({:.2f}x T_init)".format(T_current, T_current / T1_init))
        setup_model_for_dynamic(save_file_path(last_done))

    import_base_tables()

    log_path = os.path.join(OUT_DIR, "strategy_C_log.csv")
    log_is_new = (resume_from == 1) or not os.path.exists(log_path)
    log_f = open(log_path, "w" if log_is_new else "a", newline="\n")
    LOG_COLS = ["run", "record", "scale", "table_file", "record_dur_s",
                "peak_rel_OOP_mm", "residual_OOP_mm", "peak_baseshear_exp_kN",
                "T_end", "T_end_over_Tinit"]
    if log_is_new:
        log_f.write("# case {} | start {} | large-strain {} | {} | T1_init {}\n".format(
            CASE_ID, START_SAVE, LARGE_STRAIN, MAXWELL_CMD or "undamped", T1_init))
        log_f.write(",".join(LOG_COLS) + "\n")

    for idx in range(resume_from - 1, len(PROTOCOL)):
        entry = PROTOCOL[idx]
        run_summary = execute_run(entry)
        summary.append(run_summary)
        s = run_summary
        log_f.write(",".join(str(s.get(c, "")) for c in LOG_COLS) + "\n")
        log_f.flush()
        save_checkpoint(entry["run"], run_summary["T_end"], summary)
        print("  Checkpoint saved after run {:02d}.".format(entry["run"]))
    log_f.close()

    print("\n" + "=" * 70)
    print("Strategy F (full records) complete: {} runs   [{}]".format(len(PROTOCOL), CASE_ID))
    print("=" * 70)
    print("  {:>4} {:>6} {:>7} {:>9} {:>10} {:>9}".format(
        "run", "rec", "scale", "dur_s", "peak_mm", "T_end_s"))
    for s in summary:
        print("  {:>4} {:>6} {:>7.2f} {:>9.2f} {:>10} {:>9.4f}".format(
            s["run"], s["record"], s["scale"], s["record_dur_s"], s["peak_rel_OOP_mm"], s["T_end"]))

    csv_path = os.path.join(OUT_DIR, "strategy_C_summary.csv")
    if summary:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=summary[0].keys())
            w.writeheader()
            w.writerows(summary)
        print("Summary -> {}".format(csv_path))
    return summary

# =====================================================================
# ENTRY POINT
# =====================================================================
run_strategy_C()