# -*- coding: ascii -*-
"""
Strategy F (FULL) -- significant-duration-truncated records, runs 1-25.

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
    PROTOCOL_TABLE2. The velocity-z coefficient multiplies the base table, so
    the record shape is fixed and only its amplitude scales run to run.

TRUNCATION (done upstream in integrate_table_accel.py)
    Each base table is the channel-12 acceleration integrated to velocity, kept
    from t = 0 and cut at the velocity zero-crossing nearest its 95% Arias time,
    so the applied motion starts and ends at rest with no velocity step.

WHY IT EXISTS
    The diagnostic profile figure showed the model reaching about a sixth of
    the experimental out-of-plane displacement while reproducing the mechanism
    shape. Two causes are possible: intensity starvation (a loading problem the
    pulse strategies can fix) or the model itself under-responding. The single
    thing that separates them is applying the REAL motion and seeing whether
    the wall then reaches the experimental displacement. This driver does that.

    Run it against a pulse strategy at the same damage state:
      - if the record reaches the experimental displacement and the pulse does
        not, the gap is the pulse's missing cycles (a loading limitation), and
        the fidelity-vs-cost trade is real and quantified.
      - if the record ALSO falls short, the gap is in the model (boundary
        conditions, large-strain, damping), and no loading change fixes it.

TWO WAYS TO RUN IT
    Full cumulative sequence: START_SAVE at the base save, RUN_FROM/RUN_TO
        spanning all runs. About 180 s of model time for US-1's 24 records,
        roughly 11 h at ~230 s wall per s model. This is the complete fidelity
        model, directly comparable run-for-run to the pulse strategies.
    Targeted benchmark (cheaper, and the decisive test): START_SAVE at a
        damaged save from a pulse run, RUN_FROM/RUN_TO at just the FR76 runs
        22-24. One to three records, ~1 h, answers the loading-vs-model
        question without the full sequence.

    NB a warm start from another strategy's save inherits that strategy's
    damage history. For a clean run-for-run comparison, start both the record
    and the pulse analysis from the same state.

INPUTS
    vel_HU.txt, vel_EC.txt, vel_FR.txt in TABLE_DIR (the 100%-intensity base
    tables), written by integrate_table_accel.py. Generate them first:
        python integrate_table_accel.py

REPORTED
    T_end from the ring-down after each record, as a damage measure. It drives
    nothing. Run folders are Run{NN}_{record}_s{scale}, so postprocess_stratC.py
    parses them unchanged and the comparison figures work as-is.

LIMITATIONS
    Same model caveats as the pulse strategies: the top beam is driven with the
    base motion, and large-strain is off, so a large-displacement or collapse
    result needs large-strain on to be trusted. What this driver removes is the
    loading approximation, not the model approximations.

GENERATED FROM strategy_F_3dec.py: only the config block (manifest name, table
prefix, RUN_FROM/RUN_TO = 1..25, OUT_DIR) differs. The pulse-synthesis helpers
(pulse_arrays, calibrate_amplitude_pgd, newmark_sd, ...) remain in the file but
are unused here. Preserve the 3DEC-embedded Python style. run_strategy_C() is
kept as the entry-point name only, for diff-ability against the siblings.

OUTPUT
    stratF_full_results_US{n}/
"""

import itasca as it
import numpy as np
import os, csv, math, json, sys

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
T1_init     = 0.092
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
USE_RATCHETING = False   # OFF for this variant: symmetric pulse, ASYM_K inert
ASYM_K         = 1.0     # global asymmetry sharpening (>=1 sharpens); calibrate out-of-sample
# Damping is exposed here so it is explicit and printed at startup.
# The conference paper used FULL Rayleigh at 3%. The pasted driver used
# mass-proportional 0.06, which over-damps the slow rocking and is the wrong
# choice for a collapse study. Default below is full Rayleigh 3%; change if needed.
DAMP_RATIO = 0.0         # fraction of critical at the centre frequency
DAMP_TYPE  = ""          # "" = full Rayleigh (mass+stiffness); "mass"; "stiffness"

# >>> STRATEGY D : the protocol is a table, not a feedback loop.
#
# Each run is defined by two numbers measured from the corresponding
# shake-table record, and by nothing about the model:
#
#   PGD_amp_mm  amplitude of the table displacement (half its peak-to-peak)
#   T_eq_s      equivalent period, pi * PGD / PGV
#
# Both are properties of the ground motion, so the intensity measure is
# record-invariant. That is the defect that sank the Strategy C variants:
# Sd(record, T_current) is a function of the damaged structure, so it created
# an intensity-damage feedback loop with no stable regime, and its values were
# not comparable between runs or against any other study.
#
# The pulse amplitude is calibrated so the applied table displacement equals
# PGD_amp_mm exactly. The pulse period is T_eq_s. Nothing here adapts.
# T_end is still identified from the ring-down and logged, purely as a damage
# measure. It does not influence the excitation.
# >>> STRATEGY F (FULL) : significant-duration-truncated records, runs 1-25.
RECORD_TEST   = 9          # only used for the output-folder suffix
TABLE_DIR     = "velocity_output"        # folder holding vel_HU.txt / vel_EC.txt / vel_FR.txt
# record class -> (3DEC table name, table file). These are the 100%-intensity
# truncated velocity tables written by integrate_table_accel.py. Each run
# applies its record's base table SCALED by the Table-2 factor, via the
# velocity-z multiplier (e.g. run 1 -> velocity-z 0.50 table 'vel_HU').
BASE_TABLES = {
    "HU12": ("vel_HU", "vel_HU.txt"),
    "EC40": ("vel_EC", "vel_EC.txt"),
    "FR76": ("vel_FR", "vel_FR.txt"),
}
RUN_FROM      = 1          # first run to apply (inclusive)
RUN_TO        = 25         # last run to apply (inclusive)
START_SAVE    = "Part_I_MASON_v7.sav"   # undamaged base state (full sequence)
TAIL_SEC_RECORD = 2.5      # ring-down appended after each record (period ID)
# <<< STRATEGY F (FULL)

OUT_DIR = "stratF_full_results_US{}".format(1 if RECORD_TEST == 9 else 2)
# <<< STRATEGY D

# >>> RETIRED : the Strategy C feedback controls. None is read in Strategy D;
# they are kept so the stratC_results_* folders stay reproducible from history.
#
# THE PROBLEM THIS FIXES
#   Sd_target = scale * Sd(record, T). Until now T was the period identified
#   from the post-run free-vibration ring-down, which is a SMALL-AMPLITUDE
#   measurement. For a wall that has cracked into a rocking mechanism that
#   period stays short (it saturated at ~1.30 * T1 in stratC_results_GI_NORATCH)
#   while the effective response period migrates long. Because all three target
#   spectra are near their minimum at 0.09-0.12 s and rise steeply beyond it,
#   the identified period pinned the whole IDA to the flat part of the spectrum:
#
#     Sd(HU12) = 0.53 mm at 0.092 s ... 4.23 mm at 0.40 s
#     Sd(FR76) = 1.30 mm at 0.092 s ... 28.78 mm at 0.40 s
#
#   GI_NORATCH therefore topped out at Sd_target = 6.77 mm, against ~37 mm in
#   the calibration study and 61 mm of measured table displacement in the run
#   that failed the specimen. The IDA never reached collapse intensity.
#
# THE FIX
#   Look the spectrum up at the SECANT period instead, obtained by equivalent
#   linearisation at the peak response of the preceding run:
#
#       k_sec = V_peak / d_peak                (kN/m, from that run's histories)
#       T_sec = 2*pi*sqrt(M_EFF / k_sec)       (s)
#
#   M_EFF and the initial stiffness are the equivalent-SDOF values from the
#   calibration study: 1.47 t and 6875 kN/m, which reproduce T1 = 0.092 s. This
#   is standard equivalent linearisation, as used in displacement-based
#   assessment, and it is the period the rocking wall actually responds at.
#
#   The ring-down period is still identified and logged as T_end. It remains a
#   damage indicator. It no longer sets the intensity.
USE_SECANT_PERIOD = False  # OFF. The V/d estimator is not sound -- see the
                           # header. True restores it for comparison only.
M_EFF_T           = 1.47   # tonnes, equivalent SDOF mass (calibration study)
K_INIT_KNPM       = 6875.0 # kN/m, equivalent SDOF initial stiffness
T_SEC_MAX         = 0.60   # s. Hard clamp on the lookup period.
SECANT_MONOTONIC  = True   # never let the lookup period shorten. Damage is
                           # irreversible; a smaller excursion in a low-intensity
                           # run legitimately has a stiffer secant, but using it
                           # would drop the demand and make the ladder oscillate.
SD_GROWTH_CAP     = 3.0    # None to disable. Sd_target may not exceed this
                           # multiple of the previous run's, so a bad period
                           # estimate cannot reproduce the v7 runaway. Every
                           # time it binds it is printed and logged.
SCALE_MULT        = 1.0    # RETIRED in Strategy D. Was: multiplier on the
                           # what primes the loop. Calibrated on the observed
                           # stratC_results_SECANT response: the wall cracks
                           # open at Sd_target ~ 7.8 mm, and with the ring-down
                           # lookup nothing below 1.5 ever gets there (1.00 and
                           # 1.25 both stall at T1 for all 25 runs, exactly as
                           # stratC_results_GI_NORATCH did). 1.5 primes at
                           # run 20, 2.0 at run 18 but overshoots from run 19.
                           # Does NOT affect run folder names or the reported
                           # 'scale' column; only Sd_target. See scale_eff.
# <<< SECANT PERIOD

# >>> GI VARIANT : joint property overrides applied after every model restore.
# Set G_I_NEW to None to leave the properties exactly as the base save has them.
G_I_NEW      = None      # None = leave the base save's properties alone.
                         # Deliberately OFF here. G_I was raised to 20 in the
                         # GI variant to suppress the NODAMP_v7 period runaway,
                         # but that runaway is now understood as an artefact of
                         # the intensity bug this driver fixes. Keeping G_I at
                         # the base 13.2 makes this a ONE-VARIABLE change from
                         # stratC_results_NODAMP_v7. Set 20.0 to re-test G_I on
                         # top of the corrected intensity ladder.
G_II_RATIO   = 10.0      # G_II = G_II_RATIO * G_I, as in ANALYSIS_PART_I_MASON.dat
COH_RESIDUAL = 0.0       # Pa. Base model uses 0 -- cracked joints retain nothing.
                         # This is very likely the bigger lever; change it ONLY
                         # after seeing what G_I alone does, so the runs stay
                         # interpretable one variable at a time.
# <<< GI VARIANT
# <<< RATCHETING

# =====================================================================
# 2.  PROTOCOL (Table 2)
# =====================================================================
# Table 2: sequence of applied acceleration records (run, input, scale).
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
    """(duration_s, n_samples) of a 3DEC velocity table file (line 1 = name,
    line 2 = 'N <tab> 0', then 'time <tab> velocity'). Duration = last time."""
    last_t, n = 0.0, 0
    with open(path) as f:
        f.readline(); f.readline()                 # name line, count line
        for line in f:
            q = line.split()
            if len(q) >= 2:
                try:
                    last_t = float(q[0]); n += 1
                except ValueError:
                    pass
    return last_t, n


def build_protocol():
    """Table 2 restricted to [RUN_FROM, RUN_TO]; each run carries its record's
    base-table name/file and solve duration (read from the table on disk)."""
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
    """Import the three 100% base velocity tables into the model (once)."""
    for rec, (name, fname) in BASE_TABLES.items():
        pth = os.path.join(TABLE_DIR, fname)
        if os.path.isfile(pth):
            it.command("table '{}' import '{}'".format(name, cmd_path(pth)))
            print("  imported base table '{}' <- {}".format(name, fname))
        else:
            print("  ! base table missing: {}".format(pth))


PROTOCOL = build_protocol()

# Reference only (the driver reads CSVs by filename, never by index).
# Indices match instrument_history_new.dat + instrument_tilt_v2.dat, and are
# exported by instrument_history_export_v2.dat.
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
    # 3DEC accepts forward slashes on both Windows and Linux. Forcing
    # backslashes corrupts paths on Linux (where '\' is a literal filename
    # character) and was a cause of period-ID failures from CSVs that could
    # not be located. Normalise and emit forward slashes for portability.
    return os.path.normpath(path).replace("\\", "/")

def save_file_path(run_no):
    return os.path.join(OUT_DIR, "stratC_run_{:02d}.sav".format(run_no))

def _find_channel_csv(folder, needle="channel_19"):
    """First CSV in `folder` whose name contains `needle` (case-insensitive).
    Locating the file by scanning the folder makes period identification
    independent of the exact filename or path-separator style the FISH export
    used, which is what previously broke period ID across OS conventions."""
    if not os.path.isdir(folder):
        return None
    hits = [os.path.join(folder, fn) for fn in os.listdir(folder)
            if fn.lower().endswith(".csv") and needle in fn.lower()]
    return sorted(hits)[0] if hits else None

def ch19_csv_path(run_no, record=None, scale=None):
    # Primary: expected per-run folder, located by scan (robust to filename).
    if record is not None and scale is not None:
        scale_str = "{:.2f}".format(scale).replace(".", "p")
        run_folder = os.path.join(OUT_DIR, "Run{:02d}_{}_s{}".format(run_no, record, scale_str))
        hit = _find_channel_csv(run_folder)
        if hit:
            return hit
    # Fallback: any folder starting with Run{nn}_ that holds a channel-19 CSV.
    prefix = "Run{:02d}_".format(run_no)
    if os.path.isdir(OUT_DIR):
        for d in sorted(os.listdir(OUT_DIR)):
            if d.startswith(prefix) and os.path.isdir(os.path.join(OUT_DIR, d)):
                hit = _find_channel_csv(os.path.join(OUT_DIR, d))
                if hit:
                    return hit
    # Last resort: flat path (the caller already handles a missing file).
    return os.path.join(OUT_DIR,
        "Channel_19_DispTopQRight_run{:02d}.csv".format(run_no))

def run_folder_path(run_no, record, scale):
    scale_str = "{:.2f}".format(scale).replace(".", "p")
    return os.path.join(OUT_DIR,
                        "Run{:02d}_{}_s{}".format(run_no, record, scale_str))


def _read_hist(path):
    """3DEC history export: 2 header lines, whitespace separated."""
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
    """Equivalent-linearised period from the run just completed.

        k_sec = V_peak / d_peak,   T_sec = 2*pi*sqrt(M_EFF / k_sec)

    V_peak comes from 'total_shear' (kN, base reaction plus both joist
    reactions) and d_peak from 'rel_disp_top_mm' (mm, top centreline relative
    to the table). Both are exported for every run.

    Returns (T_sec, info) or (None, info) when the histories are unusable, in
    which case the caller falls back to the ring-down period.
    """
    info = {"V_peak_kN": None, "disp_peak_mm": None, "k_sec_kNpm": None}
    folder = run_folder_path(run_no, record, scale)
    fv = _find_channel_csv(folder, "total_shear")
    fd = _find_channel_csv(folder, "rel_disp_top_mm")
    if not (fv and fd):
        print("  ! secant: total_shear / rel_disp_top_mm not found in {}"
              .format(folder))
        return None, info
    av, ad = _read_hist(fv), _read_hist(fd)
    if av is None or ad is None:
        print("  ! secant: history unreadable, falling back to ring-down period")
        return None, info

    v = av[:, 1]
    u = ad[:, 1]
    V_peak = float(np.max(np.abs(v - v[0])))          # kN
    d_peak = float(np.max(np.abs(u - u[0]))) / 1000.0  # mm -> m
    info["V_peak_kN"] = round(V_peak, 4)
    info["disp_peak_mm"] = round(d_peak * 1000.0, 4)
    if d_peak <= 0 or V_peak <= 0:
        print("  ! secant: zero peak response, falling back to ring-down period")
        return None, info

    k_sec = V_peak / d_peak                            # kN/m
    info["k_sec_kNpm"] = round(k_sec, 2)
    T_sec = 2.0 * math.pi * math.sqrt(M_EFF_T / k_sec)
    return T_sec, info


def state_file_path():
    return os.path.join(OUT_DIR, STATE_FILE_NAME)

# =====================================================================
# 4.  CHECKPOINT / RESUME
# =====================================================================
def save_checkpoint(run_no, T_current, summary):
    state = {"last_completed_run": run_no, "T_current": T_current, "summary": summary}
    with open(state_file_path(), "w", newline="\n") as f:
        json.dump(state, f, indent=2)

def load_checkpoint():
    sf = state_file_path()
    if os.path.exists(sf):
        try:
            with open(sf, "r") as f:
                state = json.load(f)
            last = state["last_completed_run"]
            T_cur = state["T_current"]
            summary = state.get("summary", [])
            print("  Checkpoint JSON found: last completed = run {:02d}, "
                  "T_current = {:.4f} s".format(last, T_cur))
            return last + 1, T_cur, summary
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
# 5.  LOAD PRECOMPUTED Sd SPECTRA
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
# 6.  SINUSOID GENERATION  (original symmetric routines, kept for fallback)
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
    """(t, v) for the cosine-acceleration pulse plus its ramped zero tail.

    3DEC is given velocity, and v(0) = 0 exactly because the acceleration is a
    cosine, so the model starts from rest. Factored out of build_velocity_file
    so the PGD calibration integrates the same array that gets written."""
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
    """Displacement the base will trace out, by trapezoidal integration."""
    return np.concatenate(([0.0],
                           np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(t))))


def calibrate_amplitude_pgd(T, PGD_target):
    """Acceleration amplitude giving a table displacement of amplitude
    PGD_target, defined as half the peak-to-peak.

    Analytically A = PGD * w^2 for a pure sinusoid, but the ramped tail shifts
    that slightly, so this calibrates on the actual array. One trial pass is
    enough because the table scales linearly with A."""
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

# >>> RATCHETING : asymmetric pulse generation (record-derived)
def generate_pulse(record, T, Sd_target, run_no, out_dir):
    """
    Returns (vel_path, pulse_dur, v_peak, info_dict).
    Symmetric when USE_RATCHETING is False; otherwise asymmetric with the
    two-regime activation and the record's measured directional bias.
    """
    if not USE_RATCHETING:
        A = calibrate_amplitude(T, Sd_target)
        vp, dur, vpk = build_velocity_file(A, T, run_no, out_dir)
        return vp, dur, vpk, {"A": A, "beta": 1.0, "s": 0}

    rec = RECORD_ASYM[record]
    s = rec["s"]
    beta_target = rec["beta"] ** ASYM_K          # sharpen once
    beta_now = beta_eff(T, T1_init, beta_target) # 1.0 until rocking onset
    A, _ = calibrate_amplitude_asym(T, Sd_target, beta_now, s,
                                    n_cycles, delta_t, xi, k=1.0)
    vp, dur, vpk = build_velocity_asym(A, T, beta_now, s, run_no, out_dir,
                                       n_cycles, delta_t, tail_sec)
    return vp, dur, vpk, {"A": A, "beta": beta_now, "s": s}
# <<< RATCHETING

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
def setup_model_for_dynamic(save_file):
    it.command("model restore '{}'".format(cmd_path(save_file)))
    it.command("model dynamic active on")
    it.command("""
    block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
    block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
    """)
    it.command("call 'instrument_history_new.dat'")
    it.command("""
    block free velocity-z range group 'S'
    """)
    it.command("""
    block free rotation-y range group 'T_B'
    block free rotation-z range group 'T_B'
    block mech damp local 0.0
    block mech damp global 0.0
    """)
    # >>> GI VARIANT : working damping branch. The base driver had the command
    # commented out inside this branch with `pass` as the body, so setting
    # DAMP_RATIO non-zero there silently produced an undamped run while
    # preflight_checks() reported damping as active.
    if DAMP_RATIO > 0:
        cmd = "block mechanical damping rayleigh {ratio} {freq} {dtype}".format(
            ratio=DAMP_RATIO, freq=1.0 / T1_init, dtype=DAMP_TYPE).strip()
        it.command(cmd)
        print("  Rayleigh damping applied: {}".format(cmd))
    else:
        print("  (no Rayleigh damping applied; contact dissipation only)")

    # Joint property overrides. Applied after the restore so the base save can be
    # shared between variants without rebuilding the model.
    if G_I_NEW is not None:
        it.command("block contact property G_I {:g} G_II {:g}".format(
            G_I_NEW, G_II_RATIO * G_I_NEW))
        print("  G_I -> {:g} J/m2 , G_II -> {:g} J/m2 (base model computes 13.2)".format(
            G_I_NEW, G_II_RATIO * G_I_NEW))
    if COH_RESIDUAL:
        it.command("block contact property cohesion-residual {:g}".format(COH_RESIDUAL))
        print("  cohesion-residual -> {:g} Pa".format(COH_RESIDUAL))
    # <<< GI VARIANT
    print("  Model setup complete (BCs + damping + joist contact groups applied).")

# =====================================================================
# 10.  EXECUTE ONE RUN
# =====================================================================
def execute_run(entry):
    """Apply one Table-2 run: the record's 100% base velocity table, scaled by
    the run's factor through the velocity-z multiplier (e.g. velocity-z 0.50
    table 'vel_HU'). The base tables are imported once up front and reused;
    nothing is synthesised, nothing adapts. T_end is logged as a damage measure
    only."""
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

    # T_end is a damage measure only.
    T_end = identify_Tend_from_csv(
        ch19_csv_path(run_no, record, scale), dur)
    print("  T_end = {:.4f} s ({:.2f}x T1_init)   [damage measure]".format(
        T_end, T_end / T1_init))

    # NB: do NOT delete the base table -- it is reused every run.
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")

    return {
        "run": run_no, "record": record, "scale": round(scale, 4),
        "table_file": entry["table_file"],
        "record_dur_s": round(dur, 4),
        "T_end": round(T_end, 6),
        "T_end_over_Tinit": round(T_end / T1_init, 4),
    }

# =====================================================================
# 11.  MAIN DRIVER
# =====================================================================
SPECTRA_FILES = [os.path.join(TABLE_DIR, f) for (_n, f) in BASE_TABLES.values()]
DAT_FILES     = ["instrument_history_new.dat", 
                 "instrument_history_export_v2.dat"]
BASE_SAVE     = START_SAVE

def preflight_checks():
    """Fail loudly, before any solving, if inputs or config are not in order."""
    print("\n" + "=" * 70)
    print("  PRE-FLIGHT CHECKS  (working dir: {})".format(os.getcwd()))
    print("=" * 70)
    # A fresh start needs the base save; a resume restores a stratC_run save.
    fresh = (not os.path.exists(state_file_path()) and
             not any(os.path.exists(save_file_path(r)) for r in range(1, len(PROTOCOL) + 1)))
    required = list(SPECTRA_FILES) + list(DAT_FILES) + ([BASE_SAVE] if fresh else [])
    missing = []
    for f in required:
        ok = os.path.isfile(f)
        if not ok:
            missing.append(f)
        print("   [{:^7}] {}".format("OK" if ok else "MISSING", f))
    dtype = {"": "full Rayleigh (mass+stiffness)", "mass": "mass-proportional",
             "stiffness": "stiffness-proportional"}.get(DAMP_TYPE.strip(), DAMP_TYPE)
    print("   start mode  : {}".format("FRESH" if fresh else "RESUME"))
    if DAMP_RATIO > 0:
        print("   damping     : Rayleigh {:.1f}% at {:.2f} Hz  [{}]".format(
            DAMP_RATIO * 100.0, 1.0 / T1_init, dtype))
    else:
        print("   damping     : NONE (no Rayleigh; contact dissipation only)")
    print("   excitation  : {}   ASYM_K = {}".format(
        "ASYMMETRIC (ratcheting)" if USE_RATCHETING else "symmetric (baseline)", ASYM_K))
    print("   output dir  : {}".format(OUT_DIR))
    if missing:
        raise RuntimeError("Missing required input file(s) in the working "
                           "directory: {}. Place them next to this script or "
                           "fix the working directory.".format(", ".join(missing)))
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
    print("\n  Strategy F: {} full records, {:.0f} s of model time total.".format(
        len(PROTOCOL), _tot))
    print("  This is the fidelity benchmark and it is SLOW. On a machine that "
          "solves")
    print("  at ~230 s wall per s model, that is roughly {:.0f} h. For the "
          "decisive".format(_tot * 230.0 / 3600.0))
    print("  'does the record reach the experimental displacement' test, set "
          "START_SAVE")
    print("  to a damaged state and RUN_FROM/RUN_TO to just the FR76 runs.")
    if resume_from == 1:
        print("\n--- Starting fresh from {} ---".format(START_SAVE))
        setup_model_for_dynamic(START_SAVE)
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
    LOG_COLS = ["run", "record", "scale", "table_file",
                "record_dur_s", "T_end", "T_end_over_Tinit"]
    if log_is_new:
        log_f.write(",".join(LOG_COLS) + "\n")

    for idx in range(resume_from - 1, len(PROTOCOL)):
        entry = PROTOCOL[idx]
        run_summary = execute_run(entry)
        summary.append(run_summary)
        s = run_summary
        log_f.write(",".join(str(s.get(c, "")) for c in LOG_COLS) + "\n")
        log_f.flush()
        # stored only so a resume can report it; it drives nothing
        save_checkpoint(entry["run"], run_summary["T_end"], summary)
        print("  Checkpoint saved after run {:02d}.".format(entry["run"]))
    log_f.close()

    print("\n" + "=" * 70)
    print("Strategy F (full records) complete: {} runs".format(len(PROTOCOL)))
    print("=" * 70)
    print("\nRecords applied (100% base tables, scaled per Table 2):")
    print("  {:>4} {:>6} {:>7} {:>9} {:>9}".format(
        "run", "rec", "scale", "dur_s", "T_end_s"))
    for s in summary:
        print("  {:>4} {:>6} {:>7.2f} {:>9.2f} {:>9.4f}".format(
            s["run"], s["record"], s["scale"], s["record_dur_s"], s["T_end"]))

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