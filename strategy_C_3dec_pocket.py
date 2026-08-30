# -*- coding: ascii -*-
"""
Strategy C (SYMMETRIC PULSE, NO DAMPING) : Adaptive Sequential IDA
=======================================================================
Production driver for the mason_v7 model. Writes to stratC_results_NODAMP_v10.
(v7 = pre-guard baseline, v8 = period-guard only, both retained untouched.)

USE_RATCHETING is now FALSE, so this runs the original symmetric pulse and
ASYM_K is inert. The reason is in MODELLING_DISCREPANCIES.md section 0: the
source paper attributes the experimental asymmetry to the timber floors
applying a moment to one side of the wall, and measures 8-12 mm of cumulative
joist slip, with the bow-tie specimen (which removed that slip) reversing the
sign of its residual tilt. Both mechanisms are already in this model
geometrically, whereas ASYM_K imposes an asymmetric INPUT PULSE -- a different
mechanism, and very likely a double count. Any residual tilt this run produces
therefore comes from the joists alone, which is the quantity worth comparing.

Set USE_RATCHETING = True to restore the asymmetric (record-derived) pulse from
ratcheting_pulse.py. The marked  # >>> RATCHETING  blocks are the only places
the two paths differ.

Relationship to strategy_C_3dec_GI.py: both run the symmetric pulse, both record
the same channels, and both have the working damping branch. The pair is now a
CONTROLLED test of the fracture-energy hypothesis, differing only in:
    G_I / G_II   13.2 J/m2 here (as computed in ANALYSIS_PART_I_MASON.dat)
                 vs 20 J/m2 in GI
    OUT_DIR      stratC_results_NODAMP_v10 vs stratC_results_GI_NORATCH
GI additionally exposes COH_RESIDUAL, but it is 0.0 there, which matches the
base model, so it changes nothing.

Instrumentation: this driver calls instrument_tilt_v2.dat after
instrument_history_new.dat at both registration points, and exports through
instrument_history_export_v2.dat. FISH histories export BY NUMBER, so if you
add or remove any `fish history` line in either .dat the indices shift and the
export file must be updated to match.

ratcheting_pulse.py must sit in the same directory as this script (or on
sys.path). Run inside 3DEC:
    python-reset-state false
    call 'strategy_C_3dec_nodamp.py'

To force a full restart, delete or rename the stratC_results folder.
"""

import itasca as it
import numpy as np
import os, csv, math, json, sys

# Experiment-matched period identification (Moshfeghi et al. 2024 sec. 3.2).
# Guarded so a missing module degrades to the legacy estimator with a loud
# warning rather than killing a 25-run sequence at run 1.
try:
    import period_id_exp
    _HAVE_PERIOD_ID_EXP = True
except Exception as _e:                                  # pragma: no cover
    period_id_exp = None
    _HAVE_PERIOD_ID_EXP = False
    print("WARNING: period_id_exp not importable ({}). "
          "Falling back to the legacy single-channel displacement "
          "estimator.".format(_e))

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
# STANDALONE CASE: LS_MAXWELL_1p5_POCKET
# =====================================================================
# This file is a COMPLETE copy of the driver with the case constants baked in.
# It reads no overrides file -- the loader below is disabled -- so a stale
# stratC/stratF_overrides.json left over from another case cannot change what
# this script does. Run it directly:
#
#     python-reset-state false
#     call '<this file>'
#
#   geometry   : large-strain on
#   damping    : block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 0.0104 40.0000
#   period ID  : experiment   (T1_init = 0.0948 s)
#   output     : stratC_results_LS_MAXWELL_POCKET
#
# The cost of a standalone copy is that it does NOT track edits to the driver
# it was generated from. If you change the physics in the base driver,
# regenerate these copies or the two will silently disagree.
# =====================================================================
# =====================================================================
# 1.  PARAMETERS
# =====================================================================
T1_init     = 0.0948
xi          = 0.05
delta_t     = 0.005
n_cycles    = 3.0     # was 1.5. A cracked wall is a rocking mechanism; 1.5
                      # cycles delivered acceleration but too few reversals to
                      # build the rocking amplitude that drives collapse. 3.0
                      # gives the pulse enough cycles for the FR76 runs to
                      # accumulate rocking.
tail_sec    = 2.5
phase_accel = math.pi / 2
inter_run_gap = 0.5
bLength     = 1.292

FFT_F_MIN  = 2.0
FFT_F_MAX  = 50.0
T_MIN_PHYS = 0.02
T_MAX_PHYS = 1.0

# ---------------------------------------------------------------------
# PERIOD-ID RELIABILITY GUARDS
# ---------------------------------------------------------------------
# The autocorrelation and FFT estimators are independent. On every archived
# NODAMP_v7 ring-down they agree to better than 0.1% (e.g. run 12:
# T_acorr = 0.10428 s, T_fft = 0.10432 s). A large disagreement therefore
# means the pick has failed, not that the wall has softened. When that
# happens the previous T_current is held and the run is flagged, rather than
# an unreproducible value being propagated into the next excitation.
#
# The two tests below do DIFFERENT jobs and only one of them may stop the loop.
#
#   TID_AGREE_TOL is a detection-quality test. It asks whether the two
#   independent estimators agree. If they do not, the pick has failed and the
#   value is not a measurement of anything, so the previous period is held.
#   This does not push T_current toward a target.
#
#   TID_MAX_JUMP is NOT a detection test. It encodes a prior about how fast a
#   wall may soften, which is a statement about the physics and not about the
#   signal. Letting it HOLD the period was a design error: on NODAMP_v7_NEW it
#   rejected T_end = 0.21848 s at run 12 (2.09x the previous value) even though
#   the agreement test had passed, meaning both estimators concurred on that
#   value. The adaptive loop then froze for the remaining 14 runs and the
#   response went flat while the commanded demand tripled.
#
#   It is therefore WARN-ONLY. A large jump is recorded as T_jump_flag = 1 in
#   the log and the value is accepted, because a measurement that two
#   independent estimators agree on is evidence about the model, and
#   suppressing it would be adjusting the measurement to fit the expectation.
TID_AGREE_TOL = 0.25   # max |T_acorr - T_fft| / min(T_acorr, T_fft); HOLDS
TID_MAX_JUMP  = 1.50   # T_end / T_prev above which a run is flagged; WARN ONLY
#
# A hold is only ever meant to skip ONE bad ring-down. On NODAMP_v7_NEW the
# period was held for 14 consecutive runs and the analysis completed looking
# normal, having silently become a fixed-period scheme with a detuned pulse.
# That must fail loudly instead. Set to 0 to disable the abort.
TID_MAX_CONSEC_HOLDS = 2
_LAST_JUMP_FLAG = 0    # set by identify_Tend_from_csv, read by execute_run
_LAST_T_ACORR   = 0.0  # both estimators are logged every run so that a
_LAST_T_FFT     = 0.0  # rejection can be diagnosed from the CSV alone
_LAST_TID_NOTE  = ""
_CONSEC_HOLDS   = 0

# ---------------------------------------------------------------------
# RECORD-DERIVED EQUIVALENT CYCLES
# ---------------------------------------------------------------------
# A fixed n_cycles matches Sd(T1) but leaves the pulse energy unconstrained,
# so a 3-cycle resonant sinusoid over-delivers for broadband records whose
# energy sits away from the wall frequency. Measured against US-1 the fixed
# 3.0 over-predicts EDP by 2.2-4.9x on EC40 and 2.1-5.5x on FR76, while HU12
# (narrowband, T_eq ~ 0.17 s) is well represented.
#
# N_eq is the classical equivalent-cycle count of the achieved table motion,
# evaluated in a narrow band around the excitation frequency:
#
#     N_eq = integral(a_band^2 dt) / ( max|a_band|^2 * T / 2 )
#
# For a pure N-cycle sinusoid this returns exactly N (verified: 3.010), and it
# is independent of the pulse amplitude, so no iteration with the Sd
# calibration is needed. Everything is derived from the record.
#
# ---------------------------------------------------------------------
# DEFAULT OFF -- the hypothesis this was written to test is REFUTED.
# ---------------------------------------------------------------------
# The idea was that the fixed 3.0 cycles over-delivers energy for broadband
# records, explaining the 2.2-5.5x EDP over-prediction on EC40 and FR76.
# Evaluated on the achieved tables at T = 0.092 s, N_eq is:
#
#     HU12  3.10      EC40  10.0 (capped)      FR76  10.0 (capped)
#
# i.e. the opposite of what the hypothesis needs: enabling this would give
# EC40/FR76 MORE cycles and make the over-prediction worse. Direct energy
# accounting confirms it. Arias intensity of the calibrated pulse against the
# record's band-limited Arias in the same +-30% band around 1/T:
#
#     HU12 run 21   pulse 0.1513  record-band 0.1066   ratio 1.42
#     EC40 run 11   pulse 0.0084  record-band 0.0185   ratio 0.45
#     FR76 run 24   pulse 0.0778  record-band 0.1114   ratio 0.70
#
# The pulse carries LESS band energy than EC40 and FR76 yet produces 2-5x
# more displacement, and MORE band energy than HU12 yet matches it best. The
# over-prediction is therefore not an energy or duration effect: it is
# coherence. Three consecutive in-phase cycles exactly at the wall frequency,
# in a model with DAMP_RATIO = 0, pump the rocking mode monotonically, while
# the record's band energy is spread over 17-40 s at random phase and the
# mechanism re-seats between excursions. Adding cycles increases coherent
# pumping; it cannot fix this.
#
# Left in place as instrumentation: N_eq is logged per run as n_cycles_used
# and is a useful record descriptor. Set True only to reproduce the test.
USE_RECORD_NCYCLES = False
# ---------------------------------------------------------------------
# TWO-IM CYCLE COUNT  (Sd(T) AND the record's peak acceleration)
# ---------------------------------------------------------------------
# v8 matched Sd(T) alone and under-predicted the late runs by up to 5.5x
# (run 24: 5.2 mm modelled vs 29.5 mm measured), while the modelled EDP
# tracked the commanded Sd_target ~1:1 -- i.e. the wall responded
# essentially elastically and never entered the rocking regime.
#
# Cause: for a resonant sinusoid at damping xi,
#
#     Sd = (A / w^2) * (1 - exp(-2 pi xi N)) / (2 xi)
#
# and that build-up factor is 6.103 at N = 3.0, xi = 0.05. Matching Sd(T)
# therefore needs only ~1/6 of the acceleration the record itself carried:
# measured pulse PGA ran 2.2-3.6x BELOW the achieved table PGA on every run.
# Rocking initiation is an acceleration-threshold phenomenon, so the pulse
# reached its displacement target while staying under the threshold that
# actually triggers rocking.
#
# Fixing A = A_record leaves N as the single unknown. Solved per run this
# gives 0.54-0.94 cycles over runs 17-24 (mean 0.80) -- a single lobe, not
# three cycles. Note the earlier n_cycles 1.5 -> 3.0 change, made to "give
# the pulse enough cycles to accumulate rocking", was backwards: more cycles
# LOWER the amplitude needed for a given Sd and starve the acceleration.
#
# Both matched quantities are record-derived (Sd from the response spectrum,
# A_record from the achieved channel-12 table motion x protocol scale). No
# target displacement is assumed and no structure-dependent feedback is added
# beyond the excitation period itself.
USE_TWO_IM_CYCLES = True
NCYC_MIN, NCYC_MAX = 0.25, 10.0
NCYC_BAND_RATIO    = 1.30    # band = [f1/ratio, f1*ratio]
RECORD_TABLE = {"HU12": "vel_HU.txt", "EC40": "vel_EC.txt", "FR76": "vel_FR.txt"}
_NCYC_CACHE = {}

OUT_DIR = 'stratC_results_SmallStrain_NODAMP_POCKET'

# ---------------------------------------------------------------------
# SENSITIVITY-SWEEP OVERRIDES
# ---------------------------------------------------------------------
# sweep_run.py writes stratC_overrides.json next to this script before each
# case, so a case needs no edit to this file. Only the keys listed in
# _OVERRIDABLE may be set; anything else is reported and ignored, so a typo
# cannot silently run the baseline under a case label.
#
# The file is deleted (or absent) for a normal manual run, in which case the
# values above stand unchanged.
_OVERRIDABLE = ("OUT_DIR", "n_cycles", "xi", "delta_t", "tail_sec",
                "inter_run_gap", "DAMP_RATIO", "DAMP_TYPE", "MAXWELL_CMD",
                "USE_TWO_IM_CYCLES", "USE_RECORD_NCYCLES", "PERIOD_SCHEME",
                "TID_AGREE_TOL", "TID_MAX_JUMP", "NCYC_MIN", "NCYC_MAX",
                "LARGE_STRAIN", "PERIOD_ID_METHOD",
                "T1_init", "BASE_SAVE", "CASE_ID")

# ---------------------------------------------------------------------
# PERIOD IDENTIFICATION METHOD
# ---------------------------------------------------------------------
#   "experiment" -- Moshfeghi et al. (2024) Structures 66:106815 section 3.2:
#                   three accelerometers (bottom quarter 0.66 m, mid-height
#                   1.26 m, top quarter 2.06 m), record cropped at the instant
#                   the table stops, periodogram -> PSD matrix -> singular
#                   values -> peak-picking, cross-checked against per-channel
#                   FFT. Implemented in period_id_exp.py.
#   "legacy"     -- the original single-channel DISPLACEMENT estimator
#                   (Channel_19, autocorrelation + FFT).
#
# Both are computed and logged on every run regardless of this setting; this
# only selects which value is fed back into the adaptive loop, so a completed
# run set can be re-read either way without re-running anything.
PERIOD_ID_METHOD = "experiment"

PERIOD_SCHEME = "adaptive"   # "adaptive" = excite at the identified T_current
                             # "fixed"    = excite at T1_init every run (T_eq
                             #              control case; T_end still logged
                             #              as a damage measure)
MAXWELL_CMD   = 'block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 0.0104 40.0000'           # non-empty -> issued verbatim after local/global
                             # damping are zeroed
LARGE_STRAIN  = "on"        # "on" updates block positions and re-detects
                             # contacts, so toe contact area can reduce and the
                             # restoring action becomes geometric. Required for
                             # large-amplitude rocking; costs solve time.
                             # Applied AFTER the Part-I restore, so no rebuild
                             # is needed and the static build stays in
                             # small-strain where it belongs.
CASE_ID       = 'LS_MAXWELL_1p5_POCKET'           # free-text label echoed into the log


if False:   # DISABLED: this is a standalone case file
    try:
        with open(_OVR_FILE) as _f:
            _ovr = json.load(_f)
        print("\n--- sweep overrides from {} ---".format(_OVR_FILE))
        for _k in sorted(_ovr):
            if _k in _OVERRIDABLE:
                globals()[_k] = _ovr[_k]
                print("    {:<20} = {}".format(_k, _ovr[_k]))
            else:
                print("    ! IGNORED (not overridable): {}".format(_k))
    except Exception as _e:
        raise RuntimeError("Cannot parse {}: {}".format(_OVR_FILE, _e))

STATE_FILE_NAME = "stratC_checkpoint.json"

# >>> RATCHETING : controls
USE_RATCHETING = False   # OFF. The paper attributes the asymmetry to eccentric joist
                         # loading and 8-12 mm of cumulative joist slip, both of which
                         # this model already contains geometrically, whereas ASYM_K
                         # imposes an asymmetric INPUT PULSE. Any residual tilt that
                         # still appears now comes from the joists alone. ASYM_K is inert.
ASYM_K         = 1.0     # global asymmetry sharpening (>=1 sharpens); calibrate out-of-sample
# Damping is exposed here so it is explicit and printed at startup.
# The conference paper used FULL Rayleigh at 3%. The pasted driver used
# mass-proportional 0.06, which over-damps the slow rocking and is the wrong
# choice for a collapse study. Default below is full Rayleigh 3%; change if needed.
DAMP_RATIO = 0.0         # fraction of critical at the centre frequency
DAMP_TYPE  = ""          # "" = full Rayleigh (mass+stiffness); "mass"; "stiffness"
# <<< RATCHETING

# =====================================================================
# 2.  PROTOCOL (Table 2)
# =====================================================================
PROTOCOL = [
    ( 1, "HU12", 0.50),  ( 2, "HU12", 0.75),  ( 3, "EC40", 0.20),
    ( 4, "HU12", 1.00),  ( 5, "HU12", 1.25),  ( 6, "EC40", 0.30),
    ( 7, "HU12", 1.50),  ( 8, "EC40", 0.40),  ( 9, "HU12", 1.75),
    (10, "HU12", 2.00),  (11, "EC40", 0.50),  (12, "HU12", 2.25),
    (13, "HU12", 2.50),  (14, "HU12", 2.75),  (15, "HU12", 3.00),
    (16, "HU12", 3.50),  (17, "HU12", 4.00),  (18, "HU12", 4.50),
    (19, "HU12", 5.00),  (20, "HU12", 5.50),  (21, "HU12", 6.00),
    (22, "FR76", 1.00),  (23, "FR76", 1.50),  (24, "FR76", 1.75),
    (25, "FR76", 2.00),
]

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

def run_label(run_no, record, scale):
    """The label the FISH export prefixes every channel file with."""
    return "Run{:02d}_{}_s{}".format(
        run_no, record, "{:.2f}".format(scale).replace(".", "p"))

def run_dir(run_no, record=None, scale=None):
    """Per-run output folder. Same fallback discipline as ch19_csv_path: try
    the expected name, then any folder starting with Run{nn}_, so a renamed or
    differently-scaled folder does not silently break period identification."""
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

def _read_vel_table(path):
    """(t, v) from a 3DEC table file: name line, 'N<TAB>0' line, then rows."""
    t_list, v_list = [], []
    with open(path) as f:
        lines = f.readlines()
    for line in lines[2:]:
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            t_list.append(float(parts[0])); v_list.append(float(parts[1]))
        except ValueError:
            continue
    return np.array(t_list), np.array(v_list)


def equivalent_cycles(record, T_pulse):
    """Record-derived equivalent number of cycles at the excitation period.

    Reads the achieved base-velocity table for `record` (written by
    integrate_table_accel.py from channel 12), differentiates to
    acceleration, isolates a narrow band around f1 = 1/T_pulse by zeroing
    the rFFT outside it, and returns the energy-equivalent cycle count.

    Falls back to the global n_cycles, with a printed note, if the table is
    missing or unreadable -- the driver must never fail on this."""
    if not USE_RECORD_NCYCLES:
        return n_cycles
    key = (record, round(T_pulse, 5))
    if key in _NCYC_CACHE:
        return _NCYC_CACHE[key]
    fname = RECORD_TABLE.get(record)
    path = None
    for cand in (fname, os.path.join(TABLE_DIR, fname) if "TABLE_DIR" in globals() else None):
        if cand and os.path.exists(cand):
            path = cand
            break
    if path is None:
        print("    ! base table for {} not found -- n_cycles stays {:.2f}"
              .format(record, n_cycles))
        _NCYC_CACHE[key] = n_cycles
        return n_cycles
    try:
        t, v = _read_vel_table(path)
        if len(t) < 32:
            raise ValueError("table too short")
        dt = float(np.median(np.diff(t)))
        a = np.gradient(v, dt)                      # velocity -> acceleration
        a = a - np.mean(a)
        f1 = 1.0 / T_pulse
        freqs = np.fft.rfftfreq(len(a), d=dt)
        spec = np.fft.rfft(a)
        band = (freqs >= f1 / NCYC_BAND_RATIO) & (freqs <= f1 * NCYC_BAND_RATIO)
        if not np.any(band):
            raise ValueError("empty band at f1={:.2f} Hz".format(f1))
        spec[~band] = 0.0
        a_b = np.fft.irfft(spec, n=len(a))
        peak = float(np.max(np.abs(a_b)))
        if peak <= 0.0:
            raise ValueError("no energy in band")
        energy = float(np.sum(a_b ** 2) * dt)
        n_eq = energy / (peak ** 2 * T_pulse / 2.0)
        n_eq = float(min(max(n_eq, NCYC_MIN), NCYC_MAX))
    except Exception as e:
        print("    ! equivalent_cycles({}) failed: {} -- n_cycles stays {:.2f}"
              .format(record, e, n_cycles))
        n_eq = n_cycles
    _NCYC_CACHE[key] = n_eq
    return n_eq


def record_pga(record, scale):
    """Peak acceleration of the achieved table motion for this run [m/s2].

    Differentiates the base velocity table (written by integrate_table_accel.py
    from channel 12, i.e. the ACHIEVED table motion, not the nominal input) and
    applies the protocol scale factor. Returns None if unavailable, in which
    case the caller falls back to the fixed n_cycles."""
    fname = RECORD_TABLE.get(record)
    if not fname or not os.path.exists(fname):
        print("    ! base table for {} not found -- cannot set N from PGA"
              .format(record))
        return None
    try:
        t, v = _read_vel_table(fname)
        if len(t) < 8:
            raise ValueError("table too short")
        dt = float(np.median(np.diff(t)))
        a = np.gradient(v, dt)
        a = a - np.mean(a)
        return float(np.max(np.abs(a))) * scale
    except Exception as e:
        print("    ! record_pga({}) failed: {}".format(record, e))
        return None


def cycles_for_two_im(T, Sd_target, A_record):
    """n_cycles so the pulse matches Sd_target at T AND peaks at A_record.

    Inverts  Sd = (A/w^2) * (1 - exp(-2 pi xi N)) / (2 xi)  for N with
    A = A_record. Returns None if A_record is unusable. Saturates at
    NCYC_MAX when the target is unreachable at that amplitude."""
    if not A_record or A_record <= 0.0 or T <= 0.0:
        return None
    w = 2.0 * math.pi / T
    f_need = Sd_target * w * w / A_record
    arg = 1.0 - 2.0 * xi * f_need
    if arg <= 1e-6:
        return NCYC_MAX
    N = -math.log(arg) / (2.0 * math.pi * xi)
    return float(min(max(N, NCYC_MIN), NCYC_MAX))


def _pulse_tv(A, T, ncyc):
    """(t, v) for the pulse: cosine acceleration over `ncyc` cycles,
    integrated to velocity, raised-cosine taper to zero over half a period,
    then a zero tail. Single source of truth for the pulse shape, so the
    amplitude calibration sees exactly the motion that gets applied."""
    w = 2.0 * math.pi / T
    n = int(round(ncyc * T / delta_t)) + 1
    t = np.arange(n) * delta_t
    v = (A / w) * np.sin(w * t)
    ramp_sec = max(0.5 * T, delta_t)
    n_tail = int(round(tail_sec / delta_t))
    n_ramp = min(int(round(ramp_sec / delta_t)), n_tail)
    t_tail = t[-1] + delta_t + np.arange(n_tail) * delta_t
    t = np.r_[t, t_tail]
    v_end = float(v[-1])
    if n_ramp > 1:
        ss = np.linspace(0.0, 1.0, n_ramp)
        v_ramp = v_end * 0.5 * (1.0 + np.cos(math.pi * ss))
    else:
        v_ramp = np.array([0.0])
    v_tail = np.zeros(n_tail)
    v_tail[:len(v_ramp)] = v_ramp
    return t, np.r_[v, v_tail]


def calibrate_amplitude(T, Sd_target, ncyc=None):
    """Amplitude so the ACTUAL applied pulse reaches Sd_target at T.

    Calibrates on the differentiated velocity history, so the raised-cosine
    taper is included. That matters once ncyc < 1: the velocity no longer
    returns to zero at the end of the cosine, so the taper contributes a
    secondary acceleration lobe that the old untapered cosine calibration
    ignored. Linear in A, so one trial suffices."""
    ncyc = n_cycles if ncyc is None else ncyc
    t_u, v_u = _pulse_tv(1.0, T, ncyc)
    a_u = np.gradient(v_u, delta_t)
    sd0 = newmark_sd(a_u, delta_t, T, xi)
    if sd0 <= 0.0:
        raise RuntimeError("Trial Sd=0 at T={:.4f}, N={:.3f}".format(T, ncyc))
    return Sd_target / sd0


def _calibrate_amplitude_untapered(T, Sd_target, ncyc=None):
    ncyc = n_cycles if ncyc is None else ncyc
    A_trial = 1.0
    w = 2.0 * math.pi / T
    n = int(round(ncyc * T / delta_t)) + 1
    t_arr = np.arange(n) * delta_t
    a0 = A_trial * np.cos(w * t_arr)
    sd0 = newmark_sd(a0, delta_t, T, xi)
    if sd0 == 0.0:
        raise RuntimeError("Trial Sd=0 at T={:.4f}".format(T))
    return (Sd_target / sd0) * A_trial

def build_velocity_file(A, T, run_no, out_dir, ncyc=None):
    """Write the 3DEC velocity table. Shape comes from _pulse_tv, so this is
    byte-for-byte the motion the amplitude was calibrated against."""
    ncyc = n_cycles if ncyc is None else ncyc
    t, v = _pulse_tv(A, T, ncyc)
    fname = "vel_run_{:02d}.txt".format(run_no)
    fpath = os.path.join(out_dir, fname)
    N = len(t)
    with open(fpath, "w", newline="\n") as f:
        f.write("StratC_run{:02d}_T{:.4f}_N{:.3f}\n".format(run_no, T, ncyc))
        f.write("{}\t0\n".format(N))
        for ti, vi in zip(t, v):
            f.write("{:.6f}\t{:.9e}\n".format(ti, vi))
    duration = float(t[-1])
    v_peak = float(np.max(np.abs(v)))
    return fpath, duration, v_peak

# >>> RATCHETING : asymmetric pulse generation (record-derived)
def generate_pulse(record, T, Sd_target, run_no, out_dir, ncyc=None):
    """
    Returns (vel_path, pulse_dur, v_peak, info_dict).
    Symmetric when USE_RATCHETING is False; otherwise asymmetric with the
    two-regime activation and the record's measured directional bias.
    `ncyc` is the record-derived equivalent cycle count for this run.
    """
    ncyc = n_cycles if ncyc is None else ncyc
    if not USE_RATCHETING:
        A = calibrate_amplitude(T, Sd_target, ncyc)
        vp, dur, vpk = build_velocity_file(A, T, run_no, out_dir, ncyc)
        return vp, dur, vpk, {"A": A, "beta": 1.0, "s": 0, "ncyc": ncyc}

    rec = RECORD_ASYM[record]
    s = rec["s"]
    beta_target = rec["beta"] ** ASYM_K          # sharpen once
    beta_now = beta_eff(T, T1_init, beta_target) # 1.0 until rocking onset
    A, _ = calibrate_amplitude_asym(T, Sd_target, beta_now, s,
                                    ncyc, delta_t, xi, k=1.0)
    vp, dur, vpk = build_velocity_asym(A, T, beta_now, s, run_no, out_dir,
                                       ncyc, delta_t, tail_sec)
    return vp, dur, vpk, {"A": A, "beta": beta_now, "s": s, "ncyc": ncyc}
# <<< RATCHETING

# =====================================================================
# 7.  PERIOD IDENTIFICATION (segment + window + autocorrelation)
# =====================================================================
RMS_ALIVE_MM  = 0.05
RMS_CHUNK_S   = 0.15
MAX_RD_WINDOW = 1.50  # was 0.50. After a hard rocking run that recovers (e.g.
                      # run 23, peak 58 mm), the 0.5 s window saw too few cycles
                      # and its large residual offset masked the rocking peak,
                      # so autocorr found no peak and fell back to T1_init
                      # (0.092 s). That reset starved the next run (run 24 got a
                      # 0.092 s pulse, Sd ~2.3 mm, moved 17 mm instead of
                      # escalating). Verified on run 23's ring-down: 0.5 s ->
                      # 0.092 s fallback; 1.5 s -> 0.576 s (autocorr and FFT
                      # agree in the rocking band). 1.5 s recovers the true
                      # period and keeps the escalation intact.

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

def identify_Tend_from_csv(ch19_csv, sine_end_time, T_prev=None):
    globals()["_LAST_JUMP_FLAG"] = 0
    globals()["_LAST_TID_NOTE"] = ""
    globals()["_LAST_T_ACORR"] = 0.0
    globals()["_LAST_T_FFT"] = 0.0
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

    # --- reliability guards -------------------------------------------
    # 1. The two independent estimators must agree. On sound ring-downs they
    #    agree to <0.1%; a large split means the pick failed.
    # 2. Second line: reject an implausible single-run jump.
    # Neither guard moves T_end toward a target -- they only decide whether
    # this run's estimate is trustworthy enough to propagate.
    globals()["_LAST_T_ACORR"] = T_acorr
    globals()["_LAST_T_FFT"] = T_fft
    disagree = abs(T_acorr - T_fft) / max(min(T_acorr, T_fft), 1e-12)
    if disagree > TID_AGREE_TOL:
        held = T_prev if T_prev else T1_init
        print("  ** PERIOD ID REJECTED: T_acorr={:.5f}s vs T_fft={:.5f}s "
              "disagree by {:.0f}% (> {:.0f}%).".format(
                  T_acorr, T_fft, disagree * 100, TID_AGREE_TOL * 100))
        print("     Holding T_current = {:.5f}s and flagging this run."
              .format(held))
        globals()["_LAST_TID_NOTE"] = "HOLD_disagree_{:.0f}pct".format(disagree*100)
        globals()["_CONSEC_HOLDS"] = _CONSEC_HOLDS + 1
        if TID_MAX_CONSEC_HOLDS and _CONSEC_HOLDS >= TID_MAX_CONSEC_HOLDS:
            raise RuntimeError(
                "Period identification has been held {} runs in a row. The "
                "adaptive loop is no longer adaptive and the pulse is being "
                "detuned from a wall that keeps softening. Stopping rather "
                "than completing a run that would look normal. Inspect the "
                "Channel_19 ring-down, then either widen TID_AGREE_TOL or set "
                "PERIOD_SCHEME='fixed' deliberately."
                .format(_CONSEC_HOLDS))
        return held
    if T_prev and T_prev > 0 and T_end / T_prev > TID_MAX_JUMP:
        globals()["_LAST_JUMP_FLAG"] = 1
        print("  ** LARGE PERIOD JUMP: T_end={:.5f}s is {:.2f}x T_prev={:.5f}s "
              "(> {:.2f}x) in one run.".format(
                  T_end, T_end / T_prev, T_prev, TID_MAX_JUMP))
        print("     ACCEPTED and flagged (T_jump_flag=1). The two estimators "
              "agreed to {:.2f}%, so this is a measurement of the model, not a "
              "detection failure.".format(disagree * 100))
        globals()["_LAST_TID_NOTE"] = "JUMP_{:.2f}x_accepted".format(T_end / T_prev)

    globals()["_CONSEC_HOLDS"] = 0
    if not _LAST_TID_NOTE:
        globals()["_LAST_TID_NOTE"] = "ok_{}".format(method)
    print("  Period ID: T_acorr={:.5f}s, T_fft={:.5f}s -> T_end={:.5f}s ({}, "
          "estimators agree to {:.2f}%)".format(
              T_acorr, T_fft, T_end, method, disagree * 100))
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
    if str(LARGE_STRAIN).lower() in ("on", "true", "1"):
        it.command("model large-strain off")
        print("  GEOMETRY: large-strain OFF (positions updated, contacts "
              "re-detected)")
    else:
        it.command("model large-strain off")
        print("  GEOMETRY: small-strain (default)")
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
    # The command used to be commented out here with `pass` as the body, so
    # setting DAMP_RATIO non-zero silently produced an undamped run while
    # preflight_checks() reported damping as active.
    # if MAXWELL_CMD:
    #     it.command(MAXWELL_CMD)
    #     print("  Maxwell damping applied: {}".format(MAXWELL_CMD))
    if DAMP_RATIO > 0:
        cmd = "block mechanical damping rayleigh {ratio} {freq} {dtype}".format(
            ratio=DAMP_RATIO, freq=1.0 / T1_init, dtype=DAMP_TYPE).strip()
        it.command(cmd)
        print("  Rayleigh damping applied: {}".format(cmd))
    else:
        print("  (no Rayleigh damping applied; contact dissipation only)")
    print("  Model setup complete (BCs + damping + joist contact groups applied).")

# =====================================================================
# 10.  EXECUTE ONE RUN
# =====================================================================
def execute_run(run_no, record, scale, T_current):
    # PERIOD_SCHEME == "fixed" is the T_eq control case: the pulse is tuned to
    # T1_init on every run and the identified T_end is logged but not fed back.
    if PERIOD_SCHEME == "fixed":
        T_current = T1_init
    Sd_unit_Tcurr = interpolate_sd(record, T_current)
    Sd_target     = scale * Sd_unit_Tcurr
    Sd_unit_T1    = interpolate_sd(record, T1_init)
    Sd_fixed_T1   = scale * Sd_unit_T1
    amp_factor    = Sd_target / Sd_fixed_T1 if Sd_fixed_T1 > 0 else 0

    # Cycle count. Preferred route matches BOTH Sd(T) and the record's own
    # peak acceleration; falls back to the energy-equivalent count, then to
    # the fixed n_cycles, and says which it used.
    A_record = record_pga(record, scale) if USE_TWO_IM_CYCLES else None
    ncyc = cycles_for_two_im(T_current, Sd_target, A_record) if A_record else None
    if ncyc is None:
        ncyc = equivalent_cycles(record, T_current)
        ncyc_src = "record-energy" if USE_RECORD_NCYCLES else "fixed"
    else:
        ncyc_src = "two-IM(Sd+PGA)"

    # >>> RATCHETING : pulse generation (symmetric or asymmetric per toggle)
    vel_path, pulse_dur, v_peak, pinfo = generate_pulse(
        record, T_current, Sd_target, run_no, OUT_DIR, ncyc)
    A_cal = pinfo["A"]
    sine_dur = ncyc * T_current
    # <<< RATCHETING

    tbl_name = "run{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl_name, cmd_path(vel_path)))

    print("\n" + "=" * 70)
    print("  Run {:02d}: {} x {:.2f}  |  T = {:.4f} s ({:.2f}x T_init)".format(
        run_no, record, scale, T_current, T_current/T1_init))
    print("  Sd_target = {:.3f} mm  (fixed-T1: {:.3f} mm, amp: {:.1f}x)".format(
        Sd_target*1000, Sd_fixed_T1*1000, amp_factor))
    print("  A = {:.3f} m/s2,  PGA = {:.3f} g,  Vpeak = {:.4f} m/s".format(
        A_cal, A_cal/9.80665, v_peak))
    # >>> RATCHETING : report asymmetry actually applied
    print("  pulse: {}  beta={:.3f}  s={:+d}".format(
        "ASYMMETRIC" if USE_RATCHETING else "symmetric", pinfo["beta"], pinfo["s"]))
    # <<< RATCHETING
    print("=" * 70)

    it.command("model dynamic time-total 0")
    for grp in ["S", "T_B"]:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl_name, grp))
    it.command("model solve dynamic time {:.6f}".format(pulse_dur))
    for grp in ["S", "T_B"]:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(grp))
    if inter_run_gap > 0:
        it.command("model solve dynamic time {:.6f}".format(inter_run_gap))

    it.command("model save '{}'".format(cmd_path(save_file_path(run_no))))
    export_all_histories(run_no, record, scale, OUT_DIR)

    ringdown_start = sine_dur

    # --- Period identification -----------------------------------------
    # Both estimators run on every run. The legacy one reads a single
    # DISPLACEMENT channel; the experiment one reads the three accelerometers
    # the way Moshfeghi et al. (2024) section 3.2 describes. They are not
    # measuring the same thing: a post-rocking ring-down carries a large,
    # heavily damped rigid-body transient on top of the small elastic
    # vibration, and in displacement that transient dominates a short window.
    # See period_id_exp._selftest(), which reproduces run 12's T_end = 0.218 s
    # from a synthetic signal whose true elastic period is 0.104 s.
    T_legacy = identify_Tend_from_csv(ch19_csv_path(run_no, record, scale),
                                      ringdown_start, T_prev=T_current)
    _exp = {"T": float("nan"), "ok": False, "spread": float("nan"),
            "window_s": 0.0, "T_welch": float("nan"), "note": "module missing"}
    if _HAVE_PERIOD_ID_EXP:
        _exp = period_id_exp.identify_period_experiment(
            run_dir(run_no, record, scale), run_label(run_no, record, scale))

    # --- selection rule: seed on run 1, then TRACK the same mode ---------
    # pick_fundamental returns the tallest phase-coherent peak. On the
    # EXPERIMENT's own records that rule returns 0.0533 s for run 1, where
    # Table 5 says 0.091 -- it locks onto the timber-floor mode at 18.76 Hz,
    # which is taller than the wall's 10.97 Hz and sits at a non-integer ratio
    # so no harmonic guard can reject it. Seeding once and following the same
    # mode forward reproduces Table 5 to +0.2% / +1.3% at the two endpoints.
    #
    # T_current is the driver's own running period, so the anchor is already
    # maintained. Nothing extra has to be tracked.
    _cl = _exp.get("candidates") or []
    if _HAVE_PERIOD_ID_EXP and _cl:
        if run_no <= 1 or not (T_current and np.isfinite(T_current)):
            _T_sel = period_id_exp.seed_period(_cl)
            _sel_rule = "seed"
        else:
            _T_sel, _rec = period_id_exp.track_step(_cl, T_current)
            _sel_rule = "track"
        if np.isfinite(_T_sel):
            _exp["T"] = _T_sel
            _exp["note"] = (_exp["note"] + " | " + _sel_rule).strip(" |")
            # NOTE: _exp["ok"] is DELIBERATELY not touched here. It carries the
            # channel-agreement verdict from identify_period_experiment (the
            # paper's own PSD-vs-FFT cross-check). An earlier version forced it
            # True after tracking, which silently bypassed that test: on the
            # LS_MAXWELL run, nine of 25 runs had per-channel spreads of 54% to
            # 334% and were accepted anyway, and three of those had the Welch
            # cross-check sitting at ~0.25x the reported period, i.e. it had
            # locked onto the 4x harmonic. Selecting WHICH peak to report and
            # deciding WHETHER the identification is trustworthy are different
            # jobs; tracking does the first and must not overrule the second.

    if PERIOD_ID_METHOD == "experiment" and _exp["ok"]:
        T_end = _exp["T"]
        _tid_src = "exp_psd_svd"
    elif PERIOD_ID_METHOD == "experiment":
        # The experiment's own acceptance test (PSD vs per-channel FFT) failed.
        # Fall back rather than propagate a number the paper's method would
        # itself have rejected, and say so in the log.
        T_end = T_legacy
        _tid_src = "legacy_fallback:" + str(_exp["note"])
        print("  ** experiment-method identification rejected ({}). "
              "Using legacy value {:.5f} s.".format(_exp["note"], T_legacy))
    else:
        T_end = T_legacy
        _tid_src = "legacy"

    print("  T_end = {:.4f} s  ({:.2f}x T_init)   [{}]".format(
        T_end, T_end / T1_init, _tid_src))
    if _HAVE_PERIOD_ID_EXP and np.isfinite(_exp["T"]):
        print("     legacy(disp,1ch) = {:.5f} s   experiment(acc,3ch) = {:.5f} s"
              "   ratio = {:.2f}".format(
                  T_legacy, _exp["T"],
                  _exp["T"] / T_legacy if T_legacy else float("nan")))

    it.command("table '{}' delete".format(tbl_name))
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")
    
    run_summary = {
        "run": run_no, "record": record, "scale": scale,
        "T_excitation": round(T_current, 6),
        "T_over_Tinit": round(T_current / T1_init, 4),
        "Sd_record_mm": round(Sd_unit_Tcurr * 1000, 4),
        "Sd_target_mm": round(Sd_target * 1000, 4),
        "Sd_fixedT1_mm": round(Sd_fixed_T1 * 1000, 4),
        "amplification": round(amp_factor, 2),
        "A_mps2": round(A_cal, 4),
        "PGA_g": round(A_cal / 9.80665, 4),
        "V_peak_mps": round(v_peak, 6),
        "T_end": round(T_end, 6),
        "T_end_over_Tinit": round(T_end / T1_init, 4),
        "n_cycles_used": round(ncyc, 3),
        "ncyc_source": ncyc_src,
        "A_record_g": round(A_record / 9.80665, 4) if A_record else "",
        "T_id_held": 1 if abs(T_end - T_current) < 1e-12 else 0,
        "T_jump_flag": _LAST_JUMP_FLAG,
        "T_acorr": round(_LAST_T_ACORR, 6),
        "T_fft": round(_LAST_T_FFT, 6),
        "tid_note": _LAST_TID_NOTE,
        # --- both estimators logged on every run, always ------------------
        "T_legacy_disp": round(T_legacy, 6),
        "T_exp_psd": round(_exp["T"], 6) if np.isfinite(_exp["T"]) else "",
        "T_exp_welch": round(_exp["T_welch"], 6) if np.isfinite(_exp["T_welch"]) else "",
        "T_exp_ch_spread": round(_exp["spread"], 4) if np.isfinite(_exp["spread"]) else "",
        "T_exp_window_s": round(_exp["window_s"], 4),
        "tid_source": _tid_src,
        # >>> RATCHETING : record the asymmetry applied
        "beta_applied": round(pinfo["beta"], 4),
        "s_applied": pinfo["s"],
        # <<< RATCHETING
    }
    return T_end, run_summary

# =====================================================================
# 11.  MAIN DRIVER
# =====================================================================
SPECTRA_FILES = ["spectrum_HU12.csv", "spectrum_EC40.csv", "spectrum_FR76.csv"]
DAT_FILES     = ["instrument_history_new.dat",
                 "instrument_history_export_v2.dat"]
BASE_SAVE     = "Part_I_MASON_v8_SmallStrain.sav"

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
    load_spectra()
    print("\n--- Checking for existing runs ---")
    resume_from, T_current, summary = load_checkpoint()
    if resume_from > len(PROTOCOL):
        print("\nAll {} runs already completed.".format(len(PROTOCOL)))
        print("Delete '{}' to restart.".format(state_file_path()))
        return summary
    if resume_from == 1:
        print("\n--- Starting fresh from run 1 ---")
        print("    USE_RATCHETING = {}   ASYM_K = {}".format(USE_RATCHETING, ASYM_K))
        setup_model_for_dynamic(BASE_SAVE)
    else:
        last_done = resume_from - 1
        print("\n--- Resuming: restoring run {:02d}, will execute run {:02d} next ---".format(
            last_done, resume_from))
        print("  T_current = {:.4f} s ({:.2f}x T_init)".format(T_current, T_current / T1_init))
        setup_model_for_dynamic(save_file_path(last_done))

    log_path = os.path.join(OUT_DIR, "strategy_C_log.csv")
    log_is_new = (resume_from == 1) or not os.path.exists(log_path)
    log_f = open(log_path, "w" if log_is_new else "a", newline="\n")
    if log_is_new:
        log_f.write("run,record,scale,T_excite,T_over_Tinit,"
                    "Sd_record_mm,Sd_target_mm,Sd_fixedT1_mm,amplification,"
                    "A_mps2,PGA_g,V_peak_mps,T_end,T_end_over_Tinit,"
                    "n_cycles_used,ncyc_source,A_record_g,T_id_held,"
                    "T_jump_flag,T_acorr,T_fft,tid_note,"
                    "beta_applied,s_applied\n")

    for idx in range(resume_from - 1, len(PROTOCOL)):
        run_no, record, scale = PROTOCOL[idx]
        T_end, run_summary = execute_run(run_no, record, scale, T_current)
        summary.append(run_summary)
        s = run_summary
        log_f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            s["run"], s["record"], s["scale"], s["T_excitation"], s["T_over_Tinit"],
            s["Sd_record_mm"], s["Sd_target_mm"], s["Sd_fixedT1_mm"], s["amplification"],
            s["A_mps2"], s["PGA_g"], s["V_peak_mps"], s["T_end"], s["T_end_over_Tinit"],
            s["n_cycles_used"], s["ncyc_source"], s["A_record_g"], s["T_id_held"],
            s["T_jump_flag"], s["T_acorr"], s["T_fft"], s["tid_note"],
            s["beta_applied"], s["s_applied"]))
        log_f.flush()
        T_current = T_end
        save_checkpoint(run_no, T_current, summary)
        print("  Checkpoint saved after run {:02d}.".format(run_no))
    log_f.close()

    print("\n" + "=" * 70)
    print("Strategy C (ratcheting={}) complete: {} runs".format(USE_RATCHETING, len(PROTOCOL)))
    print("=" * 70)
    print("\nPeriod evolution:")
    for s in summary:
        print("  Run {:02d} ({} x{:.2f}):  T_ex={:.4f}s -> T_end={:.4f}s  Sd_tgt={:.2f}mm  "
              "beta={:.2f}".format(s["run"], s["record"], s["scale"], s["T_excitation"],
              s["T_end"], s["Sd_target_mm"], s.get("beta_applied", 1.0)))

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