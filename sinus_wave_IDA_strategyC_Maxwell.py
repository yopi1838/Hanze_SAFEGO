# -*- coding: ascii -*-
"""
Strategy C: Adaptive Sequential IDA with Period Tracking
==========================================================
3DEC-internal Python script with CHECKPOINT / RESUME support.

On launch the script scans for existing save files and exported CSVs.
If prior runs are found it restores the last completed state, recovers
T_end from the Ch19 ring-down CSV, and continues from the next run.

After each run:
  1. Export all instrument histories
  2. FFT on Ch19 ring-down -> T_end
  3. Save state + checkpoint JSON
  4. Generate next velocity table at T_end
  5. Apply and continue (cumulative damage)

Run inside 3DEC:
    python-reset-state false
    call 'strategy_C_3dec_python.py'

To force a full restart, delete or rename the stratC_results folder.
"""

import itasca as it
import numpy as np
import os, csv, math, json

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

OUT_DIR = "stratC_results_MP_NEW"
STATE_FILE_NAME = "stratC_checkpoint.json"

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

INSTRUMENT_CHANNELS = [
    "Channel_1_DispBot", "Channel_2_DispMid",
    "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight",
    "Channel_5_DispTable",
    "Channel_8_LeftJoistBelowSlab1", "Channel_9_LeftJoistAboveSlab1",
    "Channel_10_LeftJoistBelowSlab2", "Channel_11_LeftJoistBelowSlab2",
    "Channel_12_AccTable",
    "Channel_13_AccSlab1", "Channel_14_AccSlab2",
    "Channel_15_AccBot", "Channel_16_AccMid", "Channel_17_AccTop",
    "Channel_19_DispTopQRight", "Channel_20_TopBeam",
]

# FISH histories are stored by numeric index, not label.
# Map: index -> descriptive filename for export
# Indices 1-7: defined in ANALYSIS_PART_I_MASON.dat
# Indices 8-13: defined in instrument_history.dat (tilt angles)
FISH_HISTORIES = {
    1: "Record_Disp",
    2: "Bot_Quarter_Disp",
    3: "Mid_Disp",
    4: "Top_Quarter_A_Disp",
    5: "Top_Quarter_B_Disp",
    6: "cstav",
    7: "ncstav",
    8: "tilt_bot_seg",
    9: "tilt_low_seg",
    10: "tilt_up_seg",
    11: "tilt_beam_seg",
    12: "tilt_full_wall",
    13: "rel_disp_top_mm",
}

# =====================================================================
# 3.  PATH HELPERS
# =====================================================================

def cmd_path(path):
    """Convert path for use in it.command() -- 3DEC on Windows needs backslashes."""
    return path.replace("/", "\\")

def save_file_path(run_no):
    return os.path.join(OUT_DIR, "stratC_run_{:02d}.sav".format(run_no))

def ch19_csv_path(run_no, record=None, scale=None):
    """
    Path to Ch19 CSV. If record/scale given, use per-run subfolder.
    Falls back to scanning for the folder if record/scale unknown (resume).
    """
    if record is not None and scale is not None:
        scale_str = "{:.2f}".format(scale).replace(".", "p")
        run_folder = "Run{:02d}_{}_s{}".format(run_no, record, scale_str)
        return os.path.join(OUT_DIR, run_folder,
            "{}_Channel_19_DispTopQRight.csv".format(run_folder))
    else:
        # Scan: find any folder starting with Run{run_no:02d}_
        prefix = "Run{:02d}_".format(run_no)
        if os.path.exists(OUT_DIR):
            for d in os.listdir(OUT_DIR):
                if d.startswith(prefix) and os.path.isdir(os.path.join(OUT_DIR, d)):
                    candidate = os.path.join(OUT_DIR, d,
                        "{}_Channel_19_DispTopQRight.csv".format(d))
                    if os.path.exists(candidate):
                        return candidate
        # Fallback flat path
        return os.path.join(OUT_DIR,
            "Channel_19_DispTopQRight_run{:02d}.csv".format(run_no))

def state_file_path():
    return os.path.join(OUT_DIR, STATE_FILE_NAME)

# =====================================================================
# 4.  CHECKPOINT / RESUME
# =====================================================================

def save_checkpoint(run_no, T_current, summary):
    state = {
        "last_completed_run": run_no,
        "T_current": T_current,
        "summary": summary,
    }
    with open(state_file_path(), "w", newline="\n") as f:
        json.dump(state, f, indent=2)


def load_checkpoint():
    """
    Determine where to resume.
    Returns (resume_from_run, T_current, summary).
    resume_from_run = next run to execute (1 = fresh start).
    """
    sf = state_file_path()

    # --- Try JSON checkpoint first (most reliable) ---
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
            print("  WARNING: checkpoint corrupt ({}), "
                  "scanning save files...".format(e))

    # --- Fallback: scan for save files from highest to lowest ---
    last_found = 0
    for run_no in range(len(PROTOCOL), 0, -1):
        if os.path.exists(save_file_path(run_no)):
            last_found = run_no
            break

    if last_found == 0:
        print("  No prior runs found. Starting fresh.")
        return 1, T1_init, []

    # --- Recover T_end from Ch19 CSV of last completed run ---
    ch19_path = ch19_csv_path(last_found)
    if os.path.exists(ch19_path):
        sine_dur_approx = n_cycles * T1_init  # best guess
        T_recovered = identify_Tend_from_csv(ch19_path, sine_dur_approx)
        print("  Recovered: last run = {:02d}, "
              "T_end = {:.4f} s (from Ch19 FFT)".format(
                  last_found, T_recovered))
        return last_found + 1, T_recovered, []

    # --- Save exists but no CSV -> check run before ---
    if last_found > 1:
        prev_ch19 = ch19_csv_path(last_found - 1)
        if os.path.exists(prev_ch19):
            sine_dur_approx = n_cycles * T1_init
            T_recovered = identify_Tend_from_csv(prev_ch19, sine_dur_approx)
            print("  Run {:02d} CSV missing; recovered T from run {:02d}: "
                  "{:.4f} s".format(last_found, last_found-1, T_recovered))
            # Re-run the incomplete one
            return last_found, T_recovered, []

    print("  WARNING: cannot recover T_current. Using T1_init.")
    return last_found, T1_init, []


def check_run_complete(run_no):
    """True if both save file and Ch19 CSV exist for this run."""
    return (os.path.exists(save_file_path(run_no)) and
            os.path.exists(ch19_csv_path(run_no)))

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
# 6.  SINUSOID GENERATION
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


def build_velocity_file(A, T, run_no, out_dir):
    w = 2.0 * math.pi / T
    n = int(round(n_cycles * T / delta_t)) + 1
    t = np.arange(n) * delta_t
    V0 = A / w
    v = V0 * np.sin(w * t)

    ramp_sec = max(0.5 * T, delta_t)
    n_tail = int(round(tail_sec / delta_t))
    n_ramp = min(int(round(ramp_sec / delta_t)), n_tail)

    t_tail = t[-1] + delta_t + np.arange(n_tail) * delta_t
    t = np.r_[t, t_tail]

    v_end = float(v[-1])
    if n_ramp > 1:
        s = np.linspace(0.0, 1.0, n_ramp)
        wcos = 0.5 * (1.0 + np.cos(math.pi * s))
        v_ramp = v_end * wcos
    else:
        v_ramp = np.array([0.0])

    v_tail = np.zeros(n_tail)
    v_tail[:len(v_ramp)] = v_ramp
    v = np.r_[v, v_tail]

    fname = "vel_run_{:02d}.txt".format(run_no)
    fpath = os.path.join(out_dir, fname)
    N = len(t)
    with open(fpath, "w", newline="\n") as f:
        f.write("StratC_run{:02d}_T{:.4f}\n".format(run_no, T))
        f.write("{}\t0\n".format(N))
        for ti, vi in zip(t, v):
            f.write("{:.6f}\t{:.9e}\n".format(ti, vi))

    duration = float(t[-1])
    v_peak = float(np.max(np.abs(v)))
    return fpath, duration, v_peak

# =====================================================================
# 7.  PERIOD IDENTIFICATION (v2: segment + window + autocorrelation)
# =====================================================================

# Adaptive window parameters
RMS_ALIVE_MM  = 0.05    # minimum RMS (mm) to consider signal alive
RMS_CHUNK_S   = 0.15    # window chunk size for RMS check (seconds)
MAX_RD_WINDOW = 0.50    # max ring-down window to use (seconds)


def extract_last_segment(t_all, u_all):
    """
    3DEC 'history export' dumps ALL accumulated history.
    Since 'model dynamic time-total 0' resets the clock each run,
    the CSV contains multiple segments separated by backward time jumps.
    Extract only the LAST segment.
    """
    dt = np.diff(t_all)
    neg_idx = np.where(dt < -0.5)[0]  # jump back by > 0.5 s
    n_segments = len(neg_idx) + 1

    if len(neg_idx) > 0:
        last_start = neg_idx[-1] + 1
    else:
        last_start = 0

    return t_all[last_start:], u_all[last_start:], n_segments


def find_live_window(t, u, sine_end_time):
    """
    Find the portion of ring-down where oscillation is still visible.
    Scans in chunks and stops when RMS drops below threshold.
    """
    mask = t >= sine_end_time
    t_rd = t[mask]
    u_rd = u[mask]

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
    max_idx = max(max_idx, int(0.2 / dt_hist))  # at least 0.2 s

    t_live = t_rd[:max_idx]
    u_live = u_rd[:max_idx]
    window_dur = float(t_live[-1] - t_live[0]) if len(t_live) > 1 else 0

    return t_live, u_live, window_dur


def identify_Tend_from_csv(ch19_csv, sine_end_time):
    """
    Full pipeline:
      1. Read CSV, extract LAST segment (fixes concatenated history)
      2. Find LIVE ring-down window (ignores dead tail)
      3. 16x zero-padded FFT + parabolic interpolation
      4. Autocorrelation (primary method for decaying masonry signals)
    """
    # --- Read CSV ---
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

    t_all = data[:, 0]
    u_all = data[:, 1]

    # --- Extract last segment ---
    t_seg, u_seg, n_segments = extract_last_segment(t_all, u_all)
    if n_segments > 1:
        print("    CSV has {} segments (previous runs), using last ({} pts)".format(
            n_segments, len(t_seg)))

    # --- Find live ring-down window ---
    t_live, u_live, win_dur = find_live_window(t_seg, u_seg, sine_end_time)

    if len(u_live) < 16:
        print("  WARNING: live ring-down too short ({} pts)".format(len(u_live)))
        return T1_init

    dt_hist = float(np.median(np.diff(t_live)))
    if dt_hist <= 0:
        dt_hist = delta_t

    # --- De-mean + Hanning ---
    u_dm = u_live - np.mean(u_live)
    u_win = u_dm * np.hanning(len(u_dm))

    N_orig = len(u_win)

    # ----- FFT with 16x zero-padding -----
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

    # ----- Autocorrelation (primary for decaying masonry signals) -----
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
                y0 = float(acorr[best_lag - 1])
                y1 = float(acorr[best_lag])
                y2 = float(acorr[best_lag + 1])
                d2 = y0 - 2*y1 + y2
                if abs(d2) > 1e-12:
                    lag_ref = best_lag + 0.5*(y0 - y2)/d2
                else:
                    lag_ref = float(best_lag)
            else:
                lag_ref = float(best_lag)
            T_acorr = lag_ref * dt_hist

    # ----- Combine: autocorrelation primary, FFT fallback -----
    if T_acorr > T_MIN_PHYS and T_acorr < T_MAX_PHYS:
        T_end = T_acorr
        method = "autocorr"
    else:
        T_end = T_fft
        method = "FFT"

    T_end = max(T_MIN_PHYS, min(T_MAX_PHYS, T_end))

    print("  Period ID: T_acorr={:.5f}s, T_fft={:.5f}s -> T_end={:.5f}s ({})".format(
        T_acorr, T_fft, T_end, method))
    print("    (live window={:.3f}s, N_live={}, segments={})".format(
        win_dur, len(u_live), n_segments))

    return T_end

# =====================================================================
# 8.  EXPORT HISTORIES
# =====================================================================

def export_all_histories(run_no, record, scale, out_dir):
    # Set FISH variables for the export .dat file
    # exportdir = full path to subfolder
    # runlabel  = just the subfolder name (used as filename prefix)
    scale_str = "{:.2f}".format(scale).replace(".", "p")
    run_label = "Run{:02d}_{}_s{}".format(run_no, record, scale_str)
    full_folder = os.path.join(out_dir, run_label)
    os.makedirs(full_folder, exist_ok=True)

    it.command("[exportdir='{}']".format(full_folder))
    it.command("[runlabel='{}']".format(run_label))
    it.command("call 'instrument_history_export_new.dat'")

# =====================================================================
# 9.  MODEL SETUP (fresh or resumed)
# =====================================================================

def setup_model_for_dynamic(save_file):
    """
    Restore a save file and apply dynamic-phase BCs.
    Called both for fresh start (Part_I_MASON.sav) and resume
    (stratC_run_XX.sav). Applied BCs are not in the save file
    so they must be re-applied every time.
    """
    it.command("model restore '{}'".format(cmd_path(save_file)))
    it.command("model dynamic active on")

    # --- Group joist contacts for reaction measurement ---
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
    it.command("""
    block mechanical damping rayleigh 0.06 {freq} mass
    """.format(freq=1.0/T1_init))
    print("  Model setup complete (BCs + damping + joist contact groups applied).")

# =====================================================================
# 10.  EXECUTE ONE RUN
# =====================================================================

def execute_run(run_no, record, scale, T_current):
    """
    Full single-run cycle: generate pulse, apply, solve, export, FFT.
    Returns (T_end, summary_dict).
    """
    Sd_unit_Tcurr = interpolate_sd(record, T_current)
    Sd_target     = scale * Sd_unit_Tcurr

    Sd_unit_T1    = interpolate_sd(record, T1_init)
    Sd_fixed_T1   = scale * Sd_unit_T1
    amp_factor    = Sd_target / Sd_fixed_T1 if Sd_fixed_T1 > 0 else 0

    A_cal = calibrate_amplitude(T_current, Sd_target)
    vel_path, pulse_dur, v_peak = build_velocity_file(
        A_cal, T_current, run_no, OUT_DIR)
    sine_dur = n_cycles * T_current

    tbl_name = "run{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl_name, cmd_path(vel_path)))

    print("\n" + "=" * 70)
    print("  Run {:02d}: {} x {:.2f}  |  T = {:.4f} s ({:.2f}x T_init)".format(
        run_no, record, scale, T_current, T_current/T1_init))
    print("  Sd_target = {:.3f} mm  (fixed-T1: {:.3f} mm, amp: {:.1f}x)".format(
        Sd_target*1000, Sd_fixed_T1*1000, amp_factor))
    print("  A = {:.3f} m/s2,  PGA = {:.3f} g,  Vpeak = {:.4f} m/s".format(
        A_cal, A_cal/9.80665, v_peak))
    print("=" * 70)

    # Reset dynamic time so each run starts from t=0
    it.command("model dynamic time-total 0")

    # Apply velocity-z to base (table) and top beam only
    # Joists (Joist_S1, Joist_S2) remain fixed -- they are not on the table
    for grp in ["S", "T_B"]:
        it.command(
            "block apply velocity-z 1.0 "
            "table '{}' "
            "range group '{}'".format(tbl_name, grp))

    it.command("model solve dynamic time {:.6f}".format(pulse_dur))

    for grp in ["S", "T_B"]:
        it.command(
            "block gridpoint apply-remove velocity-z "
            "range group '{}'".format(grp))

    if inter_run_gap > 0:
        it.command("model solve dynamic time {:.6f}".format(inter_run_gap))

    it.command("model save '{}'".format(cmd_path(save_file_path(run_no))))
    export_all_histories(run_no, record, scale, OUT_DIR)

    # Ring-down starts at sine_dur (since time was reset to 0)
    ringdown_start = sine_dur
    T_end = identify_Tend_from_csv(ch19_csv_path(run_no, record, scale), ringdown_start)

    print("  T_end = {:.4f} s  ({:.2f}x T_init)".format(
        T_end, T_end/T1_init))

    it.command("table '{}' delete".format(tbl_name))

    # Delete all histories so next run's CSV is clean (not concatenated)
    # Then re-register instrument histories for the next run
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")

    run_summary = {
        "run":              run_no,
        "record":           record,
        "scale":            scale,
        "T_excitation":     round(T_current, 6),
        "T_over_Tinit":     round(T_current / T1_init, 4),
        "Sd_record_mm":     round(Sd_unit_Tcurr * 1000, 4),
        "Sd_target_mm":     round(Sd_target * 1000, 4),
        "Sd_fixedT1_mm":    round(Sd_fixed_T1 * 1000, 4),
        "amplification":    round(amp_factor, 2),
        "A_mps2":           round(A_cal, 4),
        "PGA_g":            round(A_cal / 9.80665, 4),
        "V_peak_mps":       round(v_peak, 6),
        "T_end":            round(T_end, 6),
        "T_end_over_Tinit": round(T_end / T1_init, 4),
    }

    return T_end, run_summary

# =====================================================================
# 11.  MAIN DRIVER
# =====================================================================

def run_strategy_C():

    os.makedirs(OUT_DIR, exist_ok=True)
    load_spectra()

    # --- Determine resume point ---
    print("\n--- Checking for existing runs ---")
    resume_from, T_current, summary = load_checkpoint()

    if resume_from > len(PROTOCOL):
        print("\nAll {} runs already completed.".format(len(PROTOCOL)))
        print("Delete '{}' to restart.".format(state_file_path()))
        return summary

    if resume_from == 1:
        print("\n--- Starting fresh from run 1 ---")
        setup_model_for_dynamic("Part_I_MASON.sav")
    else:
        last_done = resume_from - 1
        print("\n--- Resuming: restoring run {:02d}, "
              "will execute run {:02d} next ---".format(last_done, resume_from))
        print("  T_current = {:.4f} s ({:.2f}x T_init)".format(
            T_current, T_current / T1_init))
        setup_model_for_dynamic(save_file_path(last_done))

    # --- Open / append log ---
    log_path = os.path.join(OUT_DIR, "strategy_C_log.csv")
    log_is_new = (resume_from == 1) or not os.path.exists(log_path)
    log_f = open(log_path, "w" if log_is_new else "a", newline="\n")
    if log_is_new:
        log_f.write("run,record,scale,"
                    "T_excite,T_over_Tinit,"
                    "Sd_record_mm,Sd_target_mm,Sd_fixedT1_mm,amplification,"
                    "A_mps2,PGA_g,V_peak_mps,"
                    "T_end,T_end_over_Tinit\n")

    # --- Run loop ---
    for idx in range(resume_from - 1, len(PROTOCOL)):
        run_no, record, scale = PROTOCOL[idx]

        T_end, run_summary = execute_run(run_no, record, scale, T_current)

        summary.append(run_summary)

        s = run_summary
        log_f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            s["run"], s["record"], s["scale"],
            s["T_excitation"], s["T_over_Tinit"],
            s["Sd_record_mm"], s["Sd_target_mm"],
            s["Sd_fixedT1_mm"], s["amplification"],
            s["A_mps2"], s["PGA_g"], s["V_peak_mps"],
            s["T_end"], s["T_end_over_Tinit"]))
        log_f.flush()

        T_current = T_end
        save_checkpoint(run_no, T_current, summary)

        print("  Checkpoint saved after run {:02d}.".format(run_no))

    log_f.close()

    # --- Final report ---
    print("\n" + "=" * 70)
    print("Strategy C complete: {} runs".format(len(PROTOCOL)))
    print("=" * 70)
    print("\nPeriod evolution:")
    for s in summary:
        print("  Run {:02d} ({} x{:.2f}):  T_ex={:.4f}s -> T_end={:.4f}s  "
              "Sd_tgt={:.2f}mm  amp={:.1f}x".format(
                  s["run"], s["record"], s["scale"],
                  s["T_excitation"], s["T_end"],
                  s["Sd_target_mm"], s["amplification"]))

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
