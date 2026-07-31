# -*- coding: ascii -*-
"""
Strategy C (RATCHETING, NO DAMPING) : Adaptive Sequential IDA
=======================================================================
This is the ratcheting-enabled version of strategy_C_3dec_python.py.
It is identical to the original except for the marked  # >>> RATCHETING
blocks, which swap the symmetric pulse generation for the asymmetric
(record-derived) pulse implemented in ratcheting_pulse.py.

Toggle behaviour with USE_RATCHETING:
    USE_RATCHETING = False  -> exactly the original symmetric method
    USE_RATCHETING = True   -> asymmetric ratcheting pulse (Modification A)

ratcheting_pulse.py must sit in the same directory as this script (or on
sys.path). Run inside 3DEC:
    python-reset-state false
    call 'strategy_C_3dec_nodamp.py'

To force a full restart, delete or rename the stratC_results folder.
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

OUT_DIR = "stratC_results_NODAMP_v6_NEW"
STATE_FILE_NAME = "stratC_checkpoint.json"

# >>> RATCHETING : controls
USE_RATCHETING = True    # False reproduces the original symmetric method exactly
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
# Indices match instrument_history_new.dat / instrument_history_export_new.dat.
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

def state_file_path():
    return os.path.join(OUT_DIR, STATE_FILE_NAME)

# =====================================================================
# 4.  CHECKPOINT / RESUME
# =====================================================================
def save_checkpoint(run_no, T_current, summary):
    state = {"last_completed_run": run_no, "T_current": T_current, "summary": summary}
    with open(state_file_path(), "w") as f:
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
    with open(fpath, "w") as f:
        f.write("StratC_run{:02d}_T{:.4f}\n".format(run_no, T))
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
    it.command("call 'instrument_history_export_new.dat'")

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
    if DAMP_RATIO > 0:
        # NO-DAMP variant: Rayleigh damping intentionally disabled (commented out).
        # it.command("block mechanical damping rayleigh {ratio} {freq} {dtype}".format(
        #     ratio=DAMP_RATIO, freq=1.0/T1_init, dtype=DAMP_TYPE).strip())
        pass
    else:
        print("  (no Rayleigh damping applied; contact dissipation only)")
    print("  Model setup complete (BCs + damping + joist contact groups applied).")

# =====================================================================
# 10.  EXECUTE ONE RUN
# =====================================================================
def execute_run(run_no, record, scale, T_current):
    Sd_unit_Tcurr = interpolate_sd(record, T_current)
    Sd_target     = scale * Sd_unit_Tcurr
    Sd_unit_T1    = interpolate_sd(record, T1_init)
    Sd_fixed_T1   = scale * Sd_unit_T1
    amp_factor    = Sd_target / Sd_fixed_T1 if Sd_fixed_T1 > 0 else 0

    # >>> RATCHETING : pulse generation (symmetric or asymmetric per toggle)
    vel_path, pulse_dur, v_peak, pinfo = generate_pulse(
        record, T_current, Sd_target, run_no, OUT_DIR)
    A_cal = pinfo["A"]
    sine_dur = n_cycles * T_current
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
    T_end = identify_Tend_from_csv(ch19_csv_path(run_no, record, scale), ringdown_start)
    print("  T_end = {:.4f} s  ({:.2f}x T_init)".format(T_end, T_end/T1_init))

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
DAT_FILES     = ["instrument_history_new.dat", "instrument_history_export_new.dat"]
BASE_SAVE     = "Part_I_MASON_v6.sav"

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
        setup_model_for_dynamic("Part_I_MASON_v6.sav")
    else:
        last_done = resume_from - 1
        print("\n--- Resuming: restoring run {:02d}, will execute run {:02d} next ---".format(
            last_done, resume_from))
        print("  T_current = {:.4f} s ({:.2f}x T_init)".format(T_current, T_current / T1_init))
        setup_model_for_dynamic(save_file_path(last_done))

    log_path = os.path.join(OUT_DIR, "strategy_C_log.csv")
    log_is_new = (resume_from == 1) or not os.path.exists(log_path)
    log_f = open(log_path, "w" if log_is_new else "a")
    if log_is_new:
        log_f.write("run,record,scale,T_excite,T_over_Tinit,"
                    "Sd_record_mm,Sd_target_mm,Sd_fixedT1_mm,amplification,"
                    "A_mps2,PGA_g,V_peak_mps,T_end,T_end_over_Tinit,"
                    "beta_applied,s_applied\n")

    for idx in range(resume_from - 1, len(PROTOCOL)):
        run_no, record, scale = PROTOCOL[idx]
        T_end, run_summary = execute_run(run_no, record, scale, T_current)
        summary.append(run_summary)
        s = run_summary
        log_f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            s["run"], s["record"], s["scale"], s["T_excitation"], s["T_over_Tinit"],
            s["Sd_record_mm"], s["Sd_target_mm"], s["Sd_fixedT1_mm"], s["amplification"],
            s["A_mps2"], s["PGA_g"], s["V_peak_mps"], s["T_end"], s["T_end_over_Tinit"],
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
