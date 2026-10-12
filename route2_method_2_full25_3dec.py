# -*- coding: utf-8 -*-
"""Run fixed-target, record-informed Method 2 over the US-1 Runs 1-25.

The initial undamaged small-strain model is restored once, then damage is
carried from run to run. Wave A is carried at the fixed wall period T1 and
targets scale*Sd_record(T1). Wave B targets scale*Sd_record(Teff), while its
carrier is the record's predeclared 5%-damped pseudo-Sa peak period. Neither
target period nor carrier is updated from wall response: this is
record-informed, not damage-adaptive.

Wave B keeps the record's predominant period. Before the sequence, the script
selects its smallest feasible cycle count and relative polarity from a fixed
candidate order, then fits both amplitudes to the target spectrum. This choice
uses only the input record and the predeclared targets; it never uses the prior
run's wall response. A failed fit stops before model state changes.

Each completed run is checkpointed. Re-execution resumes from its saved state.
Use a new OUTPUT_ROOT for a clean sequence; keep this separate from the four-
method and experimental result folders.
"""

import csv
import json
import math
import os

import numpy as np
import itasca as it


# ================================ CONFIG ================================
BASE_SAVE = "part_I_mason_SS_LATEST.sav"
OUTPUT_ROOT = "route2_method2_full25_v11"

T1_JSON = os.path.join("bilinear_idealisation", "bilinear_periods.json")
SPECTRA = {"HU12": "spectrum_HU12.csv", "EC40": "spectrum_EC40.csv", "FR76": "spectrum_FR76.csv"}

HISTORY_SETUP = "instrument_history_new.dat"
BASESHEAR_SETUP = "instrument_baseshear_exp.dat"
HISTORY_EXPORT = "instrument_history_export_v2.dat"

RECORD_TP = {"HU12": 0.120, "EC40": 0.255, "FR76": 0.260}  # 5%-damped pseudo-Sa peak from record_properties.csv
RUNS = [
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
# Each standalone copy has exactly one method tuple here:
# (output_name, target_mode, carrier_mode)
METHOD = ("2_fixed_record", "fixed", "record")

XI = 0.05
DELTA_T = 0.005
TAIL_SEC = 2.5
PGA_WARN_G = 0.90
FIT_TOL = 0.02        # relative spectral target tolerance for the joint solver
MAX_WAVE_B_CYCLES = 4
WAVE_B_PHASE_ORDER = (1.0, -1.0)  # prefer same polarity; invert only if needed
FIT_RATIO_GRID_POINTS = 81
FIT_BISECT_ITER = 32
FIT_ABORT_ON_FAILURE = True
FIT_PREFLIGHT_NAME = "method2_fit_preflight.csv"
CHECKPOINT_NAME = "route2_method2_checkpoint.json"

DRIVE_GROUPS = ["S", "T_B"]


# ================================ INPUTS =================================
def cmd_path(path):
    """Use portable forward slashes in Itasca command strings."""
    return os.path.normpath(path).replace("\\", "/")


def require_files(paths):
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        lines = ["Required inputs are missing:"] + ["  " + p for p in missing]
        if BASE_SAVE in missing:
            lines.append(
                "The initial small-strain save is not in this workspace. "
                "Copy the intended virgin save here or set BASE_SAVE to its full path."
            )
        raise RuntimeError("\n".join(lines))


def load_periods(path):
    with open(path, "r") as f:
        data = json.load(f)
    model = data["model"]
    return float(model["T1_s"]), float(model["Teff_s"])


def load_spectra():
    result = {}
    for record, path in SPECTRA.items():
        data = np.genfromtxt(path, delimiter=",", skip_header=1)
        if data.ndim != 2 or data.shape[1] < 2 or len(data) < 2:
            raise RuntimeError("Could not read two-column Sd spectrum: " + path)
        result[record] = (data[:, 0].astype(float), data[:, 1].astype(float))
    return result


T1, T_EFF = load_periods(T1_JSON) if os.path.isfile(T1_JSON) else (None, None)


def sd_record(record, period_s):
    periods, sd_m = SPEC[record]
    if period_s < periods[0] or period_s > periods[-1]:
        raise ValueError(
            "Target period {:.5f}s is outside {} spectrum [{:.5f}, {:.5f}]s".format(
                period_s, record, periods[0], periods[-1]
            )
        )
    return float(np.interp(period_s, periods, sd_m))


# ============================== INPUT SIGNAL =============================
def grid_period(period_s):
    n = max(2, int(round(period_s / DELTA_T)))
    return n * DELTA_T, n


def one_velocity_cycle(period_s, amplitude_mps):
    """One zero-to-zero sine cycle, sampled on the 3DEC input-table grid."""
    period_grid, n = grid_period(period_s)
    t = np.arange(n + 1, dtype=float) * DELTA_T
    v = amplitude_mps * np.sin(2.0 * math.pi * t / period_grid)
    return t, v, period_grid


def one_velocity_train(period_s, amplitude_mps, cycles=1, phase=1.0):
    """Integer-cycle sine packet at a grid-compatible carrier period."""
    if int(cycles) < 1:
        raise ValueError("cycles must be a positive integer")
    period_grid, n = grid_period(period_s)
    t = np.arange(int(cycles) * n + 1, dtype=float) * DELTA_T
    v = float(phase) * amplitude_mps * np.sin(2.0 * math.pi * t / period_grid)
    return t, v, period_grid


def combine_cycles(period_a_s, amp_a, period_b_s, amp_b,
                   cycles_b=1, phase_b=1.0):
    """One T1 cycle followed by a record-period Wave B packet."""
    _, va, period_a_grid = one_velocity_cycle(period_a_s, amp_a)
    _, vb, period_b_grid = one_velocity_train(
        period_b_s, amp_b, cycles=cycles_b, phase=phase_b
    )
    # Both packets begin/end at zero velocity; omit the duplicate join sample.
    velocity = np.concatenate((va, vb[1:]))
    time = np.arange(len(velocity), dtype=float) * DELTA_T
    return time, velocity, period_a_grid, period_b_grid


def sd_spectrum(acceleration, dt, periods_s, damping=XI):
    """5%-damped relative-displacement spectrum, Newmark average acceleration.

    The acceleration record includes a zero-acceleration free-vibration tail
    matching the model's TAIL_SEC.
    """
    ag = np.asarray(acceleration, dtype=float)
    periods = np.atleast_1d(np.asarray(periods_s, dtype=float))
    if len(ag) == 0 or np.any(periods <= 0.0):
        raise ValueError("Non-empty acceleration and positive periods are required")

    omega = 2.0 * math.pi / periods
    c = 2.0 * damping * omega
    omega2 = omega * omega
    u = np.zeros_like(omega)
    v = np.zeros_like(omega)
    a = -ag[0] - c * v - omega2 * u
    umax = np.abs(u)

    for i in range(1, len(ag)):
        u_predict = u + dt * v + 0.25 * dt * dt * a
        v_predict = v + 0.5 * dt * a
        a_new = (-ag[i] - c * v_predict - omega2 * u_predict) / (
            1.0 + 0.5 * dt * c + 0.25 * dt * dt * omega2
        )
        u = u_predict + 0.25 * dt * dt * a_new
        v = v_predict + 0.5 * dt * a_new
        a = a_new
        umax = np.maximum(umax, np.abs(u))
    return umax


def signal_sd(velocity, periods_s):
    acceleration = np.gradient(velocity, DELTA_T)
    tail_n = int(round(TAIL_SEC / DELTA_T))
    if tail_n > 0:
        acceleration = np.concatenate((acceleration, np.zeros(tail_n)))
    return sd_spectrum(acceleration, DELTA_T, periods_s, XI)


def fit_waveform_amplitudes(period_a_s, period_b_s, cycles_b, phase_b,
                            target_a_m, target_b_m,
                            target_period_a_s, target_period_b_s):
    """Fit two positive amplitudes with a bounded 1-D spectral-ratio solve.

    Scaling both amplitudes scales both Sd ordinates linearly. The unknown
    amplitude ratio is therefore solved first; amplitude A then sets the
    absolute scale. A bracket scan plus bisection avoids the fragile local
    2-D Newton update used by the earlier version.
    """
    target = np.asarray([target_a_m, target_b_m], dtype=float)
    response_periods = np.asarray([target_period_a_s, target_period_b_s], dtype=float)
    if np.any(target <= 0.0):
        raise ValueError("Both spectral displacement targets must be positive")

    _, va_unit, _ = one_velocity_cycle(period_a_s, 1.0)
    _, vb_unit, _ = one_velocity_train(
        period_b_s, 1.0, cycles=cycles_b, phase=phase_b
    )
    unit_a = float(signal_sd(va_unit, [target_period_a_s])[0])
    unit_b = float(signal_sd(vb_unit, [target_period_b_s])[0])
    if unit_a <= 0.0 or unit_b <= 0.0:
        raise RuntimeError("A unit-cycle calibration produced zero spectral displacement")
    q0 = (target[1] / unit_b) / (target[0] / unit_a)
    if not math.isfinite(q0) or q0 <= 0.0:
        raise RuntimeError("Could not initialize the positive amplitude ratio")
    target_ratio = target[1] / target[0]
    log_q = np.linspace(
        math.log(q0) - 8.0,
        math.log(q0) + 8.0,
        FIT_RATIO_GRID_POINTS,
    )
    evaluations = 0

    def evaluate(log_ratio):
        nonlocal evaluations
        q = math.exp(float(log_ratio))
        _, velocity_unit, _, _ = combine_cycles(
            period_a_s, 1.0, period_b_s, q,
            cycles_b=cycles_b, phase_b=phase_b,
        )
        achieved_unit = signal_sd(velocity_unit, response_periods)
        achieved_unit = np.maximum(achieved_unit, 1.0e-15)
        residual = math.log((achieved_unit[1] / achieved_unit[0]) / target_ratio)
        evaluations += 1
        return residual, velocity_unit, achieved_unit

    samples = [evaluate(q)[0] for q in log_q]
    root_logs = []
    for i in range(len(log_q) - 1):
        if abs(samples[i]) <= 1.0e-12:
            root_logs.append(float(log_q[i]))
        if samples[i] * samples[i + 1] < 0.0:
            lo = float(log_q[i])
            hi = float(log_q[i + 1])
            f_lo = float(samples[i])
            for _ in range(FIT_BISECT_ITER):
                mid = 0.5 * (lo + hi)
                f_mid = evaluate(mid)[0]
                if f_lo * f_mid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    f_lo = f_mid
            root_logs.append(0.5 * (lo + hi))
    if abs(samples[-1]) <= 1.0e-12:
        root_logs.append(float(log_q[-1]))

    if root_logs:
        selected_log_q = min(root_logs, key=lambda x: abs(x - math.log(q0)))
    else:
        # A sampled ratio can still meet the declared tolerance without an
        # exact crossing (for example, a tangent response-ratio minimum).
        best_index = int(np.argmin(np.abs(samples)))
        if abs(math.exp(float(samples[best_index])) - 1.0) > FIT_TOL:
            best_residual, best_velocity, best_sd = evaluate(float(log_q[best_index]))
            amp_a = target[0] / best_sd[0]
            amp_b = amp_a * math.exp(float(log_q[best_index]))
            velocity = best_velocity * amp_a
            achieved = signal_sd(velocity, response_periods)
            return {
                "amp_a_mps": float(amp_a),
                "amp_b_mps": float(amp_b),
                "velocity": velocity,
                "achieved_a_m": float(achieved[0]),
                "achieved_b_m": float(achieved[1]),
                "error_a_pct": float((achieved[0] / target[0] - 1.0) * 100.0),
                "error_b_pct": float((achieved[1] / target[1] - 1.0) * 100.0),
                "fit_converged": False,
                "fit_iterations": int(evaluations),
            }
        selected_log_q = float(log_q[best_index])

    _, velocity_unit, achieved_unit = evaluate(selected_log_q)
    amp_a = target[0] / achieved_unit[0]
    amp_b = amp_a * math.exp(selected_log_q)
    velocity = velocity_unit * amp_a
    achieved = signal_sd(velocity, response_periods)
    error_a = float((achieved[0] / target[0] - 1.0) * 100.0)
    error_b = float((achieved[1] / target[1] - 1.0) * 100.0)
    converged = max(abs(error_a), abs(error_b)) <= FIT_TOL * 100.0
    return {
        "amp_a_mps": float(amp_a),
        "amp_b_mps": float(amp_b),
        "velocity": velocity,
        "achieved_a_m": float(achieved[0]),
        "achieved_b_m": float(achieved[1]),
        "error_a_pct": error_a,
        "error_b_pct": error_b,
        "fit_converged": bool(converged),
        "fit_iterations": int(evaluations),
    }


def choose_record_waveform(record, target_period_a_s, target_period_b_s):
    """Choose the simplest feasible Tp-carrier packet for one record.

    Search order is fixed before the IDA sequence: fewest Wave B cycles first,
    then same relative polarity, then inverted polarity. This only uses record
    spectra and the fixed target periods.
    """
    tp = RECORD_TP[record]
    target_a = sd_record(record, target_period_a_s)
    target_b = sd_record(record, target_period_b_s)
    for cycles_b in range(1, MAX_WAVE_B_CYCLES + 1):
        for phase_b in WAVE_B_PHASE_ORDER:
            fit = fit_waveform_amplitudes(
                period_a_s=target_period_a_s,
                period_b_s=tp,
                cycles_b=cycles_b,
                phase_b=phase_b,
                target_a_m=target_a,
                target_b_m=target_b,
                target_period_a_s=target_period_a_s,
                target_period_b_s=target_period_b_s,
            )
            if fit["fit_converged"]:
                fit["cycles_b"] = int(cycles_b)
                fit["phase_b"] = float(phase_b)
                fit["phase_label"] = "same" if phase_b > 0.0 else "inverted"
                fit["carrier_b_s"] = float(tp)
                return fit
    raise RuntimeError(
        "No Wave B cycle/phase candidate met {:.1f}% for {} with carrier Tp={:.3f}s "
        "and at most {} cycles. Increase the candidate set only with a documented rule.".format(
            FIT_TOL * 100.0, record, tp, MAX_WAVE_B_CYCLES
        )
    )

def write_velocity_table(path, label, velocity):
    time = np.arange(len(velocity), dtype=float) * DELTA_T
    with open(path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(label, len(time)))
        for t, v in zip(time, velocity):
            f.write("{:.6f}\t{:.9e}\n".format(t, v))
    return float(time[-1])


def pga_g(velocity):
    return float(np.max(np.abs(np.gradient(velocity, DELTA_T))) / 9.81)


# =============================== 3DEC MODEL ===============================
def restore_model_state(save_file, clear_saved_histories=False):
    """Restore a sequence state and apply the standard Route 2 setup."""
    save_no_ext = os.path.splitext(save_file)[0]
    it.command("model restore '{}'".format(cmd_path(save_no_ext)))
    it.command("python-reset-state false")
    it.command("model large-strain off")
    it.command("model dynamic active on")
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
    if clear_saved_histories:
        it.command("history delete")
    it.command("call '{}'".format(cmd_path(HISTORY_SETUP)))
    it.command("call '{}'".format(cmd_path(BASESHEAR_SETUP)))
    it.command("block free velocity-z range group 'S'")
    it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")


def read_history_values(path):
    data = np.genfromtxt(path, skip_header=2)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError("Could not read exported history: " + path)
    valid = np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])
    return data[valid, 1]


def peak_relative_displacement_mm(csv_path):
    disp_m = read_history_values(csv_path)
    if len(disp_m) == 0:
        return float("nan")
    return float(np.max(np.abs(disp_m - disp_m[0])))


def preflight_method2(t1_s, teff_s):
    """Choose each record's packet, then verify all 25 scaled inputs."""
    designs = {}
    for record in SPECTRA:
        designs[record] = choose_record_waveform(record, t1_s, teff_s)
        design = designs[record]
        print("Waveform choice for {}: Tp={:.3f}s, Wave B cycles={}, phase={}; "
              "unit-scale errors A/B={:+.3f}%/{:+.3f}%".format(
                  record, design["carrier_b_s"], design["cycles_b"],
                  design["phase_label"], design["error_a_pct"],
                  design["error_b_pct"]))

    prepared = {}
    rows = []
    for run_no, record, scale in RUNS:
        design = designs[record]
        target_a_m = scale * sd_record(record, t1_s)
        target_b_m = scale * sd_record(record, teff_s)
        velocity = design["velocity"] * scale
        achieved = signal_sd(velocity, [t1_s, teff_s])
        error_a_pct = float((achieved[0] / target_a_m - 1.0) * 100.0)
        error_b_pct = float((achieved[1] / target_b_m - 1.0) * 100.0)
        fit_ok = max(abs(error_a_pct), abs(error_b_pct)) <= FIT_TOL * 100.0
        signal = {
            "amp_a_mps": design["amp_a_mps"] * scale,
            "amp_b_mps": design["amp_b_mps"] * scale,
            "velocity": velocity,
            "achieved_a_m": float(achieved[0]),
            "achieved_b_m": float(achieved[1]),
            "error_a_pct": error_a_pct,
            "error_b_pct": error_b_pct,
            "fit_converged": bool(fit_ok),
            "fit_iterations": design["fit_iterations"],
        }
        prepared[run_no] = {
            "target_a_m": target_a_m,
            "target_b_m": target_b_m,
            "target_period_b_s": teff_s,
            "carrier_b_s": design["carrier_b_s"],
            "cycles_b": design["cycles_b"],
            "phase_b": design["phase_b"],
            "signal": signal,
        }
        duration_s = (len(velocity) - 1) * DELTA_T
        rows.append({
            "run": run_no,
            "record": record,
            "scale": scale,
            "T1_target_s": t1_s,
            "Teff_target_s": teff_s,
            "WaveB_carrier_Tp_s": design["carrier_b_s"],
            "WaveB_cycles": design["cycles_b"],
            "WaveB_phase": design["phase_label"],
            "target_A_mm": target_a_m * 1000.0,
            "achieved_A_mm": float(achieved[0]) * 1000.0,
            "error_A_pct": error_a_pct,
            "target_B_mm": target_b_m * 1000.0,
            "achieved_B_mm": float(achieved[1]) * 1000.0,
            "error_B_pct": error_b_pct,
            "fit_converged": bool(fit_ok),
            "fit_search_evaluations": design["fit_iterations"],
            "amp_A_mps": signal["amp_a_mps"],
            "amp_B_mps": signal["amp_b_mps"],
            "combined_PGA_g": pga_g(velocity),
            "input_duration_s": duration_s,
        })

    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    preflight_path = os.path.join(OUTPUT_ROOT, FIT_PREFLIGHT_NAME)
    with open(preflight_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    failures = [row for row in rows if not row["fit_converged"]]
    print("Method 2 fit preflight: {}/{} runs meet the {:.1f}% tolerance.".format(
        len(rows) - len(failures), len(rows), FIT_TOL * 100.0))
    print("  full fit table: {}".format(preflight_path))
    pga_exceedances = [row for row in rows if row["combined_PGA_g"] > PGA_WARN_G]
    if pga_exceedances:
        print("  WARNING: predicted input PGA exceeds {:.2f} g in {} run(s): {}".format(
            PGA_WARN_G, len(pga_exceedances),
            ", ".join("{:02d}={:.3f}g".format(r["run"], r["combined_PGA_g"])
                      for r in pga_exceedances)))
    if failures and FIT_ABORT_ON_FAILURE:
        raise RuntimeError(
            "Fit preflight failed for {} run(s); no model state was restored. Review {}.".format(
                len(failures), preflight_path
            )
        )
    return prepared

def protocol_signature():
    return [[int(n), str(record), float(scale)] for n, record, scale in RUNS]


def write_summary_log(path, rows):
    """Rewrite the small run log from checkpoint data, preventing duplicates."""
    if not rows:
        return
    temp_path = path + ".tmp"
    with open(temp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp_path, path)


def save_checkpoint(path, method_name, last_run, summary_rows):
    temp_path = path + ".tmp"
    state = {
        "method": method_name,
        "protocol": protocol_signature(),
        "last_completed_run": int(last_run),
        "summary": summary_rows,
    }
    with open(temp_path, "w") as f:
        json.dump(state, f, indent=2, allow_nan=False)
    os.replace(temp_path, path)


def load_checkpoint(path, method_name):
    if not os.path.isfile(path):
        return 0, []
    with open(path, "r") as f:
        state = json.load(f)
    if state.get("method") != method_name or state.get("protocol") != protocol_signature():
        raise RuntimeError(
            "Checkpoint does not match this method/protocol. Choose a clean OUTPUT_ROOT: " + path
        )
    last_run = int(state.get("last_completed_run", 0))
    completed = [n for n, _, _ in RUNS if n <= last_run]
    summary_rows = state.get("summary", [])
    if (not completed or completed[-1] != last_run or
            len(summary_rows) != len(completed)):
        raise RuntimeError("Checkpoint is incomplete or inconsistent: " + path)
    return last_run, summary_rows


def run_method(method_name, t1_s, teff_s, prepared):
    method_dir = os.path.join(OUTPUT_ROOT, method_name)
    os.makedirs(method_dir, exist_ok=True)
    log_path = os.path.join(method_dir, "route2_method_log.csv")
    checkpoint_path = os.path.join(method_dir, CHECKPOINT_NAME)
    summary_rows = []
    last_done = 0

    print("\n" + "=" * 76)
    print("Method {} | target=fixed | B carrier=record 5%-damped Sa peak".format(method_name))
    print("Initial state: {}".format(BASE_SAVE))
    print("Outputs: {}".format(method_dir))
    print("=" * 76)

    if os.path.isfile(checkpoint_path):
        last_done, summary_rows = load_checkpoint(checkpoint_path, method_name)
        if last_done == RUNS[-1][0]:
            write_summary_log(log_path, summary_rows)
            print("Checkpoint already contains all {} runs; no model solve needed.".format(last_done))
            return summary_rows
        resume_save = os.path.join(
            method_dir, "route2_{}_run_{:02d}.sav".format(method_name, last_done)
        )
        if not os.path.isfile(resume_save):
            raise RuntimeError("Checkpoint state save is missing: " + resume_save)
        restore_model_state(resume_save, clear_saved_histories=True)
        print("Resuming after run {} from {}".format(last_done, resume_save))
        write_summary_log(log_path, summary_rows)
    else:
        stale_outputs = [
            name for name in os.listdir(method_dir)
            if name.endswith(".sav") or name.endswith("_vel.txt")
        ]
        if os.path.isfile(log_path) or stale_outputs:
            raise RuntimeError(
                "Partial outputs exist without a checkpoint. Choose a new OUTPUT_ROOT "
                "or inspect the directory before restarting: " + method_dir
            )
        restore_model_state(BASE_SAVE, clear_saved_histories=False)
        print("Fresh sequential start from the initial model save.")

    for run_no, record, scale in RUNS:
        if run_no <= last_done:
            continue

        run_input = prepared[run_no]
        target_a_m = run_input["target_a_m"]
        target_b_m = run_input["target_b_m"]
        target_period_b_s = run_input["target_period_b_s"]
        carrier_b_s = run_input["carrier_b_s"]
        cycles_b = run_input["cycles_b"]
        phase_b = run_input["phase_b"]
        signal = run_input["signal"]

        label = "Run{:02d}_{}_s{}_{}".format(
            run_no, record, "{:.2f}".format(scale).replace(".", "p"), method_name
        )
        table_name = "r2_{}_{}".format(method_name, run_no)
        table_path = os.path.join(method_dir, label + "_vel.txt")
        duration_s = write_velocity_table(table_path, label, signal["velocity"])
        pga = pga_g(signal["velocity"])

        print("\nRun {:02d}: {} x {:.2f}".format(run_no, record, scale))
        print("  A: carrier {:.5f}s; Sd target {:.3f} mm".format(
            grid_period(t1_s)[0], target_a_m * 1000.0))
        print("  B: carrier {:.5f}s (record Sa peak), {} cycles, {} phase; "
              "target oscillator {:.5f}s; Sd target {:.3f} mm".format(
                  grid_period(carrier_b_s)[0], cycles_b,
                  "same" if phase_b > 0.0 else "inverted",
                  target_period_b_s, target_b_m * 1000.0))
        print("  amplitudes A/B: {:.6g} / {:.6g} m/s".format(
            signal["amp_a_mps"], signal["amp_b_mps"]))
        print("  combined fit errors A/B: {:+.2f}% / {:+.2f}% ({})".format(
            signal["error_a_pct"], signal["error_b_pct"],
            "converged" if signal["fit_converged"] else "NOT converged"))
        print("  combined PGA: {:.3f} g; input {:.3f}s + {:.1f}s tail".format(
            pga, duration_s, TAIL_SEC))
        if pga > PGA_WARN_G:
            print("  WARNING: PGA exceeds {:.2f} g.".format(PGA_WARN_G))

        # These are the same Itasca Python API commands used by the existing
        # Route 2 driver. The wall state is retained; only dynamic time resets.
        it.command("table '{}' import '{}'".format(table_name, cmd_path(table_path)))
        it.command("model dynamic time-total 0")
        for group in DRIVE_GROUPS:
            it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(
                table_name, group
            ))
        it.command("model solve dynamic time {:.6f}".format(duration_s))
        for group in DRIVE_GROUPS:
            it.command("block gridpoint apply-remove velocity-z range group '{}'".format(group))
        if TAIL_SEC > 0.0:
            it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

        # Export this run's histories, then clear/re-register them before the
        # checkpoint save so a resumed run starts with a clean history set.
        it.command("[exportdir='{}']".format(cmd_path(method_dir)))
        it.command("[runlabel='{}']".format(label))
        it.command("call '{}'".format(cmd_path(HISTORY_EXPORT)))
        disp_path = os.path.join(method_dir, label + "_rel_disp_top_exp_mm.csv")
        peak_mm = peak_relative_displacement_mm(disp_path)
        peak_value = peak_mm if math.isfinite(peak_mm) else None

        it.command("table '{}' delete".format(table_name))
        it.command("history delete")
        it.command("call '{}'".format(cmd_path(HISTORY_SETUP)))
        it.command("call '{}'".format(cmd_path(BASESHEAR_SETUP)))

        save_path = os.path.join(method_dir, "route2_{}_run_{:02d}.sav".format(
            method_name, run_no
        ))
        it.command("model save '{}'".format(cmd_path(save_path)))

        row = {
            "method": method_name,
            "target_mode": "fixed",
            "carrier_mode": "record_Sa_peak_feasibility_selected",
            "wave_B_cycles": cycles_b,
            "wave_B_phase": "same" if phase_b > 0.0 else "inverted",
            "run": run_no,
            "record": record,
            "scale": scale,
            "carrier_A_requested_s": t1_s,
            "carrier_A_grid_s": grid_period(t1_s)[0],
            "carrier_B_requested_s": carrier_b_s,
            "carrier_B_grid_s": grid_period(carrier_b_s)[0],
            "target_A_period_s": t1_s,
            "target_B_period_s": target_period_b_s,
            "target_A_mm": target_a_m * 1000.0,
            "target_B_mm": target_b_m * 1000.0,
            "achieved_A_mm": signal["achieved_a_m"] * 1000.0,
            "achieved_B_mm": signal["achieved_b_m"] * 1000.0,
            "error_A_pct": signal["error_a_pct"],
            "error_B_pct": signal["error_b_pct"],
            "fit_converged": signal["fit_converged"],
            "fit_iterations": signal["fit_iterations"],
            "amp_A_mps": signal["amp_a_mps"],
            "amp_B_mps": signal["amp_b_mps"],
            "combined_PGA_g": pga,
            "table_duration_s": duration_s,
            "peak_relative_OOP_mm": peak_value,
            "save_file": save_path,
        }
        summary_rows.append(row)
        write_summary_log(log_path, summary_rows)
        save_checkpoint(checkpoint_path, method_name, run_no, summary_rows)
        print("  peak relative OOP displacement: {}".format(
            "{:.3f} mm".format(peak_value) if peak_value is not None else "unavailable"))
        print("  saved state and checkpoint: {}".format(save_path))

    return summary_rows

def main():
    method_name, _, _ = METHOD
    method_dir = os.path.join(OUTPUT_ROOT, method_name)
    checkpoint_path = os.path.join(method_dir, CHECKPOINT_NAME)

    required = [
        T1_JSON,
        HISTORY_SETUP,
        BASESHEAR_SETUP,
        HISTORY_EXPORT,
    ] + list(SPECTRA.values())
    require_files(required)

    t1_s, teff_s = load_periods(T1_JSON)
    if set(RECORD_TP) != set(SPECTRA):
        raise RuntimeError("RECORD_TP must define exactly the records with spectra")
    for record, period_s in RECORD_TP.items():
        if not math.isfinite(period_s) or period_s <= 0.0:
            raise RuntimeError("Invalid predominant period for {}: {}".format(record, period_s))

    global SPEC
    SPEC = load_spectra()
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    prepared = preflight_method2(t1_s, teff_s)

    # Validate the initial save only after the signal preflight. This lets the
    # fit table be generated on a workstation that has the spectra but not the
    # large 3DEC save, and guarantees a bad fit cannot touch model state.
    if not os.path.isfile(checkpoint_path):
        require_files([BASE_SAVE])

    # Preserve the Python namespace across model restores, as used by the
    # existing SAFEGo scripts in the embedded Itasca Python environment.
    it.command("python-reset-state false")
    it.command("program automatic-model-save active off")
    summary = run_method(method_name, t1_s, teff_s, prepared)

    print("\nMethod 2 full-sequence output: {}".format(OUTPUT_ROOT))
    print("Run log: {}".format(os.path.join(method_dir, "route2_method_log.csv")))
    print("Fit preflight: {}".format(os.path.join(OUTPUT_ROOT, FIT_PREFLIGHT_NAME)))
    print("Base-shear/displacement postprocess after completion:")
    print('  python "Claude outputs/baseshear_exp_from_model.py" "{}" {}'.format(
        method_dir, " ".join(str(run[0]) for run in RUNS)
    ))
    print("Completed/logged runs: {} of {}".format(len(summary), len(RUNS)))


main()
