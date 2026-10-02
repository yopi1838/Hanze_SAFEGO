# -*- coding: utf-8 -*-
"""Run the four Route 2 carrier/target combinations from Strategy F run 21.

This is a self-contained per-method driver. The four distributed copies differ
only in the METHOD tuple below. Run each copy in its own 3DEC Python console,
from the Hanze_SAFEGO project root.
It uses the embedded Itasca API (``import itasca as it`` and ``it.command``)
and runs four independent FR76 sequences, each restored from the same saved
Strategy F run-21 state:

    run 22: FR76 x1.00
    run 23: FR76 x1.50
    run 24: FR76 x1.75
    run 25: FR76 x2.00

The four methods form a 2 x 2 comparison:
    fixed_wall       B carrier Teff; target scale * Sd(Teff)
    fixed_record     B carrier Tp;   target scale * Sd(Teff)
    adaptive_wall    B carrier T*;   target d* from casm_step1_targets.csv
    adaptive_record  B carrier Tp;   target d* from casm_step1_targets.csv

Wave A is always carried at T1 and targets scale * Sd(T1). The two wave
amplitudes are solved jointly against the combined velocity history. This
corrects the separate-amplitude leakage in the original Route 2. It means
fixed_wall is the joint-calibrated reference, while the old as-run Route 2
results remain a separate historical control.

All history CSVs are written flat into each method directory so that
Claude outputs/baseshear_exp_from_model.py can read them directly. Each method
restores the same baseline save before run 22, then carries its own damage
state through runs 23-25.

IMPORTANT: the expected Strategy F run-21 save is
    stratF_full_results_US1/stratC_run_21.sav
The repository copy of stratF_full_results_US1 may contain exports without
the save itself. Override BASE_SAVE below if the save is elsewhere.
"""

import csv
import json
import math
import os
import sys

import numpy as np
import itasca as it


# ================================ CONFIG ================================
BASE_SAVE = os.path.join("stratF_full_results_US1", "stratC_run_21.sav")
OUTPUT_ROOT = "route2_four_methods"

T1_JSON = os.path.join("bilinear_idealisation", "bilinear_periods.json")
TARGETS_CSV = os.path.join("Claude outputs", "casm_step1_targets.csv")
SPECTRA = {"FR76": "spectrum_FR76.csv"}

HISTORY_SETUP = "instrument_history_new.dat"
BASESHEAR_SETUP = "instrument_baseshear_exp.dat"
HISTORY_EXPORT = "instrument_history_export_v2.dat"

RECORD_TP = {"FR76": 0.260}  # 5%-damped Sa-peak period used in the handover
RUNS = [
    (22, "FR76", 1.00),
    (23, "FR76", 1.50),
    (24, "FR76", 1.75),
    (25, "FR76", 2.00),
]
# Each standalone copy has exactly one method tuple here:
# (output_name, target_mode, carrier_mode)
METHOD = ("4_adaptive_record", "adaptive", "record")
METHODS = [METHOD]

XI = 0.05
DELTA_T = 0.005
TAIL_SEC = 2.5
STAGE_B_PHASE = 1.0   # same-phase cycles, matching the current Route 2 setup
PGA_WARN_G = 0.90
FIT_TOL = 0.02        # relative spectral target tolerance for the joint solver
FIT_MAX_ITER = 30
FIT_JAC_STEP = 0.03

DRIVE_GROUPS = ["S", "T_B"]
KEY_REL_DISP = "rel_disp_top_exp_mm"


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
                "The Strategy F run-21 save is not in this workspace. "
                "Copy it here or set BASE_SAVE to its full path."
            )
        raise RuntimeError("\n".join(lines))


def load_periods(path):
    with open(path, "r") as f:
        data = json.load(f)
    model = data["model"]
    return float(model["T1_s"]), float(model["Teff_s"])


def load_fixed_point_targets(path):
    rows = {}
    with open(path, "r", newline="") as f:
        for row in csv.DictReader(f):
            run_no = int(row["run"])
            if run_no in [x[0] for x in RUNS]:
                d_star = float(row["d_star_mm"]) * 1.0e-3
                t_star = float(row["T_star_s"])
                if not (math.isfinite(d_star) and d_star > 0.0 and
                        math.isfinite(t_star) and t_star > 0.0):
                    raise RuntimeError(
                        "Invalid capacity-adaptive target for run {}".format(run_no)
                    )
                rows[run_no] = (d_star, t_star)
    absent = [run_no for run_no, _, _ in RUNS if run_no not in rows]
    if absent:
        raise RuntimeError("No finite d*/T* target for runs {} in {}".format(absent, path))
    return rows


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


def combine_cycles(period_a_s, amp_a, period_b_s, amp_b):
    ta, va, period_a_grid = one_velocity_cycle(period_a_s, amp_a)
    tb, vb, period_b_grid = one_velocity_cycle(
        period_b_s, STAGE_B_PHASE * amp_b
    )
    # Both cycles start/end at zero velocity; omit the duplicate join sample.
    velocity = np.concatenate((va, vb[1:]))
    time = np.arange(len(velocity), dtype=float) * DELTA_T
    return time, velocity, period_a_grid, period_b_grid


def sd_spectrum(acceleration, dt, periods_s, damping=XI):
    """5%-damped relative-displacement spectrum, Newmark average acceleration.

    This is the same integration scheme used by the current Route 2 driver.
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


def solve_two_amplitudes(period_a_s, period_b_s, target_a_m, target_b_m,
                         target_period_a_s, target_period_b_s):
    """Fit both oscillator targets to the *combined* two-cycle signal.

    Uses a small finite-difference Newton solve in log-amplitude coordinates.
    The solve is dependency-free beyond NumPy, which is present in the 3DEC
    scientific Python distribution. A non-converged fit is returned and
    recorded so an infeasible target is visible rather than silently hidden.
    """
    target = np.asarray([target_a_m, target_b_m], dtype=float)
    response_periods = np.asarray([target_period_a_s, target_period_b_s], dtype=float)
    if np.any(target <= 0.0):
        raise ValueError("Both spectral displacement targets must be positive")

    # Isolated-wave amplitudes give a stable initial point for the joint solve.
    _, va_unit, _ = one_velocity_cycle(period_a_s, 1.0)
    _, vb_unit, _ = one_velocity_cycle(period_b_s, 1.0)
    unit_a = signal_sd(va_unit, response_periods)[0]
    unit_b = signal_sd(vb_unit, response_periods)[1]
    if unit_a <= 0.0 or unit_b <= 0.0:
        raise RuntimeError("A unit-cycle calibration produced zero spectral displacement")
    x = np.log(np.asarray([target_a_m / unit_a, target_b_m / unit_b]))

    def evaluate(log_amplitudes):
        amp_a, amp_b = np.exp(log_amplitudes)
        _, velocity, _, _ = combine_cycles(period_a_s, amp_a, period_b_s, amp_b)
        achieved = signal_sd(velocity, response_periods)
        achieved = np.maximum(achieved, 1.0e-15)
        residual = np.log(achieved / target)
        return residual, achieved, velocity, amp_a, amp_b

    converged = False
    iterations = 0
    for iterations in range(1, FIT_MAX_ITER + 1):
        residual, achieved, velocity, amp_a, amp_b = evaluate(x)
        base_norm = float(np.linalg.norm(residual))
        if float(np.max(np.abs(residual))) <= FIT_TOL:
            converged = True
            break

        jac = np.zeros((2, 2), dtype=float)
        for j in range(2):
            trial = x.copy()
            trial[j] += FIT_JAC_STEP
            trial_residual, _, _, _, _ = evaluate(trial)
            jac[:, j] = (trial_residual - residual) / FIT_JAC_STEP

        try:
            step = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            step = -0.35 * residual
        step = np.clip(step, -1.25, 1.25)

        accepted = False
        for shrink in (1.0, 0.5, 0.25, 0.125, 0.0625):
            trial = x + shrink * step
            trial_residual, _, _, _, _ = evaluate(trial)
            if float(np.linalg.norm(trial_residual)) < base_norm:
                x = trial
                accepted = True
                break

        if not accepted:
            # Damped diagonal correction is a safe fallback when peak switching
            # makes the local finite-difference Jacobian unreliable.
            trial = x - 0.15 * residual
            trial_residual, _, _, _, _ = evaluate(trial)
            if float(np.linalg.norm(trial_residual)) < base_norm:
                x = trial
            else:
                break

    residual, achieved, velocity, amp_a, amp_b = evaluate(x)
    if float(np.max(np.abs(residual))) <= FIT_TOL:
        converged = True
    return {
        "amp_a_mps": float(amp_a),
        "amp_b_mps": float(amp_b),
        "velocity": velocity,
        "achieved_a_m": float(achieved[0]),
        "achieved_b_m": float(achieved[1]),
        "error_a_pct": float((achieved[0] / target[0] - 1.0) * 100.0),
        "error_b_pct": float((achieved[1] / target[1] - 1.0) * 100.0),
        "fit_converged": bool(converged),
        "fit_iterations": int(iterations),
    }


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
def restore_strategy_f_run21():
    """Restore the common state and apply the same dynamic setup as Route 2."""
    save_no_ext = os.path.splitext(BASE_SAVE)[0]
    it.command("model restore '{}'".format(cmd_path(save_no_ext)))
    it.command("model large-strain off")
    it.command("model dynamic active on")
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
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


def run_method(method_name, target_mode, carrier_mode, fixed_targets,
               spectra, t1_s, teff_s):
    method_dir = os.path.join(OUTPUT_ROOT, method_name)
    os.makedirs(method_dir, exist_ok=True)
    log_path = os.path.join(method_dir, "route2_method_log.csv")
    save_no_ext = os.path.splitext(BASE_SAVE)[0]
    summary_rows = []

    print("\n" + "=" * 76)
    print("Method {} | target={} | B carrier={}".format(
        method_name, target_mode, carrier_mode))
    print("Restore common Strategy F state: {}".format(BASE_SAVE))
    print("Outputs: {}".format(method_dir))
    print("=" * 76)

    restore_strategy_f_run21()

    for run_no, record, scale in RUNS:
        if target_mode == "fixed":
            target_period_b_s = teff_s
            target_b_m = scale * sd_record(record, teff_s)
        else:
            target_b_m, target_period_b_s = fixed_targets[run_no]

        if carrier_mode == "wall":
            carrier_b_s = teff_s if target_mode == "fixed" else target_period_b_s
        else:
            carrier_b_s = RECORD_TP[record]

        target_a_m = scale * sd_record(record, t1_s)
        signal = solve_two_amplitudes(
            period_a_s=t1_s,
            period_b_s=carrier_b_s,
            target_a_m=target_a_m,
            target_b_m=target_b_m,
            target_period_a_s=t1_s,
            target_period_b_s=target_period_b_s,
        )

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
        print("  B: carrier {:.5f}s; target oscillator {:.5f}s; Sd target {:.3f} mm".format(
            grid_period(carrier_b_s)[0], target_period_b_s, target_b_m * 1000.0))
        print("  amplitudes A/B: {:.6g} / {:.6g} m/s".format(
            signal["amp_a_mps"], signal["amp_b_mps"]))
        print("  combined fit errors A/B: {:+.2f}% / {:+.2f}% ({})".format(
            signal["error_a_pct"], signal["error_b_pct"],
            "converged" if signal["fit_converged"] else "NOT converged"))
        print("  combined PGA: {:.3f} g; input {:.3f}s + {:.1f}s tail".format(
            pga, duration_s, TAIL_SEC))
        if pga > PGA_WARN_G:
            print("  WARNING: PGA exceeds {:.2f} g.".format(PGA_WARN_G))

        # The table and solve are entirely in the Itasca API. Names are unique
        # per method/run; model time is reset so every table starts at t = 0.
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

        save_path = os.path.join(method_dir, "route2_{}_run_{:02d}.sav".format(
            method_name, run_no
        ))
        it.command("model save '{}'".format(cmd_path(save_path)))

        # Keep the files flat in method_dir: the base-shear script searches
        # this directory for RunNN_* histories.
        exportdir = cmd_path(method_dir)
        runlabel = label
        it.command("[exportdir='{}']".format(exportdir))
        it.command("[runlabel='{}']".format(runlabel))
        it.command("call '{}'".format(cmd_path(HISTORY_EXPORT)))

        disp_path = os.path.join(method_dir, label + "_rel_disp_top_exp_mm.csv")
        peak_mm = peak_relative_displacement_mm(disp_path)
        row = {
            "method": method_name,
            "target_mode": target_mode,
            "carrier_mode": carrier_mode,
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
            "peak_relative_OOP_mm": peak_mm,
            "save_file": save_path,
        }
        summary_rows.append(row)
        with open(log_path, "w" if run_no == RUNS[0][0] else "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            if run_no == RUNS[0][0]:
                writer.writeheader()
            writer.writerow(row)
        print("  peak relative OOP displacement: {:.3f} mm".format(peak_mm))
        print("  saved state: {}".format(save_path))

        # Clear all run histories, then re-register displacement, tilt and
        # experiment-equivalent base-shear channels before the next run.
        it.command("table '{}' delete".format(table_name))
        it.command("history delete")
        it.command("call '{}'".format(cmd_path(HISTORY_SETUP)))
        it.command("call '{}'".format(cmd_path(BASESHEAR_SETUP)))

    return summary_rows


def main():
    required = [
        BASE_SAVE,
        T1_JSON,
        TARGETS_CSV,
        HISTORY_SETUP,
        BASESHEAR_SETUP,
        HISTORY_EXPORT,
    ] + list(SPECTRA.values())
    require_files(required)

    t1_s, teff_s = load_periods(T1_JSON)
    fixed_targets = load_fixed_point_targets(TARGETS_CSV)
    global SPEC
    SPEC = load_spectra()
    if "FR76" not in RECORD_TP:
        raise RuntimeError("Missing record predominant period for FR76")

    # Preserve the Python namespace across model restores, as recommended by
    # the 3DEC Python documentation and used by the existing SAFEGo drivers.
    it.command("python-reset-state false")
    it.command("program automatic-model-save active off")
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    for method_name, target_mode, carrier_mode in METHODS:
        run_method(
            method_name, target_mode, carrier_mode,
            fixed_targets, SPEC, t1_s, teff_s
        )

    print("\nCompleted method {}.".format(METHOD[0]))
    print("Run summary: {}".format(os.path.join(OUTPUT_ROOT, METHOD[0], "route2_method_log.csv")))
    print("For base-shear versus displacement loops, run the existing postprocessor")
    print("once for each method directory, for example:")
    for method_name, _, _ in METHODS:
        print('  python "Claude outputs/baseshear_exp_from_model.py" {} 22 23 24 25'.format(
            '"{}"'.format(os.path.join(OUTPUT_ROOT, method_name))
        ))


require_files([
    T1_JSON,
    TARGETS_CSV,
    HISTORY_SETUP,
    BASESHEAR_SETUP,
    HISTORY_EXPORT,
])

main()
