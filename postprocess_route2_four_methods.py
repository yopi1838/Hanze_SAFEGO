# -*- coding: utf-8 -*-
"""Summarize and compare the four Route 2 method runs.

Run from the Hanze_SAFEGO project root with the installed Itasca Python, e.g.:
    python postprocess_route2_four_methods.py

The base-shear histories are computed by importing the existing
Claude outputs/baseshear_exp_from_model.py functions, preserving its
50-Hz filter, sign convention, masses, and mass scale.
"""

import csv
import datetime as dt
import importlib.util
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path("route2_four_methods")
OUT_DIR = ROOT / "postproc"
METHODS = [
    ("1_fixed_wall", "Method 1"),
    ("2_fixed_record", "Method 2"),
    ("3_adaptive_wall", "Method 3"),
    ("4_adaptive_record", "Method 4"),
]
RUNS = (22, 23, 24, 25)
EXP_TILT_PATH = Path("Claude outputs") / "exp_Test9_incl_tilt.csv"
SUMMARY_CSV = OUT_DIR / "route2_four_method_metrics.csv"
US2_SUMMARY_CSV = OUT_DIR / "route2_four_method_us2_summary.csv"
TILT_INPUT_DIR = Path("EXP_DATA")
TILT_SPECIMEN_RUNS = {"US1": 24, "US2": 25}
US2_TILT_CSV = OUT_DIR / "exp_US2_tiltmeter.csv"
US1_TILT_CSV = OUT_DIR / "exp_US1_tiltmeter.csv"
SUMMARY_PNG = OUT_DIR / "fig_route_2_four_method_comparison.png"
FOCUS_PNG = OUT_DIR / "fig_route_2_four_method_comparison_method2_focus.png"
TILT_PNG = OUT_DIR / "fig_route2_four_method_tiltmeter_residual.png"
TAIL_FRAC = 0.05


def load_baseshear_module():
    path = Path("Claude outputs") / "baseshear_exp_from_model.py"
    spec = importlib.util.spec_from_file_location("baseshear_exp_from_model", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Keep this postprocessor pinned to the project copy of the test data.
    module.EXP_DIR = str(Path("EXP_DATA") / "processed_globalzero")
    return module


def read_history(path):
    data = np.genfromtxt(str(path), skip_header=2)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("Could not parse history: {}".format(path))
    valid = np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])
    data = data[valid]
    if len(data) < 5:
        raise ValueError("History has fewer than five samples: {}".format(path))
    return data[:, 0], data[:, 1]


def peak_from_start(values):
    values = np.asarray(values, dtype=float)
    return float(np.max(np.abs(values - values[0])))


def tail_mean(values):
    values = np.asarray(values, dtype=float)
    n = max(1, int(math.ceil(TAIL_FRAC * len(values))))
    return float(np.mean(values[-n:]))


def load_method_log(method_dir):
    path = method_dir / "route2_method_log.csv"
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return {int(row["run"]): row for row in rows}


def load_exp_tilt():
    if not EXP_TILT_PATH.is_file():
        raise FileNotFoundError("Experimental inclinometer summary missing: {}".format(EXP_TILT_PATH))
    with EXP_TILT_PATH.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return {int(row["run"]): row for row in rows}


def load_experimental_tilt(specimen):
    """Extract settled Y-tilt using one raw-log event rule for both specimens.

    The logger is sampled every 5 s with decimal commas. Shaking windows are
    detected from its accelerometer swing; the cumulative residual is the
    median settled value after each event relative to the first five samples.
    This reports the instrument itself, not a wall chord-angle proxy.
    """
    if specimen not in TILT_SPECIMEN_RUNS:
        raise ValueError("Unknown tilt specimen: {}".format(specimen))
    input_path = TILT_INPUT_DIR / "{}_Tilt_values.csv".format(specimen)
    if not input_path.is_file():
        raise FileNotFoundError("{} inclinometer log missing: {}".format(specimen, input_path))

    with input_path.open("r", newline="", encoding="utf-8-sig") as stream:
        raw = list(csv.DictReader(stream, delimiter=";"))

    def number(row, key):
        return float(row[key].strip().replace(",", "."))

    timestamps = [dt.datetime.strptime(row["time"], "%d-%m-%Y_%H:%M:%S") for row in raw]
    seconds = np.asarray([(value - timestamps[0]).total_seconds() for value in timestamps], dtype=float)
    tilt_y = np.asarray([number(row, "sensor1 tvaly") for row in raw], dtype=float)
    swing_x = np.asarray([number(row, "sensor1 amaxx1") - number(row, "sensor1 aminx1") for row in raw])
    swing_y = np.asarray([number(row, "sensor1 amaxy1") - number(row, "sensor1 aminy1") for row in raw])
    swing = np.maximum(swing_x, swing_y)
    threshold = float(np.median(swing) * 5.0)
    active = np.flatnonzero(swing > threshold)
    if len(active) == 0:
        raise ValueError("No shaking events detected in {}".format(input_path))

    events = []
    start = previous = int(active[0])
    for index in active[1:]:
        index = int(index)
        if seconds[index] - seconds[previous] > 30.0:
            events.append((start, previous))
            start = index
        previous = index
    events.append((start, previous))
    expected = TILT_SPECIMEN_RUNS[specimen]
    if len(events) != expected:
        raise ValueError("{} tilt log produced {} event windows; expected {}".format(specimen, len(events), expected))

    initial = float(np.median(tilt_y[:5]))
    rows = []
    results = {}
    for idx, (first, last) in enumerate(events):
        next_start = events[idx + 1][0] if idx + 1 < len(events) else len(tilt_y)
        settled = tilt_y[last + 1:next_start]
        residual = ((float(np.median(settled)) if len(settled) else tilt_y[last]) - initial) / 1000.0
        run = idx + 1
        rows.append({
            "run": run,
            "t_start_s": round(float(seconds[first]), 1),
            "duration_s": round(float(seconds[last] - seconds[first]), 1),
            "resid_tilt_y_deg": round(residual, 4),
        })
        results[run] = residual

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = US1_TILT_CSV if specimen == "US1" else US2_TILT_CSV
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return results


def percent_error(model_value, experiment_value):
    if experiment_value is None or abs(experiment_value) < 1.0e-12:
        return None
    return 100.0 * (model_value / experiment_value - 1.0)


def round_or_blank(value, digits=4):
    return "" if value is None or not np.isfinite(value) else round(float(value), digits)


def collect(baseshear):
    exp_tilt = load_exp_tilt()
    if 21 not in exp_tilt:
        raise ValueError("Experimental tilt summary lacks run 21 baseline")
    exp_tilt_base = float(exp_tilt[21]["resid_tilt_y_deg"])
    # Keep the two specimens as independent benchmarks. US-1 has no Run 25;
    # US-2 contains all four FR76 runs used in the model sequences.
    exp_runs = {
        "US1": {run: baseshear.exp_run("Test9", run) for run in RUNS},
        "US2": {run: baseshear.exp_run("Test12", run) for run in RUNS},
    }
    all_rows = []
    model_tilt_curves = {}
    exp_tilt_curve = {}

    for method, label in METHODS:
        method_dir = ROOT / method
        log = load_method_log(method_dir)
        displacement_path_by_run = {}
        local_tilt_path_by_run = {}
        for run in RUNS:
            disp_hits = sorted(method_dir.glob("Run{:02d}_*rel_disp_top_exp_mm.csv".format(run)))
            tilt_hits = sorted(method_dir.glob("Run{:02d}_*tilt_local_incl.csv".format(run)))
            if not disp_hits or not tilt_hits:
                raise FileNotFoundError("Missing displacement or local-tilt history for {} run {}".format(method, run))
            displacement_path_by_run[run] = disp_hits[0]
            local_tilt_path_by_run[run] = tilt_hits[0]

        _, first_tilt = read_history(local_tilt_path_by_run[22])
        model_tilt_base = float(first_tilt[0])
        for run in RUNS:
            row_log = log[run]
            _, disp = read_history(displacement_path_by_run[run])
            disp_peak = peak_from_start(disp)  # rel_disp_top_exp_mm is already top-minus-table, in mm
            logged_peak = float(row_log["peak_relative_OOP_mm"])
            if not math.isclose(disp_peak, logged_peak, rel_tol=0.0, abs_tol=0.02):
                raise ValueError("Recomputed peak {:.3f} mm disagrees with run log {:.3f} mm for {} run {}".format(
                    disp_peak, logged_peak, method, run))

            t_model_shear, model_shear = baseshear.model_base_shear(str(method_dir), run)
            if model_shear is None:
                raise FileNotFoundError("Missing three velocity exports needed for {} run {} base shear".format(method, run))
            model_shear = model_shear - model_shear[0]
            model_shear_peak = float(np.max(np.abs(model_shear)))

            exp_peaks = {}
            for specimen, histories in exp_runs.items():
                exp_u, exp_shear = histories[run]
                if exp_u is None:
                    exp_peaks[specimen] = (None, None)
                else:
                    exp_peaks[specimen] = (
                        float(np.max(np.abs(exp_u))),
                        float(np.max(np.abs(exp_shear))),
                    )

            _, local_tilt = read_history(local_tilt_path_by_run[run])
            cumulative_local_tilt = tail_mean(local_tilt) - model_tilt_base
            if run in exp_tilt:
                exp_residual_change = float(exp_tilt[run]["resid_tilt_y_deg"]) - exp_tilt_base
            else:
                exp_residual_change = None

            exp_disp_peak_us1, exp_shear_peak_us1 = exp_peaks["US1"]
            exp_disp_peak_us2, exp_shear_peak_us2 = exp_peaks["US2"]
            all_rows.append({
                "method": method,
                "method_label": label,
                "run": run,
                "record": row_log["record"],
                "scale": float(row_log["scale"]),
                "target_A_period_s": float(row_log["target_A_period_s"]),
                "target_B_period_s": float(row_log["target_B_period_s"]),
                "carrier_B_requested_s": float(row_log["carrier_B_requested_s"]),
                "target_A_mm": float(row_log["target_A_mm"]),
                "target_B_mm": float(row_log["target_B_mm"]),
                "achieved_A_mm": float(row_log["achieved_A_mm"]),
                "achieved_B_mm": float(row_log["achieved_B_mm"]),
                "error_A_pct": float(row_log["error_A_pct"]),
                "error_B_pct": float(row_log["error_B_pct"]),
                "fit_converged": row_log["fit_converged"].lower() == "true",
                "peak_relative_OOP_mm_model": disp_peak,
                "peak_relative_OOP_mm_US1": exp_disp_peak_us1,
                "OOP_error_mm": None if exp_disp_peak_us1 is None else disp_peak - exp_disp_peak_us1,
                "OOP_error_pct": percent_error(disp_peak, exp_disp_peak_us1),
                "OOP_error_mm_US1": None if exp_disp_peak_us1 is None else disp_peak - exp_disp_peak_us1,
                "OOP_error_pct_US1": percent_error(disp_peak, exp_disp_peak_us1),
                "peak_relative_OOP_mm_US2": exp_disp_peak_us2,
                "OOP_error_mm_US2": None if exp_disp_peak_us2 is None else disp_peak - exp_disp_peak_us2,
                "OOP_error_pct_US2": percent_error(disp_peak, exp_disp_peak_us2),
                "peak_base_shear_kN_model": model_shear_peak,
                "peak_base_shear_kN_US1": exp_shear_peak_us1,
                "base_shear_error_kN": None if exp_shear_peak_us1 is None else model_shear_peak - exp_shear_peak_us1,
                "base_shear_error_pct": percent_error(model_shear_peak, exp_shear_peak_us1),
                "base_shear_error_kN_US1": None if exp_shear_peak_us1 is None else model_shear_peak - exp_shear_peak_us1,
                "base_shear_error_pct_US1": percent_error(model_shear_peak, exp_shear_peak_us1),
                "peak_base_shear_kN_US2": exp_shear_peak_us2,
                "base_shear_error_kN_US2": None if exp_shear_peak_us2 is None else model_shear_peak - exp_shear_peak_us2,
                "base_shear_error_pct_US2": percent_error(model_shear_peak, exp_shear_peak_us2),
                "cumulative_local_tilt_deg_model": cumulative_local_tilt,
                "cumulative_tilt_y_change_deg_US1": exp_residual_change,
                "local_tilt_difference_deg": None if exp_residual_change is None else cumulative_local_tilt - exp_residual_change,
            })
            model_tilt_curves.setdefault(method, {})[run] = cumulative_local_tilt
            if exp_residual_change is not None:
                exp_tilt_curve[run] = exp_residual_change

    return all_rows, model_tilt_curves, exp_tilt_curve


def write_metrics(rows):
    columns = list(rows[0].keys())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: round_or_blank(value) if isinstance(value, (float, np.floating)) else value
                             for key, value in row.items()})

    # A compact, specimen-specific scorecard. US-2 has four measured runs;
    # US-1 has three. Do not pool the observations because these are two walls.
    summary_rows = []
    for method, label in METHODS:
        subset = [r for r in rows if r["method"] == method and r["peak_relative_OOP_mm_US2"] is not None]
        summary_rows.append({
            "method": method,
            "method_label": label,
            "runs_compared_US2": len(subset),
            "mean_absolute_OOP_error_pct_US2": float(np.mean([abs(r["OOP_error_pct_US2"]) for r in subset])),
            "mean_absolute_base_shear_error_pct_US2": float(np.mean([abs(r["base_shear_error_pct_US2"]) for r in subset])),
        })
    with US2_SUMMARY_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({key: round_or_blank(value) if isinstance(value, (float, np.floating)) else value
                             for key, value in row.items()})


def make_comparison_figure(rows, focus_method=None, output_path=SUMMARY_PNG):
    colors = {
        "1_fixed_wall": "#276B8E",
        "2_fixed_record": "#E07A2D",
        "3_adaptive_wall": "#3E8C61",
        "4_adaptive_record": "#A14D8D",
    }
    if focus_method is not None and focus_method not in colors:
        raise ValueError("Unknown focus method: {}".format(focus_method))
    muted = "#AEB4B9"
    fade = 0.26
    fig, axs = plt.subplots(1, 3, figsize=(15.8, 5.0), dpi=170)
    ax = axs[0]
    ordered_methods = [entry for entry in METHODS if entry[0] != focus_method]
    if focus_method is not None:
        ordered_methods.append(next(entry for entry in METHODS if entry[0] == focus_method))
    else:
        ordered_methods = METHODS
    for method, label in ordered_methods:
        subset = [r for r in rows if r["method"] == method]
        highlighted = focus_method is None or method == focus_method
        ax.plot([r["run"] for r in subset], [r["peak_relative_OOP_mm_model"] for r in subset],
                marker="o", lw=1.8 if focus_method is None else (2.4 if focus_method == method else 1.6), ms=4,
                color=colors[method] if highlighted else muted,
                alpha=1.0 if highlighted else fade, zorder=3 if highlighted else 1,
                label=label)
    exp_rows_us1 = [r for r in rows if r["method"] == METHODS[0][0] and r["peak_relative_OOP_mm_US1"] is not None]
    exp_rows_us2 = [r for r in rows if r["method"] == METHODS[0][0] and r["peak_relative_OOP_mm_US2"] is not None]
    ax.plot([r["run"] for r in exp_rows_us1], [r["peak_relative_OOP_mm_US1"] for r in exp_rows_us1],
            marker="D", ls="--", color="black", lw=1.6, ms=5, zorder=4,
            label="US-1 experiment (Test 9)")
    ax.plot([r["run"] for r in exp_rows_us2], [r["peak_relative_OOP_mm_US2"] for r in exp_rows_us2],
            marker="^", ls=":", color="#8A5A44", lw=1.8, ms=5, zorder=4,
            label="US-2 experiment (Test 12)")
    ax.set_title("Peak relative OOP displacement", loc="left", fontweight="bold")
    ax.set_ylabel("Peak excursion from run start (mm)")
    ax.set_xticks(RUNS)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.3, frameon=False)

    ax = axs[1]
    for method, label in ordered_methods:
        subset = [r for r in rows if r["method"] == method]
        highlighted = focus_method is None or method == focus_method
        ax.plot([r["run"] for r in subset], [r["peak_base_shear_kN_model"] for r in subset],
                marker="o", lw=1.8 if focus_method is None else (2.4 if focus_method == method else 1.6), ms=4,
                color=colors[method] if highlighted else muted,
                alpha=1.0 if highlighted else fade, zorder=3 if highlighted else 1,
                label=label)
    exp_rows_us1 = [r for r in rows if r["method"] == METHODS[0][0] and r["peak_base_shear_kN_US1"] is not None]
    exp_rows_us2 = [r for r in rows if r["method"] == METHODS[0][0] and r["peak_base_shear_kN_US2"] is not None]
    ax.plot([r["run"] for r in exp_rows_us1], [r["peak_base_shear_kN_US1"] for r in exp_rows_us1],
            marker="D", ls="--", color="black", lw=1.6, ms=5, zorder=4,
            label="US-1 experiment (Test 9)")
    ax.plot([r["run"] for r in exp_rows_us2], [r["peak_base_shear_kN_US2"] for r in exp_rows_us2],
            marker="^", ls=":", color="#8A5A44", lw=1.8, ms=5, zorder=4,
            label="US-2 experiment (Test 12)")
    ax.set_title("Peak experiment-equivalent base shear", loc="left", fontweight="bold")
    ax.set_ylabel("Peak excursion from run start (kN)")
    ax.set_xlabel("FR76 sequence run")
    ax.set_xticks(RUNS)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.3, frameon=False)

    ax = axs[2]
    positions = np.arange(len(METHODS))
    width = 0.34
    for offset, key, title, color in (
        (-width / 2, "error_A_pct", "Wave A target (T1)", "#537EA5"),
        (width / 2, "error_B_pct", "Wave B target", "#CF8152"),
    ):
        for index, (method, _label) in enumerate(METHODS):
            vals = [abs(r[key]) for r in rows if r["method"] == method]
            highlighted = focus_method is None or method == focus_method
            if focus_method is None:
                legend_label = title if index == 0 else "_nolegend_"
            elif highlighted:
                legend_label = "{} ({})".format(title, dict(METHODS)[method])
            elif method == METHODS[0][0] and key == "error_A_pct":
                legend_label = "Methods 1, 3, 4 (both targets)"
            else:
                legend_label = "_nolegend_"
            ax.bar(positions[index] + offset, float(np.mean(vals)), width=width,
                   label=legend_label, color=color if highlighted else muted,
                   alpha=1.0 if highlighted else fade,
                   zorder=3 if highlighted else 1)
    ax.axhline(2.0, color="#666", lw=1, ls=":", label="2% solver tolerance")
    ax.set_xticks(positions)
    ax.set_xticklabels([label for _, label in METHODS], fontsize=8)
    if focus_method is not None:
        for index, tick in enumerate(ax.get_xticklabels()):
            method = METHODS[index][0]
            tick.set_color(colors[method] if method == focus_method else "#858B90")
            if method == focus_method:
                tick.set_fontweight("bold")
    ax.set_ylabel("Mean absolute spectral target error (%)")
    ax.set_title("Input signal calibration" + (" · Method 2 highlighted" if focus_method else ""),
                 loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=6.8 if focus_method else 7.2, frameon=False)

    suptitle = ("Route 2 four-method comparison from the common Strategy F run-21 state"
                if focus_method is None else
                "Route 2 comparison · Method 2 highlighted")
    fig.suptitle(suptitle,
                 x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.01,
             "US-1 Test 9 has no run-25 experiment; US-2 Test 12 spans runs 22–25. The specimens are scored separately. Model displacement is top-minus-table; base shear uses the supplied 50-Hz, mass-scaled inertia-force procedure.",
             fontsize=8, color="#555")
    fig.tight_layout(rect=(0, 0.045, 1, 0.93))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def make_tilt_figure(model_curves, exp_curve):
    colors = {
        "1_fixed_wall": "#276B8E",
        "2_fixed_record": "#E07A2D",
        "3_adaptive_wall": "#3E8C61",
        "4_adaptive_record": "#A14D8D",
    }
    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=170)
    for method, label in METHODS:
        runs = sorted(model_curves[method])
        ax.plot(runs, [model_curves[method][r] for r in runs], marker="o", lw=1.8,
                ms=4, color=colors[method], label=label)
    runs = sorted(exp_curve)
    ax.plot(runs, [exp_curve[r] for r in runs], marker="D", ls="--", lw=1.6,
            ms=5, color="black", label="US-1 tiltmeter settled change (from run 21)")
    ax.axhline(0, color="#555", lw=0.7)
    ax.set_title("Local tiltmeter residual after the FR76 runs", loc="left", fontweight="bold")
    ax.set_xlabel("FR76 sequence run")
    ax.set_ylabel("Cumulative change from run-21 baseline (degrees)")
    ax.set_xticks(RUNS)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, frameon=False)
    fig.text(0.01, 0.01,
             "Experiment uses quiet residuals from the inclinometer log; shaking peaks clip. Model uses tilt_local_incl with the assumed 2.33–2.53 m sensor chord, so this trend is provisional until mounting height is confirmed.",
             fontsize=8, color="#555", wrap=True)
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    fig.savefig(TILT_PNG, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    baseshear = load_baseshear_module()
    rows, model_tilt, exp_tilt = collect(baseshear)
    exp_us1_tilt = load_experimental_tilt("US1")
    exp_us2_tilt = load_experimental_tilt("US2")
    write_metrics(rows)
    make_comparison_figure(rows)
    make_comparison_figure(rows, focus_method="2_fixed_record", output_path=FOCUS_PNG)
    make_tilt_figure(model_tilt, exp_tilt)
    print("Metrics: {} ({} rows)".format(SUMMARY_CSV, len(rows)))
    print("US-2 scorecard: {}".format(US2_SUMMARY_CSV))
    print("US-1/US-2 tiltmeter residuals: {} ({} runs), {} ({} runs)".format(
        US1_TILT_CSV, len(exp_us1_tilt), US2_TILT_CSV, len(exp_us2_tilt)))
    print("Comparison figure: {}".format(SUMMARY_PNG))
    print("Method 2 focus figure: {}".format(FOCUS_PNG))
    print("Tiltmeter residual figure: {}".format(TILT_PNG))


if __name__ == "__main__":
    main()
