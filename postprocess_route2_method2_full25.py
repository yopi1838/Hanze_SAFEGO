# -*- coding: utf-8 -*-
"""Post-process the full 25-run Method 2 sequence.

Run from the Hanze_SAFEGO project root after route2_method_2_full25_3dec.py
has created its histories and route2_method_log.csv:

    python postprocess_route2_method2_full25.py

Optional paths can be supplied when the outputs live elsewhere:

    python postprocess_route2_method2_full25.py --method-dir route2_method2_full25/2_fixed_record

The script compares within-run peak displacement and cumulative residual
displacement for Method 2, the Strategy F record-on-model control, US-1, and
US-2. It also plots the fitted spectral-target errors recorded by the driver.
It deliberately does not compare base shear: the experimental mass convention
is still under review, and support reaction is not the same quantity as the
experiment-derived inertia force.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROTOCOL = {
    1: ("HU12", 0.50), 2: ("HU12", 0.75), 3: ("EC40", 0.20),
    4: ("HU12", 1.00), 5: ("HU12", 1.25), 6: ("EC40", 0.30),
    7: ("HU12", 1.50), 8: ("EC40", 0.40), 9: ("HU12", 1.75),
    10: ("HU12", 2.00), 11: ("EC40", 0.50), 12: ("HU12", 2.25),
    13: ("HU12", 2.50), 14: ("HU12", 2.75), 15: ("HU12", 3.00),
    16: ("HU12", 3.50), 17: ("HU12", 4.00), 18: ("HU12", 4.50),
    19: ("HU12", 5.00), 20: ("HU12", 5.50), 21: ("HU12", 6.00),
    22: ("FR76", 1.00), 23: ("FR76", 1.50), 24: ("FR76", 1.75),
    25: ("FR76", 2.00),
}
TAIL_FRAC = 0.05
COLORS = {
    "Method 2": "#E07A2D",
    "Strategy F control": "#5B7F74",
    "US-1 (Test 9)": "#2B5F8A",
    "US-2 (Test 12)": "#F0A020",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-dir", default="route2_method2_full25/2_fixed_record",
                        help="Method 2 folder containing flat RunNN_* history CSVs")
    parser.add_argument("--control-dir", default="stratF_full_results_US1",
                        help="Strategy F run-folder root")
    parser.add_argument("--exp-dir", default="EXP_DATA/processed_globalzero",
                        help="Folder with processed Test9/Test12 workbooks")
    return parser.parse_args()


def read_history(path):
    data = np.genfromtxt(str(path), skip_header=2)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("Could not parse history: {}".format(path))
    data = data[np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])]
    if len(data) < 2:
        raise ValueError("History has fewer than two valid samples: {}".format(path))
    return data[:, 0], data[:, 1]


def history_path(root, run, key, flat):
    root = Path(root)
    pattern = "Run{:02d}_*{}*.csv".format(run, key)
    candidates = list(root.glob(pattern)) if flat else []
    if not candidates:
        run_dirs = sorted(p for p in root.glob("Run{:02d}_*".format(run)) if p.is_dir())
        for run_dir in run_dirs:
            candidates.extend(run_dir.glob("*{}*.csv".format(key)))
    return sorted(candidates)[0] if candidates else None


def read_model_runs(root, flat):
    result = {}
    for run in PROTOCOL:
        path = history_path(root, run, "rel_disp_top_exp_mm", flat)
        if path is None:
            path = history_path(root, run, "rel_disp_top_mm", flat)
        if path is None:
            continue
        time, displacement_mm = read_history(path)
        result[run] = {"time": time, "disp": displacement_mm, "path": path}
    return result


def workbook_path(exp_dir, test, run):
    exp_dir = Path(exp_dir)
    for name in ("{}Run{}_processed_globalzero.xlsx".format(test, run),
                 "{}Run{:02d}_processed_globalzero.xlsx".format(test, run)):
        path = exp_dir / name
        if path.is_file():
            return path
    return None


def read_experiment_runs(exp_dir, test):
    try:
        import openpyxl
    except ImportError:
        print("! openpyxl is unavailable; experimental workbook comparisons skipped")
        return {}

    result = {}
    run_limit = 24 if test == "Test9" else 25
    for run in range(1, run_limit + 1):
        path = workbook_path(exp_dir, test, run)
        if path is None:
            continue
        ws = openpyxl.load_workbook(str(path), read_only=True, data_only=True).worksheets[0]
        rows = ws.iter_rows(values_only=True)
        headers = [str(value).strip() if value is not None else "" for value in next(rows)]
        if "U_avg" not in headers:
            print("! {} has no U_avg column; skipped".format(path))
            continue
        col = headers.index("U_avg")
        time_col = headers.index("Time") if "Time" in headers else None
        times, values = [], []
        for row in rows:
            if row and len(row) > col and row[col] is not None:
                try:
                    value = float(row[col])
                except (TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    values.append(value * 1000.0)  # workbook U_avg is metres; plot in mm
                    if time_col is not None and len(row) > time_col and row[time_col] is not None:
                        try:
                            times.append(float(row[time_col]))
                        except (TypeError, ValueError):
                            times.append(float(len(times)))
        if len(values) >= 2:
            if len(times) != len(values):
                times = list(np.arange(len(values), dtype=float))
            result[run] = {"time": np.asarray(times, dtype=float),
                           "time_units": "s" if time_col is not None else "samples",
                           "disp": np.asarray(values, dtype=float), "path": path}
    return result


def summarize(runs):
    if not runs:
        return {}
    baseline_run = min(runs)
    sequence_zero = float(runs[baseline_run]["disp"][0])
    metrics = {}
    for run, item in sorted(runs.items()):
        disp = np.asarray(item["disp"], dtype=float)
        tail_n = max(1, int(math.ceil(TAIL_FRAC * len(disp))))
        tail = float(np.mean(disp[-tail_n:]))
        metrics[run] = {
            "peak_incremental_mm": float(np.max(np.abs(disp - disp[0]))),
            "residual_incremental_mm": tail - float(disp[0]),
            "cumulative_residual_mm": tail - sequence_zero,
            "baseline_run": baseline_run,
        }
    return metrics


def load_log(method_dir):
    path = Path(method_dir) / "route2_method_log.csv"
    if not path.is_file():
        return {}
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        return {int(row["run"]): row for row in csv.DictReader(stream) if row.get("run")}


def write_metrics(path, datasets, log):
    fields = ["run", "record", "scale", "target_A_mm", "achieved_A_mm", "error_A_pct",
              "target_B_mm", "achieved_B_mm", "error_B_pct"]
    for name in datasets:
        key = name.lower().replace(" ", "_").replace("(", "").replace(")", "")
        fields.extend([key + "_peak_incremental_mm", key + "_residual_incremental_mm",
                       key + "_cumulative_residual_mm"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for run in PROTOCOL:
            record, scale = PROTOCOL[run]
            row = {"run": run, "record": record, "scale": scale}
            source = log.get(run, {})
            for key in fields:
                if key in source:
                    row[key] = source[key]
            for name, metrics in datasets.items():
                if run in metrics:
                    key = name.lower().replace(" ", "_").replace("(", "").replace(")", "")
                    for metric in ("peak_incremental_mm", "residual_incremental_mm", "cumulative_residual_mm"):
                        row[key + "_" + metric] = "{:.6f}".format(metrics[run][metric])
            writer.writerow(row)


def plot_sequence(path, datasets):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True, dpi=160)
    for name, metrics in datasets.items():
        runs = sorted(metrics)
        if not runs:
            continue
        color = COLORS[name]
        axes[0].plot(runs, [metrics[r]["peak_incremental_mm"] for r in runs],
                     marker="o", ms=3.5, lw=1.4, color=color, label=name)
        axes[1].plot(runs, [metrics[r]["cumulative_residual_mm"] for r in runs],
                     marker="o", ms=3.5, lw=1.4, color=color, label=name)
    for ax in axes:
        ax.set_xlim(1, 25)
        ax.set_xticks(range(1, 26))
        ax.grid(color="#E8E8E8", lw=0.7)
        ax.axvspan(21.5, 25.5, color="#F3E6D8", alpha=0.55, lw=0)
        ax.legend(frameon=False, ncol=2)
    axes[0].set_ylabel("Within-run peak displacement (mm)")
    axes[0].set_title("Peak OOP displacement by run", loc="left", fontweight="bold")
    axes[1].set_ylabel("Residual displacement from run 1 start (mm)")
    axes[1].set_xlabel("Run number (shaded: FR76 sequence)")
    axes[1].set_title("Cumulative residual displacement", loc="left", fontweight="bold")
    fig.suptitle("Method 2 full sequence vs Strategy F and experiments", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_fit(path, log):
    runs = sorted(run for run, row in log.items()
                  if row.get("error_A_pct") not in (None, "") and row.get("error_B_pct") not in (None, ""))
    if not runs:
        return False
    err_a = [float(log[run]["error_A_pct"]) for run in runs]
    err_b = [float(log[run]["error_B_pct"]) for run in runs]
    fig, ax = plt.subplots(figsize=(10, 4.4), dpi=160)
    ax.plot(runs, err_a, marker="o", ms=3.5, lw=1.3, color="#537EA5", label="Wave A target at T1")
    ax.plot(runs, err_b, marker="o", ms=3.5, lw=1.3, color="#CF8152", label="Wave B target at Teff")
    ax.axhspan(-2, 2, color="#6DAD7B", alpha=0.12, label="±2% fit tolerance")
    ax.axhline(0, color="#555555", lw=0.7)
    ax.set_xlim(1, 25)
    ax.set_xticks(range(1, 26))
    ax.set_xlabel("Run number")
    ax.set_ylabel("Spectral displacement error (%)")
    ax.set_title("Input fit check (this is not a structural-response error)", loc="left", fontweight="bold")
    ax.grid(color="#E8E8E8", lw=0.7)
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_histories(path, sources, runs=(22, 23, 24, 25)):
    available = [run for run in runs if any(run in data for data in sources.values())]
    if not available:
        return False
    fig, axes = plt.subplots(len(available), 1, figsize=(10, 2.3 * len(available)),
                             squeeze=False, dpi=150)
    for ax, run in zip(axes[:, 0], available):
        for name, data in sources.items():
            if run not in data:
                continue
            item = data[run]
            t = np.asarray(item["time"], dtype=float)
            u = np.asarray(item["disp"], dtype=float)
            ax.plot(t - t[0], u - u[0], lw=1.0, color=COLORS[name], label=name)
        record, scale = PROTOCOL[run]
        ax.set_title("Run {} · {} × {:.2f}".format(run, record, scale), loc="left", fontsize=9)
        ax.set_ylabel("Δ OOP disp. (mm)")
        ax.grid(color="#E8E8E8", lw=0.7)
        ax.legend(frameon=False, ncol=4, fontsize=8)
    axes[-1, 0].set_xlabel("Time within run (s; sample index only if a workbook lacks Time)")
    fig.suptitle("Late-run displacement histories: Method 2 vs Strategy F and tests", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    args = parse_args()
    method_dir = Path(args.method_dir)
    control_dir = Path(args.control_dir)
    exp_dir = Path(args.exp_dir)
    if not method_dir.is_dir():
        raise SystemExit("Method 2 output folder not found: {}".format(method_dir))

    log = load_log(method_dir)
    model = read_model_runs(method_dir, flat=True)
    if not model:
        raise SystemExit("No Method 2 displacement histories in {}. Run the 3DEC driver first.".format(method_dir))
    control = read_model_runs(control_dir, flat=False)
    us1 = read_experiment_runs(exp_dir, "Test9")
    us2 = read_experiment_runs(exp_dir, "Test12")

    source_runs = {
        "Method 2": model,
        "Strategy F control": control,
        "US-1 (Test 9)": us1,
        "US-2 (Test 12)": us2,
    }
    datasets = {name: summarize(runs) for name, runs in source_runs.items()}
    out_dir = method_dir / "postproc"
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = out_dir / "route2_method2_full25_metrics.csv"
    sequence_path = out_dir / "fig_method2_full25_sequence.png"
    fit_path = out_dir / "fig_method2_full25_input_fit.png"
    histories_path = out_dir / "fig_method2_full25_late_histories.png"

    write_metrics(metrics_path, datasets, log)
    plot_sequence(sequence_path, datasets)
    fit_made = plot_fit(fit_path, log)
    histories_made = plot_histories(histories_path, source_runs)

    print("Processed Method 2 runs: {} of 25".format(len(model)))
    print("Strategy F runs found: {}".format(len(control)))
    print("US-1 / US-2 runs found: {} / {}".format(len(us1), len(us2)))
    print("-> {}".format(metrics_path))
    print("-> {}".format(sequence_path))
    if fit_made:
        print("-> {}".format(fit_path))
    if histories_made:
        print("-> {}".format(histories_path))
    if len(model) < 25:
        print("! Some Method 2 runs are missing; metrics include available runs only.")


if __name__ == "__main__":
    main()
