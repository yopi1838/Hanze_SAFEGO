# -*- coding: ascii -*-
"""
Post-process a Strategy-C adaptive IDA results folder (runs OUTSIDE 3DEC).
==========================================================================
Full rewrite (2026-07). Scrapes EVERY exported instrument CSV in every
Run folder, computes per-channel peak/residual metrics, joins the driver
summary (strategy_C_summary.csv), and writes:

  <SIM_DIR>/postproc/postproc_all_channels.csv   long format, every channel
  <SIM_DIR>/postproc/postproc_summary.csv        wide format, key metrics
  <SIM_DIR>/postproc/fig1_tilt_compare.png       peak + residual tilt vs run
  <SIM_DIR>/postproc/fig2_disp.png               peak rel top OOP disp vs run
  <SIM_DIR>/postproc/fig3_activation.png         T_end/T_init and beta vs run
  <SIM_DIR>/postproc/fig4_tilt_segments.png      all tilt channels vs run
  <SIM_DIR>/postproc/fig5_channel_grid.png       small multiples, all channels

Runs and their record/scale are auto-discovered from folder names
(RunNN_REC_sXpYY); nothing is hardcoded. Comparison folders are always
re-extracted from raw CSVs (key channels only); no cached postproc files
are reused.

Metric definitions (unchanged from the previous postproc):
  peak      = max |x(t) - x(t0)|            within the run
  residual  = mean of last 5% of x(t), minus x at start of Run 1

Usage:
    python postprocess_stratC.py [SIM_DIR]
        [--label LBL] [--compare DIR ...] [--compare-labels LBL ...]
Defaults: SIM_DIR = stratC_results_NODAMP_v6_NEW,
          compare = stratC_results_DAMP1p5, stratC_results_NODAMP_v6
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, re, sys, argparse
from pathlib import Path

RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")
REC_COLOR = {"HU12": "tab:blue", "EC40": "tab:orange", "FR76": "tab:red"}
KEY_TILT = "tilt_full_wall"
KEY_DISP = "rel_disp_top_mm"
TAIL_FRAC = 0.05  # last 5% of samples -> residual


# ----------------------------------------------------------------- reading
def read_hist_csv(fpath):
    """3DEC history export: 2 header lines, whitespace-separated."""
    try:
        d = np.genfromtxt(str(fpath), skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2:
        return None, None
    m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    d = d[m]
    if len(d) < 5:
        return None, None
    return d[:, 0], d[:, 1]


def discover_runs(sim_dir):
    """[(run_no, record, scale, Path), ...] from RunNN_REC_sXpYY folders."""
    runs = []
    for d in sorted(Path(sim_dir).iterdir()):
        if not d.is_dir():
            continue
        m = RUN_RE.match(d.name)
        if m:
            runs.append((int(m.group(1)), m.group(2),
                         float(m.group(3)) + float(m.group(4)) / 100.0, d))
    return sorted(runs)


def channel_name(fname, folder_name):
    """Strip the 'RunNN_REC_sXpYY_' prefix and '.csv' suffix."""
    stem = Path(fname).name
    if stem.lower().endswith(".csv"):
        stem = stem[:-4]
    prefix = folder_name + "_"
    return stem[len(prefix):] if stem.startswith(prefix) else stem


# --------------------------------------------------------------- extraction
def extract_case(sim_dir, channels=None, verbose=False):
    """
    Scrape a results folder.
    channels=None -> every CSV; otherwise only channels whose name contains
    one of the given substrings (used for the comparison folders).
    Returns {run_no: {"record","scale","channels":{ch:{metrics}}}}.
    """
    out = {}
    run1_start = {}  # channel -> x(t0) of the first run that has it
    for run_no, record, scale, folder in discover_runs(sim_dir):
        entry = {"record": record, "scale": scale, "channels": {}}
        for f in sorted(folder.glob("*.csv")):
            ch = channel_name(f, folder.name)
            if channels is not None and not any(s in ch for s in channels):
                continue
            t, x = read_hist_csv(f)
            if x is None:
                if verbose:
                    print("    ! unreadable: {}".format(f.name))
                continue
            entry["channels"][ch] = compute_metrics(t, x, ch, run1_start)
        if entry["channels"]:
            out[run_no] = entry
    return out


def compute_metrics(t, x, ch, run1_start):
    x0 = float(x[0])
    if ch not in run1_start:
        run1_start[ch] = x0
    tail = x[int((1.0 - TAIL_FRAC) * len(x)):]
    end_mean = float(np.mean(tail))
    return {"n": int(len(x)),
            "dur_s": float(t[-1] - t[0]),
            "start": x0,
            "peak": float(np.max(np.abs(x - x0))),
            "vmin": float(np.min(x)),
            "vmax": float(np.max(x)),
            "end_mean": end_mean,
            "residual": end_mean - run1_start[ch]}


def join_driver_summary(data, sim_dir):
    """Merge every column of strategy_C_summary.csv into the run dicts."""
    sumf = Path(sim_dir) / "strategy_C_summary.csv"
    if not sumf.exists():
        print("  ! strategy_C_summary.csv not found, driver columns skipped")
        return []
    with open(str(sumf)) as f:
        rdr = csv.DictReader(f)
        cols = [c for c in rdr.fieldnames if c not in ("run", "record")]
        for row in rdr:
            rn = int(row["run"])
            if rn in data:
                for k in cols:
                    try:
                        data[rn][k] = float(row[k])
                    except (ValueError, TypeError):
                        data[rn][k] = row[k]
    return cols


# ------------------------------------------------------------------ helpers
def chan_series(data, ch, metric):
    rn = sorted(r for r in data if ch in data[r]["channels"])
    return (np.array(rn),
            np.array([data[r]["channels"][ch][metric] for r in rn]))


def run_series(data, key):
    rn = sorted(r for r in data if key in data[r] and data[r][key] != "")
    return np.array(rn), np.array([float(data[r][key]) for r in rn])


def all_channels(data):
    chs = set()
    for r in data.values():
        chs.update(r["channels"])
    return sorted(chs)


# ------------------------------------------------------------------ figures
def fig_tilt(main, main_lbl, cmps, out_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, metric, ttl in ((ax1, "peak", "Peak tilt (deg)"),
                            (ax2, "residual", "Residual (permanent) tilt (deg)")):
        x, y = chan_series(main, KEY_TILT, metric)
        ax.plot(x, y, "rs-", ms=5, lw=0.8, label=main_lbl, zorder=3)
        for lbl, d, style in cmps:
            xc, yc = chan_series(d, KEY_TILT, metric)
            if len(xc):
                ax.plot(xc, yc, style, ms=4, lw=0.8, alpha=0.7, label=lbl)
        ax.set_xlabel("Run"); ax.set_title(ttl)
        ax.grid(True, alpha=0.15); ax.legend(fontsize=9)
    fig.suptitle("Full-wall tilt: peak and residual per run",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig1_tilt_compare.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig1_tilt_compare.png")


def fig_disp(main, main_lbl, cmps, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    x, y = chan_series(main, KEY_DISP, "peak")
    ax.plot(x, y, "-", color="gray", lw=0.8, zorder=1)
    seen = set()
    for rn, v in zip(x, y):
        rec = main[rn]["record"]
        ax.plot(rn, v, "o", color=REC_COLOR.get(rec, "k"), ms=7,
                label=rec if rec not in seen else None, zorder=3)
        seen.add(rec)
    for lbl, d, style in cmps:
        xc, yc = chan_series(d, KEY_DISP, "peak")
        if len(xc):
            ax.plot(xc, yc, style, ms=4, lw=0.8, alpha=0.6, label=lbl)
    ax.set_xlabel("Run"); ax.set_ylabel("Peak relative top OOP disp (mm)")
    ax.set_title("Peak relative top OOP displacement per run ({})".format(main_lbl))
    ax.grid(True, alpha=0.15); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig2_disp.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig2_disp.png")


def fig_activation(main, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    x1, y1 = run_series(main, "T_end_over_Tinit")
    if len(x1):
        ax.plot(x1, y1, "^-", color="tab:purple", ms=5, lw=1.2)
    ax.set_xlabel("Run")
    ax.set_ylabel("T_end / T_init (purple)", color="tab:purple")
    ax.axhline(1.05, color="tab:purple", ls=":", lw=0.8, alpha=0.6)
    ax.axhline(1.20, color="tab:purple", ls=":", lw=0.8, alpha=0.6)
    ax2 = ax.twinx()
    x2, y2 = run_series(main, "beta_applied")
    if len(x2):
        ax2.plot(x2, y2, "s-", color="tab:green", ms=5, lw=1.2)
    ax2.set_ylabel("beta applied  (1=symmetric)", color="tab:green")
    ax.set_title("Period elongation and asymmetry activation per run",
                 fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.15)
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig3_activation.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig3_activation.png")


def fig_tilt_segments(main, out_dir):
    tilt_chs = [c for c in all_channels(main) if c.startswith("tilt")]
    if not tilt_chs:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, metric, ttl in ((ax1, "peak", "Peak tilt (deg)"),
                            (ax2, "residual", "Residual tilt (deg)")):
        for ch in tilt_chs:
            x, y = chan_series(main, ch, metric)
            ax.plot(x, y, "o-", ms=4, lw=0.9, label=ch)
        ax.set_xlabel("Run"); ax.set_title(ttl)
        ax.grid(True, alpha=0.15); ax.legend(fontsize=8)
    fig.suptitle("Tilt channels: segment-by-segment evolution",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig4_tilt_segments.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig4_tilt_segments.png")


def fig_channel_grid(main, out_dir):
    chs = all_channels(main)
    if not chs:
        return
    ncol = 4
    nrow = int(np.ceil(len(chs) / float(ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.6 * nrow),
                             sharex=True)
    axes = np.atleast_2d(axes)
    for i, ch in enumerate(chs):
        ax = axes[i // ncol, i % ncol]
        x, y = chan_series(main, ch, "peak")
        ax.plot(x, y, "o-", ms=3, lw=0.8, color="tab:red")
        ax.set_title(ch, fontsize=8)
        ax.grid(True, alpha=0.15)
        ax.tick_params(labelsize=7)
    for j in range(len(chs), nrow * ncol):
        axes[j // ncol, j % ncol].axis("off")
    fig.suptitle("Peak |x - x(t0)| per run, every exported channel",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(str(Path(out_dir) / "fig5_channel_grid.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig5_channel_grid.png")


# ------------------------------------------------------------------ output
def write_long_csv(main_data, out_dir):
    cols = ["run", "record", "scale", "channel", "n", "dur_s", "start",
            "peak", "vmin", "vmax", "end_mean", "residual"]
    with open(str(Path(out_dir) / "postproc_all_channels.csv"), "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for rn in sorted(main_data):
            r = main_data[rn]
            for ch in sorted(r["channels"]):
                c = r["channels"][ch]
                w.writerow([rn, r["record"], r["scale"], ch, c["n"],
                            round(c["dur_s"], 4), c["start"], c["peak"],
                            c["vmin"], c["vmax"], c["end_mean"],
                            c["residual"]])
    print("  -> postproc_all_channels.csv")


def write_wide_csv(main_data, driver_cols, out_dir):
    base_cols = ["run", "record", "peak_tilt_deg", "residual_tilt_deg",
                 "peak_rel_disp_mm", "T_end_over_Tinit", "beta_applied",
                 "Sd_target_mm"]
    extra_cols = [c for c in driver_cols if c not in base_cols]
    with open(str(Path(out_dir) / "postproc_summary.csv"), "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(base_cols + extra_cols)
        for rn in sorted(main_data):
            r = main_data[rn]
            tilt = r["channels"].get(KEY_TILT, {})
            disp = r["channels"].get(KEY_DISP, {})
            row = [rn, r["record"],
                   round(tilt["peak"], 4) if tilt else "",
                   round(tilt["residual"], 4) if tilt else "",
                   round(disp["peak"], 4) if disp else ""]
            for k in base_cols[5:] + extra_cols:
                v = r.get(k, "")
                row.append(round(v, 6) if isinstance(v, float) else v)
            w.writerow(row)
    print("  -> postproc_summary.csv")


def console_table(main_data):
    print("\n{:>4} {:>5} {:>10} {:>11} {:>10} {:>8} {:>6}".format(
        "Run", "Rec", "PeakTilt", "ResidTilt", "PeakDisp", "T/Tinit", "beta"))
    for rn in sorted(main_data):
        r = main_data[rn]
        tilt = r["channels"].get(KEY_TILT, {})
        disp = r["channels"].get(KEY_DISP, {})
        print("{:>4} {:>5} {:>10.4f} {:>11.4f} {:>10.2f} {:>8.4f} {:>6.4f}"
              .format(rn, r["record"],
                      tilt.get("peak", float("nan")),
                      tilt.get("residual", float("nan")),
                      disp.get("peak", float("nan")),
                      float(r.get("T_end_over_Tinit", float("nan"))),
                      float(r.get("beta_applied", float("nan")))))


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sim_dir", nargs="?",
                    default="stratC_results_GI_NORATCH")
    ap.add_argument("--label", default="No viscous damping (hysteretic only), asym. pulse")
    ap.add_argument("--compare", nargs="*",
                    default=["stratC_results_DAMP1p5",
                             "stratC_results_NODAMP_v6"])
    ap.add_argument("--compare-labels", nargs="*",
                    default=["1.5% Rayleigh", "NODAMP v6 (previous)"])
    args = ap.parse_args()

    sim_dir = Path(args.sim_dir)
    out_dir = sim_dir / "postproc"
    out_dir.mkdir(exist_ok=True)

    print("Extracting ALL channels: {}".format(sim_dir))
    main_data = extract_case(sim_dir, channels=None, verbose=True)
    print("  runs found: {}".format(len(main_data)))
    print("  channels found: {}".format(len(all_channels(main_data))))
    driver_cols = join_driver_summary(main_data, sim_dir)

    styles = ["o--", "d--", "v--"]
    cmps = []
    for i, cdir in enumerate(args.compare):
        if not Path(cdir).exists():
            print("  ! comparison folder missing, skipped: {}".format(cdir))
            continue
        lbl = (args.compare_labels[i] if i < len(args.compare_labels)
               else Path(cdir).name)
        print("Extracting comparison (key channels): {}".format(cdir))
        d = extract_case(cdir, channels=[KEY_TILT, KEY_DISP])
        if d:
            cmps.append((lbl, d, styles[i % len(styles)]))
            print("  comparison '{}': {} runs".format(lbl, len(d)))

    write_long_csv(main_data, out_dir)
    write_wide_csv(main_data, driver_cols, out_dir)
    console_table(main_data)

    fig_tilt(main_data, args.label, cmps, out_dir)
    fig_disp(main_data, args.label, cmps, out_dir)
    fig_activation(main_data, out_dir)
    fig_tilt_segments(main_data, out_dir)
    fig_channel_grid(main_data, out_dir)
    print("\nDone. Output: {}".format(out_dir))


if __name__ == "__main__":
    main()
