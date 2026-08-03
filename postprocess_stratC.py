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
  <SIM_DIR>/postproc/fig6_top_quarter_oop.png    top-quarter peak OOP disp vs run

ONE SIMULATION vs THE EXPERIMENTS
    This script compares a single results folder against the shake-table
    specimens US-1 (Test 9) and US-2 (Test 12). The former sim-vs-sim
    --compare / --compare-labels options have been removed; to put two
    simulations side by side, post-process each and overlay the resulting
    postproc_summary.csv files.

    Experimental series are overlaid on fig1-fig4 and come from the derived
    CSVs, resolved via safego_paths.exp_derived() (this folder's postproc/
    first, then the canonical copies):

      exp_Test9_tilt.csv  / exp_Test12_tilt.csv        fig1, fig4
      exp_Test9_metrics.csv                            fig2   (US-1 only --
                                                       there is no US-2 set)
      exp_Test9_period_psd.csv / exp_Test12_..._psd    fig3
      xval_ch19_peaks.csv (US1_peak_mm/US2_peak_mm)    fig6
        -- despite the file name those two columns are the EXPERIMENTAL
           top-quarter-right sensor, so they pair with simulated Channel_4

    Any missing file is skipped with a note; the figure still renders.

Runs and their record/scale are auto-discovered from folder names
(RunNN_REC_sXpYY); nothing is hardcoded.

Metric definitions (unchanged from the previous postproc):
  peak      = max |x(t) - x(t0)|            within the run
  residual  = mean of last 5% of x(t), minus x at start of Run 1

  The experimental resid_* columns are likewise referenced to each
  specimen's Run-1 baseline (exp_tilt_from_raw.py sets the baselines from
  Run 1), so residuals are compared like for like. Peaks are per-run in both.

  fig3 normalisation differs by necessity: the simulation reports
  T_end/T_init against the fixed T1_init = 0.092 s, while each specimen is
  normalised by its own Run-1 identified period (US-1 0.09153 s,
  US-2 0.08818 s). Read fig3 as relative period growth, not absolute period.

Usage:
    python postprocess_stratC.py [SIM_DIR] [--label LBL]
Default SIM_DIR: stratC_results_NODAMP_v6_NEW
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, re, sys, argparse
from pathlib import Path

try:
    import safego_paths as sp
except ImportError:          # allow running from a copy outside the repo
    sp = None

RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")
REC_COLOR = {"HU12": "tab:blue", "EC40": "tab:orange", "FR76": "tab:red"}
KEY_TILT = "tilt_full_wall"
KEY_DISP = "rel_disp_top_mm"
CH4_NAME = "Channel_4_DispTopQRight"       # genuinely at x = 1.186 (quarter)
CH19_NAME = "Channel_19_DispTopQRight"     # NB: on the wall centreline
CH_TABLE = "Channel_5_DispTable"
CH19_XVAL = "xval_ch19_peaks.csv"
TAIL_FRAC = 0.05  # last 5% of samples -> residual

# Experimental specimens. Colours follow the project figure convention:
# red squares = simulation, blue circles = US-1, orange circles = US-2.
EXP_SPECS = [
    {"label": "US-1 (Test 9)", "color": "royalblue", "marker": "o",
     "tilt": "exp_Test9_tilt.csv", "metrics": "exp_Test9_metrics.csv",
     "period": "exp_Test9_period_psd.csv", "ls": "--",
     "ch19_col": "US1_peak_mm"},
    {"label": "US-2 (Test 12)", "color": "darkorange", "marker": "^",
     "tilt": "exp_Test12_tilt.csv", "metrics": None,
     "period": "exp_Test12_period_psd.csv", "ls": ":",
     "ch19_col": "US2_peak_mm"},
]

# sim tilt channel -> (experimental peak column, experimental residual column)
TILT_SEG_MAP = {
    "tilt_bot_seg":   ("peak_tilt_bot_seg",   "resid_tilt_bot_seg"),
    "tilt_low_seg":   ("peak_tilt_low_seg",   "resid_tilt_low_seg"),
    "tilt_up_seg":    ("peak_tilt_up_seg",    "resid_tilt_up_seg"),
    "tilt_full_wall": ("peak_tilt_full_wall", "resid_tilt_full_wall"),
}


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


def extract_rel_to_table(sim_dir, channel):
    """Peak table-referenced OOP displacement at `channel`, per run, in mm.

    Convention copied verbatim from ch19_xval.py so these numbers reconcile
    with xval_ch19_peaks.csv and figP6:

        rel(t) = ch(t) - Ch5(t)
        peak   = max |rel(t) - rel(t0)| * 1000

    Block histories are absolute displacements, while the experimental sensor
    is mounted on the shake-table frame and therefore already reads relative,
    so the simulation has to be referenced to the table before comparing.

    Which channel to pass, and why it matters:

      CH4_NAME  = Channel_4_DispTopQRight, at x = 1.186 m -- the actual
                  top-quarter-right position, so it is the geometrically
                  correct partner for the experimental 'Disp - Top quarter
                  right' sensor. This is what fig6 plots.
      CH19_NAME = Channel_19_DispTopQRight, which despite its name sits on
                  the wall CENTRELINE at x = 0.646 m (CLAUDE.md section 5).
                  ch19_xval.py and the published figP6 compare this against
                  the quarter-right sensor, i.e. across a position mismatch.
                  Kept as a secondary reference curve so the size of that
                  mismatch is visible rather than assumed.

    Note, verified numerically rather than assumed: the Channel_19 variant
    reproduces the FISH channel rel_disp_top_mm exactly (< 1e-3 mm on every
    run tested), because instrument_history_new.dat builds that channel from
    the very same two gridpoints. The Channel_4 variant is a genuinely
    independent series."""
    out = {}
    for run_no, record, scale, folder in discover_runs(sim_dir):
        fch = sorted(folder.glob("*" + channel + "*.csv"))
        f05 = sorted(folder.glob("*" + CH_TABLE + "*.csv"))
        if not (fch and f05):
            continue
        tch, cch = read_hist_csv(fch[0])
        t05, c05 = read_hist_csv(f05[0])
        if cch is None or c05 is None:
            continue
        n = min(len(cch), len(c05))
        rel = cch[:n] - c05[:n]
        out[run_no] = float(np.max(np.abs(rel - rel[0]))) * 1000.0
    return out


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


# ------------------------------------------------------------ experimental
def exp_csv(name, out_dir):
    """Locate a derived experimental CSV: this run's own postproc/ first,
    then whatever safego_paths.exp_derived() resolves (active sim, then the
    canonical copies under stratC_results_NODAMP_v6_NEW/postproc/).
    Returns None if it cannot be found anywhere."""
    local = Path(out_dir) / name
    if local.is_file():
        return local
    if sp is not None:
        try:
            p = Path(sp.exp_derived(name))
            if p.is_file():
                return p
        except Exception:
            pass
    return None


def read_exp_csv(path):
    """-> {column_name: {run_no: value}}, skipping blanks and non-numerics."""
    out = {}
    with open(str(path)) as f:
        for row in csv.DictReader(f):
            try:
                rn = int(float(row["run"]))
            except (TypeError, ValueError, KeyError):
                continue
            for k, v in row.items():
                if k == "run" or v in ("", None):
                    continue
                try:
                    out.setdefault(k, {})[rn] = float(v)
                except (TypeError, ValueError):
                    pass
    return out


def load_exp(out_dir):
    """Load every available experimental specimen. Missing files are skipped
    with a note rather than raising -- the figures degrade gracefully."""
    exps = []
    ch19_tab = None
    p = exp_csv(CH19_XVAL, out_dir)
    if p is not None:
        ch19_tab = read_exp_csv(p)
    else:
        print("  ! {} not found -- fig6 will show simulation only"
              .format(CH19_XVAL))
    for spec in EXP_SPECS:
        e = {"label": spec["label"], "color": spec["color"],
             "marker": spec["marker"], "ls": spec["ls"]}
        found = []
        for group in ("tilt", "metrics", "period"):
            name = spec[group]
            e[group] = None
            if not name:
                continue
            p = exp_csv(name, out_dir)
            if p is None:
                print("  ! experimental file not found, skipped: {}".format(name))
                continue
            e[group] = read_exp_csv(p)
            found.append(name)
        e["ch19"] = ch19_tab
        e["ch19_col"] = spec["ch19_col"]
        if ch19_tab and spec["ch19_col"] in ch19_tab:
            found.append("{}[{}]".format(CH19_XVAL, spec["ch19_col"]))
        if found:
            exps.append(e)
            print("  experiment '{}': {}".format(spec["label"], ", ".join(found)))
    if not exps:
        print("  ! no experimental data found -- figures will show simulation only")
    return exps


def exp_series(e, group, col):
    """(runs, values) for one experimental column; empty arrays if absent."""
    d = e.get(group) or {}
    d = d.get(col) or {}
    if not d:
        return np.array([]), np.array([])
    rn = sorted(d)
    return np.array(rn), np.array([d[r] for r in rn])


def exp_period_ratio(e):
    """Measured T1 normalised by that specimen's own Run-1 period, so it is
    comparable to the simulation's T_end/T_init."""
    rn, T = exp_series(e, "period", "T1_best_s")
    if not len(rn) or T[0] <= 0:
        return np.array([]), np.array([])
    return rn, T / T[0]


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
def fig_tilt(main, main_lbl, exps, out_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    pk_col, rs_col = TILT_SEG_MAP[KEY_TILT]
    for ax, metric, ttl, ecol in (
            (ax1, "peak", "Peak tilt (deg)", pk_col),
            (ax2, "residual", "Residual (permanent) tilt (deg)", rs_col)):
        x, y = chan_series(main, KEY_TILT, metric)
        ax.plot(x, y, "s-", color="red", ms=5, lw=1.2, label=main_lbl, zorder=3)
        for e in exps:
            xe, ye = exp_series(e, "tilt", ecol)
            if len(xe):
                ax.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                        ms=5, lw=1.0, mfc="none", alpha=0.9, label=e["label"],
                        zorder=2)
        ax.set_xlabel("Run"); ax.set_title(ttl)
        ax.grid(True, alpha=0.15); ax.legend(fontsize=9)
    fig.suptitle("Full-wall tilt: simulation vs experiment, peak and residual",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig1_tilt_compare.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig1_tilt_compare.png")


def fig_disp(main, main_lbl, exps, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    x, y = chan_series(main, KEY_DISP, "peak")
    ax.plot(x, y, "-", color="gray", lw=0.8, zorder=1)
    seen = set()
    for rn, v in zip(x, y):
        rec = main[rn]["record"]
        ax.plot(rn, v, "o", color=REC_COLOR.get(rec, "k"), ms=7,
                label=rec if rec not in seen else None, zorder=3)
        seen.add(rec)
    for e in exps:
        xe, ye = exp_series(e, "metrics", "peak_rel_mm")
        if len(xe):
            ax.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                    ms=5, lw=1.0, mfc="none", alpha=0.9, label=e["label"],
                    zorder=2)
    ax.set_xlabel("Run"); ax.set_ylabel("Peak relative top OOP disp (mm)")
    ax.set_title("Peak relative top OOP displacement: simulation vs experiment\n({})"
                 .format(main_lbl))
    ax.grid(True, alpha=0.15); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig2_disp.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig2_disp.png")


def fig_activation(main, exps, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    x1, y1 = run_series(main, "T_end_over_Tinit")
    if len(x1):
        ax.plot(x1, y1, "^-", color="tab:purple", ms=5, lw=1.2,
                label="NUM  T_end/T_init")
    for e in exps:
        xe, ye = exp_period_ratio(e)
        if len(xe):
            ax.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                    ms=5, lw=1.0, mfc="none", alpha=0.9,
                    label="{}  T1/T1(run 1)".format(e["label"]), zorder=2)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlabel("Run")
    ax.set_ylabel("Period ratio (see docstring)", color="tab:purple")
    ax.axhline(1.05, color="tab:purple", ls=":", lw=0.8, alpha=0.6)
    ax.axhline(1.20, color="tab:purple", ls=":", lw=0.8, alpha=0.6)
    ax2 = ax.twinx()
    x2, y2 = run_series(main, "beta_applied")
    if len(x2):
        ax2.plot(x2, y2, "s-", color="tab:green", ms=5, lw=1.2)
    ax2.set_ylabel("beta applied  (1=symmetric)", color="tab:green")
    ax.set_title("Period elongation vs experiment, and asymmetry activation",
                 fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.15)
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig3_activation.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig3_activation.png")


def fig_tilt_segments(main, exps, out_dir):
    """Colour encodes the segment, line style encodes the source, so the
    simulated and measured curves for one segment sit on the same colour.
    Only the segments the experiment actually measured get an overlay."""
    tilt_chs = [c for c in all_channels(main) if c.startswith("tilt")]
    if not tilt_chs:
        return
    cmap = plt.get_cmap("tab10")
    seg_color = {ch: cmap(i % 10) for i, ch in enumerate(tilt_chs)}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.8))
    for ax, metric, ttl, col_idx in ((ax1, "peak", "Peak tilt (deg)", 0),
                                     (ax2, "residual", "Residual tilt (deg)", 1)):
        for ch in tilt_chs:
            x, y = chan_series(main, ch, metric)
            ax.plot(x, y, "o-", ms=4, lw=1.1, color=seg_color[ch], zorder=3)
            ecol = TILT_SEG_MAP.get(ch, (None, None))[col_idx]
            if ecol is None:
                continue
            for e in exps:
                xe, ye = exp_series(e, "tilt", ecol)
                if len(xe):
                    ax.plot(xe, ye, marker=e["marker"], ls=e["ls"],
                            color=seg_color[ch], ms=4, lw=0.9, mfc="none",
                            alpha=0.75, zorder=2)
        ax.set_xlabel("Run"); ax.set_title(ttl)
        ax.grid(True, alpha=0.15)

    # two legends: one for the segment colours, one for the source styles
    seg_handles = [plt.Line2D([], [], color=seg_color[ch], lw=2, label=ch)
                   for ch in tilt_chs]
    src_handles = [plt.Line2D([], [], color="0.35", marker="o", ls="-", lw=1.1,
                              label="NUM (simulation)")]
    for e in exps:
        src_handles.append(plt.Line2D([], [], color="0.35", marker=e["marker"],
                                      ls=e["ls"], lw=0.9, mfc="none",
                                      label=e["label"]))
    leg1 = ax2.legend(handles=seg_handles, fontsize=8, loc="upper left",
                      title="segment", title_fontsize=8)
    ax2.add_artist(leg1)
    ax2.legend(handles=src_handles, fontsize=8, loc="lower right",
               title="source", title_fontsize=8)

    fig.suptitle("Tilt segments: simulation vs experiment "
                 "(colour = segment, style = source)",
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


def fig_top_quarter(ch4, ch19, main_lbl, exps, out_dir, ycap=None):
    """Peak table-referenced OOP displacement at the top-quarter-right
    position vs run.

    The simulated series is Channel_4 (x = 1.186 m), which is the same
    location as the experimental sensor. Channel_19 (wall centreline) is
    drawn thin and grey as a reference, because ch19_xval.py / figP6 use it
    and it is useful to see how much the position mismatch costs.

    Off-scale simulated points are clipped to the cap and annotated, the same
    treatment ch19_xval.py uses on figP6, so one runaway run cannot flatten
    the experimental curves."""
    if not ch4:
        print("  ! no Channel-4/Channel-5 pairs found, fig6 skipped")
        return
    rs = np.array(sorted(ch4))
    ys = np.array([ch4[r] for r in rs])

    exp_vals = []
    for e in exps:
        _, ye = exp_series(e, "ch19", e.get("ch19_col"))
        exp_vals.extend(list(ye))
    if ycap is None and exp_vals:
        exp_max = max(exp_vals)
        # only clip when the simulation genuinely dwarfs the experiment
        if ys.size and ys.max() > 3.0 * exp_max:
            ycap = 1.5 * max(exp_max, float(np.median(ys)))
    clipped = np.minimum(ys, ycap) if ycap else ys

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    for e in exps:
        xe, ye = exp_series(e, "ch19", e.get("ch19_col"))
        if len(xe):
            ax.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                    ms=6, lw=1.2, mfc="none", label=e["label"], zorder=2)
    if ch19:
        r19 = np.array(sorted(ch19))
        y19 = np.array([ch19[r] for r in r19])
        ax.plot(r19, np.minimum(y19, ycap) if ycap else y19, "-",
                color="0.6", lw=1.1, zorder=3,
                label="NUM Ch 19 (centreline, for reference)")
    ax.plot(rs, clipped, "s-", color="red", ms=6, lw=1.6,
            label=main_lbl + " -- Ch 4 (quarter right)", zorder=4)
    if ycap:
        off = [(r, v) for r, v in zip(rs, ys) if v > ycap]
        for k, (rr, v) in enumerate(off):
            ax.plot(rr, ycap, "s", ms=8, mfc="white", mec="red", mew=1.6,
                    zorder=5)
            ax.annotate("Run {}: {:.0f} mm".format(int(rr), v), xy=(rr, ycap),
                        xytext=(rr - 9.5, ycap * (0.92 - 0.09 * k)),
                        fontsize=9, color="red",
                        arrowprops=dict(arrowstyle="->", color="red", lw=1.0))
        ax.set_ylim(-0.02 * ycap, ycap * 1.05)
        print("     (fig6 y-axis capped at {:.1f} mm; {} run(s) off-scale)"
              .format(ycap, len(off)))
    ax.set_xlabel("Run")
    ax.set_ylabel("Peak OOP displacement, table-referenced (mm)")
    ax.set_title("Top quarter right: peak out-of-plane displacement, "
                 "simulation vs experiment")
    ax.grid(True, alpha=0.15)
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig6_top_quarter_oop.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig6_top_quarter_oop.png")


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
                 "peak_rel_disp_mm", "peak_ch4_rel_mm", "peak_ch19_rel_mm",
                 "T_end_over_Tinit",
                 "beta_applied", "Sd_target_mm"]
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
    print("\n{:>4} {:>5} {:>10} {:>11} {:>10} {:>9} {:>8} {:>6}".format(
        "Run", "Rec", "PeakTilt", "ResidTilt", "PeakDisp", "Ch4Disp",
        "T/Tinit", "beta"))
    for rn in sorted(main_data):
        r = main_data[rn]
        tilt = r["channels"].get(KEY_TILT, {})
        disp = r["channels"].get(KEY_DISP, {})
        print("{:>4} {:>5} {:>10.4f} {:>11.4f} {:>10.2f} {:>9.2f} {:>8.4f} {:>6.4f}"
              .format(rn, r["record"],
                      tilt.get("peak", float("nan")),
                      tilt.get("residual", float("nan")),
                      disp.get("peak", float("nan")),
                      float(r.get("peak_ch4_rel_mm", float("nan"))),
                      float(r.get("T_end_over_Tinit", float("nan"))),
                      float(r.get("beta_applied", float("nan")))))


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sim_dir", nargs="?",
                    default="stratC_results_GI_NORATCH")
    ap.add_argument("--ycap", type=float, default=None,
                    help="fig6 y-axis cap in mm (default: automatic, only "
                         "applied when the simulation dwarfs the experiment)")
    ap.add_argument("--label", default=None,
                    help="legend label for the simulation "
                         "(default: the results folder name)")
    args = ap.parse_args()

    sim_dir = Path(args.sim_dir)
    out_dir = sim_dir / "postproc"
    out_dir.mkdir(exist_ok=True)
    label = args.label if args.label else "NUM: {}".format(sim_dir.name)

    print("Extracting ALL channels: {}".format(sim_dir))
    main_data = extract_case(sim_dir, channels=None, verbose=True)
    print("  runs found: {}".format(len(main_data)))
    print("  channels found: {}".format(len(all_channels(main_data))))
    driver_cols = join_driver_summary(main_data, sim_dir)

    print("Extracting top-quarter OOP displacement (referenced to Channel 5)")
    ch4 = extract_rel_to_table(sim_dir, CH4_NAME)
    ch19 = extract_rel_to_table(sim_dir, CH19_NAME)
    print("  runs with a Ch4/Ch5 pair: {}   Ch19/Ch5: {}"
          .format(len(ch4), len(ch19)))
    for rn, v in ch4.items():
        if rn in main_data:
            main_data[rn]["peak_ch4_rel_mm"] = round(v, 4)
    for rn, v in ch19.items():
        if rn in main_data:
            main_data[rn]["peak_ch19_rel_mm"] = round(v, 4)

    print("Loading experimental data")
    exps = load_exp(out_dir)

    write_long_csv(main_data, out_dir)
    write_wide_csv(main_data, driver_cols, out_dir)
    console_table(main_data)

    fig_tilt(main_data, label, exps, out_dir)
    fig_disp(main_data, label, exps, out_dir)
    fig_activation(main_data, exps, out_dir)
    fig_tilt_segments(main_data, exps, out_dir)
    fig_top_quarter(ch4, ch19, label, exps, out_dir, ycap=args.ycap)
    fig_channel_grid(main_data, out_dir)
    print("\nDone. Output: {}".format(out_dir))


if __name__ == "__main__":
    main()
