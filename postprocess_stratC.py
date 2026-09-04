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
  <SIM_DIR>/postproc/fig7_period_compare.png     period elongation vs exp (abs + log)

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
      US1_fig14_digitised.csv / US2_fig14_digitised.csv   fig3, fig7
        -- Fig 14 of Moshfeghi et al. (2024), digitised. PREFERRED.
           Falls back to exp_Test9/12_period_psd.csv with a printed warning.
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
  T_end/T_init against its own T1_init, while each specimen is normalised by
  its own Run-1 period -- US-1 0.09110 s and US-2 0.08747 s from the digitised
  Fig 14 (or 0.09153 / 0.08818 from the re-derived PSD series if that is what
  resolved). Read fig3 as relative period growth, not absolute period.

  A note on T1_init: it was 0.092 s historically, which is stale. Measured by
  the experiment's own method the as-built model reads 0.0948 s, so a results
  folder written with the old value carries a T_end/T_init inflated by ~3%.
  The ratio is only comparable between folders that used the same T1_init.

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
# Wall centreline at Z=2.06 referenced to the table (instrument_history_new.dat).
KEY_DISP = "rel_disp_top_mm"
# Experiment-matched EDP: mean(Ch3,Ch4) at Z=2.06 minus table, already in mm
# (FISH history 22, instrument_tilt_v2.dat). Preferred over KEY_DISP wherever
# it exists, because it is the definition used by Moshfeghi et al. Eq. (1).
KEY_DISP_EXP = "rel_disp_top_exp_mm"
CH3_NAME = "Channel_3_DispTopQLeft"        # x = 0.106 m  (quarter left)
CH4_NAME = "Channel_4_DispTopQRight"       # genuinely at x = 1.186 (quarter)
CH19_NAME = "Channel_19_DispTopQRight"     # NB: on the wall centreline
CH_TABLE = "Channel_5_DispTable"
CH19_XVAL = "xval_ch19_peaks.csv"
TAIL_FRAC = 0.05  # last 5% of samples -> residual

# Experimental specimens. Colours follow the project figure convention:
# red squares = simulation, blue circles = US-1, orange circles = US-2.
# PERIOD SOURCE, in order of preference.
#
#   1. US{1,2}_fig14_digitised.csv, column T_exp_s -- Fig 14 of Moshfeghi et
#      al. (2024) digitised from the figure's own gridlines. This is the
#      PUBLISHED experimental curve. Its calibration validates against Table 5
#      at both ends: US-1 run 1 = 0.09110 vs 0.091 and run 24 = 0.10749 vs
#      0.107 (+18.0% against the stated +17.5%); US-2 gives +13.7% against the
#      stated +13.9%.
#
#   2. exp_Test{9,12}_period_psd.csv, column T1_best_s -- re-derived from the
#      raw records. Kept as a fallback, but it is NOT equivalent: on US-1 it
#      holds runs 12-15 at a single frozen value of 0.10059, and its run-24
#      value of 0.11210 exceeds anything in the paper by 4.8% (nothing in
#      Fig 14 reaches 0.11). Mean |error| against Fig 14 is 1.51%.
#
# If the digitised file is missing the fallback is used and a note is printed,
# so a figure is never silently drawn against the weaker series.
EXP_SPECS = [
    {"label": "US-1 (Test 9)", "color": "royalblue", "marker": "o",
     "tilt": "exp_Test9_tilt.csv", "metrics": "exp_Test9_metrics.csv",
     "period_sources": [("US1_fig14_digitised.csv", "T_exp_s", "Fig 14 (published)"),
                        ("exp_Test9_period_psd.csv", "T1_best_s", "PSD re-derived")],
     "ls": "--", "ch19_col": "US1_peak_mm"},
    {"label": "US-2 (Test 12)", "color": "darkorange", "marker": "^",
     "tilt": "exp_Test12_tilt.csv", "metrics": None,
     "period_sources": [("US2_fig14_digitised.csv", "T_exp_s", "Fig 14 (published)"),
                        ("exp_Test12_period_psd.csv", "T1_best_s", "PSD re-derived")],
     "ls": ":", "ch19_col": "US2_peak_mm"},
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


def extract_topq_avg_to_table(sim_dir):
    """Primary EDP: the experiment's own definition, per run, in mm.

    Moshfeghi et al. (2024), Structures 66:106815, section 3.1, Eq. (1):

        U(t) = average( Disp-Top quarter left , Disp-Top quarter right )

    both at Z = 2.06 m (the 2nd-floor level, where the largest deformations
    were recorded), i.e. channels 3 and 4 of Table 3. The experimental
    potentiometers are carried on frame E, which is mounted on the shake
    table, so the recorded series are ALREADY table-referenced -- verified
    against the raw Test 9 data, where the bottom-quarter channel reads
    0.19 mm on run 11 while the table strokes 69.15 mm.

    3DEC block histories are absolute (they contain the prescribed base
    motion), so the simulation must be referenced to Channel 5 before it can
    be compared. The table stroke reaches 35-69 mm on the EC40 and FR76
    runs, so this is not a small correction.

        rel(t) = 0.5 * (Ch3(t) + Ch4(t)) - Ch5(t)
        peak   = max |rel(t) - rel(t0)| * 1000

    Reference values for US-1 (Test 9), from the raw data: peak of +29.48 /
    -25.52 mm at run 24, against Table 4's reported 25.0 / -30.0 mm (the
    paper's positive direction is the opposite sign convention).

    Two sources, in order of preference:

      1. The FISH channel `rel_disp_top_exp_mm` (history 22, defined in
         instrument_tilt_v2.dat), which already evaluates
         0.5*(Ch3 + Ch4) - Ch5 in mm at the gridpoint level. Preferred
         because it is computed inside 3DEC at every timestep rather than
         reconstructed from three separately-exported CSVs.
      2. Arithmetic on the Ch3 / Ch4 / Ch5 CSVs, for result folders written
         before instrument_tilt_v2.dat was in use.

    Returns {} if neither is available -- older folders (e.g. the NODAMP_v7
    set) export Channel_4 and Channel_5 but not Channel_3, so this EDP
    cannot be reconstructed for them and fig6 falls back to Ch4 alone."""
    out = {}
    for run_no, record, scale, folder in discover_runs(sim_dir):
        fexp = sorted(folder.glob("*" + KEY_DISP_EXP + "*.csv"))
        if fexp:
            texp, cexp = read_hist_csv(fexp[0])
            if cexp is not None and len(cexp):
                out[run_no] = float(np.max(np.abs(cexp - cexp[0])))
                continue
        f03 = sorted(folder.glob("*" + CH3_NAME + "*.csv"))
        f04 = sorted(folder.glob("*" + CH4_NAME + "*.csv"))
        f05 = sorted(folder.glob("*" + CH_TABLE + "*.csv"))
        if not (f03 and f04 and f05):
            continue
        t03, c03 = read_hist_csv(f03[0])
        t04, c04 = read_hist_csv(f04[0])
        t05, c05 = read_hist_csv(f05[0])
        if c03 is None or c04 is None or c05 is None:
            continue
        n = min(len(c03), len(c04), len(c05))
        rel = 0.5 * (c03[:n] + c04[:n]) - c05[:n]
        out[run_no] = float(np.max(np.abs(rel - rel[0]))) * 1000.0
    return out


def _is_missing(v):
    """Blank, None, or a non-finite float -- i.e. needs the log fallback."""
    if v == "" or v is None:
        return True
    try:
        return not np.isfinite(float(v))
    except (TypeError, ValueError):
        return False          # a non-numeric string (e.g. tid_source) counts as present


def join_driver_summary(data, sim_dir):
    """Merge every column of strategy_C_summary.csv into the run dicts, then
    backfill anything missing (absent run, blank, or NaN) from
    strategy_C_log.csv.

    WHY THE FALLBACK EXISTS: the summary CSV is rewritten at the end of a
    driver session from the in-memory summary list, so a resume that restarted
    that list truncates it (a stratF folder was seen with 12 summary rows
    against 25 log rows). The log is opened in append mode and flushed after
    every run, so it is the complete record. Backfilled T_end values are
    tagged T_end_src = "log" and drawn as OPEN markers in fig7, because for a
    run whose identification was rejected the log's T_end is the legacy
    single-channel estimator, not the three-channel PSD method -- a different
    measurement, kept visibly distinct rather than blended into the series."""
    cols = []
    sumf = Path(sim_dir) / "strategy_C_summary.csv"
    if sumf.exists():
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
                    if not _is_missing(row.get("T_end")):
                        data[rn]["T_end_src"] = "summary"
    else:
        print("  ! strategy_C_summary.csv not found, trying the log alone")

    logf = Path(sim_dir) / "strategy_C_log.csv"
    filled = []
    if logf.exists():
        with open(str(logf)) as f:
            for row in csv.DictReader(f):
                try:
                    rn = int(row["run"])
                except (TypeError, ValueError, KeyError):
                    continue
                if rn not in data:
                    continue
                for k, v in row.items():
                    if k in ("run", "record") or _is_missing(v):
                        continue
                    if k not in cols:
                        cols.append(k)
                    if _is_missing(data[rn].get(k, "")):
                        try:
                            data[rn][k] = float(v)
                        except (ValueError, TypeError):
                            data[rn][k] = v
                        if k == "T_end":
                            data[rn]["T_end_src"] = "log"
                            filled.append(rn)
    if filled:
        print("  T_end backfilled from strategy_C_log.csv for run(s): {}"
              .format(", ".join(str(r) for r in sorted(set(filled)))))
        print("    (open markers in fig7 -- the log value for a rejected "
              "identification is the legacy estimator)")
        if "T_end_src" not in cols:
            cols.append("T_end_src")
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
    # The digitised Fig 14 CSVs live with the analysis scripts, not inside a
    # results folder, so search there too before giving up.
    for cand in (Path.cwd() / name, Path(__file__).resolve().parent / name):
        if cand.is_file():
            return cand
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
        # period: try each source in order, keep the first that resolves
        e["period"] = None
        e["period_col"] = None
        e["period_src"] = None
        for fname, col, desc in spec["period_sources"]:
            pth = exp_csv(fname, out_dir)
            if pth is None:
                continue
            tab = read_exp_csv(pth)
            if col not in tab:
                print("  ! {} has no column '{}' -- skipped".format(fname, col))
                continue
            e["period"] = tab
            e["period_col"] = col
            e["period_src"] = desc
            found.append("{} [{}]".format(fname, desc))
            break
        if e["period"] is None:
            print("  ! no period source found for {} -- period figures will "
                  "show simulation only".format(spec["label"]))
        elif e["period_src"] != "Fig 14 (published)":
            print("  ! {} falling back to the RE-DERIVED period series ({}). "
                  "It holds runs 12-15 at one value and overshoots the tail; "
                  "generate US{}_fig14_digitised.csv with digitize_fig14.py "
                  "for the published curve."
                  .format(spec["label"], e["period_src"],
                          1 if "Test 9" in spec["label"] else 2))
        for group in ("tilt", "metrics"):
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


def exp_period(e):
    """(runs, T) from whichever period source resolved for this specimen."""
    col = e.get("period_col")
    if not col:
        return np.array([]), np.array([])
    return exp_series(e, "period", col)


def exp_period_ratio(e):
    """Measured T normalised by that specimen's own Run-1 period, so it is
    comparable to the simulation's T_end/T_init."""
    rn, T = exp_period(e)
    if not len(rn) or T[0] <= 0:
        return np.array([]), np.array([])
    return rn, T / T[0]


# ------------------------------------------------------------------ helpers
def chan_series(data, ch, metric):
    rn = sorted(r for r in data if ch in data[r]["channels"])
    return (np.array(rn),
            np.array([data[r]["channels"][ch][metric] for r in rn]))


def run_series(data, key):
    rn, vals = [], []
    for r in sorted(data):
        v = data[r].get(key, "")
        if v == "" or v is None:
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(fv):
            continue
        rn.append(r); vals.append(fv)
    return np.array(rn), np.array(vals)


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


def resolve_disp_key(main):
    """Prefer the experiment-matched EDP where the run exported it.

    KEY_DISP_EXP is mean(Ch3,Ch4) - Ch5 (Moshfeghi et al. Eq. 1);
    KEY_DISP is the wall-centreline equivalent, kept as the fallback for
    result folders written before instrument_tilt_v2.dat was in use."""
    for rn in main:
        if KEY_DISP_EXP in main[rn].get("channels", {}):
            return KEY_DISP_EXP
    return KEY_DISP


def fig_disp(main, main_lbl, exps, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    key = resolve_disp_key(main)
    if key != KEY_DISP:
        print("     (fig2 EDP = {}, experiment-matched)".format(key))
    x, y = chan_series(main, key, "peak")
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


def fig_period_compare(main, main_lbl, exps, out_dir):
    """Dedicated period-elongation comparison: the simulation ring-down period
    (T_end) against each specimen's small-amplitude PSD period (T1), shown both
    in absolute seconds (left) and normalised to each series' own run-1 value on
    a log axis (right).

    IMPORTANT caveat, annotated on the figure: the two are NOT the same
    measurement once the wall rocks. T_end is the LARGE-amplitude rocking
    ring-down period; the experimental T1 is the SMALL-amplitude PSD period of
    the re-closed cracked wall between runs. They are comparable only in the
    pre-rocking (elastic/cracking) regime -- the divergence at rocking onset is
    expected, not a model error."""
    xs, ys_abs = run_series(main, "T_end")
    _, ys_rat = run_series(main, "T_end_over_Tinit")
    if not len(xs):
        print("  ! no T_end column -- fig7 skipped")
        return

    # rocking-onset run: first run whose ring-down period exceeds 2x initial
    onset = None
    xr, yr = run_series(main, "T_end_over_Tinit")
    for xi, yi in zip(xr, yr):
        if yi >= 2.0:
            onset = xi
            break

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    bf = set(r for r in main
             if str(main[r].get("T_end_src", "")) == "log")

    # ---- panel A: absolute period ----
    ax1.plot(xs, ys_abs, "-o", color="tab:purple", lw=2, ms=6,
             label=r"Sim  ring-down $T_{end}$", zorder=3)
    if bf:
        m = np.isin(xs, sorted(bf))
        ax1.plot(xs[m], ys_abs[m], "o", ms=6, mfc="white", mec="tab:purple",
                 mew=1.6, zorder=4, label="backfilled from log (legacy est.)")
    for e in exps:
        xe, ye = exp_period(e)
        if len(xe):
            ax1.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                     ms=5, lw=1.2, mfc="none", alpha=0.9,
                     label="{}  $T_1$ [{}]".format(e["label"], e["period_src"]),
                     zorder=2)
    ax1.set_xlabel("Run"); ax1.set_ylabel("Identified period  T  (s)")
    ax1.set_title("Absolute period vs run")
    ax1.grid(True, alpha=0.2); ax1.legend(fontsize=9, loc="upper left")

    # ---- panel B: normalised elongation, log-y ----
    ax2.plot(xs, ys_rat, "-o", color="tab:purple", lw=2, ms=6,
             label=r"Sim  $T_{end}/T_{init}$", zorder=3)
    if bf:
        m2 = np.isin(xs, sorted(bf))
        ax2.plot(xs[m2], ys_rat[m2], "o", ms=6, mfc="white", mec="tab:purple",
                 mew=1.6, zorder=4)
    for e in exps:
        xe, ye = exp_period_ratio(e)
        if len(xe):
            ax2.plot(xe, ye, marker=e["marker"], ls=e["ls"], color=e["color"],
                     ms=5, lw=1.2, mfc="none", alpha=0.9,
                     label=r"{}  $T_1/T_1$(run 1)".format(e["label"]), zorder=2)
    ax2.set_yscale("log")
    # y-limits from the data, not a fixed 1-8 span: with elongations of only
    # 1.0-1.25 the old fixed axis flattened every series onto the baseline.
    ymax = float(np.nanmax(ys_rat)) if len(ys_rat) else 1.3
    for e in exps:
        _, ye_r = exp_period_ratio(e)
        if len(ye_r):
            ymax = max(ymax, float(np.nanmax(ye_r)))
    ytop = max(1.3, 1.08 * ymax)
    yticks = [t for t in (1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0)
              if t <= ytop * 1.001]
    ax2.set_yticks(yticks)
    ax2.set_yticklabels(["{:.2f}".format(t) for t in yticks])
    ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax2.set_ylim(0.98, 1.5)
    ax2.axhline(1.0, color="0.7", lw=0.8)
    ax2.set_xlabel("Run")
    ax2.set_ylabel("Period elongation factor  (log scale)")
    ax2.set_title("Period elongation, normalised to run 1")
    ax2.grid(True, alpha=0.2, which="both")
    ax2.legend(fontsize=9, loc="upper left")

    # rocking-regime shading + measurement-basis note
    if onset is not None:
        for ax in (ax1, ax2):
            ax.axvspan(onset - 0.5, float(max(xs)) + 0.5, color="0.9", zorder=0)
            ax.axvline(onset - 0.5, color="0.6", ls="--", lw=1)
        ax1.annotate("rocking regime: sim ring-down is the large-\n"
                     "amplitude rocking period, exp $T_1$ the small-\n"
                     "amplitude PSD period -- not directly comparable",
                     xy=(0.97, 0.03), xycoords="axes fraction",
                     ha="right", va="bottom", fontsize=8, color="0.35")

    fig.suptitle("Period elongation: simulation vs experiment\n({})".format(main_lbl),
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(str(Path(out_dir) / "fig7_period_compare.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig); print("  -> fig7_period_compare.png")


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


def fig_top_quarter(topq, ch4, ch19, main_lbl, exps, out_dir, ycap=None):
    """Peak table-referenced OOP displacement at Z = 2.06 m vs run.

    PRIMARY series is `topq` = 0.5*(Ch3 + Ch4) - Ch5, which is the
    experiment's own EDP definition (Moshfeghi et al. 2024, Eq. 1). This is
    the curve to compare against the measured data.

    Channel_4 alone (x = 1.186 m) and Channel_19 (wall centreline) are drawn
    thin for reference only: the former is a single-sensor version of the
    same quantity, the latter is what ch19_xval.py / figP6 historically used
    and is retained so the position mismatch stays visible.

    Off-scale simulated points are clipped to the cap and annotated, the same
    treatment ch19_xval.py uses on figP6, so one runaway run cannot flatten
    the experimental curves."""
    primary = topq if topq else ch4
    plabel = ("mean(Ch3,Ch4) at Z=2.06 m" if topq
              else "Ch 4 only at Z=2.06 m (Ch3 not exported)")
    if not primary:
        print("  ! no top-quarter/Channel-5 pairs found, fig6 skipped")
        return
    rs = np.array(sorted(primary))
    ys = np.array([primary[r] for r in rs])

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
    if topq and ch4:
        r4 = np.array(sorted(ch4))
        y4 = np.array([ch4[r] for r in r4])
        ax.plot(r4, np.minimum(y4, ycap) if ycap else y4, "-",
                color="salmon", lw=1.1, zorder=3,
                label="NUM Ch 4 alone (for reference)")
    ax.plot(rs, clipped, "s-", color="red", ms=6, lw=1.6,
            label=main_lbl + " -- " + plabel, zorder=4)
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
    ax.set_title("Top quarter (Z = 2.06 m): peak out-of-plane displacement, "
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
                 "peak_topq_avg_rel_mm",
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
    default="stratF_full_results_US1")
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

    print("Extracting primary EDP: mean(Ch3,Ch4) at Z=2.06 m, "
          "referenced to Channel 5")
    topq = extract_topq_avg_to_table(sim_dir)
    print("  runs with a Ch3/Ch4/Ch5 triple: {}".format(len(topq)))
    if not topq:
        print("  ! none found -- fig6 falls back to Ch4 alone. Check that "
              "Channel_3_DispTopQLeft is being exported.")
    for rn, v in topq.items():
        if rn in main_data:
            main_data[rn]["peak_topq_avg_rel_mm"] = round(v, 4)

    print("Loading experimental data")
    exps = load_exp(out_dir)

    write_long_csv(main_data, out_dir)
    write_wide_csv(main_data, driver_cols, out_dir)
    console_table(main_data)

    fig_tilt(main_data, label, exps, out_dir)
    fig_disp(main_data, label, exps, out_dir)
    fig_activation(main_data, exps, out_dir)
    fig_period_compare(main_data, label, exps, out_dir)
    fig_tilt_segments(main_data, exps, out_dir)
    fig_top_quarter(topq, ch4, ch19, label, exps, out_dir, ycap=args.ycap)
    fig_channel_grid(main_data, out_dir)
    print("\nDone. Output: {}".format(out_dir))


if __name__ == "__main__":
    main()