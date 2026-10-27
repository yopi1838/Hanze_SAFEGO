# -*- coding: utf-8 -*-
"""Post-process Strategy F runs (the record on the model: strategy_F_full_3dec_LS.py and its variants) against the
tests US-1 (Test 9) and US-2 (Test 12), for one or several Strategy F output folders side by side.

Run from the Hanze_SAFEGO project root (needs postprocess_r2tm_full25.py next to it -- its readers are reused):

    python postprocess_stratF_full25.py                                   # every stratF_full_* folder with run folders
    python postprocess_stratF_full25.py --runs "f_t 0, 1.5 %=stratF_full_LS_MAXWELL15_v12" \\
                                           --runs "f_t 0, 1.5 %, _bc tables=stratF_full_LS_MAXWELL15_bc"
    python postprocess_stratF_full25.py --loop-runs 1,15,18,21,23,24

READ, per Strategy F folder: RunNN_*/ channel CSVs (missing runs are skipped) and strategy_C_log.csv (its '#' line
gives damping and start save, its rows the table file of each run).
  * wall displacement: rel_disp_top_mm (history 16: 2.06 m minus table), else 0.5(Ch3+Ch4)-Ch5;
  * force: F = sum m_j a_j with the SAME formula as the test side (--mass-set, default v3 = the group's formula of
    26 Oct 2026: wall 12.94 kN split by sensor height on vz_acc15/16/17, slabs 165 kg on Channel_13/14), 25 Hz
    low-pass zero-phase. "old" (408.52 / 320.66 / 404.13 kg) is only like-for-like with the old processed files;

  * table: Channel_12 (table velocity) -> achieved PGA, compared with scale x the PGA of the table file the run used;
    Channel_5 (table displacement) -> the net table displacement each run leaves, and its running sum. This is the
    drift check: vel_HU/EC/FR.txt integrate to 0.19 / 18.4 / 17.8 mm net at 100 %, the _bc tables to zero.
Tests: --exp-dir, default exp_processed_v3 (make_processed_v3.py: the group's formula, channels absolute, no factor).
  The old EXP_DATA/processed_globalzero files add the table acceleration to absolute channels (peaks about 35 % high);
  use them only with --mass-set old to reproduce earlier figures. Read as in postprocess_r2tm_full25.py (offset,
  US-1 run-18 slip and 25 Hz low-pass handled there; cached in OUT/_exp_cache.pkl).

OUTPUT (--out, default stratF_postproc/)
  fig_stratF_sequence.png   peak displacement per run, cumulative residual displacement, peak base shear per run
  fig_stratF_ratios.png     model / test ratios run by run (peak displacement, peak force) and model-vs-test scatter
  fig_stratF_loops.png      F-d loops at --loop-runs, each Strategy F folder over US-1 and US-2, pushover backbone
  fig_stratF_table.png      table input: achieved vs target PGA, net table displacement per run and cumulative
  fig_stratF_damping.png    equivalent viscous damping, sum E_D / (4 pi sum E_S) over the large cycles per run, and
                            every cycle against its amplitude, model vs tests, with the measured small-amplitude band
  stratF_metrics.csv
  stratF_damping.csv        per run and series: number of large cycles, energy-weighted xi_eq, median half-amplitude
"""
import argparse, csv, glob, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
try:
    import postprocess_r2tm_full25 as R          # readers, test processing, protocol, masses
except ImportError as e:
    raise SystemExit("postprocess_r2tm_full25.py must sit next to this script ({})".format(e))

US1, US2 = R.US1, R.US2
# Experiments in front, dark and thin, told apart by line style; numerical runs behind, wide and translucent.
TEST_STYLE = {US1: dict(color="#111827", ls="-", marker="o", mfc="white", label="US-1"),
              US2: dict(color="#1D4ED8", ls="--", marker="^", mfc="white", label="US-2")}
TEST_COLORS = {k: v["color"] for k, v in TEST_STYLE.items()}
VAR_COLORS = ["#0F766E", "#B91C1C", "#7C3AED", "#A16207", "#6B7280"]
NUM_ALPHA = 0.45
PUSH_COLOR = "#9CA3AF"


def short_label(label):
    """'LS_MAXWELL15_v12' -> 'NUM 1.5% v12', 'LS_MAXWELL_1p5' -> 'NUM 1.5%', 'LS_MAXWELL15_ft0p15_bc' -> 'NUM 1.5% ft0.15 bc'"""
    import re
    s = re.sub(r"LS_MAXWELL_?(\d)p(\d)", r"NUM \1.\2%", label)
    s = re.sub(r"LS_MAXWELL(\d)(\d)", r"NUM \1.\2%", s)
    s = s.replace("LS_MAXWELL", "NUM").replace("ft0p15", "ft0.15").replace("_", " ")
    return " ".join(s.split())
GRID = dict(color="#E8E8E8", lw=0.7)
DPI = 160
# the group's formula (26 Oct 2026; make_processed_v3.py), model side: wall 12.94 kN split by the sensor heights
# 0.66 / 1.26 / 2.06 m, slabs 165 kg on the model's slab channels. Defined here too, so this script also works next to
# an older postprocess_r2tm_full25.py that has no "v3" mass set (its model_run only knows the three wall channels).
MASS_V3 = {"vz_acc15": 490.80, "vz_acc16": 357.91, "vz_acc17": 470.34,
           "Channel_13_AccSlab1": 165.0, "Channel_14_AccSlab2": 165.0}
MASS_SETS = dict(R.MASS_SETS); MASS_SETS.setdefault("v3", MASS_V3)


def force_channels(folder, n, masses):
    """F = -sum m_j a_j (kN, + with + d) over any set of velocity channels; None if one is missing."""
    rf = R.run_folder(folder, n)
    if rf is None:
        return None, None, []
    vs = {k: R.find_hist(rf, k) for k in masses}
    miss = [k for k, v in vs.items() if v[0] is None]
    if miss:
        return None, None, miss
    k = min(len(v[1]) for v in vs.values()); tv = next(iter(vs.values()))[0][:k]
    F = -sum(m * R.accel_from_vel(tv, vs[c][1][:k])[0] for c, m in masses.items()) / 1e3
    F = F - np.median(F[: max(3, len(F) // 50)])
    return tv - tv[0], F, []


def args_():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--runs", action="append", default=[],
                   help='"label=folder" (repeat for several); default: every stratF_full_* folder with Run?? folders')
    p.add_argument("--exp-dir", default="exp_processed_v3")
    p.add_argument("--pushover-pos", default="pushover_results_v5_ft0_masonrypush_TopMason/pushover_pos.csv")
    p.add_argument("--pushover-neg", default="pushover_results_v5_ft0_masonrypush_TopMason/pushover_neg.csv")
    p.add_argument("--mass-set", default="v3", choices=sorted(MASS_SETS))
    p.add_argument("--loop-runs", default="1,15,18,21,23,24")
    p.add_argument("--table-dir", default=".", help="folder with the vel_*.txt tables named in strategy_C_log.csv")
    p.add_argument("--out", default="stratF_postproc")
    return p.parse_args()


# ============================ inputs =================================
def discover(specs):
    out = []
    for s in specs:
        lab, _, path = s.partition("=")
        if not path:
            path, lab = lab, os.path.basename(lab.rstrip("/\\"))
        out.append((lab.strip(), path.strip()))
    if not out:
        for d in sorted(glob.glob("stratF_full_*")):
            if os.path.isdir(d) and glob.glob(os.path.join(d, "Run[0-9][0-9]_*")):
                out.append((d.replace("stratF_full_", ""), d))
    return out


def read_log(folder):
    p = os.path.join(folder, "strategy_C_log.csv")
    info, tables = "", {}
    if not os.path.isfile(p):
        return info, tables
    lines = open(p).read().splitlines()
    head = [l for l in lines if l.startswith("#")]
    if head:
        info = head[0].lstrip("# ").strip()
    for r in csv.DictReader(l for l in lines if not l.startswith("#")):
        try:
            tables[int(float(r["run"]))] = r.get("table_file", "")
        except (KeyError, ValueError):
            continue
    return info, tables


_REC = {}
def record_props(path):
    """PGA (g) and net displacement (mm) of a table velocity file at scale 1"""
    if path in _REC:
        return _REC[path]
    out = None
    if os.path.isfile(path):
        d = np.loadtxt(path, skiprows=2); t, v = d[:, 0], d[:, 1]
        a = np.gradient(v, t)
        trap = getattr(np, "trapezoid", None) or np.trapz
        out = dict(pga_g=float(np.max(np.abs(a)) / 9.81), net_mm=float(trap(v, t) * 1e3))
    _REC[path] = out
    return out


def variant(label, folder, masses, table_dir):
    info, tables = read_log(folder)
    runs = {}
    for n, rec, sc in R.PROTOCOL:
        r = R.model_run(folder, n, (0.0, 0.0, 0.0) if isinstance(masses, dict) else masses, None, None)
        if r is None:
            continue
        if isinstance(masses, dict):                 # force computed here, whatever the R version
            r.pop("F", None); r.pop("tF", None)
            tF, F, miss = force_channels(folder, n, masses)
            if F is not None:
                r["tF"], r["F"] = tF, F
            elif n == 1 or miss:
                print("  ! {} run {}: no force, missing channels {}".format(label, n, miss))
        rf = R.run_folder(folder, n)
        tt, dt_ = R.find_hist(rf, "Channel_5_DispTable")
        if tt is not None:
            r["tab_start_mm"], r["tab_end_mm"] = float(dt_[0] * 1e3), float(np.mean(dt_[-max(3, len(dt_) // 50):]) * 1e3)
        tf = tables.get(n, "")
        rp = record_props(os.path.join(table_dir, tf)) if tf else None
        if rp:
            r["pga_target_g"] = sc * rp["pga_g"]; r["net_target_mm"] = sc * rp["net_mm"]
        r["table_file"] = tf
        runs[n] = r
    R.model_metrics(runs)
    if runs:
        z = runs[min(runs)].get("tab_start_mm")
        for r in runs.values():
            if "tab_end_mm" in r:
                r["tab_net_mm"] = r["tab_end_mm"] - r["tab_start_mm"]
                r["tab_cum_mm"] = r["tab_end_mm"] - z if z is not None else np.nan
    return dict(label=short_label(label), folder=folder, info=info, runs=runs)


# ============================ figures ================================
def finish(ax, xl=None, yl=None, title=None, runs_axis=False, legend=True):
    ax.grid(**GRID); ax.set_axisbelow(True)
    if runs_axis:
        ax.set_xlim(0.5, 25.5); ax.set_xticks(range(1, 26)); ax.axvspan(21.5, 25.5, color="#F3E6D8", alpha=0.55, lw=0)
    if xl: ax.set_xlabel(xl)
    if yl: ax.set_ylabel(yl)
    if title: ax.set_title(title, loc="left", fontweight="bold")
    if legend and ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, ncol=2, fontsize=8)


def series(runs, key):
    ns = sorted(n for n in runs if key in runs[n] and runs[n][key] is not None and np.isfinite(runs[n][key]))
    return ns, [runs[n][key] for n in ns]


def plot_sequence(path, variants, tests):
    fig, ax = plt.subplots(3, 1, figsize=(12, 11), dpi=DPI, sharex=True)
    for key, a, yl, ttl in (("peak_incr_mm", ax[0], "peak |d| within the run (mm)", "Peak displacement at 2.06 m, per run"),
                            ("cum_resid_mm", ax[1], "residual displacement (mm)", "Cumulative residual displacement (wall relative to table)"),
                            ("peak_F_kN", ax[2], "peak |F| (kN)", "Peak base shear sum m a, per run")):
        for k, V in enumerate(variants):
            ns, v = series(V["runs"], key)
            a.plot(ns, v, "-", marker="s", ms=6, mec="none", color=VAR_COLORS[k % len(VAR_COLORS)], lw=4.0,
                   alpha=NUM_ALPHA, solid_capstyle="round", zorder=2, label=V["label"])
        for name, runs in tests.items():
            ns, v = series(runs, key); st = TEST_STYLE[name]
            a.plot(ns, v, ls=st["ls"], marker=st["marker"], color=st["color"], mfc=st["mfc"], lw=1.3, ms=4.5,
                   zorder=5, label=st["label"])
        if key == "cum_resid_mm":
            a.axhline(0, color="#555555", lw=0.6)
        finish(a, yl=yl, title=ttl, runs_axis=True, legend=(a is ax[0]))
    ax[2].set_xlabel("run")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def plot_ratios(path, variants, tests):
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), dpi=DPI)
    nv = max(1, len(variants)); w = 0.8 / nv
    for row, (key, lab) in enumerate((("peak_incr_mm", "peak displacement"), ("peak_F_kN", "peak base shear"))):
        a = ax[row, 0]
        for k, V in enumerate(variants):
            ns, rr = [], []
            for n, r in sorted(V["runs"].items()):
                t = tests[US1].get(n)
                if t is not None and key in r and key in t and t[key] > 0:
                    ns.append(n); rr.append(r[key] / t[key])
            a.bar(np.array(ns) + (k - (nv - 1) / 2) * w, rr, width=w, color=VAR_COLORS[k % len(VAR_COLORS)],
                  alpha=0.75, label=V["label"])
        a.axhline(1.0, color="#333333", lw=0.9, zorder=4); a.axhspan(0.7, 1.3, color="#DDEEDD", alpha=0.5, lw=0, zorder=0)
        finish(a, yl="model / US-1", title="{}: model / US-1 per run (band +-30 %)".format(lab.capitalize()), runs_axis=True)
        a = ax[row, 1]; mx = 0
        for k, V in enumerate(variants):
            for name, mk in ((US1, "o"), (US2, "^")):
                xs, ys = [], []
                for n, r in V["runs"].items():
                    t = tests[name].get(n)
                    if t is not None and key in r and key in t:
                        xs.append(t[key]); ys.append(r[key])
                if xs:
                    a.plot(xs, ys, mk, color=VAR_COLORS[k % len(VAR_COLORS)], ms=5, mfc="none" if name == US2 else None,
                           alpha=0.8, label="{} / {}".format(V["label"], TEST_STYLE[name]["label"]))
                    mx = max(mx, max(xs), max(ys))
        if mx > 0:
            a.plot([0, mx * 1.05], [0, mx * 1.05], color="#333333", lw=0.9)
            a.set_xlim(0, mx * 1.05); a.set_ylim(0, mx * 1.05)
        u = "mm" if key == "peak_incr_mm" else "kN"
        finish(a, xl="test ({})".format(u), yl="Strategy F ({})".format(u), title="{}: model vs test".format(lab.capitalize()))
    ax[1, 0].set_xlabel("run")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def plot_loops(path, variants, tests, loop_runs, bb):
    runs = [n for n in loop_runs if any(n in V["runs"] for V in variants)]
    if not runs:
        return False
    nc = 3; nr = int(math.ceil(len(runs) / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(5.0 * nc, 4.3 * nr), dpi=DPI, squeeze=False)
    for k, n in enumerate(runs):
        a = ax[k // nc, k % nc]
        peaks = []
        for j, V in enumerate(variants):
            r = V["runs"].get(n)
            if r is None or "F" not in r:
                continue
            x = np.interp(r["tF"], r["t"], r["x"]) - r["x"][0]; F = r["F"]
            a.plot(x, F, lw=1.6, color=VAR_COLORS[j % len(VAR_COLORS)], alpha=NUM_ALPHA, zorder=2, label=V["label"])
            peaks.append((V["label"], np.max(np.abs(x)), np.max(np.abs(F)), VAR_COLORS[j % len(VAR_COLORS)]))
        for name in (US1, US2):
            r = tests[name].get(n)
            if r is None:
                continue
            i0, i1 = r["win"]; x = r["x"][i0:i1] - r["d_start"]; F = r["F"][i0:i1]; st = TEST_STYLE[name]
            a.plot(x, F, lw=1.0, ls=st["ls"], color=st["color"], zorder=5, label=st["label"])
            peaks.append((st["label"], np.max(np.abs(x)), np.max(np.abs(F)), st["color"]))
        lim_x = max([np.max(np.abs(l.get_xdata())) for l in a.lines] + [1.0]) * 1.15
        lim_y = max([np.max(np.abs(l.get_ydata())) for l in a.lines] + [1.0]) * 1.15
        for j, (d, F) in enumerate(bb):
            if len(d):
                a.plot(d, F, color=PUSH_COLOR, lw=1.0, ls=":", zorder=1, label="pushover" if j == 0 else None)
        # peak |d| / |F| as a compact table in the lower-right corner, instead of in the legend
        for i, (lab, dm, fm, c) in enumerate(reversed(peaks)):
            a.text(0.98, 0.03 + 0.065 * i, "{:.1f} mm  {:.1f} kN".format(dm, fm), transform=a.transAxes, ha="right",
                   va="bottom", fontsize=7, color=c)
        a.set_xlim(-lim_x, lim_x); a.set_ylim(-lim_y, lim_y)
        a.axhline(0, color="#555555", lw=0.5); a.axvline(0, color="#555555", lw=0.5)
        rec, sc = R.PROT[n]
        finish(a, xl="displacement from run start (mm)", yl="base shear sum m a (kN)",
               title="Run {} - {} x {:.2f}".format(n, rec, sc), legend=False)
        pass
    hl = {}                                  # one legend, with every series that appears in any panel
    for a_ in ax.flat:
        for h, l in zip(*a_.get_legend_handles_labels()):
            hl.setdefault(l, h)
    if hl:
        ax[0, 0].legend(list(hl.values()), list(hl.keys()), frameon=False, fontsize=8, loc="upper left")
    for k in range(len(runs), nr * nc):
        ax[k // nc, k % nc].axis("off")
    fig.suptitle("Hysteresis loops at selected runs: numerical (NUM, Strategy F) vs US-1 and US-2 "
                 "(corner: peak |d|, peak |F|, in the series colour)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    return True


def plot_table(path, variants):
    fig, ax = plt.subplots(3, 1, figsize=(12, 10), dpi=DPI, sharex=True)
    for k, V in enumerate(variants):
        c = VAR_COLORS[k % len(VAR_COLORS)]; runs = V["runs"]
        ns, v = series(runs, "pga_ach_g"); ax[0].plot(ns, v, "s-", color=c, lw=1.5, ms=4, alpha=0.8, label="achieved: " + V["label"])
        ns, v = series(runs, "pga_target_g")
        if ns:
            ax[0].plot(ns, v, "o", mfc="none", color=c, ms=6, label="scale x record PGA")
        ns, v = series(runs, "tab_net_mm"); ax[1].plot(ns, v, "s-", color=c, lw=1.5, ms=4, label=V["label"])
        ns, v = series(runs, "net_target_mm")
        if ns:
            ax[1].plot(ns, v, "o", mfc="none", color=c, ms=6, label="scale x record net displacement")
        ns, v = series(runs, "tab_cum_mm"); ax[2].plot(ns, v, "s-", color=c, lw=1.5, ms=4, label=V["label"])
    for a in ax[1:]:
        a.axhline(0, color="#555555", lw=0.6)
    finish(ax[0], yl="table PGA (g)", title="Table input: achieved PGA vs the scaled record", runs_axis=True)
    finish(ax[1], yl="net table displacement (mm)", title="Net table displacement left by each run (baseline drift)", runs_axis=True)
    finish(ax[2], xl="run", yl="cumulative table displacement (mm)", title="Cumulative table displacement since run 1", runs_axis=True)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


# ============================ equivalent damping =====================
XI_FRAC = 0.3                 # "large" cycles: half amplitude >= XI_FRAC x the run's largest
XI_MIN_AMP_MM = 0.2           # cycles smaller than this are dominated by sensor noise and filter phase (both sides);
                              # below it the loop method is not meaningful -- the ring-down (exp_damping_id.py) is
XI_MEASURED = (0.009, 0.017)  # small-amplitude damping of the tests (exp_damping_id.py: ring-down / half-power)


def cycles_xi(t, d, F, frac=0.0, a_min=XI_MIN_AMP_MM):
    """[(half amplitude mm, xi_eq)] for every displacement cycle: ED / (4 pi ES), ED = |loop area|,
    ES = 0.5 x (half p-p force) x (half p-p displacement). Cycles run between upward zero crossings of d about its
    1-s running mean, so a drifting or residual offset does not break the cycles."""
    from scipy.ndimage import uniform_filter1d
    t, d, F = np.asarray(t, float), np.asarray(d, float), np.asarray(F, float)
    if len(t) < 10:
        return []
    fs = 1.0 / float(np.median(np.diff(t)))
    dc = d - uniform_filter1d(d, max(3, int(fs * 1.0)))
    zc = np.where((dc[:-1] < 0) & (dc[1:] >= 0))[0]
    trap = getattr(np, "trapezoid", None) or np.trapz
    out = []
    for i0, i1 in zip(zc[:-1], zc[1:]):
        if i1 - i0 < 5:
            continue
        dd, FF = d[i0:i1 + 1], F[i0:i1 + 1]
        A = 0.5 * (dd.max() - dd.min()); Fa = 0.5 * (FF.max() - FF.min())
        if A < a_min or Fa <= 0:
            continue
        ED, ES = abs(trap(FF, dd)), 0.5 * Fa * A
        out.append((A, ED / (4 * np.pi * ES), ED, ES))
    if out and frac > 0:
        amax = max(A for A, _ in out); out = [c for c in out if c[0] >= frac * amax]
    return out


def _xi_summary(r, cyc_all):
    """per run: energy-weighted xi_eq = sum E_D / (4 pi sum E_S) over the large cycles (robust when a run has only one
    or two big cycles, where a median of single-cycle values is noise); the per-cycle values go to the scatter."""
    big = [c for c in cyc_all if c[0] >= XI_FRAC * max(c[0] for c in cyc_all)] if cyc_all else []
    r["xi_cycles"] = [(c[0], c[1]) for c in cyc_all]
    r["xi_med"] = float(sum(c[2] for c in big) / (4 * np.pi * sum(c[3] for c in big))) if big else np.nan
    r["xi_n"] = len(big)
    r["xi_amp_mm"] = float(np.median([c[0] for c in big])) if big else np.nan


def damping_metrics(variants, tests):
    for V in variants:
        for n, r in V["runs"].items():
            if "F" not in r:
                continue
            x = np.interp(r["tF"], r["t"], r["x"])
            _xi_summary(r, cycles_xi(r["tF"], x, r["F"]))
    for name, runs in tests.items():
        for n, r in runs.items():
            i0, i1 = r["win"]
            _xi_summary(r, cycles_xi(r["t"][i0:i1], r["x"][i0:i1], r["F"][i0:i1]))


def plot_damping(path, variants, tests):
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.6), dpi=DPI, gridspec_kw=dict(width_ratios=[1.35, 1]))
    a = ax[0]
    for k, V in enumerate(variants):
        ns, v = series(V["runs"], "xi_med")
        a.plot(ns, [100 * x for x in v], "-", marker="s", ms=6, mec="none", color=VAR_COLORS[k % len(VAR_COLORS)],
               lw=4.0, alpha=NUM_ALPHA, zorder=2, label=V["label"])
    for name, runs in tests.items():
        ns, v = series(runs, "xi_med"); st = TEST_STYLE[name]
        a.plot(ns, [100 * x for x in v], ls=st["ls"], marker=st["marker"], color=st["color"], mfc=st["mfc"], lw=1.3,
               ms=4.5, zorder=5, label=st["label"])
    a.axhspan(100 * XI_MEASURED[0], 100 * XI_MEASURED[1], color="#D1D5DB", alpha=0.6, lw=0, zorder=0,
              label="measured, small amplitude")
    finish(a, xl="run", yl="xi_eq, energy-weighted (%)",
           title="Equivalent viscous damping per run (cycles >= {:.0%} of the run's largest)".format(XI_FRAC), runs_axis=True)
    a.set_ylim(bottom=0)
    a.legend(frameon=False, ncol=2, fontsize=8, loc="upper left")
    b = ax[1]
    for k, V in enumerate(variants):
        A = [c[0] for r in V["runs"].values() for c in r.get("xi_cycles", [])]
        X = [100 * c[1] for r in V["runs"].values() for c in r.get("xi_cycles", [])]
        b.plot(A, X, "s", ms=3.5, mec="none", color=VAR_COLORS[k % len(VAR_COLORS)], alpha=0.35, zorder=2, label=V["label"])
    for name, runs in tests.items():
        st = TEST_STYLE[name]
        A = [c[0] for r in runs.values() for c in r.get("xi_cycles", [])]
        X = [100 * c[1] for r in runs.values() for c in r.get("xi_cycles", [])]
        b.plot(A, X, st["marker"], ms=4, mfc="white", mec=st["color"], mew=0.9, zorder=5, label=st["label"])
    b.axhspan(100 * XI_MEASURED[0], 100 * XI_MEASURED[1], color="#D1D5DB", alpha=0.6, lw=0, zorder=0)
    b.set_xscale("log"); b.set_ylim(0, 40)
    finish(b, xl="cycle half-amplitude at 2.06 m (mm)", yl="xi_eq (%)",
           title="Every cycle of every run: damping against amplitude", legend=False)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


# ============================ main ===================================
def main():
    a = args_()
    os.makedirs(a.out, exist_ok=True)
    masses = MASS_SETS[a.mass_set]
    print("force: mass set '{}': {}".format(a.mass_set, masses))
    specs = discover(a.runs)
    if not specs:
        raise SystemExit("no Strategy F folder found (give --runs \"label=folder\")")
    variants = []
    for lab, folder in specs:
        if not os.path.isdir(folder):
            print("! {}: {} not found".format(lab, folder)); continue
        V = variant(lab, folder, masses, a.table_dir)
        print("{}: {} runs from {}{}".format(lab, len(V["runs"]), folder, "  [{}]".format(V["info"]) if V["info"] else ""))
        bad = [n for n, r in V["runs"].items() if r.get("expkey_dev_mm", 0) > R.EXPKEY_TOL_MM]
        if bad:
            print("  note: rel_disp_top_exp_mm (history 22) differs from rel_disp_top_mm on runs {} -- not used".format(bad))
        variants.append(V)
    if not variants:
        raise SystemExit("no Strategy F runs read")
    print("reading the test workbooks ...")
    tests = R.read_tests(a.exp_dir, os.path.join(a.out, "_exp_cache.pkl"))
    for name in tests:
        R.test_metrics(tests[name])
    print("US-1 runs {}, US-2 runs {}".format(len(tests[US1]), len(tests[US2])))
    damping_metrics(variants, tests)
    R.PUSHOVER_FORCE = "applied"
    bb = [R.backbone(p) if os.path.isfile(p) else (np.array([]), np.array([])) for p in (a.pushover_pos, a.pushover_neg)]

    # ---- metrics CSV
    keys = ["peak_incr_mm", "peak_pos_mm", "peak_neg_mm", "resid_incr_mm", "cum_resid_mm", "peak_F_kN",
            "pga_ach_g", "pga_target_g", "tab_net_mm", "net_target_mm", "tab_cum_mm"]
    with open(os.path.join(a.out, "stratF_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "run", "record", "scale", "table_file"] + keys +
                   ["US1_peak_mm", "US1_peak_F_kN", "US2_peak_mm", "US2_peak_F_kN", "ratio_d_US1", "ratio_F_US1"])
        for V in variants:
            for n, r in sorted(V["runs"].items()):
                rec, sc = R.PROT[n]; t1, t2 = tests[US1].get(n, {}), tests[US2].get(n, {})
                row = [V["label"], n, rec, sc, r.get("table_file", "")] + [R.fmt(r.get(k)) for k in keys]
                row += [R.fmt(t1.get("peak_incr_mm")), R.fmt(t1.get("peak_F_kN")), R.fmt(t2.get("peak_incr_mm")), R.fmt(t2.get("peak_F_kN"))]
                row += [R.fmt(r["peak_incr_mm"] / t1["peak_incr_mm"]) if t1.get("peak_incr_mm") else "",
                        R.fmt(r.get("peak_F_kN", np.nan) / t1["peak_F_kN"]) if t1.get("peak_F_kN") else ""]
                w.writerow(row)
    print("-> " + os.path.join(a.out, "stratF_metrics.csv"))
    with open(os.path.join(a.out, "stratF_damping.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "run", "record", "scale", "n_large_cycles", "xi_eq_energy_weighted", "half_amp_median_mm",
                    "n_all_cycles"])
        for lab, runs in [(V["label"], V["runs"]) for V in variants] + [(TEST_STYLE[k]["label"], v) for k, v in tests.items()]:
            for n, r in sorted(runs.items()):
                if "xi_cycles" in r and n in R.PROT:
                    rec, sc = R.PROT[n]
                    w.writerow([lab, n, rec, sc, r["xi_n"], R.fmt(r["xi_med"]), R.fmt(r["xi_amp_mm"]), len(r["xi_cycles"])])
    print("-> " + os.path.join(a.out, "stratF_damping.csv"))

    # ---- printed summary
    print("\n{:24s} {:>3s} {:>8s} {:>8s} {:>8s} {:>8s} {:>9s} {:>9s} {:>8s} {:>8s} {:>8s}".format(
        "variant", "run", "d mod", "d US-1", "F mod", "F US-1", "resid cum", "tab cum", "xi mod", "xi US-1", "xi US-2"))
    for V in variants:
        for n, r in sorted(V["runs"].items()):
            t1 = tests[US1].get(n, {}); t2 = tests[US2].get(n, {})
            print("{:24s} {:3d} {:8.2f} {:>8s} {:>8s} {:>8s} {:9.2f} {:>9s} {:>8s} {:>8s} {:>8s}".format(
                V["label"][:24], n, r["peak_incr_mm"], R.fmt(t1.get("peak_incr_mm"), "{:.2f}"),
                R.fmt(r.get("peak_F_kN"), "{:.2f}"), R.fmt(t1.get("peak_F_kN"), "{:.2f}"), r["cum_resid_mm"],
                R.fmt(r.get("tab_cum_mm"), "{:.2f}"), R.fmt(r.get("xi_med"), "{:.3f}"),
                R.fmt(t1.get("xi_med"), "{:.3f}"), R.fmt(t2.get("xi_med"), "{:.3f}")))

    # ---- figures
    loop_runs = [int(s) for s in a.loop_runs.split(",") if s.strip()]
    for fn, func in (("fig_stratF_sequence.png", lambda p: plot_sequence(p, variants, tests)),
                     ("fig_stratF_ratios.png", lambda p: plot_ratios(p, variants, tests)),
                     ("fig_stratF_loops.png", lambda p: plot_loops(p, variants, tests, loop_runs, bb)),
                     ("fig_stratF_table.png", lambda p: plot_table(p, variants)),
                     ("fig_stratF_damping.png", lambda p: plot_damping(p, variants, tests))):
        p = os.path.join(a.out, fn)
        if func(p) is not False:
            print("-> " + p)


if __name__ == "__main__":
    main()
