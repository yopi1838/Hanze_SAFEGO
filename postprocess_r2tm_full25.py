# -*- coding: utf-8 -*-
"""Post-process the full 25-run Route 2 sequence in Ihsan's form (route2_Tm_ihsan_full25_3dec.py: wave A at the
record's T_p, wave B at its T_m, amplitudes on Sd(T1) and Sd(T_eff)).

Run from the Hanze_SAFEGO project root once the 3DEC driver has written its run folders (missing runs are skipped):

    python postprocess_r2tm_full25.py
    python postprocess_r2tm_full25.py --run-dir r2tm_TpTm_full25_LS_maxwell --control-dir stratF_full_LS_MAXWELL15_v12

Compared, run by run, with the Strategy F record-on-model control (default: stratF_full_LS_MAXWELL15_v12, f_t = 0,
1.5 % Maxwell) and the tests US-1 (Test 9) and US-2 (Test 12), read from exp_processed_v3 (make_processed_v3.py):
  * within-run peak displacement and cumulative residual displacement at the 2.06 m sensor, relative to the table;
  * peak base shear F = sum m_j a_j with the group's formula (26 Oct 2026) on BOTH sides (--mass-set v3, default):
    wall 12.94 kN split by the sensor heights 0.66 / 1.26 / 2.06 m, slabs 165 kg on the slab accelerometers (model:
    Channel_13_AccSlab1 / Channel_14_AccSlab2), no top beam. "old" (408.52 / 320.66 / 404.13 kg on the wall sensors)
    only reproduces the earlier figures together with --exp-dir EXP_DATA/processed_globalzero, whose files add the
    table acceleration to absolute channels (peaks about 35 % high);
  * model / test ratios of both peaks, and model vs test scatter;
  * F-d hysteresis loops of selected runs, model over US-1 and US-2, each zeroed at its own run start.
Plus an input check (the table motion 3DEC produced against the designed signal and the two spectral targets) and
table-velocity / displacement histories of the selected runs.

MODEL DISPLACEMENT: rel_disp_top_mm (FISH history 16: dz(0.646, 2.06) - dz(table), mm), else 0.5(Ch3+Ch4)-Ch5.
  rel_disp_top_exp_mm is exported by NUMBER (history 22) and is only the tilt_v2 channel if instrument_tilt_v2.dat
  was called before instrument_baseshear_exp.dat; it is cross-checked here, never used.
MODEL FORCE: from the exported absolute velocities (vz_acc15/16/17, else Channel_15/16/17, plus Channel_13/14 for the
  slabs), differentiated, 25 Hz low-pass zero-phase.
FIGURES: experiments in front (US-1 black solid, US-2 blue dashed), the numerical runs (NUM) behind as wide translucent
  lines; one legend per figure; on the loops the peak |d| / |F| of each series sit in the corner, in its colour.
TEST: U_avg (m, global zero) and base_shear_kN / mass_scale with sign -1 (restoring force positive with positive d),
  25 Hz low-pass on both, d-F channel slip corrected by cross-correlation (US-1 run 18), force offset removed by the
  pre-event median -- the processing of exp_envelope_peaks.py. Parsed workbooks are cached in postproc/_exp_cache.pkl.
"""
import argparse, csv, glob, math, os, pickle, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROTOCOL = [(1, "HU12", 0.50), (2, "HU12", 0.75), (3, "EC40", 0.20), (4, "HU12", 1.00), (5, "HU12", 1.25),
            (6, "EC40", 0.30), (7, "HU12", 1.50), (8, "EC40", 0.40), (9, "HU12", 1.75), (10, "HU12", 2.00),
            (11, "EC40", 0.50), (12, "HU12", 2.25), (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00),
            (16, "HU12", 3.50), (17, "HU12", 4.00), (18, "HU12", 4.50), (19, "HU12", 5.00), (20, "HU12", 5.50),
            (21, "HU12", 6.00), (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75), (25, "FR76", 2.00)]
PROT = {n: (r, s) for n, r, s in PROTOCOL}
MASS_SETS = {"old": (408.52, 320.66, 404.13), "v2": (455.0, 670.0, 606.0),
             # v3 = the group's formula (26 Oct 2026; make_processed_v3.py): wall 12.94 kN split by the sensor heights
             # 0.66 / 1.26 / 2.06 m (US-1), slabs 165 kg on their own accelerometers (ch13 / ch14), no top beam.
             # The model exports the slab channels as Channel_13_AccSlab1 / Channel_14_AccSlab2 (velocity).
             "v3": {"vz_acc15": 490.80, "vz_acc16": 357.91, "vz_acc17": 470.34,
                    "Channel_13_AccSlab1": 165.0, "Channel_14_AccSlab2": 165.0}}
LP_HZ = 25.0
XI = 0.05
TAIL_S = 2.0
TAIL_FRAC = 0.05
EXPKEY_TOL_MM = 1.0
SIGN_TEST = -1.0
MODEL = "Route 2 (T_p / T_m)"
CONTROL = "Strategy F control"
US1, US2 = "US-1 (Test 9)", "US-2 (Test 12)"
COLORS = {MODEL: "#B91C1C", CONTROL: "#0F766E", US1: "#111827", US2: "#1D4ED8",
          "target": "#9A9A9A", "A": "#537EA5", "B": "#CF8152", "pushover": "#9CA3AF"}
LABEL = {MODEL: "NUM Route 2", CONTROL: "NUM Strategy F", US1: "US-1", US2: "US-2"}     # --model-label / --control-label
TEST_LS = {US1: ("-", "o"), US2: ("--", "^")}
NUM_ALPHA = 0.45
ALIAS = {"vz_acc15": "Channel_15_AccBot", "vz_acc16": "Channel_16_AccMid", "vz_acc17": "Channel_17_AccTop"}
GRID = dict(color="#E8E8E8", lw=0.7)
DPI = 160


def args_():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", default="r2tm_TpTm_full25_LS_maxwell")
    p.add_argument("--control-dir", default="stratF_full_LS_MAXWELL15_v12", help="Strategy F run-folder root (optional)")
    p.add_argument("--model-label", default="NUM Route 2")
    p.add_argument("--control-label", default="NUM Strategy F 1.5%")
    p.add_argument("--exp-dir", default="exp_processed_v3")
    p.add_argument("--pushover-pos", default="pushover_results_v5_ft0_masonrypush_TopMason/pushover_pos.csv")
    p.add_argument("--pushover-neg", default="pushover_results_v5_ft0_masonrypush_TopMason/pushover_neg.csv")
    p.add_argument("--pushover-force", choices=("applied", "reaction"), default="applied",
                   help="pushover force: |a| x m_tot (default) or F_tot_kN as logged")
    p.add_argument("--mass-set", default="v3", choices=sorted(MASS_SETS))
    p.add_argument("--loop-runs", default="1,15,18,21,23,24", help="runs for the hysteresis comparison")
    p.add_argument("--hist-runs", default="", help="runs for the history figure (default: the loop runs)")
    return p.parse_args()


def read_csv(path):
    if not os.path.isfile(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as f:
        return [r for r in csv.DictReader(l for l in f if l.strip() and not l.lstrip().startswith("#"))]


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def read_hist(path):
    """3DEC history export: 2 header lines, then time and value (whitespace- or comma-separated)."""
    best = None
    for delim in (None, ","):
        try:
            d = np.genfromtxt(path, skip_header=2, delimiter=delim)
        except Exception:
            continue
        if d is None or d.ndim != 2 or d.shape[1] < 2:
            continue
        d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
        if best is None or len(d) > len(best):
            best = d
    return (best[:, 0], best[:, 1]) if best is not None and len(best) >= 5 else (None, None)


def find_hist(folder, key):
    hits = sorted(h for h in glob.glob(os.path.join(folder, "*_" + key + ".csv")))
    if not hits:
        hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else (None, None)


def run_folder(root, n):
    hits = sorted(p for p in glob.glob(os.path.join(root, "Run{:02d}_*".format(n))) if os.path.isdir(p))
    return hits[0] if hits else None


def read_vel_table(root, n):
    hits = sorted(glob.glob(os.path.join(root, "Run{:02d}_*_vel.txt".format(n))))
    if not hits:
        return None, None
    d = np.loadtxt(hits[0], skiprows=2)
    return d[:, 0], d[:, 1]


# ============================ signal tools ===========================
def lowpass(x, fs, fc=LP_HZ):
    fc = min(fc, 0.45 * fs)
    try:
        from scipy.signal import butter, filtfilt
        b, a = butter(4, fc / (0.5 * fs))
        return filtfilt(b, a, x) if len(x) > 3 * max(len(a), len(b)) else x
    except ImportError:
        k = max(1, int(round(fs / fc / 2)))
        return np.convolve(x, np.ones(2 * k + 1) / (2 * k + 1), mode="same")


def accel_from_vel(t, v):
    fs = 1.0 / float(np.median(np.diff(t)))
    return lowpass(np.gradient(v, t), fs), fs


def sd_spectrum(ag, dt, Ts, xi=XI):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2 * np.pi / Ts; c = 2 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0] - c * v - w2 * u; um = np.abs(u)
    for i in range(1, len(ag)):
        un = u + dt * v + 0.25 * dt * dt * a; vn = v + 0.5 * dt * a
        an = (-ag[i] - c * vn - w2 * un) / (1 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = un + 0.25 * dt * dt * an; v = vn + 0.5 * dt * an; a = an; um = np.maximum(um, np.abs(u))
    return um


def sd_of_motion(t, a, Ts):
    """Sd of a table acceleration history, resampled to a uniform step, with a free tail."""
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt); au = np.interp(tu, t, a)
    return sd_spectrum(np.r_[au, np.zeros(int(TAIL_S / dt))], dt, Ts)


# ============================ model runs =============================
def model_run(root, n, masses, T1, Teff):
    folder = run_folder(root, n)
    if folder is None:
        return None
    out = dict(folder=folder)
    t, x = find_hist(folder, "rel_disp_top_mm")
    src = "rel_disp_top_mm"
    if t is None:
        t3, c3 = find_hist(folder, "Channel_3_DispTopQLeft"); _, c4 = find_hist(folder, "Channel_4_DispTopQRight")
        _, c5 = find_hist(folder, "Channel_5_DispTable")
        if t3 is None or c4 is None or c5 is None:
            return None
        k = min(len(c3), len(c4), len(c5)); t = t3[:k]; x = (0.5 * (c3[:k] + c4[:k]) - c5[:k]) * 1e3
        src = "0.5(Ch3+Ch4)-Ch5"
    out.update(t=t - t[0], x=x, edp_src=src)
    # cross-check the exported 'experiment-matched' key (history 22 by number)
    te, xe = find_hist(folder, "rel_disp_top_exp_mm")
    if te is not None:
        k = min(len(xe), len(x)); out["expkey_dev_mm"] = float(np.max(np.abs((xe[:k] - xe[0]) - (x[:k] - x[0]))))
        out["expkey_peak"] = float(np.max(np.abs(xe - xe[0])))
    # force from the absolute velocities
    mk = list(masses.items()) if isinstance(masses, dict) else list(zip(("vz_acc15", "vz_acc16", "vz_acc17"), masses))
    vs = []
    for k, _ in mk:
        v = find_hist(folder, k)
        if v[0] is None and k in ALIAS:
            v = find_hist(folder, ALIAS[k])
        vs.append(v)
    if all(v[0] is not None for v in vs):
        k = min(len(v[1]) for v in vs); tv = vs[0][0][:k]
        acc = [accel_from_vel(tv, v[1][:k])[0] for v in vs]
        F = -sum(m * a for (_, m), a in zip(mk, acc)) / 1e3          # restoring force, + with + d
    else:
        out["F_missing"] = [k for (k, _), v in zip(mk, vs) if v[0] is None]
    if all(v[0] is not None for v in vs):
        F = F - np.median(F[: max(3, len(F) // 50)])
        out.update(tF=tv - tv[0], F=F)
    # achieved table motion
    tt, vt = find_hist(folder, "Channel_12_AccTable")
    if tt is not None:
        at, fs = accel_from_vel(tt, vt)
        out["pga_ach_g"] = float(np.max(np.abs(at)) / 9.81)
        out["pgv_ach"] = float(np.max(np.abs(vt)))
        if T1 and Teff:
            sd = sd_of_motion(tt, at, [T1, Teff]) * 1e3; out["SdT1_ach_mm"], out["SdTeff_ach_mm"] = float(sd[0]), float(sd[1])
        out["t_tab"], out["v_tab"] = tt - tt[0], vt
    return out


def model_metrics(runs):
    if not runs:
        return {}
    x0 = runs[min(runs)]["x"][0]
    for n, r in runs.items():
        x = r["x"]; tail = float(np.mean(x[-max(3, len(x) // 20):]))
        r["peak_incr_mm"] = float(np.max(np.abs(x - x[0])))
        r["peak_pos_mm"], r["peak_neg_mm"] = float(np.max(x - x[0])), float(np.min(x - x[0]))
        r["resid_incr_mm"] = tail - float(x[0]); r["cum_resid_mm"] = tail - float(x0)
        if "F" in r:
            r["peak_F_kN"] = float(np.max(np.abs(r["F"])))
    return runs



# ============================ tests ==================================
def _bp(x, fs):
    from scipy.signal import butter, filtfilt
    b, a = butter(4, [0.5 / (0.5 * fs), 25.0 / (0.5 * fs)], "band"); return filtfilt(b, a, x - np.mean(x))


def _shift(x, k):
    if k == 0:
        return x
    return np.concatenate([np.full(k, x[0]), x[:-k]]) if k > 0 else np.concatenate([x[-k:], np.full(-k, x[-1])])


def load_test_run(path):
    import openpyxl
    from scipy.signal import correlate
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    rows = ws.iter_rows(values_only=True); hdr = [str(h).strip() if h is not None else "" for h in next(rows)]
    c = {h: j for j, h in enumerate(hdr)}
    rows = [r for r in rows if r and r[c["Time"]] is not None and r[c["U_avg"]] is not None]
    t = np.array([r[c["Time"]] for r in rows], float)
    d = np.array([r[c["U_avg"]] for r in rows], float) * 1e3
    ms = float(rows[0][c["mass_scale"]]) if "mass_scale" in c else 1.0
    F = SIGN_TEST * np.array([r[c["base_shear_kN"]] for r in rows], float) / ms
    fs = 1.0 / float(np.median(np.diff(t)))
    uu, ff = _bp(d, fs), _bp(F, fs)
    cc = correlate(ff, uu, "full", method="fft"); L = np.arange(-len(uu) + 1, len(uu)); m = np.abs(L) <= int(4 * fs)
    kb = int(L[m][np.argmax(cc[m])]); c0 = cc[L == 0][0] / cc[m].max()
    rho_best = cc[m].max() / (np.linalg.norm(uu) * np.linalg.norm(ff))
    k = kb if (abs(kb) > 3 and c0 < 0.5 and rho_best > 0.5) else 0
    d = _shift(d, k)
    E = np.cumsum(ff ** 2); E /= E[-1]
    i0, i1 = int(np.searchsorted(E, 0.01)), int(np.searchsorted(E, 0.99))
    pre = slice(max(i0 - int(2 * fs), 0), max(i0 - int(0.2 * fs), 1))
    F = F - np.median(F[pre])
    d, F = lowpass(d, fs), lowpass(F, fs)
    return dict(t=t - t[0], x=d, tF=t - t[0], F=F, lag=k, d_start=float(np.median(d[pre])),
                d_end=float(np.median(d[-int(2 * fs):])), win=(i0, i1), mass_scale=ms)


def read_tests(exp_dir, cache_path):
    cache = {}
    if os.path.isfile(cache_path):
        try:
            with open(cache_path, "rb") as f:
                cache = pickle.load(f)
        except Exception:
            cache = {}
    out = {US1: {}, US2: {}}
    try:
        import openpyxl, scipy  # noqa: F401
    except ImportError as e:
        print("! {}: test workbooks not read".format(e)); return out
    changed = False
    for lab, test, nmax in ((US1, "Test9", 24), (US2, "Test12", 25)):
        for n in range(1, nmax + 1):
            p = None
            for name in ("{}Run{}_processed_globalzero.xlsx".format(test, n), "{}Run{:02d}_processed_globalzero.xlsx".format(test, n)):
                if os.path.isfile(os.path.join(exp_dir, name)):
                    p = os.path.join(exp_dir, name); break
            if p is None:
                continue
            key = (os.path.abspath(p), os.path.getmtime(p))
            if key not in cache:
                try:
                    cache[key] = load_test_run(p); changed = True
                    print("  read {} run {}{}".format(test, n, "  (d-F slip {} samples corrected)".format(cache[key]["lag"]) if cache[key]["lag"] else ""))
                except Exception as e:
                    print("! {}: {}".format(p, e)); continue
            out[lab][n] = dict(cache[key])
    if changed:
        with open(cache_path, "wb") as f:
            pickle.dump(cache, f)
    return out


def test_metrics(runs):
    if not runs:
        return runs
    z = runs[min(runs)]["d_start"]
    for n, r in runs.items():
        i0, i1 = r["win"]; x = r["x"][i0:i1]
        r["peak_incr_mm"] = float(np.max(np.abs(x - r["d_start"])))
        r["resid_incr_mm"] = r["d_end"] - r["d_start"]; r["cum_resid_mm"] = r["d_end"] - z
        r["peak_F_kN"] = float(np.max(np.abs(r["F"][i0:i1])))
    return runs


PUSHOVER_FORCE = "applied"


def backbone(path):
    """pushover backbone. PUSHOVER_FORCE "applied": |a| x m_tot (m_tot from the CSV's '#' line) -- the horizontal
    force the wall carries at a settled servo point; "reaction": F_tot_kN as logged, which misses |N| sin(tilt) of the
    top-joint normal force once the top beam rotates (deformed_shape_3dec.py, 23 Oct)."""
    m_tot = None
    try:
        head = open(path).readline()
        m_tot = float(head.split("m_tot")[1].split("kg")[0])
    except Exception:
        pass
    d, F = [], []
    for r in read_csv(path):
        if r.get("converged", "1") == "0":
            continue
        dd = fnum(r.get("d_ctrl_mm"))
        if PUSHOVER_FORCE == "applied" and m_tot and "a_mps2" in r:
            ff = abs(fnum(r.get("a_mps2"))) * m_tot / 1e3
        else:
            ff = fnum(r.get("F_tot_kN"))
        if math.isfinite(dd) and math.isfinite(ff):
            d.append(dd); F.append(ff)
    d, F = np.array(d), np.array(F)
    if len(d) and np.mean(d) > 0 and "neg" in os.path.basename(path):
        d = -d
    return d, np.sign(d) * np.abs(F)


# ============================ figures ================================
def finish(ax, xl=None, yl=None, title=None, runs_axis=False, legend=True):
    ax.grid(**GRID)
    if runs_axis:
        ax.set_xlim(0.5, 25.5); ax.set_xticks(range(1, 26)); ax.axvspan(21.5, 25.5, color="#F3E6D8", alpha=0.55, lw=0)
    if xl: ax.set_xlabel(xl)
    if yl: ax.set_ylabel(yl)
    if title: ax.set_title(title, loc="left", fontweight="bold")
    if legend and ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, ncol=2, fontsize=8)


def series(d, key):
    ns = sorted(n for n in d if key in d[n] and d[n][key] is not None and math.isfinite(d[n][key]))
    return ns, [d[n][key] for n in ns]


def plot_sequence(path, sets, tg):
    fig, ax = plt.subplots(3, 1, figsize=(11, 11), sharex=True, dpi=DPI)
    spec = (("peak_incr_mm", "Within-run peak displacement (mm)", "Peak OOP displacement by run"),
            ("cum_resid_mm", "Residual displacement from run 1 start (mm)", "Cumulative residual displacement"),
            ("peak_F_kN", "Peak base shear sum m a (kN)", "Peak base shear by run (group's formula, both sides)"))
    for a, (key, yl, ttl) in zip(ax, spec):
        if key == "peak_incr_mm" and tg:
            ns = sorted(tg); a.plot(ns, [tg[n]["target_B_mm"] for n in ns], ls=":", lw=1.4, color=COLORS["target"],
                                    zorder=1, label="wave-B target Sd(T_eff)")
        for name in (CONTROL, MODEL):
            ns, v = series(sets.get(name, {}), key)
            if ns:
                a.plot(ns, v, "-", marker="s", ms=6, mec="none", lw=4.0, alpha=NUM_ALPHA, color=COLORS[name], zorder=2,
                       label=LABEL[name])
        for name in (US1, US2):
            ns, v = series(sets.get(name, {}), key)
            if ns:
                ls, mk = TEST_LS[name]
                a.plot(ns, v, ls=ls, marker=mk, mfc="white", ms=4.5, lw=1.3, color=COLORS[name], zorder=5, label=LABEL[name])
        finish(a, yl=yl, title=ttl, runs_axis=True, legend=(a is ax[0]))
    ax[-1].set_xlabel("Run number (shaded: FR76 sequence)")
    fig.suptitle("{} and {} vs US-1 and US-2".format(LABEL[MODEL], LABEL[CONTROL]), fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def plot_peak_comparison(path, sets):
    model = sets[MODEL]
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), dpi=DPI, gridspec_kw=dict(width_ratios=[1.7, 1]))
    for row, (key, lab) in enumerate((("peak_incr_mm", "peak displacement (mm)"), ("peak_F_kN", "peak base shear (kN)"))):
        a = ax[row, 0]
        ctrl = sets.get(CONTROL, {})
        for name in (US1, US2):
            d = sets.get(name, {}); ls, mk = TEST_LS[name]
            ns = [n for n in sorted(ctrl) if n in d and key in ctrl[n] and key in d[n] and d[n][key] > 0]
            if ns:
                a.plot(ns, [ctrl[n][key] / d[n][key] for n in ns], ls=ls, lw=4.0, alpha=0.3, color=COLORS[CONTROL],
                       zorder=1, label="{} / {}".format(LABEL[CONTROL], LABEL[name]))
            ns = [n for n in sorted(model) if n in d and key in model[n] and key in d[n] and d[n][key] > 0]
            if ns:
                a.plot(ns, [model[n][key] / d[n][key] for n in ns], ls=ls, marker=mk, mfc="white", ms=4.5, lw=1.4,
                       color=COLORS[MODEL], zorder=3, label="{} / {}".format(LABEL[MODEL], LABEL[name]))
        a.axhline(1, color="#555555", lw=0.8); a.axhspan(0.8, 1.25, color="#6DAD7B", alpha=0.12, label="0.8 - 1.25")
        finish(a, yl="model / test", title="Ratio of the " + lab, runs_axis=True)
        b = ax[row, 1]; hi = 0
        for name in (US1, US2):
            d = sets.get(name, {})
            ns = [n for n in sorted(model) if n in d and key in model[n] and key in d[n]]
            if not ns:
                continue
            x = [d[n][key] for n in ns]; y = [model[n][key] for n in ns]; hi = max(hi, max(x + y))
            b.plot(x, y, TEST_LS[name][1], ms=5, mfc="white", mec=COLORS[MODEL], color=COLORS[MODEL],
                   label="{} vs {}".format(LABEL[MODEL], LABEL[name]))
            for n, xi, yi in zip(ns, x, y):
                if n in (15, 21, 22, 23, 24, 25):
                    b.annotate(str(n), (xi, yi), fontsize=7, xytext=(3, 3), textcoords="offset points", color=COLORS[MODEL])
        if hi:
            g = np.array([0, hi * 1.08]); b.plot(g, g, color="#555555", lw=0.8, label="1 : 1")
            b.fill_between(g, 0.8 * g, 1.25 * g, color="#6DAD7B", alpha=0.12, lw=0)
            b.set_xlim(0, hi * 1.08); b.set_ylim(0, hi * 1.08); b.set_aspect("equal")
        finish(b, xl="test " + lab, yl="model " + lab, title="Model vs test, " + lab.split(" (")[0])
    ax[1, 0].set_xlabel("Run number (shaded: FR76 sequence)")
    fig.suptitle("{}: peak response against US-1 and US-2 ({} as translucent reference)".format(LABEL[MODEL], LABEL[CONTROL]),
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def plot_loops(path, sets, loop_runs, bb):
    runs = [n for n in loop_runs if n in sets[MODEL] or n in sets.get(CONTROL, {})]
    if not runs:
        return False
    nc = 3; nr = int(math.ceil(len(runs) / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(5.0 * nc, 4.3 * nr), dpi=DPI, squeeze=False)
    for k, n in enumerate(runs):
        a = ax[k // nc, k % nc]; peaks = []
        for name in (CONTROL, MODEL, US1, US2):
            r = sets.get(name, {}).get(n)
            if r is None or "F" not in r:
                continue
            if name in (MODEL, CONTROL):
                x = np.interp(r["tF"], r["t"], r["x"]) - r["x"][0]; F = r["F"]
                a.plot(x, F, lw=1.6, color=COLORS[name], alpha=NUM_ALPHA, zorder=2, label=LABEL[name])
            else:
                i0, i1 = r["win"]; x = r["x"][i0:i1] - r["d_start"]; F = r["F"][i0:i1]
                a.plot(x, F, lw=1.0, ls=TEST_LS[name][0], color=COLORS[name], zorder=5, label=LABEL[name])
            peaks.append((np.max(np.abs(x)), np.max(np.abs(F)), COLORS[name]))
        xs = [l.get_xdata() for l in a.lines]; ys = [l.get_ydata() for l in a.lines]
        if xs:
            a.set_xlim(*(1.15 * max(np.max(np.abs(v)) for v in xs) * np.array([-1, 1])))
            a.set_ylim(*(1.15 * max(np.max(np.abs(v)) for v in ys) * np.array([-1, 1])))
        if bb is not None:
            for j, (d, F) in enumerate(bb):
                if len(d):
                    a.plot(d, F, color=COLORS["pushover"], lw=1.0, ls=":", zorder=1, label="pushover" if j == 0 else None)
        for i, (dm, fm, c) in enumerate(reversed(peaks)):
            a.text(0.98, 0.03 + 0.065 * i, "{:.1f} mm  {:.1f} kN".format(dm, fm), transform=a.transAxes, ha="right",
                   va="bottom", fontsize=7, color=c)
        a.axhline(0, color="#555555", lw=0.5); a.axvline(0, color="#555555", lw=0.5)
        rec, sc = PROT[n]
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
    fig.suptitle("Hysteresis loops at selected runs: {} and {} vs US-1 and US-2 (corner: peak |d|, peak |F|)".format(
        LABEL[MODEL], LABEL[CONTROL]), fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    return True


def plot_input(path, model, tg):
    ns = [n for n in sorted(model) if n in tg]
    if not ns:
        return False
    fig, ax = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True, dpi=DPI)
    rec_pga = [tg[n]["PGA_tot_g"] / tg[n]["PGA_over_rec"] if tg[n].get("PGA_over_rec") else float("nan") for n in ns]
    ax[0].plot(ns, rec_pga, ls="--", lw=1.0, color="#555555", label="scaled record PGA")
    ax[0].plot(ns, [tg[n]["PGA_tot_g"] for n in ns], marker="o", ms=3.5, lw=1.4, color=COLORS["target"], label="designed signal PGA")
    ax[0].plot(ns, [model[n].get("pga_ach_g", float("nan")) for n in ns], marker="s", ms=4.5, lw=0, mfc="none",
               color=COLORS[MODEL], label="achieved table PGA (3DEC)")
    finish(ax[0], yl="PGA (g)", title="Table motion: designed vs achieved", runs_axis=True)
    for key, tk, col, lab in (("SdT1_ach_mm", "target_A_mm", COLORS["A"], "Wave A target at T1"),
                              ("SdTeff_ach_mm", "target_B_mm", COLORS["B"], "Wave B target at T_eff")):
        ax[1].plot(ns, [100 * (model[n].get(key, float("nan")) / tg[n][tk] - 1) for n in ns], marker="o", ms=3.5, lw=1.3, color=col, label=lab)
    ax[1].axhspan(-2, 2, color="#6DAD7B", alpha=0.12, label="+-2 %"); ax[1].axhline(0, color="#555555", lw=0.7)
    finish(ax[1], xl="Run number", yl="Sd error of the achieved table motion (%)",
           title="Input fit check (this is not a structural-response error)", runs_axis=True)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    return True


def plot_histories(path, root, model, sel, tg):
    sel = [n for n in sel if n in model]
    if not sel:
        return False
    fig, ax = plt.subplots(len(sel), 2, figsize=(13, 2.4 * len(sel)), squeeze=False, dpi=150)
    for i, n in enumerate(sel):
        r = model[n]; ti, vi = read_vel_table(root, n)
        if ti is not None:
            ax[i, 0].plot(ti, vi, color=COLORS["target"], lw=2.2, label="imported table velocity")
        if "t_tab" in r:
            ax[i, 0].plot(r["t_tab"], r["v_tab"], color=COLORS[MODEL], lw=1.0, label="achieved (Ch12)")
        ax[i, 1].plot(r["t"], r["x"] - r["x"][0], color=COLORS[MODEL], lw=1.1, label="2.06 m - table")
        if n in tg:
            for s in (1, -1):
                ax[i, 1].axhline(s * tg[n]["target_B_mm"], color=COLORS["B"], ls="--", lw=0.9, label="+-Sd(T_eff) target" if s == 1 else None)
        rec, sc = PROT[n]
        finish(ax[i, 0], yl="v (m/s)", title="Run {} - {} x {:.2f}".format(n, rec, sc)); finish(ax[i, 1], yl="d (mm)")
    ax[-1, 0].set_xlabel("time within run (s)"); ax[-1, 1].set_xlabel("time within run (s)")
    fig.suptitle("Route 2 (T_p / T_m): table input and response at selected runs", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    return True


# ============================ main ===================================
def collect(root, masses, T1, Teff):
    runs = {}
    if not root or not os.path.isdir(root):
        return runs
    for n, _, _ in PROTOCOL:
        r = model_run(root, n, masses, T1, Teff)
        if r is not None:
            runs[n] = r
    return model_metrics(runs)


def fmt(v, f="{:.5g}"):
    return f.format(v) if isinstance(v, (int, float)) and math.isfinite(v) else ""


def main():
    a = args_()
    root = a.run_dir
    if not os.path.isdir(root):
        raise SystemExit("run folder not found: " + root)
    out = os.path.join(root, "postproc"); os.makedirs(out, exist_ok=True)
    tg = {int(fnum(r["run"])): {k: (fnum(v) if k not in ("record", "flags", "stages") else v) for k, v in r.items()}
          for r in read_csv(os.path.join(root, "r2tm_targets.csv"))}
    log = {int(fnum(r["run"])): r for r in read_csv(os.path.join(root, "r2tm_log.csv"))}
    T1 = next((tg[n]["T1_s"] for n in tg), None); Teff = next((tg[n]["Teff_s"] for n in tg), None)
    masses = MASS_SETS[a.mass_set]
    LABEL[MODEL], LABEL[CONTROL] = a.model_label, a.control_label
    print("force: mass set '{}'; tests from {}; control {}".format(a.mass_set, a.exp_dir, a.control_dir))
    model = collect(root, masses, T1, Teff)
    if not model:
        raise SystemExit("no RunNN_* history folders in {}. Run the 3DEC driver first.".format(root))
    control = collect(a.control_dir, masses, T1, Teff)
    tests = read_tests(a.exp_dir, os.path.join(out, "_exp_cache.pkl"))
    sets = {MODEL: model, CONTROL: control, US1: test_metrics(tests[US1]), US2: test_metrics(tests[US2])}

    # ---- metrics table ----
    fields = ["run", "record", "scale", "flags", "target_A_mm", "target_B_mm", "PGA_design_g", "PGA_achieved_g",
              "SdT1_err_pct", "SdTeff_err_pct", "edp_source", "expkey_dev_mm", "driver_peak_edp_mm"]
    for name in sets:
        k = name.split(" (")[0].lower().replace(" ", "_").replace("-", "")
        fields += [k + "_peak_incr_mm", k + "_resid_incr_mm", k + "_cum_resid_mm", k + "_peak_F_kN"]
    fields += ["model_over_us1_peak_d", "model_over_us2_peak_d", "model_over_us1_peak_F", "model_over_us2_peak_F"]
    rows = []
    for n, rec, sc in PROTOCOL:
        if n not in model:
            continue
        m = model[n]; t = tg.get(n, {}); L = log.get(n, {})
        row = dict(run=n, record=rec, scale=sc, flags=t.get("flags", ""), target_A_mm=fmt(t.get("target_A_mm")),
                   target_B_mm=fmt(t.get("target_B_mm")), PGA_design_g=fmt(t.get("PGA_tot_g")), PGA_achieved_g=fmt(m.get("pga_ach_g")),
                   SdT1_err_pct=fmt(100 * (m["SdT1_ach_mm"] / t["target_A_mm"] - 1)) if "SdT1_ach_mm" in m and t else "",
                   SdTeff_err_pct=fmt(100 * (m["SdTeff_ach_mm"] / t["target_B_mm"] - 1)) if "SdTeff_ach_mm" in m and t else "",
                   edp_source=m["edp_src"], expkey_dev_mm=fmt(m.get("expkey_dev_mm")), driver_peak_edp_mm=L.get("peak_edp_mm", ""))
        for name, d in sets.items():
            k = name.split(" (")[0].lower().replace(" ", "_").replace("-", "")
            for q in ("peak_incr_mm", "resid_incr_mm", "cum_resid_mm", "peak_F_kN"):
                row[k + "_" + q] = fmt(d.get(n, {}).get(q))
        for tl, tn in (("us1", US1), ("us2", US2)):
            tr = sets[tn].get(n, {})
            row["model_over_{}_peak_d".format(tl)] = fmt(m["peak_incr_mm"] / tr["peak_incr_mm"]) if tr.get("peak_incr_mm") else ""
            row["model_over_{}_peak_F".format(tl)] = fmt(m.get("peak_F_kN", float("nan")) / tr["peak_F_kN"]) if tr.get("peak_F_kN") else ""
        rows.append(row)
    mpath = os.path.join(out, "r2tm_full25_metrics.csv")
    with open(mpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    print("{:>3} {:4} {:>5} | {:>7} {:>7} {:>7} | {:>7} {:>7} {:>7} | {:>6} {:>6}".format(
        "run", "rec", "s", "d mod", "d US1", "d US2", "F mod", "F US1", "F US2", "d/US1", "F/US1"))
    for r in rows:
        print("{:>3} {:4} {:5.2f} | {:>7} {:>7} {:>7} | {:>7} {:>7} {:>7} | {:>6} {:>6}".format(
            r["run"], r["record"], r["scale"], *[(r[k][:7] if r[k] else "-") for k in (
                "route_2_peak_incr_mm", "us1_peak_incr_mm", "us2_peak_incr_mm", "route_2_peak_F_kN", "us1_peak_F_kN",
                "us2_peak_F_kN", "model_over_us1_peak_d", "model_over_us1_peak_F")]))
    bad = [n for n in model if model[n].get("expkey_dev_mm", 0) > EXPKEY_TOL_MM]
    if bad:
        print("! rel_disp_top_exp_mm differs from rel_disp_top_mm by > {} mm on runs {}: history 22 is not the EDP "
              "(instrument_tilt_v2.dat not called); the driver's peak_edp_mm is unreliable.".format(EXPKEY_TOL_MM, sorted(bad)))
    print("model runs {}, control runs {}, US-1 runs {}, US-2 runs {}".format(len(model), len(control), len(sets[US1]), len(sets[US2])))

    # ---- figures ----
    loop_runs = [int(s) for s in a.loop_runs.split(",") if s.strip()]
    hist_runs = [int(s) for s in a.hist_runs.split(",") if s.strip()] or loop_runs
    global PUSHOVER_FORCE
    PUSHOVER_FORCE = a.pushover_force
    bb = (backbone(a.pushover_pos), backbone(a.pushover_neg))
    paths = dict(seq="fig_r2tm_full25_sequence.png", peak="fig_r2tm_full25_peak_comparison.png",
                 loops="fig_r2tm_full25_loops.png", inp="fig_r2tm_full25_input_fit.png", hist="fig_r2tm_full25_histories.png")
    paths = {k: os.path.join(out, v) for k, v in paths.items()}
    plot_sequence(paths["seq"], sets, tg)
    plot_peak_comparison(paths["peak"], sets)
    made = dict(loops=plot_loops(paths["loops"], sets, loop_runs, bb), inp=plot_input(paths["inp"], model, tg),
                hist=plot_histories(paths["hist"], root, model, hist_runs, tg))
    print("-> " + mpath)
    for k, p in paths.items():
        if made.get(k, True):
            print("-> " + p)


if __name__ == "__main__":
    main()
