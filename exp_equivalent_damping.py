"""exp_equivalent_damping.py -- evidence for the damping ratio: the equivalent linear SDOF of every test run.

For every run of US-1 (Test 9) and US-2 (Test 12) the measured response is fitted by a linear oscillator driven by
the MEASURED table acceleration:
        u'' + 2 xi w u' + w^2 u = -a_table(t),        d_meas(t) ~ G * u(t)
    a_table : raw channel 12 (g)          d_meas : relative displacement at 2.06 m, 0.5 (ch3 + ch4) (m)
Both band-passed 0.3-30 Hz; fitted over the strong motion (1-99 % of the table-acceleration energy) plus a 2 s tail,
least squares on the whole time history (phase and decay, not only the peak), three unknowns:
    T_eq = 2 pi / w   (equivalent period),  xi_eq (equivalent viscous damping),  G (participation x mode shape).
Several starting periods are tried; the best fit is kept, with its R^2.

READING: for the elastic runs T_eq ~ T1 and xi_eq is the wall's own viscous-type damping -- the evidence for the
Maxwell / Rayleigh ratio of the DEM. For the rocking runs xi_eq also contains the hysteretic and impact dissipation
that the DEM must produce through its contacts; giving the DEM that much VISCOUS damping would count it twice.
Raw files are only read.

USAGE   python exp_equivalent_damping.py <raw_dir> <out_dir>                    (the tests)
        python exp_equivalent_damping.py EXP_DATA exp_damping
        python exp_equivalent_damping.py --model <stratF_folder> <out_dir> [--tests <tests csv>] [--runs 1-6]
        python exp_equivalent_damping.py --model stratF_full_LS_MAXWELL15_char_tr exp_damping \
                                         --tests exp_damping/equivalent_damping_runs_acc.csv --runs 1-6
MODEL MODE fits the same oscillator to a Strategy F run folder, with the SAME channels: table velocity
(Channel_12_AccTable, exported as Z velocity) differentiated to a_table, and the relative acceleration at 2.06 m from
d/dt (vz_acc17 - table velocity). The identification bias then cancels between model and test, so the viscous damping
of the DEM can be set where the model's xi_eq meets the test's in the runs before cracking.
"""
import sys, os, glob, re, csv
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, lsim
from scipy.optimize import least_squares

TESTS = (("US-1", "Test9"), ("US-2", "Test12"))
BAND = (0.3, 30.0)
TAIL_S = 2.0
T_STARTS = (0.08, 0.10, 0.13, 0.18, 0.25, 0.35, 0.5)


RESP = "acc"     # "acc": relative acceleration a17 - a12 (accelerometers, good signal at every amplitude)
                 # "disp": relative displacement 0.5 (ch3 + ch4) (LVDTs, near their noise floor in the elastic runs)


def load(path):
    df = pd.read_excel(path); df.columns = df.columns.map(str).str.strip()
    t = pd.to_numeric(df["Time"], errors="coerce").values
    a = pd.to_numeric(df["12"], errors="coerce").values * 9.81
    if RESP == "acc":
        d = pd.to_numeric(df["17"], errors="coerce").values * 9.81 - a
    else:
        d = 0.5 * (pd.to_numeric(df["3"], errors="coerce").values + pd.to_numeric(df["4"], errors="coerce").values)
    m = np.isfinite(t) & np.isfinite(a) & np.isfinite(d)
    return t[m], a[m], d[m]


def band(x, fs):
    b, a = butter(4, [BAND[0] / (0.5 * fs), min(BAND[1], 0.45 * fs) / (0.5 * fs)], "band")
    return filtfilt(b, a, x - np.mean(x))


def sdof(t, ag, T, xi):
    """relative displacement (RESP 'disp') or relative acceleration (RESP 'acc') of the oscillator"""
    w = 2 * np.pi / T
    A = [[0.0, 1.0], [-w * w, -2 * xi * w]]; B = [[0.0], [-1.0]]
    C, D = ([[1.0, 0.0]], [[0.0]]) if RESP == "disp" else ([[-w * w, -2 * xi * w]], [[-1.0]])
    _, y, _ = lsim((A, B, C, D), ag, t)
    return y


def fit_run(t, a, d):
    fs = 1.0 / float(np.median(np.diff(t)))
    a = band(a, fs); d = band(d, fs)
    E = np.cumsum(a ** 2); E /= E[-1]
    i0 = max(int(np.searchsorted(E, 0.01)) - int(1.0 * fs), 0)
    i1 = min(int(np.searchsorted(E, 0.99)) + int(TAIL_S * fs), len(t))
    tt = t[i0:i1] - t[i0]; aa = a[i0:i1]; dd = d[i0:i1]
    scale = np.max(np.abs(dd)) or 1.0

    def res(p):
        T, xi, G = np.exp(p[0]), np.exp(p[1]), p[2]
        return (G * sdof(tt, aa, T, xi) - dd) / scale

    best = None
    for T0 in T_STARTS:
        u0 = sdof(tt, aa, T0, 0.05)
        G0 = float(np.dot(u0, dd) / max(np.dot(u0, u0), 1e-30))
        try:
            r = least_squares(res, [np.log(T0), np.log(0.05), G0],
                              bounds=([np.log(0.04), np.log(0.001), -5.0], [np.log(1.5), np.log(0.6), 5.0]),
                              x_scale=[0.3, 0.5, 0.5], max_nfev=400)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    T, xi, G = np.exp(best.x[0]), np.exp(best.x[1]), best.x[2]
    r2 = 1.0 - np.sum((best.fun * scale) ** 2) / np.sum((dd - dd.mean()) ** 2)
    return dict(T_eq_s=T, xi_eq=xi, G=G, R2=r2, peak_resp=float(np.max(np.abs(dd))),
                pga_g=float(np.max(np.abs(aa)) / 9.81), window_s=float(tt[-1]))


def _read_hist(path):
    with open(path) as f:
        f.readline(); f.readline()
        d = np.loadtxt(f)
    return d[:, 0], d[:, 1]


def load_model_run(rf, fs_out=200.0):
    """table acceleration and relative acceleration at 2.06 m of one model run, resampled to fs_out"""
    def one(key):
        hits = sorted(glob.glob(os.path.join(rf, "*_" + key + ".csv")))
        if not hits:
            raise FileNotFoundError("{}: no *_{}.csv".format(rf, key))
        return _read_hist(hits[0])
    tv, vt = one("Channel_12_AccTable")
    t7, v7 = one("vz_acc17")
    t = np.arange(max(tv[0], t7[0]), min(tv[-1], t7[-1]), 1.0 / fs_out)
    vt = np.interp(t, tv, vt); v7 = np.interp(t, t7, v7)
    a_tab = np.gradient(vt, t); a_rel = np.gradient(v7 - vt, t)
    return t, a_tab, a_rel


def main_model():
    global RESP
    RESP = "acc"
    i = sys.argv.index("--model"); folder = sys.argv[i + 1]
    rest = [a for k, a in enumerate(sys.argv[1:], 1) if k not in (i, i + 1) and not a.startswith("--")
            and sys.argv[k - 1] not in ("--tests", "--runs")]
    out = rest[0] if rest else "exp_damping"; os.makedirs(out, exist_ok=True)
    runs = None
    if "--runs" in sys.argv:
        spec = sys.argv[sys.argv.index("--runs") + 1]
        runs = set()
        for part in spec.split(","):
            a, _, b = part.partition("-"); runs.update(range(int(a), int(b or a) + 1))
    tests = pd.read_csv(sys.argv[sys.argv.index("--tests") + 1]) if "--tests" in sys.argv else None
    rows = []
    for rf in sorted(glob.glob(os.path.join(folder, "Run[0-9][0-9]_*"))):
        if not os.path.isdir(rf):
            continue
        n = int(re.search(r"Run(\d+)", os.path.basename(rf)).group(1))
        if runs and n not in runs:
            continue
        try:
            r = fit_run(*load_model_run(rf))
        except Exception as e:
            print("! run {}: {}".format(n, e)); continue
        r.update(run=n); rows.append(r)
        line = "model run {:2d}: T_eq {:.3f} s  xi_eq {:5.1f} %  G {:5.2f}  R2 {:.2f}".format(
            n, r["T_eq_s"], 100 * r["xi_eq"], r["G"], r["R2"])
        if tests is not None:
            for lab in ("US-1", "US-2"):
                q = tests[(tests.test == lab) & (tests.run == n)]
                if len(q):
                    line += "  | {} T {:.3f} xi {:5.1f} %".format(lab, q.T_eq_s.iloc[0], 100 * q.xi_eq.iloc[0])
        print(line)
    if not rows:
        raise SystemExit("no run folders read in " + folder)
    df = pd.DataFrame(rows)[["run", "T_eq_s", "xi_eq", "G", "R2", "peak_resp", "pga_g", "window_s"]]
    dst = os.path.join(out, "equivalent_damping_model_{}.csv".format(os.path.basename(os.path.normpath(folder))))
    df.to_csv(dst, index=False); print("->", dst)


def main():
    global RESP
    if "--model" in sys.argv:
        return main_model()
    raw, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
    if "--disp" in sys.argv:
        RESP = "disp"
    pk = {}
    rows = []
    for lab, test in TESTS:
        files = sorted(glob.glob(os.path.join(raw, test + "Run*.xlsx")),
                       key=lambda p: int(re.search(r"Run(\d+)", p).group(1)))
        for p in files:
            if "processed" in p.lower():
                continue
            n = int(re.search(r"Run(\d+)", p).group(1))
            try:
                r = fit_run(*load(p))
            except Exception as e:
                print("! {} run {}: {}".format(lab, n, e)); continue
            # peak relative displacement (LVDTs) as the amplitude measure, whatever is fitted
            r.update(test=lab, run=n); rows.append(r)
            print("{} run {:2d}: T_eq {:.3f} s  xi_eq {:5.1f} %  G {:5.2f}  R2 {:.2f}  peak resp {:.4g}  PGA {:.2f} g".format(
                lab, n, r["T_eq_s"], 100 * r["xi_eq"], r["G"], r["R2"], r["peak_resp"], r["pga_g"]))
    df = pd.DataFrame(rows)
    # amplitude axis: the peak displacement of each run from the disp file if present (same run numbering)
    df["peak_d_mm"] = np.nan
    for i, r in df.iterrows():
        p = os.path.join(raw, ("Test9" if r.test == "US-1" else "Test12") + "Run{}.xlsx".format(int(r.run)))
        try:
            q = pd.read_excel(p); q.columns = q.columns.map(str).str.strip()
            tq = pd.to_numeric(q["Time"], errors="coerce").values
            dq = 0.5 * (pd.to_numeric(q["3"], errors="coerce").values + pd.to_numeric(q["4"], errors="coerce").values)
            m = np.isfinite(tq) & np.isfinite(dq); fs = 1 / np.median(np.diff(tq[m]))
            df.loc[i, "peak_d_mm"] = float(np.max(np.abs(band(dq[m], fs)))) * 1e3
        except Exception:
            pass
    df = df[["test", "run", "T_eq_s", "xi_eq", "G", "R2", "peak_d_mm", "peak_resp", "pga_g", "window_s"]]
    df.to_csv(os.path.join(out, "equivalent_damping_runs_{}.csv".format(RESP)), index=False)
    print("->", os.path.join(out, "equivalent_damping_runs_{}.csv".format(RESP)))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        sty = {"US-1": dict(color="#1C2733", marker="o"), "US-2": dict(color="#2A5D8F", marker="^")}
        for lab in ("US-1", "US-2"):
            q = df[(df.test == lab) & (df.R2 >= 0.5)]; w = df[(df.test == lab) & (df.R2 < 0.5)]
            ax[0].scatter(q.peak_d_mm, 100 * q.xi_eq, label=lab, s=36, **sty[lab])
            ax[0].scatter(w.peak_d_mm, 100 * w.xi_eq, s=36, facecolors="none", edgecolors=sty[lab]["color"],
                          marker=sty[lab]["marker"], label="{} (R2 < 0.5)".format(lab) if len(w) else None)
            ax[1].scatter(q.peak_d_mm, q.T_eq_s, label=lab, s=36, **sty[lab])
            for _, r in q.iterrows():
                ax[0].annotate(str(int(r.run)), (r.peak_d_mm, 100 * r.xi_eq), xytext=(3, 3), textcoords="offset points",
                               fontsize=7, color=sty[lab]["color"])
        ax[0].axhline(1.5, color="#B5571A", ls="--", lw=1.2, label="DEM Maxwell 1.5 %")
        for a_ in ax:
            a_.set_xscale("log"); a_.set_xlabel("peak relative displacement at 2.06 m (mm)"); a_.grid(alpha=0.3, which="both")
        ax[0].set_ylabel("equivalent viscous damping (%)"); ax[1].set_ylabel("equivalent period (s)")
        ax[0].set_title("Equivalent linear SDOF fitted to each run: damping", loc="left", fontweight="bold")
        ax[1].set_title("... and period", loc="left", fontweight="bold")
        ax[0].legend(frameon=False, fontsize=8); ax[1].legend(frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(out, "fig_equivalent_damping_{}.png".format(RESP)), dpi=150)
        print("->", os.path.join(out, "fig_equivalent_damping_{}.png".format(RESP)))
    except Exception as e:
        print("! figure not written:", e)


if __name__ == "__main__":
    main()
