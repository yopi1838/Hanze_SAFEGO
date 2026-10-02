"""baseshear_exp_from_model.py -- the TEST's base shear (Moshfeghi et al. 2024, Eq. 2)
built from the model, and like-for-like loops against US-1 / US-2.

    F_exp = mass_scale * sum_j m_j * a_j
        a_j : d/dt of the ABSOLUTE out-of-plane velocity at the three accelerometer
              points (0.60 / 1.26 / 2.06 m), 50 Hz low-passed like an accelerometer chain
        m_j : tributary wall masses 408.52 / 320.66 / 404.13 kg
        mass_scale : 2.317348 -- the factor in the processed test files, applied to
                     BOTH sides. Its origin is unknown (question for Moshfeghi); it
                     scales both loops identically, so the comparison does not depend on it.

INPUTS, per run, in ROUTE_DIR (3DEC history export: 2 header lines, whitespace):
    Run<NN>_*_vz_acc15.csv, _vz_acc16.csv, _vz_acc17.csv     from instrument_baseshear_exp.dat
    Run<NN>_*_rel_disp_top_exp_mm.csv                        displacement, experiment definition
    Run<NN>_*_cstav*.csv                                     optional, drawn for reference
EXPERIMENT:  EXP_DIR/Test9Run<NN>_processed_globalzero.xlsx  (US-1)
             EXP_DIR/Test12Run<NN>_processed_globalzero.xlsx (US-2)   columns: Time, U_avg, base_shear_N, base_shear_kN

USAGE
    python baseshear_exp_from_model.py <ROUTE_DIR> [run ...]          default runs 21 22 23 24 25
    python baseshear_exp_from_model.py route2_full25_SS_NODAMP 22 24
OUTPUT
    <ROUTE_DIR>/postproc/fig_baseshear_exp_loops.png
    <ROUTE_DIR>/postproc/baseshear_exp_metrics.csv   (peak force, peak disp, secant slope per run and series)
"""
import sys, os, glob, csv
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = {"15": 408.52, "16": 320.66, "17": 404.13}      # kg, tributary wall masses (rho 1619, sensors 0.60/1.26/2.06)
MASS_SCALE = 2.317348                               # the test's factor, applied to both sides
LP_HZ = 50.0                                        # low-pass on the differentiated velocities
DT = 0.0005                                         # resampling step for the model histories
EXP_DIR = os.environ.get("EXP_DIR", "EXPRAW")          # folder with TestXXRunYY_processed_globalzero.xlsx
EXP = {"US-1 (Test 9)": "Test9", "US-2 (Test 12)": "Test12"}
EXP_SHEAR_SIGN, MODEL_SHEAR_SIGN = -1.0, -1.0       # same conventions as postprocess_route2_full.py
COL = {"US-1 (Test 9)": "#2b5f8a", "US-2 (Test 12)": "#f0a020", "model": "#b4491a", "cstav": "#8a8f97"}

def rd(p):
    d = np.genfromtxt(p, skip_header=2)
    m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1]); return d[m, 0], d[m, 1]
def find(folder, run, key):
    """channel CSV for a run: in the run's own subfolder (driver layout) or directly in the route folder."""
    pats = [os.path.join(folder, "Run{:02d}_*".format(run), "*{}*.csv".format(key)),
            os.path.join(folder, "Run{:02d}_*{}*.csv".format(run, key))]
    h = sorted(sum((glob.glob(p) for p in pats), []))
    return rd(h[0]) if h else (None, None)

def model_base_shear(folder, run):
    """t (uniform), F_exp in kN with the test's formula. None if the velocity exports are missing."""
    t = None; F = None
    for k in ("15", "16", "17"):
        tk, v = find(folder, run, "vz_acc" + k)
        if tk is None:
            return None, None
        tu = np.arange(tk[0], tk[-1], DT); vu = np.interp(tu, tk, v)
        a = np.gradient(vu, DT)
        b, c = butter(4, LP_HZ / (0.5 / DT), "low"); a = filtfilt(b, c, a)
        F = M[k] * a if F is None else F + M[k] * a; t = tu
    return t, MODEL_SHEAR_SIGN * MASS_SCALE * F / 1000.0

def exp_run(test, run):
    import openpyxl
    p = os.path.join(EXP_DIR, "{}Run{}_processed_globalzero.xlsx".format(test, run))
    if not os.path.isfile(p):
        print("  experiment file not found: {}   (set EXP_DIR=<folder> if it lives elsewhere)".format(p))
        return None, None
    ws = openpyxl.load_workbook(p, read_only=True).worksheets[0]
    r = np.array([x[:4] for x in ws.iter_rows(min_row=2, values_only=True)], float)
    return (r[:, 1] - r[0, 1]) * 1e3, EXP_SHEAR_SIGN * (r[:, 3] - r[0, 3])       # mm, kN, zeroed at run start

def secant(u, F, frac=0.9):
    """slope of the loop's outer envelope: median F/u over the samples beyond frac*|u|max, both signs."""
    m = np.abs(u) > frac * np.abs(u).max()
    return float(np.median(F[m] / u[m])) if m.sum() >= 3 else float("nan")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    folder = sys.argv[1]; runs = [int(x) for x in sys.argv[2:]] or [21, 22, 23, 24, 25]
    out_dir = os.path.join(folder, "postproc"); os.makedirs(out_dir, exist_ok=True)
    rows = []
    fig, ax = plt.subplots(1, len(runs), figsize=(4.4 * len(runs), 4.6)); ax = np.atleast_1d(ax)
    for a, run in zip(ax, runs):
        for lab, test in EXP.items():
            u, F = exp_run(test, run)
            if u is not None:
                a.plot(u, F, lw=0.7, alpha=0.85, color=COL[lab], label=lab)
                rows.append([run, lab, np.abs(u).max(), np.abs(F).max(), secant(u, F)])
        t, Fm = model_base_shear(folder, run)
        td, d = find(folder, run, "rel_disp_top_exp_mm")
        if Fm is not None and d is not None:
            dm = np.interp(t, td, d) - d[0]
            a.plot(dm, Fm, lw=0.9, color=COL["model"], label="model, sum(m a) x {:.3f}".format(MASS_SCALE))
            rows.append([run, "model sum(m a)", np.abs(dm).max(), np.abs(Fm).max(), secant(dm, Fm)])
            tc, c = find(folder, run, "cstav")
            if c is not None:
                dc = np.interp(tc, td, d) - d[0]
                a.plot(dc, MODEL_SHEAR_SIGN * (c - c[0]), lw=0.6, color=COL["cstav"], alpha=0.7, label="model, cstav (support only)")
                rows.append([run, "model cstav", np.abs(dc).max(), np.abs(c - c[0]).max(), secant(dc, MODEL_SHEAR_SIGN * (c - c[0]))])
        else:
            missing = [k for k in ("vz_acc15", "vz_acc16", "vz_acc17", "rel_disp_top_exp_mm") if find(folder, run, k)[0] is None]
            print("  run {}: missing channels {} under {}".format(run, missing, folder))
            a.text(0.5, 0.5, "missing: {}".format(", ".join(missing)), ha="center", transform=a.transAxes, fontsize=8)
        a.axhline(0, color="k", lw=0.6); a.axvline(0, color="k", lw=0.6); a.grid(alpha=0.25)
        a.set_title("run {}".format(run), loc="left", fontweight="bold"); a.set_xlabel("OOP rel. displacement at 2.06 m (mm)")
    ax[0].set_ylabel("base shear, sum(m_j a_j) x {:.3f}  (kN)".format(MASS_SCALE)); ax[0].legend(fontsize=7)
    fig.suptitle("Like-for-like base shear: the test's formula applied to both the model and the experiment", fontsize=11, fontweight="bold")
    fig.tight_layout()
    png = os.path.join(out_dir, "fig_baseshear_exp_loops.png"); fig.savefig(png, dpi=150); print("->", png)
    with open(os.path.join(out_dir, "baseshear_exp_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["run", "series", "peak_disp_mm", "peak_force_kN", "envelope_slope_kN_per_mm"])
        for r in rows: w.writerow([r[0], r[1]] + ["{:.3f}".format(x) for x in r[2:]])
    print("->", os.path.join(out_dir, "baseshear_exp_metrics.csv"))
