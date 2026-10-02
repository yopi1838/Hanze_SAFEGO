"""Build the TEST's base shear (Moshfeghi 2024 Eq. 2) from the model's sensor-point
velocities, and plot like-for-like loops against US-1 / US-2.

    F_exp = mass_scale * sum_j m_j * a_j        a_j = d/dt of absolute OOP velocity at the
                                                 three accelerometer points (0.60, 1.26, 2.06 m)

Inputs (per run, from instrument_baseshear_exp.dat exports in ROUTE_DIR):
    <label>_vz_acc15.csv, _vz_acc16.csv, _vz_acc17.csv     (3DEC export: 2 header lines, whitespace)
    <label>_rel_disp_top_exp_mm.csv                          (displacement, experiment definition)
Experiment: TestXXRunYY_processed_globalzero.xlsx  (U_avg, base_shear_kN)

Usage:  python baseshear_exp_from_model.py <ROUTE_DIR> 21 22 23 24 25
"""
import sys, os, glob
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = {"15": 408.52, "16": 320.66, "17": 404.13}     # tributary wall masses, kg (rho 1619, sensors 0.60/1.26/2.06)
MASS_SCALE = 2.317348                              # the test's factor -- applied to both sides
LP_HZ = 50.0                                       # low-pass on the differentiated velocities
# Processed experimental workbooks live here in the current SAFEGo workspace.
# Override with SAFEGO_EXP_DIR when running against a different copy.
EXP_DIR = os.environ.get("SAFEGO_EXP_DIR", os.path.join("EXP_DATA", "processed_globalzero"))
EXP = {"US-1 (Test 9)": "Test9", "US-2 (Test 12)": "Test12"}
EXP_SHEAR_SIGN, MODEL_SHEAR_SIGN = -1.0, -1.0      # same conventions as postprocess_route2_full.py

def rd(p):
    d = np.genfromtxt(p, skip_header=2); m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1]); return d[m, 0], d[m, 1]
def find(folder, run, key):
    h = sorted(glob.glob(os.path.join(folder, "Run{:02d}_*{}*.csv".format(run, key))))
    return rd(h[0]) if h else (None, None)
def model_base_shear(folder, run):
    t = None; F = None
    for k in ("15", "16", "17"):
        tk, v = find(folder, run, "vz_acc" + k)
        if tk is None:
            return None, None
        dt = 0.0005; tu = np.arange(tk[0], tk[-1], dt); vu = np.interp(tu, tk, v)
        a = np.gradient(vu, dt)
        b, c = butter(4, LP_HZ / (0.5 / dt), "low"); a = filtfilt(b, c, a)
        F = M[k] * a if F is None else F + M[k] * a; t = tu
    return t, MODEL_SHEAR_SIGN * MASS_SCALE * F / 1000.0
def exp_run(test, run):
    import openpyxl
    p = os.path.join(EXP_DIR, "{}Run{}_processed_globalzero.xlsx".format(test, run))
    if not os.path.isfile(p):
        return None, None
    ws = openpyxl.load_workbook(p, read_only=True).worksheets[0]
    r = np.array([x[:4] for x in ws.iter_rows(min_row=2, values_only=True)], float)
    return r[:, 1] * 1e3 - r[0, 1] * 1e3, EXP_SHEAR_SIGN * (r[:, 3] - r[0, 3])

if __name__ == "__main__":
    folder = sys.argv[1]; runs = [int(x) for x in sys.argv[2:]] or [21, 22, 23, 24, 25]
    fig, ax = plt.subplots(1, len(runs), figsize=(4.2 * len(runs), 4.4))
    ax = np.atleast_1d(ax)
    for a, run in zip(ax, runs):
        for lab, test in EXP.items():
            u, F = exp_run(test, run)
            if u is not None:
                a.plot(u, F, lw=0.7, alpha=0.8, label=lab)
        t, Fm = model_base_shear(folder, run)
        td, d = find(folder, run, "rel_disp_top_exp_mm")
        if Fm is not None and d is not None:
            dm = np.interp(t, td, d) - d[0]
            a.plot(dm, Fm, lw=0.9, color="#b4491a", label="model, sum(m a) same formula")
            tc, c = find(folder, run, "cstav")
            if c is not None:
                a.plot(np.interp(tc, td, d) - d[0], MODEL_SHEAR_SIGN * (c - c[0]), lw=0.6, color="#8a8f97", alpha=0.7, label="model, cstav (support)")
        else:
            a.text(0.5, 0.5, "no vz_acc exports\nfor this run", ha="center", transform=a.transAxes)
        a.axhline(0, color="k", lw=0.6); a.axvline(0, color="k", lw=0.6); a.grid(alpha=0.25)
        a.set_title("run {}".format(run), loc="left", fontweight="bold"); a.set_xlabel("OOP rel. displacement at 2.06 m (mm)")
    ax[0].set_ylabel("base shear, sum(m a) x {:.3f} (kN)".format(MASS_SCALE)); ax[0].legend(fontsize=7)
    fig.suptitle("Like-for-like base shear: test formula applied to model and experiment", fontsize=11, fontweight="bold")
    fig.tight_layout(); out = os.path.join(folder, "postproc", "fig_baseshear_exp_loops.png"); os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150); print("->", out)
