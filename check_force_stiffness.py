# -*- coding: utf-8 -*-
"""check_force_stiffness.py -- why is the model's peak force ~1.4 x the tests' at about the same displacement?
(16 Nov 2026). Two checks, per run, model (Strategy F folder) and tests (raw EXP_DATA, read only):

  1. SELF-CONSISTENCY (amplitude ratio = rms ratio in the strong part, signed by the correlation r) of each side: the relative acceleration at 2.06 m measured by the ACCELEROMETERS
     (a17 - a12; model: d/dt of vz_acc17 - Channel_12) against the one from the DISPLACEMENT (d'' of the LVDT mean
     0.5(ch3+ch4); model: rel_disp_top_mm). Both band-passed 1-25 Hz; slope by least squares (accelerometer on d'').
     Slope 1 = the accelerations that build the force agree with the displacements we compare.
  2. EFFECTIVE STIFFNESS: omega_eff^2 = -(ABSOLUTE a17) / d (rms ratio; a_abs = -w^2 d - 2 zeta w d') -> f_eff = sqrt / 2 pi, and the
     secant peak|F| / peak|d| from the metrics CSV. Same f_eff on both sides = same stiffness.
USAGE  python check_force_stiffness.py [stratF folder] [metrics csv]
Output check_force_stiffness/check_force_stiffness.csv and .png (linear axes).
"""
import os, sys, csv, glob
import numpy as np
from scipy.signal import butter, sosfiltfilt, resample_poly
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, ".")
import postprocess_r2tm_full25 as R

FOLDER = sys.argv[1] if len(sys.argv) > 1 else "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr"
METRICS = sys.argv[2] if len(sys.argv) > 2 else "stratF_postproc_unclamp_vs_glued_mean/stratF_metrics.csv"
RUNS = [1, 4, 7, 10, 13, 15, 17, 19, 21]
G = 9.8104
OUT = "check_force_stiffness"


def band(t, x, lo=1.0, hi=25.0):
    fs = 1.0 / np.median(np.diff(t))
    sos = butter(4, [lo / (0.5 * fs), min(hi, 0.45 * fs) / (0.5 * fs)], "band", output="sos")
    return sosfiltfilt(sos, x - np.mean(x))


def d2(t, x):
    return np.gradient(np.gradient(x, t), t)


def fit(y, x):
    """amplitude ratio y / x: rms ratio signed by the correlation (unbiased by noise in x, unlike least squares)"""
    c = float(np.corrcoef(x, y)[0, 1])
    return float(np.sign(c) * np.std(y) / np.std(x)), c


def uniform(t, *xs, fs=500.0):
    """resample to a uniform 500 Hz grid (the model histories are on the solver step)"""
    tu = np.arange(t[0], t[-1], 1.0 / fs)
    return (tu,) + tuple(np.interp(tu, t, x) for x in xs)


def strong(x, q=0.05):            # samples where |x| > q * max: drops the quiet lead-in / tail
    return np.abs(x) > q * np.max(np.abs(x))


def test_run(test, n):
    import pandas as pd
    p = os.path.join("EXP_DATA", "{}Run{}.xlsx".format(test, n))
    if not os.path.isfile(p):
        return None
    df = pd.read_excel(p); df.columns = df.columns.map(str).str.strip()
    t = df["Time"].values.astype(float)
    d = 0.5 * (df["3"].values + df["4"].values).astype(float)          # m
    a = (df["17"].values - df["12"].values).astype(float) * G           # m/s2, relative
    aa = df["17"].values.astype(float) * G                              # m/s2, absolute
    tab = df["5"].values.astype(float) if "5" in df else None
    db, ab, aab = band(t, d), band(t, a), band(t, aa); adb = band(t, d2(t, d))
    m = strong(db)
    (sl, c1), (w2, c2) = fit(ab[m], adb[m]), fit(aab[m], db[m])
    out = dict(slope_acc_vs_d2=sl, corr1=c1, w2=-w2, corr2=c2,
               d_pk_mm=1e3 * np.max(np.abs(db)), a_pk=np.max(np.abs(ab)), a_d2_pk=np.max(np.abs(adb)))
    if tab is not None:
        out["table_pk_mm"] = 1e3 * np.max(np.abs(band(t, tab, 0.2, 25)))
    return out


def model_run(n):
    rf = R.run_folder(FOLDER, n)
    t, d = R.find_hist(rf, "rel_disp_top_mm")
    t17, v17 = R.find_hist(rf, "vz_acc17")
    t12, v12 = R.find_hist(rf, "Channel_12_AccTable")
    if t is None or t17 is None or t12 is None:
        return None
    t, d = np.asarray(t, float), np.asarray(d, float) * 1e-3
    t, d, va, vt = uniform(t, d, np.interp(t, t17, v17), np.interp(t, t12, v12))
    a, aa = np.gradient(va - vt, t), np.gradient(va, t)
    db, ab, aab, adb = band(t, d), band(t, a), band(t, aa), band(t, d2(t, d))
    m = strong(db)
    (sl, c1), (w2, c2) = fit(ab[m], adb[m]), fit(aab[m], db[m])
    return dict(slope_acc_vs_d2=sl, corr1=c1, w2=-w2, corr2=c2,
                d_pk_mm=1e3 * np.max(np.abs(db)), a_pk=np.max(np.abs(ab)), a_d2_pk=np.max(np.abs(adb)))


def main():
    os.makedirs(OUT, exist_ok=True)
    met = {}
    for r in csv.DictReader(open(METRICS)):
        if r["variant"].startswith("unclamped"):
            met[int(r["run"])] = r
    rows = []
    for n in RUNS:
        for lab, res in (("model", model_run(n)), ("US-1", test_run("Test9", n)), ("US-2", test_run("Test12", n))):
            if res is None:
                continue
            r = met.get(n, {})
            key = {"model": ("peak_incr_mm", "peak_F_kN"), "US-1": ("US1_peak_mm", "US1_peak_F_kN"),
                   "US-2": ("US2_peak_mm", "US2_peak_F_kN")}[lab]
            try:
                ks = float(r[key[1]]) / float(r[key[0]])
            except (KeyError, ValueError, ZeroDivisionError):
                ks = np.nan
            res.update(run=n, series=lab, f_eff_Hz=np.sqrt(max(res["w2"], 0)) / (2 * np.pi), K_sec_kN_mm=ks)
            rows.append(res)
            print("run {:2d} {:6s} acc/d'' {:5.2f} (r {:.2f})  f_eff {:5.2f} Hz  d_pk {:6.2f} mm  a_pk {:5.2f} m/s2  "
                  "K_sec {:5.2f} kN/mm{}".format(n, lab, res["slope_acc_vs_d2"], res["corr1"], res["f_eff_Hz"], res["d_pk_mm"],
                                                  res["a_pk"], ks, "  table {:.1f} mm".format(res["table_pk_mm"])
                                                  if "table_pk_mm" in res else ""))
    keys = ["run", "series", "slope_acc_vs_d2", "corr1", "f_eff_Hz", "corr2", "w2", "d_pk_mm", "a_pk", "a_d2_pk", "K_sec_kN_mm", "table_pk_mm"]
    with open(os.path.join(OUT, "check_force_stiffness.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r.get(k), float) else r.get(k, "")) for k in keys})
    st = {"model": dict(color="#C2410C", marker="s", lw=2.2, label="DEM unclamped, 1.5 %"),
          "US-1": dict(color="#111827", marker="o", mfc="white", lw=1.2, label="US-1"),
          "US-2": dict(color="#1D4ED8", marker="^", mfc="white", ls="--", lw=1.2, label="US-2")}
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), dpi=160)
    for k, (key, yl, ttl) in enumerate((("slope_acc_vs_d2", "accelerometer / d'' (-)", "1. Accelerations agree with displacements?"),
                                        ("f_eff_Hz", "f_eff (Hz)", "2. Effective frequency, sqrt(a_rel / d)"),
                                        ("K_sec_kN_mm", "peak F / peak d (kN/mm)", "Secant: peak force / peak displacement"))):
        for s in ("model", "US-1", "US-2"):
            rr = [r for r in rows if r["series"] == s]
            ax[k].plot([r["run"] for r in rr], [r[key] for r in rr], **st[s])
        if key == "slope_acc_vs_d2":
            ax[k].axhline(1.0, color="#555555", lw=0.8)
        ax[k].set_title(ttl, loc="left", fontweight="bold", fontsize=10); ax[k].set_ylabel(yl); ax[k].set_xlabel("run")
        ax[k].set_xticks(RUNS); ax[k].grid(color="#E8E8E8", lw=0.7); ax[k].set_ylim(bottom=0)
    ax[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "check_force_stiffness.png"), bbox_inches="tight")
    print("->", OUT)


if __name__ == "__main__":
    main()
