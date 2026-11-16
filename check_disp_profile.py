# -*- coding: utf-8 -*-
"""check_disp_profile.py -- displacement profile over the height, model vs tests, and the model's top beam
(Channel_20_TopBeam, 2.68 m) relative to the table. 16 Nov 2026.
Test displacement channels (TestX_Info.xlsx): 1 bottom quarter 0.66 m, 2 mid 1.26 m, 3/4 top quarter 2.06 m (left/right),
5 shake table; the tests have NO top-beam channel. Model: block histories at the same points (absolute, table removed
here with Channel_5) and Channel_20 on the top beam.
Shape phi_k = regression of d_k on d(2.06 m) (0.5-25 Hz band-pass, relative to the table): the in-phase share of each
point's motion, so the shape is not mixed from peaks at different instants.
USAGE python check_disp_profile.py [stratF folder] [runs, e.g. 1,7,15,19]
      -> check_force_stiffness/check_disp_profile_<folder tag>.png / .csv"""
import os, sys, csv
import numpy as np
from scipy.signal import butter, sosfiltfilt
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, ".")
import postprocess_r2tm_full25 as R
FOLDER = sys.argv[1] if len(sys.argv) > 1 else "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr"
RUNS = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [1, 7, 15, 19]
TAG = FOLDER.rstrip("/\\").replace("stratF_full_LS_", "")
Z = {"1": 0.66, "2": 1.26, "34": 2.06, "20": 2.68}


def band(t, x, lo=0.5, hi=25.0):
    fs = 1.0 / np.median(np.diff(t))
    return sosfiltfilt(butter(4, [lo / (0.5 * fs), min(hi, 0.45 * fs) / (0.5 * fs)], "band", output="sos"), x - np.mean(x))


def shape(t, d):
    f = {k: band(t, v) for k, v in d.items()}
    ref = f["34"]; m = np.abs(ref) > 0.05 * np.max(np.abs(ref))
    phi = {k: float(np.dot(f[k][m], ref[m]) / np.dot(ref[m], ref[m])) for k in f}
    pk = {k: 1e3 * float(np.max(np.abs(f[k]))) for k in f}
    return phi, pk


def model(n):
    rf = R.run_folder(FOLDER, n)
    t5, d5 = R.find_hist(rf, "Channel_5_DispTable")
    if t5 is None:
        return None
    tu = np.arange(t5[0], t5[-1], 1 / 500.0); tab = np.interp(tu, t5, d5)
    d = {}
    for k, name in (("1", "Channel_1_DispBot"), ("2", "Channel_2_DispMid"), ("20", "Channel_20_TopBeam")):
        t, v = R.find_hist(rf, name)
        if t is not None:
            d[k] = np.interp(tu, t, v) - tab
    t, v = R.find_hist(rf, "rel_disp_top_mm"); d["34"] = np.interp(tu, t, v) * 1e-3   # already relative, mm -> m
    return shape(tu, d)


def test(nm, n):
    import pandas as pd
    df = pd.read_excel("EXP_DATA/{}Run{}.xlsx".format(nm, n)); df.columns = df.columns.map(str).str.strip()
    t = df["Time"].values.astype(float)
    d = {"1": df["1"].values.astype(float), "2": df["2"].values.astype(float),
         "34": 0.5 * (df["3"].values + df["4"].values).astype(float)}
    return shape(t, d)


st = {"model": dict(color="#C2410C", marker="s", lw=2.2, label="DEM unclamped, 1.5 %"),
      "US-1": dict(color="#111827", marker="o", mfc="white", lw=1.2, label="US-1"),
      "US-2": dict(color="#1D4ED8", marker="^", mfc="white", ls="--", lw=1.2, label="US-2")}
fig, ax = plt.subplots(1, len(RUNS), figsize=(4.2 * len(RUNS), 5.4), dpi=160, sharey=True)
rows = []
for a, n in zip(ax, RUNS):
    for lab, res in (("model", model(n)), ("US-1", test("Test9", n)), ("US-2", test("Test12", n))):
        if res is None:
            continue
        phi, pk = res
        ks = [k for k in ("1", "2", "34", "20") if k in phi]
        rows.append(dict(run=n, series=lab, **{"phi_" + k: round(phi[k], 3) for k in ks},
                         **{"peak_mm_" + k: round(pk[k], 3) for k in ks}))
        a.plot([0] + [phi[k] for k in ks], [0] + [Z[k] for k in ks], **st[lab])
        if "20" in phi:
            a.annotate("top beam: {:.2f} x d(2.06 m)\n(peak {:.3f} mm rel. to table)".format(phi["20"], pk["20"]),
                       (phi["20"], 2.68), xytext=(0.25, 2.45), fontsize=7, color=st[lab]["color"])
    a.axhline(2.58, color="#9CA3AF", lw=0.8, ls=":")
    a.set_title("Run {}".format(n), loc="left", fontweight="bold"); a.set_xlabel("d / d(2.06 m), in phase")
    a.set_xlim(-0.1, 1.5); a.set_ylim(0, 2.8); a.grid(color="#E8E8E8", lw=0.7)
ax[0].set_ylabel("height (m)"); ax[0].legend(frameon=False, fontsize=7, loc="center right")
fig.suptitle(TAG + ": displacement profile relative to the table (LVDT points 0.66 / 1.26 / 2.06 m; model also top beam 2.68 m)",
             fontweight="bold", fontsize=11)
os.makedirs("check_force_stiffness", exist_ok=True)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig("check_force_stiffness/check_disp_profile_{}.png".format(TAG), bbox_inches="tight")
keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("run", "series"), k))
with open("check_force_stiffness/check_disp_profile_{}.csv".format(TAG), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
