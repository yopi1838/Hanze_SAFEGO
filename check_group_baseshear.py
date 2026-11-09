# check_group_baseshear.py -- the group's base-shear script (Ihsan / Eleni, 26 Oct 2026), verbatim formula, run on
# EXP_DATA, against (1) exp_processed_v3 and (2) the test peaks the Strategy F post-processor reports.
import os, csv, glob
import numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
G, H, WALL_WEIGHT, S1, S2 = 9.8104, 2.58, 12.94, 165.0, 165.0

def read_heights(test):
    info = pd.read_excel("EXP_DATA/Test{}_Info.xlsx".format(test), header=None)
    return tuple(float(info.iloc[r, 5]) for r in (15, 16, 17))

def group_V(df, h):                                  # their calculate_base_shear, unchanged
    h1, h2, h3 = h
    f1 = (h1 + (h2 - h1) / 2) / H; f2 = ((h2 - h1) / 2 + (h3 - h2) / 2) / H; f3 = (H - h3 + (h3 - h2) / 2) / H
    a = {k: pd.to_numeric(df[k], errors="coerce") for k in ("13", "14", "15", "16", "17")}
    return G / 1000 * (S1 * a["13"] + S2 * a["14"]) + WALL_WEIGHT * (f1 * a["15"] + f2 * a["16"] + f3 * a["17"])

met = {}
for r in csv.DictReader(open("stratF_postproc_pK11_1p5/stratF_metrics.csv")):
    met[int(r["run"])] = r
rows = []
for test, lab, nmax in ((9, "US1", 24), (12, "US2", 25)):
    h = read_heights(test); print(lab, "heights", h)
    for n in range(1, nmax + 1):
        df = pd.read_excel("EXP_DATA/Test{}Run{}.xlsx".format(test, n)); df.columns = df.columns.map(str).str.strip()
        V = group_V(df, h).values.astype(float); t = df["Time"].values.astype(float); ok = np.isfinite(V)
        V, t = V[ok], t[ok]
        raw = float(np.nanmax(np.abs(V)))
        fs = 1 / np.median(np.diff(t)); b, a = butter(4, 25.0 / (0.5 * fs))
        Vf = filtfilt(b, a, V - np.median(V[: int(1.5 * fs)]))
        filt = float(np.max(np.abs(Vf)))
        p = "exp_processed_v3/Test{}Run{}_processed_globalzero.xlsx".format(test, n)
        v3 = float(np.nanmax(np.abs(pd.read_excel(p)["base_shear_kN"].values.astype(float)))) if os.path.isfile(p) else np.nan
        pp = float(met.get(n, {}).get(lab + "_peak_F_kN") or np.nan)
        a17 = pd.to_numeric(df["17"], errors="coerce").values.astype(float)[ok]
        Vb = filtfilt(b, a, (V + G / 1000 * 61.4 * a17) - np.median((V + G / 1000 * 61.4 * a17)[: int(1.5 * fs)]))
        rows.append((lab, n, raw, v3, filt, pp, float(np.max(np.abs(Vb)))))
print("test run | group raw | v3 file raw | group, offset + 25 Hz | post-processor | + beam 61.4 kg, filtered")
for r in rows:
    print("{} {:2d} | {:6.2f} | {:6.2f} | {:6.2f} | {:6.2f} | {:6.2f}".format(*r))
