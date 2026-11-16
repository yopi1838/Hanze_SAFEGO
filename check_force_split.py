# -*- coding: utf-8 -*-
"""check_force_split.py -- split the peak base shear into the TABLE part (sum m x a_table) and the WALL RESPONSE part
(sum m_j x (a_j - a_table)), model vs tests, and compare the table accelerations themselves. 16 Nov 2026.
Same masses on both sides (group formula v3: 490.8/357.9/470.3 kg on 15/16/17, 165 kg on 13/14; US-2 475.5/373.2/470.3).
All signals 25 Hz low-pass (as the force); table PGA from the same filtered signal.
USAGE python check_force_split.py [stratF folder]   -> check_force_stiffness/check_force_split.csv"""
import os, sys, csv
import numpy as np
from scipy.signal import butter, sosfiltfilt
sys.path.insert(0, ".")
import postprocess_r2tm_full25 as R
FOLDER = sys.argv[1] if len(sys.argv) > 1 else "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr"
RUNS = [1, 4, 7, 15, 19, 21]
G = 9.8104
M_MODEL = {"vz_acc15": 490.80, "vz_acc16": 357.91, "vz_acc17": 470.34, "Channel_13_AccSlab1": 165.0, "Channel_14_AccSlab2": 165.0}
M_TEST = {"Test9": {15: 490.80, 16: 357.91, 17: 470.34, 13: 165.0, 14: 165.0},
          "Test12": {15: 475.5, 16: 373.2, 17: 470.3, 13: 165.0, 14: 165.0}}


def lp(t, x, fc=25.0):
    fs = 1.0 / np.median(np.diff(t))
    return sosfiltfilt(butter(4, fc / (0.5 * fs), output="sos"), x)


def split(t, acc, a_tab, masses):
    """acc: {key: absolute accel}, returns peak |F|, peak |F_table|, peak |F_rel| (kN) and PGA (g)"""
    a_tab = lp(t, a_tab - np.median(a_tab))
    Ft = sum(masses.values()) * a_tab / 1e3
    Fr = sum(m * (lp(t, acc[k] - np.median(acc[k])) - a_tab) for k, m in masses.items()) / 1e3
    F = Ft + Fr
    return np.max(np.abs(F)), np.max(np.abs(Ft)), np.max(np.abs(Fr)), np.max(np.abs(a_tab)) / G


def model(n):
    rf = R.run_folder(FOLDER, n)
    tt, vt = R.find_hist(rf, "Channel_12_AccTable")
    if tt is None:
        return None
    tu = np.arange(tt[0], tt[-1], 1 / 500.0)
    a_tab = np.gradient(np.interp(tu, tt, vt), tu)
    acc = {}
    for k in M_MODEL:
        t, v = R.find_hist(rf, k)
        acc[k] = np.gradient(np.interp(tu, t, v), tu)
    return split(tu, acc, a_tab, M_MODEL)


def test(name, n):
    import pandas as pd
    p = os.path.join("EXP_DATA", "{}Run{}.xlsx".format(name, n))
    if not os.path.isfile(p):
        return None
    df = pd.read_excel(p); df.columns = df.columns.map(str).str.strip()
    t = df["Time"].values.astype(float)
    acc = {k: df[str(k)].values.astype(float) * G for k in M_TEST[name]}
    return split(t, acc, df["12"].values.astype(float) * G, M_TEST[name])


rows = []
print("run series   peak F   F_table  F_wall   PGA(g)   [F_table = sum m x a_table, F_wall = sum m x (a - a_table)]")
for n in RUNS:
    for lab, r in (("model", model(n)), ("US-1", test("Test9", n)), ("US-2", test("Test12", n))):
        if r is None:
            continue
        rows.append(dict(run=n, series=lab, F=r[0], F_table=r[1], F_wall=r[2], pga_g=r[3]))
        print("{:3d} {:6s} {:8.2f} {:8.2f} {:8.2f} {:8.3f}".format(n, lab, *r))
os.makedirs("check_force_stiffness", exist_ok=True)
with open("check_force_stiffness/check_force_split.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
    for r in rows:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
