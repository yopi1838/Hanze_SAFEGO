# -*- coding: utf-8 -*-
"""check_mode_shape.py -- shape of the wall's response over the height, model vs tests (16 Nov 2026).
Peak RELATIVE acceleration (a - a_table, 15 Hz low-pass) at the wall sensors 15 / 16 / 17 (0.60 / 1.26 / 2.06 m)
and the slabs 13 / 14 (0.49 / 2.28 m), normalised by the 2.06 m sensor, at a few runs. Same masses x these shapes
give the force per unit top acceleration: sum m_j phi_j (the 'participating mass' of the formula).
USAGE python check_mode_shape.py [stratF folder] [runs, e.g. 1,4,7]
      -> check_force_stiffness/check_mode_shape_<folder tag>.png / .csv"""
import os, sys, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ".")
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_force_split.py")).read().split("rows = []")[0]
exec(src)
import pandas as pd
RUNS_S = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [1, 7, 15, 19]
TAG = FOLDER.rstrip("/\\").replace("stratF_full_LS_", "")
H = {"15": 0.60, "16": 1.26, "17": 2.06}
M = {"15": 490.8, "16": 357.9, "17": 470.3, "13": 165.0, "14": 165.0}


def shape(t, c):
    f = {k: lp(t, v - np.median(v), 15.0) for k, v in c.items()}
    r = {k: np.max(np.abs(f[k] - f["12"])) for k in ("13", "14", "15", "16", "17")}
    return {k: r[k] / r["17"] for k in r}


def cm(n):
    rf = R.run_folder(FOLDER, n); tt, vt = R.find_hist(rf, "Channel_12_AccTable")
    tu = np.arange(tt[0], tt[-1], 1 / 500.0); out = {"12": np.gradient(np.interp(tu, tt, vt), tu)}
    for k, c in zip(M_MODEL, ("15", "16", "17", "13", "14")):
        t, v = R.find_hist(rf, k); out[c] = np.gradient(np.interp(tu, t, v), tu)
    return tu, out


def ct(nm, n):
    df = pd.read_excel("EXP_DATA/{}Run{}.xlsx".format(nm, n)); df.columns = df.columns.map(str).str.strip()
    return df["Time"].values.astype(float), {c: df[c].values.astype(float) * G for c in ("12", "13", "14", "15", "16", "17")}


st = {"model": dict(color="#C2410C", marker="s", lw=2.2, label="DEM " + TAG),
      "US-1": dict(color="#111827", marker="o", mfc="white", lw=1.2, label="US-1"),
      "US-2": dict(color="#1D4ED8", marker="^", mfc="white", ls="--", lw=1.2, label="US-2")}
fig, ax = plt.subplots(1, len(RUNS_S), figsize=(4.2 * len(RUNS_S), 5.2), dpi=160, sharey=True)
rows = []
for a, n in zip(ax, RUNS_S):
    for lab, tc in (("model", cm(n)), ("US-1", ct("Test9", n)), ("US-2", ct("Test12", n))):
        s = shape(*tc); pm = sum(M[k] * s[k] for k in M)
        rows.append(dict(run=n, series=lab, **{"phi" + k: round(s[k], 3) for k in s}, participating_kg=round(pm, 1)))
        a.plot([0] + [s[k] for k in ("15", "16", "17")], [0] + [H[k] for k in ("15", "16", "17")],
               **dict(st[lab], label="{} (sum m phi {:.0f} kg)".format(st[lab]["label"], pm)))
    a.axhline(2.58, color="#9CA3AF", lw=0.8, ls=":"); a.text(0.02, 2.6, "top beam 2.58 m", fontsize=7, color="#6B7280")
    a.set_title("Run {}".format(n), loc="left", fontweight="bold"); a.set_xlabel("peak a_rel / peak a_rel at 2.06 m")
    a.set_xlim(0, 1.5); a.set_ylim(0, 2.75); a.grid(color="#E8E8E8", lw=0.7); a.legend(frameon=False, fontsize=7, loc="lower right")
ax[0].set_ylabel("height (m)")
fig.suptitle(TAG + ": shape of the wall response over the height (relative acceleration, 15 Hz low-pass), normalised at 2.06 m",
             fontweight="bold", fontsize=11)
os.makedirs("check_force_stiffness", exist_ok=True)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig("check_force_stiffness/check_mode_shape_{}.png".format(TAG), bbox_inches="tight")
with open("check_force_stiffness/check_mode_shape_{}.csv".format(TAG), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
