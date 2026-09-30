# -*- coding: ascii -*-
"""
fig_chord_tilt_validation.py -- ONE figure: chord tilt per run, model vs the
US-1 potentiometers. This is the DEFENSIBLE tilt comparison; the inclinometer
(the paper's Fig. 16) measures a local rotation of the instrumented unit and
is a separate quantity -- see instrument_tilt_v2.dat.

Both sides are the same definition:
    chord tilt = atan( u_top / 2.06 ),  u_top = mean(Ch3, Ch4) - Ch5   [deg]
Experiment: exp_Test9_tilt.csv, made by exp_tilt_from_raw.py from the raw
Test9RunNN.xlsx channels. (exp_Test12_tilt.csv for US-2 is plotted if present.)
Model: computed from the displacement channel rather than the tilt channel, so
it works with whatever the export actually wrote --
    rel_disp_top_exp_mm  (index 22, experiment-matched, preferred)
    rel_disp_top_mm      (index 16, centreline; differs by <0.6%)
    0.5(Ch3+Ch4) - Ch5   rebuilt from the raw channel CSVs
Where a tilt channel (tiltx_full_wall / tilt_full_wall) is also present it is
printed as a cross-check; the figure uses the displacement-derived value.

    peak     = max |tilt(t) - tilt(t0)|                within the run
    residual = mean(last 5%) - value at the start of RUN 1   (cumulative)

Usage: python fig_chord_tilt_validation.py
Output: fig_chord_tilt_validation.png, chord_tilt_validation.csv
"""
import os, re, csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
MODELS = [  # (label, results dir, colour, marker)
    ("3DEC model, record (Strategy F)", "stratF_full_results_US1", "#5b7f74", "s"),
    ("3DEC model, Route 2",             "route2_full25_revised",   "#b85042", "^"),
]
EXPS = [    # (label, csv, colour, marker)   columns: run, peak_tilt_full_wall, resid_tilt_full_wall
    ("US-1 (Test 9), potentiometers",  "exp_Test9_tilt.csv",  "royalblue", "o"),
    ("US-2 (Test 12), potentiometers", "exp_Test12_tilt.csv", "#f0a020",   "D"),
]
H_M        = 2.06
DISP_KEYS  = ["rel_disp_top_exp_mm", "rel_disp_top_mm"]      # mm, absolute (not reset between runs)
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
CHECK_KEYS = ["tiltx_full_wall", "tilt_full_wall"]           # deg, cross-check only
TAIL_FRAC  = 0.05
FR76_FROM  = 22
OUT_PNG, OUT_CSV = "fig_chord_tilt_validation.png", "chord_tilt_validation.csv"
RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")

# ============================ READERS ================================
def read_hist(p):
    try:
        d = np.genfromtxt(str(p), skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None

def find(folder, key):
    hits = sorted(Path(folder).glob("*" + key + "*.csv"))
    return read_hist(hits[0]) if hits else None

def top_disp_mm(folder):
    """mean top-quarter displacement relative to the table, mm, absolute."""
    for k in DISP_KEYS:
        u = find(folder, k)
        if u is not None:
            return u, k
    c3, c4, c5 = find(folder, CH3), find(folder, CH4), find(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, None
    n = min(len(c3), len(c4), len(c5))
    return (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0, "0.5(Ch3+Ch4)-Ch5"

def model_series(results_dir):
    """{run: (peak_deg, resid_deg)}, source name, and the tilt-channel cross-check."""
    out, base, src, chk = {}, None, set(), {}
    for f in sorted(Path(results_dir).iterdir()) if Path(results_dir).is_dir() else []:
        m = RUN_RE.match(f.name)
        if not (f.is_dir() and m):
            continue
        u, k = top_disp_mm(f)
        if u is None:
            continue
        src.add(k)
        if base is None:
            base = float(u[0])
        th = np.degrees(np.arctan(u / 1000.0 / H_M))
        th0 = np.degrees(np.arctan(base / 1000.0 / H_M))
        out[int(m.group(1))] = (float(np.max(np.abs(th - th[0]))),
                                float(np.mean(th[int((1 - TAIL_FRAC) * len(th)):]) - th0))
        for ck in CHECK_KEYS:
            x = find(f, ck)
            if x is not None:
                chk.setdefault(ck, {})[int(m.group(1))] = float(np.max(np.abs(x - x[0])))
                break
    return out, "/".join(sorted(src)) or "-", chk

def exp_series(path):
    if not os.path.isfile(path):
        return {}
    out = {}
    for row in csv.DictReader(open(path)):
        try:
            rn = int(float(row["run"]))
            out[rn] = (float(row["peak_tilt_full_wall"]), float(row["resid_tilt_full_wall"]))
        except Exception:
            pass
    return out

# ============================ COLLECT ================================
print("chord tilt = atan(u_top / {:.2f} m), same definition both sides\n".format(H_M))
series = []
for lab, d, c, mk in MODELS:
    v, src, chk = model_series(d)
    print("  {:34s} {:26s} {:2d} runs  from {}".format(lab, d, len(v), src))
    if v and chk:
        ck = list(chk)[0]; common = [r for r in v if r in chk[ck] and v[r][0] > 0.01]
        if common:
            ratio = np.median([chk[ck][r] / v[r][0] for r in common])
            print("      cross-check {} / displacement-derived peak = {:.3f}  ({})".format(
                ck, ratio, "consistent" if 0.95 < ratio < 1.05 else "DIFFERENT -- the channel is not this chord"))
    series.append((lab, v, c, mk, False))
for lab, p, c, mk in EXPS:
    v = exp_series(p)
    print("  {:34s} {:26s} {:2d} runs".format(lab, p, len(v)))
    if v:
        series.append((lab, v, c, mk, True))
if not any(v for _, v, _, _, _ in series):
    raise SystemExit("nothing found -- check the folder names in MODELS and that exp_Test9_tilt.csv exists "
                     "(make it with exp_tilt_from_raw.py)")

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["series", "run", "peak_chord_tilt_deg", "resid_chord_tilt_deg"])
    for lab, v, c, mk, is_exp in series:
        for rn in sorted(v):
            w.writerow([lab, rn, round(v[rn][0], 5), round(v[rn][1], 5)])
print("-> " + OUT_CSV)

# ============================ FIGURE =================================
fig, axs = plt.subplots(1, 2, figsize=(15, 5.6), dpi=150)
for ax, idx, ttl in ((axs[0], 0, "(a) peak chord tilt within each run"),
                     (axs[1], 1, "(b) cumulative chord tilt (residual, from the run-1 zero)")):
    for lab, v, c, mk, is_exp in series:
        if not v:
            continue
        ks = sorted(v)
        ax.plot(ks, [v[r][idx] for r in ks], marker=mk, ms=5, lw=1.5,
                ls="-" if is_exp else "-", alpha=1.0 if is_exp else 0.95, color=c, label=lab)
    ax.axhline(0, color="k", lw=0.6)
    ax.axvspan(FR76_FROM - 0.5, 25.5, color="#f3e6d8", alpha=0.55, lw=0)
    ax.set_xlabel("Test Run No"); ax.set_ylabel("chord tilt at 2.06 m (deg)")
    ax.set_xticks(range(1, 26)); ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(color="#eeeeee"); ax.legend(fontsize=8.5, frameon=False, loc="upper left")
    ax.set_title(ttl, loc="left", fontsize=10.5, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axs[0].text(FR76_FROM + 1.5, 0.02, "FR76", ha="center", fontsize=8, color="#8a6a4a",
            transform=axs[0].get_xaxis_transform())
fig.suptitle("Chord tilt validation: model vs the potentiometers, same definition atan(u_top / 2.06 m)",
             fontsize=12, fontweight="bold")
fig.text(0.01, 0.01, "Not comparable with the paper's Fig. 16: the tiltmeter at Z = 2.43 m reads the local rotation of "
                     "the instrumented unit, a different quantity (instrument_tilt_v2.dat).", fontsize=8, color="#666")
plt.tight_layout(rect=(0, 0.03, 1, 0.95)); plt.savefig(OUT_PNG); print("-> " + OUT_PNG)