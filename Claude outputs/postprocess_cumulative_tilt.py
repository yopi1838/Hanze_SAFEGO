# -*- coding: ascii -*-
"""
postprocess_cumulative_tilt.py -- the model's CUMULATIVE (residual) tilt over
the 25-run sequence, every tilt channel, for the record-driven control and
for Route 2, beside US-1 from both of its instruments.

Model channels (instrument_tilt_v2.dat exports, whichever exist in a run
folder): tilt_bot_seg, tilt_low_seg, tilt_up_seg, tilt_beam_seg,
tilt_full_wall. They are ABSOLUTE (each run starts where the previous one
ended), so the cumulative tilt after run N is
    resid_N = mean of the last 5% of run N  -  first sample of run 1
and the peak within run N is max |tilt - tilt(t0 of run N)|.
A proxy from the displacement channel is computed alongside for the model,
atan(u_2.06m / 2.06), and the ratio tilt_full_wall / proxy is printed: if it
is ~1 the channel is the base-to-2.06 m angle; if not, the channel is
something else and the label says so.

US-1:
    potentiometers  exp_Test9_tilt.csv      (exp_tilt_from_raw.py; column
                    resid_tilt_full_wall = atan(u_top/2.06) from run-1 zero)
    inclinometer    exp_Test9_incl_tilt.csv (exp_tilt_from_inclinometer.py;
                    resid_tilt_y_deg, the paper's curve; ~100x larger than the
                    potentiometers say, see that script's docstring)
Either file may be absent; the panel then shows what is there.

OUTPUT (OUT_DIR = cumulative_tilt_postproc/)
    cumulative_tilt.csv          run x {series} x {channel}: resid_deg, peak_deg
    fig_cumtilt_sequence.png     (a) cumulative full-wall tilt, model vs US-1
                                 potentiometers; (b) the same on the paper's tilt
                                 figure: US-1 from the log, US-2 / ST-1-FDM / ST-2-HB
                                 digitised from the figure (embedded, no file needed);
                                 (c)(d) per-segment cumulative tilt, control / Route 2
    fig_cumtilt_history.png      tilt_full_wall stitched over the whole sequence,
                                 control and Route 2, run boundaries marked
Usage: python postprocess_cumulative_tilt.py
"""
import os, re, csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
SERIES = [  # (label, results dir, colour, marker)
    ("record on the model (Strategy F)", "stratF_full_results_US1", "#5b7f74", "s"),
    ("Route 2 (full sequence)",          "route2_full25_revised",   "#b85042", "^"),
]
TILT_KEYS   = ["tilt_bot_seg", "tilt_low_seg", "tilt_up_seg", "tilt_beam_seg", "tilt_full_wall"]
EXP_POT_CSV = "exp_Test9_tilt.csv"        # potentiometer-based (resid_tilt_full_wall, peak_tilt_full_wall)
EXP_INC_CSV = "exp_Test9_incl_tilt.csv"   # inclinometer-based (resid_tilt_y_deg)
KEY_DISP_EXP, KEY_DISP = "rel_disp_top_exp_mm", "rel_disp_top_mm"
H_EDP_M     = 2.06
TAIL_FRAC   = 0.05
OUT_DIR     = "cumulative_tilt_postproc"; os.makedirs(OUT_DIR, exist_ok=True)
RUN_RE      = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")
SEG_COL     = {"tilt_bot_seg": "#d98c21", "tilt_low_seg": "#5b7f74", "tilt_up_seg": "#b85042", "tilt_beam_seg": "#7a5c99", "tilt_full_wall": "k", "proxy_2.06m": "#888888"}
COL_EXP     = "royalblue"; COL_INC = "#e0201b"

# ---- the paper's tilt figure (Moshfeghi et al. 2024), for panel (b) when the CSV is absent
# US-1: exact, from the raw inclinometer log (exp_tilt_from_inclinometer.py). US-2 and the paper's own
# numerical curves ST-1-FDM / ST-2-HB: digitised from the figure by colour at each run (about +-0.05 deg;
# runs hidden behind other markers interpolated). Residual tilt in degrees, run-1 zero.
PAPER_CURVES = {
    "US-1 (inclinometer, from the raw log)": {"1": 0.05, "2": 0.05, "3": 0.07, "4": 0.05, "5": 0.08, "6": 0.09, "7": 0.06, "8": 0.09, "9": 0.09, "10": 0.21, "11": 0.23, "12": 0.28, "13": 0.37, "14": 0.47, "15": 0.58, "16": 0.74, "17": 1.03, "18": 1.51, "19": 2.07, "20": 2.39, "21": 2.76, "22": 2.53, "23": 2.4, "24": 6.04},
    "US-2 (paper figure, digitised)": {"1": 0.0, "2": 0.1, "3": 0.1, "4": 0.09, "5": 0.12, "6": -0.05, "7": 0.11, "8": 0.12, "9": 0.17, "10": 0.26, "11": 0.12, "12": 0.16, "13": 0.21, "14": 0.21, "15": 0.33, "16": 0.37, "17": 0.37, "18": 0.47, "19": 0.7, "20": 1.01, "21": 1.24, "22": 1.0, "23": 0.36, "24": 0.87, "25": 2.3},
    "ST-1-FDM (paper's model, digitised)": {"1": -0.01, "2": -0.01, "3": 0.0, "4": -0.03, "5": -0.03, "6": -0.03, "7": -0.02, "8": -0.02, "9": -0.02, "10": -0.02, "11": 0.02, "12": 0.02, "13": 0.02, "14": 0.02, "15": 0.0, "16": 0.04, "17": 0.05, "18": 0.1, "19": 0.19, "20": 0.19, "21": 0.4, "22": 0.44, "23": 0.59, "24": 1.0, "25": 1.61},
    "ST-2-HB (paper's model, digitised)": {"1": 0.07, "2": -0.1, "3": -0.13, "4": -0.3, "5": -0.35, "6": -0.51, "7": -0.75, "8": -1.15, "9": -1.38, "10": -1.53, "11": -1.64, "12": -1.87, "13": -2.07, "14": -2.28, "15": -2.42, "16": -2.72, "17": -3.06, "18": -3.53, "19": -3.91, "20": -4.15, "21": -4.49, "22": -4.87, "23": -5.13, "24": -4.72, "25": -3.68},
}
PAPER_CURVES = {k: {int(r): float(v) for r, v in d.items()} for k, d in PAPER_CURVES.items()}   # integer run keys
PAPER_STYLE = {"US-1 (inclinometer, from the raw log)": ("#e0201b", "x", "-"), "US-2 (paper figure, digitised)": ("#f0a020", "x", "-"),
               "ST-1-FDM (paper's model, digitised)": ("#3a5a9a", "x", "--"), "ST-2-HB (paper's model, digitised)": ("#2060ff", "x", "--")}

# ============================ READERS ================================
def read_hist(p):
    try:
        d = np.genfromtxt(str(p), skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2:
        return None, None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return (d[:, 0], d[:, 1]) if len(d) >= 5 else (None, None)
def find_channel(folder, key):
    hits = sorted(Path(folder).glob("*" + key + "*.csv"))
    return read_hist(hits[0]) if hits else (None, None)
def discover(d):
    out = {}
    if not Path(d).is_dir():
        return out
    for f in sorted(Path(d).iterdir()):
        m = RUN_RE.match(f.name)
        if f.is_dir() and m:
            out[int(m.group(1))] = f
    return out

def series_metrics(results_dir):
    """{channel: {run: (resid_deg, peak_deg)}} + stitched histories {channel: (t, x)}; zero = first sample of the first run."""
    folders = discover(results_dir)
    res = {k: {} for k in TILT_KEYS + ["proxy_2.06m"]}; base = {}; hist = {k: ([], []) for k in TILT_KEYS}; t_off = 0.0
    for rn in sorted(folders):
        f = folders[rn]; dur = 0.0
        for k in TILT_KEYS:
            t, x = find_channel(f, k)
            if x is None:
                continue
            base.setdefault(k, float(x[0]))
            res[k][rn] = (float(np.mean(x[int((1 - TAIL_FRAC) * len(x)):]) - base[k]), float(np.max(np.abs(x - x[0]))))
            hist[k][0].append(t - t[0] + t_off); hist[k][1].append(x - base[k]); dur = max(dur, float(t[-1] - t[0]))
        for key in (KEY_DISP_EXP, KEY_DISP):
            t, u = find_channel(f, key)
            if u is not None:
                base.setdefault("u", float(u[0]))
                r = float(np.mean(u[int((1 - TAIL_FRAC) * len(u)):]) - base["u"]); pk = float(np.max(np.abs(u - u[0])))
                res["proxy_2.06m"][rn] = (np.degrees(np.arctan(r / 1000.0 / H_EDP_M)), np.degrees(np.arctan(pk / 1000.0 / H_EDP_M)))
                dur = max(dur, float(t[-1] - t[0])); break
        t_off += dur
    hist = {k: (np.concatenate(v[0]), np.concatenate(v[1])) for k, v in hist.items() if v[0]}
    bounds = []
    t_off = 0.0
    for rn in sorted(folders):
        f = folders[rn]; dur = 0.0
        for k in TILT_KEYS + [KEY_DISP_EXP, KEY_DISP]:
            t, x = find_channel(f, k)
            if x is not None:
                dur = max(dur, float(t[-1] - t[0]))
        bounds.append((rn, t_off)); t_off += dur
    return res, hist, bounds

def read_table(path, run_col="run"):
    if not os.path.isfile(path):
        return {}
    out = {}
    for row in csv.DictReader(open(path)):
        try:
            out[int(float(row[run_col]))] = {k: float(v) for k, v in row.items() if k != run_col and v not in ("", None)}
        except Exception:
            pass
    return out

# ============================ COLLECT ================================
print("discovery:")
data = []
for lab, d, c, mk in SERIES:
    res, hist, bounds = series_metrics(d)
    n = {k: len(v) for k, v in res.items() if v}
    print("  {:34s} {}  ->  {}".format(lab, d, ", ".join("{} {}".format(k, n[k]) for k in n) or "NOT FOUND"))
    if res["tilt_full_wall"] and res["proxy_2.06m"]:
        common = [rn for rn in res["tilt_full_wall"] if rn in res["proxy_2.06m"] and abs(res["proxy_2.06m"][rn][1]) > 0.01]
        if common:
            ratio = np.median([res["tilt_full_wall"][rn][1] / res["proxy_2.06m"][rn][1] for rn in common])
            print("      tilt_full_wall / atan(u_2.06/2.06) (peaks, median over runs) = {:.2f}  ->  {}".format(
                ratio, "the channel IS the base-to-2.06 m angle" if 0.9 < ratio < 1.1 else "the channel is NOT the base-to-2.06 m angle: compare with care"))
    data.append((lab, d, c, mk, res, hist, bounds))
exp_pot = read_table(EXP_POT_CSV); exp_inc = read_table(EXP_INC_CSV)
print("  US-1 potentiometers: {}".format("{} runs from {}".format(len(exp_pot), EXP_POT_CSV) if exp_pot else "no " + EXP_POT_CSV))
print("  US-1 inclinometer  : {}".format("{} runs from {}".format(len(exp_inc), EXP_INC_CSV) if exp_inc else "no " + EXP_INC_CSV))

# ============================ CSV ====================================
with open(os.path.join(OUT_DIR, "cumulative_tilt.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["series", "channel", "run", "resid_deg", "peak_deg"])
    for lab, d, c, mk, res, hist, bounds in data:
        for k, v in res.items():
            for rn in sorted(v):
                w.writerow([lab, k, rn, round(v[rn][0], 5), round(v[rn][1], 5)])
    for rn in sorted(exp_pot):
        w.writerow(["US-1 potentiometers", "tilt_full_wall", rn, round(exp_pot[rn].get("resid_tilt_full_wall", np.nan), 5), round(exp_pot[rn].get("peak_tilt_full_wall", np.nan), 5)])
    for rn in sorted(exp_inc):
        w.writerow(["US-1 inclinometer", "tilt_y", rn, round(exp_inc[rn].get("resid_tilt_y_deg", np.nan), 5), ""])
print("-> " + os.path.join(OUT_DIR, "cumulative_tilt.csv"))

# ============================ FIG 1: per-run =========================
fig, axs = plt.subplots(2, 2, figsize=(15, 9.5), dpi=140)
def plot_series(ax, key, resid=True):
    for lab, d, c, mk, res, hist, bounds in data:
        v = res.get(key, {})
        if v:
            ks = sorted(v); ax.plot(ks, [v[rn][0 if resid else 1] for rn in ks], marker=mk, ms=5, lw=1.4, color=c, label=lab)
ax = axs[0, 0]
plot_series(ax, "tilt_full_wall")
if exp_pot:
    ks = sorted(exp_pot); ax.plot(ks, [exp_pot[rn].get("resid_tilt_full_wall", np.nan) for rn in ks], "o-", ms=5, lw=1.4, color=COL_EXP, label="US-1, potentiometers: atan(u_top/2.06) from run-1 zero")
ax.set_title("(a) cumulative full-wall tilt: model channel vs US-1 potentiometers", loc="left", fontsize=10.5, fontweight="bold")
ax = axs[0, 1]
plot_series(ax, "tilt_full_wall")
for name, vals in PAPER_CURVES.items():
    c, mk, ls = PAPER_STYLE[name]
    if name.startswith("US-1") and exp_inc:      # prefer the CSV when it exists (same numbers, from the log)
        ks = sorted(exp_inc); ys = [exp_inc[rn].get("resid_tilt_y_deg", np.nan) for rn in ks]; name = "US-1 (inclinometer, " + EXP_INC_CSV + ")"
    else:
        ks = sorted(vals); ys = [vals[rn] for rn in ks]
    ax.plot(ks, ys, marker=mk, ms=5, lw=1.4 if name.startswith("US") else 1.0, ls=ls, color=c, label=name)
ax.set_title("(b) the same model curves on the paper's tilt figure (inclinometer scale)", loc="left", fontsize=10.5, fontweight="bold")
for j, (lab, d, c, mk, res, hist, bounds) in enumerate(data[:2]):
    ax = axs[1, j]
    for k in TILT_KEYS + ["proxy_2.06m"]:
        v = res.get(k, {})
        if v:
            ks = sorted(v); ax.plot(ks, [v[rn][0] for rn in ks], marker="o" if k != "proxy_2.06m" else ".", ms=4, lw=1.3, ls="-" if k != "proxy_2.06m" else ":", color=SEG_COL[k], label=k)
    ax.set_title("({}) {}: cumulative tilt per segment".format("cd"[j], lab), loc="left", fontsize=10.5, fontweight="bold")
for ax in axs.flat:
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("run"); ax.set_ylabel("residual tilt from run-1 zero (deg)")
    ax.set_xticks(range(1, 26)); ax.tick_params(axis="x", labelsize=7.5); ax.grid(color="#eee"); ax.legend(fontsize=8, frameon=False)
    ax.axvspan(21.5, 25.5, color="#f3e6d8", alpha=0.5, lw=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.suptitle("Cumulative (residual) tilt over the sequence: model vs US-1  (residual = mean of the last 5% of each run, from the run-1 start)", fontsize=11.5, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.95)); p = os.path.join(OUT_DIR, "fig_cumtilt_sequence.png"); plt.savefig(p); plt.close(fig); print("-> " + p)

# ============================ FIG 2: stitched history ================
have = [(lab, c, hist, bounds) for lab, d, c, mk, res, hist, bounds in data if "tilt_full_wall" in hist]
if have:
    fig, axs = plt.subplots(len(have), 1, figsize=(16, 3.6 * len(have)), dpi=140, squeeze=False)
    for ax, (lab, c, hist, bounds) in zip(axs[:, 0], have):
        t, x = hist["tilt_full_wall"]; ax.plot(t, x, color=c, lw=0.6)
        for rn, t0 in bounds:
            ax.axvline(t0, color="#ccc", lw=0.6); ax.text(t0, 0.97, " {}".format(rn), fontsize=6.5, va="top", color="#8a6a4a", transform=ax.get_xaxis_transform())
        ax.axhline(0, color="k", lw=0.5); ax.set_ylabel("tilt_full_wall from run-1 zero (deg)"); ax.set_xlabel("stitched time over the sequence (s)")
        ax.set_title("{}: full-wall tilt over all runs (run numbers at the top)".format(lab), loc="left", fontsize=10, fontweight="bold"); ax.grid(color="#eee")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.suptitle("Tilt history stitched across the sequence: the lean accumulates run by run", fontsize=11.5, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.95)); p = os.path.join(OUT_DIR, "fig_cumtilt_history.png"); plt.savefig(p); plt.close(fig); print("-> " + p)
