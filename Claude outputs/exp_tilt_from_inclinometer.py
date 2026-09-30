# -*- coding: ascii -*-
"""
exp_tilt_from_inclinometer.py -- replicate the paper's "Tilt Angle vs Test
Run No" curve for US-1 from the raw inclinometer log (EXP_DATA/
US1_Tilt_values.csv) and put the model's residual tilt channels beside it.

THE LOG: one row every 5 s, ';' separated, decimal comma. Columns used:
    sensor1 tvaly / tminy / tmaxy : tilt about y, value / min / max within
                                    the 5 s, in 0.001 deg (tvaly rises from
                                    985 to 7019 over the test = 6.03 deg,
                                    the paper's run-24 value)
    sensor1 tvalx                 : tilt about x, same units
Shaking is found from the within-5-s range (tmaxy - tminy): quiet samples
have a range of ~0.03 deg, runs show 18-51 deg (the sensor CLIPS at
+-25.75 deg, so the PEAK tilt from this sensor is meaningless; only the
residual after each run is usable). Events shorter than MIN_EVENT_S or with
a range below MIN_RANGE_DEG are dropped (one 0.4 deg blip at 550 s is not a
run). The 24 remaining events are runs 1-24 in order; their durations
(15-20 s HU12, 35-40 s EC40, 40-45 s FR76) match the protocol.
Residual after run N = median of tvaly over the 30 s after the event, minus
the median before the first event, in degrees.

WHAT IT IS NOT: the wall's global lean. The potentiometers on the same wall
show an accumulated residual displacement of 1.7 mm at 2.06 m after run 24
(0.05 deg), against 6.0 deg here. The inclinometer therefore reports the
rotation of whatever it is fixed to (a local block/course, or the top beam),
or drift after saturation -- a question for the authors. The model's
channels are of the potentiometer order, so on this figure the model sits
30-100x below US-1 for reasons that are not the model's.

OUTPUT
    exp_Test9_incl_tilt.csv  : run, t_start_s, t_end_s, duration_s,
                               resid_tilt_y_deg, resid_tilt_x_deg, clipped
    fig_paper_tilt_replica.png : the paper's figure: US-1 residual tilt from
                               the log, plus the model's residual tilt
                               channels (tilt_full_wall, tilt_up_seg,
                               tilt_beam_seg if exported) from CONTROL_DIR
                               and ROUTE_DIR, referenced to their run-1 start
Usage: python exp_tilt_from_inclinometer.py
"""
import os, csv, glob, re, datetime as dt
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
LOG_CANDIDATES = ["EXP_DATA/US1_Tilt_values.csv", "US1_Tilt_values.csv", "EXP_DATA/raw/US1_Tilt_values.csv"]
UNIT_DEG       = 0.001          # log units -> degrees
MIN_RANGE_DEG  = 1.0            # a run shows a within-5-s range far above this (quiet ~0.03; a 0.4 deg blip at 550 s is not a run)
MIN_EVENT_S    = 5.0
POST_WINDOW_S  = 30.0           # residual = median over this window after the event
CLIP_DEG       = 25.7           # sensor saturation (tmaxy - tminy ~ 51.5 deg)
N_RUNS         = 24
CONTROL_DIR, CONTROL_LABEL = "stratF_full_results_US1", "record on the model (Strategy F)"
ROUTE_DIR, ROUTE_LABEL     = "route2_full25_revised", "Route 2 (full sequence)"
MODEL_TILT_KEYS = ["tilt_full_wall", "tilt_up_seg", "tilt_beam_seg"]
TAIL_FRAC = 0.05
OUT_CSV, OUT_PNG = "exp_Test9_incl_tilt.csv", "fig_paper_tilt_replica.png"
RUN_RE = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")

# ============================ LOG ====================================
log = next((p for p in LOG_CANDIDATES if os.path.isfile(p)), None)
if log is None:
    raise SystemExit("inclinometer log not found: " + ", ".join(LOG_CANDIDATES))
lines = open(log).read().strip().splitlines(); hdr = [h.strip() for h in lines[0].split(";")]
ci = {h: i - 1 for i, h in enumerate(hdr)}           # data columns start after 'time'
t, X = [], []
for l in lines[1:]:
    p = l.replace(",", ".").split(";")
    if len(p) < len(hdr):
        continue
    t.append(dt.datetime.strptime(p[0], "%d-%m-%Y_%H:%M:%S")); X.append([float(v) for v in p[1:len(hdr)]])
X = np.array(X); ts = np.array([(x - t[0]).total_seconds() for x in t])
ty, tmn, tmx = X[:, ci["sensor1 tvaly"]], X[:, ci["sensor1 tminy"]], X[:, ci["sensor1 tmaxy"]]
tx = X[:, ci["sensor1 tvalx"]]
rng = (tmx - tmn) * UNIT_DEG
print("log: {}  {} samples over {:.0f} min, quiet within-5s range {:.3f} deg".format(log, len(ts), ts[-1] / 60, np.median(rng)))

active = rng > MIN_RANGE_DEG
events, i = [], 0
while i < len(active):
    if active[i]:
        j = i
        while j + 1 < len(active) and (active[j + 1] or (j + 2 < len(active) and active[j + 2])):
            j += 1
        if ts[j] - ts[i] + 5 >= MIN_EVENT_S:
            events.append((i, j))
        i = j + 1
    else:
        i += 1
print("{} shaking events found (expected {})".format(len(events), N_RUNS))
if len(events) != N_RUNS:
    print("  ! event count differs from the protocol; check MIN_RANGE_DEG / MIN_EVENT_S -- runs are numbered in order of occurrence")
base_y = float(np.median(ty[:events[0][0]])); base_x = float(np.median(tx[:events[0][0]]))
rows = []
for k, (i, j) in enumerate(events, 1):
    post = (ts > ts[j]) & (ts <= ts[j] + POST_WINDOW_S)
    if events[k - 1:k] and k < len(events):
        post &= ts < ts[events[k][0]]
    ry = (float(np.median(ty[post])) - base_y) * UNIT_DEG if post.any() else np.nan
    rx = (float(np.median(tx[post])) - base_x) * UNIT_DEG if post.any() else np.nan
    clipped = bool(rng[i:j + 1].max() >= 2 * CLIP_DEG - 1)
    rows.append(dict(run=k, t_start_s=float(ts[i]), t_end_s=float(ts[j] + 5), duration_s=float(ts[j] + 5 - ts[i]),
                     resid_tilt_y_deg=ry, resid_tilt_x_deg=rx, peak_range_deg=float(rng[i:j + 1].max()), clipped=int(clipped)))
    print("  run {:2d}: {:5.0f}-{:5.0f} s ({:2.0f} s)  residual tilt-y {:+6.2f} deg  tilt-x {:+6.2f}  within-run range {:5.1f} deg{}".format(
        k, ts[i], ts[j] + 5, ts[j] + 5 - ts[i], ry, rx, rng[i:j + 1].max(), "  [clipped]" if clipped else ""))
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
print("-> " + OUT_CSV)

# ============================ MODEL residual tilts ==================
def read_hist(p):
    try:
        d = np.genfromtxt(str(p), skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None
def model_resid(results_dir):
    out = {k: {} for k in MODEL_TILT_KEYS}; base = {}
    if not os.path.isdir(results_dir):
        return out
    for d in sorted(Path(results_dir).iterdir()):
        m = RUN_RE.match(d.name)
        if not (d.is_dir() and m):
            continue
        rn = int(m.group(1))
        for k in MODEL_TILT_KEYS:
            hits = sorted(d.glob("*" + k + "*.csv"))
            if not hits:
                continue
            x = read_hist(hits[0])
            if x is None:
                continue
            base.setdefault(k, float(x[0]))
            out[k][rn] = float(np.mean(x[int((1 - TAIL_FRAC) * len(x)):]) - base[k])
    return out
models = [(CONTROL_LABEL, "#5b7f74", "s", model_resid(CONTROL_DIR)), (ROUTE_LABEL, "#b85042", "^", model_resid(ROUTE_DIR))]
for lab, _, _, mr in models:
    print("  {}: {}".format(lab, ", ".join("{} ({} runs)".format(k, len(v)) for k, v in mr.items() if v) or "no tilt channels found"))

# ============================ FIGURE ================================
fig, axs = plt.subplots(1, 2, figsize=(15, 5.4), dpi=140)
ax = axs[0]
ax.plot([r["run"] for r in rows], [r["resid_tilt_y_deg"] for r in rows], "x-", color="#e0201b", lw=1.6, ms=6, label="US-1, inclinometer (tilt about y), from the raw log")
ax.plot([r["run"] for r in rows], [r["resid_tilt_x_deg"] for r in rows], "x:", color="#e0201b", lw=1.0, ms=4, alpha=0.6, label="US-1, inclinometer, tilt about x")
ls = {"tilt_full_wall": "-", "tilt_up_seg": "--", "tilt_beam_seg": ":"}
for lab, c, mk, mr in models:
    for k in MODEL_TILT_KEYS:
        if mr[k]:
            ks = sorted(mr[k]); ax.plot(ks, [mr[k][r] for r in ks], marker=mk, ms=4, lw=1.2, ls=ls[k], color=c, label="{}, {}".format(lab, k))
ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("test run no"); ax.set_ylabel("residual tilt angle (deg)")
ax.set_title("Paper's figure replicated: residual tilt per run", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(fontsize=8, frameon=False)
ax = axs[1]
ax.plot(ts / 60, (ty - base_y) * UNIT_DEG, color="#e0201b", lw=1.0, label="tvaly (5 s value)")
ax.fill_between(ts / 60, (tmn - base_y) * UNIT_DEG, (tmx - base_y) * UNIT_DEG, color="#e0201b", alpha=0.15, lw=0, label="min-max within 5 s (clips at +-{:.1f} deg)".format(CLIP_DEG))
for r in rows:
    ax.axvspan(r["t_start_s"] / 60, r["t_end_s"] / 60, color="#f3e6d8", alpha=0.6, lw=0)
    ax.text(0.5 * (r["t_start_s"] + r["t_end_s"]) / 60, ax.get_ylim()[1] if False else 27, str(r["run"]), fontsize=6.5, ha="center", color="#8a6a4a")
ax.set_ylim(-30, 30); ax.set_xlabel("time in the log (min)"); ax.set_ylabel("inclinometer tilt about y (deg)")
ax.set_title("The raw log: runs found from the within-5 s range; peaks are saturated", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(fontsize=8, frameon=False, loc="lower left")
for ax in axs:
    ax.grid(color="#eee")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.suptitle("US-1 residual tilt: inclinometer (paper) vs model channels  --  NB potentiometers give 1.7 mm residual at 2.06 m after run 24 (0.05 deg)", fontsize=10.5, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.95)); plt.savefig(OUT_PNG); print("-> " + OUT_PNG)
