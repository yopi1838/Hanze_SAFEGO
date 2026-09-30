# -*- coding: ascii -*-
"""
exp_tilt_from_raw.py -- reconstruct US-1's tilt per run from its own
displacement channels (raw Test9RunNN.xlsx), so tilt can be validated like
displacement and shear.

Channels (Test9_Info.xlsx): 1 = bottom quarter, z = 0.66 m; 2 = mid,
z = 1.26 m; 3 / 4 = top quarter left / right, z = 2.06 m; 5 = shake table.
Checked on run 24: channels 1-4 are already RELATIVE to the table (their
correlation with channel 5 is 0.01-0.03, and subtracting channel 5 makes
them track the 121 mm table stroke), so they are used as they are. The
raw workbook has column 0 = sample index, column 1 = time, channel k in
column k + 1.

Tilts (deg), positive = the wall leans in the +displacement direction:
    tilt_bot_seg   : atan( u1 / 0.66 )                 base   -> 0.66 m
    tilt_low_seg   : atan( (u2 - u1) / 0.60 )          0.66   -> 1.26 m
    tilt_up_seg    : atan( (u_top - u2) / 0.80 )       1.26   -> 2.06 m
    tilt_full_wall : least-squares slope of the displaced shape through
                     (0, 0), (0.66, u1), (1.26, u2), (2.06, u_top)
    tilt_top_proxy : atan( u_top / 2.06 )              (the rigid-rocking proxy)
with u_top = 0.5 (u3 + u4).
Per run: peak = max |tilt(t) - tilt(t0)| within the run;
         residual = mean of the last 5% minus the value at the start of RUN 1
         (the potentiometers keep their offsets between runs, so run 1's
         first sample is the common zero; check RESID_OK below).
NOTE the model's tilt_full_wall channel (instrument_tilt_v2.dat) must use the
same construction for the comparison to be like for like; if it is the
base-to-top-of-wall angle instead, compare against tilt_top_proxy or
recompute one of the two.

OUTPUT
    exp_Test9_tilt.csv  : run, peak_/resid_ for the five tilts, plus the
                          three peak relative displacements (mm) -- read by
                          postprocess_route2_full.py (peak_tilt_full_wall,
                          resid_tilt_full_wall) and postprocess_stratC.py
    fig_exp_tilt_segments.png : the five tilts per run, and the displaced
                          shape at the run peak for a few runs
Usage: python exp_tilt_from_raw.py   (from the working directory; raw files
       are searched in EXP_DATA/, EXP_DATA/raw/, EXP_DATA/Test9/ and here)
"""
import os, csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
TEST      = "Test9"
RUNS      = range(1, 25)
SEARCH    = ["EXP_DATA", "EXP_DATA/raw", "EXP_DATA/Test9", "EXP_DATA/RAW", "."]
Z         = {"bot": 0.66, "mid": 1.26, "top": 2.06}         # gauge heights (m)
COL_T, COL = 1, {"bot": 2, "mid": 3, "topL": 4, "topR": 5, "table": 6}   # workbook columns (0-based)
TAIL_FRAC = 0.05
SHAPE_RUNS = (12, 18, 21, 22, 23, 24)                       # displaced shapes drawn for these
OUT_CSV   = "exp_{}_tilt.csv".format(TEST)
OUT_PNG   = "fig_exp_tilt_segments.png"

def find_raw(run):
    for d in SEARCH:
        for nm in ("{}Run{:02d}.xlsx", "{}Run{}.xlsx"):
            p = Path(d) / nm.format(TEST, run)
            if p.is_file():
                return p
    return None

def read_raw(p):
    import openpyxl
    ws = openpyxl.load_workbook(str(p), read_only=True, data_only=True).worksheets[0]
    rows = ws.iter_rows(min_row=2, values_only=True)
    d = np.array([r for r in rows if r and r[0] is not None and r[COL_T] is not None], dtype=float)
    t = d[:, COL_T]
    u = {k: d[:, c] for k, c in COL.items()}
    u["top"] = 0.5 * (u["topL"] + u["topR"])
    return t, u

def tilts_deg(u):
    """dict of tilt time histories (deg) from RELATIVE displacements (m)."""
    zb, zm, zt = Z["bot"], Z["mid"], Z["top"]
    out = {"tilt_bot_seg": np.arctan(u["bot"] / zb),
           "tilt_low_seg": np.arctan((u["mid"] - u["bot"]) / (zm - zb)),
           "tilt_up_seg":  np.arctan((u["top"] - u["mid"]) / (zt - zm)),
           "tilt_top_proxy": np.arctan(u["top"] / zt)}
    zz = np.array([0.0, zb, zm, zt]); szz = float(np.sum(zz * zz))
    out["tilt_full_wall"] = np.arctan((zb * u["bot"] + zm * u["mid"] + zt * u["top"]) / szz)   # LSQ slope through the origin
    return {k: np.degrees(v) for k, v in out.items()}

KEYS = ["tilt_bot_seg", "tilt_low_seg", "tilt_up_seg", "tilt_full_wall", "tilt_top_proxy"]
rows, base, shapes = [], None, {}
print("reconstructing US-1 tilt from raw channels 1-4:")
for rn in RUNS:
    p = find_raw(rn)
    if p is None:
        print("  run {:2d}: raw workbook not found".format(rn)); continue
    t, u = read_raw(p)
    th = tilts_deg(u)
    if base is None:
        base = {k: float(v[0]) for k, v in th.items()}; base_run = rn
        if rn != 1:
            print("  ! residual baseline taken from run {} (run 1 missing)".format(rn))
    row = {"run": rn}
    for k in KEYS:
        v = th[k]; row["peak_" + k] = float(np.max(np.abs(v - v[0])))
        row["resid_" + k] = float(np.mean(v[int((1 - TAIL_FRAC) * len(v)):]) - base[k])
    for k in ("bot", "mid", "top"):
        row["peak_disp_{}_mm".format(k)] = float(np.max(np.abs(u[k] - u[k][0]))) * 1000.0
    i = int(np.argmax(np.abs(u["top"] - u["top"][0])))
    shapes[rn] = np.array([0.0, (u["bot"][i] - u["bot"][0]), (u["mid"][i] - u["mid"][0]), (u["top"][i] - u["top"][0])]) * 1000.0
    row["table_stroke_mm"] = float(np.ptp(u["table"])) * 1000.0
    rows.append(row)
    print("  run {:2d}: peak full-wall {:.3f} deg (bot {:.3f}, low {:.3f}, up {:.3f}), residual full-wall {:+.3f} deg | disp bot/mid/top {:.1f}/{:.1f}/{:.1f} mm".format(
        rn, row["peak_tilt_full_wall"], row["peak_tilt_bot_seg"], row["peak_tilt_low_seg"], row["peak_tilt_up_seg"], row["resid_tilt_full_wall"],
        row["peak_disp_bot_mm"], row["peak_disp_mid_mm"], row["peak_disp_top_mm"]))
if not rows:
    raise SystemExit("no raw Test9RunNN.xlsx found in " + ", ".join(SEARCH))

# ---- sanity check on the residual baseline: run-to-run start values should be continuous
RESID_OK = True
if len(rows) > 1:
    jumps = []
    prev = None
    for rn in RUNS:
        p = find_raw(rn)
        if p is None:
            continue
        t, u = read_raw(p); th = tilts_deg(u)
        s = float(th["tilt_full_wall"][0]); e = float(np.mean(th["tilt_full_wall"][int((1 - TAIL_FRAC) * len(t)):]))
        if prev is not None:
            jumps.append((rn, s - prev))
        prev = e
    big = [(rn, j) for rn, j in jumps if abs(j) > 0.02]
    if big:
        RESID_OK = False
        print("  ! start-of-run tilt differs from the previous run's end by > 0.02 deg at runs {} -- sensors re-zeroed? residuals from run 1 are then not valid".format(
            ", ".join("{} ({:+.3f})".format(rn, j) for rn, j in big)))

cols = ["run"] + ["{}_{}".format(a, k) for k in KEYS for a in ("peak", "resid")] + ["peak_disp_bot_mm", "peak_disp_mid_mm", "peak_disp_top_mm", "table_stroke_mm"]
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows:
        w.writerow({k: (round(r[k], 5) if isinstance(r[k], float) else r[k]) for k in cols})
print("-> " + OUT_CSV + ("" if RESID_OK else "   (residual columns flagged, see above)"))

# ---- figure
fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=140)
rr = [r["run"] for r in rows]
cl = {"tilt_bot_seg": "#d98c21", "tilt_low_seg": "#5b7f74", "tilt_up_seg": "#b85042", "tilt_full_wall": "royalblue", "tilt_top_proxy": "#888888"}
ax = axs[0]
for k in KEYS:
    ax.plot(rr, [r["peak_" + k] for r in rows], "o-" if k != "tilt_top_proxy" else "o:", ms=4, lw=1.3, color=cl[k], label=k)
ax.set_title("US-1 peak tilt per run, reconstructed from channels 1-4", loc="left", fontsize=10.5, fontweight="bold"); ax.set_ylabel("peak tilt within run (deg)")
ax = axs[1]
for k in KEYS:
    ax.plot(rr, [r["resid_" + k] for r in rows], "o-" if k != "tilt_top_proxy" else "o:", ms=4, lw=1.3, color=cl[k], label=k)
ax.axhline(0, color="k", lw=0.5)
ax.set_title("US-1 residual tilt (from run-{} start){}".format(base_run, "" if RESID_OK else "  [baseline suspect]"), loc="left", fontsize=10.5, fontweight="bold"); ax.set_ylabel("residual tilt (deg)")
for ax in axs[:2]:
    ax.set_xlabel("run"); ax.set_xticks(range(1, 26)); ax.tick_params(axis="x", labelsize=7.5); ax.grid(color="#eee"); ax.legend(fontsize=8, frameon=False)
    ax.axvspan(21.5, 25.5, color="#f3e6d8", alpha=0.5, lw=0)
ax = axs[2]
zz = [0.0, Z["bot"], Z["mid"], Z["top"]]
for rn in SHAPE_RUNS:
    if rn in shapes:
        ax.plot(shapes[rn], zz, "o-", ms=4, lw=1.3, label="run {}".format(rn))
ax.axvline(0, color="k", lw=0.5); ax.set_xlabel("relative displacement at the top-quarter peak (mm)"); ax.set_ylabel("height (m)"); ax.set_ylim(0, 2.6)
ax.axhline(2.58, color="#bbb", lw=0.8, ls=":"); ax.text(ax.get_xlim()[0], 2.58, " wall top 2.58 m", fontsize=8, va="bottom", color="#888")
ax.set_title("Displaced shape at the run peak: rigid rocking = straight line", loc="left", fontsize=10.5, fontweight="bold"); ax.grid(color="#eee"); ax.legend(fontsize=8, frameon=False)
for ax in axs:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(OUT_PNG); print("-> " + OUT_PNG)