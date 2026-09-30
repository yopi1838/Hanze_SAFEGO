# -*- coding: ascii -*-
"""
fig_paper_tilt_compare.py -- ONE figure: the paper's "Tilt Angle vs Test Run
No" (Moshfeghi et al. 2024) redrawn, with our model's readings on it.

Paper curves (embedded, no files needed):
    US-1      exact, from the raw inclinometer log (tvaly, 0.001 deg units); the
              2.2-6 deg range for the unstrengthened walls in the paper's text is
              reproduced (US-2 2.3 deg, US-1 6.04 deg), so log and figure agree
    US-2, ST-1-FDM, ST-2-HB   digitised from the published figure (+-0.05 deg)
Our model: residual LOCAL rotation at the tiltmeter height per run, = mean of
the last 5% of the run minus the first sample of run 1, from the run folders
below. Channel chosen by TILT_KEYS (see the note there).

Usage: python fig_paper_tilt_compare.py
Output: fig_paper_tilt_compare.png
"""
import os, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
OURS = [  # (label, results dir, colour, marker)
    ("3DEC model, record (Strategy F)", "stratF_full_results_US1", "#5b7f74", "s"),
    ("3DEC model, Route 2",             "route2_full25_revised",   "#b85042", "^"),
]
# WHICH MODEL CHANNEL (settled by instrument_tilt_v2.dat).
#   tilt_local_incl  chord 2.33 -> 2.53 m, centred on the tiltmeter at Z = 2.43 m (Table 3,
#                    sensor 18). The LOCAL rotation of the instrumented unit -- the counterpart
#                    of the inclinometer, and the only one of these that is like for like.
#   tilt_beam_seg    2.06 -> 2.68 m; also brackets 2.43 m, a usable coarser proxy.
#   tilt_full_wall   atan(u_2.06 / 2.06), the wall's CHORD lean. This is NOT what the tiltmeter
#                    reads: 2.77 deg over the 2.06 m chord would be 99 mm of permanent top
#                    displacement against the -1.18 mm the transducer records.
# The first of these found in the run folders is used.
# CAVEAT carried over from instrument_tilt_v2.dat: incl_y_lo / incl_y_hi (2.33 / 2.53) are
# PLACEHOLDERS. The paper gives the tiltmeter's position but not the extent of the unit it is
# bonded to, so the chord length is a guess; a wrong height gives a plausible but wrong number.
# EXPORT CHECK: tilt_local_incl is FISH history index 23. instrument_history_export_v2.dat
# exports by numeric index -- if it stops at 16, this channel is simply never written and the
# script will fall back with a message.
TILT_KEYS = ["tilt_local_incl", "tilt_beam_seg", "tilt_full_wall"]
TAIL_FRAC = 0.05
OUT_PNG   = "fig_paper_tilt_compare.png"
RUN_RE    = re.compile(r"^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$")

PAPER = {
    "US-1":     (dict(zip(range(1, 25), [0.05, 0.05, 0.07, 0.05, 0.08, 0.09, 0.06, 0.09, 0.09, 0.21, 0.23, 0.28,
                                         0.37, 0.47, 0.58, 0.74, 1.03, 1.51, 2.07, 2.39, 2.76, 2.53, 2.40, 6.04])), "#e0201b", "-"),
    "US-2":     (dict(zip(range(1, 26), [0.00, 0.10, 0.10, 0.09, 0.12, -0.05, 0.11, 0.12, 0.17, 0.26, 0.12, 0.16,
                                         0.21, 0.21, 0.33, 0.37, 0.37, 0.47, 0.70, 1.01, 1.24, 1.00, 0.36, 0.87, 2.30])), "#f0a020", "-"),
    "ST-1-FDM": (dict(zip(range(1, 26), [-0.01, -0.01, 0.00, -0.03, -0.03, -0.03, -0.02, -0.02, -0.02, -0.02, 0.02, 0.02,
                                         0.02, 0.02, 0.00, 0.04, 0.05, 0.10, 0.19, 0.19, 0.40, 0.44, 0.59, 1.00, 1.61])), "#3a5a9a", "-"),
    "ST-2-HB":  (dict(zip(range(1, 26), [0.07, -0.10, -0.13, -0.30, -0.35, -0.51, -0.75, -1.15, -1.38, -1.53, -1.64, -1.87,
                                         -2.07, -2.28, -2.42, -2.72, -3.06, -3.53, -3.91, -4.15, -4.49, -4.87, -5.13, -4.72, -3.68])), "#2060ff", "-"),
}

# ============================ OUR MODEL ==============================
def read_hist(p):
    try:
        d = np.genfromtxt(str(p), skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None

def resid_per_run(results_dir, used):
    out, base = {}, None
    if not Path(results_dir).is_dir():
        return out
    for f in sorted(Path(results_dir).iterdir()):
        m = RUN_RE.match(f.name)
        if not (f.is_dir() and m):
            continue
        hits, key = [], None
        for k in TILT_KEYS:
            hits = sorted(f.glob("*" + k + "*.csv"))
            if hits:
                key = k; break
        if not hits:
            continue
        used.add(key)
        x = read_hist(hits[0])
        if x is None:
            continue
        if base is None:
            base = float(x[0])
        out[int(m.group(1))] = float(np.mean(x[int((1 - TAIL_FRAC) * len(x)):]) - base)
    return out

print("our model:")
ours = []
for lab, d, c, mk in OURS:
    used = set(); v = resid_per_run(d, used)
    print("  {:34s} {:26s} {:2d} runs  channel {:16s}{}".format(
        lab, d, len(v), "/".join(sorted(used)) or "-", "  (last: {:+.3f} deg)".format(v[max(v)]) if v else "  -- NOT FOUND"))
    if used and "tilt_local_incl" not in used:
        print("      ! not tilt_local_incl -- the tiltmeter counterpart is missing from the export "
              "(FISH index 23; check instrument_history_export_v2.dat). This is a coarser chord.")
    ours.append((lab, v, c, mk))

# ============================ FIGURE =================================
fig, ax = plt.subplots(figsize=(11, 6.4), dpi=150)
for name, (vals, col, ls) in PAPER.items():
    ks = sorted(vals)
    ax.plot(ks, [vals[r] for r in ks], marker="x", ms=6, lw=1.6, ls=ls, color=col, label=name)
for lab, v, c, mk in ours:
    if v:
        ks = sorted(v)
        ax.plot(ks, [v[r] for r in ks], marker=mk, ms=5, lw=1.8, color=c, label=lab)
        ax.annotate("{} ends at {:+.2f} deg".format(lab, v[max(v)]), (max(v), v[max(v)]),
                    xytext=(12, -14 - 18 * ours.index((lab, v, c, mk))), textcoords="offset points",
                    fontsize=8.5, color=c, arrowprops=dict(arrowstyle="-", color=c, lw=0.8))
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(0, 30); ax.set_xticks(range(0, 31, 5)); ax.set_ylim(-6, 8)
ax.set_xlabel("Test Run No", fontsize=11, fontweight="bold"); ax.set_ylabel("Tilt Angle (deg)", fontsize=11, fontweight="bold")
ax.grid(color="#e8e8e8"); ax.legend(fontsize=9.5, frameon=False, loc="upper left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_title("Residual tilt per run: the paper's figure, with the 3DEC model on it", loc="left", fontsize=12, fontweight="bold")
fig.text(0.01, 0.01, "US-1 from the raw inclinometer log; US-2 / ST-1-FDM / ST-2-HB digitised from the published figure. "
                     "Model = residual local rotation at the tiltmeter height, same run-1 zero (instrument_tilt_v2.dat).", fontsize=8, color="#666")
plt.tight_layout(rect=(0, 0.03, 1, 1)); plt.savefig(OUT_PNG); print("-> " + OUT_PNG)
