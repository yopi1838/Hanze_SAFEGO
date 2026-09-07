# -*- coding: ascii -*-
"""
postprocess_routes.py -- Route 1 / Route 2 / Strategy-F control, run by run.

Metrics per run (same definitions as postprocess_stratC.py):
    peak OOP rel. disp  : max |rel_disp_top_mm - start|        [mm]
    peak base shear     : max |cstav|  (base + joist reaction; the top-joint
                          share is NOT in the dynamic exports, so this is
                          ~55% of the total -- consistent across cases)  [kN]
    max tilt            : max |tilt_full_wall|                 [deg]
    residual tilt       : mean of the last 5% of tilt_full_wall [deg]
Plus base shear vs displacement loops for the runs in HYST_RUNS.

Experiment overlay: exp_Test9_metrics.csv (US-1) if present -- columns are
matched loosely on 'peak'/'resid' + 'tilt'/'disp'.

Usage: python postprocess_routes.py   (edit CASES below)
"""
import os, re, glob, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
CASES = [   # (label, results dir, style)
    ("Control (Strategy F)", "stratF_results",  dict(color="#5b7f74", marker="s")),
    ("Route 2 (two-stage)",  "route2_full25",   dict(color="#b85042", marker="o")),
    ("Route 1 (velocity)",   "route1_results",  dict(color="#3f6fb5", marker="^")),
]
EXP_CSV   = "exp_Test9_metrics.csv"     # optional US-1 overlay
HYST_RUNS = (21, 24, 25)
OUT_PNG_M = "fig_routes_metrics.png"
OUT_PNG_H = "fig_routes_hysteresis.png"
OUT_CSV   = "routes_metrics.csv"
TAIL_FRAC = 0.05
KEY_DISP, KEY_TILT, KEY_SHEAR = "rel_disp_top_mm", "tilt_full_wall", "cstav"
RUN_RE = re.compile(r"Run(\d+)_|route\d_[A-Za-z0-9]+_s")   # RunNN_... or route1_FR76_s2p00_N1p0

# ============================ READERS ================================
def read_hist(path):
    try:
        d = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2 or len(d) < 5:
        return None, None
    return d[:, 0], d[:, 1]

def channel(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else (None, None)

def discover(results_dir):
    """{run_no: folder}. A route1 single-run folder is mapped to run 25."""
    out = {}
    if not os.path.isdir(results_dir):
        return out
    for d in sorted(os.listdir(results_dir)):
        p = os.path.join(results_dir, d)
        if not os.path.isdir(p):
            continue
        m = re.match(r"Run(\d+)_", d)
        if m:
            out[int(m.group(1))] = p
        elif re.match(r"route\d_[A-Za-z0-9]+_s", d):
            m2 = re.search(r"_s(\d+)p(\d+)", d)
            sc = float("{}.{}".format(m2.group(1), m2.group(2))) if m2 else 0
            out[25 if abs(sc - 2.0) < 1e-6 else 24] = p     # FR76 x2.0 = run 25
    return out

def metrics(folder):
    t, u = channel(folder, KEY_DISP); _, th = channel(folder, KEY_TILT); _, F = channel(folder, KEY_SHEAR)
    m = {}
    if u is not None:
        m["peak_disp_mm"] = float(np.max(np.abs(u - u[0])))
    if th is not None:
        m["max_tilt_deg"] = float(np.max(np.abs(th)))
        n = max(5, int(len(th) * TAIL_FRAC))
        m["resid_tilt_deg"] = float(np.mean(th[-n:]))
    if F is not None:
        m["peak_shear_kN"] = float(np.max(np.abs(F - F[0])))
    return m

def hysteresis(folder):
    _, u = channel(folder, KEY_DISP); _, F = channel(folder, KEY_SHEAR)
    if u is None or F is None:
        return None
    n = min(len(u), len(F))
    return u[:n] - u[0], F[:n] - F[0]

def load_exp(path):
    if not os.path.isfile(path):
        return None
    rows = list(csv.DictReader(open(path)))
    def col(*words):
        for k in rows[0].keys():
            kl = k.lower()
            if all(w in kl for w in words):
                return k
        return None
    keys = {"peak_disp_mm": col("peak", "disp"), "max_tilt_deg": col("peak", "tilt", "full"),
            "resid_tilt_deg": col("resid", "tilt", "full"), "peak_shear_kN": col("peak", "shear")}
    runcol = col("run")
    out = {}
    for r in rows:
        try:
            rn = int(float(r[runcol]))
        except Exception:
            continue
        out[rn] = {k: float(r[c]) for k, c in keys.items() if c and r.get(c, "") not in ("", None)}
    return out

# ============================ COLLECT ================================
table = {}   # (label, run) -> metrics
folders = {}
for label, d, _ in CASES:
    folders[label] = discover(d)
    for rn, f in folders[label].items():
        table[(label, rn)] = metrics(f)
exp = load_exp(EXP_CSV)

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "run", "peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg"])
    for (label, rn), m in sorted(table.items(), key=lambda x: (x[0][0], x[0][1])):
        w.writerow([label, rn] + [round(m[k], 4) if k in m else "" for k in
                                  ("peak_disp_mm", "peak_shear_kN", "max_tilt_deg", "resid_tilt_deg")])
print("-> " + OUT_CSV)

# ============================ FIG 1: metrics vs run ==================
panels = [("peak_disp_mm", "Peak OOP rel. displacement (mm)"),
          ("peak_shear_kN", "Peak base shear, cstav (kN)"),
          ("max_tilt_deg", "Max tilt, full wall (deg)"),
          ("resid_tilt_deg", "Residual tilt, full wall (deg)")]
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5), dpi=140)
for ax, (key, ttl) in zip(axs.ravel(), panels):
    if exp:
        xs = sorted(r for r in exp if key in exp[r])
        if xs:
            ax.plot(xs, [exp[r][key] for r in xs], "o--", color="royalblue", mfc="none",
                    ms=6, lw=1.0, label="US-1 (experiment)")
    for label, _, st in CASES:
        xs = sorted(r for (l, r) in table if l == label and key in table[(l, r)])
        if not xs:
            continue
        ys = [table[(label, r)][key] for r in xs]
        if len(xs) == 1:
            ax.plot(xs, ys, linestyle="none", ms=11, mew=2, label=label, **st)
        else:
            ax.plot(xs, ys, "-", lw=1.6, ms=5, label=label, **st)
    if key == "peak_disp_mm":
        ax.axhline(29.5, color="#444", lw=0.9, ls=":"); ax.text(0.6, 30.3, "29.5 mm", fontsize=8, color="#444")
    ax.axvline(21.5, color="#bbb", lw=0.8); ax.set_xlim(0, 26)
    ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left")
    ax.set_xlabel("run"); ax.grid(axis="y", color="#e5e5e5")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0, 0].legend(fontsize=9, frameon=False)
fig.suptitle("Route 1 / Route 2 / Control vs US-1, run by run  (base shear = base+joist reaction only)",
             fontsize=12, fontweight="bold")
plt.tight_layout(); plt.savefig(OUT_PNG_M, dpi=140); print("-> " + OUT_PNG_M)

# ============================ FIG 2: base shear vs displacement =======
fig, axs = plt.subplots(1, len(HYST_RUNS), figsize=(5.2 * len(HYST_RUNS), 4.8), dpi=140)
axs = np.atleast_1d(axs)
for ax, rn in zip(axs, HYST_RUNS):
    for label, _, st in CASES:
        f = folders[label].get(rn)
        h = hysteresis(f) if f else None
        if h is None:
            continue
        ax.plot(h[0], h[1], lw=0.9, color=st["color"], label=label, alpha=0.9)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_title("Run {}".format(rn), fontsize=11, fontweight="bold", loc="left")
    ax.set_xlabel("OOP rel. displacement (mm)"); ax.set_ylabel("base shear, cstav (kN)")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0].legend(fontsize=9, frameon=False)
plt.tight_layout(); plt.savefig(OUT_PNG_H, dpi=140); print("-> " + OUT_PNG_H)