# -*- coding: ascii -*-
"""
postprocess_pushover.py -- OOP pushover curves of the US-1 model vs the test.

Reads any of the pushover driver outputs (format detected from the header):
    v1        : lam_g, d_ctrl_mm, F_base_kN
    v2.x / v3 : phase, a_mps2, d_ctrl_mm, F_base_kN[, F_top_kN, F_tot_kN], conv,
                converged, m_part_kg[, a_pred_mps2, settle_chunks]
    disp-ctrl : d_ctrl_mm, F_base_kN, F_appl_kN, KE_frac, maxvel_mms
and overlays the digitised Fig 13 envelope of US-1 (US1_fig13_digitised.csv,
positive and negative branches drawn on |d|, |F|) and the Table 4 peaks as
read from that file (largest |F| on each branch) -- nothing typed in.

Force column: F_tot_kN when present (base + joists + top joint = the whole
lateral reaction, comparable to Fig 13), else F_base_kN (base + joists only,
~55% of the total on this wall -- stated in the legend). Sign: plotted as
|F| against |d| so both push directions and the test share one quadrant.

Per curve it reports: peak |F| and the displacement there, initial stiffness
K1 (least-squares slope through the origin over d <= D_LINEAR_MM), secant
stiffness and force at D_SEC_MM if reached, and the periods
T = 2 pi sqrt(M_EFF / K) for both -- same law as capacity_law.py.

Outputs (OUT_DIR):
    fig_pushover_curves.png   : |F| vs |d|, all curves + envelope + Table 4 peaks
    fig_pushover_control.png  : control variable and quality per curve
                                (a vs d with unconverged points hollow; or
                                F_appl/F_base and KE_frac for disp-control)
    pushover_summary.csv
Usage: python postprocess_pushover.py
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
CURVES = [   # (label, csv path, colour)
    ("force control (tilted gravity)",        "pushover_results/pushover_pos.csv",      "#b85042"),
    ("displacement control (top beam driven)", "pushover_results/pushover_disp_pos.csv", "#3f6fb5"),
    # ("force control, negative",             "pushover_results/pushover_neg.csv",      "#d98c21"),
]
ENVELOPE_CSV = "US1_fig13_digitised.csv"   # drift_pct, F_kN, d_mm (paper sign)
M_EFF        = 1635.0     # kg, effective mass used for the periods (same as Route 2)
D_LINEAR_MM  = 2.0        # K1 fitted over 0 < d <= D_LINEAR_MM
D_SEC_MM     = 29.5       # secant evaluated here (US-1 run 24)
OUT_DIR      = "pushover_postproc"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================ READERS ================================
def read_curve(path):
    """-> dict(d, F, Fsrc, kind, extra...) with d, F as positive arrays, or None"""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    hdr = list(rows[0].keys())
    def col(name, default=np.nan):
        out = []
        for r in rows:
            v = r.get(name, "")
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                out.append(default)
        return np.array(out)
    d = np.abs(col("d_ctrl_mm"))
    if "F_tot_kN" in hdr and np.isfinite(col("F_tot_kN")).any():
        F, Fsrc = np.abs(col("F_tot_kN")), "F_tot (base+joists+top)"
    else:
        F, Fsrc = np.abs(col("F_base_kN")), "F_base (base+joists only)"
    out = {"d": d, "F": F, "Fsrc": Fsrc, "n": len(rows)}
    if "phase" in hdr:
        out["kind"] = "load"
        out["phase"] = np.array([r.get("phase", "") for r in rows])
        out["a"] = np.abs(col("a_mps2"))
        out["conv_ok"] = col("converged", 1.0) >= 0.5
        out["m_part"] = col("m_part_kg")
    elif "F_appl_kN" in hdr:
        out["kind"] = "disp"
        out["F_appl"] = np.abs(col("F_appl_kN")); out["KE_frac"] = col("KE_frac")
    elif "lam_g" in hdr:
        out["kind"] = "load"; out["a"] = np.abs(col("lam_g")) * 9.81
        out["phase"] = np.array(["asc"] * len(rows)); out["conv_ok"] = np.ones(len(rows), bool)
    else:
        out["kind"] = "unknown"
    return out

def read_envelope(path):
    if not os.path.isfile(path):
        print("  ! {} not found -- no test envelope".format(path)); return None
    e = np.genfromtxt(path, delimiter=",", names=True)
    if "d_mm" not in e.dtype.names or "F_kN" not in e.dtype.names:
        print("  ! {} lacks d_mm / F_kN".format(path)); return None
    br = {}
    for nm, m in (("test, positive branch", e["d_mm"] > 0), ("test, negative branch", e["d_mm"] < 0)):
        d = np.abs(e["d_mm"][m]); F = np.abs(e["F_kN"][m]); o = np.argsort(d)
        br[nm] = (d[o], F[o])
    return br

# ============================ METRICS ================================
def k1_fit(d, F):
    m = (d > 0) & (d <= D_LINEAR_MM)
    if m.sum() < 2:
        m = np.argsort(d)[:3]
    return float(np.sum(d[m] * F[m]) / np.sum(d[m] * d[m])) if np.sum(d[m] * d[m]) > 0 else np.nan

def secant(d, F, dq):
    o = np.argsort(d); d, F = d[o], F[o]
    if dq > d.max():
        return np.nan, np.nan
    Fq = float(np.interp(dq, d, F)); return Fq, Fq / dq

def period(K_kNmm):
    return 2 * np.pi * np.sqrt(M_EFF / (K_kNmm * 1e6)) if K_kNmm and K_kNmm > 0 else np.nan

def metrics(name, d, F):
    i = int(np.argmax(F))
    K1 = k1_fit(d, F); Fs, Ks = secant(d, F, D_SEC_MM)
    return {"curve": name, "peak_F_kN": float(F[i]), "d_at_peak_mm": float(d[i]), "d_max_mm": float(d.max()),
            "K1_kNmm": K1, "T1_s": period(K1), "F_at_{:g}mm_kN".format(D_SEC_MM): Fs,
            "Ksec_{:g}mm_kNmm".format(D_SEC_MM): Ks, "Teff_{:g}mm_s".format(D_SEC_MM): period(Ks)}

# ============================ COLLECT ================================
print("discovery:")
curves = []
for label, path, colr in CURVES:
    c = read_curve(path)
    if c is None:
        print("  {:40s} NOT FOUND [{}]".format(label, path)); continue
    print("  {:40s} {} rows, kind {}, force = {}, d_max {:.1f} mm, peak {:.2f} kN".format(
        label, c["n"], c["kind"], c["Fsrc"], c["d"].max(), c["F"].max()))
    curves.append((label, colr, c))
env = read_envelope(ENVELOPE_CSV)

rowsout = []
for label, colr, c in curves:
    rowsout.append(metrics(label, c["d"], c["F"]) | {"force": c["Fsrc"]})
if env:
    for nm, (d, F) in env.items():
        rowsout.append(metrics(nm, d, F) | {"force": "Fig 13 (total)"})
keys = ["curve", "force"] + [k for k in rowsout[0] if k not in ("curve", "force")] if rowsout else []
with open(os.path.join(OUT_DIR, "pushover_summary.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(keys)
    for r in rowsout:
        w.writerow([r.get(k, "") if not isinstance(r.get(k), float) else round(r[k], 4) for k in keys])
print("-> " + os.path.join(OUT_DIR, "pushover_summary.csv"))
for r in rowsout:
    print("  {:40s} peak {:6.2f} kN at {:5.1f} mm | K1 {:.2f} kN/mm -> T1 {:.4f} s | at {:g} mm: F {} K {} -> T {}".format(
        r["curve"], r["peak_F_kN"], r["d_at_peak_mm"], r["K1_kNmm"], r["T1_s"], D_SEC_MM,
        "{:.1f}".format(r["F_at_{:g}mm_kN".format(D_SEC_MM)]) if np.isfinite(r["F_at_{:g}mm_kN".format(D_SEC_MM)]) else "n/a",
        "{:.3f}".format(r["Ksec_{:g}mm_kNmm".format(D_SEC_MM)]) if np.isfinite(r["Ksec_{:g}mm_kNmm".format(D_SEC_MM)]) else "n/a",
        "{:.3f}".format(r["Teff_{:g}mm_s".format(D_SEC_MM)]) if np.isfinite(r["Teff_{:g}mm_s".format(D_SEC_MM)]) else "n/a"))

# ============================ FIG 1: curves =========================
fig, ax = plt.subplots(figsize=(9.5, 5.8), dpi=150)
if env:
    for (nm, (d, F)), ls in zip(env.items(), ("-", "--")):
        ax.plot(d, F, ls, color="#444444", lw=1.8, label=nm.replace("test, ", "US-1 test, "))
        i = int(np.argmax(F)); ax.plot(d[i], F[i], "o", color="#444444", ms=5)
        ax.annotate("{:.1f} kN".format(F[i]), (d[i], F[i]), xytext=(4, 5), textcoords="offset points", fontsize=8, color="#444444")
for label, colr, c in curves:
    d, F = c["d"], c["F"]
    short = "total reaction" if c["Fsrc"].startswith("F_tot") else "base + joists only"
    if c["kind"] == "load":
        for ph, mk, ls in (("asc", "o", "-"), ("desc", "s", "-")):
            m = c["phase"] == ph
            if not m.any():
                continue
            ok = m & c["conv_ok"]; bad = m & ~c["conv_ok"]
            ax.plot(d[ok], F[ok], ls, color=colr, lw=1.4, marker=mk, ms=3.5,
                    label="model, {} ({})".format(label, short) if ph == "asc" else "model, descending (servo-held)")
            if bad.any():
                ax.plot(d[bad], F[bad], mk, mfc="white", mec=colr, ms=5, label="not converged" if ph == "asc" else None)
    else:
        ax.plot(d, F, "-", color=colr, lw=1.4, label="model, {} ({})".format(label, short))
ax.axvline(D_SEC_MM, color="#bbbbbb", lw=0.8, ls=":"); ax.text(D_SEC_MM + 0.4, 0.8, "{:g} mm".format(D_SEC_MM), fontsize=8, color="#888888")
ax.set_xlabel("|displacement| at 2.06 m (mm)"); ax.set_ylabel("|base shear| (kN)")
ax.set_title("OOP pushover of the model vs the US-1 envelope (Fig 13)", loc="left", fontweight="bold")
ax.grid(color="#eeeeee"); ax.set_xlim(left=0); ax.set_ylim(0, 36)
# legend below the axes, two columns -- keeps the curves clear
ax.legend(fontsize=8.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "fig_pushover_curves.png"), bbox_inches="tight"); print("-> fig_pushover_curves.png")

# ============================ FIG 2: control + quality ==============
if curves:
    fig, axs = plt.subplots(1, len(curves), figsize=(5.4 * len(curves), 4.4), dpi=150, squeeze=False)
    for ax, (label, colr, c) in zip(axs[0], curves):
        d = c["d"]
        if c["kind"] == "load":
            ok = c["conv_ok"]
            ax.plot(d[ok], c["a"][ok] / 9.81, "o-", color=colr, ms=4, lw=1, label="converged / settled")
            if (~ok).any():
                ax.plot(d[~ok], c["a"][~ok] / 9.81, "o", mfc="white", mec=colr, ms=5, label="unconverged")
            ax.set_ylabel("applied lateral acceleration (g)")
            if "m_part" in c and np.isfinite(c["m_part"]).any():
                ax2 = ax.twinx(); ax2.plot(d, c["m_part"], ":", color="#888888", lw=1); ax2.set_ylabel("participating mass (kg)", color="#888888")
        elif c["kind"] == "disp":
            ax.plot(d, c["F"], "-", color=colr, lw=1.2, label="F_base")
            ax.plot(d, c["F_appl"], "--", color=colr, lw=1.0, label="F_appl (top beam)")
            ax.set_ylabel("force (kN)")
            ax2 = ax.twinx(); ax2.plot(d, c["KE_frac"], ":", color="#888888", lw=1); ax2.set_ylabel("KE_frac (quasi-static check)", color="#888888")
        ax.set_xlabel("|displacement| at 2.06 m (mm)"); ax.set_title(label, loc="left", fontsize=10, fontweight="bold")
        ax.grid(color="#eeeeee"); ax.legend(fontsize=8, frameon=False)
        for sp in ("top",): ax.spines[sp].set_visible(False)
    plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "fig_pushover_control.png")); print("-> fig_pushover_control.png")
