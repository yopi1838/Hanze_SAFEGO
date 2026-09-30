# -*- coding: ascii -*-
"""
postprocess_focus_run24_25.py -- ONE case, told in full: the wall at the end
of run 24, then run 25 as the record (Strategy F) or as the revised Route 2
pulse, with US-1 run 24 (29.5 mm) as the target.

Three players
    record   : stratF_full_results_US1/Run25_*   (FR76 x2.0 on the run-24 state)
    Route 2  : route2rev_from24_bil/route2_FR76_s2p00   (the pulse on the SAME state)
    US-1     : EXP_DATA/processed_globalzero/Test9Run24_processed_globalzero.xlsx
               (run 24, FR76 x1.75 -- the last run the wall received)

Figure 1  fig_focus_inputs.png   "what the table did"
    (a) record velocity, full length, strongest 1 s shaded
    (b) that 1 s of the record next to the 0.35 s pulse, same scale
    (c) 5%-damped Sd of record x2.0, record x1.75 (US-1's run 24 input), and
        the pulse, with T1 / T_eff marked: matched at two periods, not elsewhere
Figure 2  fig_focus_response.png "what the wall did"
    (d) displacement histories aligned at the compared POSITIVE peak (the
        first large one, or the one nearest ALIGN_T_HINT), with the trough
        that follows it marked -- the same shape as the pulse's push-then-pull;
        free-vibration period after the peak from zero crossings
    (e) base shear vs displacement, one panel each, Fig 13 envelope faint
    (f) scorecard table: peak disp, peak shear, residual disp, ring-down period,
        table stroke, ratio to US-1
CSV   focus_run24_25_metrics.csv

Conventions as the other postprocess scripts: EDP = rel_disp_top_exp_mm (else
rel_disp_top_mm, else 0.5(Ch3+Ch4)-Ch5); model shear = cstav (+topj_shear),
sign flipped (MODEL_SHEAR_SIGN); US-1 shear = sum(m a) flipped (EXP_SHEAR_SIGN).
Usage: python postprocess_focus_run24_25.py   (from the working directory)
"""
import os, csv, glob, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG =================================
CONTROL_DIR  = "stratF_full_results_US1"
CONTROL_RUN  = 25
R2_FOLDER    = "route2rev_from24_bil/route2_FR76_s2p00"     # channel exports of the pulse run
R2_VEL       = R2_FOLDER + "_vel.txt"                        # the applied table velocity
R2_SUMMARY   = "route2rev_from24_bil/route2_revised_summary.csv"
REC_TABLE    = "vel_FR.txt"                                  # driver's FR76 table (x scale)
REC_SCALE    = 2.0                                           # run 25
EXP_SCALE    = 1.75                                          # US-1 run 24
SPECTRUM     = "spectrum_FR76.csv"
FIG13_CSV    = "US1_fig13_digitised.csv"
BILINEAR_JSON, BILINEAR_KEY = "bilinear_idealisation/bilinear_periods.json", "model"
EXP_ROOT, EXP_GZ_DIR, EXP_TEST, EXP_RUN = "EXP_DATA", "processed_globalzero", "Test9", 24
EXP_SHEAR_SIGN, MODEL_SHEAR_SIGN, EXP_FIG13_SIGN = -1.0, -1.0, -1.0
OUT_DIR = "focus_run24_25"; os.makedirs(OUT_DIR, exist_ok=True)
XI, DT = 0.05, 0.005
WIN_S = 1.0            # length of the record window shown next to the pulse
ALIGN_BEFORE, ALIGN_AFTER = 0.4, 1.6   # seconds shown before / after the positive peak in (d)
# The swing that is compared: the FIRST large POSITIVE peak (>= SWING_FRAC x the largest positive value)
# and the negative trough that follows it -- the same shape as the pulse (positive push, then negative).
# ALIGN_T_HINT overrides the choice per case: the positive local maximum nearest to that time is used.
SWING_FRAC   = 0.5
ALIGN_T_HINT = {"rec": 12.5, "r2": None, "exp": None}     # seconds on each case's own time axis
KEY_DISP_EXP, KEY_DISP, KEY_TILT, KEY_SHEAR, KEY_SHEAR_TOP = \
    "rel_disp_top_exp_mm", "rel_disp_top_mm", "tilt_full_wall", "cstav", "topj_shear"
CH3, CH4, CH5 = "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
C = {"rec": "#5b7f74", "r2": "#b85042", "exp": "royalblue", "A": "#5b7f74", "B": "#b85042"}
LBL = {"rec": "record, run 25 (FR76 x2.0)", "r2": "Route 2 pulse (x2.0)", "exp": "US-1, run 24 (FR76 x1.75)"}

# ============================ READERS ================================
def read_hist_csv(p):
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
    return read_hist_csv(hits[0]) if hits else (None, None)
def edp_series(folder):
    for key in (KEY_DISP_EXP, KEY_DISP):
        t, u = find_channel(folder, key)
        if u is not None:
            return t, u, key
    t3, c3 = find_channel(folder, CH3); _, c4 = find_channel(folder, CH4); _, c5 = find_channel(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, None, None
    n = min(len(c3), len(c4), len(c5))
    return t3[:n], (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0, "0.5(Ch3+Ch4)-Ch5"
def shear_series(folder):
    t, F = find_channel(folder, KEY_SHEAR)
    if F is None:
        return None, None, None
    _, Ft = find_channel(folder, KEY_SHEAR_TOP)
    if Ft is not None:
        n = min(len(F), len(Ft)); return t[:n], MODEL_SHEAR_SIGN * (F[:n] + Ft[:n]), "cstav+topj_shear (total)"
    return t, MODEL_SHEAR_SIGN * F, "cstav (base+joist only)"
def model_case(folder):
    out = {"folder": folder}
    t, u, src = edp_series(folder)
    if u is None:
        return None
    out.update(t=t, u=u - u[0], edp_src=src)
    t, F, ssrc = shear_series(folder)
    if F is not None:
        out.update(tF=t, F=F - F[0], shear_src=ssrc)
    t, th = find_channel(folder, KEY_TILT)
    if th is not None:
        out.update(tilt=th - th[0])
    t, c5 = find_channel(folder, CH5)
    if c5 is not None:
        out.update(table_mm=(c5 - c5[0]) * 1000.0)
    return out
def control_folder(run_no):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    return hits[0] if hits else None
def exp_case():
    for root in (Path(EXP_ROOT), Path.cwd() / EXP_ROOT, Path(__file__).resolve().parent / EXP_ROOT):
        for nm in ("{}Run{:02d}_processed_globalzero.xlsx", "{}Run{}_processed_globalzero.xlsx"):
            p = root / EXP_GZ_DIR / nm.format(EXP_TEST, EXP_RUN)
            if p.is_file():
                import openpyxl
                ws = openpyxl.load_workbook(str(p), read_only=True, data_only=True).worksheets[0]
                rows = ws.iter_rows(values_only=True); hdr = [str(h) for h in next(rows)]
                d = np.array([r for r in rows if r and r[0] is not None], dtype=float)
                col = {h: d[:, j] for j, h in enumerate(hdr)}
                if "U_avg" not in col or "base_shear_kN" not in col:
                    print("  ! {} lacks U_avg / base_shear_kN".format(p.name)); return None
                u = col["U_avg"] * 1000.0; F = col["base_shear_kN"]
                t = col["Time"] if "Time" in col else np.arange(len(u)) * (col["Time"][1] - col["Time"][0] if "Time" in col else 0.005)
                print("  experiment: {}".format(p))
                return dict(t=t - t[0], u=u - u[0], tF=t - t[0], F=EXP_SHEAR_SIGN * (F - F[0]), edp_src="U_avg (globalzero)", shear_src="sum(m a)")
    print("  ! US-1 run 24 workbook not found"); return None
def read_vel(p):
    d = np.loadtxt(p, skiprows=2); return d[:, 0], d[:, 1]

# ============================ ANALYSIS HELPERS ======================
def sd_spectrum(ag, dt, Ts, xi=0.05):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2 * np.pi / Ts; c = 2 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0] - c * v - w2 * u; umax = np.abs(u)
    for i in range(1, len(ag)):
        un = u + dt * v + 0.25 * dt * dt * a; vn = v + 0.5 * dt * a
        an = (-ag[i] - c * vn - w2 * un) / (1 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = un + 0.25 * dt * dt * an; v = vn + 0.5 * dt * an; a = an; umax = np.maximum(umax, np.abs(u))
    return umax
def strongest_window(t, v, win):
    """start time of the win-long window with the largest integral of a^2 (velocity-derived)."""
    dt = float(np.median(np.diff(t))); a2 = np.gradient(v, dt) ** 2; n = int(round(win / dt))
    cs = np.concatenate([[0], np.cumsum(a2)]); e = cs[n:] - cs[:-n]
    i = int(np.argmax(e)); return t[i], t[i + n]
def ringdown_period(t, u, t_peak, n_cross=6):
    """free-vibration period after the peak: 2 x mean zero-crossing interval of (u - tail mean)."""
    m = t > t_peak; tt, uu = t[m], u[m]
    if len(uu) < 10:
        return np.nan
    uu = uu - np.mean(uu[int(0.9 * len(uu)):])
    s = np.sign(uu); idx = np.where(s[1:] * s[:-1] < 0)[0]
    if len(idx) < 3:
        return np.nan
    tc = tt[idx[:n_cross + 1]]
    return float(2 * np.mean(np.diff(tc)))
def swing_pair(t, u, t_hint=None, frac=SWING_FRAC):
    """(i_pos, i_neg): the positive peak compared and the negative trough after it."""
    lm = np.where((u[1:-1] > u[:-2]) & (u[1:-1] >= u[2:]) & (u[1:-1] > 0))[0] + 1
    if len(lm) == 0:
        i = int(np.argmax(u)); lm = np.array([i])
    if t_hint is not None:
        i_pos = int(lm[np.argmin(np.abs(t[lm] - t_hint))])
    else:
        big = lm[u[lm] >= frac * u.max()]; i_pos = int(big[0]) if len(big) else int(lm[np.argmax(u[lm])])
    after = np.arange(i_pos + 1, len(u))
    neg = after[u[after] < 0]
    if len(neg) == 0:
        return i_pos, None
    j0 = int(neg[0]); back = np.where(u[j0:] >= 0)[0]
    j1 = j0 + int(back[0]) if len(back) else len(u)
    i_neg = j0 + int(np.argmin(u[j0:j1]))
    return i_pos, i_neg
def residual(u, frac=0.05):
    return float(np.mean(u[int((1 - frac) * len(u)):]))

# ============================ COLLECT ================================
print("discovery:")
cases = {}
cf = control_folder(CONTROL_RUN)
cases["rec"] = model_case(cf) if cf else None
print("  record run {}: {}".format(CONTROL_RUN, cf or "NOT FOUND"))
cases["r2"] = model_case(R2_FOLDER) if os.path.isdir(R2_FOLDER) else None
print("  Route 2: {}".format(R2_FOLDER if cases["r2"] else "NOT FOUND [{}]".format(R2_FOLDER)))
cases["exp"] = exp_case()
for k, c in cases.items():
    if c:
        ip, ineg = swing_pair(c["t"], c["u"], ALIGN_T_HINT.get(k))
        tp = float(c["t"][ip]); up = float(c["u"][ip]); un = float(c["u"][ineg]) if ineg is not None else np.nan
        c.update(i_pos=ip, i_neg=ineg, t_peak=tp, u_pos=up, u_neg=un, t_neg=float(c["t"][ineg]) if ineg is not None else np.nan,
                 swing_mm=up - un if np.isfinite(un) else np.nan, absmax_mm=float(np.max(np.abs(c["u"]))),
                 peak_mm=max(abs(up), abs(un)) if np.isfinite(un) else abs(up), resid_mm=residual(c["u"]),
                 T_ring=ringdown_period(c["t"], c["u"], tp),
                 peak_kN=float(np.max(np.abs(c["F"]))) if "F" in c else np.nan,
                 stroke_mm=float(c["table_mm"].max() - c["table_mm"].min()) if "table_mm" in c else np.nan)
        print("    {:8s} +peak {:+.2f} mm at {:.2f} s, following trough {:+.2f} mm at {:.2f} s (swing {:.1f} mm; |max| over the run {:.1f} mm), "
              "ring-down period {:.3f} s, residual {:.2f} mm, shear {:.1f} kN ({})".format(
            k, up, tp, un, c["t_neg"], c["swing_mm"], c["absmax_mm"], c["T_ring"], c["resid_mm"], c["peak_kN"], c.get("shear_src", "-")))
if cases["rec"] is None and cases["r2"] is None:
    raise SystemExit("neither the control run nor the Route 2 folder was found")

# inputs
t_rec, v_rec = read_vel(REC_TABLE); v25 = REC_SCALE * v_rec; v24 = EXP_SCALE * v_rec
t_r2, v_r2 = read_vel(R2_VEL) if os.path.isfile(R2_VEL) else (None, None)
if t_r2 is not None and cases["r2"] is not None and np.isnan(cases["r2"]["stroke_mm"]):
    cases["r2"]["stroke_mm"] = float(np.ptp(np.concatenate([[0], np.cumsum(0.5 * (v_r2[1:] + v_r2[:-1]) * np.diff(t_r2))]))) * 1000
if cases["rec"] is not None and np.isnan(cases["rec"]["stroke_mm"]):
    cases["rec"]["stroke_mm"] = float(np.ptp(np.concatenate([[0], np.cumsum(0.5 * (v25[1:] + v25[:-1]) * np.diff(t_rec))]))) * 1000
stage = {}
if os.path.isfile(R2_SUMMARY):
    rows = list(csv.DictReader(open(R2_SUMMARY)))
    if rows:
        stage = rows[-1]
T1 = float(stage.get("T_A_s", 0)) or None; TE = float(stage.get("T_B_s", 0)) or None
if (T1 is None or TE is None) and os.path.isfile(BILINEAR_JSON):
    b = json.load(open(BILINEAR_JSON)).get(BILINEAR_KEY, {}); T1 = T1 or b.get("T1_s"); TE = TE or b.get("Teff_s")
Ts = np.arange(0.02, 0.8001, 0.005)
dtr = float(np.median(np.diff(t_rec)))
Sd25 = sd_spectrum(np.gradient(v25, dtr), dtr, Ts, XI) * 1000; Sd24 = sd_spectrum(np.gradient(v24, dtr), dtr, Ts, XI) * 1000
Sd_r2 = sd_spectrum(np.gradient(v_r2, DT), DT, Ts, XI) * 1000 if t_r2 is not None else None
fig13 = None
if os.path.isfile(FIG13_CSV):
    d = np.genfromtxt(FIG13_CSV, delimiter=",", names=True); fig13 = (EXP_FIG13_SIGN * d["d_mm"], EXP_FIG13_SIGN * d["F_kN"])

# ============================ FIGURE 1: inputs ======================
fig, axs = plt.subplots(1, 3, figsize=(18, 4.9), dpi=150, gridspec_kw={"width_ratios": [1.5, 1.1, 1.1]})
ax = axs[0]
w0, w1 = strongest_window(t_rec, v25, WIN_S)
ax.plot(t_rec, v25, color=C["rec"], lw=0.8); ax.axvspan(w0, w1, color="#f3e6d8", alpha=0.8, lw=0)
ax.text(w1 + 0.2, -0.95 * abs(v25).max(), "strongest {:.0f} s\n(shown in b)".format(WIN_S), fontsize=9, va="bottom", color="#8a6a4a")
ax.set_title("(a) the record: FR76 x{:g}, {:.1f} s long, PGV {:.2f} m/s, PGA {:.2f} g".format(REC_SCALE, t_rec[-1], abs(v25).max(), abs(np.gradient(v25, dtr)).max() / 9.81), loc="left", fontsize=10.5, fontweight="bold")
ax.set_xlabel("time (s)"); ax.set_ylabel("table velocity (m/s)")
ax = axs[1]
m = (t_rec >= w0) & (t_rec <= w1)
ax.plot(t_rec[m] - w0, v25[m], color=C["rec"], lw=1.3, label="record, its strongest {:.0f} s".format(WIN_S))
if t_r2 is not None:
    mA = t_r2 <= (T1 or 0.11) + 1e-6
    ax.fill_between(t_r2[mA], 0, v_r2[mA], color=C["A"], alpha=0.35); ax.fill_between(t_r2[~mA | (t_r2 == t_r2[mA][-1])], 0, v_r2[~mA | (t_r2 == t_r2[mA][-1])], color=C["B"], alpha=0.3)
    ax.plot(t_r2, v_r2, color="k", lw=2.0, label="Route 2 pulse, {:.2f} s in total".format(t_r2[-1]))
    ax.text(t_r2[-1] + 0.02, abs(v_r2).max() * 0.95, "pulse: PGV {:.2f} m/s\nPGA {:.2f} g, stroke {:.0f} mm".format(abs(v_r2).max(), abs(np.gradient(v_r2, DT)).max() / 9.81, cases["r2"]["stroke_mm"] if cases["r2"] else np.nan), fontsize=8.5, va="top")
ax.axhline(0, color="k", lw=0.5); ax.set_xlim(-0.02, WIN_S); ax.set_xlabel("time within the window (s)"); ax.set_ylabel("table velocity (m/s)")
ax.set_title("(b) same time scale: many cycles vs one push", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(fontsize=8.5, frameon=False, loc="lower right")
ax = axs[2]
ax.plot(Ts, Sd25, color=C["rec"], lw=1.6, label="record x2.0 (run 25 input)")
ax.plot(Ts, Sd24, color=C["exp"], lw=1.2, ls="--", label="record x1.75 (US-1 run 24 input)")
if Sd_r2 is not None:
    ax.plot(Ts, Sd_r2, color="k", lw=1.6, label="Route 2 pulse")
note = []
for T, c, lab in ((T1, C["A"], "T1"), (TE, C["B"], "T_eff")):
    if T:
        ax.axvline(T, color=c, lw=0.9, ls=":"); ax.plot(T, np.interp(T, Ts, Sd25), "o", color=c, ms=8, mec="white")
        if Sd_r2 is not None:
            ax.plot(T, np.interp(T, Ts, Sd_r2), "s", color="k", ms=6, mec="white")
        ax.text(T, -4, "{}\n{:.3f} s".format(lab, T), fontsize=8.5, color=c, ha="center", va="top")
        note.append("at {}: record {:.1f} mm, pulse {}".format(lab, np.interp(T, Ts, Sd25), "{:.1f} mm".format(np.interp(T, Ts, Sd_r2)) if Sd_r2 is not None else "-"))
if note:
    ax.text(0.98, 0.04, "\n".join(note), transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, bbox=dict(fc="white", ec="#ddd"))
ax.set_xlim(0, 0.8); ax.set_ylim(-16, None); ax.set_xlabel("period (s)"); ax.set_ylabel("Sd, 5% damped (mm)")
ax.set_title("(c) demand by period: equal at T_eff, not elsewhere", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(fontsize=8.5, frameon=False, loc="upper left")
for ax in axs:
    ax.grid(color="#eee")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.suptitle("Run 24 -> 25 on the same damaged wall: what the table did", fontsize=12.5, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.94)); p1 = os.path.join(OUT_DIR, "fig_focus_inputs.png"); plt.savefig(p1); plt.close(fig); print("-> " + p1)

# ============================ FIGURE 2: response ====================
fig = plt.figure(figsize=(18, 10), dpi=150)
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.38, wspace=0.28)
ax = fig.add_subplot(gs[0, :2])
for k in ("exp", "rec", "r2"):
    c = cases.get(k)
    if not c:
        continue
    tt = c["t"] - c["t_peak"]; m = (tt >= -ALIGN_BEFORE) & (tt <= ALIGN_AFTER)
    ax.plot(tt[m], c["u"][m], color=C[k], lw=1.6 if k != "exp" else 1.2,
            label="{}: {:+.1f} then {:+.1f} mm, period after {:.2f} s".format(LBL[k], c["u_pos"], c["u_neg"], c["T_ring"]))
    ax.plot(0, c["u_pos"], "o", color=C[k], ms=7, mec="white"); ax.text(0.02, c["u_pos"], "{:+.1f}".format(c["u_pos"]), color=C[k], fontsize=9, va="bottom")
    if c["i_neg"] is not None:
        tn = c["t_neg"] - c["t_peak"]; ax.plot(tn, c["u_neg"], "v", color=C[k], ms=7, mec="white"); ax.text(tn + 0.02, c["u_neg"], "{:+.1f}".format(c["u_neg"]), color=C[k], fontsize=9, va="top")
ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="#999", lw=0.7, ls=":")
ax.set_xlim(-ALIGN_BEFORE, ALIGN_AFTER); ax.set_xlabel("time from the positive peak (s)"); ax.set_ylabel("OOP displacement at 2.06 m (mm)")
ax.set_title("(d) the swing: first large positive peak and the trough after it, aligned at the positive peak", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(fontsize=9, frameon=False, loc="upper right")
# scorecard
ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
rows_tbl = [("positive peak (mm)", "u_pos", "{:+.1f}"), ("trough after it (mm)", "u_neg", "{:+.1f}"), ("peak-to-peak swing (mm)", "swing_mm", "{:.1f}"),
            ("swing / US-1 swing", None, None), ("largest |u| in the run (mm)", "absmax_mm", "{:.1f}"), ("peak base shear (kN)", "peak_kN", "{:.1f}"),
            ("residual displ. (mm)", "resid_mm", "{:+.1f}"), ("period after peak (s)", "T_ring", "{:.2f}"), ("table stroke (mm)", "stroke_mm", "{:.0f}")]
cols = [k for k in ("exp", "rec", "r2") if cases.get(k)]
cell = []
for name, key, fmt in rows_tbl:
    r = [name]
    for k in cols:
        c = cases[k]
        if key is None:
            r.append("{:.2f}".format(c["swing_mm"] / cases["exp"]["swing_mm"]) if cases.get("exp") and np.isfinite(c["swing_mm"]) else "-")
        else:
            v = c.get(key, np.nan); r.append(fmt.format(v) if np.isfinite(v) else "-")
    cell.append(r)
tb = ax.table(cellText=cell, colLabels=["", *["US-1\nrun 24" if k == "exp" else "record\nrun 25" if k == "rec" else "Route 2\nx2.0" for k in cols]],
              loc="center", cellLoc="center", colLoc="center", colWidths=[0.46] + [0.18] * len(cols))
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1.0, 1.55)
for (r, cidx), cl in tb.get_celld().items():
    cl.set_edgecolor("#dddddd")
    if r == 0:
        cl.set_facecolor("#f4f4f4"); cl.set_text_props(fontweight="bold")
    if cidx == 0:
        cl.set_text_props(ha="left")
ax.set_title("(f) scorecard  (model shear = {})".format(next((cases[k]["shear_src"] for k in ("rec", "r2") if cases.get(k) and "shear_src" in cases[k]), "cstav")), loc="left", fontsize=10.5, fontweight="bold")
# loops
lim_u = max([cases[k]["absmax_mm"] for k in cols] + [30]) * 1.1
lim_F = max([cases[k]["peak_kN"] for k in cols if np.isfinite(cases[k]["peak_kN"])] + [30]) * 1.1
for j, k in enumerate(("exp", "rec", "r2")):
    ax = fig.add_subplot(gs[1, j]); c = cases.get(k)
    if fig13 is not None:
        ax.plot(fig13[0], fig13[1], color="#bbbbbb", lw=1.0, label="US-1 Fig 13 envelope")
    if c and "F" in c:
        n = min(len(c["u"]), len(c["F"])); ax.plot(c["u"][:n], c["F"][:n], color=C[k], lw=0.7, alpha=0.9, label=LBL[k])
        for idx, mk in ((c["i_pos"], "o"), (c["i_neg"], "v")):
            if idx is not None and idx < n:
                ax.plot(c["u"][idx], c["F"][idx], mk, color=C[k], ms=8, mec="white")
                ax.text(c["u"][idx], c["F"][idx], "  {:+.1f} mm".format(c["u"][idx]), fontsize=9, color=C[k], va="center")
    else:
        ax.text(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", color="#888")
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.set_xlim(-lim_u, lim_u); ax.set_ylim(-lim_F, lim_F)
    ax.set_xlabel("OOP displacement (mm)"); ax.set_ylabel("base shear (kN)" if j == 0 else "")
    ax.set_title("(e{}) {}".format(j + 1, LBL[k]), loc="left", fontsize=10, fontweight="bold"); ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(color="#eee")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.axes[0].grid(color="#eee")
for s in ("top", "right"): fig.axes[0].spines[s].set_visible(False)
fig.suptitle("Run 24 -> 25 on the same damaged wall: what the wall did", fontsize=12.5, fontweight="bold")
p2 = os.path.join(OUT_DIR, "fig_focus_response.png"); plt.savefig(p2, bbox_inches="tight"); plt.close(fig); print("-> " + p2)

# ============================ CSV ====================================
with open(os.path.join(OUT_DIR, "focus_run24_25_metrics.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["case", "label", "t_pos_peak_s", "pos_peak_mm", "t_trough_s", "trough_mm", "swing_mm", "swing_over_US1", "absmax_mm",
                                   "peak_shear_kN", "shear_source", "resid_disp_mm", "ringdown_period_s", "table_stroke_mm", "edp_source"])
    for k in cols:
        c = cases[k]
        w.writerow([k, LBL[k], round(c["t_peak"], 4), round(c["u_pos"], 3), round(c["t_neg"], 4) if np.isfinite(c["t_neg"]) else "", round(c["u_neg"], 3) if np.isfinite(c["u_neg"]) else "",
                    round(c["swing_mm"], 3) if np.isfinite(c["swing_mm"]) else "", round(c["swing_mm"] / cases["exp"]["swing_mm"], 3) if cases.get("exp") and np.isfinite(c["swing_mm"]) else "", round(c["absmax_mm"], 3),
                    round(c["peak_kN"], 3) if np.isfinite(c["peak_kN"]) else "", c.get("shear_src", ""), round(c["resid_mm"], 3),
                    round(c["T_ring"], 4) if np.isfinite(c["T_ring"]) else "", round(c["stroke_mm"], 2) if np.isfinite(c["stroke_mm"]) else "", c["edp_src"]])
print("-> " + os.path.join(OUT_DIR, "focus_run24_25_metrics.csv"))
