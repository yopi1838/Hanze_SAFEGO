# -*- coding: ascii -*-
"""
route2_revised_3dec.py -- ROUTE 2, REVISED (Ihsan's sketch, bottom panel):
ONE cycle at the elastic period T1 followed IMMEDIATELY by ONE cycle at the
effective period T_eff, no gap, same phase (both start positive), so the
T1 peak and the T_eff peak sit on one straight peak-to-peak line.

PERIODS (PERIOD_SOURCE):
    "bilinear" (default): T1 and T_eff from bilinear_idealisation/
        bilinear_periods.json, block BILINEAR_KEY ("model" = the v3 pushover,
        "fig13" = the envelope). Rule: K1 = 0.5Vu/d50, dy = Vu/K1, du at 0.8Vu
        on the descending branch, mu = du/dy, T1 = 2 pi sqrt(m/K1),
        T_eff = T1 sqrt(mu/0.8). Model, m_eff 1635: T1 0.112 s, T_eff 0.241 s.
        Output folders get the suffix _bil.
    "secant": the previous rule (T1_A fixed, T_eff = secant to TARGET_D_MM).

    v(t) = A1 sin(2 pi t / T1)            0     <= t <= T1
         = A2 sin(2 pi (t-T1) / Teff)     T1    <  t <= T1 + Teff
    A1 : Sd(T1)   of the cycle = SCALE x record Sd(T1)     (5% damped)
    A2 : Sd(Teff) of the cycle = SCALE x record Sd(Teff)
    "secant" periods (T1 0.0948 s, T_eff 0.361 s), FR76 x 2.0: A1 0.094 m/s,
    A2 0.421 m/s, total 0.455 s, PGA 0.75 g, join step 0.13 g in phase /
    1.37 g inverted. "bilinear" model periods (0.112 / 0.241 s) give a
    shorter, smaller stage B (Sd target 27 mm instead of 51 mm at x2.0).

What changed versus route2_from_save_3dec.py: GAP_S = 0 (was 0.5 s), stage B
phase selectable (was always positive-first after a quiet gap), and the
virgin wall is a case alongside the two restored states.

CASES (each writes its own folder, postprocess_route2_from_save.py-compatible)
    virgin (Part_I_MASON_v8_SmallStrain.sav), x2.0   -> compare control run 25
    end of run 24 (stratC_run_24.sav),        x2.0   -> compare control run 25
    end of run 23 (stratC_run_23.sav),        x1.75  -> compare control run 24
Reference peaks are read from the control folders' channel exports
(rel_disp_top_exp_mm, else 0.5(Ch3+Ch4)-Ch5), never typed in.

SETUP is the control's resume path: restore; large-strain off; dynamic on;
joist contact groups; instrument_history_new.dat; free velocity-z on S; free
rotation-y/z on T_B; local/global damping 0. No joint overrides.

OUTPUT per case (OUT_DIR = route2rev_<case>[_bil][_inv])
    route2_<REC>_s<scale>/          channel CSVs (same layout as the other routes)
    route2_<REC>_s<scale>.sav
    route2_<REC>_s<scale>_vel.txt   the table that was applied
    route2_revised_summary.csv      case, from_save, from_run, scale, phase,
        T_A, SdA_mm, VampA, PGA_A_g, T_B, SdB_mm, VampB, PGA_B_g,
        join_step_g, pulse_dur_s, table_end_mm, SdA_total_mm, SdB_total_mm,
        peak_edp_mm, edp_source, control_ref_run, control_ref_mm
    SdA_total / SdB_total = Sd of the WHOLE two-cycle signal at T1 / Teff
    (the per-stage calibration ignores the other cycle; this shows how much
    the neighbour changes the delivered Sd).
Requires capacity_law.py, spectrum_FR76.csv, US1_fig13_digitised.csv,
instrument_history_new.dat, instrument_history_export_v2.dat next to this
script. Run inside 3DEC: program call 'route2_revised_3dec.py'
"""

import itasca as it
import os, csv, math, sys, glob
import numpy as np

for _p in (os.getcwd(),
           os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)
try:
    import capacity_law as cap
except Exception as _e:
    raise RuntimeError("capacity_law.py must sit next to this script: {}".format(_e))

# ============================ CONFIG =================================
PERIOD_SOURCE = "bilinear"    # "bilinear": T1 and T_eff from bilinear_idealisation/bilinear_periods.json
                              #             (the sketch's rule: K1 at 0.5Vu, du at 0.8Vu, mu = du/dy,
                              #              T_eff = T1 sqrt(mu/0.8); run bilinear_idealisation.py first)
                              # "secant"  : T1 = T1_A below, T_eff = secant to TARGET_D_MM on Fig 13 (previous rule)
BILINEAR_JSON = "bilinear_idealisation/bilinear_periods.json"
BILINEAR_KEY  = "model"       # "model" (pushover_pos_v3, prediction mode) or "fig13" (envelope, validation mode)
STAGE_B_PHASE = "same"        # "same" (bottom panel: both cycles start positive)
                              # or "inverted" (stage B starts negative)
CASES = [
    # (case name, save to restore,                          from_run, SCALE, control run, )
    ("virgin", "Part_I_MASON_v8_SmallStrain.sav",           0,        2.0,   25),
    ("from24", "stratF_full_results_US1/stratC_run_24.sav", 24,       2.0,   25),
    ("from23", "stratF_full_results_US1/stratC_run_23.sav", 23,       1.75,  24),
]
CONTROL_DIR = "stratF_full_results_US1"
LABEL       = "FR76"
SPECTRUM    = "spectrum_FR76.csv"
T1_A        = 0.0948          # s, model elastic period (stage A) -- used only when PERIOD_SOURCE = "secant"
XI          = 0.05
DELTA_T     = 0.005
N_A, N_B    = 1.0, 1.0        # one cycle each -- the revised scheme
GAP_S       = 0.0             # no gap -- the revised scheme
TAIL_SEC    = 2.5
DRIVE_GROUPS = ["S", "T_B"]
CAP_CSV     = "US1_fig13_digitised.csv"
CAP_DIR     = "neg"
M_EFF       = 1635.0
TARGET_D_MM = 29.5
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
if STAGE_B_PHASE not in ("same", "inverted"):
    raise RuntimeError("STAGE_B_PHASE must be 'same' or 'inverted'")
PHASE_SIGN = 1.0 if STAGE_B_PHASE == "same" else -1.0
DIR_SUFFIX = ("_bil" if PERIOD_SOURCE == "bilinear" else "") + ("" if STAGE_B_PHASE == "same" else "_inv")

# ============================ Sd INTEGRATOR (inlined, self-tested) ====
def sd_spectrum(ag, dt, Ts, xi=0.05):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2.0 * np.pi / Ts; c = 2.0 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0] - c * v - w2 * u; umax = np.abs(u)
    for i in range(1, len(ag)):
        u_new = u + dt * v + 0.25 * dt * dt * a
        v_new = v + 0.5 * dt * a
        a_new = (-ag[i] - c * v_new - w2 * u_new) / (1.0 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = u_new + 0.25 * dt * dt * a_new; v = v_new + 0.5 * dt * a_new; a = a_new
        umax = np.maximum(umax, np.abs(u))
    return umax

def _sd_selftest():
    T, xi, A, dt = 0.4, 0.05, 1.0, 0.001
    t = np.arange(int(60 * T / dt)) * dt
    got = sd_spectrum(A * np.sin(2 * np.pi / T * t), dt, [T], xi)[0]
    want = A / ((2 * np.pi / T) ** 2 * 2 * xi)
    if abs(got - want) / want > 0.01:
        raise RuntimeError("Sd integrator self-test failed: {} vs {}".format(got, want))
_sd_selftest()

# ============================ HELPERS ================================
def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")

def read_hist(path):
    try:
        d = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None

def find_ch(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else None

def peak_edp_mm(folder):
    """(peak, source): rel_disp_top_exp_mm channel, else 0.5(Ch3+Ch4)-Ch5 in mm."""
    x = find_ch(folder, KEY_DISP_EXP)
    if x is not None:
        return float(np.max(np.abs(x - x[0]))), KEY_DISP_EXP
    c3, c4, c5 = find_ch(folder, CH3), find_ch(folder, CH4), find_ch(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, "n/a"
    n = min(len(c3), len(c4), len(c5)); rel = 0.5 * (c3[:n] + c4[:n]) - c5[:n]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0, "0.5(Ch3+Ch4)-Ch5"

def control_ref(run_no):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    if not hits:
        return None, "control folder Run{:02d}_* not found".format(run_no)
    return peak_edp_mm(hits[0])

def grid_period(T):
    spc = max(1, int(round(T / DELTA_T))); return spc * DELTA_T, spc

def one_cycle(Vamp, T, ncyc=1.0):
    """velocity of ncyc sine cycles at grid period Tg, samples 0..ncyc*Tg inclusive."""
    Tg, spc = grid_period(T); npl = int(round(ncyc * spc)); w = 2.0 * math.pi / Tg
    tt = np.arange(npl + 1) * DELTA_T
    return tt, Vamp * np.sin(w * tt), Tg

def calibrate_vel(T, Sd_target, ncyc=1.0):
    _, v, Tg = one_cycle(1.0, T, ncyc)
    sd0 = float(sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [Tg], XI)[0])
    if sd0 <= 0:
        raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg

def build_signal(Vamp_A, Vamp_B):
    """Stage A cycle, then stage B cycle starting one step after A's last
    (zero) sample: velocity continuous, no gap, no duplicated zero."""
    tA, vA, TgA = one_cycle(Vamp_A, T1_A, N_A)
    tB, vB, TgB = one_cycle(PHASE_SIGN * Vamp_B, T_eff, N_B)
    gap = np.zeros(int(round(GAP_S / DELTA_T)))
    v = np.r_[vA, gap, vB[1:]]
    t = np.arange(len(v)) * DELTA_T
    return t, v, vA, vB, TgA, TgB

def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81)

def join_step_g(vA, vB):
    """acceleration just before / just after the join (finite-difference slopes)."""
    dA = (vA[-1] - vA[-2]) / DELTA_T
    dB = (vB[1] - vB[0]) / DELTA_T
    return abs(dB - dA) / 9.81

# ============================ STAGE PERIODS (fixed law) ==============
if PERIOD_SOURCE == "bilinear":
    import json
    if not os.path.isfile(BILINEAR_JSON):
        raise RuntimeError("run bilinear_idealisation.py first: {} not found".format(BILINEAR_JSON))
    _b = json.load(open(BILINEAR_JSON))
    if BILINEAR_KEY not in _b:
        raise RuntimeError("{} has no '{}' block (blocks present: {}). For 'model', bilinear_idealisation.py must find the v3 "
                           "pushover_pos.csv (next to it or in pushover_results/) -- rerun it and check its console output.".format(
                               BILINEAR_JSON, BILINEAR_KEY, ", ".join(k for k in _b if k != "rule") or "none"))
    _b = _b[BILINEAR_KEY]
    T1_A, T_eff = float(_b["T1_s"]), float(_b["Teff_s"])
    period_note = "bilinear ({}; K1 {:.2f} kN/mm, mu {:.2f}, branches {})".format(BILINEAR_KEY, _b["K1_kNmm"], _b["mu"], "+".join(_b["branches_used"]))
elif PERIOD_SOURCE == "secant":
    law = cap.CapacityLaw(CAP_CSV, m_eff=M_EFF, direction=CAP_DIR, label="route2rev")
    T_eff, cap_note = law.t_eff(TARGET_D_MM)
    period_note = "secant to {:.1f} mm on Fig 13 {} {}".format(TARGET_D_MM, CAP_DIR, cap_note or "")
else:
    raise RuntimeError("PERIOD_SOURCE must be 'bilinear' or 'secant'")
FR = np.genfromtxt(SPECTRUM, delimiter=",", skip_header=1)
def sd_record(T):
    return float(np.interp(T, FR[:, 0], FR[:, 1]))
print("ROUTE 2 REVISED: one cycle T1 = {:.4f} s then one cycle T_eff = {:.4f} s  [{}; m_eff {:.0f}] ; "
      "no gap ; stage B phase = {}".format(T1_A, T_eff, period_note, M_EFF, STAGE_B_PHASE))

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ RUN THE CASES ==========================
for (case, save, from_run, SCALE, ctrl_run) in CASES:
    OUT_DIR = "route2rev_" + case + DIR_SUFFIX
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.isfile(save):
        print("!! save not found, case skipped: {}".format(save)); continue

    # ---- signal for this scale
    Sd_A = SCALE * sd_record(T1_A); Sd_B = SCALE * sd_record(T_eff)
    Vamp_A, TgA = calibrate_vel(T1_A, Sd_A, N_A)
    Vamp_B, TgB = calibrate_vel(T_eff, Sd_B, N_B)
    t_all, v_all, vA, vB, _, _ = build_signal(Vamp_A, Vamp_B)
    pulse_dur = float(t_all[-1])
    ag_all = np.gradient(v_all, DELTA_T)
    Sd_tot = sd_spectrum(ag_all, DELTA_T, [TgA, TgB], XI)
    step_g = join_step_g(vA, vB)
    end_mm = float(np.sum(0.5 * (v_all[1:] + v_all[:-1]) * np.diff(t_all))) * 1000.0   # trapezoid, numpy-version agnostic
    scale_str = "{:.2f}".format(SCALE).replace(".", "p")
    run_label = "route2_{}_s{}".format(LABEL, scale_str)
    vel_path = os.path.join(OUT_DIR, run_label + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(run_label, len(t_all)))
        for a, b in zip(t_all, v_all):
            f.write("{:.6f}\t{:.9e}\n".format(a, b))

    print("=" * 70)
    print("ROUTE 2 REVISED [{}] from {} (end of run {})  scale {:.2f}  -> compare control run {}".format(
        case, save, from_run, SCALE, ctrl_run))
    print("  stage A: T {:.4f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.3f} g".format(TgA, Sd_A * 1e3, Vamp_A, pga_g(vA)))
    print("  stage B: T {:.4f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.3f} g  phase {}".format(TgB, Sd_B * 1e3, Vamp_B, pga_g(vB), STAGE_B_PHASE))
    print("  join acceleration step {:.2f} g ; whole-signal Sd(T1) {:.2f} mm, Sd(Teff) {:.2f} mm (targets {:.2f}, {:.2f})".format(
        step_g, Sd_tot[0] * 1e3, Sd_tot[1] * 1e3, Sd_A * 1e3, Sd_B * 1e3))
    print("  table {:.3f}s (ends at {:+.2f} mm) + tail {:.1f}s".format(pulse_dur, end_mm, TAIL_SEC))

    # ---- restore + setup exactly as the control's resume path
    it.command("model restore '{}'".format(cmd_path(save)))
    it.command("python-reset-state false")
    it.command("model large-strain off")
    it.command("model dynamic active on")
    it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
    it.command("call 'instrument_history_new.dat'")
    it.command("block free velocity-z range group 'S'")
    it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
block mech damp local 0.0
block mech damp global 0.0
""")
    it.command("@Record_Disp")
    d0 = 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1000.0
    print("  restored state: top displacement {:.3f} mm (absolute)".format(d0))

    # ---- apply the signal
    it.command("table 'r2' import '{}'".format(cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for grp in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table 'r2' range group '{}'".format(grp))
    it.command("model solve dynamic time {:.6f}".format(pulse_dur))
    for grp in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(grp))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    # ---- save + export
    it.command("model save '{}'".format(cmd_path(os.path.join(OUT_DIR, run_label + ".sav"))))
    full_folder = os.path.join(OUT_DIR, run_label)
    os.makedirs(full_folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(full_folder)))
    it.command("[runlabel='{}']".format(run_label))
    it.command("call 'instrument_history_export_v2.dat'")

    peak, src = peak_edp_mm(full_folder)
    ref, ref_src = control_ref(ctrl_run)
    print("\n  PEAK EDP (Route 2 revised, {}) = {}".format(src, "{:.2f} mm".format(peak) if peak is not None else "n/a"))
    print("  control run {} ({}) = {}".format(ctrl_run, ref_src, "{:.2f} mm".format(ref) if ref is not None else "n/a"))
    if peak is not None and ref:
        print("  ratio revised Route 2 / record = {:.2f}".format(peak / ref))

    sum_path = os.path.join(OUT_DIR, "route2_revised_summary.csv")
    new = not os.path.isfile(sum_path)
    with open(sum_path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["case", "from_save", "from_run", "scale", "stage_B_phase", "period_source",
                        "T_A_s", "SdA_mm", "VampA_mps", "PGA_A_g",
                        "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g",
                        "join_step_g", "pulse_dur_s", "table_end_mm", "SdA_total_mm", "SdB_total_mm",
                        "N_A", "N_B", "gap_s", "m_eff",
                        "peak_edp_mm", "edp_source", "control_ref_run", "control_ref_mm"])
        w.writerow([case, save, from_run, SCALE, STAGE_B_PHASE, PERIOD_SOURCE + ("/" + BILINEAR_KEY if PERIOD_SOURCE == "bilinear" else ""),
                    round(TgA, 5), round(Sd_A * 1e3, 3), round(Vamp_A, 5), round(pga_g(vA), 4),
                    round(TgB, 5), round(Sd_B * 1e3, 3), round(Vamp_B, 5), round(pga_g(vB), 4),
                    round(step_g, 3), round(pulse_dur, 4), round(end_mm, 3), round(Sd_tot[0] * 1e3, 3), round(Sd_tot[1] * 1e3, 3),
                    N_A, N_B, GAP_S, M_EFF,
                    round(peak, 3) if peak is not None else "", src, ctrl_run,
                    round(ref, 3) if ref else ""])
    print("  -> {}".format(sum_path))

print("\ndone. Compare peak_edp_mm with control_ref_mm per case; the virgin case is the one-shot")
print("counterpart of the earlier route2_results run (4.12 mm with the 0.5 s gap).")