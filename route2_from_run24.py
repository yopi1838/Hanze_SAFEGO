# -*- coding: ascii -*-
"""
route2_from_save_3dec.py -- ROUTE 2 (two sequential sinusoids) applied to a
RESTORED Strategy F state, not to the virgin wall.

PURPOSE
    Test whether the two-stage pulse reproduces what the RECORD does to the
    MODEL in the same state:
      case 1: restore stratC_run_24.sav (end of FR76 x1.75), apply Route 2 at
              scale 2.0  -> compare with the control's run 25 (FR76 x2.0).
      case 2: restore stratC_run_23.sav (end of FR76 x1.5),  apply Route 2 at
              scale 1.75 -> compare with the control's run 24 (FR76 x1.75).
    Same state, same intensity, record replaced by the sinusoids. The
    reference values are read from the control folder's own channel exports
    (rel_disp_top_exp_mm, else 0.5(Ch3+Ch4)-Ch5), never typed in.

THE RULE (unchanged from route2_sequential_pulse_3dec.py)
    Stage A: period = the model's ELASTIC period T1_A (0.0948 s, as agreed),
             amplitude so that the pulse's 5%-damped Sd(T1_A) = SCALE x the
             record's Sd(T1_A).
    Stage B: period = T_eff from the digitised Fig 13 envelope (neg branch,
             ENVELOPE force at TARGET_D_MM, m_eff = M_EFF) -> 0.361 s;
             amplitude so that Sd(T_eff) = SCALE x record Sd(T_eff).
    N_A = N_B = 1 cycle, GAP_S quiet between stages.
    NB the restored wall is cracked: its ring-down period at run 24 is
    0.1024 s, not 0.0948 s. Stage A is still tuned to the elastic period
    because that is the rule as agreed; the number is printed so the
    difference is visible.

SETUP mirrors strategy_F_full_3dec.setup_model_for_dynamic() on resume:
    model restore <save>; large-strain off; dynamic on; joist contact
    groups; instrument_history_new.dat; free velocity-z on S; free
    rotation-y/z on T_B; local/global damping 0. No Maxwell/Rayleigh, no
    joint overrides (the control used none). Everything else comes back
    with the save.

PEAK EDP is 0.5(Ch3+Ch4) - Ch5 (paper Eq. 1) via the FISH channel
rel_disp_top_exp_mm when exported, else rebuilt from the Channel_3/4/5
CSVs. NOT Record_Disp: that is the integrated input table and differs from
the base block's actual displacement (0.5 mm on the one-shot, more on long
records) -- the cause of the earlier 4.61 vs 4.12 mm discrepancy.

OUTPUT
    <OUT_DIR>/route2_<REC>_s<scale>/   channel CSVs   (postprocess_routes.py
        reads this folder; set SINGLE_SHOT_RUN to the compared control run)
    <OUT_DIR>/route2_<REC>_s<scale>.sav
    <OUT_DIR>/route2_from_save_summary.csv : from_save, from_run, scale,
        T_A, SdA_mm, VampA, PGA_A_g, T_B, SdB_mm, VampB, PGA_B_g, N_A, N_B,
        gap_s, m_eff, peak_edp_mm, edp_source, control_ref_run,
        control_ref_mm
Requires capacity_law.py, spectrum_FR76.csv, US1_fig13_digitised.csv,
instrument_history_new.dat, instrument_history_export_v2.dat in the working
directory (same as the other drivers).
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
CASES = [
    # (save to restore,                              from_run, SCALE, control run to compare, OUT_DIR)
    ("stratF_full_results_US1/stratC_run_24.sav",   24,       2.0,   25, "route2_from24"),
    ("stratF_full_results_US1/stratC_run_23.sav",   23,       1.75,  24, "route2_from23"),
]
CONTROL_DIR = "stratF_full_results_US1"   # RunNN_REC_sXpYY folders with channel CSVs
LABEL       = "FR76"
SPECTRUM    = "spectrum_FR76.csv"
T1_A        = 0.0948          # s, model elastic period for stage A (as agreed)
XI          = 0.05
DELTA_T     = 0.005
N_A, N_B    = 1.0, 1.0
GAP_S       = 0.5
TAIL_SEC    = 2.5             # ring-down after the table is removed (control uses 2.5)
DRIVE_GROUPS = ["S", "T_B"]
CAP_CSV     = "US1_fig13_digitised.csv"
CAP_DIR     = "neg"
M_EFF       = 1635.0
TARGET_D_MM = 29.5
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"

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

def pulse_vel(Vamp, T, ncyc):
    spc = max(1, int(round(T / DELTA_T))); Tg = spc * DELTA_T
    npl = int(round(ncyc * spc)); w = 2.0 * math.pi / Tg
    tt = np.arange(npl + 1) * DELTA_T
    return tt, Vamp * np.sin(w * tt), Tg

def calibrate_vel(T, Sd_target, ncyc):
    _, v, Tg = pulse_vel(1.0, T, ncyc)
    sd0 = float(sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [Tg], XI)[0])
    if sd0 <= 0:
        raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg

def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81)

# ============================ STAGE PERIODS (fixed law) ==============
law = cap.CapacityLaw(CAP_CSV, m_eff=M_EFF, direction=CAP_DIR, label="route2")
T_eff, cap_note = law.t_eff(TARGET_D_MM)
FR = np.genfromtxt(SPECTRUM, delimiter=",", skip_header=1)
def sd_record(T):
    return float(np.interp(T, FR[:, 0], FR[:, 1]))
print("stage A period (model elastic) = {:.4f} s ; stage B period T_eff({:.1f} mm, {}, m_eff {:.0f}) = {:.4f} s {}".format(
    T1_A, TARGET_D_MM, CAP_DIR, M_EFF, T_eff, cap_note or ""))

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ RUN THE CASES ==========================
for (save, from_run, SCALE, ctrl_run, OUT_DIR) in CASES:
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.isfile(save):
        print("!! save not found, case skipped: {}".format(save)); continue

    # ---- pulse for this scale
    Sd_A = SCALE * sd_record(T1_A); Sd_B = SCALE * sd_record(T_eff)
    Vamp_A, TgA = calibrate_vel(T1_A, Sd_A, N_A)
    Vamp_B, TgB = calibrate_vel(T_eff, Sd_B, N_B)
    tA, vA, _ = pulse_vel(Vamp_A, T1_A, N_A)
    tB, vB, _ = pulse_vel(Vamp_B, T_eff, N_B)
    v_all = np.r_[vA, np.zeros(int(round(GAP_S / DELTA_T))), vB]
    t_all = np.arange(len(v_all)) * DELTA_T
    pulse_dur = float(t_all[-1])
    scale_str = "{:.2f}".format(SCALE).replace(".", "p")
    run_label = "route2_{}_s{}".format(LABEL, scale_str)
    vel_path = os.path.join(OUT_DIR, run_label + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(run_label, len(t_all)))
        for a, b in zip(t_all, v_all):
            f.write("{:.6f}\t{:.9e}\n".format(a, b))

    print("=" * 70)
    print("ROUTE 2 from {} (end of run {})  scale {:.2f}  -> compare control run {}".format(
        save, from_run, SCALE, ctrl_run))
    print("  stage A: T {:.4f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.3f} g".format(TgA, Sd_A * 1e3, Vamp_A, pga_g(vA)))
    print("  stage B: T {:.4f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.3f} g".format(TgB, Sd_B * 1e3, Vamp_B, pga_g(vB)))
    print("  gap {:.2f}s  table {:.2f}s + tail {:.1f}s".format(GAP_S, pulse_dur, TAIL_SEC))

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
    print("  restored state: top displacement {:.3f} mm (absolute, carries the run history)".format(d0))

    # ---- apply the two-stage pulse (same commands as the control's execute_run)
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
    print("\n  PEAK EDP (Route 2, {}) = {}".format(src, "{:.2f} mm".format(peak) if peak is not None else "n/a"))
    print("  control run {} ({}) = {}".format(ctrl_run, ref_src, "{:.2f} mm".format(ref) if ref is not None else "n/a"))
    if peak is not None and ref:
        print("  ratio Route 2 / record = {:.2f}".format(peak / ref))

    sum_path = os.path.join(OUT_DIR, "route2_from_save_summary.csv")
    new = not os.path.isfile(sum_path)
    with open(sum_path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["from_save", "from_run", "scale", "T_A_s", "SdA_mm", "VampA_mps", "PGA_A_g",
                        "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g", "N_A", "N_B", "gap_s", "m_eff",
                        "peak_edp_mm", "edp_source", "control_ref_run", "control_ref_mm"])
        w.writerow([save, from_run, SCALE, round(TgA, 5), round(Sd_A * 1e3, 3), round(Vamp_A, 5), round(pga_g(vA), 4),
                    round(TgB, 5), round(Sd_B * 1e3, 3), round(Vamp_B, 5), round(pga_g(vB), 4), N_A, N_B, GAP_S, M_EFF,
                    round(peak, 3) if peak is not None else "", src, ctrl_run,
                    round(ref, 3) if ref else ""])
    print("  -> {}".format(sum_path))

print("\ndone. Read: Route 2 / record ratio near 1 in BOTH cases = the scheme reproduces the record on this model;")
print("only case 2 (from run 23, x1.75) near 1 = stage B is the shortfall; neither = stage A.")