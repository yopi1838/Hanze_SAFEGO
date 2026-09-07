# -*- coding: ascii -*-
"""
route2_sequential_pulse_3dec.py -- ROUTE 2 (Ihsan): the sequential
two-stage pulse, SMALL STRAIN.

THE IDEA (Ihsan, meeting 2026-08-31)
    Apply two sinusoids back to back in ONE run:
      Stage A -- at the elastic period T1, amplitude matched to the
                 record's spectral displacement Sd(record, scale, T1).
      Stage B -- at the EFFECTIVE period T_eff from the bilinear capacity
                 fit, amplitude matched to Sd(record, scale, T_eff).
    Stage A represents the near-elastic onset; stage B puts demand at the
    period the softened wall actually occupies. This is the paper-cleanest
    hybrid of "fixed law (capacity curve + secant + m_eff), adaptive state
    (the displacement the wall has reached)".

  ** WHAT STAGE B ACTUALLY COMMANDS -- read before interpreting. **
    At scale 2.0 (run 25) and T_eff = 0.36 s, FR76's own 5%-damped Sd is
    ~51 mm. So stage B commands a pulse whose LINEAR-SDOF spectral
    displacement is 51 mm -- NOT 29.5. Route 2 is therefore not "aim for
    29.5"; it applies the scheme and we see what the wall does. If the
    wall returns ~29.5, the 51-vs-29.5 gap is the linear-SDOF-vs-rocking
    difference (a real, reportable finding). If it returns ~51, the scheme
    over-predicts. Decide which BEFORE looking, as with Route 1.

  ** SMALL STRAIN -- the caveat you asked for, stated plainly. **
    model large-strain OFF freezes block positions and contact detection:
    no toe migration, no P-Delta geometric softening -- the very mechanism
    large-amplitude rocking uses. So a shortfall here CANNOT be cleanly
    blamed on the scheme vs the frozen geometry; that confound must be in
    any conclusion. The reason to use it anyway: it matches the existing
    small-strain runs (Strategy F control, ROCKPERIOD), so Route 1 (rerun
    small strain -- see note), Route 2 and the control become comparable.
    NOTE: the Route 1 driver as delivered used large-strain ON. Flip it to
    OFF for a clean 3-way comparison, or run all three large strain.

CAPACITY SOURCE decides the claim (as capacity_law documents): the
digitised Fig 13 envelope makes this VALIDATION against US-1; the model
pushover makes it prediction. Default here is the envelope, neg branch
(US-1's large excursions were negative), m_eff = 1635 kg (the mid-height
mechanism mass that reproduced Table 4's secant period). T_eff(29.5 mm,
neg, 1635) = 0.36 s.

SETUP mirrors Route 1 and the Strategy C/F drivers exactly, except
large-strain is OFF: restore Part_I_MASON_v8_SmallStrain, no damping,
joist groups, instrument_history_new.dat, drive S and T_B with a
velocity-z table (multiplier 1.0), ring down, export.

Requires capacity_law.py and spectrum_FR76.csv next to this script
(the Sd integrator is inlined and self-tests at start-up).

OUTPUT
    <OUT_DIR>/route2_<record>_s<scale>/  channel CSVs
    <OUT_DIR>/route2_summary.csv : record, scale, T_A, SdA_mm, T_B,
        SdB_mm, N_A, N_B, gap_s, peak_edp_mm
"""

import itasca as it
import os, csv, math, sys
import numpy as np

# 3DEC's embedded Python does not put the working directory on sys.path,
# so sibling modules are invisible unless added explicitly (same fix your
# strategy drivers use for ratcheting_pulse.py).
for _p in (os.getcwd(),
           os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import capacity_law as cap   # bilinear T_eff (self-tests on import)
except Exception as _e:
    raise RuntimeError("capacity_law.py must sit next to this script "
                       "(or in the 3DEC working directory). Error: {}".format(_e))

# ---- Newmark average-acceleration Sd spectrum, INLINED so this driver
# ---- has no sdof.py dependency. Self-tested below before any use.
def sd_spectrum(ag, dt, Ts, xi=0.05):
    ag = np.asarray(ag, float)
    Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2.0 * np.pi / Ts
    c = 2.0 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w)
    a = -ag[0] - c * v - w2 * u
    umax = np.abs(u)
    for i in range(1, len(ag)):
        u_new = u + dt * v + 0.25 * dt * dt * a
        v_new = v + 0.5 * dt * a
        a_new = (-ag[i] - c * v_new - w2 * u_new) / (1.0 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = u_new + 0.25 * dt * dt * a_new
        v = v_new + 0.5 * dt * a_new
        a = a_new
        umax = np.maximum(umax, np.abs(u))
    return umax

def _sd_selftest():
    # resonant sinusoid: steady-state Sd = A / (w^2 * 2 xi), must match to 1%
    T, xi, A, dt = 0.4, 0.05, 1.0, 0.001
    t = np.arange(int(60 * T / dt)) * dt
    got = sd_spectrum(A * np.sin(2 * np.pi / T * t), dt, [T], xi)[0]
    want = A / ((2 * np.pi / T) ** 2 * 2 * xi)
    if abs(got - want) / want > 0.01:
        raise RuntimeError("inlined Sd integrator failed self-test: {} vs {}".format(got, want))
_sd_selftest()

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ CONFIG =================================
BASE_SAVE   = "part_I_mason_LS.sav"
SPECTRUM    = "spectrum_FR76.csv"     # record Sd(T), the amplitude target
OUT_DIR     = "route2_results"
LABEL       = "FR76"
SCALE       = 2.0            # run 25 = FR76 x 2.0
T1_INIT     = 0.0948         # s, elastic period (stage A)
XI          = 0.05
DELTA_T     = 0.005
N_A         = 1.0            # stage A cycles (honest single lobe, as Route 1)
N_B         = 1.0            # stage B cycles
GAP_S       = 0.5           # quiet time between stages (PHASE-SENSITIVE;
                            # rerun 0.25 / 1.0 for the band)
TAIL_SEC    = 3.0
DRIVE_GROUPS = ["S", "T_B"]

# --- capacity law (stage B period) ---
CAP_CSV      = "US1_fig13_digitised.csv"   # envelope = VALIDATION mode
CAP_DIR      = "neg"        # US-1's large excursions were negative
M_EFF        = 1635.0       # mid-height mechanism mass (reproduces Table 4 T_sec)
TARGET_D_MM  = 29.5         # the observed run-25 displacement; sets T_eff

os.makedirs(OUT_DIR, exist_ok=True)

# ============================ STAGE PERIODS + AMPLITUDES =============
law = cap.CapacityLaw(CAP_CSV, m_eff=M_EFF, direction=CAP_DIR, label="route2")
print(law.report())
T_eff, cap_note = law.t_eff(TARGET_D_MM)
print("stage B period: T_eff({:.1f} mm) = {:.4f} s   {}".format(
    TARGET_D_MM, T_eff, cap_note or ""))

FR = np.genfromtxt(SPECTRUM, delimiter=",", skip_header=1)
def sd_record(T):
    return float(np.interp(T, FR[:, 0], FR[:, 1]))

def pulse_vel(Vamp, T, ncyc):
    # snap T to the table grid so N cycles end exactly at zero (Route 1 fix)
    spc = max(1, int(round(T / DELTA_T)))
    Tg = spc * DELTA_T
    npl = int(round(ncyc * spc))
    w = 2.0 * math.pi / Tg
    tt = np.arange(npl + 1) * DELTA_T
    return tt, Vamp * np.sin(w * tt), Tg

def calibrate_vel(T, Sd_target, ncyc):
    """Velocity amplitude so the pulse's 5%-damped Sd at T = Sd_target.
    Linear in Vamp; calibrated on the differentiated velocity (the actual
    base acceleration), using the verified sdof spectrum."""
    _, v, Tg = pulse_vel(1.0, T, ncyc)
    a = np.gradient(v, DELTA_T)
    sd0 = float(sd_spectrum(a, DELTA_T, [Tg], XI)[0])
    if sd0 <= 0:
        raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg

Sd_A = SCALE * sd_record(T1_INIT)
Sd_B = SCALE * sd_record(T_eff)
Vamp_A, TgA = calibrate_vel(T1_INIT, Sd_A, N_A)
Vamp_B, TgB = calibrate_vel(T_eff,  Sd_B, N_B)

# build the two stages and splice: A -> gap -> B -> tail
tA, vA, _ = pulse_vel(Vamp_A, T1_INIT, N_A)
tB, vB, _ = pulse_vel(Vamp_B, T_eff, N_B)
gap_n = int(round(GAP_S / DELTA_T))
tail_n = int(round(TAIL_SEC / DELTA_T))
v_all = np.r_[vA, np.zeros(gap_n), vB, np.zeros(tail_n)]
t_all = np.arange(len(v_all)) * DELTA_T
pulse_dur = float(tA[-1] + GAP_S + tB[-1])   # end of stage B (before tail)

print("=" * 66)
print("ROUTE 2 -- sequential two-stage pulse (SMALL STRAIN)")
print("  stage A: T {:.4f}s  Sd {:.2f}mm  N {:.1f}  Vamp {:.4f} m/s  "
      "PGA {:.3f} g".format(TgA, Sd_A * 1000, N_A, Vamp_A,
                            np.max(np.abs(np.gradient(vA, DELTA_T))) / 9.81))
print("  stage B: T {:.4f}s  Sd {:.2f}mm  N {:.1f}  Vamp {:.4f} m/s  "
      "PGA {:.3f} g".format(TgB, Sd_B * 1000, N_B, Vamp_B,
                            np.max(np.abs(np.gradient(vB, DELTA_T))) / 9.81))
print("  gap {:.2f}s   stage-B Sd target {:.1f}mm >> observed 29.5mm "
      "(see header)".format(GAP_S, Sd_B * 1000))
print("=" * 66)

scale_str = "{:.2f}".format(SCALE).replace(".", "p")
run_label = "route2_{}_s{}".format(LABEL, scale_str)
vel_path = os.path.join(OUT_DIR, run_label + "_vel.txt")
with open(vel_path, "w", newline="\n") as f:
    f.write("{}\n{}\t0\n".format(run_label, len(t_all)))
    for a, b in zip(t_all, v_all):
        f.write("{:.6f}\t{:.9e}\n".format(a, b))

# ============================ MODEL SETUP (SMALL STRAIN) ============
if not os.path.isfile(BASE_SAVE):
    raise RuntimeError("base save not found: {}".format(BASE_SAVE))
it.command("model restore '{}'".format(BASE_SAVE.replace(".sav", "")))
it.command("python-reset-state false")
it.command("model large-strain on")     # <<< SMALL STRAIN, as requested
it.command("model dynamic active on")
it.command("block mech damp local 0.0")
it.command("block mech damp global 0.0")
it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
it.command("call 'instrument_history_new.dat'")
it.command("block free velocity-z range group 'S'")
it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")

it.command("@Record_Disp")
d0 = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
            it.fish.get("Top_Quarter_B_Disp")) * 1000.0
print("baseline top displacement = {:.4f} mm".format(d0))

# ============================ APPLY THE TWO-STAGE PULSE ==============
it.command("table 'r2' import '{}'".format(vel_path.replace("\\", "/")))
it.command("model dynamic time-total 0")
for grp in DRIVE_GROUPS:
    it.command("block apply velocity-z 1.0 table 'r2' range group '{}'".format(grp))
# solve through both stages + a little of the tail so the peak is captured
it.command("model solve dynamic time {:.6f}".format(pulse_dur + TAIL_SEC))
for grp in DRIVE_GROUPS:
    it.command("block gridpoint apply-remove velocity-z range group '{}'".format(grp))

it.command("@Record_Disp")
d_end = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
               it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0
print("end-of-run top displacement (rel. baseline) = {:.3f} mm".format(d_end))

# ============================ EXPORT + PEAK =========================
full_folder = os.path.join(OUT_DIR, run_label)
os.makedirs(full_folder, exist_ok=True)
it.command("[exportdir='{}']".format(full_folder.replace("\\", "/")))
it.command("[runlabel='{}']".format(run_label))
it.command("call 'instrument_history_export_v2.dat'")
it.command("model save '{}'".format(
    os.path.join(OUT_DIR, run_label + ".sav").replace("\\", "/")))

def _peak_from_channels(folder):
    import glob
    def rd(pat):
        hits = sorted(glob.glob(os.path.join(folder, "*" + pat + "*.csv")))
        if not hits:
            return None
        try:
            a = np.genfromtxt(hits[0], skip_header=2)
            return a[:, 1] if a.ndim == 2 and a.shape[1] >= 2 else None
        except Exception:
            return None
    a = rd("Top_Quarter_A_Disp"); b = rd("Top_Quarter_B_Disp")
    tab = rd("Record_Disp")
    if a is None or b is None:
        return None
    top = 0.5 * (a + b)
    if tab is not None:
        n = min(len(top), len(tab)); rel = top[:n] - tab[:n]
    else:
        rel = top - top[0]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0

peak_edp = _peak_from_channels(full_folder)
print("\nPEAK EDP (table-referenced, paper Eq.1) = {}".format(
    "{:.2f} mm".format(peak_edp) if peak_edp is not None else "n/a"))
print("stage-B Sd target = {:.1f} mm | observed run-25 = 29.5 mm | "
      "strat-F control = 41.6 mm".format(Sd_B * 1000))

sum_path = os.path.join(OUT_DIR, "route2_summary.csv")
new = not os.path.isfile(sum_path)
with open(sum_path, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["record", "scale", "T_A_s", "SdA_mm", "T_B_s", "SdB_mm",
                    "N_A", "N_B", "gap_s", "m_eff", "peak_edp_mm"])
    w.writerow([LABEL, SCALE, round(TgA, 5), round(Sd_A * 1000, 3),
                round(TgB, 5), round(Sd_B * 1000, 3), N_A, N_B, GAP_S,
                M_EFF, round(peak_edp, 3) if peak_edp is not None else ""])
print("-> {}".format(sum_path))
print("\ninterpretation (decide BEFORE looking):")
print("  ~29-35 mm : wall reaches the observed response; the 51mm Sd target")
print("              vs 29.5mm response IS the linear-SDOF-vs-rocking gap.")
print("  ~45-55 mm : wall tracks the linear Sd target -> scheme over-predicts")
print("              (and in SMALL strain, may be missing rocking limit).")
print("  <15 mm    : two-stage pulse under-delivers -- but small strain")
print("              confounds this; rerun large strain before concluding.")