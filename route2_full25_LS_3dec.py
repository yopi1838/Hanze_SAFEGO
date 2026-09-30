# -*- coding: ascii -*-
"""
route2_full25_cases_3dec.py -- the REVISED Route 2 full 25-run protocol, run
as a 2x2 of geometry (small / large strain) x damping (none / Maxwell 1.5%).
Identical signal, protocol, instrumentation, checkpointing and log in every
case; one CASE constant selects which of the four this invocation is.

Supersedes route2_full25_revised_3dec.py (= case SS_NODAMP) and the earlier
route2_full25_revised_damped_3dec.py (= case SS_MAXWELL_1p5).

WHY: the undamped model has no energy sink except contact friction and
impact. That flatters the late runs (the wall keeps whatever a pulse gives
it) and it also makes the "residual" tilt/displacement meaningless, because
the wall is still ringing when the 2.5 s tail ends. Adding damping tests
whether the run-24 agreement (Route 2 ~32 mm vs US-1 29.5 mm) is a property
of the wall or a property of having no dissipation.

DAMPING AND GEOMETRY ARE TAKEN FROM THE CONTROL, NOT CHOSEN HERE.
    The matching control is the standalone Strategy F case F_LS_MAXWELL_1p5
    (output stratF_full_LS_MAXWELL_1p5), which runs:
        block mech damp local 0.0
        block mech damp global 0.0
        block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 0.0104 40.0000
        DAMP_RATIO = 0.0   (no Rayleigh)
        LARGE_STRAIN = "on",  base save Part_I_MASON_v8.sav
    The Maxwell command below is that string verbatim: 1.5% over 1-40 Hz,
    the same three (weight, frequency) pairs as strategy_C_3dec_maxwell.py
    and sweep_cases.py, so a damped Route 2 is directly comparable to the
    damped record run.

    NOT THE SAME CONTROL AS BEFORE: stratF_full_results_US1 is the UNDAMPED
    record sequence. Compare this run against stratF_full_LS_MAXWELL_1p5.
    Comparing it against stratF_full_results_US1 changes damping twice over.

    SECOND DIFFERENCE, AND IT IS NOT DAMPING: that control runs LARGE-STRAIN
    ON from Part_I_MASON_v8.sav, while every Route 2 run so far has been
    small-strain from Part_I_MASON_v8_SmallStrain.sav. LARGE_STRAIN and
    BASE_SAVE are exposed below. Defaults keep small strain, so this run is
    a ONE-VARIABLE change from route2_full25_revised/ (damping only); set
    them to the control's values to make it a one-variable change from
    stratF_full_LS_MAXWELL_1p5 instead. You cannot have both at once --
    decide which comparison the figure is making before running.

    MAXWELL CAUTION, carried over verbatim from the Strategy F case: 3DEC's
    getMaxShearStiffness() returns a constant ks_ regardless of contact
    state, so if the Maxwell shear dashpot is applied to open subcontacts a
    separated joint transmits viscous shear c*ks*v_rel across the gap. Route
    2 is only 0.35 s of table per run against a full record, so the exposure
    is far smaller here than in Strategy F -- but the 2.5 s ring-down is
    where the wall sits cracked and open, and that is the window the settled
    displacement is read from. Treat the residual column with suspicion
    until energy-shear has been checked.

DAMPING CHOICE -- this is a base-driven model, which rules some options out:
  "maxwell"  (DEFAULT, and what the control uses) 3DEC's frequency-band
             damping, 1.5% over 1-40 Hz. Near-constant across the band, does
             not over-damp the slow rocking mode the way mass-proportional
             Rayleigh does, and does not touch the timestep.
  "rayleigh_stiff"  stiffness-proportional only. Correct physics for a model
             whose base moves, but it CUTS THE TIMESTEP hard (cost scales
             with 1/beta); expect the run to take several times longer.
  "rayleigh_mass"   mass-proportional. CHEAP AND WRONG HERE: it damps
             ABSOLUTE velocity, so while the table is driving the S and T_B
             groups it fights the input itself and it drags the whole wall
             during rigid-body rocking, which is exactly the motion being
             measured. Provided for comparison only -- do not report it.
  "local"    3DEC local damping, alpha = pi*xi. Amplitude-proportional, no
             timestep penalty, but it is not viscous and it also removes
             energy from rigid-body rocking. A blunt instrument; useful to
             bracket how sensitive the answer is to dissipation at all.
  "none"     reproduces route2_full25_revised_3dec.py exactly (sanity check).
Only "maxwell" matches the control. The other four are diagnostics; if one of
them is used, the Strategy F comparison is off the table.

The command actually issued is printed at startup and written to
OUT_DIR/damping.txt. There is no silent path: an unknown DAMP_MODE raises,
and there is no overrides file -- the constants in this file are the run.

OUTPUT: same layout and columns as the undamped driver, plus damping_mode /
damping_cmd / large_strain columns in the log, under
    route2_full25_MAXWELL_1p5_SS/            (defaults)
    route2_full25_MAXWELL_1p5_LS/            (LARGE_STRAIN = "on")
so nothing collides with the undamped results. Postprocess with
postprocess_route2_full.py (set ROUTE_DIR to that folder), or put both
folders in fig_chord_tilt_validation.py's MODELS list to see them together.
"""

import itasca as it
import os, csv, math, sys, json, glob
import numpy as np

# ============================ THE CASE ===============================
# ONE line selects the run. The 2x2 is geometry x damping:
#
#   CASE              geometry      damping        base save
#   SS_NODAMP         small-strain  none           Part_I_MASON_v8_SmallStrain.sav
#   SS_MAXWELL_1p5    small-strain  Maxwell 1.5%   Part_I_MASON_v8_SmallStrain.sav
#   LS_NODAMP         large-strain  none           part_I_mason_LS_strong.sav
#   LS_MAXWELL_1p5    large-strain  Maxwell 1.5%   part_I_mason_LS_strong.sav
#
# SS_NODAMP reproduces route2_full25_revised/ -- run it only to confirm this
# file and the old driver agree before trusting the other three. Each case
# writes to its own folder, so the four can run in any order, in parallel
# 3DEC sessions, and resume independently.
#
# Override from the shell without editing: ROUTE2_CASE=LS_NODAMP 3dec ...
CASE = "SS_NODAMP"
CASE = os.environ.get("ROUTE2_CASE", CASE).strip()

# Damping string verbatim from the Strategy F standalone case
# F_LS_MAXWELL_1p5. Do not retune it here: the point is that both sequences
# dissipate the same way. 1.5% over 1-40 Hz, three (weight, frequency) pairs.
MAXWELL_CMD = "block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 0.0104 40.0000"

# base save per geometry. The small-strain save is the one every Route 2 run
# so far has used. The large-strain save is the one named for this study.
SAVE_SS = "Part_I_MASON_v8_SmallStrain.sav"
SAVE_LS = "part_I_mason_LS_strong.sav"

# NOT THE SAME FILE AS THE STRATEGY F LARGE-STRAIN CONTROL, which starts from
# Part_I_MASON_v8.sav. If part_I_mason_LS_strong.sav differs by anything other
# than having been saved in large-strain mode -- the "_strong" in the name
# suggests joint properties may differ -- then LS Route 2 and
# stratF_full_LS_MAXWELL_1p5 differ by loading AND by model, and the pair is
# not a controlled comparison. Check cstav at t = 0 in run 1 of each before
# putting them on the same axes. The small-strain pair is unaffected.
CASES = {
    "SS_NODAMP":      dict(ls=False, damp="none",    save=SAVE_SS),
    "SS_MAXWELL_1p5": dict(ls=False, damp="maxwell", save=SAVE_SS),
    "LS_NODAMP":      dict(ls=True,  damp="none",    save=SAVE_LS),
    "LS_MAXWELL_1p5": dict(ls=True,  damp="maxwell", save=SAVE_LS),
}
if CASE not in CASES:
    raise RuntimeError("CASE must be one of {} (got '{}')".format(", ".join(sorted(CASES)), CASE))
_C = CASES[CASE]
LS_ON, DAMP_MODE, BASE_SAVE = _C["ls"], _C["damp"], _C["save"]
LARGE_STRAIN = "on" if LS_ON else "off"
DAMP_RATIO = 0.0          # Rayleigh fraction of critical -- 0.0 in the control.
                          # Read ONLY by the rayleigh_* / local diagnostic modes,
                          # which the four cases above never select; "maxwell"
                          # ignores it (its 1.5% is set by MAXWELL_CMD's weights).

# ============================ CONFIG =================================
BILINEAR_JSON = "bilinear_idealisation/bilinear_periods.json"
BILINEAR_KEY  = "model"        # "model" (v3 pushover, prediction mode) or "fig13" (envelope)
STAGE_B_PHASE = "same"         # "same" or "inverted"
XI            = 0.05           # damping of the Sd TARGET spectrum -- unrelated to DAMP_RATIO
DELTA_T       = 0.005
TAIL_SEC      = 2.5
DRIVE_GROUPS  = ["S", "T_B"]
PGA_WARN_G    = 0.90
SPECTRA = {"HU12": "spectrum_HU12.csv", "EC40": "spectrum_EC40.csv", "FR76": "spectrum_FR76.csv"}
PROTOCOL = [
    ( 1, "HU12", 0.50), ( 2, "HU12", 0.75), ( 3, "EC40", 0.20),
    ( 4, "HU12", 1.00), ( 5, "HU12", 1.25), ( 6, "EC40", 0.30),
    ( 7, "HU12", 1.50), ( 8, "EC40", 0.40), ( 9, "HU12", 1.75),
    (10, "HU12", 2.00), (11, "EC40", 0.50), (12, "HU12", 2.25),
    (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00),
    (16, "HU12", 3.50), (17, "HU12", 4.00), (18, "HU12", 4.50),
    (19, "HU12", 5.00), (20, "HU12", 5.50), (21, "HU12", 6.00),
    (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75),
    (25, "FR76", 2.00),
]
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
if STAGE_B_PHASE not in ("same", "inverted"):
    raise RuntimeError("STAGE_B_PHASE must be 'same' or 'inverted'")
PHASE_SIGN = 1.0 if STAGE_B_PHASE == "same" else -1.0

_MODES = ("maxwell", "rayleigh_stiff", "rayleigh_mass", "local", "none")
if DAMP_MODE not in _MODES:
    raise RuntimeError("DAMP_MODE must be one of {}".format(", ".join(_MODES)))
CASE_ID = "R2_" + CASE
OUT_DIR = "route2_full25_" + CASE
os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "route2_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "route2_full_log.csv")
LOG_COLS = ["run", "record", "scale", "d_hist_mm", "two_stage", "T_A_s", "SdA_mm", "VampA_mps", "PGA_A_g",
            "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g", "gap_s", "cap_note", "peak_edp_mm",
            "join_step_g", "SdA_total_mm", "SdB_total_mm", "period_source", "edp_source",
            "damping_mode", "damping_cmd", "large_strain", "base_save", "case_id"]

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

# ============================ PERIODS (fixed) + SPECTRA ==============
if not os.path.isfile(BILINEAR_JSON):
    raise RuntimeError("run bilinear_idealisation.py first: {} not found".format(BILINEAR_JSON))
_bj = json.load(open(BILINEAR_JSON))
if BILINEAR_KEY not in _bj:
    raise RuntimeError("{} has no '{}' block (present: {})".format(
        BILINEAR_JSON, BILINEAR_KEY, ", ".join(k for k in _bj if k != "rule")))
_b = _bj[BILINEAR_KEY]
T1, T_EFF = float(_b["T1_s"]), float(_b["Teff_s"])
PERIOD_SOURCE = "bilinear/{} ({}; K1 {:.2f} kN/mm, mu {:.2f}, m_eff {:.0f})".format(
    BILINEAR_KEY, "+".join(_b["branches_used"]), _b["K1_kNmm"], _b["mu"], _bj["rule"]["m_eff"])
SPEC = {}
for k, f in SPECTRA.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing spectrum file: {}".format(f))
    d = np.genfromtxt(f, delimiter=",", skip_header=1); SPEC[k] = (d[:, 0], d[:, 1])
def sd_record(rec, T):
    Ts, Sd = SPEC[rec]; return float(np.interp(T, Ts, Sd))

# centre frequency for Rayleigh: geometric mean of the two stage frequencies,
# so neither cycle sits at the over/under-damped end of the Rayleigh curve.
F_CENTRE = math.sqrt((1.0 / T1) * (1.0 / T_EFF))
if DAMP_MODE == "maxwell":
    DAMP_CMD = MAXWELL_CMD
elif DAMP_MODE == "rayleigh_stiff":
    DAMP_CMD = "block mechanical damping rayleigh {:.6f} {:.4f} stiffness".format(DAMP_RATIO, F_CENTRE)
elif DAMP_MODE == "rayleigh_mass":
    DAMP_CMD = "block mechanical damping rayleigh {:.6f} {:.4f} mass".format(DAMP_RATIO, F_CENTRE)
elif DAMP_MODE == "local":
    DAMP_CMD = "block mechanical damping local {:.6f}".format(math.pi * DAMP_RATIO)
else:
    DAMP_CMD = ""

# The record (Strategy F) control that pairs with this case. Damping and
# geometry must match for the pair to isolate LOADING.
#   F_SS_MAXWELL_1p5  being generated from the standalone driver with
#                     LARGE_STRAIN="off", START_SAVE="Part_I_MASON_v8_SmallStrain.sav",
#                     OUT_DIR='stratF_full_SS_MAXWELL_1p5'
#   F_LS_MAXWELL_1p5  exists; LARGE_STRAIN="on", START_SAVE="Part_I_MASON_v8.sav"
# A NODAMP case has no record control at all unless an undamped Strategy F is
# also run at the same geometry -- stratF_full_results_US1 is undamped but its
# geometry and save are not documented here, so it is not assumed to match.
CONTROLS = {
    "SS_MAXWELL_1p5": ("stratF_full_SS_MAXWELL_1p5", "Part_I_MASON_v8_SmallStrain.sav"),
    "LS_MAXWELL_1p5": ("stratF_full_LS_MAXWELL_1p5", "Part_I_MASON_v8.sav"),
}
CONTROL_DIR, CONTROL_SAVE = CONTROLS.get(CASE, (None, None))
MATCHES_CONTROL = CONTROL_DIR is not None and BASE_SAVE == CONTROL_SAVE
# the other axis: the same geometry with damping removed. This is always a
# sibling case of this file, so it always exists once both have been run.
TWIN_CASE = CASE.replace("MAXWELL_1p5", "NODAMP") if "MAXWELL" in CASE else CASE.replace("NODAMP", "MAXWELL_1p5")
TWIN_DIR  = "route2_full25_" + TWIN_CASE
# and the other geometry at the same damping
GEOM_TWIN = ("SS_" if LS_ON else "LS_") + CASE.split("_", 1)[1]
GEOM_DIR  = "route2_full25_" + GEOM_TWIN

print("=" * 72)
print("REVISED ROUTE 2, full protocol   [case {}]".format(CASE_ID))
print("  T1 {:.4f} s ({:.2f} Hz), T_eff {:.4f} s ({:.2f} Hz) on every run  [{}]".format(
    T1, 1.0 / T1, T_EFF, 1.0 / T_EFF, PERIOD_SOURCE))
print("  stage B phase : {}".format(STAGE_B_PHASE))
print("  DAMP_MODE     : {}".format(DAMP_MODE))
print("  command       : {}".format(DAMP_CMD if DAMP_CMD else "(none -- contact dissipation only)"))
print("  geometry      : large-strain {}".format("on" if LS_ON else "off"))
print("  base save     : {}".format(BASE_SAVE))
if DAMP_MODE in ("rayleigh_stiff", "rayleigh_mass"):
    print("  centre freq   : {:.3f} Hz (geometric mean of {:.2f} and {:.2f} Hz)".format(F_CENTRE, 1.0 / T1, 1.0 / T_EFF))
if DAMP_MODE == "rayleigh_mass":
    print("  ** mass-proportional damping resists ABSOLUTE velocity: it fights the")
    print("     driven table motion and drags rigid-body rocking. Diagnostic only.")
if DAMP_MODE == "rayleigh_stiff":
    print("  ** stiffness-proportional damping reduces the critical timestep; this run")
    print("     will be substantially slower than the undamped one.")
print("  output        : {}".format(OUT_DIR))
print("  " + "-" * 68)
print("  WHICH COMPARISON THIS RUN IS VALID FOR:")
print("    vs {}  -- damping is the only difference".format(TWIN_DIR))
print("    vs {}  -- geometry is the only difference".format(GEOM_DIR))
if CONTROL_DIR is None:
    print("    record control: NONE. {} is a damped Strategy F case; an".format(
        CONTROLS.get(CASE.replace("NODAMP", "MAXWELL_1p5"), ("stratF_full_*",))[0]))
    print("       undamped Route 2 has no matching record run unless one is made.")
elif MATCHES_CONTROL:
    print("    vs {}  -- LOADING is the only difference".format(CONTROL_DIR))
    print("       (damping and geometry identical; that control must exist)")
else:
    print("    vs {}  -- ** SUSPECT **".format(CONTROL_DIR))
    print("       damping and geometry match, but the base save does not:")
    print("         this run  {}".format(BASE_SAVE))
    print("         control   {}".format(CONTROL_SAVE))
    print("       If those two saves differ by more than large-strain mode -- the")
    print("       '_strong' in the name suggests joint properties may -- then this")
    print("       pair differs by LOADING AND MODEL and proves nothing. Check cstav")
    print("       at t=0 in run 1 of each before plotting them together.")
print("    stratF_full_results_US1 is the UNDAMPED record sequence. Its geometry and")
print("    base save are not pinned here; do not treat it as a control.")
print("  " + "-" * 68)
if not LS_ON:
    print("  SMALL-STRAIN CAVEAT AT THE RUNS YOU CARE ABOUT: block positions are not")
    print("  updated, so the restoring lever arm does not shorten as the wall leans.")
    print("  At 30 mm on a 210 mm wall that is ~14% of the thickness -- P-delta is")
    print("  left out, and left out the SAME way in every small-strain run (all")
    print("  over-predict stability late). Comparisons stay internally valid; an")
    print("  absolute run-24 capacity claim does not. That is what the LS_* cases")
    print("  are for.")
else:
    print("  LARGE-STRAIN NOTE: contacts CREATED during the run inherit the material-")
    print("  table default, not the mason assignment made at build time. If the base")
    print("  save was built without that default, new toe contacts are ungoverned.")
    print("  Expect this case to be slower and to need more restarts than the SS pair.")
print("=" * 72)
with open(os.path.join(OUT_DIR, "case.txt"), "w") as f:
    f.write("CASE={}\nCASE_ID={}\nDAMP_MODE={}\nDAMP_RATIO={}\nF_CENTRE_Hz={:.6f}\ncommand={}\n"
            "LARGE_STRAIN={}\nBASE_SAVE={}\ndamping_twin={}\ngeometry_twin={}\n"
            "record_control={}\nrecord_control_save={}\nmatches_control={}\n".format(
                CASE, CASE_ID, DAMP_MODE, DAMP_RATIO, F_CENTRE, DAMP_CMD,
                "on" if LS_ON else "off", BASE_SAVE, TWIN_DIR, GEOM_DIR,
                CONTROL_DIR, CONTROL_SAVE, MATCHES_CONTROL))

# ============================ SIGNAL BUILDERS =======================
def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")
def grid_period(T):
    spc = max(1, int(round(T / DELTA_T))); return spc * DELTA_T, spc
def one_cycle(Vamp, T):
    Tg, spc = grid_period(T); tt = np.arange(spc + 1) * DELTA_T
    return tt, Vamp * np.sin(2.0 * math.pi * tt / Tg), Tg
def calibrate_vel(T, Sd_target):
    _, v, Tg = one_cycle(1.0, T)
    sd0 = float(sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [Tg], XI)[0])
    if sd0 <= 0:
        raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg
def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81)
def build_signal(rec, scale):
    SdA, SdB = scale * sd_record(rec, T1), scale * sd_record(rec, T_EFF)
    VA, TgA = calibrate_vel(T1, SdA); VB, TgB = calibrate_vel(T_EFF, SdB)
    _, vA, _ = one_cycle(VA, TgA); _, vB, _ = one_cycle(PHASE_SIGN * VB, TgB)
    v = np.r_[vA, vB[1:]]; t = np.arange(len(v)) * DELTA_T
    step = abs((vB[1] - vB[0]) - (vA[-1] - vA[-2])) / DELTA_T / 9.81
    tot = sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [TgA, TgB], XI)
    return dict(t=t, v=v, TgA=TgA, TgB=TgB, SdA=SdA, SdB=SdB, VA=VA, VB=VB,
                pgaA=pga_g(vA), pgaB=pga_g(vB), step=step, SdA_tot=float(tot[0]), SdB_tot=float(tot[1]))

# ============================ MODEL SETUP ============================
def setup_dynamic(save_file):
    it.command("model restore '{}'".format(cmd_path(save_file).replace(".sav", "")))
    it.command("python-reset-state false")
    if LS_ON:
        it.command("model large-strain on")
        print("  GEOMETRY: large-strain ON (positions updated, contacts re-detected).")
        print("           NB contacts CREATED during the run inherit the material-table")
        print("           default, not the mason assignment made at build time.")
    else:
        it.command("model large-strain off")
        print("  GEOMETRY: small-strain")
    it.command("model dynamic active on")
    # zero first, unconditionally: a restored save can carry damping from
    # whatever wrote it (the pushover driver saves with local 0.9).
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    if DAMP_CMD:
        it.command(DAMP_CMD)
        print("  damping applied: {}".format(DAMP_CMD))
    else:
        print("  damping: none (contact dissipation only)")
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

def run_label(n, rec, sc):
    return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def save_path(n):
    return os.path.join(OUT_DIR, "route2_run_{:02d}.sav".format(n))

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
    x = find_ch(folder, KEY_DISP_EXP)
    if x is not None:
        return float(np.max(np.abs(x - x[0]))), KEY_DISP_EXP
    c3, c4, c5 = find_ch(folder, CH3), find_ch(folder, CH4), find_ch(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, "n/a"
    n = min(len(c3), len(c4), len(c5)); rel = 0.5 * (c3[:n] + c4[:n]) - c5[:n]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0, "0.5(Ch3+Ch4)-Ch5"
def resid_edp_mm(folder, tail_frac=0.05):
    """settled displacement at the end of the tail -- only meaningful WITH damping."""
    x = find_ch(folder, KEY_DISP_EXP)
    if x is None:
        return None
    n = max(1, int(tail_frac * len(x)))
    return float(np.mean(x[-n:]) - x[0])

def load_ckpt():
    if os.path.isfile(CKPT):
        with open(CKPT) as f:
            s = json.load(f)
        return s["last_run"], s["summary"]
    return 0, []
def save_ckpt(last, summary):
    with open(CKPT, "w") as f:
        json.dump({"last_run": last, "summary": summary, "damping": DAMP_CMD}, f, indent=1)

# ============================ CHECKPOINT =============================
it.command("python-reset-state false")
it.command("program automatic-model-save active off")
last_done, summary = load_ckpt()
if last_done == 0:
    if not os.path.isfile(BASE_SAVE):
        raise RuntimeError("missing " + BASE_SAVE)
    setup_dynamic(BASE_SAVE); print("fresh start from", BASE_SAVE)
else:
    setup_dynamic(save_path(last_done)); print("resuming after run", last_done)

log_new = (last_done == 0) or not os.path.isfile(LOG)
logf = open(LOG, "w" if log_new else "a", newline="")
logw = csv.writer(logf)
if log_new:
    logw.writerow(LOG_COLS)

# ============================ MAIN LOOP ==============================
for run_no, rec, sc in PROTOCOL[last_done:]:
    d_hist = max([r["peak_edp_mm"] for r in summary if r.get("peak_edp_mm") is not None] + [0.0])   # logged only, not used
    S = build_signal(rec, sc)
    lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(S["t"])))
        for a_, b_ in zip(S["t"], S["v"]):
            f.write("{:.6f}\t{:.9e}\n".format(a_, b_))
    pulse_dur = float(S["t"][-1])

    print("\n" + "=" * 66)
    print("Run {:02d}  {} x {:.2f}   [damping: {}]".format(run_no, rec, sc, DAMP_MODE))
    print("  stage A: T {:.3f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.2f}g".format(S["TgA"], S["SdA"] * 1e3, S["VA"], S["pgaA"]))
    print("  stage B: T {:.3f}s  Sd {:.2f}mm  Vamp {:.4f} m/s  PGA {:.2f}g   join step {:.2f} g".format(S["TgB"], S["SdB"] * 1e3, S["VB"], S["pgaB"], S["step"]))
    print("  whole signal: Sd(T1) {:.2f} mm, Sd(T_eff) {:.2f} mm ; table {:.3f} s + tail {:.1f} s".format(S["SdA_tot"] * 1e3, S["SdB_tot"] * 1e3, pulse_dur, TAIL_SEC))
    for tag, p in (("A", S["pgaA"]), ("B", S["pgaB"])):
        if p > PGA_WARN_G:
            print("  ** stage {} PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(tag, p, PGA_WARN_G))
    print("=" * 66)

    tbl = "r2d_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(pulse_dur))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("model save '{}'".format(cmd_path(save_path(run_no))))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(folder)))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call 'instrument_history_export_v2.dat'")
    peak, src = peak_edp_mm(folder)
    resid = resid_edp_mm(folder)
    print("  PEAK EDP ({}) = {}{}".format(
        src, "{:.2f} mm".format(peak) if peak is not None else "n/a",
        "   settled {:+.2f} mm".format(resid) if resid is not None else ""))

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")

    row = dict(run=run_no, record=rec, scale=sc, d_hist_mm=round(d_hist, 3), two_stage=1,
               T_A_s=round(S["TgA"], 5), SdA_mm=round(S["SdA"] * 1e3, 3), VampA_mps=round(S["VA"], 5), PGA_A_g=round(S["pgaA"], 4),
               T_B_s=round(S["TgB"], 5), SdB_mm=round(S["SdB"] * 1e3, 3), VampB_mps=round(S["VB"], 5), PGA_B_g=round(S["pgaB"], 4),
               gap_s=0.0, cap_note="fixed_periods", peak_edp_mm=round(peak, 3) if peak is not None else None,
               join_step_g=round(S["step"], 3), SdA_total_mm=round(S["SdA_tot"] * 1e3, 3), SdB_total_mm=round(S["SdB_tot"] * 1e3, 3),
               period_source=PERIOD_SOURCE.replace(",", ";"), edp_source=src,
               damping_mode=DAMP_MODE, damping_cmd=DAMP_CMD.replace(",", ";"),
               large_strain="on" if LS_ON else "off", base_save=BASE_SAVE, case_id=CASE_ID)
    summary.append(row)
    logw.writerow([row[k] if row[k] is not None else "" for k in LOG_COLS])
    logf.flush()
    save_ckpt(run_no, summary)

logf.close()
print("\ndone: {} runs in {}  [case {}].".format(len(summary), OUT_DIR, CASE_ID))
print("Postprocess with postprocess_route2_full.py (ROUTE_DIR = '{}'), or add".format(OUT_DIR))
print("  ('3DEC Route 2, {}', '{}', '#8a6a4a', 'v')".format(CASE.replace("_", " "), OUT_DIR))
print("to MODELS in fig_chord_tilt_validation.py.")
print("Remaining cases: {}".format(", ".join(
    c for c in ("SS_NODAMP", "SS_MAXWELL_1p5", "LS_NODAMP", "LS_MAXWELL_1p5")
    if not os.path.isfile(os.path.join("route2_full25_" + c, "case.txt"))) or "none -- all four done"))