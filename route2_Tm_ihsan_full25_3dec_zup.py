# -*- coding: ascii -*-
"""
route2_Tm_ihsan_full25_3dec_zup.py -- route2_Tm_ihsan_full25_3dec.py for the Z-UP frame (29 Oct 2026).
  Frame (Geo_Prep_zup*.dat / ANALYSIS_PART_I_MASON_zup*.dat): X = 1.292 - x in-plane, Y = z out-of-plane (+Y towards
  the joists), Z = y vertical. Changes against the old-frame driver, nothing else:
    * BUILD picks the save, the pushover, the joist-support contact ranges and the joist-probe heights:
        "zup_pk2" : Part_I_MASON_v8_LargeStrain_springforce_zup_pk2 (pockets at the test heights, f_t = 0)
        "zup"     : Part_I_MASON_v8_LargeStrain_springforce_zup     (old pocket heights, f_t = 0)
        "zup_pk2_char" : Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char (test pocket heights, characterisation
                    strengths on every joint: tension 0.344 MPa, cohesion 0.75 MPa; the material of the char Strategy F)
    * table and top-beam drive along Y (was the old z); S freed in Y
    * spring-force check: the vertical (z) component of the T_B applied force (was y)
    * Joist_S1/S2_contact ranges: pos-z (height) / pos-y (OOP)
    * instrumentation: instrument_history_new_zup.dat (calls instrument_tilt_v2_zup.dat), instrument_baseshear_exp_zup.dat;
      instrument_history_export_v2.dat unchanged (channel names and indices are the same)
    * frame guard: stops on an old-frame save (top beam must sit at Z > 2.5 m)
  The signal, the targets, the fit, the gate, the checkpoints and the outputs are those of the old-frame driver, so
  the results compare channel by channel with r2tm_TpTm_full25_LS_maxwell (old frame, old pockets, f_t = 0).
  T1 and T_eff come from the z-up pushover of the same BUILD (PUSHOVER_CSV / _NEG below), same recipe as before
  (BACKBONE_COMBINE "mean", F_SOURCE "column", D_U_MM 27.4).

route2_Tm_ihsan_full25_3dec.py -- FULL PROTOCOL (runs 1-25). Route 2 in Ihsan's form: the AMPLITUDE targets come
from the wall (T1, T_eff), the PERIOD of both sine waves comes from the RECORD. The record period is Rathje's mean
period T_m instead of the peak-Sa period T_p. LARGE STRAIN, MAXWELL damping (MAXWELL_CMD).
BUILD INTEGRITY: joint properties, densities, spring force and boundary conditions are those of BASE_SAVE, i.e. of
    ANALYSIS_PART_I_MASON(...).dat + mat_case.dat. This driver assigns NO joint property. After setup it audits the
    masonry-masonry subcontacts against the build's own FISH values (tension, cohesion, kn, friction) and stops on a
    mismatch. The only boundary changes are the ones driving needs (S freed in z and driven, T_B driven in z);
    T_B's rotation-y/z fixity is kept unless FREE_TB_ROTATIONS = True (the pre-23-Oct behaviour).

WHY T_m AND NOT T_p (preflight of 20 Oct, route2_Tp_preflight/): with one cycle at T_p the wave that has to deliver
    Sd(T_eff) is much shorter than T_eff for HU12 (0.121 vs 0.357 s). A long-period oscillator then only follows the
    ground displacement of the pulse, so the acceleration needed grows as Sd / T_p^2: wave B needed 3.6 x the
    record's PGA on every HU12 run (2.85 g at run 21), and because 0.121 s sits next to T1 it also drove Sd(T1) to
    5 x wave A's target. T_m (Rathje, Abrahamson & Bray 1998; Rathje, Faraj, Russell & Bray 2004) is the Fourier-
    amplitude-weighted mean period of the record, a standard, wall-independent frequency-content parameter, and it
    sits closer to T_eff:   HU12 0.206 s, EC40 0.544 s, FR76 0.374 s   (T_p: 0.121, 0.251, 0.261 s).
        T_m = sum(C_i^2 / f_i) / sum(C_i^2),  C_i = Fourier amplitudes of the acceleration, 0.25 <= f_i <= 20 Hz
    It is computed here from the record velocity files (vel_*.txt), not hard-coded.

THE SIGNAL, every run
    wave A : CYCLES_A sine cycles of table VELOCITY at PERIOD_A (default T_p of the record)
    wave B : CYCLES_B sine cycles at PERIOD_B (default T_m of the record), phase PHASE_B (+1 same, -1 inverted)
    targets (XI-damped Sd of the WHOLE signal, a 2 s free tail included):
             Sd(T1)    = s * Sd_record(T1)          wave A's job
             Sd(T_eff) = s * Sd_record(T_eff)       wave B's job
    The two amplitudes are solved JOINTLY on the combined signal (Newton on log-amplitudes, as in the fixed-point
    driver). If one wave alone already meets or exceeds the other wave's target, the run collapses to that single
    wave and the overshoot of the other target is logged (SdA_over_tgtA / SdB_over_tgtB).
    PERIOD_A / PERIOD_B accept "Tm", "Tp", "T1", "Teff" or a number in s, so Ihsan's original T_p form ("Tp", "Tp"),
    the all-T_m form ("Tm", "Tm") and the original Route 2 ("T1", "Teff") run from the same file for comparison.

WHY ("Tp", "Tm") IS THE DEFAULT (preflight, 22 Oct, route2_Tm_preflight/):
    ("Tm","Tm")  two cycles at the same period cannot meet two independent targets -- the pulse's Sd(T_eff)/Sd(T1)
                 ratio is fixed by T_m -- so the fit drops wave B; one T_m cycle meets Sd(T1) and overshoots Sd(T_eff)
                 by 14 % (HU12), 33 % (FR76), 76 % (EC40). 8 runs flagged.
    ("Tp","Tp")  wave B alone overshoots Sd(T1) 4x and needs 3x the record PGA on HU12. 21 runs flagged.
    ("Tp","Tm")  both periods from the record, two different spectral shapes: both targets met exactly on all 25
                 runs, PGA 0.96-1.04 x the scaled record, max 0.74 g (run 21). 0 runs flagged.

T1, T_EFF (fixed for all 25 runs, the point of the method)
    from the pushover in PUSHOVER_CSV with the measured effective mass L (M_EFF_KG = 1200 kg, L = sum m_j phi_j from
    the walls' free-vibration mode shape; exp_effective_mass.py):
        K_el  = least-squares slope through the origin over K_EL_RANGE_MM
        F_y   = peak force of the curve (before D_CUT_MM if given),   d_u = D_U_MM
        T1    = 2 pi sqrt(L / K_el)          T_eff = 2 pi sqrt(L d_u / F_y)  (= T1 sqrt(K_el d_u / F_y))
    D_U_MM = 27.4 mm is the end of the US-1 test envelope, NOT a model output: the ft0 gravity-tilt push stops at
    19.4 mm (joist out at 17.4). T_eff moves with it (16.4 mm -> 0.276 s, 27.4 -> 0.357, 76 -> 0.565 on v5).
    T1_S_OVERRIDE / TEFF_S_OVERRIDE bypass the computation.

  TWO BRANCHES (PUSHOVER_CSV = positive push, PUSHOVER_CSV_NEG = negative push; signs are irrelevant, |d| and |F|
    are used). BACKBONE_COMBINE decides how they enter -- the same recipe as the test envelopes
    (t1_teff_envelopes.py), so model and test periods stay comparable:
        "mean" : K_el = least-squares slope over K_EL_RANGE_MM with the points of BOTH branches pooled;
                 F_y, D_Y = mean of the two branch peaks (each before D_CUT_MM); T1, T_eff from those
        "pos"  : positive branch only (the behaviour before 22 Oct)
        "neg"  : negative branch only
    Whatever the mode, T1 / T_eff of each branch alone are printed and kept in the checkpoint, so the asymmetry is
    visible. With no PUSHOVER_CSV_NEG the mode falls back to "pos".
  FORCE SOURCE (F_SOURCE): "column" reads F_COL (F_tot_kN, the sum of the support reactions); "applied" uses
    |a_mps2| x M_PUSH_KG, the load the servo actually applies. On the negative ft0 masonry push (22 Oct) F_tot
    exceeds the applied load by a term growing ~0.03 kN/mm after the runaway (F_tot/a = 1308 kg on the ascent,
    1432 kg at 31 mm): the applied load peaks at 20.8 mm and softens, F_tot keeps creeping up. Until that excess
    is explained, "applied" is the safer backbone for a masonry-only push (M_PUSH_KG = masonry mass, 1308 kg).
    For a gravity-tilt push the applied load acts on every block -- set the mass to the total mass, or use "column".
    M_PUSH_KG / M_PUSH_KG_NEG are per branch, because the two pushes need not be of the same kind: on 22 Oct the
    positive file was a gravity-tilt push (F_tot / a = 1750 kg) and the negative one a masonry-only push (1308 kg).
    Check: on the elastic ascent F_tot / |a| equals the mass the push loads.

FEASIBILITY GATE (preflight, every run)
    'pga'    : PGA of the signal > PGA_RATIO_MAX x (s x record PGA)
    'table'  : PGA of the signal > PGA_WARN_G (the table reached 0.78 g)
    'leakA'  : wave B alone overshoots Sd(T1) target by more than LEAK_MAX (the run collapses to wave B)
    'leakB'  : wave A alone overshoots Sd(T_eff) target by more than LEAK_MAX (the run collapses to wave A)
    GATE = "flag" runs everything and logs the flags; "stop" refuses to start 3DEC if any run is flagged.

PREFLIGHT: set PREFLIGHT_ONLY = True and run with plain python (no 3DEC): targets, amplitudes, PGA/PGV, flags go to
    OUT_DIR/r2tm_targets.csv and the velocity tables are written.

OUTPUT   OUT_DIR/RunNN_REC_sXpYY/ (channel CSVs), r2tm_run_NN.sav, RunNN_..._vel.txt, r2tm_targets.csv,
         r2tm_log.csv, r2tm_checkpoint.json (delete to start the folder over)
REQUIRES next to this script: PUSHOVER_CSV and PUSHOVER_CSV_NEG (unless both overrides are set), spectrum_HU12/EC40/FR76.csv,
    vel_HU.txt / vel_EC.txt / vel_FR.txt, the base save, instrument_history_new.dat, instrument_baseshear_exp.dat,
    instrument_history_export_v2.dat
"""
import os, csv, math, json, glob
import numpy as np
try:
    import itasca as it
except ImportError:
    it = None

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False            # True: targets + tables only, no 3DEC

BUILD        = "zup_pk2"          # "zup_pk2" | "zup_pk2_char" | "zup"   (see the header)
_BUILDS = {   # save, pushover folder, joist bottoms Z (support contacts), joist centres Z (slip probes)
    "zup_pk2": ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2.sav", "pushover_results_v5zup_pk2_ft0_TopMason",
                (0.23455, 2.11091), (0.32250, 2.19886)),
    "zup":     ("Part_I_MASON_v8_LargeStrain_springforce_zup.sav", "pushover_results_v5zup_ft0_TopMason",
                (0.29318, 2.05227), (0.38114, 2.14023)),
    "zup_pk2_char": ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char.sav", "pushover_results_v5zup_pk2_char_TopMason",
                (0.23455, 2.11091), (0.32250, 2.19886)),
}
OUT_DIR      = "r2tm_TpTm_full25_LS_maxwell_" + BUILD
RESUME_AFTER = 0
RESUME_SAVE  = ""
STOP_AFTER   = 25
BASE_SAVE    = _BUILDS[BUILD][0]
JOIST_BOTTOM_Z, JOIST_CENTRE_Z = _BUILDS[BUILD][2], _BUILDS[BUILD][3]

# ---- wall periods ----------------------------------------------------------------------
PUSHOVER_CSV     = _BUILDS[BUILD][1] + "/pushover_pos.csv"
PUSHOVER_CSV_NEG = _BUILDS[BUILD][1] + "/pushover_neg.csv"   # "" for positive only
BACKBONE_COMBINE = "mean"         # "mean" | "pos" | "neg"   (see the docstring)
D_COL, F_COL     = "d_ctrl_mm", "F_tot_kN"
F_SOURCE         = "column"       # "column": F_COL | "applied": |A_COL| x M_PUSH_KG
A_COL            = "a_mps2"
M_PUSH_KG        = 1308.0         # mass the positive push loads: masonry 1308 kg (masonry-only push) / all blocks (gravity tilt)
M_PUSH_KG_NEG    = None           # same for the negative push; None = M_PUSH_KG. Check: F_tot / a on the ascent = this mass
M_EFF_KG       = 1200.0           # L = sum m_j phi_j (US-1 1202, US-2 1195 kg)
K_EL_RANGE_MM  = (0.1, 1.5)
D_U_MM         = 27.4             # end of the US-1 envelope; see the docstring
D_CUT_MM       = 17.0             # F_y taken before this d: the gravity-tilt push loses the joist at 17.4 mm (harmless on v5: peak at 12 mm)
T1_S_OVERRIDE  = None             # s, or None to compute
TEFF_S_OVERRIDE = None            # s, or None to compute

# ---- signal ----------------------------------------------------------------------------
PERIOD_A   = "Tp"                 # "Tm" | "Tp" | "T1" | "Teff" | number (s)   -- default: record T_p (preflight 22 Oct)
PERIOD_B   = "Tm"                 #   default: record T_m. ("Tm","Tm") collapses to one wave; ("Tp","Tp") needs 3x PGA on HU12
CYCLES_A, CYCLES_B = 1, 1
PHASE_B    = +1.0
XI         = 0.05                 # damping of the target spectra AND of the signal calibration (same, deliberately)
TM_BAND_HZ = (0.25, 20.0)         # Rathje et al. (2004)
DELTA_T    = 0.0025
TAIL_SEC   = 3.0
FIT_TOL, FIT_MAXIT = 1e-8, 80
RECORD_VEL = {"HU12": "vel_HU.txt", "EC40": "vel_EC.txt", "FR76": "vel_FR.txt"}
SPECTRA    = {"HU12": "spectrum_HU12.csv", "EC40": "spectrum_EC40.csv", "FR76": "spectrum_FR76.csv"}

# ---- feasibility gate ------------------------------------------------------------------
GATE          = "flag"            # "flag" | "stop"
PGA_RATIO_MAX = 1.5
PGA_WARN_G    = 0.90
LEAK_MAX      = 1.25

# ---- model -----------------------------------------------------------------------------
LARGE_STRAIN   = True
MAXWELL_CMD    = "block mech damp maxwell 0.0112 1.0005 0.0081 7.4315 0.0096 25.0004"   # 1.5 % flat 1-25 Hz
# MAXWELL_CMD  = "block mech damp maxwell 0.0224 1.0005 0.0162 7.4315 0.0192 25.0004"   # 3 %  (fractions x 2, check on the free tails)
# MAXWELL_CMD  = "block mech damp maxwell 0.0373 1.0005 0.0270 7.4315 0.0320 25.0004"   # 5 %  (fractions x 3.33, check on the free tails)
DRIVE_GROUPS   = ["S", "T_B"]
SPRING_FORCE_N = 10350.0
# ---- build integrity: this driver never assigns joint properties -----------------------------------------------
VERIFY_BUILD      = True          # audit masonry-masonry subcontact properties against the build's FISH values; stop on mismatch
PROP_TOL          = 0.01          # relative tolerance of the audit
PROP_MISMATCH_MAX = 0.02          # fraction of subcontacts allowed off the build values (material-table defaults of contacts
                                  #   created after the build); above it the save is not the build you think it is
FREE_TB_ROTATIONS = False         # False: keep the build's T_B rotation-y/z fixity. True: free them (behaviour before 23 Oct)
INSTRUMENT_FILES = ["instrument_history_new_zup.dat", "instrument_baseshear_exp_zup.dat"]   # the first calls instrument_tilt_v2_zup.dat
EXPORT_FILE      = "instrument_history_export_v2.dat"
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
# EDP channel: rel_disp_top_mm = FISH history 16 of instrument_history_new.dat (2.06 m minus table). The previous key,
# rel_disp_top_exp_mm (history 22), left the wall's curve by 36-81 mm on runs 23-25 of the 23 Oct set.
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"
KEY_BS = "base_shear_exp_kN"
ACH_FFT_BAND_HZ = (1.0, 20.0)
# ======================================================================

os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "r2tm_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "r2tm_log.csv")
TARGETS_CSV = os.path.join(OUT_DIR, "r2tm_targets.csv")
TARGET_COLS = ["run", "record", "scale", "T1_s", "Teff_s", "Tm_s", "Tp_s", "flags", "stages",
               "target_A_mm", "target_B_mm", "T_A_s", "cyclesA", "VampA_mps", "PGA_A_g",
               "T_B_s", "cyclesB", "phaseB", "VampB_mps", "PGA_B_g", "PGA_tot_g", "PGV_tot_mps",
               "SdA_tot_mm", "SdB_tot_mm", "SdA_over_tgtA", "SdB_over_tgtB", "PGA_over_rec", "PGV_over_rec",
               "fit_iter", "fit_resid", "table_s"]
LOG_COLS = TARGET_COLS + ["peak_edp_mm", "resid_edp_mm", "peak_bs_exp_kN", "T_ach_halfcycle_s", "T_ach_fft_s",
                          "edp_source", "state_source"]

# ============================ Sd INTEGRATOR (self-tested) =============
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

# ============================ WALL PERIODS ============================
def read_backbone(path, m_push=None):
    """|d|, |F| of a pushover CSV ('#' lines skipped), monotone in |d| (unloading / repeated rows dropped).
    Works for either push direction: the signs of d and F are discarded. Rows with converged == 0 (the runaway
    step that ends the force-controlled ascent) are skipped."""
    rows = [r for r in open(path) if r.strip() and not r.lstrip().startswith("#")]
    d, F = [], []
    for r in csv.DictReader(rows):
        if str(r.get("converged", "1")).strip() in ("0", "0.0", "False", "false"):
            continue                                   # the runaway step of the force-controlled ascent: not an equilibrium
        try:
            dd = abs(float(r[D_COL]))
            ff = abs(float(r[A_COL])) * m_push / 1e3 if F_SOURCE == "applied" else abs(float(r[F_COL]))
        except (KeyError, ValueError, TypeError):
            continue
        if dd > 0 and (not d or dd > d[-1]):
            d.append(dd); F.append(ff)
    if len(d) < 5:
        raise RuntimeError("backbone {}: fewer than 5 usable rows (columns {} / {})".format(
            path, D_COL, A_COL if F_SOURCE == "applied" else F_COL))
    return np.array(d), np.array(F)
def branch_props(d, F):
    """elastic points (for the pooled fit), K_el of this branch alone, and the peak before D_CUT_MM."""
    m = (d >= K_EL_RANGE_MM[0]) & (d <= K_EL_RANGE_MM[1])
    if m.sum() < 2:
        m = np.zeros(len(d), bool); m[0] = True
    mc = d <= D_CUT_MM if D_CUT_MM else np.ones(len(d), bool)
    iy = int(np.argmax(F[mc]))
    return dict(d_el=d[m], F_el=F[m], K=float(np.sum(d[m] * F[m]) / np.sum(d[m] ** 2)),
                Fy=float(F[mc][iy]), Dy=float(d[mc][iy]), dmax=float(d[-1]))
def periods_of(K, Fy):
    T1 = 2 * math.pi * math.sqrt(M_EFF_KG / (K * 1e6))
    Te = 2 * math.pi * math.sqrt(M_EFF_KG * D_U_MM * 1e-3 / (Fy * 1e3))
    return T1, Te
BRANCHES = {}
if T1_S_OVERRIDE is not None and TEFF_S_OVERRIDE is not None:
    T1_S, TEFF_S = float(T1_S_OVERRIDE), float(TEFF_S_OVERRIDE)
    PERIOD_NOTE = "overrides: T1 {:.4f} s, T_eff {:.4f} s".format(T1_S, TEFF_S)
else:
    if BACKBONE_COMBINE not in ("mean", "pos", "neg"):
        raise RuntimeError("BACKBONE_COMBINE must be 'mean', 'pos' or 'neg'")
    for lab, path, mp in (("pos", PUSHOVER_CSV, M_PUSH_KG), ("neg", PUSHOVER_CSV_NEG, M_PUSH_KG_NEG if M_PUSH_KG_NEG else M_PUSH_KG)):
        if path:
            if not os.path.isfile(path):
                raise RuntimeError("missing pushover CSV ({} branch): {}".format(lab, path))
            BRANCHES[lab] = branch_props(*read_backbone(path, mp))
            BRANCHES[lab]["m_push"] = mp
            BRANCHES[lab]["T1"], BRANCHES[lab]["Teff"] = periods_of(BRANCHES[lab]["K"], BRANCHES[lab]["Fy"])
    mode = BACKBONE_COMBINE if (BACKBONE_COMBINE != "mean" or len(BRANCHES) == 2) else "pos"
    if mode != "mean" and mode not in BRANCHES:
        raise RuntimeError("BACKBONE_COMBINE='{}' but that branch was not given".format(BACKBONE_COMBINE))
    if mode == "mean":
        dd = np.concatenate([BRANCHES[b]["d_el"] for b in ("pos", "neg")])
        ff = np.concatenate([BRANCHES[b]["F_el"] for b in ("pos", "neg")])
        K_EL = float(np.sum(dd * ff) / np.sum(dd * dd))
        F_Y = 0.5 * (BRANCHES["pos"]["Fy"] + BRANCHES["neg"]["Fy"]); D_Y = 0.5 * (BRANCHES["pos"]["Dy"] + BRANCHES["neg"]["Dy"])
    else:
        K_EL, F_Y, D_Y = BRANCHES[mode]["K"], BRANCHES[mode]["Fy"], BRANCHES[mode]["Dy"]
    T1c, Tec = periods_of(K_EL, F_Y)
    T1_S = float(T1_S_OVERRIDE) if T1_S_OVERRIDE is not None else T1c
    TEFF_S = float(TEFF_S_OVERRIDE) if TEFF_S_OVERRIDE is not None else Tec
    PERIOD_NOTE = "{} [{}{}]: K_el {:.3f} kN/mm, F_y {:.2f} kN at {:.1f} mm, d_u {:.1f} mm, L {:.0f} kg -> T1 {:.4f} s, T_eff {:.4f} s (mu {:.2f})".format(
        "+".join(sorted(BRANCHES)), mode, ", applied a*m_push" if F_SOURCE == "applied" else "",
        K_EL, F_Y, D_Y, D_U_MM, M_EFF_KG, T1_S, TEFF_S, D_U_MM / D_Y)
    for lab in sorted(BRANCHES):
        b = BRANCHES[lab]
        print("  branch {}: K_el {:.3f} kN/mm, F_y {:.2f} kN at {:.1f} mm, curve to {:.1f} mm -> T1 {:.4f} s, T_eff {:.4f} s{}".format(
            lab, b["K"], b["Fy"], b["Dy"], b["dmax"], b["T1"], b["Teff"],
            "   [force = |a| x {:.0f} kg]".format(b["m_push"]) if F_SOURCE == "applied" else ""))
        if D_U_MM > b["dmax"]:
            print("    (d_u {:.1f} mm is beyond this branch's last point {:.1f} mm: F_y is taken as the plateau)".format(D_U_MM, b["dmax"]))

# ============================ RECORDS: spectra, T_m, T_p, PGA/PGV =====
SPEC, REC = {}, {}
for k, f in SPECTRA.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing spectrum file: {}".format(f))
    s = np.genfromtxt(f, delimiter=",", skip_header=1); SPEC[k] = (s[:, 0], s[:, 1] * 1e3)   # mm
def sd_record_mm(rec, T):
    Ts, Sd = SPEC[rec]; return float(np.interp(T, Ts, Sd))
for k, f in RECORD_VEL.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing record velocity file (needed for T_m): {}".format(f))
    d = np.loadtxt(f, skiprows=2); dtr = float(d[1, 0] - d[0, 0]); ag = np.gradient(d[:, 1], dtr)
    C = np.abs(np.fft.rfft(ag)); fr = np.fft.rfftfreq(len(ag), dtr)
    b = (fr >= TM_BAND_HZ[0]) & (fr <= TM_BAND_HZ[1])
    Tm = float(np.sum(C[b] ** 2 / fr[b]) / np.sum(C[b] ** 2))
    Ts, Sd = SPEC[k]; Sa = Sd * (2 * np.pi / Ts) ** 2; Tp = float(Ts[int(np.argmax(Sa))])
    REC[k] = dict(Tm=Tm, Tp=Tp, pga=float(np.max(np.abs(ag)) / 9.81), pgv=float(np.max(np.abs(d[:, 1]))))
def carrier(spec, rec):
    if isinstance(spec, (int, float)):
        return float(spec)
    return {"Tm": REC[rec]["Tm"], "Tp": REC[rec]["Tp"], "T1": T1_S, "Teff": TEFF_S}[spec]

print("=" * 78)
print("ROUTE 2, Ihsan's form with record periods   [{}]".format(OUT_DIR))
print("  wall : " + PERIOD_NOTE)
for k in sorted(REC):
    print("  {}  T_m {:.3f} s  T_p {:.3f} s  PGA {:.3f} g  PGV {:.3f} m/s   wave A at {:.3f} s, wave B at {:.3f} s".format(
        k, REC[k]["Tm"], REC[k]["Tp"], REC[k]["pga"], REC[k]["pgv"], carrier(PERIOD_A, k), carrier(PERIOD_B, k)))
print("  targets: Sd(T1 {:.4f}) and Sd(T_eff {:.4f}) at xi {:.2f}; PERIOD_A {}, PERIOD_B {}, cycles {}/{}, phase B {:+.0f}".format(
    T1_S, TEFF_S, XI, PERIOD_A, PERIOD_B, CYCLES_A, CYCLES_B, PHASE_B))
print("=" * 78)

# ============================ SIGNAL BUILDERS =======================
def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")
def grid_period(T):
    spc = max(1, int(round(T / DELTA_T))); return spc * DELTA_T, spc
def cycles(Vamp, T, n):
    Tg, spc = grid_period(T); tt = np.arange(n * spc + 1) * DELTA_T
    return Vamp * np.sin(2.0 * math.pi * tt / Tg), Tg
def join(vA, vB):
    return np.r_[vA, vB[1:]] if len(vA) else vB
_TAIL = np.zeros(int(round(2.0 / DELTA_T)))
def sd_signal(v, Ts):
    return sd_spectrum(np.gradient(np.r_[v, _TAIL], DELTA_T), DELTA_T, Ts, XI)
def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81) if len(v) > 1 else 0.0

def build_signal(rec, scale):
    Ts = [T1_S, TEFF_S]
    tA = scale * sd_record_mm(rec, T1_S); tB = scale * sd_record_mm(rec, TEFF_S)
    uA, TgA = cycles(1.0, carrier(PERIOD_A, rec), CYCLES_A)
    uB, TgB = cycles(PHASE_B, carrier(PERIOD_B, rec), CYCLES_B)
    sA = sd_signal(uA, Ts) * 1e3; sB = sd_signal(uB, Ts) * 1e3          # mm per unit amplitude, at (T1, T_eff)
    VA0 = tA / sA[0]; VB0 = tB / sB[1]
    stages, it_n, resid = "A+B", 0, 0.0
    if VA0 * sA[1] >= tB * (1 - 1e-6):          # wave A alone covers the T_eff target
        VA, VB, stages = VA0, 0.0, "A"
    elif VB0 * sB[0] >= tA * (1 - 1e-6):        # wave B alone covers the T1 target
        VA, VB, stages = 0.0, VB0, "B"
    else:
        x = np.log([VA0, VB0])
        def resid_of(x):
            s = sd_signal(join(np.exp(x[0]) * uA, np.exp(x[1]) * uB), Ts) * 1e3
            return np.array([s[0] / tA - 1.0, s[1] / tB - 1.0])
        r = resid_of(x)
        for it_n in range(1, FIT_MAXIT + 1):
            if np.max(np.abs(r)) < FIT_TOL:
                break
            J = np.zeros((2, 2)); h = 1e-4
            for k in range(2):
                xp = x.copy(); xp[k] += h; J[:, k] = (resid_of(xp) - r) / h
            try:
                dx = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                dx = -r
            dx = np.clip(dx, -1.0, 1.0); lam = 1.0
            for _ in range(8):
                xn = x + lam * dx; rn = resid_of(xn)
                if np.max(np.abs(rn)) < np.max(np.abs(r)):
                    break
                lam *= 0.5
            x, r = xn, rn
        resid = float(np.max(np.abs(r)))
        if resid > 1e-3:
            print("  ** joint fit did not converge (resid {:.2e}); amplitudes are the last iterate".format(resid))
        VA, VB = float(np.exp(x[0])), float(np.exp(x[1]))
    vA = VA * uA if VA > 0 else np.zeros(0); vB = VB * uB if VB > 0 else np.zeros(0)
    v = join(vA, vB) if len(vA) else vB
    t = np.arange(len(v)) * DELTA_T
    tot = sd_signal(v, Ts) * 1e3
    pga = pga_g(v); R = REC[rec]
    flags = []
    if pga > PGA_RATIO_MAX * scale * R["pga"]: flags.append("pga")
    if pga > PGA_WARN_G: flags.append("table")
    if tot[0] / tA > LEAK_MAX: flags.append("leakA")
    if tot[1] / tB > LEAK_MAX: flags.append("leakB")
    return dict(t=t, v=v, tA=tA, tB=tB, stages=stages, TgA=TgA, VA=VA, pgaA=pga_g(vA), TgB=TgB, VB=VB, pgaB=pga_g(vB),
                pga=pga, pgv=float(np.max(np.abs(v))), SdA_tot=float(tot[0]), SdB_tot=float(tot[1]),
                flags="|".join(flags) if flags else "ok", fit_iter=it_n, fit_resid=resid, table_s=float(t[-1]))

def target_row(n, rec, sc, S):
    R = REC[rec]
    return dict(run=n, record=rec, scale=sc, T1_s=round(T1_S, 4), Teff_s=round(TEFF_S, 4), Tm_s=round(R["Tm"], 4), Tp_s=round(R["Tp"], 4),
                flags=S["flags"], stages=S["stages"], target_A_mm=round(S["tA"], 3), target_B_mm=round(S["tB"], 3),
                T_A_s=round(S["TgA"], 5), cyclesA=CYCLES_A, VampA_mps=round(S["VA"], 5), PGA_A_g=round(S["pgaA"], 4),
                T_B_s=round(S["TgB"], 5), cyclesB=CYCLES_B, phaseB=int(PHASE_B), VampB_mps=round(S["VB"], 5), PGA_B_g=round(S["pgaB"], 4),
                PGA_tot_g=round(S["pga"], 4), PGV_tot_mps=round(S["pgv"], 5),
                SdA_tot_mm=round(S["SdA_tot"], 3), SdB_tot_mm=round(S["SdB_tot"], 3),
                SdA_over_tgtA=round(S["SdA_tot"] / S["tA"], 3), SdB_over_tgtB=round(S["SdB_tot"] / S["tB"], 3),
                PGA_over_rec=round(S["pga"] / (sc * R["pga"]), 3), PGV_over_rec=round(S["pgv"] / (sc * R["pgv"]), 3),
                fit_iter=S["fit_iter"], fit_resid="{:.1e}".format(S["fit_resid"]), table_s=round(S["table_s"], 4))
def run_label(n, rec, sc):
    return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def write_table(path, lbl, S):
    with open(path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(S["t"])))
        for a_, b_ in zip(S["t"], S["v"]):
            f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

# ============================ PREFLIGHT: all 25 targets ==============
print("{:>3} {:>5} {:>5} | {:>7} {:>7} {:>4} | {:>6} {:>7} {:>6} | {:>6} {:>7} {:>6} | {:>6} {:>7} {:>7} {:>7} {:>7} | flags".format(
    "run", "rec", "s", "tgtA", "tgtB", "stg", "T_A", "VA", "pgaA", "T_B", "VB", "pgaB", "pga", "pga/rec", "SdA/tA", "SdB/tB", "pgv/rec"))
SIG = {}; FLAGGED = []
with open(TARGETS_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=TARGET_COLS); w.writeheader()
    for n, rec, sc in PROTOCOL:
        S = build_signal(rec, sc); SIG[n] = S; row = target_row(n, rec, sc, S); w.writerow(row)
        if S["flags"] != "ok":
            FLAGGED.append(n)
        print("{:>3} {:>5} {:>5.2f} | {:>7.2f} {:>7.2f} {:>4} | {:>6.3f} {:>7.4f} {:>6.3f} | {:>6.3f} {:>7.4f} {:>6.3f} | {:>6.3f} {:>7.2f} {:>7.2f} {:>7.2f} {:>7.2f} | {}".format(
            n, rec, sc, S["tA"], S["tB"], S["stages"], S["TgA"], S["VA"], S["pgaA"], S["TgB"], S["VB"], S["pgaB"], S["pga"],
            row["PGA_over_rec"], row["SdA_over_tgtA"], row["SdB_over_tgtB"], row["PGV_over_rec"], S["flags"]))
        write_table(os.path.join(OUT_DIR, run_label(n, rec, sc) + "_vel.txt"), run_label(n, rec, sc), S)
print("-> {}  ({} runs, {} flagged: {})".format(TARGETS_CSV, len(SIG), len(FLAGGED), FLAGGED or "none"))
if PREFLIGHT_ONLY:
    print("PREFLIGHT_ONLY: stopping before 3DEC.")
    raise SystemExit(0)
if GATE == "stop" and FLAGGED:
    raise RuntimeError("GATE='stop': runs {} failed the feasibility gate -- see {}".format(FLAGGED, TARGETS_CSV))
if it is None:
    raise RuntimeError("itasca module not available: run inside 3DEC, or set PREFLIGHT_ONLY = True")

# ============================ MODEL SETUP ============================
def call_instruments():
    for f in INSTRUMENT_FILES:
        if not os.path.isfile(f):
            raise RuntimeError("missing instrumentation file: " + f)
        it.command("call '{}'".format(f))
    # joist slip probes of instrument_tilt_v2_zup.dat at this build's joist centres (Z)
    it.command("[jy_lower = {:.5f}]".format(JOIST_CENTRE_Z[0]))
    it.command("[jy_upper = {:.5f}]".format(JOIST_CENTRE_Z[1]))

CHK_GROUP, CHK_SLOT = "MasJchk", "chk"
FISH_OK = True
FISH_SRC = """
fish define chk_force_tb
    ; a FISH function name used on the right-hand side CALLS the function (recursion) -- accumulate in a local
    local fsum = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'T_B') then
            fsum = fsum + comp.z(block.force.app(b))   ; vertical is z
        endif
    endloop
    chk_force_tb = fsum
end
fish define chk_tb_z
    ; frame guard: the top beam's highest centroid (Z-up: ~2.63 m; old frame: Z ~ 0.1 m)
    local zmax = -1e30
    loop foreach local b block.list
        if block.isgroup(b, 'T_B') then
            zmax = math.max(zmax, comp.z(block.pos(b)))
        endif
    endloop
    chk_tb_z = zmax
end
fish define chk_props
    ; masonry-masonry subcontacts (group MasJchk in slot chk): min / max of the joint properties and the number of
    ; subcontacts that differ from the build's FISH values (tension, cohesion, kn, friction defined by mat_case.dat and
    ; stored in the save). Properties live on subcontacts in 3DEC.
    local n = 0
    local nbad = 0
    local nzero = 0
    global chk_tmin = 1e30
    global chk_tmax = -1e30
    global chk_cmin = 1e30
    global chk_cmax = -1e30
    global chk_kmin = 1e30
    global chk_kmax = -1e30
    global chk_fmin = 1e30
    global chk_fmax = -1e30
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'MasJchk', 'chk') then
            loop foreach local sc block.contact.subcontactlist(cx)
                local tn = block.subcontact.prop(sc, 'tension')
                local co = block.subcontact.prop(sc, 'cohesion')
                local kk = block.subcontact.prop(sc, 'stiffness-normal')
                local fr = block.subcontact.prop(sc, 'friction')
                n = n + 1
                chk_tmin = math.min(chk_tmin, tn)
                chk_tmax = math.max(chk_tmax, tn)
                chk_cmin = math.min(chk_cmin, co)
                chk_cmax = math.max(chk_cmax, co)
                chk_kmin = math.min(chk_kmin, kk)
                chk_kmax = math.max(chk_kmax, kk)
                chk_fmin = math.min(chk_fmin, fr)
                chk_fmax = math.max(chk_fmax, fr)
                if tn = 0.0 then
                    nzero = nzero + 1
                endif
                local off = 0
                if math.abs(tn - tension) > chk_tol * math.max(math.abs(tension), 1.0) then
                    off = 1
                endif
                if math.abs(co - cohesion) > chk_tol * math.max(math.abs(cohesion), 1.0) then
                    off = 1
                endif
                if math.abs(kk - kn) > chk_tol * kn then
                    off = 1
                endif
                if math.abs(fr - friction) > chk_tol * math.max(math.abs(friction), 1.0) then
                    off = 1
                endif
                nbad = nbad + off
            endloop
        endif
    endloop
    global chk_nzero = nzero
    global chk_nbad = nbad
    chk_props = n
end
"""
def define_chk_fish():
    """model restore replaces the FISH environment with the save's own, so the read-back functions are
    (re)defined after every restore, never once at import. An intrinsic this 3DEC build does not know
    disables the check instead of killing the run."""
    global FISH_OK
    FISH_OK = True
    try:
        it.command(FISH_SRC)
    except Exception as e:
        FISH_OK = False
        print("build read-back FISH not accepted by this 3DEC ({}); the check is skipped".format(e))
def fish_value(name):
    """evaluate a FISH function and fetch the result through a plain symbol (portable across itasca python builds)."""
    it.command("[{0}_val = {0}]".format(name))
    return float(it.fish.get(name + "_val"))

def audit_props(state_source, strict):
    """compare masonry-masonry subcontact properties with the build's FISH values (tension, cohesion, kn, friction from
    mat_case.dat, stored in the save). strict: stop on a mismatch (after setup); otherwise report only (after a run --
    a damage model may store degraded strengths in the properties)."""
    try:
        b = {k: float(it.fish.get(k)) for k in ("tension", "cohesion", "kn", "friction")}
    except Exception as e:
        print("verify [{}]: build FISH values not in the save ({}); properties reported, not compared".format(state_source, e))
        b = None
    if b is None:
        return
    it.command("[global chk_tol = {:g}]".format(PROP_TOL))
    n = int(fish_value("chk_props"))
    if n == 0:
        print("  ! no masonry-masonry contacts in group '{}' -- properties not verified".format(CHK_GROUP)); return
    g = lambda k: float(it.fish.get(k))
    nbad, nzero = int(g("chk_nbad")), int(g("chk_nzero"))
    print("verify [{}]: {} masonry-masonry subcontacts; build: tension {:.4g} Pa, cohesion {:.4g} Pa, kn {:.4g} Pa/m, "
          "friction {:.4g} deg".format(state_source, n, b["tension"], b["cohesion"], b["kn"], b["friction"]))
    print("    found: tension {:.4g}-{:.4g}, cohesion {:.4g}-{:.4g}, kn {:.4g}-{:.4g}, friction {:.4g}-{:.4g};  "
          "{} off the build values ({:.2%}), {} with zero tension".format(
              g("chk_tmin"), g("chk_tmax"), g("chk_cmin"), g("chk_cmax"), g("chk_kmin"), g("chk_kmax"),
              g("chk_fmin"), g("chk_fmax"), nbad, nbad / n, nzero))
    if strict and nbad > PROP_MISMATCH_MAX * n:
        raise RuntimeError("{:.1%} of the masonry-masonry subcontacts differ from the build values in the save's FISH "
                           "(tension/cohesion/kn/friction): wrong save, or properties reassigned after the build".format(nbad / n))

def setup_dynamic(save_file, state_source):
    it.command("model restore '{}'".format(cmd_path(save_file).replace(".sav", "")))
    it.command("python-reset-state false")
    it.command("model large-strain {}".format("on" if LARGE_STRAIN else "off"))
    it.command("model dynamic active on")
    it.command("block mech damp local 0.0")
    it.command("block mech damp global 0.0")
    if MAXWELL_CMD.strip():
        it.command(MAXWELL_CMD)
    it.command("""
block contact group 'Joist_S1_contact' range pos-z {:.4f} {:.4f} pos-y 0.9 1.5
block contact group 'Joist_S2_contact' range pos-z {:.4f} {:.4f} pos-y 0.9 1.5
""".format(JOIST_BOTTOM_Z[0] - 0.03, JOIST_BOTTOM_Z[0] + 0.03, JOIST_BOTTOM_Z[1] - 0.03, JOIST_BOTTOM_Z[1] + 0.03))
    for g in DRIVE_GROUPS:                       # a control save may still carry its last apply
        try:
            it.command("block gridpoint apply-remove velocity-y range group '{}'".format(g))
        except Exception:
            pass
    it.command("history delete")
    call_instruments()
    it.command("block free velocity-y range group 'S'")
    if FREE_TB_ROTATIONS:
        it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")
    # ---- build verification: precompression as applied force, joint properties exactly as built
    define_chk_fish()
    if not FISH_OK:
        print("verify [{}]: skipped (read-back FISH unavailable)".format(state_source)); return
    it.command("block contact group '{}' slot '{}' range group-intersection 'Masonry' 'Masonry'".format(CHK_GROUP, CHK_SLOT))
    ztb = fish_value("chk_tb_z")
    print("verify [{}]: frame check, top-beam centroid Z = {:.3f} m".format(state_source, ztb))
    if ztb < 2.5:
        raise RuntimeError("top beam at Z = {:.3f} m: an OLD-frame save; use route2_Tm_ihsan_full25_3dec.py".format(ztb))
    fy = fish_value("chk_force_tb")
    print("verify [{}]: applied force on T_B z = {:.1f} N (expected {:.1f})".format(state_source, fy, -SPRING_FORCE_N))
    if abs(fy + SPRING_FORCE_N) > 0.02 * SPRING_FORCE_N:
        raise RuntimeError("applied force on T_B is {:.1f} N, expected {:.1f}: wrong save?".format(fy, -SPRING_FORCE_N))
    if VERIFY_BUILD:
        audit_props(state_source, strict=True)
    print("setup: large-strain {}, damping '{}', state from {}".format("on" if LARGE_STRAIN else "off", MAXWELL_CMD or "none", state_source))

def save_path(n):
    return os.path.join(OUT_DIR, "r2tm_run_{:02d}.sav".format(n))
def read_hist(path):
    try:
        d = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None
    if d.ndim < 2 or d.shape[1] < 2:
        return None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return d[:, 1] if len(d) >= 5 else None
def read_hist_t(path):
    try:
        d = np.genfromtxt(path, skip_header=2)
    except Exception:
        return None, None
    if d.ndim < 2 or d.shape[1] < 2:
        return None, None
    d = d[np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])]
    return (d[:, 0], d[:, 1]) if len(d) >= 5 else (None, None)
def find_ch(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist(hits[0]) if hits else None
def find_ch_t(folder, key):
    hits = sorted(glob.glob(os.path.join(folder, "*" + key + "*.csv")))
    return read_hist_t(hits[0]) if hits else (None, None)
def achieved_period(t, x, table_s):
    """dominant period of the response: (a) 2 x the time between the largest excursion and the
    nearest opposite-sign extreme, (b) FFT peak in ACH_FFT_BAND_HZ over [0, table_s + 1 s]."""
    if t is None or len(t) < 16:
        return None, None
    x = x - x[0]
    i0 = int(np.argmax(np.abs(x)))
    # (a) half-cycle: walk forward and backward from the peak to the next extreme of opposite sign
    def next_extreme(direction):
        j = i0; best = None
        while 0 < j + direction < len(x) - 1:
            j += direction
            if np.sign(x[j]) == -np.sign(x[i0]) and abs(x[j]) >= abs(x[j - 1]) and abs(x[j]) >= abs(x[j + 1]):
                best = j; break
        return best
    cands = [abs(t[j] - t[i0]) for j in (next_extreme(+1), next_extreme(-1)) if j is not None]
    T_hc = 2.0 * min(cands) if cands else None
    # (b) FFT over the forced part + 1 s
    m = t <= (t[0] + table_s + 1.0)
    xw = x[m] - np.mean(x[m]); dt = float(np.median(np.diff(t[m])))
    if len(xw) < 16 or dt <= 0:
        return T_hc, None
    xw = xw * np.hanning(len(xw)); n = len(xw) * 16
    fr = np.fft.rfftfreq(n, dt); am = np.abs(np.fft.rfft(xw, n))
    b = (fr >= ACH_FFT_BAND_HZ[0]) & (fr <= ACH_FFT_BAND_HZ[1])
    if not np.any(b):
        return T_hc, None
    fp = fr[b][np.argmax(am[b])]
    return T_hc, (1.0 / fp if fp > 0 else None)
def edp_series(folder):
    x = find_ch(folder, KEY_DISP_EXP)
    if x is not None:
        return x, KEY_DISP_EXP
    c3, c4, c5 = find_ch(folder, CH3), find_ch(folder, CH4), find_ch(folder, CH5)
    if c3 is None or c4 is None or c5 is None:
        return None, "n/a"
    n = min(len(c3), len(c4), len(c5))
    return (0.5 * (c3[:n] + c4[:n]) - c5[:n]) * 1000.0, "0.5(Ch3+Ch4)-Ch5"
def peak_resid(folder):
    x, src = edp_series(folder)
    if x is None:
        return None, None, src
    ntail = max(5, len(x) // 20)
    return float(np.max(np.abs(x - x[0]))), float(np.mean(x[-ntail:]) - x[0]), src
def peak_bs(folder):
    x = find_ch(folder, KEY_BS)
    return float(np.max(np.abs(x))) if x is not None else None

def load_ckpt():
    if os.path.isfile(CKPT):
        with open(CKPT) as f:
            s = json.load(f)
        return s["last_run"], s["summary"]
    return 0, []
def save_ckpt(last, summary):
    with open(CKPT, "w") as f:
        json.dump({"last_run": last, "summary": summary, "resume_after": RESUME_AFTER, "resume_save": RESUME_SAVE,
                   "base_save": BASE_SAVE, "large_strain": LARGE_STRAIN, "damping": MAXWELL_CMD,
                   "periods": PERIOD_NOTE, "T1_s": T1_S, "Teff_s": TEFF_S, "period_A": str(PERIOD_A), "period_B": str(PERIOD_B),
                   "branches": {k: {q: v[q] for q in ("K", "Fy", "Dy", "dmax", "T1", "Teff")} for k, v in BRANCHES.items()},
                   "backbone_combine": BACKBONE_COMBINE, "f_source": F_SOURCE,
                   "cycles": [CYCLES_A, CYCLES_B], "xi": XI, "Tm_s": {k: REC[k]["Tm"] for k in REC}}, f, indent=1)

# ============================ CHECKPOINT =============================
it.command("python-reset-state false")
it.command("program automatic-model-save active off")
last_done, summary = load_ckpt()
if last_done == 0:
    if RESUME_AFTER > 0:
        if not os.path.isfile(RESUME_SAVE):
            raise RuntimeError("RESUME_SAVE not found: " + RESUME_SAVE)
        if os.path.abspath(os.path.dirname(RESUME_SAVE)) == os.path.abspath(OUT_DIR):
            raise RuntimeError("OUT_DIR must be a NEW folder, not the one holding RESUME_SAVE")
        setup_dynamic(RESUME_SAVE, "resume:" + RESUME_SAVE); last_done = RESUME_AFTER
        print("partial run: state after run {} from {}; runs {}..{} go to {}".format(RESUME_AFTER, RESUME_SAVE, RESUME_AFTER + 1, STOP_AFTER, OUT_DIR))
        state_source = RESUME_SAVE
    else:
        if not os.path.isfile(BASE_SAVE):
            raise RuntimeError("missing " + BASE_SAVE)
        setup_dynamic(BASE_SAVE, "base:" + BASE_SAVE); print("fresh start from", BASE_SAVE)
        state_source = BASE_SAVE
else:
    if last_done < RESUME_AFTER:
        raise RuntimeError("checkpoint says last_run {} but RESUME_AFTER is {}: delete {} or fix the config".format(last_done, RESUME_AFTER, CKPT))
    setup_dynamic(save_path(last_done), "checkpoint:" + save_path(last_done)); print("resuming after run", last_done)
    state_source = save_path(last_done)

log_new = (len(summary) == 0) or not os.path.isfile(LOG)
logf = open(LOG, "w" if log_new else "a", newline="")
logw = csv.writer(logf)
if log_new:
    logw.writerow(LOG_COLS)

# ============================ MAIN LOOP ==============================
for run_no, rec, sc in PROTOCOL[last_done:]:
    if run_no > STOP_AFTER:
        print("STOP_AFTER {} reached".format(STOP_AFTER)); break
    if run_no not in SIG:
        raise RuntimeError("run {}: no signal -- see the preflight table".format(run_no))
    S = SIG[run_no]; lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    print("\n" + "=" * 70)
    print("Run {:02d}  {} x {:.2f}   targets Sd(T1 {:.3f}) {:.2f} mm, Sd(T_eff {:.3f}) {:.2f} mm   flags {}".format(
          run_no, rec, sc, T1_S, S["tA"], TEFF_S, S["tB"], S["flags"]))
    print("  wave A: T {:.4f}s x{}  Vamp {:.4f} m/s  PGA {:.2f} g".format(S["TgA"], CYCLES_A, S["VA"], S["pgaA"]))
    print("  wave B: T {:.4f}s x{}  Vamp {:.4f} m/s  PGA {:.2f} g   stages {}".format(S["TgB"], CYCLES_B, S["VB"], S["pgaB"], S["stages"]))
    print("  whole signal: Sd(T1) {:.3f} / Sd(T_eff) {:.3f} mm ; PGA {:.2f} g PGV {:.3f} m/s ; table {:.3f} s + tail {:.1f} s".format(
          S["SdA_tot"], S["SdB_tot"], S["pga"], S["pgv"], S["table_s"], TAIL_SEC))
    if S["pga"] > PGA_WARN_G:
        print("  ** PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(S["pga"], PGA_WARN_G))
    print("=" * 70)

    tbl = "r2tm_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-y 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(S["table_s"]))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-y range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("model save '{}'".format(cmd_path(save_path(run_no))))
    if VERIFY_BUILD and FISH_OK:
        try:
            audit_props("after run {}".format(run_no), strict=False)
        except Exception as e:
            print("  property audit after run {} skipped ({})".format(run_no, e))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(cmd_path(folder)))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call '{}'".format(EXPORT_FILE))
    peak, resid, src = peak_resid(folder); bs = peak_bs(folder)
    t_x, x_x = find_ch_t(folder, KEY_DISP_EXP)
    T_hc, T_ff = achieved_period(t_x, x_x, S["table_s"]) if t_x is not None else (None, None)
    print("  PEAK EDP ({}) = {}   residual {}   peak base shear (sum m a) {}".format(
        src, "{:.2f} mm".format(peak) if peak is not None else "n/a",
        "{:.2f} mm".format(resid) if resid is not None else "n/a",
        "{:.2f} kN".format(bs) if bs is not None else "n/a"))
    if peak is not None:
        print("  achieved / Sd(T_eff) target = {:.2f}".format(peak / S["tB"]))
    print("  T_achieved: half-cycle {}  FFT {}   (T1 {:.3f} s, T_eff {:.3f} s)".format(
        "{:.3f} s".format(T_hc) if T_hc else "n/a", "{:.3f} s".format(T_ff) if T_ff else "n/a", T1_S, TEFF_S))

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    call_instruments()

    row = target_row(run_no, rec, sc, S)
    row.update(peak_edp_mm=round(peak, 3) if peak is not None else "", resid_edp_mm=round(resid, 3) if resid is not None else "",
               peak_bs_exp_kN=round(bs, 3) if bs is not None else "",
               T_ach_halfcycle_s=round(T_hc, 4) if T_hc else "", T_ach_fft_s=round(T_ff, 4) if T_ff else "",
               edp_source=src, state_source=state_source)
    summary.append(row)
    logw.writerow([row[k] for k in LOG_COLS]); logf.flush()
    save_ckpt(run_no, summary)
    state_source = save_path(run_no)

logf.close()
print("\ndone: {} runs in {}.".format(len(summary), OUT_DIR))