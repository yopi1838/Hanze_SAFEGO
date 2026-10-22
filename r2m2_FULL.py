# -*- coding: ascii -*-
"""
route2_m2_fixedpoint_v13_full25_3dec.py -- FULL PROTOCOL (runs 1-25). Route 2 / Method 2 with the CAPACITY-ADAPTIVE
(fixed-point) target, LARGE STRAIN, MAXWELL 1.5 % (1-25 Hz), on the f_t = 0 build.

THE SIGNAL, every run
    wave A : one velocity cycle at the wall's T1 (T1_S, measured ring-down),
             carrier gridded to DELTA_T.
    wave B : Method 2 packets (cycles and phase from the record):
                 HU12  T_p 0.120 s, 1 cycle, inverted phase
                 EC40  T_p 0.255 s, 1 cycle, same phase
                 FR76  T_p 0.260 s, 2 cycles, same phase
             carrier period = max(T_p, T*) by default (WAVE_B_CARRIER). Preflight on the v5
             backbone: a 0.26 s carrier cannot deliver d* ~ 60 mm at T* ~ 0.6 s (FR76 x1.5-2.0)
             below 1.5-1.7 g, 2.5-2.8x the record's PGA; at T* it takes 0.33-0.46 g (0.6x record).
    targets: whole-signal 5 %-damped Sd(T1)  = s * Sd_record(T1)          (wave A target)
             whole-signal 5 %-damped Sd(T*)  = d*                          (wave B target)
             where d* solves the fixed point  d* = s * Sd_record( T_sec(d*) ),
             T_sec(d) = 2 pi sqrt( M_EFF / (F(d)/d) ) on the PUSHOVER backbone.
             No bilinear, no T_eff: the period is read where the wall sits on its own curve.
    The two amplitudes are solved JOINTLY on the combined signal (Newton on log-amplitudes),
    so the T1 leak of the old driver is gone. Degenerate cases (T* ~ T1, or one wave alone
    already exceeding the other target) collapse to a single wave and are logged as such.

TIPPING: a run with no stable fixed point (demand above the 45-deg line everywhere) is a
    collapse prediction. ON_NO_FIXED_POINT = "stop" halts before that run; "dmax" targets
    the demand at D_MAX and flags the row (the wall is then expected to run away).

PREFLIGHT: set PREFLIGHT_ONLY = False and run with plain python (no 3DEC): every run's
    fixed point, targets, amplitudes, PGA/PGV go to OUT_DIR/m2fp_targets.csv and the
    velocity tables are written. Look at that table before you burn 3DEC time.

START POINTS
    full protocol        : RESUME_AFTER = 0,  BASE_SAVE = the ft0 large-strain springforce save
    FR76 only, from the  : RESUME_AFTER = 21, RESUME_SAVE = "<stratF folder>/stratC_run_21.sav"
    control's run 21       (runs 22..25 into a NEW OUT_DIR; the control folder is never written)
    STOP_AFTER limits how far the loop goes (e.g. 23 to look at 22-23 first).

DAMPING / STRAIN: model large-strain on, block mech damp maxwell <MATLAB 1.5 % triplet>,
    local and global damping 0. The base save carries the 10.35 kN applied force on T_B
    (block.force.app) and tension = 0 on the masonry joints; both are verified after restore.

OUTPUT (same layout as the earlier drivers)
    OUT_DIR/RunNN_REC_sXpYY/            channel CSVs (instrument_history_export_v2.dat)
    OUT_DIR/m2fp_run_NN.sav
    OUT_DIR/RunNN_REC_sXpYY_vel.txt     the applied table
    OUT_DIR/m2fp_targets.csv            preflight table (all 25 runs)
    OUT_DIR/m2fp_log.csv                per-run results
    OUT_DIR/m2fp_checkpoint.json        resume state (delete to start the folder over)
REQUIRES next to this script: PUSHOVER_CSV, spectrum_HU12/EC40/FR76.csv, the save(s),
    instrument_history_new.dat, instrument_baseshear_exp.dat, instrument_history_export_v2.dat
"""
import os, csv, math, json, glob
import numpy as np
try:
    import itasca as it
except ImportError:
    it = None

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False            # True: build targets + tables only, no 3DEC

# ---- start point: FULL protocol, runs 1..25 from the f_t = 0 base save ---------------
# This file is the FULL-PROTOCOL variant. The from-save-21 variant is
# route2_m2_fixedpoint_v13_from21_3dec.py (identical apart from this block).
OUT_DIR      = "m2fp_full25_v13_LS_maxwell"
RESUME_AFTER = 0                                          # 0 = fresh start from BASE_SAVE
RESUME_SAVE  = ""
STOP_AFTER   = 25

BASE_SAVE    = "Part_I_MASON_v8_LargeStrain_springforce.sav"

# ---- wall / backbone ------------------------------------------------------------------
PUSHOVER_CSV     = "pushover_results_v5_ft0_masonrypush/pushover_pos.csv"   # v5 masonry-only push, '#' lines skipped
PUSHOVER_CSV_NEG = ""             # optional: negative branch; if given, |F| of both branches is averaged
D_COL, F_COL     = "d_ctrl_mm", "F_tot_kN"
T1_S             = 0.103          # measured undamped ring-down (9.70-9.83 Hz); Maxwell stiffens ~1.8 % -> 0.102
M_EFF            = "match_T1"     # kg, or "match_T1": M_EFF = K_el*(T1_S/2pi)^2 so the backbone's elastic
                                  #   secant period equals the measured T1 (1635 reproduces the figure numbers)
D_MAX            = 120.0          # mm, search range for the fixed point
ON_NO_FIXED_POINT = "stop"        # "stop" | "dmax"

# ---- signal ---------------------------------------------------------------------------
XI        = 0.05                  # record spectra are 5 %-damped (verified on vel_HU.txt)
DELTA_T   = 0.0025
TAIL_SEC  = 3.0
PACKETS   = {"HU12": dict(Tp=0.120, cycles=1, phase=-1.0),
             "EC40": dict(Tp=0.255, cycles=1, phase=+1.0),
             "FR76": dict(Tp=0.260, cycles=2, phase=+1.0)}
WAVE_B_CARRIER = "longer"         # "longer": carrier = max(record T_p, T*), packet cycles/phase kept  [default, see preflight]
                                  # "record_Tp": Method 2 packets as run in v11/v12 (needs 1.5-1.7 g on FR76 x1.5-2.0)
                                  # "T_star": one cycle at T* (original Route 2 form)
RECORD_VEL = {"HU12": "vel_HU.txt", "FR76": "vel_FR.txt"}   # optional; PGA/PGV of the input are reported relative to s x record
PGA_WARN_G = 0.90                 # table reached 0.78 g
FIT_TOL, FIT_MAXIT = 1e-8, 80

# ---- model ----------------------------------------------------------------------------
LARGE_STRAIN   = True
MAXWELL_CMD    = "block mech damp maxwell 0.0112 1.0005 0.0081 7.4315 0.0096 25.0004"   # 1.5 % flat 1-25 Hz (MATLAB)
DRIVE_GROUPS   = ["S", "T_B"]
SPRING_FORCE_N = 10350.0          # applied force on T_B, -y
VERIFY_TENSION = True             # read back max 'tension' on the masonry joints
TENSION_MAX_PA = 1.0
INSTRUMENT_FILES = ["instrument_history_new.dat", "instrument_baseshear_exp.dat"]
EXPORT_FILE      = "instrument_history_export_v2.dat"
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
KEY_BS = "base_shear_exp_kN"
# ======================================================================

os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "m2fp_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "m2fp_log.csv")
TARGETS_CSV = os.path.join(OUT_DIR, "m2fp_targets.csv")
TARGET_COLS = ["run", "record", "scale", "d_star_mm", "T_star_s", "tipping_mm", "fp_flag",
               "target_A_mm", "target_B_mm", "stages", "T_A_s", "VampA_mps", "PGA_A_g",
               "T_B_s", "cyclesB", "phaseB", "VampB_mps", "PGA_B_g", "PGA_tot_g", "PGV_tot_mps",
               "SdA_tot_mm", "SdB_tot_mm", "SdA_over_tgtA", "PGA_over_rec", "PGV_over_rec", "fit_iter", "fit_resid", "table_s"]
LOG_COLS = TARGET_COLS + ["peak_edp_mm", "resid_edp_mm", "peak_bs_exp_kN", "T_ach_halfcycle_s", "T_ach_fft_s",
                          "T_ach_over_Tstar", "edp_source", "state_source"]
ACH_FFT_BAND_HZ = (1.0, 20.0)

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

# ============================ BACKBONE + SPECTRA ======================
def read_backbone(path):
    rows = [r for r in open(path) if r.strip() and not r.lstrip().startswith("#")]
    rd = csv.DictReader(rows)
    d, F = [], []
    for r in rd:
        try:
            dd, ff = float(r[D_COL]), abs(float(r[F_COL]))
        except (KeyError, ValueError):
            continue
        if dd > 0 and (not d or dd > d[-1]):          # monotone in d: drop unloading / repeated rows
            d.append(dd); F.append(ff)
    if len(d) < 5:
        raise RuntimeError("backbone {}: fewer than 5 usable rows ({} / {} columns?)".format(path, D_COL, F_COL))
    return np.array(d), np.array(F)
if not os.path.isfile(PUSHOVER_CSV):
    raise RuntimeError("missing PUSHOVER_CSV: " + PUSHOVER_CSV)
d_bb, F_bb = read_backbone(PUSHOVER_CSV)
BACKBONE_NOTE = "pos"
if PUSHOVER_CSV_NEG:
    dn, Fn = read_backbone(PUSHOVER_CSV_NEG)
    dmax_c = min(d_bb[-1], dn[-1]); dg = np.linspace(d_bb[0], dmax_c, 600)
    F_bb = 0.5 * (np.interp(dg, d_bb, F_bb) + np.interp(dg, dn, Fn)); d_bb = dg; BACKBONE_NOTE = "pos+neg mean"
K_EL = F_bb[0] / d_bb[0]                                          # kN/mm, first point
if M_EFF == "match_T1":
    M_EFF_KG = K_EL * 1e6 * (T1_S / (2 * math.pi)) ** 2
    M_NOTE = "match_T1 -> {:.0f} kg".format(M_EFF_KG)
else:
    M_EFF_KG = float(M_EFF); M_NOTE = "{:.0f} kg (given)".format(M_EFF_KG)
def F_of(d):
    d = np.asarray(d, float)
    out = np.interp(d, d_bb, F_bb)
    out = np.where(d < d_bb[0], K_EL * d, out)
    return np.where(d > d_bb[-1], np.maximum(F_bb[-1] * (1.0 - (d - d_bb[-1]) / 300.0), 1e-6), out)   # mild extrapolation
def T_sec(d):
    d = np.asarray(d, float)
    return 2.0 * np.pi * np.sqrt(M_EFF_KG / (F_of(d) * 1e3 / (d * 1e-3)))
SPEC = {}
for k, f in SPECTRA.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing spectrum file: {}".format(f))
    s = np.genfromtxt(f, delimiter=",", skip_header=1); SPEC[k] = (s[:, 0], s[:, 1] * 1e3)   # mm
def sd_record_mm(rec, T):
    Ts, Sd = SPEC[rec]; return np.interp(T, Ts, Sd)
REC_PG = {}
for k, f in RECORD_VEL.items():
    if os.path.isfile(f):
        d = np.loadtxt(f, skiprows=2); dtr = d[1, 0] - d[0, 0]
        REC_PG[k] = (float(np.max(np.abs(np.gradient(d[:, 1], dtr))) / 9.81), float(np.max(np.abs(d[:, 1]))))
_DGRID = np.linspace(0.02, D_MAX, 6000)
def fixed_point(rec, scale):
    """stable d* (demand crosses the 45-deg line from above), first tipping point after it (from below)."""
    g = scale * sd_record_mm(rec, T_sec(_DGRID)) - _DGRID
    i = np.where((g[:-1] > 0) & (g[1:] <= 0))[0]
    j = np.where((g[:-1] <= 0) & (g[1:] > 0))[0]
    ds = float(np.interp(0.0, [g[i[0] + 1], g[i[0]]], [_DGRID[i[0] + 1], _DGRID[i[0]]])) if len(i) else None
    jt = [jj for jj in j if ds is None or _DGRID[jj] > ds]
    dt = float(_DGRID[jt[0]]) if jt else None
    return ds, dt
print("backbone {} ({}): K_el {:.3f} kN/mm at {:.3f} mm, {} rows to {:.1f} mm; M_EFF {}; T_sec(first) {:.4f} s vs T1_S {:.4f}".format(
    PUSHOVER_CSV, BACKBONE_NOTE, K_EL, d_bb[0], len(d_bb), d_bb[-1], M_NOTE, float(T_sec(d_bb[0])), T1_S))

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
def sd_signal(v, Ts):
    return sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, Ts, XI)
def pga_g(v):
    return float(np.max(np.abs(np.gradient(v, DELTA_T))) / 9.81) if len(v) > 1 else 0.0

def build_signal(rec, scale):
    pk = PACKETS[rec]
    ds, dt = fixed_point(rec, scale)
    flag = "ok"
    if ds is None:
        if ON_NO_FIXED_POINT == "stop":
            raise RuntimeError("run {} x{:.2f}: no stable fixed point up to {:.0f} mm (collapse prediction) -- ON_NO_FIXED_POINT='stop'".format(rec, scale, D_MAX))
        ds = D_MAX; flag = "no_fixed_point:dmax"
    T_star = float(T_sec(ds))
    tA = float(scale * sd_record_mm(rec, T1_S))                     # mm
    tB = float(scale * sd_record_mm(rec, T_star))                   # = d* when flag ok
    uA, TgA = cycles(1.0, T1_S, 1)
    if WAVE_B_CARRIER == "record_Tp":
        TpB, ncB = pk["Tp"], pk["cycles"]
    elif WAVE_B_CARRIER == "T_star":
        TpB, ncB = T_star, 1
    elif WAVE_B_CARRIER == "longer":
        TpB, ncB = max(pk["Tp"], T_star), pk["cycles"]
    else:
        raise RuntimeError("WAVE_B_CARRIER must be record_Tp | T_star | longer")
    uB, TgB = cycles(pk["phase"], TpB, ncB)
    Ts = [T1_S, T_star]
    # ---- single-wave degeneracies: calibrate each alone and see whether it already covers the other target
    sA = sd_signal(uA, Ts) * 1e3; sB = sd_signal(uB, Ts) * 1e3          # mm per unit amplitude
    VA0 = tA / sA[0]; VB0 = tB / sB[1]
    stages, it_n, resid = "A+B", 0, 0.0
    if VA0 * sA[1] >= tB * (1 - 1e-6):              # wave A alone already gives >= target B at T*  (T* ~ T1)
        VA, VB, stages = VA0, 0.0, "A"
    elif VB0 * sB[0] >= tA * (1 - 1e-6):            # wave B alone already covers target A at T1
        VA, VB, stages = 0.0, VB0, "B"
    else:
        x = np.log([VA0, VB0])                      # Newton on log-amplitudes, FD Jacobian
        def resid_of(x):
            v = join(np.exp(x[0]) * uA, np.exp(x[1]) * uB); s = sd_signal(v, Ts) * 1e3
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
            dx = np.clip(dx, -1.0, 1.0)
            lam = 1.0
            for _ in range(8):                      # backtracking
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
    if not len(vB):
        v = vA
    t = np.arange(len(v)) * DELTA_T
    tot = sd_signal(v, Ts) * 1e3                                        # mm
    return dict(t=t, v=v, d_star=ds, T_star=T_star, tipping=dt, flag=flag, tA=tA, tB=tB, stages=stages,
                TgA=TgA, VA=VA, pgaA=pga_g(vA), TgB=TgB, cyclesB=ncB, phaseB=pk["phase"], VB=VB, pgaB=pga_g(vB),
                pga=pga_g(v), pgv=float(np.max(np.abs(v))), SdA_tot=float(tot[0]), SdB_tot=float(tot[1]),
                fit_iter=it_n, fit_resid=resid, table_s=float(t[-1]))

def target_row(n, rec, sc, S):
    return dict(run=n, record=rec, scale=sc, d_star_mm=round(S["d_star"], 3), T_star_s=round(S["T_star"], 4),
                tipping_mm=round(S["tipping"], 1) if S["tipping"] is not None else "", fp_flag=S["flag"],
                target_A_mm=round(S["tA"], 3), target_B_mm=round(S["tB"], 3), stages=S["stages"],
                T_A_s=round(S["TgA"], 5), VampA_mps=round(S["VA"], 5), PGA_A_g=round(S["pgaA"], 4),
                T_B_s=round(S["TgB"], 5), cyclesB=S["cyclesB"], phaseB=int(S["phaseB"]), VampB_mps=round(S["VB"], 5), PGA_B_g=round(S["pgaB"], 4),
                PGA_tot_g=round(S["pga"], 4), PGV_tot_mps=round(S["pgv"], 5),
                SdA_tot_mm=round(S["SdA_tot"], 3), SdB_tot_mm=round(S["SdB_tot"], 3),
                SdA_over_tgtA=round(S["SdA_tot"] / S["tA"], 3) if S["tA"] > 0 else "",
                PGA_over_rec=round(S["pga"] / (sc * REC_PG[rec][0]), 3) if rec in REC_PG else "",
                PGV_over_rec=round(S["pgv"] / (sc * REC_PG[rec][1]), 3) if rec in REC_PG else "",
                fit_iter=S["fit_iter"], fit_resid="{:.1e}".format(S["fit_resid"]), table_s=round(S["table_s"], 4))
def run_label(n, rec, sc):
    return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def write_table(path, lbl, S):
    with open(path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(S["t"])))
        for a_, b_ in zip(S["t"], S["v"]):
            f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

# ============================ PREFLIGHT: all 25 targets ==============
print("\nFIXED-POINT TARGETS  (T1 {:.4f} s; wave B carrier = record T_p; targets at 5 % damping)".format(T1_S))
print("wave B carrier: {}".format(WAVE_B_CARRIER))
print("{:>3} {:>5} {:>5} | {:>7} {:>6} {:>7} | {:>7} {:>7} {:>4} | {:>7} {:>7} {:>6} {:>6} {:>6} | {:>7} {:>7} {:>7} | {:>3} {:>7}".format(
    "run", "rec", "s", "d*mm", "T*s", "tip mm", "tgtA", "tgtB", "stg", "VA m/s", "VB m/s", "pgaA", "pgaB", "pga", "SdA/tgA", "pga/rec", "pgv/rec", "it", "resid"))
SIG = {}
with open(TARGETS_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=TARGET_COLS); w.writeheader()
    for n, rec, sc in PROTOCOL:
        try:
            S = build_signal(rec, sc)
        except RuntimeError as e:
            print("{:>3} {:>5} {:>5.2f} | {}".format(n, rec, sc, e)); continue
        SIG[n] = S; row = target_row(n, rec, sc, S); w.writerow(row)
        print("{:>3} {:>5} {:>5.2f} | {:>7.2f} {:>6.3f} {:>7} | {:>7.2f} {:>7.2f} {:>4} | {:>7.4f} {:>7.4f} {:>6.3f} {:>6.3f} {:>6.3f} | {:>7} {:>7} {:>7} | {:>3d} {:>7.1e}{}".format(
            n, rec, sc, S["d_star"], S["T_star"], row["tipping_mm"], S["tA"], S["tB"], S["stages"], S["VA"], S["VB"],
            S["pgaA"], S["pgaB"], S["pga"], row["SdA_over_tgtA"], row["PGA_over_rec"], row["PGV_over_rec"], S["fit_iter"], S["fit_resid"],
            "  ** PGA > {:.2f} g".format(PGA_WARN_G) if S["pga"] > PGA_WARN_G else ("  ** " + S["flag"] if S["flag"] != "ok" else "")))
        write_table(os.path.join(OUT_DIR, run_label(n, rec, sc) + "_vel.txt"), run_label(n, rec, sc), S)
print("-> {}  ({} runs)".format(TARGETS_CSV, len(SIG)))
if PREFLIGHT_ONLY:
    print("PREFLIGHT_ONLY: stopping before 3DEC.")
    raise SystemExit(0)
if it is None:
    raise RuntimeError("itasca module not available: run inside 3DEC, or set PREFLIGHT_ONLY = False")

# ============================ MODEL SETUP ============================
def call_instruments():
    for f in INSTRUMENT_FILES:
        if not os.path.isfile(f):
            raise RuntimeError("missing instrumentation file: " + f)
        it.command("call '{}'".format(f))

CHK_GROUP, CHK_SLOT = "MasJchk", "chk"
FISH_OK = True
FISH_SRC = """
fish define chk_force_tb
    ; a FISH function name used on the right-hand side CALLS the function (recursion) -- accumulate in a local
    local fsum = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'T_B') then
            fsum = fsum + comp.y(block.force.app(b))
        endif
    endloop
    chk_force_tb = fsum
end
fish define chk_tension_max
    ; max 'tension' property over the subcontacts of masonry-masonry contacts (group MasJchk in slot chk,
    ; assigned by range group-intersection after restore). Properties live on subcontacts in 3DEC.
    local tmax = -1.0
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'MasJchk', 'chk') then
            loop foreach local sc block.contact.subcontactlist(cx)
                local tn = block.subcontact.prop(sc, 'tension')
                if tn > tmax then
                    tmax = tn
                endif
            endloop
        endif
    endloop
    chk_tension_max = tmax
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
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
    for g in DRIVE_GROUPS:                       # a control save may still carry its last apply
        try:
            it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
        except Exception:
            pass
    it.command("history delete")
    call_instruments()
    it.command("block free velocity-z range group 'S'")
    it.command("""
block free rotation-y range group 'T_B'
block free rotation-z range group 'T_B'
""")
    # ---- build verification: precompression as applied force, joints unbonded
    define_chk_fish()
    if not FISH_OK:
        print("verify [{}]: skipped (read-back FISH unavailable)".format(state_source)); return
    it.command("block contact group '{}' slot '{}' range group-intersection 'Masonry' 'Masonry'".format(CHK_GROUP, CHK_SLOT))
    fy = fish_value("chk_force_tb")
    print("verify [{}]: applied force on T_B y = {:.1f} N (expected {:.1f})".format(state_source, fy, -SPRING_FORCE_N))
    if abs(fy + SPRING_FORCE_N) > 0.02 * SPRING_FORCE_N:
        raise RuntimeError("applied force on T_B is {:.1f} N, expected {:.1f}: wrong save?".format(fy, -SPRING_FORCE_N))
    if VERIFY_TENSION:
        tmax = fish_value("chk_tension_max")
        print("verify [{}]: max joint tension on masonry-masonry subcontacts = {:.3g} Pa".format(state_source, tmax))
        if tmax < 0:
            print("  ! no masonry-masonry contacts in group '{}' -- tension not verified".format(CHK_GROUP))
        elif tmax > TENSION_MAX_PA:
            raise RuntimeError("joint tension {:.3g} Pa > {:.3g}: this is not the f_t = 0 build".format(tmax, TENSION_MAX_PA))
    print("setup: large-strain {}, damping '{}', state from {}".format("on" if LARGE_STRAIN else "off", MAXWELL_CMD or "none", state_source))

def save_path(n):
    return os.path.join(OUT_DIR, "m2fp_run_{:02d}.sav".format(n))
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
                   "backbone": PUSHOVER_CSV, "m_eff_kg": M_EFF_KG, "T1_s": T1_S}, f, indent=1)

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
        raise RuntimeError("run {}: no signal (fixed point failed) -- see the preflight table".format(run_no))
    S = SIG[run_no]; lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    print("\n" + "=" * 70)
    print("Run {:02d}  {} x {:.2f}   d* {:.2f} mm  T* {:.3f} s{}".format(run_no, rec, sc, S["d_star"], S["T_star"],
          "  (tipping at {:.0f} mm)".format(S["tipping"]) if S["tipping"] is not None else ""))
    print("  wave A: T {:.4f}s  target {:.2f} mm  Vamp {:.4f} m/s  PGA {:.2f} g".format(S["TgA"], S["tA"], S["VA"], S["pgaA"]))
    print("  wave B: T {:.4f}s x{}  target {:.2f} mm  Vamp {:.4f} m/s  PGA {:.2f} g   stages {}".format(S["TgB"], S["cyclesB"], S["tB"], S["VB"], S["pgaB"], S["stages"]))
    print("  whole signal: Sd(T1) {:.3f} / Sd(T*) {:.3f} mm ; PGA {:.2f} g PGV {:.3f} m/s ; table {:.3f} s + tail {:.1f} s".format(
          S["SdA_tot"], S["SdB_tot"], S["pga"], S["pgv"], S["table_s"], TAIL_SEC))
    if S["pga"] > PGA_WARN_G:
        print("  ** PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(S["pga"], PGA_WARN_G))
    print("=" * 70)

    tbl = "m2fp_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(S["table_s"]))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("model save '{}'".format(cmd_path(save_path(run_no))))
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
    if peak is not None and S["flag"] == "ok":
        print("  achieved / d* = {:.2f}".format(peak / S["d_star"]))
    print("  T_achieved: half-cycle {}  FFT {}   vs T* {:.3f} s{}".format(
        "{:.3f} s".format(T_hc) if T_hc else "n/a", "{:.3f} s".format(T_ff) if T_ff else "n/a", S["T_star"],
        "   (ratio {:.2f})".format(T_hc / S["T_star"]) if T_hc else ""))

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    call_instruments()

    row = target_row(run_no, rec, sc, S)
    row.update(peak_edp_mm=round(peak, 3) if peak is not None else "", resid_edp_mm=round(resid, 3) if resid is not None else "",
               peak_bs_exp_kN=round(bs, 3) if bs is not None else "",
               T_ach_halfcycle_s=round(T_hc, 4) if T_hc else "", T_ach_fft_s=round(T_ff, 4) if T_ff else "",
               T_ach_over_Tstar=round(T_hc / S["T_star"], 3) if T_hc else "", edp_source=src, state_source=state_source)
    summary.append(row)
    logw.writerow([row[k] for k in LOG_COLS]); logf.flush()
    save_ckpt(run_no, summary)
    state_source = save_path(run_no)

logf.close()
print("\ndone: {} runs in {}.".format(len(summary), OUT_DIR))