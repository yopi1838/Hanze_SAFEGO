# -*- coding: ascii -*-
"""
route2_full25_3dec.py -- ROUTE 2 over the FULL 25-run protocol, SMALL STRAIN.

Same physics and same pulse construction as route2_sequential_pulse_3dec.py
(the single-run diagnostic), applied run after run with the model's damage
state CARRIED BETWEEN RUNS, exactly as the Strategy F control did with the
real records. Comparable to that control run for run: same base save, small
strain, no damping, same drive groups, same instrumentation.

PER RUN (record, scale from Table 2):
    Stage A -- sinusoid at the elastic period T1, amplitude so its 5%-damped
               Sd at T1 equals scale * Sd(record, T1).
    Stage B -- only once the wall has moved: sinusoid at
               T_eff = 2 pi sqrt(m_eff / K_sec(d_hist)), amplitude so its Sd
               at T_eff equals scale * Sd(record, T_eff).
    d_hist = the largest peak EDP the MODEL has produced over completed runs.
    Nothing measured in the experiment enters the loop.

THE LAW (fixed a priori)  : capacity curve + secant construction + m_eff
THE STATE (adaptive)      : d_hist, from the model's own response
Stage B is gated: fires when d_hist >= ONSET_MM and T_eff/T1 >= TRATIO_MIN.
Before that a run is the plain elastic pulse -- an intact wall is never given
softened-period demand.

READ BEFORE INTERPRETING (unchanged from the single-run file):
  * stage B commands scale * Sd(record, T_eff) -- at run 25 that is ~51 mm
    of LINEAR-SDOF spectral displacement, not 29.5. The wall's response is
    the measurement; the gap to 29.5 is the linear-vs-rocking difference.
  * SMALL STRAIN freezes toe migration and P-Delta softening. A shortfall
    cannot be cleanly blamed on the scheme vs the frozen geometry. It is
    used because the control run used it; rerun large strain to separate.
  * N = 1 cycle per stage, deliberately (Route 1 lesson: constant-amplitude
    cycles over-deliver; more cycles is not "more energy", it is a bias).
  * GAP_S between stages is phase-sensitive. Sweep 0.25 / 0.5 / 1.0 on the
    FR76 runs before trusting one number.

CAPACITY SOURCE: the digitised Fig 13 envelope by default (VALIDATION mode
against US-1). Point CAP_CSV at a model pushover CSV for prediction mode --
but NOT the top-driven displacement-control file (that measured sway
resistance, not capacity). The load-controlled ascending branch is valid up
to its 17 mm peak; beyond that capacity_law flat-extrapolates and says so.

RESUME: OUT_DIR/route2_checkpoint.json holds the last completed run and the
summary; a restart restores that run's save and continues. Delete it to
start from run 1.

REQUIRES next to this script: capacity_law.py, spectrum_HU12/EC40/FR76.csv,
US1_fig13_digitised.csv, Part_I_MASON_v8_SmallStrain.sav,
instrument_history_new.dat, instrument_history_export_v2.dat.

OUTPUT: OUT_DIR/route2_full_log.csv (one row per run, appended as it goes)
        OUT_DIR/route2_run_NN.sav + Run folders with channel CSVs
"""

import itasca as it
import os, csv, math, sys, json
import numpy as np

for _p in (os.getcwd(),
           os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)
try:
    import capacity_law as cap
except Exception as _e:
    raise RuntimeError("capacity_law.py must sit next to this script. {}".format(_e))

# ---- inlined, self-tested Newmark Sd (no sdof.py dependency) ----------
def sd_spectrum(ag, dt, Ts, xi=0.05):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2.0 * np.pi / Ts; c = 2.0 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w)
    a = -ag[0] - c * v - w2 * u; umax = np.abs(u)
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
        raise RuntimeError("Sd integrator self-test failed")
_sd_selftest()

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ CONFIG =================================
BASE_SAVE    = "Part_I_MASON_v8_SmallStrain.sav"
OUT_DIR      = "route2_full25"
T1_INIT      = 0.0948
XI           = 0.05
DELTA_T      = 0.005
N_A          = 1.0
N_B          = 1.0
GAP_S        = 0.5
TAIL_SEC     = 3.0
DRIVE_GROUPS = ["S", "T_B"]
ONSET_MM     = 2.0          # stage B fires once the model has moved this far
TRATIO_MIN   = 1.05         # ... and T_eff has separated from T1
T_EFF_MAX    = 1.0
PGA_WARN_G   = 0.90         # the table reached 0.78 g; above this is extrapolation

CAP_CSV      = "US1_fig13_digitised.csv"
CAP_DIR      = "neg"
M_EFF        = 1635.0
SPECTRA      = {"HU12": "spectrum_HU12.csv", "EC40": "spectrum_EC40.csv",
                "FR76": "spectrum_FR76.csv"}

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

os.makedirs(OUT_DIR, exist_ok=True)
CKPT = os.path.join(OUT_DIR, "route2_checkpoint.json")
LOG  = os.path.join(OUT_DIR, "route2_full_log.csv")

# ============================ LAW + SPECTRA ==========================
law = cap.CapacityLaw(CAP_CSV, m_eff=M_EFF, direction=CAP_DIR,
                      t_eff_max=T_EFF_MAX, label="route2_full25")
print(law.report())
SPEC = {}
for k, f in SPECTRA.items():
    if not os.path.isfile(f):
        raise RuntimeError("missing spectrum file: {}".format(f))
    d = np.genfromtxt(f, delimiter=",", skip_header=1); SPEC[k] = (d[:, 0], d[:, 1])
def sd_record(rec, T):
    Ts, Sd = SPEC[rec]; return float(np.interp(T, Ts, Sd))

# ============================ PULSE BUILDERS =========================
def pulse_vel(Vamp, T, ncyc):
    spc = max(1, int(round(T / DELTA_T))); Tg = spc * DELTA_T
    npl = int(round(ncyc * spc)); w = 2.0 * math.pi / Tg
    tt = np.arange(npl + 1) * DELTA_T
    return tt, Vamp * np.sin(w * tt), Tg
def calibrate_vel(T, Sd_target, ncyc):
    _, v, Tg = pulse_vel(1.0, T, ncyc)
    sd0 = float(sd_spectrum(np.gradient(v, DELTA_T), DELTA_T, [Tg], XI)[0])
    if sd0 <= 0: raise RuntimeError("trial Sd=0 at T={:.4f}".format(Tg))
    return Sd_target / sd0, Tg
def stage(rec, scale, T, ncyc):
    Sd_t = scale * sd_record(rec, T)
    Vamp, Tg = calibrate_vel(T, Sd_t, ncyc)
    t, v, _ = pulse_vel(Vamp, T, ncyc)
    pga = float(np.max(np.abs(np.gradient(v, DELTA_T)))) / 9.81
    return dict(T=Tg, Sd=Sd_t, Vamp=Vamp, pga_g=pga, t=t, v=v)

# ============================ MODEL SETUP ============================
def setup_dynamic(save_file):
    it.command("model restore '{}'".format(save_file.replace(".sav", "").replace("\\", "/")))
    it.command("python-reset-state false")
    it.command("model large-strain off")          # SMALL STRAIN, as the control
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

def run_label(n, rec, sc):
    return "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(sc).replace(".", "p"))
def save_path(n):
    return os.path.join(OUT_DIR, "route2_run_{:02d}.sav".format(n))

def peak_from_channels(folder):
    import glob
    def rd(pat):
        hits = sorted(glob.glob(os.path.join(folder, "*" + pat + "*.csv")))
        if not hits: return None
        try:
            a = np.genfromtxt(hits[0], skip_header=2)
            return a[:, 1] if a.ndim == 2 and a.shape[1] >= 2 else None
        except Exception:
            return None
    a = rd("Top_Quarter_A_Disp"); b = rd("Top_Quarter_B_Disp"); tab = rd("Record_Disp")
    if a is None or b is None: return None
    top = 0.5 * (a + b)
    if tab is not None:
        n = min(len(top), len(tab)); rel = top[:n] - tab[:n]
    else:
        rel = top - top[0]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0

# ============================ CHECKPOINT =============================
def load_ckpt():
    if os.path.isfile(CKPT):
        with open(CKPT) as f: s = json.load(f)
        return s["last_run"], s["summary"]
    return 0, []
def save_ckpt(last, summary):
    with open(CKPT, "w") as f: json.dump({"last_run": last, "summary": summary}, f, indent=1)

last_done, summary = load_ckpt()
if last_done == 0:
    if not os.path.isfile(BASE_SAVE): raise RuntimeError("missing " + BASE_SAVE)
    setup_dynamic(BASE_SAVE); print("fresh start from", BASE_SAVE)
else:
    setup_dynamic(save_path(last_done)); print("resuming after run", last_done)

log_new = (last_done == 0) or not os.path.isfile(LOG)
logf = open(LOG, "w" if log_new else "a", newline="")
logw = csv.writer(logf)
if log_new:
    logw.writerow(["run", "record", "scale", "d_hist_mm", "two_stage",
                   "T_A_s", "SdA_mm", "VampA_mps", "PGA_A_g",
                   "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g", "gap_s",
                   "cap_note", "peak_edp_mm"])

# ============================ MAIN LOOP ==============================
for run_no, rec, sc in PROTOCOL[last_done:]:
    d_hist = max([r["peak_edp_mm"] for r in summary if r.get("peak_edp_mm") is not None] + [0.0])
    T_eff, note = law.t_eff(d_hist) if d_hist > 0 else (T1_INIT, "")
    two = (d_hist >= ONSET_MM) and (T_eff >= T1_INIT * TRATIO_MIN)
    if not two: note = (note + "|" if note else "") + "single_stage"

    A = stage(rec, sc, T1_INIT, N_A)
    if two:
        B = stage(rec, sc, T_eff, N_B)
        gap_n = int(round(GAP_S / DELTA_T)); tail_n = int(round(TAIL_SEC / DELTA_T))
        v_all = np.r_[A["v"], np.zeros(gap_n), B["v"], np.zeros(tail_n)]
        excite_dur = float(A["t"][-1] + GAP_S + B["t"][-1])
    else:
        B = None
        v_all = np.r_[A["v"], np.zeros(int(round(TAIL_SEC / DELTA_T)))]
        excite_dur = float(A["t"][-1])
    t_all = np.arange(len(v_all)) * DELTA_T
    total_dur = float(t_all[-1])

    lbl = run_label(run_no, rec, sc)
    vel_path = os.path.join(OUT_DIR, lbl + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(lbl, len(t_all)))
        for a_, b_ in zip(t_all, v_all): f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

    print("\n" + "=" * 66)
    print("Run {:02d}  {} x {:.2f}   d_hist {:.2f} mm   {}".format(run_no, rec, sc, d_hist, note))
    print("  stage A: T {:.4f}s Sd {:.2f}mm Vamp {:.4f} PGA {:.2f}g".format(
        A["T"], A["Sd"] * 1e3, A["Vamp"], A["pga_g"]))
    if two:
        print("  stage B: T {:.4f}s Sd {:.2f}mm Vamp {:.4f} PGA {:.2f}g  (gap {:.2f}s)".format(
            B["T"], B["Sd"] * 1e3, B["Vamp"], B["pga_g"], GAP_S))
    for tag, st in (("A", A), ("B", B)):
        if st and st["pga_g"] > PGA_WARN_G:
            print("  ** stage {} PGA {:.2f} g > {:.2f} g the table reached -- extrapolation".format(
                tag, st["pga_g"], PGA_WARN_G))
    print("=" * 66)

    tbl = "r2_{:02d}".format(run_no)
    it.command("table '{}' import '{}'".format(tbl, vel_path.replace("\\", "/")))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table '{}' range group '{}'".format(tbl, g))
    it.command("model solve dynamic time {:.6f}".format(total_dur))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))

    it.command("model save '{}'".format(save_path(run_no).replace("\\", "/")))
    folder = os.path.join(OUT_DIR, lbl); os.makedirs(folder, exist_ok=True)
    it.command("[exportdir='{}']".format(folder.replace("\\", "/")))
    it.command("[runlabel='{}']".format(lbl))
    it.command("call 'instrument_history_export_v2.dat'")
    peak = peak_from_channels(folder)
    print("  PEAK EDP = {}".format("{:.2f} mm".format(peak) if peak is not None else "n/a"))

    it.command("table '{}' delete".format(tbl))
    it.command("history delete")
    it.command("call 'instrument_history_new.dat'")

    row = dict(run=run_no, record=rec, scale=sc, d_hist_mm=round(d_hist, 3),
               two_stage=1 if two else 0,
               T_A_s=round(A["T"], 5), SdA_mm=round(A["Sd"] * 1e3, 3),
               VampA_mps=round(A["Vamp"], 5), PGA_A_g=round(A["pga_g"], 4),
               T_B_s=round(B["T"], 5) if two else "", SdB_mm=round(B["Sd"] * 1e3, 3) if two else "",
               VampB_mps=round(B["Vamp"], 5) if two else "", PGA_B_g=round(B["pga_g"], 4) if two else "",
               gap_s=GAP_S if two else "", cap_note=note.replace(",", ";"),
               peak_edp_mm=round(peak, 3) if peak is not None else None)
    summary.append(row)
    logw.writerow([row[k] if row[k] is not None else "" for k in
                   ["run", "record", "scale", "d_hist_mm", "two_stage", "T_A_s", "SdA_mm",
                    "VampA_mps", "PGA_A_g", "T_B_s", "SdB_mm", "VampB_mps", "PGA_B_g",
                    "gap_s", "cap_note", "peak_edp_mm"]])
    logf.flush()
    save_ckpt(run_no, summary)

logf.close()
print("\nRoute 2 full protocol complete. Log -> {}".format(LOG))
print("compare peak_edp_mm run by run against US-1 (run 21: 16.8, run 24: 29.5) "
      "and the Strategy F control (run 21: 10.5, run 25: 41.6).")
