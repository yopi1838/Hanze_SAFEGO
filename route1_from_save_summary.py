# -*- coding: ascii -*-
"""
route1_from_save_3dec.py -- ROUTE 1 (Eleni: velocity-matched sinusoid)
applied to a RESTORED Strategy F state, not to the virgin wall.

Same test as route2_from_save_3dec.py, for Route 1:
    case 1: restore stratC_run_24.sav, pulse at scale 2.0  -> vs control run 25
    case 2: restore stratC_run_23.sav, pulse at scale 1.75 -> vs control run 24

THE RULE (as agreed for Route 1)
    v(t) = v0 sin(2 pi t / T),  v0 = PGV x scale,  T = dominant velocity
    period of the signal (FFT peak in FFT_FMIN..FFT_FMAX), N_CYCLES cycles.
    Nothing about the wall enters.

WHICH VELOCITY SIGNAL -- VEL_SOURCE
    "table"  : vel_FR.txt, the achieved table motion the CONTROL was driven
               with (PGV 0.253 m/s at x1.0). Default here, because this test
               compares the pulse with the record ON THE MODEL and the model
               never saw Friuli_x1.0.csv.
    "record" : Friuli_x1.0.csv, the commanded input record (PGV 0.220 m/s),
               as in the original virgin-wall Route 1 run.
    The two differ by ~15% in PGV; both PGVs are printed whichever is used.

STRAIN MODE: small strain, mirroring the control (the virgin-wall Route 1
run used large strain -- that inconsistency is removed here). Setup is the
control's resume path: restore, large-strain off, dynamic on, joist contact
groups, histories, free S / T_B, zero damping.

PEAK EDP: rel_disp_top_exp_mm (0.5(Ch3+Ch4)-Ch5), else rebuilt from the
Channel_3/4/5 CSVs. Not Record_Disp.

OUTPUT
    <OUT_DIR>/route1_<REC>_s<scale>_N<n>/   channel CSVs
    <OUT_DIR>/route1_<REC>_s<scale>_N<n>.sav
    <OUT_DIR>/route1_from_save_summary.csv
"""

import itasca as it
import os, csv, math, glob
import numpy as np

# ============================ CONFIG =================================
CASES = [
    # (save to restore,                              from_run, SCALE, control run to compare, OUT_DIR)
    ("stratF_full_results_US1/stratC_run_24.sav",   24,       2.0,   25, "route1_from24"),
    ("stratF_full_results_US1/stratC_run_23.sav",   23,       1.75,  24, "route1_from23"),
]
CONTROL_DIR = "stratF_full_results_US1"
LABEL       = "FR76"
VEL_SOURCE  = "table"          # "table" (vel_FR.txt, what the control used) or "record" (Friuli_x1.0.csv)
TABLE_FILE  = "vel_FR.txt"     # 3DEC table: name line, count line, then t v
RECORD_CSV  = "Friuli_x1.0.csv"  # 3 header lines; columns pos[mm], vel[mm/s], acc, time[s]
N_CYCLES    = 1.0
T_PULSE     = None             # s; None = dominant velocity period from FFT
FFT_FMIN, FFT_FMAX = 0.5, 20.0
DELTA_T     = 0.005
TAIL_SEC    = 2.5              # ring-down after the table is removed (control uses 2.5)
DRIVE_GROUPS = ["S", "T_B"]
KEY_DISP_EXP, CH3, CH4, CH5 = "rel_disp_top_exp_mm", "Channel_3_DispTopQLeft", "Channel_4_DispTopQRight", "Channel_5_DispTable"

# ============================ HELPERS ================================
def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")

def read_table(path):
    d = np.loadtxt(path, skiprows=2)
    return d[:, 0], d[:, 1]

def read_record(path):
    d = np.genfromtxt(path, delimiter=",", skip_header=3)
    return d[:, 3], d[:, 1] / 1000.0

def dominant_period(t, vel):
    dt = float(np.median(np.diff(t)))
    v = vel - np.mean(vel)
    fr = np.fft.rfftfreq(len(v), dt); amp = np.abs(np.fft.rfft(v))
    band = (fr >= FFT_FMIN) & (fr <= FFT_FMAX)
    return 1.0 / fr[band][int(np.argmax(amp[band]))]

def build_pulse(v0, T, ncyc, dt):
    spc = max(1, int(round(T / dt))); T_grid = spc * dt
    n_pulse = int(round(ncyc * spc)); w = 2.0 * math.pi / T_grid
    tt = np.arange(n_pulse + 1) * dt
    return tt, v0 * np.sin(w * tt), T_grid

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

def control_ref(run_no):
    hits = [d for d in glob.glob(os.path.join(CONTROL_DIR, "Run{:02d}_*".format(run_no))) if os.path.isdir(d)]
    return peak_edp_mm(hits[0]) if hits else (None, "control folder Run{:02d}_* not found".format(run_no))

# ============================ VELOCITY SIGNAL -> PULSE PARAMS ========
pgv = {}
if os.path.isfile(TABLE_FILE):
    t_tab, v_tab = read_table(TABLE_FILE); pgv["table"] = (float(np.max(np.abs(v_tab))), dominant_period(t_tab, v_tab))
if os.path.isfile(RECORD_CSV):
    t_rec, v_rec = read_record(RECORD_CSV); pgv["record"] = (float(np.max(np.abs(v_rec))), dominant_period(t_rec, v_rec))
if VEL_SOURCE not in pgv:
    raise RuntimeError("VEL_SOURCE '{}' file not found ({})".format(VEL_SOURCE, TABLE_FILE if VEL_SOURCE == "table" else RECORD_CSV))
PGV, T_dom = pgv[VEL_SOURCE]
T_pulse = float(T_PULSE) if T_PULSE else T_dom
print("=" * 70)
print("ROUTE 1 on restored states -- velocity source '{}'".format(VEL_SOURCE))
for k, (p, T) in pgv.items():
    print("  {:6s}: PGV {:.4f} m/s  dominant T {:.4f} s{}".format(k, p, T, "   <- used" if k == VEL_SOURCE else ""))
print("  pulse: T {:.4f} s, N {:g} cycles, v0 = PGV x scale".format(T_pulse, N_CYCLES))
print("=" * 70)

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ RUN THE CASES ==========================
for (save, from_run, SCALE, ctrl_run, OUT_DIR) in CASES:
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.isfile(save):
        print("!! save not found, case skipped: {}".format(save)); continue
    v0 = PGV * SCALE
    t_p, v_p, T_grid = build_pulse(v0, T_pulse, N_CYCLES, DELTA_T)
    pulse_dur = float(t_p[-1])
    pga = float(np.max(np.abs(np.gradient(v_p, DELTA_T))) / 9.81)
    pgd = float(np.max(np.abs(np.cumsum(v_p) * DELTA_T))) * 1000.0
    run_label = "route1_{}_s{}_N{}".format(LABEL, "{:.2f}".format(SCALE).replace(".", "p"),
                                            "{:.1f}".format(N_CYCLES).replace(".", "p"))
    vel_path = os.path.join(OUT_DIR, run_label + "_vel.txt")
    with open(vel_path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(run_label, len(t_p)))
        for a, b in zip(t_p, v_p):
            f.write("{:.6f}\t{:.9e}\n".format(a, b))
    print("\nROUTE 1 from {} (end of run {})  scale {:.2f}  -> compare control run {}".format(save, from_run, SCALE, ctrl_run))
    print("  v0 {:.4f} m/s  T {:.4f} s  N {:g}  PGA {:.3f} g  stroke {:.1f} mm  duration {:.2f} s + tail {:.1f} s".format(
        v0, T_grid, N_CYCLES, pga, pgd, pulse_dur, TAIL_SEC))

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

    # ---- apply the pulse (same commands as the control's execute_run)
    it.command("table 'r1' import '{}'".format(cmd_path(vel_path)))
    it.command("model dynamic time-total 0")
    for grp in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table 'r1' range group '{}'".format(grp))
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
    print("  PEAK EDP (Route 1, {}) = {}".format(src, "{:.2f} mm".format(peak) if peak is not None else "n/a"))
    print("  control run {} ({}) = {}".format(ctrl_run, ref_src, "{:.2f} mm".format(ref) if ref is not None else "n/a"))
    if peak is not None and ref:
        print("  ratio Route 1 / record = {:.2f}".format(peak / ref))

    sum_path = os.path.join(OUT_DIR, "route1_from_save_summary.csv")
    new = not os.path.isfile(sum_path)
    with open(sum_path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["from_save", "from_run", "scale", "vel_source", "PGV_mps", "T_pulse_s", "N_cycles",
                        "v0_mps", "PGA_g", "stroke_mm", "peak_edp_mm", "edp_source", "control_ref_run", "control_ref_mm"])
        w.writerow([save, from_run, SCALE, VEL_SOURCE, round(PGV, 5), round(T_grid, 5), N_CYCLES, round(v0, 5),
                    round(pga, 4), round(pgd, 2), round(peak, 3) if peak is not None else "", src, ctrl_run,
                    round(ref, 3) if ref else ""])
    print("  -> {}".format(sum_path))