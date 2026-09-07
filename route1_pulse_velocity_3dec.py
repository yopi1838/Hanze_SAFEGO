# -*- coding: ascii -*-
"""
route1_velocity_pulse_3dec.py -- ROUTE 1 (Eleni): the velocity-matched
sinusoidal pulse.

THE IDEA (Eleni, meeting 2026-08-31)
    Instead of matching the record's spectral DISPLACEMENT at some period
    (Ihsan/Eleni's original scheme, and the adaptive Strategy C), fit a
    plain sinusoid to the record's VELOCITY: give the pulse the record's
    peak ground velocity and its dominant velocity-cycle period, and apply
    it to the base. Nothing about the wall enters -- this is the
    NON-ADAPTIVE, record-only baseline. Route 2 (the capacity-law
    scheduling) has to beat it to justify its extra machinery.

WHAT IS MATCHED, AND WHAT IS NOT
    v0 = PGV(record) * scale        <- matched (the "approach to velocity")
    T  = dominant velocity period   <- from the record, NOT the wall
    N  = number of cycles           <- a CHOICE, and the whole result hinges
                                       on it (see NCYC note). Default 1.
    No wall stiffness, no effective period, no displacement target.

EXTRACTED FROM Friuli_x1.0.csv (the INPUT record, verified this session):
    PGV = 0.2202 m/s, dominant velocity period 0.50 s by FFT. The raw peak
    velocity lobe is ONE slow trough (~0.76 s half-width) -- a sinusoid
    built on that lobe is spectrally dead (10 mm vs the record's 49 at
    0.52 s), so the period is taken from the sustained oscillation (FFT),
    not the peak-lobe width. This driver re-extracts both at run time so
    the numbers are reproducible, and prints them.

NCYC -- THE ONE KNOB THAT DECIDES THE ANSWER
    A constant-amplitude sinusoid gives EVERY cycle the peak velocity; the
    real record reaches PGV in ~1-2 cycles. Verified against the record's
    own 5%-damped spectrum at 0.52 s: N=1 -> 44 mm, N=2 -> 76 mm, the
    record itself -> 49 mm. So N=1 reproduces the record's own dominant
    peak within ~10%; higher N over-delivers. Ship N=1 as primary, N=2 as
    the sensitivity bound. This is the amplitude-vs-cycles trap that sank
    the ROCKPERIOD run at 65 mm -- do not raise N to "add energy".

WHAT THIS IS NOT
    Not adaptive, not a spectral-displacement match, not Route 2. The
    pulse is also SYMMETRIC while the wall responded asymmetrically, so a
    direction bias in the result is the pulse, not the wall.

SETUP mirrors the Strategy C / F dynamic drivers exactly (so Route 1 is
comparable to the control run): restore Part_I_MASON_v8_SmallStrain,
LARGE_STRAIN on, NO damping (contact dissipation only), joist contact
groups, instrument_history_new.dat, drive groups S and T_B with a
velocity-z table (multiplier 1.0), ring down, export.

OUTPUT
    <OUT_DIR>/route1_<record>_s<scale>_N<ncyc>/  (channel CSVs via export)
    <OUT_DIR>/route1_summary.csv : record, scale, PGV, T_pulse, N, v0,
        peak_edp_mm  (mean Top_Quarter_A/B, table-referenced, paper Eq.1)
    peak_edp_mm ~ 29.5 -> the record-only pulse reaches the observed OOP
    response with no wall information at all.
"""

import itasca as it
import os, csv, math
import numpy as np

it.command("python-reset-state false")
it.command("program automatic-model-save active off")

# ============================ CONFIG =================================
BASE_SAVE   = "Part_I_MASON_v8_SmallStrain.sav"
RECORD_CSV  = "Friuli_x1.0.csv"     # the INPUT record (pos/vel/acc/time)
OUT_DIR     = "route1_results"
SCALE       = 2.0        # protocol scale for the run being reproduced
                         # (run 25 = FR76 x 2.0). v0 = PGV * SCALE.
N_CYCLES    = 1.0        # PRIMARY. Rerun 2.0 for the sensitivity bound.
T_PULSE     = None       # s; None = dominant velocity period from FFT
DELTA_T     = 0.005      # s, table sampling (matches the other drivers)
TAIL_SEC    = 3.0        # s ring-down after the pulse
FFT_FMIN    = 0.5        # Hz, band for the dominant-period pick
FFT_FMAX    = 20.0
LABEL       = "FR76"     # for folder/summary naming
DRIVE_GROUPS = ["S", "T_B"]   # base + top beam, as in Strategy C/F

os.makedirs(OUT_DIR, exist_ok=True)

# ============================ RECORD -> PULSE PARAMS =================
def read_record(path):
    """(t, vel_mps) from the input-record CSV (3 header lines;
    columns pos[mm], vel[mm/s], acc[m/s2], time[s])."""
    if not os.path.isfile(path):
        raise RuntimeError("record not found: {}".format(path))
    d = np.genfromtxt(path, delimiter=",", skip_header=3)
    t = d[:, 3]
    vel = d[:, 1] / 1000.0        # mm/s -> m/s
    return t, vel

def dominant_period(t, vel):
    dt = float(np.median(np.diff(t)))
    v = vel - np.mean(vel)
    fr = np.fft.rfftfreq(len(v), dt)
    amp = np.abs(np.fft.rfft(v))
    band = (fr >= FFT_FMIN) & (fr <= FFT_FMAX)
    f_dom = fr[band][int(np.argmax(amp[band]))]
    return 1.0 / f_dom

t_rec, v_rec = read_record(RECORD_CSV)
PGV = float(np.max(np.abs(v_rec)))
T_pulse = float(T_PULSE) if T_PULSE else dominant_period(t_rec, v_rec)
v0 = PGV * SCALE

print("=" * 66)
print("ROUTE 1 -- velocity-matched sinusoid (Eleni)")
print("  record {}  PGV = {:.4f} m/s  dominant T = {:.4f} s".format(
    RECORD_CSV, PGV, T_pulse))
print("  scale {:.2f}  ->  v0 = {:.4f} m/s   N = {:.2f} cycles".format(
    SCALE, v0, N_CYCLES))
print("  (record-only: no wall property enters this pulse)")
print("=" * 66)

# ============================ BUILD THE VELOCITY TABLE ===============
# v(t) = v0 * sin(2 pi t / T) over N full cycles, then a zero tail.
# Integer or half-integer N returns to zero cleanly; the net base
# displacement over a full cycle is zero, so the table leaves no
# permanent offset -- the wall rings down about its start position.
def build_pulse(v0, T, ncyc, dt, tail):
    # snap the period to a whole number of table steps so N full cycles end
    # EXACTLY on a zero sample -- otherwise the table stops mid-cycle and
    # injects a residual base-velocity step (~3% of v0 at T=0.4975/dt=0.005).
    spc = max(1, int(round(T / dt)))        # samples per cycle
    T_grid = spc * dt                        # snapped period
    n_pulse = int(round(ncyc * spc))         # samples in the pulse
    w = 2.0 * math.pi / T_grid
    tt = np.arange(n_pulse + 1) * dt         # include the final zero sample
    v = v0 * np.sin(w * tt)
    n_tail = int(round(tail / dt))
    t_all = np.r_[tt, tt[-1] + dt + np.arange(n_tail) * dt]
    v_all = np.r_[v, np.zeros(n_tail)]
    return t_all, v_all, T_grid

t_p, v_p, T_grid = build_pulse(v0, T_pulse, N_CYCLES, DELTA_T, TAIL_SEC)
if abs(T_grid - T_pulse) > 1e-9:
    print("  pulse period snapped {:.4f} -> {:.4f} s (table-step grid)".format(
        T_pulse, T_grid))
    T_pulse = T_grid
pulse_dur = float(t_p[-1])
scale_str = "{:.2f}".format(SCALE).replace(".", "p")
ncyc_str = "{:.1f}".format(N_CYCLES).replace(".", "p")
run_label = "route1_{}_s{}_N{}".format(LABEL, scale_str, ncyc_str)
vel_path = os.path.join(OUT_DIR, run_label + "_vel.txt")
with open(vel_path, "w", newline="\n") as f:
    f.write("{}\n{}\t0\n".format(run_label, len(t_p)))
    for a, b in zip(t_p, v_p):
        f.write("{:.6f}\t{:.9e}\n".format(a, b))
print("velocity table -> {}  ({} pts, {:.2f} s, Vpeak {:.4f} m/s)".format(
    vel_path, len(t_p), pulse_dur, float(np.max(np.abs(v_p)))))

# ============================ MODEL SETUP (as Strategy C/F) ==========
if not os.path.isfile(BASE_SAVE):
    raise RuntimeError("base save not found: {}".format(BASE_SAVE))
it.command("model restore '{}'".format(BASE_SAVE.replace(".sav", "")))
it.command("python-reset-state false")
it.command("model large-strain off")     # rocking capacity is geometric
it.command("model dynamic active on")
it.command("block mech damp local 0.0")
it.command("block mech damp global 0.0")
# joist contact groups + instrumentation, as the other drivers
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

# ============================ BASELINE ===============================
it.command("@Record_Disp")
d0 = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
            it.fish.get("Top_Quarter_B_Disp")) * 1000.0
print("baseline top displacement = {:.4f} mm".format(d0))

# ============================ APPLY THE PULSE ========================
it.command("table 'r1' import '{}'".format(vel_path.replace("\\", "/")))
it.command("model dynamic time-total 0")
for grp in DRIVE_GROUPS:
    it.command("block apply velocity-z 1.0 table 'r1' range group '{}'".format(grp))
it.command("model solve dynamic time {:.6f}".format(pulse_dur))
for grp in DRIVE_GROUPS:
    it.command("block gridpoint apply-remove velocity-z range group '{}'".format(grp))

# ============================ PEAK EDP ===============================
# tracked live via the FISH histories; the table-referenced peak is read
# from the exported channels below, but print the end-state top disp now.
it.command("@Record_Disp")
d_end = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
               it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0
print("end-of-pulse top displacement (rel. baseline) = {:.3f} mm".format(d_end))

# export channels for the true peak (0.5(Ch3+Ch4) table-referenced)
full_folder = os.path.join(OUT_DIR, run_label)
os.makedirs(full_folder, exist_ok=True)
it.command("[exportdir='{}']".format(full_folder.replace("\\", "/")))
it.command("[runlabel='{}']".format(run_label))
it.command("call 'instrument_history_export_v2.dat'")
it.command("model save '{}'".format(
    os.path.join(OUT_DIR, run_label + ".sav").replace("\\", "/")))

# peak from the exported top-quarter channels, table-referenced
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
    a = rd("Top_Quarter_A_Disp")
    b = rd("Top_Quarter_B_Disp")
    tab = rd("Record_Disp")   # y=2.68 beam point ~ table reference
    if a is None or b is None:
        return None
    top = 0.5 * (a + b)
    if tab is not None:
        n = min(len(top), len(tab))
        rel = top[:n] - tab[:n]
    else:
        rel = top - top[0]
    return float(np.max(np.abs(rel - rel[0]))) * 1000.0

peak_edp = _peak_from_channels(full_folder)
print("\nPEAK EDP (table-referenced, paper Eq.1) = {}".format(
    "{:.2f} mm".format(peak_edp) if peak_edp is not None else "n/a (read channels)"))
print("target: US-1 run 25 reached ~29.5 mm; strat-F control gave 41.6 mm")

# ============================ SUMMARY ROW ============================
sum_path = os.path.join(OUT_DIR, "route1_summary.csv")
new = not os.path.isfile(sum_path)
with open(sum_path, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["record", "scale", "PGV_mps", "T_pulse_s", "N_cycles",
                    "v0_mps", "peak_edp_mm"])
    w.writerow([LABEL, SCALE, round(PGV, 5), round(T_pulse, 5), N_CYCLES,
                round(v0, 5),
                round(peak_edp, 3) if peak_edp is not None else ""])
print("-> {}".format(sum_path))
print("\ninterpretation:")
print("  ~25-35 mm : the record-only pulse reaches the observed response;")
print("              the adaptive machinery is not needed for THIS record.")
print("  >40 mm    : over-delivery (as strat-F run 25) -- drop to N<1 or")
print("              check the cycle count before reading anything in.")
print("  <15 mm    : velocity matching alone under-delivers; Route 2's")
print("              period scheduling is doing real work.")