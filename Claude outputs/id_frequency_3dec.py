# -*- coding: ascii -*-
"""
id_frequency_3dec.py -- fundamental frequency of the SIN model by impulse +
free vibration, compared with Graziotti et al. (2016) Table 6.

This script changes NOTHING in the model. It restores START_SAVE as built by
00_Master_KN.dat (or 00_Master.dat), applies one short base pulse, lets the
wall ring, exports the channels and identifies the frequency. Overburden,
joint state, joint stiffness and damping are whatever the save carries.

METHOD
    One full-sine velocity pulse, PULSE_SEC long, PULSE_G peak acceleration
    (0.05 g = the experiment's random-test amplitude), applied to the
    foundation and the steel beam in Y as 08_IDT_load_*.dat does. A full sine
    has zero net displacement, so the base ends where it started. Then the
    base is held and the wall rings for RING_SEC. Total ~1 s of dynamic time.
    Frequency from the ring-down of the mid-height acceleration:
      (1) mean period between upward zero crossings   (primary)
      (2) zero-padded FFT peak                          (check)
      (3) FFT of the mid-height relative displacement   (independent check)
    plus the log-decrement damping, so a damped and a no-damp save can be
    told apart from the output.

WHAT IT MEASURES
    At 0.05 g no joint opens (bending stress ~34 kPa against 100 kPa of
    precompression), so tension / cohesion / the pre-cracked MidCrack do not
    enter. The frequency is set by kn, ks, block E and the supports. As built
    (kn 9.14e10) the model reads 20.0 Hz; the measured wall reads 14.27 Hz
    (SIN-01-00) / 18.75 Hz (SIN-03-00). The homogenised estimate for the
    joint scaling is 0.31 / 0.75; 00_Master_KN.dat rebuilds the model with it.

MID-HEIGHT SENSORS (added here; instrument_history.dat has none)
    ACC_MID  velocity-y at z = 0.400 + 1.260 m   (paper Table 3, "1/2 Acc. CS", 1260 mm)
    WP_MID   displacement-y at z = 0.400 + 1.340 m (paper Table 3, "1/2 WP CS",  1340 mm)

Usage (3DEC):  model new ; python-reset-state false ; call 'id_frequency_3dec.py'
Output: idfreq_<TAG>/  pulse table, channel CSVs, a save, idfreq_summary.csv
"""

import itasca as it
import os, csv, math, glob
import numpy as np

# ============================ CONFIG =================================
TAG          = "kn0p31"            # matches case_tag in 00_Master_KN.dat; "asbuilt" for the original model
TAG          = os.environ.get("IDFREQ_TAG", TAG).strip()
START_SAVE   = ("SAV/SinleLeaf_onerun_ready.sav" if TAG == "asbuilt"
                else "SAV/SinleLeaf_onerun_ready_{}.sav".format(TAG))
F_TARGET, REF = 14.27, "SIN-01-00 (0.1 MPa, tested after the 0.3 MPa sequence)"
# F_TARGET, REF = 18.75, "SIN-03-00 (0.3 MPa, virgin)"       # for a 0.3 MPa build

DRIVE_GROUPS = ["foundation", "steel beam"]
DRIVE_AXIS   = "y"                 # this model shakes in Y (gravity is -Z)
PULSE_G      = 0.05                # peak base acceleration of the pulse
PULSE_SEC    = 0.03                # one full sine of velocity; energy up to ~1/PULSE_SEC = 33 Hz
RING_SEC     = 1.0                 # free vibration: 14 cycles at 14 Hz, 20 at 20 Hz
DELTA_T      = 0.0005              # pulse table step (60 rows)
ID_BAND      = (5.0, 35.0)         # where the FFT peak is looked for
# History numbering: `history delete`, then instrument_history.dat (ACC2=1, ACC6=2,
# POT11=3, WP13=4, WP14=5, WP15=6), then ACC_MID=7, WP_MID=8, then dynamic time = 9.
# Exported by NAME against TIME_HIST_ID (export_csv.dat is not used: it only writes
# the labels in hist_labels.txt and would drop the two new channels).
HIST = {"ACC2": 1, "ACC6": 2, "POT11": 3, "WP13": 4, "WP14": 5, "WP15": 6, "ACC_MID": 7, "WP_MID": 8}
TIME_HIST_ID = 9
Z_ACC_MID    = 0.400 + 1.260
Z_WP_MID     = 0.400 + 1.340
X_MID, Y_FACE = 1.05, 0.474
OUT_DIR      = "idfreq_" + TAG
os.makedirs(OUT_DIR, exist_ok=True)

def cmd_path(p):
    return os.path.normpath(p).replace("\\", "/")

# ============================ THE PULSE TABLE ========================
# v(t) = V sin(2 pi t / PULSE_SEC): peak acceleration V * 2 pi / PULSE_SEC = PULSE_G * g
V_AMP = PULSE_G * 9.81 * PULSE_SEC / (2.0 * math.pi)
t = np.arange(int(round(PULSE_SEC / DELTA_T)) + 1) * DELTA_T
v = V_AMP * np.sin(2.0 * math.pi * t / PULSE_SEC); v[0] = 0.0; v[-1] = 0.0
d = np.concatenate(([0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * DELTA_T)))
VEL_PATH = os.path.join(OUT_DIR, "idfreq_pulse_vel.txt")
with open(VEL_PATH, "w", newline="\n") as f:
    f.write("idfreq_pulse\n{}\t0\n".format(len(t)))
    for a_, b_ in zip(t, v):
        f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

print("=" * 72)
print("FREQUENCY IDENTIFICATION   [{}]".format(TAG))
print("  start save   : {}   (restored as is -- nothing in the model is changed)".format(START_SAVE))
print("  excitation   : one full-sine pulse, {:.0f} ms, {:.3f} g peak ({:.2f} mm/s), then {:.1f} s free vibration".format(
    PULSE_SEC * 1000, PULSE_G, V_AMP * 1000, RING_SEC))
print("  base stroke  : {:+.4f} mm peak, {:+.5f} mm net   |   dynamic time {:.2f} s total".format(
    np.max(np.abs(d)) * 1000, d[-1] * 1000, PULSE_SEC + RING_SEC))
print("  target       : {:.2f} Hz (T {:.4f} s) -- {}".format(F_TARGET, 1.0 / F_TARGET, REF))
print("  output       : {}".format(OUT_DIR))
print("=" * 72)

# ============================ MODEL ==================================
it.command("python-reset-state false")
it.command("program automatic-model-save active off")
if not os.path.isfile(START_SAVE):
    raise RuntimeError("start save not found: {} (build it with 00_Master_KN.dat, case_tag '{}')".format(START_SAVE, TAG))
it.command("model restore '{}'".format(cmd_path(START_SAVE).replace(".sav", "")))
it.command("python-reset-state false")

# --- histories -------------------------------------------------------
it.command("history delete")
it.command("call 'instrument_history.dat'")
it.command("block history name 'ACC_MID' velocity-{} position {:g},{:g},{:g}".format(DRIVE_AXIS, X_MID, Y_FACE, Z_ACC_MID))
it.command("block history name 'WP_MID' displacement-{} position {:g},{:g},{:g}".format(DRIVE_AXIS, X_MID, Y_FACE, Z_WP_MID))
it.command("model history dynamic time-total")
it.command("history interval 20")

# --- pulse, then hold the base and ring ------------------------------
it.command("table 'idpulse' import '{}'".format(cmd_path(VEL_PATH)))
for g in DRIVE_GROUPS:
    it.command("block gridpoint apply-remove vel-{} range group '{}'".format(DRIVE_AXIS, g))
    it.command("block gridpoint apply velocity-{} 1 table 'idpulse' range group '{}'".format(DRIVE_AXIS, g))
it.command("model dynamic time-total 0")
it.command("model update-interval 2000")
it.command("model solve time {:.6f}".format(PULSE_SEC))
for g in DRIVE_GROUPS:
    it.command("block gridpoint apply-remove velocity-{} range group '{}'".format(DRIVE_AXIS, g))
    it.command("block gridpoint apply vel-{} 0 range group '{}'".format(DRIVE_AXIS, g))
it.command("model solve time {:.6f}".format(RING_SEC))
for nm in HIST:
    it.command("history export '{}' vs {} file '{}' truncate".format(nm, TIME_HIST_ID, cmd_path(os.path.join(OUT_DIR, nm + ".csv"))))
it.command("model save '{}'".format(cmd_path(os.path.join(OUT_DIR, "idfreq_" + TAG + ".sav"))))

# ============================ IDENTIFY ===============================
def read_hist(name):
    for h in sorted(glob.glob(os.path.join(OUT_DIR, "*" + name + "*.csv"))):
        try:
            d_ = np.genfromtxt(h, skip_header=2)            # 3DEC history export: 2 header lines, whitespace
            if d_.ndim == 2 and d_.shape[1] >= 2:
                m = np.isfinite(d_[:, 0]) & np.isfinite(d_[:, 1])
                if m.sum() >= 64:
                    return d_[m, 0], d_[m, 1]
        except Exception:
            pass
    return None, None

def ring_down(tt, yy, t_from):
    m = tt >= t_from
    return tt[m], yy[m] - np.mean(yy[m])

def f_zero_cross(tt, yy):
    s_ = np.sign(yy); idx = np.where((s_[1:] > 0) & (s_[:-1] <= 0))[0]
    if len(idx) < 3:
        return None, len(idx)
    tc = tt[idx] - yy[idx] * (tt[idx + 1] - tt[idx]) / (yy[idx + 1] - yy[idx])
    return float(1.0 / np.mean(np.diff(tc))), len(idx)

def f_fft(tt, yy, band=ID_BAND):
    dt = float(np.median(np.diff(tt))); n = len(yy); npad = max(8 * n, 1 << 16)
    f = np.fft.rfftfreq(npad, dt); A_ = np.abs(np.fft.rfft(yy * np.hanning(n), npad))
    m = (f >= band[0]) & (f <= band[1])
    return float(f[m][int(np.argmax(A_[m]))])

def log_decrement(tt, yy, f0):
    pi_ = np.where((yy[1:-1] > yy[:-2]) & (yy[1:-1] >= yy[2:]) & (yy[1:-1] > 0))[0] + 1
    if len(pi_) < 4:
        return None
    pk, tp = yy[pi_], tt[pi_]; sel = pk > 0.05 * pk.max()
    if sel.sum() < 4:
        return None
    return float(-np.polyfit(tp[sel], np.log(pk[sel]), 1)[0] / (2.0 * math.pi * f0))

t_mid, acc_mid = read_hist("ACC_MID")
t_wp, wp_mid = read_hist("WP_MID")
t_pot, pot = read_hist("POT11")
f_id, f_rel, n_cyc, xi = None, None, 0, None
if acc_mid is not None:
    tr, yr = ring_down(t_mid, acc_mid, PULSE_SEC + 0.02)
    f_zc, n_cyc = f_zero_cross(tr, yr); f_ff = f_fft(tr, yr)
    f_id = f_zc if f_zc else f_ff
    if f_zc and abs(f_zc / f_ff - 1) > 0.05:
        print("  ! zero-crossing ({:.2f}) and FFT ({:.2f}) differ by {:.0f}% -- more than one mode; using the FFT peak".format(
            f_zc, f_ff, 100 * abs(f_zc / f_ff - 1)))
        f_id = f_ff
    xi = log_decrement(tr, yr, f_id)
    print("\n  ring-down, mid-height acceleration ({:.2f} s, {} cycles counted)".format(tr[-1] - tr[0], n_cyc))
    print("    zero-crossing period  : {}".format("{:.2f} Hz  (T = {:.4f} s)".format(f_zc, 1 / f_zc) if f_zc else "too few cycles"))
    print("    zero-padded FFT peak  : {:.2f} Hz  (T = {:.4f} s)".format(f_ff, 1 / f_ff))
    if xi is not None:
        print("    damping, log decrement: {:.2f}% of critical{}".format(100 * xi, "   (undamped save)" if xi < 0.005 else ""))
else:
    print("\n  ! ACC_MID not found in {}. Export is by NAME against history {}; if 3DEC numbered".format(OUT_DIR, TIME_HIST_ID))
    print("    the histories differently, fix TIME_HIST_ID (`history list`).")
if wp_mid is not None and pot is not None:
    nn = min(len(wp_mid), len(pot))
    tr, rel = ring_down(t_wp[:nn], wp_mid[:nn] - pot[:nn], PULSE_SEC + 0.02)
    f_rel = f_fft(tr, rel)
    print("  mid-height relative displacement, FFT : {:.2f} Hz   (independent check)".format(f_rel))

print("\n  MEASURED (paper Table 6)               : {:.2f} Hz   {}".format(F_TARGET, REF))
if f_id:
    err = 100.0 * (f_id / F_TARGET - 1.0)
    print("  model / measured                       : {:+.1f}%".format(err))
    print("  " + "-" * 68)
    if abs(err) <= 5:
        print("  WITHIN 5%: this save reproduces the measured small-amplitude stiffness.")
    else:
        print("  first-order correction on top of this build: kn_scale x {:.2f}".format((F_TARGET / f_id) ** 2))
        print("  (all-in-joints estimate; the blocks carry ~half the compliance, so the true")
        print("   factor is further from 1 than that -- interpolate between the sweep points)")

with open(os.path.join(OUT_DIR, "idfreq_summary.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["tag", "start_save", "f_Hz", "T_s", "f_rel_disp_Hz", "cycles", "xi_logdec",
                "f_measured_Hz", "T_measured_s", "pct_error", "reference", "pulse"])
    w.writerow([TAG, START_SAVE, round(f_id, 3) if f_id else "", round(1.0 / f_id, 5) if f_id else "",
                round(f_rel, 3) if f_rel else "", n_cyc, round(xi, 4) if xi is not None else "",
                F_TARGET, round(1.0 / F_TARGET, 5), round(100.0 * (f_id / F_TARGET - 1.0), 2) if f_id else "", REF,
                "{:.0f}ms {:.2f}g + {:.1f}s".format(PULSE_SEC * 1000, PULSE_G, RING_SEC)])
print("\n-> {}/idfreq_summary.csv".format(OUT_DIR))
