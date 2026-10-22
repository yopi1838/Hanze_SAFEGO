# -*- coding: ascii -*-
"""
ringdown_check.py (v2) -- verify the Maxwell damping level on the corrected wall.

v1 result: three units at 0.015 each (1.5 / 6 / 25 Hz) measured xi = 0.0266 at 10 Hz, because the
units superpose: xi(f) = sum_i xi_i * 2(f/f_i) / (1 + (f/f_i)^2) = 1.86 x 0.015 at 10 Hz.
v2: each unit at 0.0085 -> predicted 0.0155 at 10 Hz, 0.0153 at 3 Hz, 0.013 at 1.5 and 25 Hz.
Pass: xi 0.014-0.017, frequency within ~2 % of the undamped 9.7 Hz, no residual, no beam drift.

Pulse: one-cycle velocity sine (PULSE_T) with base PGA = PULSE_PGA_G; free vibration TAIL_SEC.
Needs the save and the three instrument .dat files the sequence drivers use.
Output: <OUT_DIR>/<case>_*.csv and <OUT_DIR>/ringdown_summary.csv
"""
import os, csv, math
import numpy as np
import itasca as it

# ================================ CONFIG ================================
BASE_SAVE   = "Part_I_MASON_v8_LargeStrain_springforce.sav"
OUT_DIR     = "ringdown_check_v2"
HISTORY_SETUP   = "instrument_history_new.dat"
BASESHEAR_SETUP = "instrument_baseshear_exp.dat"
HISTORY_EXPORT  = "instrument_history_export_v2.dat"

LARGE_STRAIN = True
DAMP_RATIO_UNIT = 0.0085            # PER UNIT; the three units superpose to ~0.015 over 1.5-25 Hz
DAMP_FREQS      = (1.5, 6.0, 25.0)
CASES = [                            # (label, mode, per-unit ratio, freqs)
    ("maxwell_0085", "maxwell", DAMP_RATIO_UNIT, DAMP_FREQS),
    # ("undamped", "none", 0.0, None),   # reference, already measured: xi 0.0002, f 9.70 Hz
]
RAYLEIGH_FC = 6.0

PULSE_T     = 0.10
PULSE_PGA_G = 0.05
TAIL_SEC    = 3.0
DELTA_T     = 0.0025
DRIVE_GROUPS = ["S", "T_B"]
TB_POS      = (0.646, 2.63, 0.105)
NOISE_MM    = 0.01

os.makedirs(OUT_DIR, exist_ok=True)

def cmd_path(p): return os.path.normpath(p).replace("\\", "/")

def maxwell_sum(f, ratio, freqs):
    """predicted damping ratio at frequency f for three equal Maxwell units"""
    return sum(ratio * 2.0 * (f / fi) / (1.0 + (f / fi) ** 2) for fi in freqs)

def damping_command(mode, ratio, freqs):
    if mode == "maxwell":
        d = ratio; f1, f2, f3 = freqs
        return "block mechanical damping maxwell {:.4f} {:.3f} {:.4f} {:.3f} {:.4f} {:.3f}".format(d, f1, d, f2, d, f3)
    if mode == "rayleigh":
        return "block mechanical damping rayleigh {:.4f} {:.3f} mass".format(ratio, RAYLEIGH_FC)
    return "block mechanical damping local 0.0"

def restore_dynamic(save_file, dcmd):
    it.command("model restore '{}'".format(cmd_path(os.path.splitext(save_file)[0])))
    it.command("python-reset-state false")
    it.command("model large-strain {}".format("on" if LARGE_STRAIN else "off"))
    it.command("model dynamic active on")
    it.command("block mechanical damping local 0.0")
    it.command(dcmd)
    it.command("""
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
""")
    it.command("call '{}'".format(cmd_path(HISTORY_SETUP)))
    it.command("call '{}'".format(cmd_path(BASESHEAR_SETUP)))
    it.command("block free velocity-z range group 'S'")
    it.command("block free rotation-y range group 'T_B'")
    it.command("block free rotation-z range group 'T_B'")
    it.command("model history name 'rd_time' dynamic time-total")
    it.command("block history name 'rd_tb_uy' displacement-y position {} {} {}".format(*TB_POS))
    it.command("block history name 'rd_tb_uz' displacement-z position {} {} {}".format(*TB_POS))

def write_pulse_table(path):
    v_amp = PULSE_PGA_G * 9.81 * PULSE_T / (2.0 * math.pi)
    n = int(round(PULSE_T / DELTA_T))
    t = np.arange(n + 1) * DELTA_T
    v = v_amp * np.sin(2.0 * math.pi * t / PULSE_T)
    with open(path, "w", newline="\n") as f:
        f.write("ringdown_pulse\n{}\t0\n".format(len(t)))
        for ti, vi in zip(t, v): f.write("{:.6f}\t{:.9e}\n".format(ti, vi))
    return float(t[-1]), v_amp

def read_hist(path):
    d = np.genfromtxt(path, skip_header=2)
    m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    return d[m, 0], d[m, 1]

def analyse_tail(t, u, t_start):
    m = t >= t_start
    tt, uu = t[m], u[m]
    uu = uu - uu[-int(0.1 * len(uu)):].mean()
    s = np.sign(uu); idx = np.where(np.diff(s) != 0)[0]
    tz = tt[idx]; per = np.diff(tz)
    per = per[per > 0.25 * np.median(per)] if len(per) > 3 else per
    f_zc = 0.5 / per.mean() if len(per) else float("nan")
    dt = float(np.median(np.diff(tt))); n = len(uu)
    spec = np.abs(np.fft.rfft(uu * np.hanning(n))); fr = np.fft.rfftfreq(n, dt)
    f_fft = float(fr[1:][np.argmax(spec[1:])]) if n > 8 else float("nan")
    pk = [i for i in range(1, n - 1) if uu[i] > uu[i - 1] and uu[i] >= uu[i + 1] and uu[i] > NOISE_MM]
    amps = np.array([uu[i] for i in pk])
    if len(amps) >= 3:
        k = np.arange(len(amps)); slope = np.polyfit(k, np.log(amps), 1)[0]
        xi = -slope / (2.0 * math.pi)
        xi_pair = float(np.mean(np.log(amps[:-1] / amps[1:]))) / (2.0 * math.pi)
    else:
        xi = xi_pair = float("nan")
    return f_zc, f_fft, xi, xi_pair, len(amps), float(np.abs(uu).max())

summary = []
pulse_path = os.path.join(OUT_DIR, "ringdown_pulse_vel.txt")
dur, v_amp = write_pulse_table(pulse_path)
print("pulse: T {} s, V {:.4f} m/s, PGA {:.3f} g, {} s tail".format(PULSE_T, v_amp, PULSE_PGA_G, TAIL_SEC))

for label, mode, ratio, freqs in CASES:
    dcmd = damping_command(mode, ratio, freqs)
    print("\n=== case {}: {} ===".format(label, dcmd))
    if mode == "maxwell":
        print("  predicted xi: {:.4f} at 3 Hz, {:.4f} at 10 Hz, {:.4f} at 1.5 Hz, {:.4f} at 25 Hz".format(
            maxwell_sum(3.0, ratio, freqs), maxwell_sum(10.0, ratio, freqs),
            maxwell_sum(1.5, ratio, freqs), maxwell_sum(25.0, ratio, freqs)))
    restore_dynamic(BASE_SAVE, dcmd)
    it.command("history interval 5")
    it.command("table 'rd' import '{}'".format(cmd_path(pulse_path)))
    it.command("model dynamic time-total 0")
    for g in DRIVE_GROUPS:
        it.command("block apply velocity-z 1.0 table 'rd' range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(dur))
    for g in DRIVE_GROUPS:
        it.command("block gridpoint apply-remove velocity-z range group '{}'".format(g))
    it.command("model solve dynamic time {:.6f}".format(TAIL_SEC))

    it.command("[exportdir='{}']".format(cmd_path(OUT_DIR)))
    it.command("[runlabel='{}']".format(label))
    it.command("call '{}'".format(cmd_path(HISTORY_EXPORT)))
    for h in ("rd_tb_uy", "rd_tb_uz"):
        it.command("history export '{}' vs 'rd_time' file '{}' truncate".format(
            h, cmd_path(os.path.join(OUT_DIR, "{}_{}.csv".format(label, h)))))
    it.command("table 'rd' delete")

    t, u = read_hist(os.path.join(OUT_DIR, label + "_rel_disp_top_exp_mm.csv"))
    f_zc, f_fft, xi, xi_pair, npk, u_tail = analyse_tail(t, u, dur + 0.05)
    u_peak = float(np.abs(u - u[0]).max()); u_res = float(u[-int(0.05 * len(u)):].mean() - u[0])
    ty, uy = read_hist(os.path.join(OUT_DIR, "{}_rd_tb_uy.csv".format(label)))
    tb_drift = float((uy[-1] - uy[0]) * 1e3)
    f_mid = 0.5 * (f_zc + f_fft)
    xi_pred = maxwell_sum(f_mid, ratio, freqs) if mode == "maxwell" else 0.0
    row = dict(case=label, damping=dcmd, ratio_per_unit=ratio, f_zc_Hz=f_zc, f_fft_Hz=f_fft,
               xi_fit=xi, xi_pairs=xi_pair, xi_predicted_at_f=xi_pred, n_peaks=npk,
               u_peak_mm=u_peak, u_tail_peak_mm=u_tail, u_residual_mm=u_res, tb_vertical_drift_mm=tb_drift)
    summary.append(row)
    print("  frequency: {:.2f} Hz (zero-crossings), {:.2f} Hz (FFT)   [undamped reference 9.70 / 9.83]".format(f_zc, f_fft))
    print("  damping  : xi = {:.4f} (envelope), {:.4f} (pairs), {} peaks   [predicted {:.4f} at {:.1f} Hz; pass 0.014-0.017]".format(
        xi, xi_pair, npk, xi_pred, f_mid))
    print("  wall     : peak {:.3f} mm, residual {:+.4f} mm".format(u_peak, u_res))
    print("  top beam : vertical drift {:+.4f} mm".format(tb_drift))

with open(os.path.join(OUT_DIR, "ringdown_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
print("\n->", os.path.join(OUT_DIR, "ringdown_summary.csv"))