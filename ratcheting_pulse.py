# -*- coding: ascii -*-
"""
Ratcheting (asymmetric) spectrum-matched pulse  --  Modification A
==================================================================
Drop-in replacement for the symmetric pulse generation in the
strategy-C 3DEC driver. Implements an asymmetric base-acceleration
waveform whose dominant-direction peak still matches the target
spectral displacement, but whose reduced reverse lobes leave a net
per-cycle drift (ratcheting), in the direction dictated by the record.

Reduces EXACTLY to the original symmetric method when beta = 1.

This module is pure NumPy (no itasca) so it can be unit-tested outside
3DEC. The driver imports calibrate_amplitude_asym and build_velocity_asym
and calls them in place of the originals.

Asymmetry is sourced from the record, not tuned:
    s    = sign of the record's dominant acceleration peak (+1 or -1)
    beta = |peak_reverse| / |peak_dominant|   in (0, 1]
A single optional global constant k (default 1.0, applied identically to
every run/record) can sharpen the asymmetry if the record's own value
under-ratchets; it is the only free parameter and must be validated
out of sample (see ratcheting_pulse_spec.md).

NOTE on sign convention: s is expressed in the record (channel-12 table
accelerometer) positive direction. Before use, confirm this maps to the
3DEC velocity-z positive direction; if the model's +z is opposite to the
sensor +x, negate s. The lean-direction check (HU12 vs FR76) fixes this.
"""

import os, math
import numpy as np


# =====================================================================
# Record-derived asymmetry (measured on the achieved table motion,
# strong-motion window; recompute with asym_from_record for exactness)
# =====================================================================
RECORD_ASYM = {
    "HU12": {"s": +1, "beta": 0.54},
    "EC40": {"s": +1, "beta": 0.70},
    "FR76": {"s": -1, "beta": 0.78},
}


def asym_from_record(a_g, pga_frac=0.02):
    """
    Compute (s, beta) from a raw acceleration array (any units).
    Trims to the strong-motion window (|a|>pga_frac*PGA), de-means,
    and returns the dominant sign and the reverse/dominant peak ratio.
    """
    a = np.asarray(a_g, float)
    a = a - np.mean(a)
    env = np.abs(a) > pga_frac * np.abs(a).max()
    idx = np.where(env)[0]
    if len(idx):
        a = a[idx[0]:idx[-1] + 1]
    ap, an = a.max(), abs(a.min())
    if ap >= an:
        return +1, float(an / ap)
    return -1, float(ap / an)


# =====================================================================
# Asymmetric unit waveform  (dominant peak ~ 1 in direction s)
# =====================================================================
def asym_shape(t, w, beta, s, demean=True):
    """
    Acceleration shape g(t): cosine with the lobes opposite to the
    dominant direction s scaled by beta. With beta==1 this is cos(w t).
    Dominant lobes keep full amplitude; reverse lobes are scaled by beta,
    so the inertial demand is biased toward direction s.

    demean=True subtracts the mean so the integral of acceleration over
    the pulse is zero, i.e. the base returns to zero net velocity (no
    runaway base translation across sequential runs). For beta==1 over a
    half-integer number of cycles the mean is already ~0, so the symmetric
    case is unchanged and the regression against the original is preserved.
    """
    base = np.cos(w * t)
    g = base.copy()
    opp = (s * base) < 0.0          # lobes opposite to dominant direction
    g[opp] = beta * base[opp]
    if demean:
        # Remove ONLY the net offset introduced by the asymmetry, so the
        # base returns to zero net velocity (no runaway translation across
        # runs) while beta==1 stays exactly cos(w t) and reproduces the
        # original symmetric method.
        delta = g - base
        g = base + (delta - np.mean(delta))
    return g


def beta_eff(T_end, T1, beta_record, ramp_lo=1.05, ramp_hi=1.20):
    """
    Two-regime activation: beta = 1 (symmetric) while the wall is
    essentially elastic, ramping linearly to beta_record as the
    period-elongation ratio T_end/T1 goes from ramp_lo to ramp_hi.
    Keeps the early, well-matched runs untouched.
    """
    r = T_end / T1
    if r <= ramp_lo:
        return 1.0
    if r >= ramp_hi:
        return beta_record
    frac = (r - ramp_lo) / (ramp_hi - ramp_lo)
    return 1.0 + frac * (beta_record - 1.0)


# =====================================================================
# Newmark linear SDOF peak relative displacement (average accel.)
# (identical to the driver's routine; included for standalone testing)
# =====================================================================
def newmark_sd(a_g, dt, T, xi):
    if T <= 0:
        return 0.0
    w = 2 * math.pi / T
    c = 2 * xi * w
    w2 = w * w
    u = 0.0
    v = 0.0
    acc = -a_g[0]
    u_max = 0.0
    for i in range(1, len(a_g)):
        u_new = u + dt * v + dt**2 * 0.25 * acc
        v_new = v + dt * 0.5 * acc
        denom = 1.0 + 0.5 * dt * c + 0.25 * dt**2 * w2
        acc_new = (-a_g[i] - c * v_new - w2 * u_new) / denom
        u_new += 0.25 * dt**2 * acc_new
        v_new += 0.5 * dt * acc_new
        u, v, acc = u_new, v_new, acc_new
        if abs(u) > u_max:
            u_max = abs(u)
    return u_max


# =====================================================================
# Calibration: amplitude so dominant-direction Sd(T) == Sd_target
# =====================================================================
def calibrate_amplitude_asym(T, Sd_target, beta, s, n_cycles, delta_t, xi, k=1.0):
    """
    Build the unit asymmetric waveform, compute its elastic Sd at T, and
    scale so the dominant-direction peak displacement equals Sd_target.
    Linear in amplitude because beta fixes the shape.
    k sharpens asymmetry: beta_used = beta**k (k>=1 => more asymmetric).
    Returns (A, beta_used).
    """
    beta_used = beta ** k
    w = 2.0 * math.pi / T
    n = int(round(n_cycles * T / delta_t)) + 1
    t = np.arange(n) * delta_t
    g = asym_shape(t, w, beta_used, s)
    sd0 = newmark_sd(g, delta_t, T, xi)
    if sd0 == 0.0:
        raise RuntimeError("Trial Sd=0 at T={:.4f}".format(T))
    return (Sd_target / sd0), beta_used


# =====================================================================
# Velocity-history file (asymmetric acceleration integrated to velocity)
# =====================================================================
def build_velocity_asym(A, T, beta, s, run_no, out_dir,
                        n_cycles, delta_t, tail_sec):
    """
    a(t) = A * asym_shape ; v = cumulative integral with v(0)=0 ;
    raised-cosine taper over half a period to v=0 ; zero tail.
    Same file format as the original build_velocity_file.
    Returns (fpath, duration, v_peak).
    """
    w = 2.0 * math.pi / T
    n = int(round(n_cycles * T / delta_t)) + 1
    t = np.arange(n) * delta_t

    a = A * asym_shape(t, w, beta, s)                 # acceleration (m/s2)
    v = np.concatenate(([0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * delta_t)))

    # raised-cosine taper from v_end to 0 over half a period
    ramp_sec = max(0.5 * T, delta_t)
    n_tail = int(round(tail_sec / delta_t))
    n_ramp = min(int(round(ramp_sec / delta_t)), n_tail)

    t_tail = t[-1] + delta_t + np.arange(n_tail) * delta_t
    t = np.r_[t, t_tail]

    v_end = float(v[-1])
    if n_ramp > 1:
        ss = np.linspace(0.0, 1.0, n_ramp)
        wcos = 0.5 * (1.0 + np.cos(math.pi * ss))
        v_ramp = v_end * wcos
    else:
        v_ramp = np.array([0.0])

    v_tail = np.zeros(n_tail)
    v_tail[:len(v_ramp)] = v_ramp
    v = np.r_[v, v_tail]

    fname = "vel_run_{:02d}.txt".format(run_no)
    fpath = os.path.join(out_dir, fname)
    N = len(t)
    with open(fpath, "w") as f:
        f.write("StratC_run{:02d}_T{:.4f}_beta{:.3f}_s{:+d}\n".format(run_no, T, beta, s))
        f.write("{}\t0\n".format(N))
        for ti, vi in zip(t, v):
            f.write("{:.6f}\t{:.9e}\n".format(ti, vi))

    return fpath, float(t[-1]), float(np.max(np.abs(v)))