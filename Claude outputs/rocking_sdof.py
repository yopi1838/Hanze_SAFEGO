# -*- coding: ascii -*-
"""
rocking_sdof.py -- a rocking single-degree-of-freedom oscillator built from a
capacity curve, and the amplitude rule that uses it.

WHY
    Matching the spectral displacement of a LINEAR oscillator at one period
    does not fix what a rocking wall does: on the cracked US-1 model the same
    Sd(0.36 s) gave 85 mm with one sine cycle and 5 mm with three. A rocking
    wall answers the size and coherence of the individual pushes relative to
    its own strength, not Sd. So the reduced model used to size the pulse
    has to be a rocking oscillator, not a linear one.

THE OSCILLATOR
    m u'' + c u' + F(u) = -m a_g(t)
    F(u)  : the wall's capacity curve (base shear vs control displacement),
            piecewise linear, taken from a pushover CSV (d_ctrl_mm, F_tot_kN
            or F_base_kN) or from the digitised Fig 13 envelope (d_mm, F_kN).
            Loading and unloading follow the same curve (nonlinear elastic --
            Housner's rocking idealisation: no hysteresis except at impact).
            Beyond the last point the curve is continued along its last slope
            down to F = 0 (the instability displacement); reaching it is
            reported as collapse.
    m     : the effective mass M_EFF (same value as the capacity law), so
            F/m is the mass-proportional acceleration the pushover applied.
    c     : small viscous damping XI on the initial stiffness (numerical
            robustness; keep <= 0.02).
    impact: at every zero crossing of u the velocity is multiplied by the
            restitution coefficient R (Housner). R is the one calibration
            parameter; estimate it from a ring-down of the 3DEC model with
            restitution_from_ringdown(), default 0.95.
    Integration: explicit central difference with sub-stepping.

THE AMPLITUDE RULE  (amplitude_for_state)
    Given the current pulse period T (the state, e.g. the secant period at
    the displacement the wall has reached), the cycle count N (fixed by the
    law) and the record a_g(t) x scale:
        d_rec  = peak |u| of the oscillator under the record
        Vamp   = velocity amplitude of  v(t) = Vamp sin(2 pi t / T), N cycles,
                 such that the oscillator's peak |u| under that pulse = d_rec
    solved by bisection. N sits inside the oscillator's response on both
    sides, so the rule does not depend on N the way an Sd match does; the
    check of THAT claim is the 3DEC run, not this module.

SELF-TESTS run on import (SELFTEST = True): linear limit reproduces the
Newmark Sd of a damped linear SDOF; restitution damps a free vibration;
the bisection recovers a known amplitude. numpy only (3DEC-Python safe).
"""
import math
import numpy as np

SELFTEST = True
G = 9.81


# ---------------------------------------------------------------- backbone
class Backbone(object):
    """Piecewise-linear F(d) in N and m, from a CSV. Symmetric in sign."""

    def __init__(self, csv_path, branch="neg", force_col=None, label=""):
        self.path = csv_path; self.label = label or csv_path
        d, F = self._load(csv_path, branch, force_col)
        self.d, self.F = d, F                      # m, N, ascending d, d[0] = 0
        self.K1 = float(F[1] / d[1]) if len(d) > 1 and d[1] > 0 else float("nan")
        # continuation past the last point along the last slope to F = 0
        k_end = (F[-1] - F[-2]) / (d[-1] - d[-2]) if len(d) > 2 else -F[-1] / d[-1]
        if k_end >= 0:                              # curve still rising/flat: hold flat, no instability
            self.d_inst = float("inf"); self.k_end = 0.0
        else:
            self.d_inst = float(d[-1] - F[-1] / k_end); self.k_end = float(k_end)
        self.F_peak = float(F.max()); self.d_peak = float(d[int(np.argmax(F))])

    @staticmethod
    def _load(path, branch, force_col):
        with open(path) as f:
            hdr = [h.strip() for h in f.readline().split(",")]
        raw = np.genfromtxt(path, delimiter=",", skip_header=1)
        if raw.ndim == 1:
            raw = raw[None, :]
        col = {h: raw[:, i] for i, h in enumerate(hdr)}
        if "d_ctrl_mm" in col:                                  # pushover driver output
            fc = force_col or ("F_tot_kN" if "F_tot_kN" in col and np.isfinite(col["F_tot_kN"]).any() else "F_base_kN")
            d = np.abs(col["d_ctrl_mm"]) / 1000.0; F = np.abs(col[fc]) * 1000.0
            if "converged" in col:                              # drop unconverged / unsettled points
                ok = col["converged"] >= 0.5; d, F = d[ok], F[ok]
        elif "d_mm" in col:                                     # digitised Fig 13 envelope
            m = col["d_mm"] < 0 if branch == "neg" else col["d_mm"] > 0
            d = np.abs(col["d_mm"][m]) / 1000.0; F = np.abs(col["F_kN"][m]) * 1000.0
        else:
            raise RuntimeError("backbone: unknown CSV format " + path)
        o = np.argsort(d); d, F = d[o], F[o]
        # make d strictly increasing (keep the max F at duplicate d)
        du, idx = np.unique(d, return_index=True); Fu = np.array([F[d == x].max() for x in du])
        if du[0] > 0:
            du = np.r_[0.0, du]; Fu = np.r_[0.0, Fu]
        return du, Fu

    def force(self, u):
        """N, odd in u; continued along the last slope to zero past the curve end."""
        a = abs(u)
        if a <= self.d[-1]:
            f = np.interp(a, self.d, self.F)
        else:
            f = max(0.0, self.F[-1] + self.k_end * (a - self.d[-1]))
        return math.copysign(f, u) if u != 0 else 0.0

    def secant_period(self, d_m, m_eff):
        """T = 2 pi sqrt(m / (F(d)/d)) at displacement d_m [m]."""
        if d_m <= 0:
            return 2 * math.pi * math.sqrt(m_eff / self.K1)
        Fd = abs(self.force(d_m))
        return float("inf") if Fd <= 0 else 2 * math.pi * math.sqrt(m_eff * d_m / Fd)


# ---------------------------------------------------------------- oscillator
def rocking_response(ag, dt, backbone, m_eff, restitution=0.95, xi=0.02, substeps=10,
                     u0=0.0, v0=0.0, collapse_factor=1.0):
    """Peak |u| [m] and the u(t) history of the rocking SDOF under ground
    acceleration ag [m/s2] sampled at dt. Returns (peak_m, u_hist, collapsed)."""
    ag = np.asarray(ag, float)
    h = dt / substeps
    w1 = math.sqrt(backbone.K1 / m_eff); c = 2.0 * xi * w1 * m_eff
    u, v = u0, v0
    acc = (-m_eff * ag[0] - c * v - backbone.force(u)) / m_eff
    peak = abs(u); hist = np.empty(len(ag)); collapsed = False
    d_col = backbone.d_inst * collapse_factor
    for i in range(len(ag)):
        a_i = ag[i]
        for _ in range(substeps):
            v_half = v + 0.5 * h * acc
            u_new = u + h * v_half
            if (u_new > 0) != (u > 0) and u != 0.0:      # zero crossing: impact
                v_half *= restitution
            u = u_new
            acc = (-m_eff * a_i - c * v_half - backbone.force(u)) / m_eff
            v = v_half + 0.5 * h * acc
            if abs(u) > peak:
                peak = abs(u)
        hist[i] = u
        if abs(u) >= d_col:
            collapsed = True; hist[i + 1:] = u; break
    return peak, hist, collapsed


def pulse_accel(Vamp, T, ncyc, dt):
    """Ground acceleration of v(t) = Vamp sin(2 pi t/T) over ncyc cycles (T snapped to dt)."""
    spc = max(1, int(round(T / dt))); Tg = spc * dt
    n = int(round(ncyc * spc)); w = 2 * math.pi / Tg
    t = np.arange(n + 1) * dt
    a = Vamp * w * np.cos(w * t)          # exact derivative of the velocity sine
    return a, Tg


def amplitude_for_state(ag_record, dt, T, ncyc, backbone, m_eff, restitution=0.95, xi=0.02,
                        tail_s=3.0, vmax=3.0, tol=0.01, itmax=40):
    """Vamp so that the oscillator's peak under the N-cycle sine at period T
    equals its peak under ag_record. Returns dict with Vamp, T (snapped), d_rec,
    d_pulse, pga_g, collapsed_record, converged."""
    d_rec, _, col_rec = rocking_response(ag_record, dt, backbone, m_eff, restitution, xi)
    ntail = int(round(tail_s / dt))
    def d_pulse(V):
        a, Tg = pulse_accel(V, T, ncyc, dt)
        a = np.r_[a, np.zeros(ntail)]
        p, _, col = rocking_response(a, dt, backbone, m_eff, restitution, xi)
        return p, col, Tg
    lo, hi = 0.0, vmax
    p_hi, col_hi, Tg = d_pulse(hi)
    if col_rec:
        return dict(Vamp=float("nan"), T=Tg, d_rec=d_rec, d_pulse=float("nan"), pga_g=float("nan"),
                    collapsed_record=True, converged=False,
                    note="oscillator collapses under the record: no finite target")
    if p_hi < d_rec and not col_hi:
        return dict(Vamp=float("nan"), T=Tg, d_rec=d_rec, d_pulse=p_hi, pga_g=float("nan"),
                    collapsed_record=False, converged=False,
                    note="pulse cannot reach d_rec below vmax={} m/s".format(vmax))
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        p, col, Tg = d_pulse(mid)
        if col or p > d_rec:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol * max(hi, 1e-6):
            break
    V = lo if lo > 0 else hi
    p, col, Tg = d_pulse(V)
    return dict(Vamp=V, T=Tg, d_rec=d_rec, d_pulse=p, pga_g=V * 2 * math.pi / Tg / G,
                collapsed_record=False, converged=abs(p - d_rec) <= 0.05 * d_rec, note="")


def restitution_from_ringdown(u, n_peaks=4):
    """Estimate R from successive absolute peaks of a free-vibration history:
    R ~ (u_{k+1}/u_k) per half cycle, averaged over n_peaks. Crude but model-based."""
    u = np.asarray(u, float) - np.mean(u[-max(10, len(u) // 20):])
    s = np.sign(u); idx = np.where(np.diff(s) != 0)[0]
    peaks = [np.max(np.abs(u[idx[k]:idx[k + 1]])) for k in range(len(idx) - 1)]
    peaks = [p for p in peaks if p > 0][:n_peaks + 1]
    if len(peaks) < 2:
        return float("nan")
    r = [peaks[k + 1] / peaks[k] for k in range(len(peaks) - 1)]
    return float(np.clip(np.mean(r), 0.5, 1.0))


# ---------------------------------------------------------------- self-tests
def _linear_sd(ag, dt, T, xi):
    w = 2 * math.pi / T; c = 2 * xi * w; w2 = w * w; u = v = 0.0; a = -ag[0]; umax = 0.0
    for i in range(1, len(ag)):
        un = u + dt * v + 0.25 * dt * dt * a; vn = v + 0.5 * dt * a
        an = (-ag[i] - c * vn - w2 * un) / (1 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = un + 0.25 * dt * dt * an; v = vn + 0.5 * dt * an; a = an; umax = max(umax, abs(u))
    return umax


def _selftest():
    import tempfile, os
    # linear backbone (K = 1 MN/m, up to 1 m, then flat -> no instability)
    p = os.path.join(tempfile.gettempdir(), "_bb_lin.csv")
    with open(p, "w") as f:
        f.write("d_mm,F_kN\n0,0\n1000,1000\n")
    bb = Backbone(p, branch="pos")
    m = 1000.0; T = 2 * math.pi * math.sqrt(m / bb.K1)          # 0.1987 s
    dt = 0.002; t = np.arange(0, 6.0, dt); ag = 2.0 * np.sin(2 * math.pi / 0.3 * t)
    pk, _, col = rocking_response(ag, dt, bb, m, restitution=1.0, xi=0.05)
    ref = _linear_sd(ag, dt, T, 0.05)
    assert not col and abs(pk - ref) / ref < 0.03, ("linear limit", pk, ref)
    # restitution damps a free vibration
    pk1, h1, _ = rocking_response(np.zeros(2000), dt, bb, m, restitution=1.0, xi=0.0, u0=0.01)
    pk2, h2, _ = rocking_response(np.zeros(2000), dt, bb, m, restitution=0.8, xi=0.0, u0=0.01)
    assert abs(h1[-1]) > 0.5 * abs(h2[-1]) and np.max(np.abs(h2[1500:])) < 0.3 * 0.01, "restitution"
    # bisection recovers a known amplitude: target = pulse with V=0.3 at T=0.3, N=1
    a_t, Tg = pulse_accel(0.3, 0.3, 1, dt); a_t = np.r_[a_t, np.zeros(1500)]
    r = amplitude_for_state(a_t, dt, 0.3, 1, bb, m, restitution=1.0, xi=0.05)
    assert r["converged"] and abs(r["Vamp"] - 0.3) / 0.3 < 0.03, ("bisection", r)
    os.remove(p)

if SELFTEST:
    _selftest()
