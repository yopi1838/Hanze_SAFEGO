"""Large-displacement virtual work for the two-block (three-hinge) strip under a mass-proportional lateral load.

Geometry and kinematics (push towards +z):
    base pivot B at the leeward toe (t, 0); mid hinge M on the windward face (0, h1); top contact at the leeward top
    corner (t, h), held sideways (the beam), free to rise. Overburden O (kN) rides on the top contact.

TOP AND BASE OPTIONS (23 Oct, from the DEM deformed shapes):
    top=None            top pivot at the leeward top corner (t, h): the wall rocks on the top beam (mason top joint);
    top=(z_p, y_p)      top pivot at a fixed point of the upper block, held sideways at z_p and free to rise -- e.g.
                        (t/2, 2.65) for the elastic (glued) top joint, where the top beam rotates with the upper block
                        about its own centroid (DEM v5 elastic-top runs: the top joint tilts 140-194 mrad, stays closed);
                        the overburden O acts at that point.
    base_y=0.0          height of the lower pivot (leeward face). base_y > 0 treats the masonry below it as not moving
                        (the DEM lower block kinks at the lower floor, 0.3-0.5 m; base_y = 0.30 bounds that restraint).

Equilibrium by potential energy: the push a*m does work through the mass-weighted lateral displacement Z, gravity and
the joists store U. F = dU/dtheta / dZ/dtheta is the total lateral force on the WALL mass (kN) at equilibrium, i.e.
F = a * M_wall. Everything is evaluated in the displaced geometry, so the curve is exact for rigid blocks.

JOISTS (optional, list of Joist) -- each one is a floor (both joists of a floor lumped):
    R_wall   bearing load the pocket puts on the wall (kN), applied at the bearing point (z_b, y) of the block it sits in
    W_lever  joist weight carried by the far support at rest (kN)
    k_r      rotational clamp of the joist in its pocket (N m/rad); the far support acts as a vertical spring
             k_v = k_r / L^2 under the far end, preloaded by W_lever, compression only:
               far end pushed DOWN -> spring force grows       -> rotational restraint (clamping)
               far end lifted       -> spring unloads, then the joist hangs off the wall -> constant restoring moment
             k_r = 0: pinned joist (rotates freely in the pocket, as on the test's hinge) -> only R_wall acts
    M_cap    largest moment the pocket can transmit (kN m), in both senses; beyond it the joist turns in the pocket
             at a constant moment (elastic-perfectly-plastic clamp). Monotonic loading assumed (pushover).
    m_floor  floor mass (kg) that the push also accelerates (inertia through the pocket), used if inertia=True
    side     "leeward" (push towards the joists) or "windward" (push away from them)
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

t, h, zs = 0.21, 2.58, 2.06
G = 9.81


def rot(p, c, a):           # rotate point p about c by a (a>0 tilts the up-vector toward +z)
    dz, dy = p[0] - c[0], p[1] - c[1]
    return np.array([c[0] + dz * np.cos(a) + dy * np.sin(a), c[1] - dz * np.sin(a) + dy * np.cos(a)])


def config(th1, h1, top=None, base_y=0.0):
    P0 = np.array([t, base_y])
    H = rot(np.array([0.0, h1]), P0, th1)
    ptop = np.array([t, h]) if top is None else np.array(top, dtype=float)

    def up(p, th2):          # upper-block point p (original coords), block rotated by -th2 about its (0,h1)
        q = rot(p, np.array([0.0, h1]), -th2)
        return q - np.array([0.0, h1]) + H
    f = lambda th2: up(ptop, th2)[0] - ptop[0]
    th2 = brentq(f, -1.0, 1.5) if th1 > 0 else 0.0
    lo = lambda p: rot(p, P0, th1) if p[1] >= base_y else np.array(p, dtype=float)
    return lo, (lambda p: up(p, th2)), th2


@dataclass
class Joist:
    y: float                 # height of the joist bearing (m)
    R_wall: float            # kN, joist weight carried by the wall through the pocket
    W_lever: float           # kN, joist weight carried by the far support at rest
    L: float = 1.02          # m, from the bearing to the far support
    k_r: float = 0.0         # N m/rad, pocket clamp (0 = pinned)
    M_cap: float = float("inf")   # kN m, largest clamp moment the pocket can transmit (pressing side)
    m_floor: float = 0.0     # kg, floor mass accelerated with the wall
    side: str = "leeward"    # "leeward": push towards the joists; "windward": push away from them
    pocket: float = 0.5      # pocket depth as a fraction of t (US-1: 105 of 210 mm)

    def zb(self):            # bearing point: middle of the pocket, on the joist side
        d = self.pocket * t / 2
        return t - d if self.side == "leeward" else d

    def zf(self):            # far support, beyond the joist-side face
        return self.zb() + self.L if self.side == "leeward" else self.zb() - self.L


def far_end_energy(v, kv, W_lever, f_cap):
    """work stored by the far support + the lever weight for a far-end rise v (m), kN m. The far end feels
    g(v) = clip(kv v, -f_cap, min(W_lever, f_cap)):
      pushed down -> clamp reaction, capped at f_cap = M_cap/L (the pocket lets the joist turn beyond M_cap);
      lifted      -> the support unloads, then the joist weight W_lever hangs off the wall -- also only as far as
                     the pocket can transmit the moment (min(W_lever, f_cap))."""
    up_cap = min(W_lever, f_cap)
    vu = up_cap / kv
    vc = f_cap / kv
    if v >= vu:
        return 0.5 * kv * vu ** 2 + up_cap * (v - vu)
    if v >= -vc:
        return 0.5 * kv * v ** 2
    return 0.5 * kv * vc ** 2 + f_cap * (-v - vc)


def state(th1, h1, W, O, kr=0.0, joists=(), inertia=False, top=None, base_y=0.0):
    lo, up, th2 = config(th1, h1, top, base_y)
    m1, m2 = (h1 - base_y) / h, (h - h1) / h
    c1, c2 = lo(np.array([t / 2, (h1 + base_y) / 2])), up(np.array([t / 2, (h1 + h) / 2]))
    ytop = up(np.array([t, h]) if top is None else np.array(top, dtype=float))[1]
    U = W * (m1 * c1[1] + m2 * c2[1]) + O * ytop + 0.5 * kr * th1 ** 2 / 1e3   # kN*m
    Z = m1 * (c1[0] - t / 2) + m2 * (c2[0] - t / 2)                            # mass-weighted sway / M_wall
    M_wall = W * 1e3 / G
    for j in joists:
        blk = lo if j.y <= h1 else up
        pb = blk(np.array([j.zb(), j.y]))
        U += j.R_wall * pb[1]                                                # (1) bearing load
        if j.k_r > 0:                                                        # (3)+(4) clamp / lift-off
            kv = j.k_r / 1e3 / j.L ** 2                                      # kN/m at the far end
            v = blk(np.array([j.zf(), j.y]))[1] - j.y                        # far-end rise (m)
            U += far_end_energy(v, kv, j.W_lever, j.M_cap / j.L)
        if inertia and j.m_floor > 0:                                        # floor mass pushed with the wall
            Z += j.m_floor / M_wall * (pb[0] - j.zb())
    ps = lo(np.array([t, zs])) if zs <= h1 else up(np.array([t, zs]))
    return U, Z, ps[0] - t


def curve(h1, W, O, kr=0.0, thmax=0.12, n=400, joists=(), inertia=False, top=None, base_y=0.0):
    th = np.linspace(1e-5, thmax, n)
    U, Z, d = np.array([state(a, h1, W, O, kr, joists, inertia, top, base_y) for a in th]).T
    F = np.gradient(U, th) / np.gradient(Z, th)
    return d * 1e3, F


def us1_joists(side, k_r=1.3e5, M_cap=float("inf"), m_floor=203.0, R_wall=0.186 * 2, W_lever=(0.996 - 0.186) * 2):
    """US-1 DEM build: two floors (2 joists each, 101.5 kg per joist); the wall carries 76 kg of the 406 kg at rest.
    side: 'leeward' for a push towards the joists (positive z in the 3DEC model), 'windward' for the negative push."""
    return [Joist(0.38, R_wall, W_lever, k_r=k_r, M_cap=M_cap, m_floor=m_floor, side=side),
            Joist(2.14, R_wall, W_lever, k_r=k_r, M_cap=M_cap, m_floor=m_floor, side=side)]
