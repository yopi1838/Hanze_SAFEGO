# -*- coding: ascii -*-
"""
slip_work_tracker.py -- frictional slip work per contact REGION, tracked subcontact by subcontact (3DEC 9, Python).

Why: 'model energy mechanical on' only keeps model-wide totals. To test the mechanical model's pocket-friction term
(rigid_two_block_vw.py, POCKET CLAMP) we need the slip work and the clamp force of the LOWER POCKET on its own, and
the same for the upper pocket, the joist supports, the base, the top joint and the hinge courses.

WHAT IT ACCUMULATES, per subcontact, between two updates (old -> new):
    dW = s * 0.5 (F_old + F_new) . (u_new - u_old)          F = shear force vector, u = shear displacement vector
    s  = the sign convention of 3DEC's shear update, found from the data: over every increment in which |F| grows
         (elastic loading) F_mid . du must be positive in the right convention; s = sign of that sum.
    W_shear (region) = sum dW                                total shear work done on the contacts (J)
    E_el (region)    = sum |F|^2 / (2 k_s A)                 shear energy stored elastically NOW (J), k_s = stiffness-shear
    W_slip           = W_shear - (E_el - E_el at install)    dissipated by slip since install (friction / softening) (J)
  and, at each update, per region: N = sum of the normal forces (kN), the number of live subcontacts, and the largest
  accumulated shear-path length of a subcontact (mm).
  Subcontacts that disappear (large strain) keep what they accumulated; their last elastic energy is counted as
  dissipated. New subcontacts are classified when first seen.

REGIONS (contact groups in slot 'swt', set once at install with the same 'range group-intersection' commands the
pushover driver uses; masonry joints are split by height):
    pocket_low   Joist - Masonry, y < 1.3 m            pocket_up   Joist - Masonry, y >= 1.3 m
    support      Joist - Joist_S1 / Joist_S2           base        y in [-0.01, 0.01]
    top          T_B - Masonry                         new         formed after install (large strain)
    mas_lowfloor masonry joints 0.25-0.55 m            mas_upfloor masonry joints 2.03-2.25 m
    mas_other    every other masonry joint
  CAVEAT for 'support': those contacts have stiffness-shear 1 Pa/m and 1 deg friction, so their SHEAR work is ~0 by
  construction; any drag they exert comes through normal forces on inclined faces, which this tracker does not count.

USE (inside 3DEC, after the model is restored and settled):
    import slip_work_tracker as swt
    swt.install("slipwork_pos.csv")          # classify, zero the accumulators, write the header
    ... after every solve chunk:  swt.update(d_ctrl_mm)
    ... at a logged pushover point: swt.update(d_ctrl_mm, tag="asc") / tag="desc"
    swt.close()
Every update appends one row to the CSV (flushed), so an interrupted run keeps its data.
Accuracy: midpoint rule per update interval; exact for monotonic slip between updates (a pushover).
"""
import os, math, time
import itasca as it

REGIONS = ("pocket_low", "pocket_up", "support", "base", "top", "mas_lowfloor", "mas_upfloor", "mas_other", "new",
           "pocket_low_top", "pocket_low_bot", "pocket_low_side", "pocket_up_top", "pocket_up_bot", "pocket_up_side")
# the last six are the pocket FACES (also counted in pocket_low / pocket_up): top = ceiling (masonry above the joist),
# bot = floor, side = the vertical faces. Pocket heights in the build:
POCKET_Y = {"pocket_low": (0.29318, 0.46909), "pocket_up": (2.05227, 2.22818)}
FACE_TOL = 0.012
SLOT = "swt"
Y_POCKET_SPLIT = 1.3
LOWFLOOR = (0.25, 0.55)
UPFLOOR = (2.03, 2.25)
TRACK_MASONRY = True            # False: only pockets, supports, base and top (faster)

_S = {}                          # sc id -> [region, Fx, Fy, Fz, ux, uy, uz, ks_A or None, path]
_acc = {}                        # region -> dict(W=..., gone=...)
_sign = {"pos": 0.0, "neg": 0.0}
_csv = None
_t0 = None
_ks_ok = {"n": 0, "fail": 0}


def _v3(v):
    try:
        return float(v.x()), float(v.y()), float(v.z())
    except Exception:
        try:
            return float(v[0]), float(v[1]), float(v[2])
        except Exception:
            return float(v.x), float(v.y), float(v.z)


def _set_groups():
    cmds = [
        "block contact group 'swt_mas' slot '{}'".format(SLOT),
        "block contact group 'swt_pocket' slot '{}' range group-intersection 'Joist' 'Masonry'".format(SLOT),
        "block contact group 'swt_support' slot '{}' range group-intersection 'Joist' 'Joist_S1'".format(SLOT),
        "block contact group 'swt_support' slot '{}' range group-intersection 'Joist' 'Joist_S2'".format(SLOT),
        "block contact group 'swt_top' slot '{}' range group-intersection 'T_B' 'Masonry'".format(SLOT),
        "block contact group 'swt_base' slot '{}' range pos-y -0.01 0.01".format(SLOT),
    ]
    for c in cmds:
        it.command(c)


def _group_of(sc):
    """contact group of the subcontact's contact in slot 'swt' ('' if none)"""
    for fn in (lambda: sc.contact().group(SLOT), lambda: sc.group(SLOT),
               lambda: it.block.contact.find(sc.contact()).group(SLOT)):
        try:
            g = fn()
            if g is not None:
                return str(g)
        except Exception:
            continue
    return ""


def _region(sc, installed):
    g = _group_of(sc)
    y = float(sc.pos_y())
    if g == "swt_pocket":
        return "pocket_low" if y < Y_POCKET_SPLIT else "pocket_up"
    if g == "swt_support":
        return "support"
    if g == "swt_base":
        return "base"
    if g == "swt_top":
        return "top"
    if g == "swt_mas":
        if not TRACK_MASONRY:
            return None
        if LOWFLOOR[0] <= y <= LOWFLOOR[1]:
            return "mas_lowfloor"
        if UPFLOOR[0] <= y <= UPFLOOR[1]:
            return "mas_upfloor"
        return "mas_other"
    return "new" if installed else None


def _face(region, y):
    lo, hi = POCKET_Y.get(region, (None, None))
    if lo is None:
        return None
    if abs(y - hi) < FACE_TOL:
        return region + "_top"
    if abs(y - lo) < FACE_TOL:
        return region + "_bot"
    return region + "_side"


def _ksA(sc):
    try:
        ks = float(sc.prop("stiffness-shear"))
        A = float(sc.area())
        _ks_ok["n"] += 1
        return ks * A if ks * A > 0 else None
    except Exception:
        _ks_ok["fail"] += 1
        return None


def _snapshot(sc):
    F = _v3(sc.force_shear())
    u = _v3(sc.disp_shear())
    return F, u


def install(csv_path, note=""):
    global _csv, _t0
    _S.clear(); _acc.clear(); _sign["pos"] = _sign["neg"] = 0.0
    for r in REGIONS:
        _acc[r] = dict(W=0.0, gone=0.0, E0=0.0)
    _set_groups()
    n = 0
    for sc in it.block.subcontact.list():
        r = _region(sc, installed=False)
        if r is None:
            continue
        F, u = _snapshot(sc)
        _S[sc.id()] = [r, F[0], F[1], F[2], u[0], u[1], u[2], _ksA(sc), 0.0, _face(r, float(sc.pos_y()))]
        n += 1
    cols = ["time_s", "tag", "d_ctrl_mm", "sign"]
    for r in REGIONS:
        cols += [r + "_Wslip_J", r + "_Wshear_J", r + "_Eel_J", r + "_N_kN", r + "_n", r + "_path_max_mm"]
    d = os.path.dirname(csv_path)
    if d:
        os.makedirs(d, exist_ok=True)
    _csv = open(csv_path, "w", newline="\n")
    _csv.write("# slip_work_tracker {}; {} subcontacts tracked at install; k_s read for {} (failed {}); {}\n".format(
        time.strftime("%Y-%m-%d %H:%M"), n, _ks_ok["n"], _ks_ok["fail"], note))
    _csv.write(",".join(cols) + "\n"); _csv.flush()
    _t0 = time.time()
    counts = {}
    for v in _S.values():
        for k in (v[0], v[9]):
            if k:
                counts[k] = counts.get(k, 0) + 1
                if v[7]:
                    _acc[k]["E0"] += (v[1] ** 2 + v[2] ** 2 + v[3] ** 2) / (2.0 * v[7])
    print("[slip_work_tracker] installed: {} subcontacts; per region {}; k_s read {} / failed {}".format(
        n, counts, _ks_ok["n"], _ks_ok["fail"]))
    if counts.get("pocket_low", 0) == 0:
        print("[slip_work_tracker] WARNING: no lower-pocket subcontacts found -- check the contact groups")
    return counts


def update(d_ctrl_mm=float("nan"), tag="chunk"):
    if _csv is None:
        return
    live = set()
    N = dict((r, 0.0) for r in REGIONS)
    E = dict((r, 0.0) for r in REGIONS)
    cnt = dict((r, 0) for r in REGIONS)
    pmax = dict((r, 0.0) for r in REGIONS)
    for sc in it.block.subcontact.list():
        i = sc.id()
        st = _S.get(i)
        if st is None:
            r = _region(sc, installed=True)
            if r is None:
                continue
            F, u = _snapshot(sc)
            st = _S[i] = [r, F[0], F[1], F[2], u[0], u[1], u[2], _ksA(sc), 0.0, _face(r, float(sc.pos_y()))]
        r = st[0]
        F, u = _snapshot(sc)
        du = (u[0] - st[4], u[1] - st[5], u[2] - st[6])
        Fm = (0.5 * (F[0] + st[1]), 0.5 * (F[1] + st[2]), 0.5 * (F[2] + st[3]))
        dw = Fm[0] * du[0] + Fm[1] * du[1] + Fm[2] * du[2]
        f_old = math.sqrt(st[1] ** 2 + st[2] ** 2 + st[3] ** 2)
        f_new = math.sqrt(F[0] ** 2 + F[1] ** 2 + F[2] ** 2)
        if f_new > 1.01 * f_old and f_new > 1.0:            # elastic loading -> sign convention
            _sign["pos" if dw > 0 else "neg"] += abs(dw)
        st[8] += math.sqrt(du[0] ** 2 + du[1] ** 2 + du[2] ** 2)
        st[1], st[2], st[3], st[4], st[5], st[6] = F[0], F[1], F[2], u[0], u[1], u[2]
        live.add(i)
        try:
            fn = float(sc.force_norm())
        except Exception:
            fn = 0.0
        for k in (r, st[9]):
            if not k:
                continue
            _acc[k]["W"] += dw
            cnt[k] += 1
            N[k] += fn
            if st[7]:
                E[k] += f_new ** 2 / (2.0 * st[7])
            if st[8] > pmax[k]:
                pmax[k] = st[8]
    for i in [k for k in _S if k not in live]:               # subcontacts that disappeared
        st = _S.pop(i)
        if st[7]:
            f = math.sqrt(st[1] ** 2 + st[2] ** 2 + st[3] ** 2)
            for k in (st[0], st[9]):
                if k:
                    _acc[k]["gone"] += f ** 2 / (2.0 * st[7])
    s = 1.0 if _sign["pos"] > _sign["neg"] else -1.0         # no evidence yet: -1 (3DEC, confirmed by the probe)
    row = ["{:.1f}".format(time.time() - _t0), tag, "{:.4f}".format(d_ctrl_mm), "{:+.0f}".format(s)]
    for r in REGIONS:
        Wsh = s * _acc[r]["W"]
        row += ["{:.4f}".format(Wsh - (E[r] - _acc[r]["E0"])), "{:.4f}".format(Wsh), "{:.4f}".format(E[r]),
                "{:.4f}".format(N[r] / 1e3), str(cnt[r]), "{:.4f}".format(pmax[r] * 1e3)]
    _csv.write(",".join(row) + "\n"); _csv.flush()


def close():
    global _csv
    if _csv is not None:
        _csv.close(); _csv = None
    print("[slip_work_tracker] closed; sign evidence +{:.3g} / -{:.3g} J".format(_sign["pos"], _sign["neg"]))
