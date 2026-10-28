# -*- coding: ascii -*-
"""
run_r2_TpTp.py -- Route 2 with BOTH sine waves at the record's T_p ("Tp", "Tp"), and T_eff taken from a MODEL,
not from the test envelope.

Runs route2_Tm_ihsan_full25_3dec.py with its config lines rewritten in memory (the driver file is not touched):
    PERIOD_A = PERIOD_B = "Tp"
    T_eff from TEFF_SOURCE:
      "dem"  : the DEM pushover itself. d_u = displacement where the DEM force has dropped to 80 % of its peak, read
               from PO_POS / PO_NEG (mean of the two branches), with the force as |a| x M_PUSH_KG (the applied load,
               F_SOURCE "applied", the 'DEM |a| m_tot' curve of the mechanical-model figure). The driver then uses
               T_eff = 2 pi sqrt(L d_u / F_y) on the same curves. Expected d_u ~ 57 mm (56.6 / 57.1 mm on 27 Oct).
      "mech" : the adjusted-hinge mechanical model: F0 = F0_MECH_KN, d_u = D80_MECH_MM (80 % of F0 on the rigid
               curve), T_eff = 2 pi sqrt(L d_u / F0) passed as TEFF_S_OVERRIDE (0.362 s with 10.1 kN, 28 mm).
    T1 always from the DEM pushover's elastic slope (the driver's own computation).
The feasibility gate is left on "flag": with T_p for both waves the HU12 runs WILL be flagged (pga / table / leakA);
the point of this run is to see what the DEM does with that input.

USAGE
    MODE = "preflight": plain python, no 3DEC. Writes targets, amplitudes, PGA and flags to OUT_DIR/r2tm_targets.csv
                        and the velocity tables. Check this first.
    MODE = "run"      : inside 3DEC (program call 'run_r2_TpTp.py'); needs the driver's usual files next to it.
"""
import os, re, math, csv

# ============================ CONFIG =================================
MODE         = "preflight"                     # "preflight" | "run"
DRIVER       = "route2_Tm_ihsan_full25_3dec.py"
TEFF_SOURCE  = "dem"                           # "dem" | "mech"
PO_POS       = "pushover_results_v5_ft0_masonrypush_TopMason/pushover_pos.csv"
PO_NEG       = "pushover_results_v5_ft0_masonrypush_TopMason/pushover_neg.csv"   # "" = positive only
M_PUSH_KG    = 1308.3                          # masonry mass the masonry-only push loads (F = |a| x this)
DROP         = 0.80                            # d_u at DROP x peak on the descending branch
M_EFF_KG     = 1200.0                          # L, as in the driver
F0_MECH_KN   = 10.1                            # adjusted-hinge rigid model, base pivot (mechanical_model_US1)
D80_MECH_MM  = 28.0                            #   its 80 % point at the 2.06 m sensor
OUT_DIR      = None                            # None: r2tm_TpTp_full25_LS_maxwell_teff<source>
# ======================================================================


def log(s):
    print("[run_r2_TpTp] " + s)


def backbone(path):
    """|d|, |a| x M_PUSH_KG (kN) of the converged rows."""
    d, F = [], []
    with open(path) as f:
        rows = csv.DictReader(line for line in f if not line.startswith("#"))
        for r in rows:
            try:
                if "converged" in r and r["converged"] not in ("", None) and int(float(r["converged"])) == 0:
                    continue
                d.append(abs(float(r["d_ctrl_mm"])))
                F.append(abs(float(r["a_mps2"])) * M_PUSH_KG / 1000.0)
            except (KeyError, ValueError):
                continue
    if len(d) < 5:
        raise RuntimeError("{}: fewer than 5 usable rows".format(path))
    return d, F


def d_drop(path):
    d, F = backbone(path)
    i = max(range(len(F)), key=lambda k: F[k])
    for k in range(i + 1, len(d)):
        if d[k] > d[i] and F[k] <= DROP * F[i]:
            # interpolate between k-1 and k
            d0, d1, f0, f1 = d[k - 1], d[k], F[k - 1], F[k]
            x = d0 + (DROP * F[i] - f0) * (d1 - d0) / (f1 - f0) if f1 != f0 else d1
            return F[i], d[i], x, max(d)
    raise RuntimeError("{}: the curve never drops to {:.0f} % of its peak {:.2f} kN (ends at {:.1f} mm)".format(
        path, DROP * 100, F[i], max(d)))


def patch(src, over):
    for name, val in over.items():
        m = re.search(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M)
        if m is None:
            raise RuntimeError("{}: config line '{} = ...' not found".format(DRIVER, name))
        src = src[:m.start()] + "{} = {!r}   # [run_r2_TpTp override]".format(name, val) + src[m.end():]
    return src


def main():
    out = OUT_DIR or "r2tm_TpTp_full25_LS_maxwell_teff{}".format(TEFF_SOURCE)
    over = dict(PERIOD_A="Tp", PERIOD_B="Tp", OUT_DIR=out, GATE="flag", PREFLIGHT_ONLY=(MODE == "preflight"),
                PUSHOVER_CSV=PO_POS, PUSHOVER_CSV_NEG=PO_NEG, F_SOURCE="applied", M_PUSH_KG=M_PUSH_KG,
                M_EFF_KG=M_EFF_KG)
    if TEFF_SOURCE == "dem":
        res = []
        for lab, p in (("pos", PO_POS), ("neg", PO_NEG)):
            if not p:
                continue
            Fy, dy, du, dmax = d_drop(p)
            res.append(du)
            log("DEM {}: peak {:.2f} kN at {:.1f} mm, {:.0f} % at d_u = {:.1f} mm (curve to {:.1f} mm)".format(
                lab, Fy, dy, DROP * 100, du, dmax))
        D_U = sum(res) / len(res)
        over.update(D_U_MM=round(D_U, 2), TEFF_S_OVERRIDE=None, T1_S_OVERRIDE=None)
        log("T_eff source: DEM pushover, d_u = {:.1f} mm (mean of {} branch(es))".format(D_U, len(res)))
    elif TEFF_SOURCE == "mech":
        Te = 2 * math.pi * math.sqrt(M_EFF_KG * D80_MECH_MM * 1e-3 / (F0_MECH_KN * 1e3))
        over.update(TEFF_S_OVERRIDE=round(Te, 4), T1_S_OVERRIDE=None, D_U_MM=D80_MECH_MM)
        log("T_eff source: mechanical model, F0 {:.1f} kN, d_u {:.1f} mm -> T_eff {:.4f} s".format(
            F0_MECH_KN, D80_MECH_MM, Te))
    else:
        raise RuntimeError("TEFF_SOURCE must be 'dem' or 'mech'")
    with open(DRIVER) as f:
        src = f.read()
    src = patch(src, over)
    compile(src, DRIVER, "exec")
    log("{} patched: {}".format(DRIVER, ", ".join("{}={!r}".format(k, v) for k, v in over.items())))
    log("mode {} -> {}".format(MODE, out))
    exec(compile(src, DRIVER, "exec"), {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})


main()
