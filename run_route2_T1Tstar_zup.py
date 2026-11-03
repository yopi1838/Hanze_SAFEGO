# -*- coding: ascii -*-
"""
run_route2_T1Tstar_zup.py -- Route 2 with the fixed-point period T* as the second target, 3DEC runs 1-25 (3 Nov 2026).

Runs route2_Tm_ihsan_full25_3dec_zup.py (NOT edited: read and patched in memory, like the Strategy F runners) with:
    wave A : one sine cycle of table velocity at T1 (the driver's T1 from the pushover, L = 1200 kg)
    wave B : one sine cycle at T* of the run
    targets: Sd(T1) = s Sd_rec(T1) and Sd(T*) = s Sd_rec(T*), xi 5 %, Sd of the WHOLE signal (driver's integrator, 2 s tail)
    fit    : joint (Newton on log-amplitudes); one wave alone covering the other's target -> that wave only (stage A or B)
    T*     : fixed point d* = s Sd_rec(T_sec(d*), 5 %) on the mean +/- pushover backbone, L = 1200 kg
             (fixed_point_check_zup.py, column *_mean_L1200_xi05; values embedded below so the workstation needs no CSV)
    rule 1 : T* < SINGLE_RATIO x T1  ->  one target, Sd(T1), wave A only            (runs 1-14, 22)
    rule 2 : d* > the driver's d_u   ->  T* = the driver's T_eff, flag 'beyond_du'  (run 25)
    band   : Sd(signal) / Sd(s x record) over [T1, 1.5 T*] is LOGGED (band_min, band_max), not used as a gate
Same build, damping (5 % Maxwell), large strain, tail, instrumentation and outputs as r2sep_TpTm_du80_full25_LS_MAXWELL50_zup_pk2_ft0p05.

USAGE (Hanze_SAFEGO folder)
    check first, plain python (no 3DEC):   set PREFLIGHT_ONLY = True,  python run_route2_T1Tstar_zup.py
                                           -> OUT_DIR/r2tm_targets.csv and the 25 velocity tables
    3DEC 9:                                set PREFLIGHT_ONLY = False, program call 'run_route2_T1Tstar_zup.py'
    a stopped run resumes from OUT_DIR/r2tm_checkpoint.json
Post-processing:  python postprocess_r2tm_full25.py with this OUT_DIR (as for the r2sep_ case)
"""
import os, re

PREFLIGHT_ONLY = False
BASE    = "route2_Tm_ihsan_full25_3dec_zup.py"
OUT_DIR = "r2_T1Tstar_xi05_full25_LS_MAXWELL50_zup_pk2_ft0p05"
SINGLE_RATIO = 1.10
BAND_HI, N_BAND = 1.5, 40
# run: (d* mm, T* s) at xi_eq 5 %, mean backbone, L = 1200 kg  (fixed_point_check_zup/fixed_point_check.csv, 3 Nov 2026)
TSTAR = {1: (0.2374, 0.0892), 2: (0.3561, 0.0892), 3: (0.2306, 0.0892), 4: (0.4748, 0.0892), 5: (0.5935, 0.0892),
         6: (0.3459, 0.0892), 7: (0.7122, 0.0892), 8: (0.4612, 0.0892), 9: (0.8309, 0.0892), 10: (0.9496, 0.0892),
         11: (0.5764, 0.0892), 12: (1.0684, 0.0892), 13: (1.2181, 0.0899), 14: (1.5392, 0.0935), 15: (2.4804, 0.1056),
         16: (4.5673, 0.1303), 17: (6.8441, 0.1554), 18: (7.8217, 0.1696), 19: (8.4201, 0.1761), 20: (9.0482, 0.1830),
         21: (9.8882, 0.1920), 22: (1.2023, 0.0898), 23: (2.5841, 0.1069), 24: (3.1264, 0.1131), 25: (63.835, 0.6120)}
TAG = "# [run_route2_T1Tstar_zup]"

NEW_BUILD = '''
# ---- injected by run_route2_T1Tstar_zup.py ------------------------------------------------------------------
_TSTAR = {tstar!r}
_SINGLE_RATIO, _BAND_HI, _N_BAND = {single!r}, {bandhi!r}, {nband!r}
def tstar_of(n):
    ds, Ts = _TSTAR[n]
    if ds > DU:
        return float(TEFF_S), "beyond_du"
    return float(Ts), ""
def build_signal_tstar(n, rec, scale):
    Tst, tflag = tstar_of(n)
    single = Tst < _SINGLE_RATIO * T1_S
    Ts = [T1_S, Tst]
    tA = scale * sd_record_mm(rec, T1_S); tB = scale * sd_record_mm(rec, Tst)
    uA, TgA = cycles(1.0, T1_S, CYCLES_A)
    uB, TgB = cycles(PHASE_B, Tst, CYCLES_B)
    sA = sd_signal(uA, Ts) * 1e3; sB = sd_signal(uB, Ts) * 1e3
    VA0 = tA / sA[0]; VB0 = tB / sB[1]
    it_n, resid = 0, 0.0
    if single:
        VA, VB, stages = VA0, 0.0, "A1"            # one target (rule 1)
    elif VA0 * sA[1] >= tB * (1 - 1e-6):
        VA, VB, stages = VA0, 0.0, "A"
    elif VB0 * sB[0] >= tA * (1 - 1e-6):
        VA, VB, stages = 0.0, VB0, "B"
    else:
        stages = "A+B"; x = np.log([VA0, VB0])
        def resid_of(x):
            s = sd_signal(join(np.exp(x[0]) * uA, np.exp(x[1]) * uB), Ts) * 1e3
            return np.array([s[0] / tA - 1.0, s[1] / tB - 1.0])
        r = resid_of(x)
        for it_n in range(1, FIT_MAXIT + 1):
            if np.max(np.abs(r)) < FIT_TOL:
                break
            J = np.zeros((2, 2)); h = 1e-4
            for k in range(2):
                xp = x.copy(); xp[k] += h; J[:, k] = (resid_of(xp) - r) / h
            try:
                dx = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                dx = -r
            dx = np.clip(dx, -1.0, 1.0); lam = 1.0
            for _ in range(8):
                xn = x + lam * dx; rn = resid_of(xn)
                if np.max(np.abs(rn)) < np.max(np.abs(r)):
                    break
                lam *= 0.5
            x, r = xn, rn
        resid = float(np.max(np.abs(r)))
        if resid > 1e-3:
            print("  ** joint fit did not converge (resid {{:.2e}}); amplitudes are the last iterate".format(resid))
        VA, VB = float(np.exp(x[0])), float(np.exp(x[1]))
    vA = VA * uA if VA > 0 else np.zeros(0); vB = VB * uB if VB > 0 else np.zeros(0)
    v = join(vA, vB) if len(vA) else vB
    t = np.arange(len(v)) * DELTA_T
    tot = sd_signal(v, Ts) * 1e3
    Tb = np.geomspace(T1_S, _BAND_HI * (T1_S if single else Tst), _N_BAND)
    rb = sd_signal(v, Tb) * 1e3 / np.array([scale * sd_record_mm(rec, T) for T in Tb])
    pga = pga_g(v); R = REC[rec]
    flags = [tflag] if tflag else []
    if pga > PGA_RATIO_MAX * scale * R["pga"]: flags.append("pga")
    if pga > PGA_WARN_G: flags.append("table")
    if tot[0] / tA > LEAK_MAX: flags.append("leakA")
    if not single and tot[1] / tB > LEAK_MAX: flags.append("leakB")
    return dict(t=t, v=v, tA=tA, tB=tB, stages=stages, TgA=TgA, VA=VA, pgaA=pga_g(vA), TgB=TgB, VB=VB, pgaB=pga_g(vB),
                pga=pga, pgv=float(np.max(np.abs(v))), SdA_tot=float(tot[0]), SdB_tot=float(tot[1]),
                flags="|".join(flags) if flags else "ok", fit_iter=it_n, fit_resid=resid, table_s=float(t[-1]),
                Tstar=Tst, band_min=float(rb.min()), band_max=float(rb.max()))
print("  T1/T* form: wave A at T1 {{:.4f}} s, wave B at T* of the run (fixed point, xi_eq 5 %); one target when "
      "T* < {{:.2f}} T1; T* = T_eff {{:.4f}} s when d* > d_u {{:.1f}} mm".format(T1_S, _SINGLE_RATIO, TEFF_S, DU))
# ----------------------------------------------------------------------------------------------------------------
'''


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def rep1(src, old, new):
    k = src.count(old)
    if k != 1:
        raise RuntimeError("{}: expected {!r} once, found {} -- update this runner".format(BASE, old, k))
    return src.replace(old, new)


with open(BASE) as f:
    src = f.read().replace("\r\n", "\n")
for k, v in (("PREFLIGHT_ONLY", PREFLIGHT_ONLY), ("OUT_DIR", OUT_DIR), ("PERIOD_A", "T1"), ("PERIOD_B", "Tstar"),
             ("FIT_MODE", "joint")):
    src = set_last(src, k, v)
src = rep1(src, '"T1": T1_S, "Teff": TEFF_S}[spec]', '"T1": T1_S, "Teff": TEFF_S, "Tstar": float("nan")}[spec]')
src = rep1(src, '"fit_iter", "fit_resid", "table_s"]', '"fit_iter", "fit_resid", "table_s", "Tstar_s", "band_min", "band_max"]')
src = rep1(src, 'table_s=round(S["table_s"], 4))',
           'table_s=round(S["table_s"], 4), Tstar_s=round(S["Tstar"], 4), band_min=round(S["band_min"], 3), '
           'band_max=round(S["band_max"], 3))')
src = rep1(src, "\ndef target_row(", NEW_BUILD.format(tstar=TSTAR, single=SINGLE_RATIO, bandhi=BAND_HI, nband=N_BAND)
           + "\ndef target_row(")
src = rep1(src, "S = build_signal(rec, sc); SIG[n] = S", "S = build_signal_tstar(n, rec, sc); SIG[n] = S")
src = rep1(src, 'Sd(T_eff {:.3f}) {:.2f} mm   flags {}".format(\n          run_no, rec, sc, T1_S, S["tA"], TEFF_S, S["tB"]',
           'Sd(T* {:.3f}) {:.2f} mm   flags {}".format(\n          run_no, rec, sc, T1_S, S["tA"], S["Tstar"], S["tB"]')
print("[run_route2_T1Tstar_zup] {} -> {}  (PREFLIGHT_ONLY {})".format(BASE, OUT_DIR, PREFLIGHT_ONLY))
exec(compile(src, "{} [T1/T*]".format(BASE), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(BASE)})
