# -*- coding: ascii -*-
"""
run_r2sep_TpTm_du30_char_mx15_zup.py -- Route 2, T_p / T_m two waves, d_u = 30 mm, on the CHARACTERISATION-strength
build (f_t 0.344 MPa, cohesion 0.75 MPa on every joint, pockets at the test heights: zup_pk2_char) with 1.5 % Maxwell
damping (16 Nov 2026).
    wave A : one cycle at the record's T_p, sized to Sd(T1)    = s x Sd_record(T1)
    wave B : one cycle at the record's T_m, sized to Sd(T_eff) = s x Sd_record(T_eff)
    T1, T_eff from THIS build's own pushover (pushover_results_v5zup_pk2_char_TopMason, mean of + and - push),
    T_eff = 2 pi sqrt(L d_u / F(d_u)), d_u = 30 mm, L = 1200 kg.
Same as run_r2sep_TpTm_du30_zup.py (glued f_t 0.05, 5 %) except BUILD "zup_pk2_char" and DAMP_PCT 1.5.

NEEDS, in the Hanze_SAFEGO folder (checked below before anything runs):
    Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char.sav   <- program call 'ANALYSIS_PART_I_MASON_zup_pk2_char.dat'
    pushover_results_v5zup_pk2_char_TopMason/pushover_pos.csv, pushover_neg.csv   <- program call
                                                                   'run_pushover_zup_pk2_char_queue.py' (29 Oct queue)
    route2_Tm_ihsan_full25_3dec_zup.py with the DU_FIXED_MM option (16 Nov 2026)
OUT_DIR (named by the driver): r2sep_TpTm_du30mm_full25_LS_MAXWELL15_zup_pk2_char   (a stopped run resumes)

USAGE (Hanze_SAFEGO folder)
    check first, plain python:  PREFLIGHT_ONLY = True,  python run_r2sep_TpTm_du30_char_mx15_zup.py
                                -> T1, T_eff of the char build and OUT_DIR/r2tm_targets.csv
    3DEC 9:                     PREFLIGHT_ONLY = False, program call 'run_r2sep_TpTm_du30_char_mx15_zup.py'
Post-processing: postprocess_r2_full25.py --runs "TpTm char 1.5 %=r2sep_TpTm_du30mm_full25_LS_MAXWELL15_zup_pk2_char"
"""
import os, re

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False
BUILD       = "zup_pk2_char"
DAMP_PCT    = 1.5
PERIOD_A    = "Tp"
PERIOD_B    = "Tm"
FIT_MODE    = "separate"
GATE        = "flag"
DU_FIXED_MM = 30.0
# ======================================================================
DRIVER = "route2_Tm_ihsan_full25_3dec_zup.py"
TAG = "# [run_r2sep_TpTm_du30_char_mx15_zup]"
NEED = [("pushover_results_v5zup_pk2_char_TopMason/pushover_pos.csv", "program call 'run_pushover_zup_pk2_char_queue.py'"),
        ("pushover_results_v5zup_pk2_char_TopMason/pushover_neg.csv", "program call 'run_pushover_zup_pk2_char_queue.py'")]
if not PREFLIGHT_ONLY:
    NEED.append(("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char.sav",
                 "program call 'ANALYSIS_PART_I_MASON_zup_pk2_char.dat'"))
_miss = [(p, how) for p, how in NEED if not os.path.isfile(p)]
if _miss:
    raise RuntimeError("[run_r2sep_TpTm_du30_char_mx15_zup] missing:\n" +
                       "\n".join("    {}   -> first: {}".format(p, how) for p, how in _miss))


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...' (old driver? it needs DU_FIXED_MM)".format(DRIVER, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


with open(DRIVER) as f:
    _src = f.read()
for _k, _v in (("PREFLIGHT_ONLY", PREFLIGHT_ONLY), ("BUILD", BUILD), ("DAMP_PCT", DAMP_PCT), ("PERIOD_A", PERIOD_A),
               ("PERIOD_B", PERIOD_B), ("FIT_MODE", FIT_MODE), ("GATE", GATE), ("DU_FIXED_MM", DU_FIXED_MM)):
    _src = set_last(_src, _k, _v)
print("[run_r2sep_TpTm_du30_char_mx15_zup] {}: BUILD {}, {} % Maxwell, waves {} / {}, fit {}, d_u {} mm, preflight {}".format(
    DRIVER, BUILD, DAMP_PCT, PERIOD_A, PERIOD_B, FIT_MODE, DU_FIXED_MM, PREFLIGHT_ONLY))
exec(compile(_src, "{} [TpTm du30 char 1.5]".format(DRIVER), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})
