# -*- coding: ascii -*-
"""
run_r2sep_TpTm_du30_zup.py -- Route 2, two waves (T_p / T_m, each fitted on its own), with the ultimate point of the
pushover fixed at d_u = 30 mm (16 Nov 2026).
    wave A : one cycle of table velocity at the record's T_p, sized to Sd(T1)    = s x Sd_record(T1)
    wave B : one cycle at the record's T_m,                   sized to Sd(T_eff) = s x Sd_record(T_eff)
    T_eff  = 2 pi sqrt(L d_u / F(d_u)), d_u = 30 mm on the mean +/- pushover, L = 1200 kg.
WHY 30 mm: the model is pushed to the experimental limit, about 30 mm (US-1 run 24). The 80 %-of-peak drop that set
T_eff before (0.4214 s) lies beyond it on the + push and on the mean curve, so the ultimate point is taken at 30 mm.
Everything else as r2sep_TpTm_du80_full25_LS_MAXWELL50_zup_pk2_ft0p05: glued pockets (zup_pk2_ft0p05), 5 % Maxwell,
large strain, T1 from the same pushover, FIT_MODE "separate", GATE "flag".

Runs route2_Tm_ihsan_full25_3dec_zup.py (needs its DU_FIXED_MM option, 16 Nov 2026), patched in memory.
OUT_DIR (named by the driver): r2sep_TpTm_du30mm_full25_LS_MAXWELL50_zup_pk2_ft0p05   (a stopped run resumes)

USAGE (Hanze_SAFEGO folder)
    check first, plain python:  PREFLIGHT_ONLY = True,  python run_r2sep_TpTm_du30_zup.py  -> OUT_DIR/r2tm_targets.csv
    3DEC 9:                     PREFLIGHT_ONLY = False, program call 'run_r2sep_TpTm_du30_zup.py'
Post-processing: postprocess_r2_full25.py --runs "TpTm d_u 30 mm=r2sep_TpTm_du30mm_full25_LS_MAXWELL50_zup_pk2_ft0p05"
"""
import os, re

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False
BUILD       = "zup_pk2_ft0p05"
DAMP_PCT    = 5.0
PERIOD_A    = "Tp"
PERIOD_B    = "Tm"
FIT_MODE    = "separate"
GATE        = "flag"
DU_FIXED_MM = 30.0
# ======================================================================
DRIVER = "route2_Tm_ihsan_full25_3dec_zup.py"
TAG = "# [run_r2sep_TpTm_du30_zup]"


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
print("[run_r2sep_TpTm_du30_zup] {}: BUILD {}, {} % Maxwell, waves {} / {}, fit {}, d_u {} mm, preflight {}".format(
    DRIVER, BUILD, DAMP_PCT, PERIOD_A, PERIOD_B, FIT_MODE, DU_FIXED_MM, PREFLIGHT_ONLY))
exec(compile(_src, "{} [TpTm du30]".format(DRIVER), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})
