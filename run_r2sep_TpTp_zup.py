# -*- coding: ascii -*-
"""
run_r2sep_TpTp_zup.py -- Route 2, both sine waves at the record's T_p, each wave fitted on its own (10 Nov 2026).
    wave A : one cycle of table velocity at T_p, sized to Sd(T1)    = s x Sd_record(T1)
    wave B : one cycle at T_p,                   sized to Sd(T_eff) = s x Sd_record(T_eff)
    FIT_MODE "separate": each amplitude set alone, then joined (A, then B); the combined signal's Sd at T1 and T_eff
    is logged (SdA_over_tgtA / SdB_over_tgtB), not corrected. Only the period of wave B differs from the T_p / T_m run.

Same build, damping, large strain, T1 / T_eff, tail and instrumentation as r2sep_TpTm_du80_full25_LS_MAXWELL50_zup_pk2_ft0p05,
so the two compare run by run: zup_pk2_ft0p05 (glued pockets), 5 % Maxwell, T1 0.0883 s, T_eff 0.4214 s (d_u at 80 %
of the peak on the mean +/- pushover, L = 1200 kg).

Runs route2_Tm_ihsan_full25_3dec_zup.py (NOT edited: read and patched in memory).
OUT_DIR (named by the driver): r2sep_TpTp_du80_full25_LS_MAXWELL50_zup_pk2_ft0p05   (a stopped run resumes from its checkpoint)

USAGE (Hanze_SAFEGO folder)
    check first, plain python:  PREFLIGHT_ONLY = True,  python run_r2sep_TpTp_zup.py  -> OUT_DIR/r2tm_targets.csv
    3DEC 9:                     PREFLIGHT_ONLY = False, program call 'run_r2sep_TpTp_zup.py'
Post-processing as for the T_p / T_m run (postprocess_r2tm_full25.py with this OUT_DIR).
"""
import os, re

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False
BUILD     = "zup_pk2_ft0p05"
DAMP_PCT  = 5.0
PERIOD_A  = "Tp"
PERIOD_B  = "Tp"
FIT_MODE  = "separate"
GATE      = "flag"          # the HU12 runs ARE flagged (pga, table, leakA): they run anyway, that is the comparison
# ======================================================================
DRIVER = "route2_Tm_ihsan_full25_3dec_zup.py"
TAG = "# [run_r2sep_TpTp_zup]"


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(DRIVER, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


with open(DRIVER) as f:
    _src = f.read()
for _k, _v in (("PREFLIGHT_ONLY", PREFLIGHT_ONLY), ("BUILD", BUILD), ("DAMP_PCT", DAMP_PCT), ("PERIOD_A", PERIOD_A),
               ("PERIOD_B", PERIOD_B), ("FIT_MODE", FIT_MODE), ("GATE", GATE)):
    _src = set_last(_src, _k, _v)
print("[run_r2sep_TpTp_zup] {}: BUILD {}, {} % Maxwell, waves {} / {}, fit {}, preflight {}".format(
    DRIVER, BUILD, DAMP_PCT, PERIOD_A, PERIOD_B, FIT_MODE, PREFLIGHT_ONLY))
exec(compile(_src, "{} [TpTp sep]".format(DRIVER), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})
