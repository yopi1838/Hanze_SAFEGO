# -*- coding: ascii -*-
"""
run_r2sep_TpTp_unclamp_mx15.py -- Route 2, both sine waves at the record's T_p, each wave fitted on its own, on the
configuration of the unclamped Strategy F run stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr (12 Nov 2026):
    build    Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_unclamp.sav (pk2 f_t 0.05, rigid blocks, joist top
             faces in the pockets soft and tensionless: ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_unclamp.dat)
    damping  1.5 % Maxwell (0.0112 1.0005 0.0081 7.4315 0.0096 25.0004), large strain
    T1, T_eff from the UNCLAMPED pushovers (pushover_results_v5zup_pk2_ft0p05_unclamp_TopMason/pushover_{pos,neg}.csv),
             same rule as before: d_u at 80 % of the peak on each branch, mean of the two, L = 1200 kg
    wave A : one cycle of table velocity at T_p, sized to Sd(T1)    = s x Sd_record(T1)
    wave B : one cycle at T_p,                   sized to Sd(T_eff) = s x Sd_record(T_eff)
    FIT_MODE "separate": each amplitude set alone, then joined (A, then B); the combined signal's Sd at T1 and T_eff
    is logged (SdA_over_tgtA / SdB_over_tgtB), not corrected.

Runs route2_Tm_ihsan_full25_3dec_zup.py (NOT edited: read and patched in memory). The driver has no entry for the
unclamped save, so one is added to its _BUILDS here (same joist heights as zup_pk2_ft0p05). The driver's build audit
(masonry-masonry subcontacts against the save's FISH tension / cohesion / kn / friction) still runs: the unclamp .dat
changes joist-masonry contacts only.
OUT_DIR (named by the driver): r2sep_TpTp_du80_full25_LS_MAXWELL15_zup_pk2_ft0p05_unclamp   (resumes from its checkpoint)

USAGE (Hanze_SAFEGO folder)
    check first, plain python:  PREFLIGHT_ONLY = True,  python run_r2sep_TpTp_unclamp_mx15.py -> OUT_DIR/r2tm_targets.csv
    3DEC 9:                     PREFLIGHT_ONLY = False, program call 'run_r2sep_TpTp_unclamp_mx15.py'
"""
import os, re

# ============================ CONFIG =================================
PREFLIGHT_ONLY = False
BUILD     = "zup_pk2_ft0p05_unclamp"
DAMP_PCT  = 1.5
PERIOD_A  = "Tp"
PERIOD_B  = "Tp"
FIT_MODE  = "separate"
GATE      = "flag"          # flagged runs (pga, table, leakA) run anyway: that is the comparison
BUILD_ENTRY = ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_unclamp.sav",
               "pushover_results_v5zup_pk2_ft0p05_unclamp_TopMason",
               (0.23455, 2.11091), (0.32250, 2.19886))      # save, pushover folder, joist bottoms Z, joist centres Z
# ======================================================================
DRIVER = "route2_Tm_ihsan_full25_3dec_zup.py"
TAG = "# [run_r2sep_TpTp_unclamp_mx15]"


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(DRIVER, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


with open(DRIVER) as f:
    _src = f.read()
if '"{}"'.format(BUILD) not in _src:
    _m = re.search(r"^_BUILDS[ \t]*=[ \t]*\{[^\n]*\n", _src, re.M)
    if not _m:
        raise RuntimeError("{}: '_BUILDS = {{' not found".format(DRIVER))
    _src = _src[:_m.end()] + "    {!r}: {!r},   {}\n".format(BUILD, BUILD_ENTRY, TAG) + _src[_m.end():]
for _k, _v in (("PREFLIGHT_ONLY", PREFLIGHT_ONLY), ("BUILD", BUILD), ("DAMP_PCT", DAMP_PCT), ("PERIOD_A", PERIOD_A),
               ("PERIOD_B", PERIOD_B), ("FIT_MODE", FIT_MODE), ("GATE", GATE)):
    _src = set_last(_src, _k, _v)
print("[run_r2sep_TpTp_unclamp_mx15] {}: BUILD {} ({}), {} % Maxwell, waves {} / {}, fit {}, preflight {}".format(
    DRIVER, BUILD, BUILD_ENTRY[0], DAMP_PCT, PERIOD_A, PERIOD_B, FIT_MODE, PREFLIGHT_ONLY))
exec(compile(_src, "{} [TpTp sep, {}, {} %]".format(DRIVER, BUILD, DAMP_PCT), "exec"),
     {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})
