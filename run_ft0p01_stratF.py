# -*- coding: ascii -*-
"""
run_ft0p01_stratF.py -- Strategy F (records, runs 1-25) on the pk2 build with joint tensile strength f_t = 0.01 MPa
(Ihsan, 9 Nov 2026), 13 Nov 2026.

BUILD  ANALYSIS_PART_I_MASON_zup_pk2_ft0p01.dat -> Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p01.sav
       (built on the fly by run_stratF_zup.py if missing). Every contact, the joist-pocket (Joist-Masonry) contacts
       included, is the MASON joint (mason_v12) with the SAME mortar set: f_t 0.01 MPa, cohesion 0.2 MPa, friction,
       kn = ks = 2.5e9 Pa/m, G_I following f_t. Only the joist-to-support contacts are Mohr (as in every pk2 build).
       I.e. the clamped (glued) pocket of the f_t 0.05 reference, one parameter changed.
RUN    run_stratF_zup.py (not edited: read and patched in memory), 1.5 % Maxwell, large strain, truncated records
       + 4 s tail, as stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_tr.
OUT    stratF_full_LS_MAXWELL15_zup_pk2_ft0p01_tr/   (resumes from its own checkpoint if stopped)

USAGE  3DEC 9, Hanze_SAFEGO folder:  program call 'run_ft0p01_stratF.py'
       check first in plain python:  DRY_RUN = True, python run_ft0p01_stratF.py
If the f_t 0.01 pushovers (run_ft0p01_unclamp_push_*.py) are building the same base save at that moment, start this
one after "save ready" appears in pocket_jobs_ft0p01_unclamp.log.
Post-processing: python postprocess_stratF_full25.py --runs "f_t 0.01, glued, 1.5 %=stratF_full_LS_MAXWELL15_zup_pk2_ft0p01_tr"
                 --out stratF_postproc_ft0p01_1p5
"""
import os, re

# ============================ CONFIG =================================
DRY_RUN  = False
BUILD    = "zup_pk2_ft0p01"
DAMP_PCT = 1.5
ENTRY    = ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p01.sav", "ANALYSIS_PART_I_MASON_zup_pk2_ft0p01.dat",
            (0.008e6, 0.012e6),            # tension check of the save (Pa)
            (0.23455, 2.11091), (0.32250, 2.19886),   # joist bottoms Z, joist centres Z
            0.091)                         # T1_init (s): reference for T_end / T_init only (the f_t 0.05 value)
# ======================================================================
BASE = "run_stratF_zup.py"
TAG = "# [run_ft0p01_stratF]"


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


with open(BASE) as f:
    _src = f.read()
if '"{}"'.format(BUILD) not in _src:
    _m = re.search(r"^_BUILDS[ \t]*=[ \t]*\{[^\n]*\n", _src, re.M)
    if not _m:
        raise RuntimeError("{}: '_BUILDS = {{' not found".format(BASE))
    _src = _src[:_m.end()] + "    {!r}: {!r},   {}\n".format(BUILD, ENTRY, TAG) + _src[_m.end():]
for _k, _v in (("BUILD", BUILD), ("DAMP_PCT", DAMP_PCT), ("DRY_RUN", DRY_RUN)):
    _src = set_last(_src, _k, _v)
print("[run_ft0p01_stratF] {}: BUILD {} ({}), {} % Maxwell, dry run {}".format(BASE, BUILD, ENTRY[0], DAMP_PCT, DRY_RUN))
exec(compile(_src, "{} [{}]".format(BASE, BUILD), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(BASE)})
