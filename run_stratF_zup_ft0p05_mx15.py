# -*- coding: ascii -*-
"""
run_stratF_zup_ft0p05_mx15.py -- Strategy F, runs 1-25, on the z-up pk2 f_t 0.05 MPa build with 1.5 % Maxwell damping
(2 Nov 2026). The 1.5 % counterpart of stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr: same build, truncated tables
(vel_*_tr.txt) and 4 s tail; only the Maxwell set differs (MAXWELL_SETS[1.5] = 0.0112 1.0005 0.0081 7.4315 0.0096
25.0004, the set of the f_t = 0 1.5 % cases).

Runs run_stratF_zup.py (NOT edited: read and patched in memory) with BUILD / DAMP_PCT / DRY_RUN set here.
OUT_DIR: stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_tr  (a stopped run resumes from its checkpoint)

USAGE (3DEC 9, Hanze_SAFEGO folder):  program call 'run_stratF_zup_ft0p05_mx15.py'
Post-processing:  python postprocess_stratF_full25.py --runs "pk2 f_t 0.05, 1.5 %=stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_tr" \
                      --runs "pk2 f_t 0.05, 5 %=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr"   (add --run-max N while it runs)
"""
import os, re

DRY_RUN  = False
BUILD    = "zup_pk2_ft0p05"
DAMP_PCT = 1.5
BASE     = "run_stratF_zup.py"
TAG = "# [run_stratF_zup_ft0p05_mx15]"


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


with open(BASE) as f:
    _src = f.read()
for _k, _v in (("BUILD", BUILD), ("DAMP_PCT", DAMP_PCT), ("DRY_RUN", DRY_RUN)):
    _src = set_last(_src, _k, _v)
print("[run_stratF_zup_ft0p05_mx15] {} with BUILD {}, DAMP {} %, DRY_RUN {}".format(BASE, BUILD, DAMP_PCT, DRY_RUN))
exec(compile(_src, "{} [{} {} %]".format(BASE, BUILD, DAMP_PCT), "exec"),
     {"__name__": "__main__", "__file__": os.path.abspath(BASE)})
