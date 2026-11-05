# -*- coding: ascii -*-
"""
run_stratF_zup_pocket.py -- Strategy F (runs 1-25, recorded table motions, 1.5 % Maxwell by default, truncated tables + 4 s tail)
on a POCKET variant of the f_t 0.05 MPa pk2 build (3 Nov 2026).

The build: ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_<POCKET>.dat = the zup_pk2_ft0p05 save with the joist-pocket contacts given
the timber-mortar set of Arslan et al. (2024): ks 1.1e8 (pK11) or 3.6e8 (pK36) Pa/m, cohesion 0.1 MPa, friction 31 deg,
tension 0. Same case as stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_tr in every other respect, so the runs differ by the
pocket only. Needs the NEW instrument_history_new_zup.dat (3 Nov): Channel_8..11 are then joist-to-wall slip in mm,
to be scored against the test's ch8-11 (permanent slip only in runs 23-24, ~11 mm at the upper joist).

HOW: reads run_stratF_zup.py (NOT edited), adds the build "zup_pk2_ft0p05_<POCKET>" to its _BUILDS table, sets
BUILD / DAMP_PCT / DRY_RUN, and runs it (building the save from the .dat if it is missing).
OUT_DIR: stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_<POCKET>_tr

Superseded for parallel runs by run_pK11_stratF.py (pocket_jobs.py), which coordinates the save build.
USAGE (3DEC 9, Hanze_SAFEGO folder): set POCKET, then
    DRY_RUN = True   -> program call 'run_stratF_zup_pocket.py'  (or plain python): patch + compile only
    DRY_RUN = False  -> program call 'run_stratF_zup_pocket.py'
Post-processing: python postprocess_stratF_full25.py --runs "pocket <POCKET>, 1.5 %=stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_<POCKET>_tr"
"""
import os, re

# ============================ CONFIG =================================
POCKET   = "pK11"                               # "pK11" (ks 1.1e8 Pa/m) | "pK36" (ks 3.6e8 Pa/m)
DRY_RUN  = False
DAMP_PCT = 1.5                                  # 1.5 % (5 Nov 2026): closer to US-1 ring-down; was 5.0
BASE     = "run_stratF_zup.py"
NEW_KEY  = "zup_pk2_ft0p05_" + POCKET
BUILD_DAT = "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_{}.dat".format(POCKET)
NEW_ENTRY = ('    "{}": ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_{}.sav", "{}",\n'
             '                       (0.04e6, 0.06e6), (0.23455, 2.11091), (0.32250, 2.19886), 0.091),\n').format(NEW_KEY, POCKET, BUILD_DAT)
# ======================================================================
TAG = "# [run_stratF_zup_pocket]"


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def main():
    if not os.path.isfile(BASE):
        raise RuntimeError("{} not found in {}".format(BASE, os.getcwd()))
    with open(BASE) as f:
        src = f.read()
    if '"{}"'.format(NEW_KEY) not in src:
        m = re.search(r"^_BUILDS[ \t]*=[ \t]*\{[^\n]*\n", src, re.M)
        if not m:
            raise RuntimeError("{}: '_BUILDS = {{' not found -- update this runner".format(BASE))
        src = src[:m.end()] + NEW_ENTRY + src[m.end():]
    src = set_last(src, "BUILD", NEW_KEY)
    src = set_last(src, "DAMP_PCT", DAMP_PCT)
    src = set_last(src, "DRY_RUN", DRY_RUN)
    if not os.path.isfile(BUILD_DAT):
        raise RuntimeError("{} not found (needed to build the save)".format(BUILD_DAT))
    print("[run_stratF_zup_pocket] {} with BUILD {}, DAMP {} %, DRY_RUN {}".format(BASE, NEW_KEY, DAMP_PCT, DRY_RUN))
    code = compile(src, "{} [{}]".format(BASE, NEW_KEY), "exec")
    exec(code, {"__name__": "__main__", "__file__": os.path.abspath(BASE)})


main()
