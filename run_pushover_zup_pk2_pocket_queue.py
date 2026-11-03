# -*- coding: ascii -*-
"""
run_pushover_zup_pk2_pocket_queue.py -- +/- pushovers of the two POCKET variants of the f_t 0.05 MPa pk2 build
(3 Nov 2026): ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11.dat / _pK36.dat = the zup_pk2_ft0p05 save with the joist-pocket
contacts given the timber-mortar set of Arslan et al. (2024): ks 1.1e8 / 3.6e8 Pa/m, cohesion 0.1 MPa, friction 31 deg,
tension 0. Compare with pushover_results_v5zup_pk2_ft0p05_TopMason (same build, brick-mortar pockets).
pushover_oop_3dec_v5_zup.py is not edited: its config lines are set in memory, as in the earlier queues.

QUEUE: pK11 +1, pK11 -1, pK36 +1, pK36 -1
       -> pushover_results_v5zup_pk2_ft0p05_pK11_TopMason/pushover_pos.csv, pushover_neg.csv (and _pK36_)
A missing save is built first (program call its .dat; needs Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05.sav).
A job whose CSV already exists is skipped, so a stopped queue resumes (delete a half-finished CSV first).
Check the build log for the POCKET lines: AFTER must read stiffness-shear 1.1e8 (or 3.6e8) and cohesion 1e5.
USAGE (3DEC 9, Hanze_SAFEGO folder):
    MODE = "dry"  -> program call 'run_pushover_zup_pk2_pocket_queue.py'   (patch + compile only, seconds)
    MODE = "run"  -> program call 'run_pushover_zup_pk2_pocket_queue.py'
Post-processing: python postprocess_pushover_oop.py   (add the two folders to its CASES)
"""
import os, re, time

# ============================ CONFIG =================================
MODE   = "run"                                  # "dry" | "run"
DRIVER = "pushover_oop_3dec_v5_zup.py"
JOBS = [   # (label, save, OUT_DIR, PUSH_DIR, build .dat)
    ("pK11 +", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_pK11", "pushover_results_v5zup_pk2_ft0p05_pK11_TopMason", +1,
     "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11.dat"),
    ("pK11 -", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_pK11", "pushover_results_v5zup_pk2_ft0p05_pK11_TopMason", -1,
     "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11.dat"),
    ("pK36 +", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_pK36", "pushover_results_v5zup_pk2_ft0p05_pK36_TopMason", +1,
     "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK36.dat"),
    ("pK36 -", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_pK36", "pushover_results_v5zup_pk2_ft0p05_pK36_TopMason", -1,
     "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK36.dat"),
]
TOP_JOINT = "mason"
TENSION_MAX_PA = 0.06e6                        # the save has f_t = 0.05 MPa on purpose
LOG = "run_pushover_zup_pk2_pocket_queue.log"
# ======================================================================
TAG = "# [run_pushover_zup_pk2_pocket_queue]"


def log(s=""):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), s)
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def patch_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("driver has no config line '{} = ...'".format(name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def main():
    log("=" * 70)
    log("run_pushover_zup_pk2_pocket_queue  MODE {}".format(MODE))
    if not os.path.isfile(DRIVER):
        raise RuntimeError("{} not found in {}".format(DRIVER, os.getcwd()))
    with open(DRIVER) as f:
        src0 = f.read()
    if "Z-UP" not in src0:
        raise RuntimeError("{} is not the z-up driver".format(DRIVER))
    for _l, sv0, _o, _d, dat in JOBS:
        if not os.path.isfile(dat):
            raise RuntimeError("{} not found".format(dat))
        if MODE != "dry" and not (os.path.isfile(sv0) or os.path.isfile(sv0 + ".sav")):
            import itasca as it
            log("save {} missing: program call '{}'".format(sv0, dat))
            it.command("python-reset-state false")
            it.command("program call '{}'".format(dat))
            it.command("python-reset-state false")
    for label, sv, out, dr, _dat in JOBS:
        src = src0
        for k, v in dict(RESTORE_SAV=sv, OUT_DIR=out, PUSH_DIR=dr, TOP_JOINT=TOP_JOINT, TENSION_MAX_PA=TENSION_MAX_PA).items():
            src = patch_last(src, k, v)
        code = compile(src, "{} [{}]".format(DRIVER, label), "exec")
        csvp = os.path.join(out, "pushover_{}.csv".format("pos" if dr > 0 else "neg"))
        log("{}: {} -> {}".format(label, sv, csvp))
        if MODE == "dry":
            continue
        if not (os.path.isfile(sv) or os.path.isfile(sv + ".sav")):
            log("  ! save {} missing -- skipped".format(sv)); continue
        if os.path.isfile(csvp):
            log("  exists -- skipped (delete it to re-run)"); continue
        t0 = time.time()
        try:
            exec(code, {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)})
            log("  done in {:.1f} min".format((time.time() - t0) / 60.0))
        except Exception as e:      # one failed job does not stop the queue
            log("  ! FAILED after {:.1f} min: {!r}".format((time.time() - t0) / 60.0, e))
    log("finished" if MODE != "dry" else "dry run: all jobs patch and compile; set MODE = 'run'")


main()
