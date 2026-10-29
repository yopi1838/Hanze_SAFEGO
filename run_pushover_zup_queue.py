# -*- coding: ascii -*-
"""
run_pushover_zup_queue.py -- runs pushover_oop_3dec_v5_zup.py for several z-up saves and both directions in one
3DEC session (29 Oct 2026). The driver file is not edited: its config lines are set in memory, as in
run_pushover_slipwork.py.

QUEUE (in this order; a job whose pushover_<pos|neg>.csv already exists is skipped, so a stopped queue resumes):
    1-2  zup_pk2 : pockets at the test heights (Geo_Prep_zup_pk2.dat)  -> pushover_results_v5zup_pk2_ft0_TopMason
    3-4  zup     : old pocket heights, frame check at large displacement -> pushover_results_v5zup_ft0_TopMason
Both builds have f_t = 0 (tension 0, c = 0.2, Gs = 100), the same as the old-frame curve "ft0_mason":
    zup  vs ft0_mason  -> the frame only (the curves must overlay)
    pk2  vs zup        -> the pocket heights only
USAGE (3DEC 9, working folder with the saves and pushover_oop_3dec_v5_zup.py):
    MODE = "dry"  -> program call 'run_pushover_zup_queue.py'   (patch + compile only, seconds)
    MODE = "run"  -> program call 'run_pushover_zup_queue.py'
Post-processing: python postprocess_pushover_oop.py   (cases ft0_zup and ft0_zup_pk2 are in its CASES)
"""
import os, re, time

# ============================ CONFIG =================================
MODE   = "run"                                  # "dry" | "run"
DRIVER = "pushover_oop_3dec_v5_zup.py"
JOBS = [   # (label, save, OUT_DIR, PUSH_DIR)
    ("zup_pk2 +", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2", "pushover_results_v5zup_pk2_ft0_TopMason", +1),
    ("zup_pk2 -", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2", "pushover_results_v5zup_pk2_ft0_TopMason", -1),
    ("zup +",     "Part_I_MASON_v8_LargeStrain_springforce_zup",     "pushover_results_v5zup_ft0_TopMason",     +1),
    ("zup -",     "Part_I_MASON_v8_LargeStrain_springforce_zup",     "pushover_results_v5zup_ft0_TopMason",     -1),
]
TOP_JOINT = "mason"
LOG = "run_pushover_zup_queue.log"
# ======================================================================
TAG = "# [run_pushover_zup_queue]"


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
    log("run_pushover_zup_queue  MODE {}".format(MODE))
    if not os.path.isfile(DRIVER):
        raise RuntimeError("{} not found in {}".format(DRIVER, os.getcwd()))
    with open(DRIVER) as f:
        src0 = f.read()
    if "Z-UP" not in src0:
        raise RuntimeError("{} is not the z-up driver".format(DRIVER))
    for label, sv, out, dr in JOBS:
        src = src0
        for k, v in dict(RESTORE_SAV=sv, OUT_DIR=out, PUSH_DIR=dr, TOP_JOINT=TOP_JOINT).items():
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
