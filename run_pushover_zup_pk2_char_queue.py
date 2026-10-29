# -*- coding: ascii -*-
"""
run_pushover_zup_pk2_char_queue.py -- run_pushover_zup_pk2_queue.py for the CHARACTERISATION-strength pk2 build
(ANALYSIS_PART_I_MASON_zup_pk2_char.dat: tension 0.344 MPa, cohesion 0.75 MPa on every joint): pushover_oop_3dec_v5_zup.py
on the zup_pk2_char save, PUSH_DIR +1 then -1, in one 3DEC session (29 Oct 2026). The driver's bonded-build guard
TENSION_MAX_PA is raised to 0.35 MPa for this save. The driver file is not edited: its config lines are set in memory, as in
run_pushover_slipwork.py.

QUEUE: zup_pk2_char +1, then -1 (pockets at the test heights, characterisation strengths)
       -> pushover_results_v5zup_pk2_char_TopMason/pushover_pos.csv, pushover_neg.csv
A job whose CSV already exists is skipped, so a stopped queue resumes (delete a half-finished CSV first).
Both builds have f_t = 0 (tension 0, c = 0.2, Gs = 100), the same as the old-frame curve "ft0_mason":
    zup  vs ft0_mason  -> the frame only (the curves must overlay)
    pk2  vs zup        -> the pocket heights only
USAGE (3DEC 9, working folder with the saves and pushover_oop_3dec_v5_zup.py):
    MODE = "dry"  -> program call 'run_pushover_zup_pk2_char_queue.py'   (patch + compile only, seconds)
    MODE = "run"  -> program call 'run_pushover_zup_pk2_char_queue.py'
Post-processing: python postprocess_pushover_oop.py   (case char_zup_pk2 is in its CASES)
"""
import os, re, time

# ============================ CONFIG =================================
MODE   = "run"                                  # "dry" | "run"
DRIVER = "pushover_oop_3dec_v5_zup.py"
JOBS = [   # (label, save, OUT_DIR, PUSH_DIR)
    ("zup_pk2_char +", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char", "pushover_results_v5zup_pk2_char_TopMason", +1),
    ("zup_pk2_char -", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char", "pushover_results_v5zup_pk2_char_TopMason", -1),
]
TOP_JOINT = "mason"
TENSION_MAX_PA = 0.35e6                        # the save is bonded (0.344 MPa) on purpose
LOG = "run_pushover_zup_pk2_char_queue.log"
# ======================================================================
TAG = "# [run_pushover_zup_pk2_char_queue]"


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
    log("run_pushover_zup_pk2_char_queue  MODE {}".format(MODE))
    if not os.path.isfile(DRIVER):
        raise RuntimeError("{} not found in {}".format(DRIVER, os.getcwd()))
    with open(DRIVER) as f:
        src0 = f.read()
    if "Z-UP" not in src0:
        raise RuntimeError("{} is not the z-up driver".format(DRIVER))
    for label, sv, out, dr in JOBS:
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
