# -*- coding: ascii -*-
"""
run_stratF_ft0_damp3.py -- Strategy F (recorded table motions, runs 1-25) on the f_t = 0 build with 3 % damping.

The counterpart of stratF_full_LS_MAXWELL15_v12 (f_t = 0, 1.5 %) with ONE change, the damping:
    MAXWELL_CMD   block mech damp maxwell 0.0224 1.0005 0.0162 7.4315 0.0192 25.0004    3 % flat, 1-25 Hz
                  (v12: 0.0112 / 0.0081 / 0.0096 = 1.5 %; the batch's 5 %: 0.0373 / 0.0270 / 0.0320 -- the three
                  coefficients scale linearly with the target ratio, 5 / 1.5 = 3.33 exactly, so 3 % = 2 x v12)
    OUT_DIR/CASE_ID stratF_full_LS_MAXWELL30_ft0 / F_LS_MAXWELL30_ft0
Kept as in v12, so the comparison isolates the damping: START_SAVE ..._springforce.sav (f_t = 0, c = 0.2 MPa),
large strain, Table-2 protocol, T1_init 0.103 s, T_B rotation-y/z freed, no joint property touched, the ORIGINAL
vel_*.txt tables (USE_BC_TABLES = True switches to the baseline-corrected ones the f_t = 0.15 wrapper uses; that is
a second change, so leave it off for the damping comparison).

HOW IT WORKS (as run_stratF_ft0p15.py)
    Reads strategy_F_full_3dec_LS.py, rewrites the config lines IN MEMORY (last assignment of each name), checks
    every rewrite took, and executes the result. The driver on disk is not changed. Before anything runs it
    restores START_SAVE once and checks FISH 'tension' <= SAVE_TENSION_MAX (the f_t = 0 build, not a bonded save).
    The driver checkpoints after every run (stratC_checkpoint.json in OUT_DIR); calling this file again resumes.
    The checkpoint stores MAXWELL_CMD, so it cannot resume a folder written with other damping.

USAGE (3DEC 9, working folder = the one holding the saves, the driver and vel_*.txt):
    DRY_RUN = True :  program call 'run_stratF_ft0_damp3.py'   (or plain python) -> files, overrides, compile
    DRY_RUN = False:  program call 'run_stratF_ft0_damp3.py'
Post-processing:  python postprocess_stratF_full25.py --runs "f_t 0, 1.5 %=stratF_full_LS_MAXWELL15_v12" \\
                                                       --runs "f_t 0, 3 %=stratF_full_LS_MAXWELL30_ft0"
"""
import os, re

# ============================ CONFIG =================================
DRY_RUN       = True
DRIVER        = "strategy_F_full_3dec_LS.py"
SAVE_FT0      = "Part_I_MASON_v8_LargeStrain_springforce.sav"
SAVE_TENSION_MAX = 0.05e6                  # Pa, the save must hold at most this joint tension (f_t = 0 build)
MAXWELL_3PCT  = "block mech damp maxwell 0.0224 1.0005 0.0162 7.4315 0.0192 25.0004"
USE_BC_TABLES = False                      # False: v12's vel_*.txt (damping is the only change)
OVERRIDES = dict(                          # NAME = value lines rewritten in the driver (last assignment of each)
    OUT_DIR="stratF_full_LS_MAXWELL30_ft0",
    CASE_ID="F_LS_MAXWELL30_ft0",
    START_SAVE=SAVE_FT0,
    TENSION_MAX_PA=SAVE_TENSION_MAX,
    MAXWELL_CMD=MAXWELL_3PCT,
    T1_init=0.103,
    RUN_FROM=1,
    RUN_TO=25,
)
BASE_TABLES_BC = {
    "HU12": ("vel_HU", "vel_HU_bc.txt"),
    "EC40": ("vel_EC", "vel_EC_bc.txt"),
    "FR76": ("vel_FR", "vel_FR_bc.txt"),
}
KEEP_AS_V12 = ("LARGE_STRAIN", "G_I_NEW", "COH_RESIDUAL", "DAMP_RATIO", "BASE_TABLES")   # printed, not changed
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
TAG = "# [run_stratF_ft0_damp3 override]"


def log(s):
    print("[run_stratF_ft0_damp3] " + s)


def patch_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...' -- the driver changed; update this wrapper".format(DRIVER, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def patch_tables(src, tables):
    hits = list(re.finditer(r"^BASE_TABLES[ \t]*=[ \t]*\{.*?^\}", src, re.M | re.S))
    if not hits:
        raise RuntimeError("{}: BASE_TABLES block not found -- update this wrapper".format(DRIVER))
    m = hits[-1]
    body = "BASE_TABLES = {   " + TAG + "\n" + "".join("    {!r}: {!r},\n".format(k, v) for k, v in tables.items()) + "}"
    return src[:m.start()] + body + src[m.end():]


def current_value(src, name):
    hits = re.findall(r"^{}[ \t]*=[ \t]*(.*)$".format(re.escape(name)), src, re.M)
    return hits[-1].strip() if hits else "(not found)"


def main():
    log("DRY RUN" if DRY_RUN else "RUN")
    tables = [fn for _, fn in BASE_TABLES_BC.values()] if USE_BC_TABLES else ["vel_HU.txt", "vel_EC.txt", "vel_FR.txt"]
    miss = [f for f in [DRIVER, SAVE_FT0] + tables if not os.path.isfile(f)]
    if miss:
        raise RuntimeError("missing in {}: {}".format(os.getcwd(), ", ".join(miss)))
    with open(DRIVER) as f:
        src = f.read()
    for name in KEEP_AS_V12:
        if not (name == "BASE_TABLES" and USE_BC_TABLES):
            log("  kept   {} = {}".format(name, current_value(src, name)))
    log("  was    MAXWELL_CMD = {}".format(current_value(src, "MAXWELL_CMD")))
    out = src
    for name, val in OVERRIDES.items():
        out = patch_last(out, name, val)
        log("  set    {} = {!r}".format(name, val))
    if USE_BC_TABLES:
        out = patch_tables(out, BASE_TABLES_BC)
        log("  set    BASE_TABLES = {}".format({k: v[1] for k, v in BASE_TABLES_BC.items()}))
    for name in OVERRIDES:
        if TAG not in current_value(out, name):
            raise RuntimeError("override of {} did not take (a later assignment wins?)".format(name))
    if USE_BC_TABLES and TAG not in current_value(out, "BASE_TABLES"):
        raise RuntimeError("override of BASE_TABLES did not take")
    code = compile(out, DRIVER + " [ft0 damp3]", "exec")
    log("  patched driver compiles")
    if os.path.isfile(os.path.join(OVERRIDES["OUT_DIR"], "stratC_checkpoint.json")):
        log("  checkpoint found in {}: the driver resumes after its last finished run".format(OVERRIDES["OUT_DIR"]))
    if DRY_RUN:
        log("DRY_RUN: stopping before 3DEC (set DRY_RUN = False to run)")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC (or set DRY_RUN = True)")
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(SAVE_FT0.replace("\\", "/")))
    it.command("python-reset-state false")
    try:
        ft = float(it.fish.get("tension"))
    except Exception as e:
        raise RuntimeError("{}: FISH 'tension' not readable ({}) -- cannot verify the build".format(SAVE_FT0, e))
    log("  {}: joint tension {:.3f} MPa (expected <= {:.2f})".format(SAVE_FT0, ft / 1e6, SAVE_TENSION_MAX / 1e6))
    if ft > SAVE_TENSION_MAX:
        raise RuntimeError("wrong save: joint tension {:.3f} MPa, not the f_t = 0 build".format(ft / 1e6))
    ns = {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)}
    exec(code, ns)
    log("done: {}".format(OVERRIDES["OUT_DIR"]))


main()
