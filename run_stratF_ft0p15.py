# -*- coding: ascii -*-
"""
run_stratF_ft0p15.py -- Strategy F (recorded table motions, runs 1-25) on the f_t = 0.15 MPa build.

The counterpart of stratF_full_LS_MAXWELL15_v12 (f_t = 0) with ONE physical change, the joint bond:
    START_SAVE      Part_I_MASON_v8_LargeStrain_springforce_ft_0p15.sav   (v12: ..._springforce.sav, f_t = 0)
    TENSION_MAX_PA  0.2e6        (the driver's guard against a bonded save; v12: 0.05e6)
    BASE_TABLES     vel_*_bc.txt (baseline-corrected: same motion, net table displacement 0 instead of
                    0.19 / 18.4 / 17.8 mm per 100 % record -- removes the table drift, not the wall response)
    T1_init         0.090 s      (period-ID windows and the T_end / T_init column only; the bonded pushover gives
                    T1 0.089 s; v12: 0.103 s. Never used for loading.)
    OUT_DIR/CASE_ID stratF_full_LS_MAXWELL15_ft0p15_bc / F_LS_MAXWELL15_ft0p15_bc
Kept as in v12, so the comparison isolates the bond: Maxwell 1.5 % (MAXWELL_CMD), large strain, Table-2 protocol,
T_B rotation-y/z freed (as v12 did), no joint property touched (G_I_NEW = None, COH_RESIDUAL = 0).

HOW IT WORKS
    Reads strategy_F_full_3dec_LS.py, rewrites the config lines above IN MEMORY (the last assignment of each name,
    because OUT_DIR is assigned twice in the driver; BASE_TABLES as a whole block), checks every rewrite took, and
    executes the result. strategy_F_full_3dec_LS.py on disk is not changed.
    Before anything runs it restores START_SAVE once and checks the joint tension the build stored (FISH 'tension'):
    it must lie in SAVE_TENSION_RANGE, so the f_t = 0 save cannot be used by mistake.
    The driver checkpoints after every run (stratC_checkpoint.json in OUT_DIR); calling this file again resumes.

USAGE (3DEC 9, working folder = the one holding the saves and the driver):
    DRY_RUN = True :  program call 'run_stratF_ft0p15.py'   (or plain python) -> files, overrides, compile; no analysis
    DRY_RUN = False:  program call 'run_stratF_ft0p15.py'
Post-processing:  python postprocess_stratF_full25.py --runs "f_t 0, 1.5 %=stratF_full_LS_MAXWELL15_v12" \\
                                                       --runs "f_t 0.15, 1.5 %, bc=stratF_full_LS_MAXWELL15_ft0p15_bc"
"""
import os, re, sys

# ============================ CONFIG =================================
DRY_RUN       = True
DRIVER        = "strategy_F_full_3dec_LS.py"
SAVE_FT015    = "Part_I_MASON_v8_LargeStrain_springforce_ft_0p15.sav"
SAVE_TENSION_RANGE = (0.10e6, 0.20e6)      # Pa, the save must hold this joint tension
OVERRIDES = dict(                          # NAME = value lines rewritten in the driver (last assignment of each)
    OUT_DIR="stratF_full_LS_MAXWELL15_ft0p15_bc",
    CASE_ID="F_LS_MAXWELL15_ft0p15_bc",
    START_SAVE=SAVE_FT015,
    TENSION_MAX_PA=0.2e6,
    T1_init=0.090,
    RUN_FROM=1,
    RUN_TO=25,
)
BASE_TABLES = {                            # replaces the driver's BASE_TABLES block (table names unchanged)
    "HU12": ("vel_HU", "vel_HU_bc.txt"),
    "EC40": ("vel_EC", "vel_EC_bc.txt"),
    "FR76": ("vel_FR", "vel_FR_bc.txt"),
}
KEEP_AS_V12 = ("MAXWELL_CMD", "LARGE_STRAIN", "G_I_NEW", "COH_RESIDUAL")   # printed, not changed
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None


def log(s):
    print("[run_stratF_ft0p15] " + s)


def patch_last(src, name, val):
    pat = re.compile(r"^{}[ \t]*=.*$".format(re.escape(name)), re.M)
    hits = list(pat.finditer(src))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...' -- the driver changed; update this wrapper".format(DRIVER, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   # [run_stratF_ft0p15 override]".format(name, val) + src[m.end():]


def patch_tables(src, tables):
    pat = re.compile(r"^BASE_TABLES[ \t]*=[ \t]*\{.*?^\}", re.M | re.S)
    hits = list(pat.finditer(src))
    if not hits:
        raise RuntimeError("{}: BASE_TABLES block not found -- update this wrapper".format(DRIVER))
    m = hits[-1]
    body = "BASE_TABLES = {   # [run_stratF_ft0p15 override]\n" + "".join(
        "    {!r}: {!r},\n".format(k, v) for k, v in tables.items()) + "}"
    return src[:m.start()] + body + src[m.end():]


def current_value(src, name):
    hits = re.findall(r"^{}[ \t]*=[ \t]*(.*)$".format(re.escape(name)), src, re.M)
    return hits[-1].strip() if hits else "(not found)"


def main():
    log("{}".format("DRY RUN" if DRY_RUN else "RUN"))
    miss = [f for f in [DRIVER, SAVE_FT015] + [fn for _, fn in BASE_TABLES.values()] if not os.path.isfile(f)]
    if miss:
        raise RuntimeError("missing in {}: {}".format(os.getcwd(), ", ".join(miss)))
    with open(DRIVER) as f:
        src = f.read()
    for name in KEEP_AS_V12:
        log("  kept   {} = {}".format(name, current_value(src, name)))
    out = src
    for name, val in OVERRIDES.items():
        out = patch_last(out, name, val)
        log("  set    {} = {!r}".format(name, val))
    out = patch_tables(out, BASE_TABLES)
    log("  set    BASE_TABLES = {}".format({k: v[1] for k, v in BASE_TABLES.items()}))
    for name in OVERRIDES:
        if "# [run_stratF_ft0p15 override]" not in current_value(out, name):
            raise RuntimeError("override of {} did not take (a later assignment wins?)".format(name))
    if "[run_stratF_ft0p15 override]" not in current_value(out, "BASE_TABLES"):
        raise RuntimeError("override of BASE_TABLES did not take")
    code = compile(out, DRIVER + " [ft0p15]", "exec")
    log("  patched driver compiles")
    if os.path.isfile(os.path.join(OVERRIDES["OUT_DIR"], "stratC_checkpoint.json")):
        log("  checkpoint found in {}: the driver resumes after its last finished run".format(OVERRIDES["OUT_DIR"]))
    if DRY_RUN:
        log("DRY_RUN: stopping before 3DEC (set DRY_RUN = False to run)")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC (or set DRY_RUN = True)")
    # build check: the save must be the bonded one
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(SAVE_FT015.replace("\\", "/")))
    it.command("python-reset-state false")
    try:
        ft = float(it.fish.get("tension"))
    except Exception as e:
        raise RuntimeError("{}: FISH 'tension' not readable ({}) -- cannot verify the build".format(SAVE_FT015, e))
    log("  {}: joint tension {:.3f} MPa (expected {:.2f}-{:.2f})".format(
        SAVE_FT015, ft / 1e6, SAVE_TENSION_RANGE[0] / 1e6, SAVE_TENSION_RANGE[1] / 1e6))
    if not (SAVE_TENSION_RANGE[0] <= ft <= SAVE_TENSION_RANGE[1]):
        raise RuntimeError("wrong save: joint tension {:.3f} MPa".format(ft / 1e6))
    ns = {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)}
    exec(code, ns)
    log("done: {}".format(OVERRIDES["OUT_DIR"]))


main()
