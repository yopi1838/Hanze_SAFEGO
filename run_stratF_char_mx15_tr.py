# -*- coding: ascii -*-
"""
run_stratF_char_mx15_tr_tr.py -- Strategy F on the bonded characterisation build (f_t = 0.344 MPa, c = 0.75 MPa),
large strain, 1.5 % Maxwell (as run_stratF_char_mx15_tr.py), with TRUNCATED record tables and a longer FREE-VIBRATION
tail after every run.

    tables     vel_HU_tr.txt / vel_EC_tr.txt / vel_FR_tr.txt  (make_truncated_tables.py: quiet lead-in cut, 0.5 s
               fade-in, 0.25 s fade-out; 1.5 %-damped Sd at T1 and T_eff within 0.5 % of the full records).
               Written first if missing (MAKE_TABLES). HU12 35.9 -> 3.3 s, EC40 39.6 -> 26.1 s, FR76 17.3 -> 9.3 s.
    tail       TAIL_SEC_RECORD = TAIL_S (4 s; driver default 2.5 s): after each record the table velocity is removed
               and the wall rings down freely -- the part to check for over-damping against the tests.
    OUT_DIR    stratF_full_LS_MAXWELL15_char_tr
Everything else as run_stratF_char_mx15_tr.py (START_SAVE ..._springforce_char.sav, T1_init 0.089 s, Table-2 protocol).
The post-processor reads the table file of each run from strategy_C_log.csv, so its input-velocity panel uses the
truncated tables automatically.

USAGE (3DEC 9, working folder = the one holding the saves, the driver and vel_*.txt):
    DRY_RUN = True :  program call 'run_stratF_char_mx15_tr_tr.py'   (or plain python) -> tables, overrides, compile
    DRY_RUN = False:  program call 'run_stratF_char_mx15_tr_tr.py'
Post-processing:  python postprocess_stratF_full25.py --runs "f_t 0, 1.5 %=stratF_full_LS_MAXWELL15_v12" \\
                      --runs "0.344 MPa, 1.5 %, truncated=stratF_full_LS_MAXWELL15_char_tr"
"""
import os, re

# ============================ CONFIG =================================
DRY_RUN       = True
DRIVER        = "strategy_F_full_3dec_LS.py"
SAVE_FT0      = "Part_I_MASON_v8_LargeStrain_springforce_char.sav"
BUILD_DAT     = "ANALYSIS_PART_I_MASON_char.dat"   # run first if SAVE_FT0 is missing (None: never)
SAVE_TENSION_RANGE = (0.30e6, 0.40e6)     # Pa, FISH tension of the save (the characterisation build)
TAIL_S        = 4.0                        # s, free vibration after each record (driver default 2.5)
MAKE_TABLES   = True                       # write the truncated tables first if they are missing
RAY_EFF       = 0.0   # no Rayleigh in this case
MAXWELL_15    = "block mech damp maxwell 0.0112 1.0005 0.0081 7.4315 0.0096 25.0004"   # v12, 1.5 % flat 1-25 Hz
RAY_RATIO     = 2.0 * RAY_EFF              # (Rayleigh mass term, unused here: RAY_EFF = 0)
T1_INIT       = 0.089                      # s, bonded wall (the driver's period tracking)
_TAGR         = "MAXWELL15" if RAY_EFF <= 0 else "RAYM{:03d}".format(int(round(RAY_EFF * 1000)))   # named by the EFFECTIVE ratio
OVERRIDES = dict(                          # NAME = value lines rewritten in the driver (last assignment of each)
    OUT_DIR="stratF_full_LS_{}_char_tr".format(_TAGR),
    CASE_ID="F_LS_{}_char_tr".format(_TAGR),
    START_SAVE=SAVE_FT0,
    TENSION_MAX_PA=0.4e6,
    MAXWELL_CMD=MAXWELL_15,
    DAMP_RATIO=float(RAY_RATIO),
    DAMP_TYPE="mass" if RAY_RATIO > 0 else "",
    LARGE_STRAIN="on",
    T1_init=T1_INIT,
    RUN_FROM=1,
    RUN_TO=25,
    TAIL_SEC_RECORD=TAIL_S,
)
BASE_TABLES_TR = {                         # truncated tables (make_truncated_tables.py)
    "HU12": ("vel_HU_tr", "vel_HU_tr.txt"),
    "EC40": ("vel_EC_tr", "vel_EC_tr.txt"),
    "FR76": ("vel_FR_tr", "vel_FR_tr.txt"),
}
KEEP_AS_V12 = ("G_I_NEW", "COH_RESIDUAL")   # printed, not changed
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
TAG = "# [run_stratF_char_mx15_tr override]"


def log(s):
    print("[run_stratF_char_mx15_tr] " + s)


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
    tables = [fn for _, fn in BASE_TABLES_TR.values()]
    if MAKE_TABLES and not all(os.path.isfile(f) for f in tables):
        if not os.path.isfile("make_truncated_tables.py"):
            raise RuntimeError("truncated tables missing and make_truncated_tables.py not found")
        log("  writing the truncated tables (make_truncated_tables.py)")
        import importlib.util
        spec = importlib.util.spec_from_file_location("make_truncated_tables", os.path.abspath("make_truncated_tables.py"))
        mtt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mtt); mtt.main()
    miss = [f for f in [DRIVER] + tables if not os.path.isfile(f)]
    if not os.path.isfile(SAVE_FT0) and not (BUILD_DAT and os.path.isfile(BUILD_DAT)):
        miss.append("{} (and no {})".format(SAVE_FT0, BUILD_DAT))
    if miss:
        raise RuntimeError("missing in {}: {}".format(os.getcwd(), ", ".join(miss)))
    with open(DRIVER) as f:
        src = f.read()
    for name in KEEP_AS_V12:
        log("  kept   {} = {}".format(name, current_value(src, name)))
    log("  was    MAXWELL_CMD = {}".format(current_value(src, "MAXWELL_CMD")))
    out = src
    for name, val in OVERRIDES.items():
        out = patch_last(out, name, val)
        log("  set    {} = {!r}".format(name, val))
    out = patch_tables(out, BASE_TABLES_TR)
    log("  set    BASE_TABLES = {}".format({k: v[1] for k, v in BASE_TABLES_TR.items()}))
    for name in OVERRIDES:
        if TAG not in current_value(out, name):
            raise RuntimeError("override of {} did not take (a later assignment wins?)".format(name))
    if TAG not in current_value(out, "BASE_TABLES"):
        raise RuntimeError("override of BASE_TABLES did not take")
    code = compile(out, DRIVER + " [char mx15 tr]", "exec")
    log("  patched driver compiles")
    if os.path.isfile(os.path.join(OVERRIDES["OUT_DIR"], "stratC_checkpoint.json")):
        log("  checkpoint found in {}: the driver resumes after its last finished run".format(OVERRIDES["OUT_DIR"]))
    if DRY_RUN:
        log("DRY_RUN: stopping before 3DEC (set DRY_RUN = False to run)")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC (or set DRY_RUN = True)")
    if not os.path.isfile(SAVE_FT0):
        log("  {} missing: program call '{}'".format(SAVE_FT0, BUILD_DAT))
        it.command("python-reset-state false")
        it.command("program call '{}'".format(BUILD_DAT))
        it.command("python-reset-state false")
        if not os.path.isfile(SAVE_FT0):
            raise RuntimeError("{} did not write {}".format(BUILD_DAT, SAVE_FT0))
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(SAVE_FT0.replace("\\", "/")))
    it.command("python-reset-state false")
    try:
        ft = float(it.fish.get("tension"))
    except Exception as e:
        raise RuntimeError("{}: FISH 'tension' not readable ({}) -- cannot verify the build".format(SAVE_FT0, e))
    lo, hi = SAVE_TENSION_RANGE
    log("  {}: joint tension {:.3f} MPa (expected {:.2f}-{:.2f})".format(SAVE_FT0, ft / 1e6, lo / 1e6, hi / 1e6))
    if not (lo <= ft <= hi):
        raise RuntimeError("wrong save: joint tension {:.3f} MPa, not the characterisation build".format(ft / 1e6))
    ns = {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)}
    exec(code, ns)
    log("done: {}".format(OVERRIDES["OUT_DIR"]))


main()
