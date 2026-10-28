# -*- coding: ascii -*-
"""
run_stratF_char_mx_tr.py -- Strategy F on the bonded characterisation build (f_t = 0.344 MPa, c = 0.75 MPa), large
strain, truncated records + 4 s free-vibration tail (as run_stratF_char_mx_tr.py), with the Maxwell damping level
set by DAMP_PCT. Levels with their own optimised (MATLAB) 1-25 Hz set are taken from MAXWELL_SETS as given:
    1.5 %  0.0112 1.0005  0.0081 7.4315  0.0096 25.0004   (v12)
    5 %    0.0388 1.0008  0.0230 5.4866  0.0473 24.9990
Any other DAMP_PCT falls back to the 1.5 % set with its ratios scaled by DAMP_PCT / 1.5 (a warning is printed).

Why (28 Oct): the same equivalent-SDOF identification on model and tests (exp_equivalent_damping.py) gave, for
runs 1-6 before cracking, xi_eq 2.2-4.7 % for the 1.5 % model against 3.4-18.6 % (US-1) and 5.4-10.8 % (US-2); the
model's free vibration rang for ~2 s where the tests died out in ~0.3 s. The model adds ~0.7 % to its Maxwell level
before cracking, so DAMP_PCT = 5 puts it at ~6 %, inside the tests' range. Re-check runs 1-6 and then 15-21 (where
the contacts add their own hysteresis) with exp_equivalent_damping.py --model, for over-damping.

    OUT_DIR    stratF_full_LS_MAXWELL<pct x10>_char_tr   (DAMP_PCT 5 -> stratF_full_LS_MAXWELL50_char_tr)
Everything else as run_stratF_char_mx_tr.py. Starts at run 1 in its own folder (a new damping cannot resume a
1.5 % checkpoint).

USAGE (3DEC 9):  DRY_RUN = True first (plain python is fine), then DRY_RUN = False and
                 program call 'run_stratF_char_mx_tr.py'
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
DAMP_PCT      = 5.0                        # target Maxwell damping (%), flat 1-25 Hz
MAXWELL_SETS  = {                          # optimised 1-25 Hz sets, (ratio, freq) x 3
    1.5: (0.0112, 1.0005, 0.0081, 7.4315, 0.0096, 25.0004),   # v12
    5.0: (0.0388, 1.0008, 0.0230, 5.4866, 0.0473, 24.9990),   # 28 Oct
}
if DAMP_PCT in MAXWELL_SETS:
    _MX = MAXWELL_SETS[DAMP_PCT]
else:
    _k  = DAMP_PCT / 1.5
    _b  = MAXWELL_SETS[1.5]
    _MX = (_b[0] * _k, _b[1], _b[2] * _k, _b[3], _b[4] * _k, _b[5])
    print("[run_stratF_char_mx_tr] WARNING: no optimised set for {} %, using the 1.5 % set x {:.3f}".format(DAMP_PCT, _k))
MAXWELL_15    = "block mech damp maxwell {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(*_MX)
RAY_RATIO     = 2.0 * RAY_EFF              # (Rayleigh mass term, unused here: RAY_EFF = 0)
T1_INIT       = 0.089                      # s, bonded wall (the driver's period tracking)
_TAGR         = "MAXWELL{:02d}".format(int(round(DAMP_PCT * 10))) if RAY_EFF <= 0 else "RAYM{:03d}".format(int(round(RAY_EFF * 1000)))   # named by the EFFECTIVE ratio
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
TAG = "# [run_stratF_char_mx_tr override]"


def log(s):
    print("[run_stratF_char_mx_tr] " + s)


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
    code = compile(out, DRIVER + " [char mx tr]", "exec")
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