# -*- coding: ascii -*-
"""
run_stratF_char_raym.py -- Strategy F (recorded table motions, runs 1-25) on the BONDED characterisation build
(f_t = 0.344 MPa, c = 0.75 MPa on every joint), LARGE STRAIN, with only a SMALL MASS-PROPORTIONAL Rayleigh damping
(RAY_EFF at the wall's first frequency 1 / T1_init; no local, no global, no Maxwell). Set RAY_EFF = 0.0 for
hysteretic damping only (contact friction and softening).

MASS-PROPORTIONAL DAMPING: xi(f) = RAY_EFF * f1 / f. 3DEC is given RAY_RATIO = 2 x RAY_EFF, because with the "mass"
keyword it applies only the mass half of the Rayleigh pair defined at (ratio, f1). It does not shrink the time step (stiffness-proportional
would), but it grows at LOW frequency: 0.5 % at f1 = 11.2 Hz is ~2 % at the rocking frequency (T_eff ~0.36 s,
2.8 Hz) and acts on the ABSOLUTE velocity, so it also drags on the wall's rigid-body motion with the table. Keep it
small; compare against RAY_EFF = 0 and the 1.5 % Maxwell run.

Changes against stratF_full_LS_MAXWELL15_v12 (f_t = 0 build, 1.5 % Maxwell):
    START_SAVE   Part_I_MASON_v8_LargeStrain_springforce_char.sav   (built by ANALYSIS_PART_I_MASON_char.dat if missing)
    MAXWELL_CMD  ''; DAMP_RATIO = RAY_RATIO (= 2 x RAY_EFF), DAMP_TYPE "mass" -> after 'block mech damp local 0.0 /
                 global 0.0' the driver issues 'block mechanical damping rayleigh <RAY_RATIO> <1/T1_init> mass'
                 (RAY_EFF = 0: nothing is added, "contact dissipation only")
    TENSION_MAX_PA 0.4e6 (the driver's own guard), T1_init 0.089 s (bonded wall; f_t = 0.15 build: 0.089 s)
    OUT_DIR / CASE_ID  stratF_full_LS_RAYM<ratio>_char / F_LS_RAYM<ratio>_char  (RAY_RATIO 0: ..._NODAMP_char)
Kept: large strain, Table-2 protocol, T_B rotation-y/z freed, the original vel_*.txt tables.
Before anything runs it restores START_SAVE and checks FISH 'tension' is within SAVE_TENSION_RANGE.
Expect: with no viscous damping the elastic runs ring down only through contact hysteresis (very little while the
wall is bonded), so free-vibration tails last long and early-run peaks can exceed the 1.5 % runs.

USAGE (3DEC 9, working folder = the one holding the saves, the driver and vel_*.txt):
    DRY_RUN = True :  program call 'run_stratF_char_raym.py'   (or plain python) -> files, overrides, compile
    DRY_RUN = False:  program call 'run_stratF_char_raym.py'   (builds the save first if it is missing)
Post-processing:  python postprocess_stratF_full25.py --runs "f_t 0, 1.5 %=stratF_full_LS_MAXWELL15_v12" \\
                                                       --runs "0.344 MPa, 0.5 % mass=stratF_full_LS_RAYM005_char"
"""
import os, re

# ============================ CONFIG =================================
DRY_RUN       = True
DRIVER        = "strategy_F_full_3dec_LS.py"
SAVE_FT0      = "Part_I_MASON_v8_LargeStrain_springforce_char.sav"
BUILD_DAT     = "ANALYSIS_PART_I_MASON_char.dat"   # run first if SAVE_FT0 is missing (None: never)
SAVE_TENSION_RANGE = (0.30e6, 0.40e6)     # Pa, FISH tension of the save (the characterisation build)
USE_BC_TABLES = False                      # False: v12's vel_*.txt
RAY_EFF       = 0.005                      # EFFECTIVE mass-proportional damping ratio at f1 = 1 / T1_init (0.0 = none)
RAY_RATIO     = 2.0 * RAY_EFF                 # what 3DEC is given: with the "mass" keyword 3DEC keeps only the mass half of the
                                              # Rayleigh pair defined at (ratio, f1), i.e. ratio / 2 at f1 -> input 2 x RAY_EFF
T1_INIT       = 0.089                      # s, bonded wall; sets f1 for the Rayleigh term
_TAGR         = "NODAMP" if RAY_EFF <= 0 else "RAYM{:03d}".format(int(round(RAY_EFF * 1000)))   # named by the EFFECTIVE ratio
OVERRIDES = dict(                          # NAME = value lines rewritten in the driver (last assignment of each)
    OUT_DIR="stratF_full_LS_{}_char".format(_TAGR),
    CASE_ID="F_LS_{}_char".format(_TAGR),
    START_SAVE=SAVE_FT0,
    TENSION_MAX_PA=0.4e6,
    MAXWELL_CMD="",
    DAMP_RATIO=float(RAY_RATIO),
    DAMP_TYPE="mass" if RAY_RATIO > 0 else "",
    LARGE_STRAIN="on",
    T1_init=T1_INIT,
    RUN_FROM=1,
    RUN_TO=25,
)
BASE_TABLES_BC = {
    "HU12": ("vel_HU", "vel_HU_bc.txt"),
    "EC40": ("vel_EC", "vel_EC_bc.txt"),
    "FR76": ("vel_FR", "vel_FR_bc.txt"),
}
KEEP_AS_V12 = ("G_I_NEW", "COH_RESIDUAL", "BASE_TABLES")   # printed, not changed
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
TAG = "# [run_stratF_char_raym override]"


def log(s):
    print("[run_stratF_char_raym] " + s)


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
    miss = [f for f in [DRIVER] + tables if not os.path.isfile(f)]
    if not os.path.isfile(SAVE_FT0) and not (BUILD_DAT and os.path.isfile(BUILD_DAT)):
        miss.append("{} (and no {})".format(SAVE_FT0, BUILD_DAT))
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
    code = compile(out, DRIVER + " [char nodamp]", "exec")
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
