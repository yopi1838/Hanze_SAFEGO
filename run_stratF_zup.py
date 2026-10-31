# -*- coding: ascii -*-
"""
run_stratF_zup.py -- Strategy F (the recorded table motions, runs 1-25) on a Z-UP build (31 Oct 2026).

Runs strategy_F_full_3dec_LS.py with its config AND its frame-dependent lines rewritten in memory (the driver file is
not edited), as run_stratf_char_mx15_tr.py does for the config. Same case set-up as the 5 % characterisation run
(stratF_full_LS_MAXWELL50_char_tr): large strain, truncated records (vel_*_tr.txt) + 4 s free tail, Maxwell set by
DAMP_PCT (optimised 1-25 Hz sets), no Rayleigh, G_I / cohesion-residual of the save untouched.

FRAME (Geo_Prep_zup*.dat / ANALYSIS_PART_I_MASON_zup*.dat): X = 1.292 - x in-plane, Y = z out-of-plane (+Y towards the
joists), Z = y vertical. The rewrites, each checked for its exact number of matches (the run stops if the driver
changed):
    table and top-beam drive          block apply / apply-remove velocity-z  ->  velocity-y
    base block S freed                block free velocity-z range group 'S'  ->  velocity-y
    spring-force check                comp.y(block.force.app(b))             ->  comp.z(...)
    joist-support contact groups      pos-y (height) / pos-z (OOP) ranges    ->  pos-z / pos-y at this build's joists
    instrumentation                   instrument_history_new.dat / instrument_baseshear_exp.dat -> the _zup versions
                                      (the first calls instrument_tilt_v2_zup.dat); joist slip probes set to this
                                      build's joist centres after every call; the export file is unchanged
    frame guard                       after the restore: stops if the top beam is not at Z > 2.5 m (old-frame save)
T_B rotation-y and rotation-z are freed as in the driver (the same pair in both frames).

BUILD picks the save, its tension check and the joist heights:
    "zup_pk2_ft0p05"  Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05   f_t 0.05 MPa, test pocket heights
    "zup_pk2"         Part_I_MASON_v8_LargeStrain_springforce_zup_pk2          f_t 0,        test pocket heights
    "zup_pk2_char"    Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char     0.344 / 0.75 MPa, test pocket heights
    "zup"             Part_I_MASON_v8_LargeStrain_springforce_zup              f_t 0,        old pocket heights
OUT_DIR  stratF_full_LS_MAXWELL<pct x10>_<BUILD>_tr   (e.g. stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr)
The Route 2 z-up driver with the same BUILD and DAMP_PCT is the surrogate scored against this run.

NEEDS next to it: strategy_F_full_3dec_LS.py, period_id_exp.py, the save of BUILD (or its ANALYSIS_..._zup_pk2*.dat,
built on the fly), vel_HU_tr.txt / vel_EC_tr.txt / vel_FR_tr.txt (or make_truncated_tables.py), instrument_history_new_zup.dat,
instrument_tilt_v2_zup.dat, instrument_baseshear_exp_zup.dat, instrument_history_export_v2.dat.
USAGE (3DEC 9):  DRY_RUN = True first (plain python is fine: patches, checks the match counts, compiles), then
                 DRY_RUN = False and  program call 'run_stratF_zup.py'
Post-processing: python postprocess_stratF_full25.py --runs "pk2 f_t 0.05, 5 %=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr"
"""
import os, re

# ============================ CONFIG =================================
DRY_RUN   = True
BUILD     = "zup_pk2_ft0p05"            # "zup_pk2_ft0p05" | "zup_pk2" | "zup_pk2_char" | "zup"
DAMP_PCT  = 5.0                         # a key of MAXWELL_SETS
TAIL_S    = 4.0                         # s, free vibration after each record (as the 5 % char run)
DRIVER    = "strategy_F_full_3dec_LS.py"
MAKE_TABLES = True                      # write the truncated tables first if they are missing
_BUILDS = {   # save, build file, tension range of the save (Pa), joist bottoms Z, joist centres Z, T1_init (s)
    "zup_pk2_ft0p05": ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05.sav", "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05.dat",
                       (0.04e6, 0.06e6), (0.23455, 2.11091), (0.32250, 2.19886), 0.091),
    "zup_pk2":        ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2.sav", "ANALYSIS_PART_I_MASON_zup_pk2.dat",
                       (0.0, 0.01e6), (0.23455, 2.11091), (0.32250, 2.19886), 0.095),
    "zup_pk2_char":   ("Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_char.sav", "ANALYSIS_PART_I_MASON_zup_pk2_char.dat",
                       (0.30e6, 0.40e6), (0.23455, 2.11091), (0.32250, 2.19886), 0.089),
    "zup":            ("Part_I_MASON_v8_LargeStrain_springforce_zup.sav", "ANALYSIS_PART_I_MASON_zup.dat",
                       (0.0, 0.01e6), (0.29318, 2.05227), (0.38114, 2.14023), 0.100),
}
MAXWELL_SETS = {                        # optimised 1-25 Hz sets (ratio, freq) x 3, as run_stratf_char_mx15_tr.py
    1.5: (0.0112, 1.0005, 0.0081, 7.4315, 0.0096, 25.0004),
    5.0: (0.0388, 1.0008, 0.0230, 5.4866, 0.0473, 24.9990),
}
BASE_TABLES_TR = {                      # truncated tables (make_truncated_tables.py)
    "HU12": ("vel_HU_tr", "vel_HU_tr.txt"),
    "EC40": ("vel_EC_tr", "vel_EC_tr.txt"),
    "FR76": ("vel_FR_tr", "vel_FR_tr.txt"),
}
# ======================================================================
SAVE, BUILD_DAT, TRANGE, JB, JC, T1_INIT = _BUILDS[BUILD]
_TAGR = "MAXWELL{:02d}".format(int(round(DAMP_PCT * 10)))
OVERRIDES = dict(
    OUT_DIR="stratF_full_LS_{}_{}_tr".format(_TAGR, BUILD),
    CASE_ID="F_LS_{}_{}_tr".format(_TAGR, BUILD),
    START_SAVE=SAVE,
    TENSION_MAX_PA=TRANGE[1],
    MAXWELL_CMD="block mech damp maxwell {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(*MAXWELL_SETS[DAMP_PCT]),
    DAMP_RATIO=0.0,
    DAMP_TYPE="",
    LARGE_STRAIN="on",
    T1_init=T1_INIT,
    RUN_FROM=1,
    RUN_TO=25,
    TAIL_SEC_RECORD=TAIL_S,
)
TAG = "# [run_stratF_zup]"

try:
    import itasca as it
except ImportError:
    it = None


def log(s):
    print("[run_stratF_zup] " + s)


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


def rep(src, old, new, count, what):
    k = src.count(old)
    if k != count:
        raise RuntimeError("frame rewrite '{}': expected {} match(es), found {} -- the driver changed; update this "
                           "wrapper".format(what, count, k))
    return src.replace(old, new)


def frame_zup(src):
    """the z-up rewrites of the driver's frame-dependent lines (see the header)"""
    src = rep(src, "block apply velocity-z {:g} table", "block apply velocity-y {:g} table", 1, "table drive")
    src = rep(src, "block gridpoint apply-remove velocity-z range group", "block gridpoint apply-remove velocity-y range group",
              1, "drive removal")
    src = rep(src, "block free velocity-z range group 'S'", "block free velocity-y range group 'S'", 1, "free S")
    src = rep(src, "comp.y(block.force.app(b))", "comp.z(block.force.app(b))", 1, "spring-force check")
    src = rep(src, "block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5",
              "block contact group 'Joist_S1_contact' range pos-z {:.4f} {:.4f} pos-y 0.9 1.5".format(JB[0] - 0.03, JB[0] + 0.03),
              1, "Joist_S1_contact")
    src = rep(src, "block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5",
              "block contact group 'Joist_S2_contact' range pos-z {:.4f} {:.4f} pos-y 0.9 1.5".format(JB[1] - 0.03, JB[1] + 0.03),
              1, "Joist_S2_contact")
    # instrumentation calls (setup and after every run's history delete), each followed by the joist probe heights
    pat = re.compile(r"^([ \t]*)it\.command\(\"call 'instrument_history_new\.dat'\"\)[ \t]*$", re.M)
    if len(pat.findall(src)) != 2:
        raise RuntimeError("frame rewrite 'instrument_history_new calls': expected 2, found {}".format(len(pat.findall(src))))
    src = pat.sub(lambda m: ("{0}it.command(\"call 'instrument_history_new_zup.dat'\")   {1}\n"
                             "{0}it.command(\"[jy_lower = {2:.5f}]\")   {1}\n"
                             "{0}it.command(\"[jy_upper = {3:.5f}]\")   {1}").format(m.group(1), TAG, JC[0], JC[1]), src)
    src = rep(src, "it.command(\"call 'instrument_baseshear_exp.dat'\")",
              "it.command(\"call 'instrument_baseshear_exp_zup.dat'\")", 2, "instrument_baseshear_exp calls")
    src = rep(src, '''DAT_FILES     = ["instrument_history_new.dat",
                 "instrument_baseshear_exp.dat",''', '''DAT_FILES     = ["instrument_history_new_zup.dat", "instrument_tilt_v2_zup.dat",   {}
                 "instrument_baseshear_exp_zup.dat",'''.format(TAG), 1, "DAT_FILES")
    # frame guard right after the restore of setup_model_for_dynamic
    old = "    it.command(\"model restore '{}'\".format(cmd_path(save_file)))\n"
    guard = old + ('''    it.command("""
fish define _zup_tb_z
    local zmax = -1e30
    loop foreach local b block.list
        if block.isgroup(b, 'T_B')
            zmax = math.max(zmax, comp.z(block.pos(b)))
        endif
    endloop
    _zup_tb_z = zmax
end
[global _zup_tbz = _zup_tb_z]
""")   {0}
    _ztb = float(it.fish.get("_zup_tbz"))   {0}
    print("  frame check (z-up): top-beam centroid Z = {{:.3f}} m".format(_ztb))   {0}
    if _ztb < 2.5:   {0}
        raise RuntimeError("top beam at Z = {{:.3f}} m: an OLD-frame save; this is the z-up runner".format(_ztb))   {0}
''').format(TAG)
    src = rep(src, old, guard, 1, "restore in setup_model_for_dynamic")
    if re.search(r"velocity-z [{]|apply-remove velocity-z|free velocity-z", src):
        raise RuntimeError("a velocity-z drive line survived the rewrite")
    return src


def current_value(src, name):
    hits = re.findall(r"^{}[ \t]*=[ \t]*(.*)$".format(re.escape(name)), src, re.M)
    return hits[-1].strip() if hits else "(not found)"


def main():
    log("{}  BUILD {}  DAMP {} %  ->  {}".format("DRY RUN" if DRY_RUN else "RUN", BUILD, DAMP_PCT, OVERRIDES["OUT_DIR"]))
    tables = [fn for _, fn in BASE_TABLES_TR.values()]
    if MAKE_TABLES and not all(os.path.isfile(f) for f in tables) and not DRY_RUN:
        if not os.path.isfile("make_truncated_tables.py"):
            raise RuntimeError("truncated tables missing and make_truncated_tables.py not found")
        log("  writing the truncated tables (make_truncated_tables.py)")
        import importlib.util
        spec = importlib.util.spec_from_file_location("make_truncated_tables", os.path.abspath("make_truncated_tables.py"))
        mtt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mtt); mtt.main()
    need = [DRIVER, "instrument_history_new_zup.dat", "instrument_tilt_v2_zup.dat", "instrument_baseshear_exp_zup.dat",
            "instrument_history_export_v2.dat"] + tables
    miss = [f for f in need if not os.path.isfile(f)]
    if not os.path.isfile(SAVE) and not os.path.isfile(BUILD_DAT):
        miss.append("{} (and no {})".format(SAVE, BUILD_DAT))
    if miss:
        if DRY_RUN:
            log("  ! missing (fine for a dry run elsewhere, not for the run): " + ", ".join(miss))
        else:
            raise RuntimeError("missing in {}: {}".format(os.getcwd(), ", ".join(miss)))
    with open(DRIVER) as f:
        src = f.read()
    out = src
    for name, val in OVERRIDES.items():
        out = patch_last(out, name, val)
        log("  set    {} = {!r}".format(name, val))
    out = patch_tables(out, BASE_TABLES_TR)
    out = frame_zup(out)
    log("  frame rewrites applied: drive/free along Y, spring check on Z, joist contacts at Z {:.4f} / {:.4f}, "
        "instrumentation _zup (joist probes Z {:.4f} / {:.4f}), frame guard".format(JB[0], JB[1], JC[0], JC[1]))
    for name in OVERRIDES:
        if TAG not in current_value(out, name):
            raise RuntimeError("override of {} did not take (a later assignment wins?)".format(name))
    code = compile(out, DRIVER + " [z-up {}]".format(BUILD), "exec")
    log("  patched driver compiles")
    if os.path.isfile(os.path.join(OVERRIDES["OUT_DIR"], "stratC_checkpoint.json")):
        log("  checkpoint found in {}: the driver resumes after its last finished run".format(OVERRIDES["OUT_DIR"]))
    if DRY_RUN:
        log("DRY_RUN: stopping before 3DEC (set DRY_RUN = False to run)")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC (or set DRY_RUN = True)")
    if not os.path.isfile(SAVE):
        log("  {} missing: program call '{}'".format(SAVE, BUILD_DAT))
        it.command("python-reset-state false")
        it.command("program call '{}'".format(BUILD_DAT))
        it.command("python-reset-state false")
        if not os.path.isfile(SAVE):
            raise RuntimeError("{} did not write {}".format(BUILD_DAT, SAVE))
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(SAVE.replace("\\", "/")))
    it.command("python-reset-state false")
    try:
        ft = float(it.fish.get("tension"))
    except Exception as e:
        raise RuntimeError("{}: FISH 'tension' not readable ({}) -- cannot verify the build".format(SAVE, e))
    lo, hi = TRANGE
    log("  {}: joint tension {:.3f} MPa (expected {:.3f}-{:.3f})".format(SAVE, ft / 1e6, lo / 1e6, hi / 1e6))
    if not (lo <= ft <= hi):
        raise RuntimeError("wrong save for BUILD {}: joint tension {:.3f} MPa".format(BUILD, ft / 1e6))
    ns = {"__name__": "__main__", "__file__": os.path.abspath(DRIVER)}
    exec(code, ns)
    log("done: {}".format(OVERRIDES["OUT_DIR"]))


main()
