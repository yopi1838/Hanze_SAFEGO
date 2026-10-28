# -*- coding: ascii -*-
"""
run_pushover_slipwork.py -- the v5 pushover with SLIP WORK TRACKED PER CONTACT REGION (slip_work_tracker.py).

Purpose (28 Oct): test the mechanical model's pocket-friction term (rigid_two_block_vw.py, POCKET CLAMP) against the
DEM directly, instead of through a model-wide energy total:
    * the lower pocket's clamp force N: the model predicts about 23.4 kN when pushed TOWARDS the joists and ~0 the
      other way;
    * its slip work: the model predicts mu N (0.088 m) x rotation, about 74 J at 80 mm in the positive push;
    * where the rest of the dissipation sits (upper pocket, base, top joint, the lower- and upper-floor courses).

What it does, per CASE and direction:
    1. reads the pushover driver (PUSHOVER_SCRIPT; the laptop copy is pushover_oop_3dec_v4.py, which holds v5);
    2. sets its config lines (RESTORE_SAV, OUT_DIR, PUSH_DIR, ...) and inserts, IN MEMORY ONLY:
         - slip_work_tracker.install() after the settled baseline (just before the ascent),
         - slip_work_tracker.update(d) after every solve chunk (ascent: CHUNK_ASC cycles, servo: CHUNK cycles)
           and at every logged pushover point (tag asc / desc),
         - slip_work_tracker.close() at the end;
    3. runs it. The pushover CSV is written as usual; the tracker writes OUT_DIR/slipwork_<pos|neg>.csv.
The driver file is not edited. Results go to their own folders (suffix _slipwork), nothing existing is overwritten.

MODES
    "dry"  : patch and compile for every case and direction; no model work (plain python is fine)
    "probe": first case, first direction, ascent stopped at A_MAX = PROBE_A_MAX (a few steps, minutes):
             check the log line '[slip_work_tracker] installed: ...' -- the lower pocket must have subcontacts,
             k_s must be read (failed = 0), and the CSV must show a sign of +1 or -1 and a lower-pocket N.
    "run"  : every case and direction (each about as long as the sweep's pushover, ~5-9 h)
USAGE (3DEC 9, working folder with the saves, the driver, slip_work_tracker.py):
    MODE = "dry" -> program call 'run_pushover_slipwork.py'; then "probe"; then "run"
Post-processing: python postprocess_slipwork.py
"""
import os, re, sys, time, importlib

# ============================ CONFIG =================================
MODE        = "dry"                           # "dry" | "probe" | "run"
PUSHOVER_SCRIPT = ("pushover_oop_3dec_v5.py", "pushover_oop_3dec_v4.py")   # first one found
CASES = [   # (label, save for PUSH_DIR +1, save for PUSH_DIR -1, OUT_DIR, TENSION_MAX_PA)
    ("ft0p05", "Part_I_MASON_v8_LargeStrain_springforce_ft_0p050", "Part_I_MASON_v8_LargeStrain_springforce_ft_0p050_n",
     "pushover_results_v5_ft_0p050_TopMason_slipwork", 60000.0),
    # ("ft0", "Part_I_MASON_v8_LargeStrain_springforce", "Part_I_MASON_v8_LargeStrain_springforce",
    #  "pushover_results_v5_ft0_TopMason_slipwork", 50000.0),
]
DIRS        = (+1, -1)
TOP_JOINT   = "mason"
PROBE_A_MAX = 2.0                             # m/s2: probe ascends 0.5 ... 2.0 and stops (elastic range)
LOG         = "run_pushover_slipwork.log"
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
TAG = "# [run_pushover_slipwork]"


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


def insert_after(src, pattern, text, count=1, what=""):
    """insert `text` (one line, indented like the matched line) after each line matching pattern"""
    hits = list(re.finditer(r"^([ \t]*)" + pattern + r".*$", src, re.M))
    if len(hits) != count:
        raise RuntimeError("{}: expected {} match(es) of {!r}, found {} -- the driver changed; update this wrapper".format(
            what, count, pattern, len(hits)))
    for m in reversed(hits):
        src = src[:m.end()] + "\n" + m.group(1) + text + "   " + TAG + src[m.end():]
    return src


def build(src, save, out_dir, push_dir, tmax, probe):
    over = dict(RESTORE_SAV=save, OUT_DIR=out_dir, PUSH_DIR=push_dir, TOP_JOINT=TOP_JOINT, TENSION_MAX_PA=tmax)
    if probe:
        over["A_MAX"] = PROBE_A_MAX
    for k, v in over.items():
        src = patch_last(src, k, v)
    # 1. install after the helpers, just before the ascent bookkeeping
    m = re.search(r"^rows = \[\][ \t]*$", src, re.M)
    if not m:
        raise RuntimeError("anchor 'rows = []' not found -- update this wrapper")
    inst = ("import sys as _sys; _sys.path.insert(0, os.getcwd()); import importlib as _il   {t}\n"
            "import slip_work_tracker as _swt; _swt = _il.reload(_swt)   {t}\n"
            "_swt.install(os.path.join(OUT_DIR, 'slipwork_{{}}.csv'.format('pos' if PUSH_DIR > 0 else 'neg')), "
            "note='save {{}}, dir {{:+d}}'.format(RESTORE_SAV, int(PUSH_DIR)))   {t}\n"
            "_swt.update(ctrl_disp(), tag='baseline')   {t}\n").format(t=TAG)
    src = src[:m.start()] + inst + src[m.start():]
    # 2. after every solve chunk
    src = insert_after(src, r'it\.command\("model solve ratio \{:g\} cycles \{:d\}"\.format\(RATIO, CHUNK_ASC\)\)',
                       "_swt.update(ctrl_disp())", 1, "ascent chunk")
    src = insert_after(src, r'it\.command\("model cycle \{:d\}"\.format\(CHUNK\)\)',
                       "_swt.update(ctrl_disp())", 1, "servo chunk")
    # 3. at every logged point
    src = insert_after(src, r'rows\.append\(row\("asc"', "_swt.update(d, tag='asc')", 1, "asc row")
    src = insert_after(src, r'rows\.append\(row\("desc"', "_swt.update(d, tag='desc')", 1, "desc row")
    # 4. close at the end
    src = insert_after(src, r'print\("-> \{\}"\.format\(csv_path\)\)', "_swt.close()", 1, "end")
    return src, over


def main():
    log("=" * 70)
    log("run_pushover_slipwork  MODE {}".format(MODE))
    drv = next((p for p in PUSHOVER_SCRIPT if os.path.isfile(p)), None)
    if drv is None:
        raise RuntimeError("none of {} found in {}".format(PUSHOVER_SCRIPT, os.getcwd()))
    if not os.path.isfile("slip_work_tracker.py"):
        raise RuntimeError("slip_work_tracker.py not found in {}".format(os.getcwd()))
    with open(drv) as f:
        src0 = f.read()
    if "masonry only" not in src0.lower() and "po_push" not in src0:
        log("  WARNING: {} does not look like the v5 (masonry-push) driver".format(drv))
    jobs = []
    for label, sp, sn, out, tmax in CASES:
        for dr in DIRS:
            jobs.append((label, sp if dr > 0 else sn, out, dr, tmax))
    if MODE == "probe":
        label, sv, out, dr, tmax = jobs[0]
        jobs = [(label, sv, out + "_probe", dr, tmax)]
    for label, sv, out, dr, tmax in jobs:
        code_src, over = build(src0, sv, out, dr, tmax, MODE == "probe")
        code = compile(code_src, "{} [slipwork {} {:+d}]".format(drv, label, dr), "exec")
        log("{} dir {:+d}: {} -> {} (driver {}, patched and compiled)".format(label, dr, sv, out, drv))
        if MODE == "dry":
            continue
        if it is None:
            raise RuntimeError("itasca not available: run inside 3DEC, or MODE = 'dry'")
        if not (os.path.isfile(sv) or os.path.isfile(sv + ".sav")):
            log("  ! save {} missing -- skipped".format(sv)); continue
        csvp = os.path.join(out, "pushover_{}.csv".format("pos" if dr > 0 else "neg"))
        if MODE == "run" and os.path.isfile(csvp):
            log("  {} exists -- skipped (delete it to re-run)".format(csvp)); continue
        t0 = time.time()
        ns = {"__name__": "__main__", "__file__": os.path.abspath(drv)}
        exec(code, ns)
        log("  done in {:.1f} min -> {}".format((time.time() - t0) / 60.0, out))
    log("finished" if MODE != "dry" else "dry run: all patches apply and compile; set MODE = 'probe'")


main()
