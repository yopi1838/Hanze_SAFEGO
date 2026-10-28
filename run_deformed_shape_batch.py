# -*- coding: ascii -*-
"""
run_deformed_shape_batch.py -- deformed_shape_3dec.py over every pushover of the batch, in one 3DEC session.

For each CASE (a build + one push direction) it
  * finds the snapshot saves the batch wrote in that pushover folder: pushover_<dir>_d05/_d10/_d20/_d40 and
    pushover_<dir>_end (whatever exists, in order of displacement);
  * runs deformed_shape_3dec.py with REF_SAV = the save the pushover started from, SAVES = those snapshots and
    OUT_DIR = its own folder -- the config lines are rewritten in memory, deformed_shape_3dec.py on disk is not changed.
Each case restores its reference save first, so the cases do not share model state. One figure per case shows the
snapshot sequence (light = early, dark = late): where the first hinge opens and how the mechanism develops.

    deformed_shape_ft0_pos/     f_t = 0,    positive push   (REF: ..._springforce)
    deformed_shape_ft0_neg/     f_t = 0,    negative push
    deformed_shape_ft0p15_pos/  f_t = 0.15, positive push   (REF: ..._springforce_ft_0p15)
    deformed_shape_ft0p15_neg/  f_t = 0.15, negative push
    deformed_shape_hinge_pos/   forced hinges (f_t = 0 at base, 2.110 m, top, joist pockets; 0.344 MPa elsewhere),
    deformed_shape_hinge_neg/     REF: ..._springforce_hingeft0 (ANALYSIS_PART_I_MASON_hingeft0.dat)
ONLY selects cases by output folder (empty = all), so finished cases need not be repeated.
Per case: <tag>_summary.txt (course profile, bed-joint map, hinges, top-joint tilt), <tag>_courses/joints/blocks.csv,
deformed_shape.png and the raw _dump_<tag>.log files.

USAGE (3DEC 9, working folder = the one with the saves and the pushover_results_* folders; needs the licence, so run
it when the batch / Strategy F is not running):
    MODE = "dry"  : lists the cases and the snapshot saves found, patches and compiles; no 3DEC model work
    MODE = "run"  : program call 'run_deformed_shape_batch.py'
    MODE = "post" : re-processes the existing dumps of every case (plain python, no licence), e.g. after a plot change
Missing folders / saves are reported and skipped; a failing case is logged and the next one runs.
"""
import os, re, glob, traceback

# ============================ CONFIG =================================
MODE   = "dry"                                   # "dry" | "run" | "post"
SCRIPT = "deformed_shape_3dec.py"
SAVE_FT0   = "Part_I_MASON_v8_LargeStrain_springforce"
SAVE_FT015 = "Part_I_MASON_v8_LargeStrain_springforce_ft_0p15"
SAVE_HINGE = "Part_I_MASON_v8_LargeStrain_springforce_hingeft0"
CASES = [   # (output folder, reference save, pushover folder, direction "pos" | "neg")
    ("deformed_shape_ft0_pos",    SAVE_FT0,   "pushover_results_v5_ft0_masonrypush_TopMason",    "pos"),
    ("deformed_shape_ft0_neg",    SAVE_FT0,   "pushover_results_v5_ft0_masonrypush_TopMason",    "neg"),
    ("deformed_shape_ft0p15_pos", SAVE_FT015, "pushover_results_v5_ft0p15_masonrypush_TopMason", "pos"),
    ("deformed_shape_ft0p15_neg", SAVE_FT015, "pushover_results_v5_ft0p15_masonrypush_TopMason", "neg"),
    ("deformed_shape_hinge_pos",  SAVE_HINGE, "pushover_results_v5_hingeft0_char_TopMason",      "pos"),
    ("deformed_shape_hinge_neg",  SAVE_HINGE, "pushover_results_v5_hingeft0_char_TopMason",      "neg"),
]
ONLY = ["deformed_shape_hinge_pos", "deformed_shape_hinge_neg"]   # output folders to process; [] = every case
LOG = "run_deformed_shape_batch.log"
# ======================================================================


def log(s):
    print("[deformed_shape batch] " + s)
    with open(LOG, "a") as f:
        f.write(s + "\n")


def snapshots(folder, direction):
    """[(tag, save path without .sav, push dir)] ordered d05, d10, ..., end"""
    out = []
    for p in glob.glob(os.path.join(folder, "pushover_{}_d[0-9]*.sav".format(direction))):
        m = re.search(r"_d(\d+)\.sav$", p)
        if m:
            out.append((int(m.group(1)), p))
    out.sort()
    end = os.path.join(folder, "pushover_{}_end.sav".format(direction))
    if os.path.isfile(end):
        out.append((10 ** 6, end))
    sign = +1 if direction == "pos" else -1
    res = []
    for k, p in out:
        stem = os.path.splitext(p)[0].replace("\\", "/")
        tag = "{}_end".format(direction) if k == 10 ** 6 else "{}_d{:02d}".format(direction, k)
        res.append((tag, stem, sign))
    return res


def patch(src, ref, po_dir, saves, out_dir):
    def one(src, name, val):
        pat = re.compile(r"^{}[ \t]*=.*$".format(name), re.M)
        if not pat.search(src):
            raise RuntimeError("{}: config line '{} = ...' not found".format(SCRIPT, name))
        return pat.sub(lambda m: "{} = {!r}   # [batch]".format(name, val), src, count=1)
    src = one(src, "REF_SAV", ref)
    src = one(src, "PO_DIR", po_dir)
    src = one(src, "OUT_DIR", out_dir)
    blk = re.compile(r"^SAVES[ \t]*=[ \t]*\[.*?^\]", re.M | re.S)
    if not blk.search(src):
        raise RuntimeError("{}: SAVES block not found".format(SCRIPT))
    body = "SAVES = [   # [batch]\n" + "".join("    ({!r}, {!r}, {:+d}),\n".format(t, p, s) for t, p, s in saves) + "]"
    return blk.sub(lambda m: body, src, count=1)


def main():
    log("=" * 70)
    log("mode {}".format(MODE))
    if not os.path.isfile(SCRIPT):
        raise RuntimeError("{} not found in {}".format(SCRIPT, os.getcwd()))
    with open(SCRIPT) as f:
        src0 = f.read()
    import sys
    argv0 = list(sys.argv)
    todo = [c for c in CASES if not ONLY or c[0] in ONLY]
    log("cases: " + ", ".join(c[0] for c in todo))
    for out_dir, ref, po_dir, direction in todo:
        log("-" * 70)
        if MODE == "post":
            dumps = sorted(glob.glob(os.path.join(out_dir, "_dump_*.log")))
            saves = [(os.path.basename(d)[6:-4], "", +1 if os.path.basename(d).startswith("_dump_pos") else -1)
                     for d in dumps if not d.endswith("_dump_ref.log")]
            if not saves:
                log("{}: no dumps to re-process -- skipped".format(out_dir)); continue
        else:
            if not os.path.isdir(po_dir):
                log("{}: pushover folder {} not found -- skipped".format(out_dir, po_dir)); continue
            if not (os.path.isfile(ref) or os.path.isfile(ref + ".sav")):
                log("{}: reference save {}(.sav) not found -- skipped".format(out_dir, ref)); continue
            saves = snapshots(po_dir, direction)
            if not saves:
                log("{}: no pushover_{}_*.sav in {} -- skipped".format(out_dir, direction, po_dir)); continue
        log("{}: REF {} | {} save(s): {}".format(out_dir, ref, len(saves), ", ".join(t for t, _, _ in saves)))
        code = compile(patch(src0, ref, po_dir, saves, out_dir), "{} [{}]".format(SCRIPT, out_dir), "exec")
        if MODE == "dry":
            log("  patched script compiles (dry run)"); continue
        sys.argv = [SCRIPT] + (["--post"] if MODE == "post" else [])
        try:
            exec(code, {"__name__": "__main__", "__file__": os.path.abspath(SCRIPT)})
            log("  done -> {}".format(out_dir))
        except SystemExit:
            log("  done -> {}".format(out_dir))
        except Exception as e:
            log("  FAILED: {}".format(e))
            with open(LOG, "a") as f:
                traceback.print_exc(file=f)
        finally:
            sys.argv = argv0
    log("finished")


main()
