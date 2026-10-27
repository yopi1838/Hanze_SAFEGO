# -*- coding: ascii -*-
"""
run_pushover_hinge_nofric.py -- the forced-hinge pushover again, with FRICTIONLESS joist pockets.

Why: the forced-hinge push (run_pushover_hinge_ft0.py) has the rigid model's mechanism exactly, yet carries 2-5 kN
more. Its energy balance (27 Oct) puts the whole gap in non-conservative work (152 J of 420 J at 80 mm), most likely
friction of the lower joist sliding in its pocket (44 kN normal force, 36 deg friction inherited from the masonry
joint). This run removes that friction and nothing else. If the curve then drops to the rigid line (plus the small
support drag), the pocket friction is identified.

WHAT IT CHANGES (everything else is the forced-hinge build: hinges at the base / 2.110 m / top, 0.344 MPa elsewhere)
    joist - masonry contacts (both floors): jmodel mohr, friction 0, cohesion 0, tension 0, same normal and shear
    stiffness as the masonry joints (kn, ks of mat_case.dat). The joists still rest in their pockets (normal contact),
    they just cannot transmit friction; the pocket back face still pushes them along on the positive push.
    Joist supports (Joist - Joist_S1/S2) unchanged.
STEPS
    1. restore SAVE_IN (the forced-hinge build), check FISH 'tension' (0.30-0.40 MPa);
    2. reassign the pocket contacts, list how many changed, re-settle (model solve ratio 1e-5), switch on the
       mechanical energy bookkeeping (ENERGY = True: 3DEC then records the slip work of the contacts in the save,
       so the remaining dissipation can be read from the snapshot saves), save SAVE_OUT;
    3. pushover_oop_3dec_v5.py from SAVE_OUT for each direction in DIRS, snapshots at SNAP_MM, into PO_DIR.
USAGE (3DEC 9, working folder with the forced-hinge save and the pushover script; needs the licence)
    MODE = "dry" : patches and compiles the pushover; no model work
    MODE = "prep": steps 1-2 only (check the log: number of pocket contacts changed, settle)
    MODE = "run" : steps 1-3
Post-processing: add PO_DIR as a case in postprocess_pushover_oop.py (CASES), and run_deformed_shape_batch.py with
REF = SAVE_OUT for the deformed shape.
"""
import os, re, time

# ============================ CONFIG =================================
MODE        = "dry"                                   # "dry" | "prep" | "run"
PUSHOVER_SCRIPT = "pushover_oop_3dec_v4.py"
SAVE_IN     = "Part_I_MASON_v8_LargeStrain_springforce_hingeft0"
SAVE_OUT    = "Part_I_MASON_v8_LargeStrain_springforce_hingeft0_nofric"
FT_RANGE    = (0.30e6, 0.40e6)                        # SAVE_IN must hold this joint tension (FISH 'tension')
PO_DIR      = "pushover_results_v5_hingeft0_nofric_TopMason"
DIRS        = (+1, -1)
SNAP_MM     = (5, 10, 20, 40)
POCKET_PROPS = "stiffness-normal @kn stiffness-shear @ks friction 0 cohesion 0 tension 0"
ENERGY      = True                                    # model energy mechanical on (slip work kept in the saves)
SETTLE      = "model solve ratio 1e-5 cycles 150000"
PO_OVER     = dict(TOP_JOINT="mason", TENSION_MAX_PA=0.4e6)
LOG         = "run_pushover_hinge_nofric.log"
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
POCKET = "group-intersection 'Joist' 'Masonry'"


def log(s=""):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), s)
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def cmdpath(p):
    return p.replace("\\", "/")


def patch(src, over):
    for name, val in over.items():
        m = re.search(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M)
        if m is None:
            raise RuntimeError("{}: config line '{} = ...' not found".format(PUSHOVER_SCRIPT, name))
        src = src[:m.start()] + "{} = {!r}   # [hinge_nofric override]".format(name, val) + src[m.end():]
    return src


SNAP_ANCHOR = ('        rows.append(row("desc", a, d, Fb, Ft, Ftot, settled, m_part, a_pred, nch, dg))\n'
               '        flush()\n')
SNAP_CODE = ('        for _s in _BATCH_SNAP_MM:\n'
             '            if _s not in _batch_snapped and abs(d) >= _s:\n'
             '                it.command("model save \'{}\'".format(os.path.join(OUT_DIR, "pushover_{}_d{:02d}".format(\n'
             '                    "pos" if PUSH_DIR > 0 else "neg", int(_s))).replace("\\\\", "/")))\n'
             '                _batch_snapped.add(_s)\n'
             '                print("  snapshot save at d = {:.1f} mm (target {} mm)".format(d, _s))\n')


def pushover_source(push_dir):
    with open(PUSHOVER_SCRIPT) as f:
        src = f.read()
    over = dict(PO_OVER, RESTORE_SAV=SAVE_OUT, OUT_DIR=PO_DIR, PUSH_DIR=int(push_dir))
    src = patch(src, over)
    if SNAP_MM:
        if src.count(SNAP_ANCHOR) == 1:
            src = "_BATCH_SNAP_MM = {!r}\n_batch_snapped = set()\n".format(tuple(SNAP_MM)) + src.replace(SNAP_ANCHOR, SNAP_ANCHOR + SNAP_CODE)
        else:
            log("  ! snapshot anchor not found -- no snapshot saves (the pushover is unaffected)")
    return src, over


FISH_COUNT = r"""
fish define hnf_count
    local n = 0
    local nf = 0
    loop foreach local cx block.contact.list
        local b1 = block.contact.b1(cx)
        local b2 = block.contact.b2(cx)
        local j = 0
        if block.isgroup(b1, 'Joist') then
            if block.isgroup(b2, 'Masonry') then
                j = 1
            endif
        endif
        if block.isgroup(b2, 'Joist') then
            if block.isgroup(b1, 'Masonry') then
                j = 1
            endif
        endif
        if j = 1 then
            n = n + 1
            if block.contact.prop(cx, 'friction') < 1e-6 then
                nf = nf + 1
            endif
        endif
    endloop
    hnf_count = string(n) + ' pocket contacts, ' + string(nf) + ' with friction 0'
end
"""


def prepare():
    os.makedirs(PO_DIR, exist_ok=True)
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(cmdpath(SAVE_IN)))
    it.command("python-reset-state false")
    ft = float(it.fish.get("tension"))
    log("{}: joint tension {:.3f} MPa".format(SAVE_IN, ft / 1e6))
    if not (FT_RANGE[0] <= ft <= FT_RANGE[1]):
        raise RuntimeError("wrong save: joint tension {:.3f} MPa (expected the forced-hinge build)".format(ft / 1e6))
    it.command("block contact jmodel assign mohr range {}".format(POCKET))
    it.command("block contact property {} range {}".format(POCKET_PROPS, POCKET))
    log("  pocket contacts -> mohr, {}".format(POCKET_PROPS))
    try:
        it.command(FISH_COUNT)
        it.command("[global hnf_msg = hnf_count]")
        log("  " + str(it.fish.get("hnf_msg")))
    except Exception as e:
        log("  (pocket contact count not available: {})".format(e))
    log("  " + SETTLE)
    it.command(SETTLE)
    if ENERGY:
        it.command("model energy mechanical on")
        log("  model energy mechanical on")
    it.command("model save '{}'".format(cmdpath(SAVE_OUT)))
    log("  saved {}.sav".format(SAVE_OUT))


def main():
    log("=" * 70)
    log("run_pushover_hinge_nofric  mode {}".format(MODE))
    if not os.path.isfile(PUSHOVER_SCRIPT):
        raise RuntimeError("{} not found in {}".format(PUSHOVER_SCRIPT, os.getcwd()))
    for d in DIRS:
        src, over = pushover_source(d)
        compile(src, PUSHOVER_SCRIPT, "exec")
        log("pushover {}: {} -- patched source compiles".format("POS" if d > 0 else "NEG",
                                                             ", ".join("{}={!r}".format(k, v) for k, v in over.items())))
    if MODE == "dry":
        log("dry run: no model work. MODE = 'prep' or 'run'.")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC")
    if not (os.path.isfile(SAVE_IN) or os.path.isfile(SAVE_IN + ".sav")):
        raise RuntimeError("missing {}.sav (build it with ANALYSIS_PART_I_MASON_hingeft0.dat)".format(SAVE_IN))
    prepare()
    if MODE == "prep":
        return
    for d in DIRS:
        tag = "pos" if d > 0 else "neg"
        csvp = os.path.join(PO_DIR, "pushover_{}.csv".format(tag))
        if os.path.isfile(csvp):
            os.rename(csvp, csvp.replace(".csv", "_prev_{}.csv".format(time.strftime("%Y%m%d_%H%M%S"))))
        log("-" * 70)
        log("pushover {} from {}".format(tag.upper(), SAVE_OUT))
        src, _ = pushover_source(d)
        t0 = time.time()
        try:
            exec(compile(src, PUSHOVER_SCRIPT, "exec"), {"__name__": "__main__", "__file__": os.path.abspath(PUSHOVER_SCRIPT)})
        except SystemExit:
            pass
        except Exception as e:
            log("  pushover {} FAILED: {}".format(tag, e)); continue
        log("  pushover {} done in {:.1f} min -> {}".format(tag, (time.time() - t0) / 60, csvp))
    log("finished")


main()