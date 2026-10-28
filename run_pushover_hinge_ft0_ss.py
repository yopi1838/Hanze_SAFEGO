# run_pushover_hinge_ft0_ss.py -- SMALL-STRAIN variant of run_pushover_hinge_ft0.py (27 Oct).
# Default: restores YOUR small-strain hinge build Part_I_MASON_v8_SmallStrain_springforce_hingeft0 (f_t = 0 at the
# base, the 2.110 m bed joint, the top joint and both pockets; 0.344 / 0.75 MPa elsewhere; built with
# ANALYSIS_PART_I_MASON_hingeft0_SS.dat if missing), checks it, and pushes both ways into
# pushover_results_v5_hingeft0_SS_TopMason. "large-strain on" in the pushover script is patched to "off".
# (SET_STRIPS / SET_REST = True: set the strips and the 0.344 MPa joints after the fact on a plain small-strain save.)
# -*- coding: ascii -*-
"""
run_pushover_hinge_ft0.py -- pushover on a bonded build with f_t = 0 only in the strips where the wall hinges.

Idea (L0 of the DEM-vs-mechanical-model ladder): give the DEM the two-block mechanism the mechanical model assumes.
Every joint has the characterisation-based properties (f_t = 0.344 MPa, c = 0.75 MPa; same stiffness, joists as
built), EXCEPT the joints in the hinge strips, which get the f_t = 0 build's joint (f_t 0, c 0.2 MPa, G_I 0):
    base joint                              y = -0.01 ... 0.01 m
    upper-floor course (the DEM's hinge)    y =  2.03 ... 2.25 m  (bed joints 2.051, 2.110, 2.168, 2.227 m, the head
                                                                   joints of the pocket course, the upper joist pocket)
    top joint (masonry - top beam)          y =  2.57 ... 2.59 m
If the DEM then lands on the rigid line (mechanical_model_US1, adjusted hinge 2.09 m, lower pivot at the base), the
two models agree once the mechanism is the same; what is left (plateau, joist clamp) is the DEM's own.

TWO WAYS TO GET THE MODEL (SET_STRIPS)
  False (default): SAVE_IN is built by ANALYSIS_PART_I_MASON_hingeft0.dat, which already sets the characterisation
         joint everywhere and the f_t = 0 joint in the strips BEFORE gravity and the springs are applied. This script
         only checks it (step 4) and pushes from it. If the save is missing and BUILD_DAT is set, it is built first.
  True:  SAVE_IN is any bonded save; steps 2-3 change the strips after the fact, then it settles and saves SAVE_OUT.

WHAT IT DOES
  1. restores SAVE_IN and checks its FISH 'tension' against FT_RANGE;
  2. (SET_STRIPS) puts the contacts of the strips into contact group 'HingeFT0' (slot 'hingeft0') by contact position, and takes
     the joist SUPPORT contacts (Joist - Joist_S1/S2, Mohr) out of it again;
  3. (SET_STRIPS) block contact property <HINGE_PROPS> range group 'HingeFT0' slot 'hingeft0';
  4. writes every strip contact (id, y, kind, tension read back) to <PO_DIR>/hinge_ft0_contacts.csv, a summary per
     joint height to the log, and the tension range of all other masonry contacts (expect 0.344 MPa); with
     SET_STRIPS it then re-equilibrates and saves SAVE_OUT;
  5. runs pushover_oop_3dec_v5.py from SAVE_OUT (= SAVE_IN without SET_STRIPS) for each direction in DIRS (config lines rewritten in memory, as the
     batch does), with snapshot saves at SNAP_MM.
The pushover's own guard reads FISH 'tension', which is the build's 0.344 MPa: TENSION_MAX_PA is set to 0.4 MPa.

JOIST POCKETS (POCKETS, SET_STRIPS only): the build assigns mason_v12 to every contact, Mohr only to the joist
supports, so the pocket contacts are bonded like the masonry; in the f_t = 0 build they have f_t = 0. The .dat build
corresponds to "strip".
    "strip" (default): pocket contacts inside a strip (the upper joist) get f_t = 0, the lower joist stays bonded
    "all":   every pocket contact gets f_t = 0 (joists behave as in the f_t = 0 build)
    "build": no pocket contact changes, not even in the upper strip

USAGE (3DEC 9, working folder = the one with the saves and pushover_oop_3dec_v5.py; needs the licence):
    MODE = "dry" : patches and compiles the pushover for each direction; no 3DEC model work (plain python is fine)
    The build (BUILD_DAT) unbonds the base, the 2.110 m bed joint, the top joint and both joist pockets.
    MODE = "prep": steps 1-4 only (the contact check, and SAVE_OUT with SET_STRIPS; check before pushing)
    MODE = "run" : steps 1-5        program call 'run_pushover_hinge_ft0.py'
Post-processing: add PO_DIR as a case in postprocess_pushover_oop.py (CASES), next to the f_t = 0 and f_t = 0.15 runs.
"""
import os, re, csv, time

# ============================ CONFIG =================================
MODE        = "dry"                                    # "dry" | "prep" | "run"
PUSHOVER_SCRIPT = "pushover_oop_3dec_v5.py"
SET_STRIPS  = False                                    # False: SAVE_IN already has the strips (small-strain hinge build)
SET_REST    = False                                    # True: first give every masonry joint the characterisation values
REST_PROPS  = {                                        #   (as ANALYSIS_PART_I_MASON_hingeft0.dat does before gravity)
    "tension": 0.344e6,
    "cohesion": 0.75e6,
    "cohesion-residual": 0.075e6,
    "G_I": 0.025 * ((2 * 0.344) ** 0.7) * 1000,        # Gt of mat_case.dat at 0.344 MPa (N/m)
}
REST_RANGES = ["group-intersection 'Masonry' 'Masonry'", "group-intersection 'Masonry' 'S'",
               "group-intersection 'Masonry' 'T_B'", "group-intersection 'Joist' 'Masonry'"]   # never the joist supports
BUILD_DAT   = "ANALYSIS_PART_I_MASON_hingeft0_SS.dat"     # run first if SAVE_IN is missing (None: never)
SAVE_IN     = "Part_I_MASON_v8_SmallStrain_springforce_hingeft0"   # your small-strain hinge build
SAVE_OUT    = "Part_I_MASON_v8_SmallStrain_springforce_hingeft0_set"   # SET_STRIPS only
FT_RANGE    = (0.30e6, 0.40e6)                         # SAVE_IN must hold this joint tension (FISH 'tension')
REST_FT_PA  = 0.344e6                                  # expected tension of every other masonry contact (check)
PO_DIR      = "pushover_results_v5_hingeft0_SS_TopMason"
DIRS        = (+1, -1)                                 # push directions, one pushover each
SNAP_MM     = (5, 10, 20, 40)                          # snapshot saves on the descending branch; () for none
BANDS = [                                              # (label, y_lo, y_hi): contacts whose position lies inside
    ("base joint",          -0.01, 0.01),              #   must match the ranges in BUILD_DAT
    ("bed joint 2.110 m",    2.105, 2.115),            #   the strip alternative: ("upper-floor course", 2.03, 2.25)
    ("top joint",            2.57, 2.59),
]
#  BANDS.append(("lower-floor course", 0.27, 0.55))    # optional: the lower-floor kink of the f_t = 0 runs
POCKETS     = "all"                                    # "strip" | "all" | "build"   (see the docstring)
HINGE_PROPS = {                                        # the f_t = 0 build's joint (mat_case.dat: tension 0,
    "tension": 0.0,                                    #   cohesion 0.2 MPa, G_I = Gt(0) = 0, residual 0.1 c)
    "G_I": 0.0,
    "cohesion": 0.2e6,
    "cohesion-residual": 0.02e6,
}
SETTLE      = "model solve ratio 1e-5 cycles 150000"   # after the property change, before SAVE_OUT
PO_OVER     = dict(TOP_JOINT="mason", TENSION_MAX_PA=0.4e6)
LOG         = "run_pushover_hinge_ft0_ss.log"
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
GRP, SLOT = "HingeFT0", "hingeft0"


def log(s=""):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), s)
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def cmdpath(p):
    return p.replace("\\", "/")


# ---------------- pushover patching (as batch_ft0p15_pipeline.py) ----------------
def patch(src, over):
    for name, val in over.items():
        pat = re.compile(r"^{}[ \t]*=.*$".format(re.escape(name)), re.M)
        m = pat.search(src)
        if m is None:
            raise RuntimeError("{}: config line '{} = ...' not found".format(PUSHOVER_SCRIPT, name))
        src = src[:m.start()] + "{} = {!r}   # [hinge_ft0 override]".format(name, val) + src[m.end():]
    return src


SNAP_ANCHOR = ('        rows.append(row("desc", a, d, Fb, Ft, Ftot, settled, m_part, a_pred, nch, dg))\n'
               '        flush()\n')
SNAP_CODE = ('        for _s in _BATCH_SNAP_MM:\n'
             '            if _s not in _batch_snapped and abs(d) >= _s:\n'
             '                it.command("model save \'{}\'".format(os.path.join(OUT_DIR, "pushover_{}_d{:02d}".format(\n'
             '                    "pos" if PUSH_DIR > 0 else "neg", int(_s))).replace("\\\\", "/")))\n'
             '                _batch_snapped.add(_s)\n'
             '                print("  snapshot save at d = {:.1f} mm (target {} mm)".format(d, _s))\n')


def add_snapshots(src):
    if not SNAP_MM:
        return src
    if src.count(SNAP_ANCHOR) != 1:
        log("  ! snapshot anchor not found in {} -- no snapshot saves (the pushover is unaffected)".format(PUSHOVER_SCRIPT))
        return src
    src = src.replace(SNAP_ANCHOR, SNAP_ANCHOR + SNAP_CODE)
    return "_BATCH_SNAP_MM = {!r}\n_batch_snapped = set()\n".format(tuple(SNAP_MM)) + src


def push_save():
    return SAVE_OUT if SET_STRIPS else SAVE_IN


def pushover_source(push_dir):
    with open(PUSHOVER_SCRIPT) as f:
        src = f.read()
    over = dict(PO_OVER, RESTORE_SAV=push_save(), OUT_DIR=PO_DIR, PUSH_DIR=int(push_dir))
    src = add_snapshots(patch(src, over))
    # SMALL STRAIN: the pushover must not switch large strain back on
    n_on = len(re.findall(r"large-strain\s+on", src))
    if n_on:
        src = re.sub(r"large-strain\s+on", "large-strain off", src)
        log("  {}: {} 'large-strain on' -> 'large-strain off'".format(PUSHOVER_SCRIPT, n_on))
    return src, over


# ---------------- model preparation ----------------
FISH_LIST = r"""
fish define hft_kind(cx)
    local b1 = block.contact.b1(cx)
    local b2 = block.contact.b2(cx)
    local k = 'other'
    if block.isgroup(b1, 'Masonry') and block.isgroup(b2, 'Masonry') then
        k = 'masonry'
    endif
    if block.isgroup(b1, 'Joist') or block.isgroup(b2, 'Joist') then
        k = 'pocket'
    endif
    if block.isgroup(b1, 'T_B') or block.isgroup(b2, 'T_B') then
        k = 'top'
    endif
    if block.isgroup(b1, 'S') or block.isgroup(b2, 'S') then
        k = 'base'
    endif
    hft_kind = k
end
fish define hft_list
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'HingeFT0', 'hingeft0') then
            local ys = 0.0
            local n = 0
            loop foreach local sc block.contact.subcontactlist(cx)
                ys = ys + comp.y(block.subcontact.pos(sc))
                n = n + 1
            endloop
            if n > 0 then
                ys = ys / n
            endif
            io.out('@@H,' + string(block.contact.id(cx)) + ',' + string(ys) + ',' + hft_kind(cx) + ',' + string(n) + ',' + string(PROP_EXPR))
        endif
    endloop
end
fish define hft_rest
    local tmin = 1e20
    local tmax = -1e20
    local n = 0
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'HingeFT0', 'hingeft0') then
            n = n
        else
            if hft_kind(cx) = 'masonry' then
                local t = PROP_EXPR
                if t < tmin then
                    tmin = t
                endif
                if t > tmax then
                    tmax = t
                endif
                n = n + 1
            endif
        endif
    endloop
    io.out('@@R,' + string(n) + ',' + string(tmin) + ',' + string(tmax))
end
"""


def prepare():
    os.makedirs(PO_DIR, exist_ok=True)
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(cmdpath(SAVE_IN)))
    it.command("python-reset-state false")
    ft = float(it.fish.get("tension"))
    log("{}: joint tension {:.3f} MPa".format(SAVE_IN, ft / 1e6))
    if SET_REST:
        props = " ".join("{} {:g}".format(k, v) for k, v in REST_PROPS.items())
        for rg in REST_RANGES:
            it.command("block contact property {} range {}".format(props, rg))
        it.command("[global tension = {:g}]".format(REST_PROPS["tension"]))
        it.command("[global cohesion = {:g}]".format(REST_PROPS["cohesion"]))
        log("  SET_REST: {} on masonry, base, top and pocket contacts (joist supports untouched)".format(props))
        ft = REST_PROPS["tension"]
    if not (FT_RANGE[0] <= ft <= FT_RANGE[1]):
        raise RuntimeError("wrong save: joint tension {:.3f} MPa, expected {:.2f}-{:.2f} MPa".format(
            ft / 1e6, FT_RANGE[0] / 1e6, FT_RANGE[1] / 1e6))

    # 1. group the strips (also without SET_STRIPS: the same ranges the .dat build used, for the check)
    it.command("block contact group 'none' slot '{}'".format(SLOT))
    for lab, lo, hi in BANDS:
        it.command("block contact group '{}' slot '{}' range pos-y {} {}".format(GRP, SLOT, lo, hi))
        log("  strip {:<20s} y {:+.3f} ... {:+.3f} m".format(lab, lo, hi))
    # 2. never the joist supports (Mohr, tension 1 Pa / cohesion 1 Pa by design)
    for sup in ("Joist_S1", "Joist_S2"):
        it.command("block contact group 'none' slot '{}' range group-intersection 'Joist' '{}'".format(SLOT, sup))
    if POCKETS == "all":
        it.command("block contact group '{}' slot '{}' range group-intersection 'Joist' 'Masonry'".format(GRP, SLOT))
    elif POCKETS == "build":
        it.command("block contact group 'none' slot '{}' range group-intersection 'Joist' 'Masonry'".format(SLOT))
    if SET_STRIPS:
        log("  joist pockets: {}".format(POCKETS))
        # 3. the f_t = 0 joint on the group
        props = " ".join("{} {:g}".format(k, v) for k, v in HINGE_PROPS.items())
        it.command("block contact property {} range group '{}' slot '{}'".format(props, GRP, SLOT))
        log("  block contact property {}  (group {})".format(props, GRP))
    else:
        log("  strips taken from the build ({}); properties not changed here".format(BUILD_DAT))

    # 4. list what was changed, read the tension back
    out = os.path.join(PO_DIR, "_hinge_ft0_list.log")
    if os.path.exists(out):
        os.remove(out)
    prop_ok = True
    for expr in ("block.contact.prop(cx, 'tension')", "-1.0"):
        try:
            it.command(FISH_LIST.replace("PROP_EXPR", expr))
            prop_ok = expr != "-1.0"
            break
        except Exception as e:
            log("  contact list with '{}' failed ({}); simpler variant".format(expr, e))
    it.command('program log-file "{}"'.format(cmdpath(out)))
    it.command("program log on")
    it.command("[hft_list]")
    if prop_ok:
        it.command("[hft_rest]")
    it.command("program log off")
    rows, rest = [], None
    with open(out) as f:
        for line in f:
            if line.startswith("@@H,"):
                p = line.strip().split(",")
                rows.append(dict(id=p[1], y=float(p[2]), kind=p[3], n_sub=int(p[4]), tension_Pa=float(p[5])))
            elif line.startswith("@@R,"):
                p = line.strip().split(",")
                rest = (int(p[1]), float(p[2]), float(p[3]))
    with open(os.path.join(PO_DIR, "hinge_ft0_contacts.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "y", "kind", "n_sub", "tension_Pa"])
        w.writeheader()
        w.writerows(rows)
    if not rows:
        raise RuntimeError("no contact in group {}: the pos-y ranges matched nothing".format(GRP))
    summ = {}
    for r in rows:
        key = (round(r["y"], 2), r["kind"])
        summ.setdefault(key, [0, 0.0])
        summ[key][0] += 1
        summ[key][1] = max(summ[key][1], r["tension_Pa"])
    log("  {} contacts set to f_t = 0 (by y, kind: count, max tension read back):".format(len(rows)))
    for (y, k), (n, t) in sorted(summ.items()):
        log("    y {:+.2f} m  {:<8s} {:4d}   {}".format(y, k, n, "{:.0f} Pa".format(t) if prop_ok else "(not read)"))
    bad = [r for r in rows if prop_ok and abs(r["tension_Pa"] - HINGE_PROPS["tension"]) > 1.0]
    if bad:
        raise RuntimeError("{} contacts did not take tension = {} (e.g. id {})".format(len(bad), HINGE_PROPS["tension"], bad[0]["id"]))
    others = [k for k in summ if k[1] not in ("masonry", "pocket", "top", "base")]
    if others:
        log("  ! contacts of an unexpected kind in the group: {}".format(others))
    if rest:
        log("  all other masonry contacts: {} with tension {:.3f} ... {:.3f} MPa (expected {:.3f})".format(
            rest[0], rest[1] / 1e6, rest[2] / 1e6, REST_FT_PA / 1e6))
        if (SET_REST or not SET_STRIPS) and (abs(rest[1] - REST_FT_PA) > 1e3 or abs(rest[2] - REST_FT_PA) > 1e3):
            raise RuntimeError("the masonry outside the strips does not have tension {:.3f} MPa".format(REST_FT_PA / 1e6))

    # 5. settle and save (only when the strips were changed here)
    if SET_STRIPS:
        log("  {}".format(SETTLE))
        it.command(SETTLE)
        it.command("model save '{}'".format(cmdpath(SAVE_OUT)))
        log("  saved {}.sav".format(SAVE_OUT))


def main():
    log("=" * 70)
    log("run_pushover_hinge_ft0  mode {}".format(MODE))
    if not os.path.isfile(PUSHOVER_SCRIPT):
        raise RuntimeError("{} not found in {}".format(PUSHOVER_SCRIPT, os.getcwd()))
    for d in DIRS:
        src, over = pushover_source(d)
        compile(src, PUSHOVER_SCRIPT, "exec")
        log("pushover {}: {} -- patched source compiles".format("POS" if d > 0 else "NEG",
                                                             ", ".join("{}={!r}".format(k, v) for k, v in over.items())))
    if MODE == "dry":
        log("dry run: no model work. Set MODE = 'prep' (strip setup only) or 'run'.")
        return
    if it is None:
        raise RuntimeError("itasca not available: run inside 3DEC")
    if not (os.path.isfile(SAVE_IN) or os.path.isfile(SAVE_IN + ".sav")):
        if BUILD_DAT and not SET_STRIPS:
            if not os.path.isfile(BUILD_DAT):
                raise RuntimeError("missing {}.sav and {}".format(SAVE_IN, BUILD_DAT))
            log("building {}.sav: program call '{}'".format(SAVE_IN, BUILD_DAT))
            t0 = time.time()
            it.command("python-reset-state false")
            it.command("program call '{}'".format(BUILD_DAT))
            it.command("python-reset-state false")
            log("  built in {:.1f} min".format((time.time() - t0) / 60))
        if not (os.path.isfile(SAVE_IN) or os.path.isfile(SAVE_IN + ".sav")):
            raise RuntimeError("missing {}.sav".format(SAVE_IN))
    prepare()
    if MODE == "prep":
        log("prep done: check {} and the log, then MODE = 'run'".format(os.path.join(PO_DIR, "hinge_ft0_contacts.csv")))
        return
    for d in DIRS:
        tag = "pos" if d > 0 else "neg"
        old = os.path.join(PO_DIR, "pushover_{}.csv".format(tag))
        if os.path.isfile(old):
            os.rename(old, old.replace(".csv", "_prev_{}.csv".format(time.strftime("%Y%m%d_%H%M%S"))))
        log("-" * 70)
        log("pushover {} from {}".format(tag.upper(), push_save()))
        src, _ = pushover_source(d)
        t0 = time.time()
        try:
            exec(compile(src, PUSHOVER_SCRIPT, "exec"), {"__name__": "__main__", "__file__": os.path.abspath(PUSHOVER_SCRIPT)})
        except SystemExit:
            pass
        except Exception as e:
            log("  pushover {} FAILED: {}".format(tag, e))
            continue
        log("  pushover {} done in {:.1f} min -> {}".format(tag, (time.time() - t0) / 60, old))
    log("finished")


main()
