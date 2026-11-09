# -*- coding: ascii -*-
"""
deformed_shape_stratF_zup.py -- residual deformed shape and joint-damage map from the Strategy F run saves,
Z-UP frame (3DEC 9).  2 Nov 2026.

PURPOSE
  Find where the f_t 0.05 MPa pk2 wall cracks under the recorded table motions, and in which run, so the hinge
  joint(s) of a hinged build (f_t = 0 at the base, the hinge and the top joint, 0.05 MPa elsewhere) can be picked
  from the dynamic response rather than from a pushover.

FRAME (Geo_Prep_zup*.dat): X = 1.292 - x in-plane, Y out-of-plane (+Y towards the joists), Z vertical.

WHAT IT DOES
  1. Restores the build save (REF_SAV) and records every block centroid and every masonry subcontact with its
     joint properties (tension, cohesion, stiffness-normal): the undamaged reference.
  2. Restores stratC_run_NN.sav of every run in RUNS (4 s after the record ends: the wall is AT REST) and records
       - every block centroid -> residual displacement; the OOP displacement is taken RELATIVE TO THE BASE SLAB 'S'
         (the table drive leaves the whole model translated);
       - every subcontact of masonry-masonry joints (MasJ), pockets (Masonry-Joist), the top joint (Masonry-T_B)
         and the base joint (Masonry-S): position, normal displacement, normal force, area, contact normal, block ids
         and the joint properties.
  3. Post-processes each run, and all runs together:
       - course profile: residual OOP displacement of each course (relative to S) against its original height,
         segment rotation and kink;
       - bed-joint map: residual largest opening and open fraction, and the DAMAGED area fraction -- subcontacts
         whose tension property has dropped below DAMAGED_RATIO x the build value (MASON softening). Contacts created
         during the run (tension 0 AND cohesion 0, the material-table default) are counted apart as "new";
       - damage history: for every bed joint, the run in which it first reaches 10 % and 50 % damaged area.
     A bed joint is flagged as a HINGE when its damaged fraction >= HINGE_DMG_FRAC or its residual opening > HINGE_MM.

WHY DAMAGE AND NOT ONLY OPENING
  A cracked bed joint under the precompression recloses once the shaking stops, so the residual opening understates
  where the wall hinged. The softened tension stays. The script checks on the first save whether the properties do
  change from the build values; if no subcontact ever differs, it says so and the hinge reading falls back to the
  residual openings and course kinks. EXTRA_PROPS can name further state variables of the joint model to dump
  (each is probed on the reference save first; a name the model does not know is dropped with a message).

OUTPUT  <OUT_DIR>/run_NN_courses.csv  course profile of run NN
        <OUT_DIR>/run_NN_joints.csv   joint map of run NN (bed joints, base, top joint, pockets per floor)
        <OUT_DIR>/joints_by_run.csv   every joint of every run in one long table
        <OUT_DIR>/damage_onset.csv    first run with >= 10 % / >= 50 % damaged area, per bed joint
        <OUT_DIR>/summary.txt         the printed tables
        <OUT_DIR>/deformed_shape_stratF.png
        <OUT_DIR>/_dump_<tag>.log     raw dumps (python deformed_shape_stratF_zup.py --post re-processes them)

USAGE   inside 3DEC (no other model open):  program call 'deformed_shape_stratF_zup.py'
        outside 3DEC:                       python deformed_shape_stratF_zup.py --post [--run-dir <folder>]
        another run folder inside 3DEC:     a wrapper setting DFM_RUN_DIR, e.g. run_deformed_shape_1p5.py
        A run still in progress is fine: only the saves that exist are read. Re-running skips dumps already written
        (delete a _dump_run_NN.log to re-dump that run, or set REDUMP = True).
"""

import os, sys, math, csv, re

# ============================ CONFIG =================================
RUN_DIR  = "stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr"                 # the Strategy F output folder (default)
# override without editing (9 Nov 2026): a wrapper that sets DFM_RUN_DIR and execs this file, or, outside 3DEC,
#   python deformed_shape_stratF_zup.py --post --run-dir <folder>
RUN_DIR  = globals().get("DFM_RUN_DIR", RUN_DIR)
if "--run-dir" in sys.argv[:-1]:
    RUN_DIR = sys.argv[sys.argv.index("--run-dir") + 1]
REF_SAV  = "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05"      # the save the run started from
RUNS     = None            # None = every stratC_run_NN.sav found in RUN_DIR; or a list, e.g. [12, 17, 18, 19, 20, 21]
OUT_DIR  = os.path.join(RUN_DIR, "deformed_shape")
DMG      = "tension"     # set in main() from DAMAGE_PROP
REDUMP   = False           # True: dump again even when _dump_<tag>.log exists
EXTRA_PROPS = ["stiffness-initial"]   # further joint-model properties / state variables to dump, probed first

T_WALL        = 0.21       # m
SENSOR_Z      = 2.06       # m, Top_Quarter_A/B height
FLOOR_Z       = (0.32250, 2.19886)    # m, joist centres of the pk2 build (drawn on the figure)
FLOOR_SPLIT_Z = 1.26       # m, pockets below = lower floor
OPEN_MM       = 0.01       # a subcontact counts as open above this residual opening
HINGE_MM      = 0.5        # residual opening that flags a hinge on its own
DAMAGE_PROP   = "auto"     # property that carries the joint damage: "tension", "stiffness-normal", or "auto" (tension
                           # if any bed joint's tension ever changes, else stiffness-normal -- MASON degrades kn, 2 Nov)
DAMAGED_RATIO = 0.5        # damaged subcontact: DAMAGE_PROP < this x its build value
HINGE_DMG_FRAC = 0.5       # damaged area fraction that flags a hinge
COURSE_GAP    = 0.02       # m, centroid-height gap that separates two courses / two joints
BED_NZ_MIN    = 0.8        # |normal_z| above this = bed joint

# ============================ FISH ===================================
FISH_BLOCKS = r"""
fish define dfm_blocks
    loop foreach local b block.list
        local g = 'other'
        if block.isgroup(b, 'Masonry') then
            g = 'Masonry'
        else
            if block.isgroup(b, 'T_B') then
                g = 'T_B'
            else
                if block.isgroup(b, 'Joist') then
                    g = 'Joist'
                else
                    if block.isgroup(b, 'S') then
                        g = 'S'
                    endif
                endif
            endif
        endif
        if g # 'other' then
            local p = block.pos(b)
            io.out('@@B,' + g + ',' + string(block.id(b)) + ',' + string(comp.x(p)) + ',' + string(comp.y(p)) + ',' + string(comp.z(p)) + ',' + string(block.mass(b)))
        endif
    endloop
end
"""
FISH_GLOBALS = r"""
fish define dfm_globals
    io.out('@@G,' + string(tension) + ',' + string(cohesion) + ',' + string(kn))
end
"""

# contact classes by elimination, as the v5 pushover driver: everything not named below is masonry-masonry
GROUP_CMDS = """
block contact group "MasJ"   slot "dfm"
block contact group "Base"   slot "dfm" range group-intersection 'S' 'Masonry'
block contact group "TopTB"  slot "dfm" range group-intersection 'T_B' 'Masonry'
block contact group "Pocket" slot "dfm" range group-intersection 'Joist' 'Masonry'
block contact group "Skip"   slot "dfm" range group-intersection 'Joist' 'Joist_S1'
block contact group "Skip"   slot "dfm" range group-intersection 'Joist' 'Joist_S2'
"""

FISH_CONTACTS = r"""
fish define dfm_cgroup(cx)
    local g = ''
    if block.contact.isgroup(cx, 'MasJ', 'dfm') then
        g = 'MasJ'
    endif
    if block.contact.isgroup(cx, 'Pocket', 'dfm') then
        g = 'Pocket'
    endif
    if block.contact.isgroup(cx, 'TopTB', 'dfm') then
        g = 'TopTB'
    endif
    if block.contact.isgroup(cx, 'Base', 'dfm') then
        g = 'Base'
    endif
    dfm_cgroup = g
end
fish define dfm_contacts
    loop foreach local cx block.contact.list
        local g = dfm_cgroup(cx)
        if g # '' then
            local nz = NORMAL_EXPR
            loop foreach local sc block.contact.subcontactlist(cx)
                local p = block.subcontact.pos(sc)
                io.out('@@C,' + g + ',' + string(block.contact.id(cx)) + ',' + string(nz) + ',' + string(comp.x(p)) + ',' + string(comp.y(p)) + ',' + string(comp.z(p)) + ',' + string(block.subcontact.disp.norm(sc)) + ',' + string(block.subcontact.force.norm(sc)) + ',' + string(block.subcontact.area(sc)) BIDS_EXPR PROPS_EXPR)
            endloop
        endif
    endloop
end
"""

FISH_PROBE = r"""
fish define dfm_probe
    dfm_probe = -1.0
    local done = 0
    loop foreach local cx block.contact.list
        if done = 0 then
            if block.contact.isgroup(cx, 'MasJ', 'dfm') then
                loop foreach local sc block.contact.subcontactlist(cx)
                    if done = 0 then
                        dfm_probe = EXPR
                        done = 1
                    endif
                endloop
            endif
        endif
    endloop
end
"""

BASE_PROPS = ["tension", "cohesion", "stiffness-normal"]


def _cmdpath(p):
    return p.replace("\\", "/")


def _restore(it, save):
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(_cmdpath(save[:-4] if save.endswith(".sav") else save)))
    it.command("python-reset-state false")
    it.command(GROUP_CMDS)
    it.command(FISH_BLOCKS)


def probe_props(it):
    """joint properties this build answers for (probed on one masonry-masonry subcontact)"""
    ok = []
    for name in BASE_PROPS + [p for p in EXTRA_PROPS if p not in BASE_PROPS]:
        try:
            it.command(FISH_PROBE.replace("EXPR", "block.subcontact.prop(sc, '{}')".format(name)))
            it.command("[dfm_probe_val = dfm_probe]")
            v = float(it.fish.get("dfm_probe_val"))
            ok.append(name)
            print("  property '{}' readable (first masonry subcontact: {:.6g})".format(name, v))
        except Exception as e:
            print("  property '{}' not readable here ({}); dropped".format(name, str(e).splitlines()[0][:80]))
    return ok


def run_dump(it, tag, save, props, restored=False):
    log = os.path.join(OUT_DIR, "_dump_{}.log".format(tag))
    if os.path.exists(log):
        os.remove(log)
    if not restored:
        _restore(it, save)
    pexpr = "".join(" + ',' + string(block.subcontact.prop(sc, '{}'))".format(p) for p in props)
    bids = "+ ',' + string(block.id(block.contact.b1(cx))) + ',' + string(block.id(block.contact.b2(cx)))"
    variants = [("comp.z(block.contact.normal(cx))", bids, pexpr), ("comp.z(block.contact.normal(cx))", "+ ',,'", pexpr),
                ("-9.0", "+ ',,'", pexpr), ("-9.0", "+ ',,'", "")]
    done = None
    for nexpr, bexpr, pe in variants:
        try:
            it.command(FISH_CONTACTS.replace("NORMAL_EXPR", nexpr).replace("BIDS_EXPR", bexpr).replace("PROPS_EXPR", pe))
            done = (nexpr, bexpr, pe)
            break
        except Exception as e:
            print("contact dump variant failed ({}); trying a simpler one".format(str(e).splitlines()[0][:100]))
    if done is None:
        raise RuntimeError("could not define the contact dump")
    used_props = props if done[2] else []
    have_g = True
    try:
        it.command(FISH_GLOBALS)
    except Exception:
        have_g = False
    it.command('program log-file "{}"'.format(_cmdpath(log)))
    it.command("program log on")
    it.command("[io.out('@@P,' + '{}')]".format("|".join(used_props)))
    if have_g:
        try:
            it.command("[dfm_globals]")
        except Exception:
            have_g = False
    it.command("[dfm_blocks]")
    it.command("[dfm_contacts]")
    it.command("program log off")
    print("dumped {} -> {}  (normals {}, block ids {}, properties {}, build globals {})".format(
        save, log, "yes" if done[0] != "-9.0" else "NO", "yes" if done[1] == bids else "NO",
        ",".join(used_props) or "NONE", "yes" if have_g else "NO"))
    return log


# ============================ PARSING ================================
def read_dump(log):
    blocks, subs, props, glob = {}, [], [], None
    with open(log, errors="replace") as f:
        for line in f:
            i = line.find("@@P,")
            if i >= 0:
                s = line[i + 4:].strip()
                props = [p for p in s.split("|") if p]
                continue
            i = line.find("@@G,")
            if i >= 0:
                s = line[i + 4:].strip().split(",")
                try:
                    glob = dict(tension=float(s[0]), cohesion=float(s[1]), kn=float(s[2]))
                except Exception:
                    pass
                continue
            i = line.find("@@B,")
            if i >= 0:
                s = line[i + 4:].strip().split(",")
                blocks[int(float(s[1]))] = (s[0], float(s[2]), float(s[3]), float(s[4]), float(s[5]))
                continue
            i = line.find("@@C,")
            if i >= 0:
                s = line[i + 4:].strip().split(",")
                d = dict(g=s[0], cid=int(float(s[1])), nz=float(s[2]), x=float(s[3]), y=float(s[4]),
                         z=float(s[5]), un=float(s[6]), fn=float(s[7]), area=float(s[8]), b1=None, b2=None)
                if len(s) >= 11 and s[9].strip() and s[10].strip():
                    d["b1"], d["b2"] = int(float(s[9])), int(float(s[10]))
                for k, p in enumerate(props):
                    try:
                        d[p] = float(s[11 + k])
                    except (IndexError, ValueError):
                        d[p] = float("nan")
                subs.append(d)
    return blocks, subs, props, glob


def cluster(values, gap):
    order = sorted(range(len(values)), key=lambda k: values[k])
    groups, cur = [], []
    for k in order:
        if cur and values[k] - values[cur[-1]] > gap:
            groups.append(cur); cur = []
        cur.append(k)
    if cur:
        groups.append(cur)
    return groups


def interp(x, xs, ys):
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


# ============================ POST ===================================
def build_values(ref_subs, glob):
    """build tension / cohesion / kn: the most common value on the reference save's masonry-masonry subcontacts
    (what the joints actually carry); the save's FISH globals only when the dump has no properties"""
    if not ref_subs or "tension" not in ref_subs[0]:
        return glob or dict(tension=float("nan"), cohesion=float("nan"), kn=float("nan"))
    out = {}
    for key, p in (("tension", "tension"), ("cohesion", "cohesion"), ("kn", "stiffness-normal")):
        vals = [round(s[p], 3) for s in ref_subs if s["g"] == "MasJ" and p in s and s[p] == s[p]]
        out[key] = max(set(vals), key=vals.count) if vals else float("nan")
    return out


def reference_heights(subs, mas, ref_blocks, blocks, crs_rot):
    """undeformed height of each subcontact: halfway between its two masonry blocks' original centroids when the dump
    has block ids; one masonry block (pocket, base, top joint): the offset along that block's rotated vertical axis
    added to its original height; no ids: nearest-block rules (see deformed_shape_3dec.py)"""
    course_of = {}
    for gi, grp in enumerate(cluster([r["z0"] for r in mas], COURSE_GAP)):
        for k in grp:
            course_of[mas[k]["id"]] = gi
    mas_ids = [r["id"] for r in mas]
    cur = {b: blocks[b] for b in mas_ids if b in blocks}

    def local(bid, s):
        _, _, yc, zc, _ = cur[bid]
        th = crs_rot[course_of[bid]]
        return ref_blocks[bid][3] + (s["z"] - zc) * math.cos(th) + (s["y"] - yc) * math.sin(th)

    d2 = lambda b, s: (cur[b][1] - s["x"]) ** 2 + (cur[b][2] - s["y"]) ** 2 + (cur[b][3] - s["z"]) ** 2
    by_course = {}
    for b in mas_ids:
        by_course.setdefault(course_of[b], []).append(b)
    nearest = lambda s, pool=None: min(pool if pool is not None else mas_ids, key=lambda b: d2(b, s))
    out = []
    for s in subs:
        ids = [b for b in (s.get("b1"), s.get("b2")) if b is not None and b in cur]
        if len(ids) == 2:
            out.append(0.5 * (ref_blocks[ids[0]][3] + ref_blocks[ids[1]][3]))
        elif len(ids) == 1:
            out.append(local(ids[0], s))
        elif s["g"] == "MasJ" and not (-2 < s["nz"] and abs(s["nz"]) < BED_NZ_MIN):
            a = nearest(s); ca = course_of[a]
            cand = [nearest(s, by_course[c]) for c in (ca - 1, ca + 1) if c in by_course]
            b = min(cand, key=lambda q: d2(q, s)) if cand else a
            out.append(0.5 * (ref_blocks[a][3] + ref_blocks[b][3]))
        else:
            out.append(local(nearest(s), s))
    return out


def post(tag, run_no, ref_blocks, bv, blocks, subs, P):
    # ---- blocks: residual displacement, OOP relative to the base slab S ----
    s_ids = [b for b, v in blocks.items() if v[0] == "S" and b in ref_blocks]
    dS = (0.0, 0.0, 0.0)
    if s_ids:
        b = s_ids[0]
        dS = tuple(blocks[b][k] - ref_blocks[b][k] for k in (1, 2, 3))
    rows = []
    for bid, (g, x, y, z, m) in blocks.items():
        if bid not in ref_blocks:
            continue
        _, x0, y0, z0, _ = ref_blocks[bid]
        rows.append(dict(id=bid, group=g, x0=x0, y0=y0, z0=z0, mass=m, dx_mm=(x - x0 - dS[0]) * 1e3,
                         dy_mm=(y - y0 - dS[1]) * 1e3, dz_mm=(z - z0 - dS[2]) * 1e3))
    mas = [r for r in rows if r["group"] == "Masonry"]
    if not mas:
        P("{}: no Masonry blocks matched the reference -- block ids differ?".format(tag)); return None
    crs = []
    for grp in cluster([r["z0"] for r in mas], COURSE_GAP):
        bl = [mas[k] for k in grp]; M = sum(b["mass"] for b in bl)
        wm = lambda key: sum(b[key] * b["mass"] for b in bl) / M
        crs.append(dict(z0=wm("z0"), n=len(bl), dy_mm=wm("dy_mm"), dy_min=min(b["dy_mm"] for b in bl),
                        dy_max=max(b["dy_mm"] for b in bl), dz_mm=wm("dz_mm")))
    for i, c in enumerate(crs):
        c["rot_mrad"] = float("nan") if i == 0 else (c["dy_mm"] - crs[i - 1]["dy_mm"]) / (c["z0"] - crs[i - 1]["z0"])
    for i, c in enumerate(crs):
        c["kink_mrad"] = (crs[i + 1]["rot_mrad"] - c["rot_mrad"]) if 0 < i < len(crs) - 1 else float("nan")
    with open(os.path.join(OUT_DIR, tag + "_courses.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["z0", "n", "dy_mm", "dy_min", "dy_max", "dz_mm", "rot_mrad", "kink_mrad"])
        w.writeheader(); [w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in c.items()}) for c in crs]
    d_sens = interp(SENSOR_Z, [c["z0"] for c in crs], [c["dy_mm"] for c in crs])
    P("=" * 96)
    P("{}: {} masonry blocks in {} courses; base slab S moved {:+.2f} mm OOP (removed below)".format(
        tag, len(mas), len(crs), dS[1] * 1e3))
    P("residual OOP displacement at {:.2f} m, relative to S: {:+.2f} mm".format(SENSOR_Z, d_sens))
    tb = [r for r in rows if r["group"] == "T_B"]
    if tb:
        P("top beam relative to S: dy {:+.2f} mm, dz {:+.2f} mm".format(tb[0]["dy_mm"], tb[0]["dz_mm"]))
    crs_rot = []
    for i in range(len(crs)):
        v = [crs[j]["rot_mrad"] for j in (i, i + 1) if 0 <= j < len(crs) and crs[j]["rot_mrad"] == crs[j]["rot_mrad"]]
        crs_rot.append(math.atan(sum(v) / len(v) / 1e3) if v else 0.0)
    for s_, zr in zip(subs, reference_heights(subs, mas, ref_blocks, blocks, crs_rot)):
        s_["zref"] = zr
    # ---- joints ----
    ft = bv.get("tension", float("nan")); kn = bv.get("kn", float("nan"))
    has_t = any("tension" in s for s in subs)
    dref = ft if DMG == "tension" else kn

    def joint_stats(sub, label):
        A = sum(s["area"] for s in sub) or 1.0
        yc = sum(s["y"] * s["area"] for s in sub) / A
        zc = sum(s["zref"] * s["area"] for s in sub) / A
        smax = max(sub, key=lambda s: s["un"])
        op = [s for s in sub if s["un"] * 1e3 > OPEN_MM]
        Aop = sum(s["area"] for s in op)
        yw = (sum((s["y"] - yc) * s["un"] * s["area"] for s in op) / sum(s["un"] * s["area"] for s in op)) if op else 0.0
        r = dict(run=run_no, kind=label, z0=zc, n_sub=len(sub), max_open_mm=smax["un"] * 1e3, open_frac=Aop / A,
                 open_face=("" if not op else ("+Y" if yw > 0 else "-Y")), N_kN=-sum(s["fn"] for s in sub) / 1e3)
        if has_t and ft == ft:
            iszero = lambda s: s.get("tension", 1) == 0.0 and s.get("cohesion", 1) == 0.0 and ft > 0
            zero = [s for s in sub if iszero(s)]
            old = [s for s in sub if not iszero(s)]
            Ao = sum(s["area"] for s in old) or 1.0
            if dref == dref and dref > 0 and old and DMG in old[0]:
                dmg = [s for s in old if s[DMG] < DAMAGED_RATIO * dref]
                r["dmg_frac"] = sum(s["area"] for s in dmg) / Ao
            else:
                r["dmg_frac"] = float("nan")
            r["t_ratio_mean"] = (sum(s["tension"] * s["area"] for s in old) / Ao / ft) if (ft > 0 and old) else float("nan")
            r["n_zero"] = len(zero)
            if "stiffness-normal" in sub[0] and kn == kn and kn > 0 and old:
                r["kn_ratio_min"] = min(s["stiffness-normal"] for s in old) / kn
                r["kn_ratio_mean"] = sum(s["stiffness-normal"] * s["area"] for s in old) / Ao / kn
        return r

    mj = [s for s in subs if s["g"] == "MasJ"]
    have_n = any(s["nz"] > -2 for s in mj)
    bed = [s for s in mj if (abs(s["nz"]) >= BED_NZ_MIN if have_n else True)]
    head = [s for s in mj if have_n and abs(s["nz"]) < BED_NZ_MIN]
    jrows = []
    for grp in cluster([s["zref"] for s in bed], COURSE_GAP):
        sub = [bed[k] for k in grp]
        if len(sub) >= 4:
            jrows.append(joint_stats(sub, "bed"))
    for g, label, sel in (("Base", "base", None), ("TopTB", "top_joint", None),
                          ("Pocket", "pocket_lower", lambda s: s["zref"] < FLOOR_SPLIT_Z),
                          ("Pocket", "pocket_upper", lambda s: s["zref"] >= FLOOR_SPLIT_Z)):
        sub = [s for s in subs if s["g"] == g and (sel is None or sel(s))]
        if sub:
            jrows.append(joint_stats(sub, label))
    if head:
        jrows.append(joint_stats(head, "head_all"))
    jrows.sort(key=lambda r: r["z0"])
    for r in jrows:
        r["hinge"] = int(r["kind"] in ("bed", "base", "top_joint") and (
            r["max_open_mm"] > HINGE_MM or r.get("dmg_frac", 0) >= HINGE_DMG_FRAC))
    with open(os.path.join(OUT_DIR, tag + "_joints.csv"), "w", newline="") as f:
        keys = sorted({k for r in jrows for k in r}, key=lambda k: list(JKEYS).index(k) if k in JKEYS else 99)
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        [w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()}) for r in jrows]
    P("")
    P("  joint         z0 (m)  res. open (mm)  open frac  face   damaged  kn_mean/kn  kn_min/kn  zero  N (kN)")
    for r in jrows:
        P("  {:12s} {:6.3f} {:10.3f}      {:6.2f}     {:3s}   {:6s}   {:6s}     {:6s}   {:4s} {:8.2f}{}".format(
            r["kind"], r["z0"], r["max_open_mm"], r["open_frac"], r["open_face"],
            _f(r.get("dmg_frac"), "{:.2f}"), _f(r.get("kn_ratio_mean"), "{:.2f}"), _f(r.get("kn_ratio_min"), "{:.2f}"),
            _f(r.get("n_zero"), "{:d}"), r["N_kN"], "  <- HINGE" if r["hinge"] else ""))
    hinges = [r for r in jrows if r["hinge"]]
    P("HINGES ({}): ".format(tag) + (", ".join("{:.3f} m ({})".format(r["z0"], r["kind"]) for r in hinges) or "none"))
    return dict(tag=tag, run=run_no, crs=crs, joints=jrows, d_sens=d_sens)


JKEYS = ("run", "kind", "z0", "n_sub", "max_open_mm", "open_frac", "open_face", "dmg_frac", "t_ratio_mean",
         "kn_ratio_mean", "kn_ratio_min", "n_zero", "N_kN", "hinge")


def _f(v, fmt):
    if v is None or (isinstance(v, float) and v != v):
        return "  -"
    return fmt.format(int(v) if fmt == "{:d}" else v)


def damage_onset(results, P):
    """per bed joint (matched across runs by its original height): first run with >= 10 % and >= 50 % damaged area"""
    bed = [(r["run"], j) for r in results for j in r["joints"] if j["kind"] in ("bed", "base", "top_joint")]
    if not bed or not any("dmg_frac" in j for _, j in bed):
        return []
    groups = cluster([j["z0"] for _, j in bed], COURSE_GAP / 2)
    rows = []
    for grp in groups:
        seq = sorted((bed[k] for k in grp), key=lambda t: t[0])
        z = sum(j["z0"] for _, j in seq) / len(seq)
        first = lambda thr: next((rn for rn, j in seq if j.get("dmg_frac", 0) >= thr), "")
        rows.append(dict(kind=seq[-1][1]["kind"], z0=round(z, 4), first_run_10pct=first(0.10),
                         first_run_50pct=first(0.50), dmg_frac_last=round(seq[-1][1].get("dmg_frac", float("nan")), 3),
                         res_open_last_mm=round(seq[-1][1]["max_open_mm"], 3)))
    with open(os.path.join(OUT_DIR, "damage_onset.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); [w.writerow(r) for r in rows]
    hit = [r for r in rows if r["first_run_10pct"] != ""]
    P("")
    P("=" * 96)
    P("DAMAGE ONSET (bed / base / top joints that reach >= 10 % damaged area; damaged = {} < {:.0%} of build)".format(
        DMG, DAMAGED_RATIO))
    P("  z0 (m)   kind        >=10 % in run   >=50 % in run   damaged, last run   res. open, last run (mm)")
    for r in sorted(hit, key=lambda r: r["z0"]):
        P("  {:6.3f}   {:10s}  {:>8}        {:>8}          {:6.2f}             {:7.3f}".format(
            r["z0"], r["kind"], r["first_run_10pct"], r["first_run_50pct"], r["dmg_frac_last"], r["res_open_last_mm"]))
    if not hit:
        P("  none: no bed joint has softened in the runs read")
    return rows


def zero_report(ref_subs, bv, P):
    """where the masonry subcontacts with tension 0 AND cohesion 0 sit in the BUILD save (material-table default:
    contacts that formed after the property assignment)"""
    if not ref_subs or "tension" not in ref_subs[0] or not (bv["tension"] > 0):
        return
    mj = [s for s in ref_subs if s["g"] in ("MasJ", "Pocket", "Base", "TopTB")]
    z = [s for s in mj if s["tension"] == 0.0 and s["cohesion"] == 0.0]
    P("")
    P("ZERO-STRENGTH SUBCONTACTS IN THE BUILD SAVE (tension 0 and cohesion 0): {} of {} ({:.1%} of the area)".format(
        len(z), len(mj), sum(s["area"] for s in z) / (sum(s["area"] for s in mj) or 1)))
    if not z:
        return
    by_g = {}
    for s in z:
        by_g[s["g"]] = by_g.get(s["g"], 0) + 1
    P("  by class: " + ", ".join("{} {}".format(k, v) for k, v in sorted(by_g.items())))
    bed = [s for s in z if s["nz"] > -2 and abs(s["nz"]) >= BED_NZ_MIN]
    head = [s for s in z if s["nz"] > -2 and abs(s["nz"]) < BED_NZ_MIN]
    P("  bed-joint orientation {}, head-joint orientation {}".format(len(bed), len(head)))
    cids = {}
    for s in mj:
        c = cids.setdefault(s["cid"], [0, 0]); c[0] += 1
    zc = {}
    for s in z:
        zc[s["cid"]] = zc.get(s["cid"], 0) + 1
    whole = sum(1 for c, n in zc.items() if n == cids[c][0])
    P("  in {} contacts: {} entirely zero, {} partly zero".format(len(zc), whole, len(zc) - whole))
    xs = sorted(s["x"] for s in z); ys = sorted(s["y"] for s in z)
    q = lambda v, f: v[min(len(v) - 1, int(f * len(v)))]
    P("  X (in-plane, 0-1.292 m): min {:.3f}  10% {:.3f}  median {:.3f}  90% {:.3f}  max {:.3f}".format(
        xs[0], q(xs, .1), q(xs, .5), q(xs, .9), xs[-1]))
    P("  Y (OOP, 0-{:.3f} m):     min {:.3f}  10% {:.3f}  median {:.3f}  90% {:.3f}  max {:.3f}".format(
        T_WALL, ys[0], q(ys, .1), q(ys, .5), q(ys, .9), ys[-1]))
    edge_x = sum(1 for s in z if s["x"] < 0.02 or s["x"] > 1.272)
    face_y = sum(1 for s in z if s["y"] < 0.01 or s["y"] > T_WALL - 0.01)
    mid_y = sum(1 for s in z if abs(s["y"] - T_WALL / 2) < 0.01)
    P("  within 20 mm of a wall end: {}; within 10 mm of a face: {}; within 10 mm of mid-thickness (leaf split): {}".format(
        edge_x, face_y, mid_y))
    hz = {}
    for s in z:
        k = round(s["z"] / 0.0586)
        hz[k] = hz.get(k, 0) + 1
    P("  per height band (Z, count): " + "  ".join("{:.2f}:{}".format(k * 0.0586, v) for k, v in sorted(hz.items())))
    with open(os.path.join(OUT_DIR, "zero_strength_subcontacts_ref.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["class", "contact_id", "x", "y", "z", "normal_z", "area", "b1", "b2"])
        [w.writerow([s["g"], s["cid"], round(s["x"], 4), round(s["y"], 4), round(s["z"], 4), round(s["nz"], 3),
                     s["area"], s["b1"], s["b2"]]) for s in z]


# ============================ PLOT ===================================
def plot(results, onset):
    try:
        import matplotlib
        matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print("matplotlib not available ({}); CSVs written, no figure".format(e)); return
    has_d = any("dmg_frac" in j for r in results for j in r["joints"])
    fig, ax = plt.subplots(1, 3, figsize=(15, 6.5), sharey=True, gridspec_kw=dict(width_ratios=[1, 1.4, 1]))
    n = len(results)
    cols = [plt.cm.viridis(0.9 - 0.85 * k / max(n - 1, 1)) for k in range(n)]
    show = set(range(n)) if n <= 8 else set([0] + list(range(max(0, n - 7), n)))
    for k, r in enumerate(results):
        if k in show:
            ax[0].plot([q["dy_mm"] for q in r["crs"]], [q["z0"] for q in r["crs"]], "o-", ms=2.5, lw=1.2, color=cols[k],
                       label="run {}".format(r["run"]))
    ax[0].axvline(0, color="grey", lw=0.6)
    ax[0].set_xlabel("residual OOP displacement rel. to base (mm)")
    ax[0].set_ylabel("original height Z (m)")
    ax[0].set_title("Residual deformed shape", loc="left", fontsize=10)
    ax[0].legend(fontsize=7, ncol=2)
    # damage map: joints x runs
    runs = [r["run"] for r in results]
    key = "dmg_frac" if has_d else "max_open_mm"
    for k, r in enumerate(results):
        for j in r["joints"]:
            if j["kind"] not in ("bed", "base", "top_joint"):
                continue
            v = j.get(key, float("nan"))
            if v != v:
                continue
            c = plt.cm.inferno_r(min(v, 1.0) if has_d else min(v / 2.0, 1.0))
            ax[1].add_patch(plt.Rectangle((r["run"] - 0.45, j["z0"] - 0.03), 0.9, 0.06, color=c, lw=0))
            if j["max_open_mm"] > HINGE_MM:
                ax[1].plot(r["run"], j["z0"], "o", mfc="none", mec="#2a78d6", ms=5)
    ax[1].set_xlim(min(runs) - 0.6, max(runs) + 0.6)
    ax[1].set_xticks(runs)
    ax[1].set_xlabel("run")
    sm = plt.cm.ScalarMappable(cmap="inferno_r", norm=plt.Normalize(0, 1 if has_d else 2.0))
    cb = fig.colorbar(sm, ax=ax[1], pad=0.01)
    cb.set_label("damaged area fraction ({} < {:.0%} build)".format(DMG, DAMAGED_RATIO) if has_d
                 else "residual opening (mm)")
    ax[1].set_title("Bed-joint damage by run (blue ring: residual opening > {} mm)".format(HINGE_MM), loc="left",
                    fontsize=10)
    # last run: damage and residual opening against height
    last = results[-1]
    bj = [j for j in last["joints"] if j["kind"] in ("bed", "base", "top_joint")]
    ax[2].plot([j["max_open_mm"] for j in bj], [j["z0"] for j in bj], "s-", ms=3, color="#2a78d6",
               label="residual opening (mm)")
    if has_d:
        ax2 = ax[2].twiny()
        ax2.plot([j.get("dmg_frac", float("nan")) for j in bj], [j["z0"] for j in bj], "o-", ms=3, color="#c2410c",
                 label="damaged fraction")
        ax2.set_xlim(0, 1.05); ax2.set_xlabel("damaged area fraction", color="#c2410c")
    ax[2].axvline(HINGE_MM, color="grey", ls=":", lw=0.8)
    ax[2].set_xlabel("residual largest opening (mm)", color="#2a78d6")
    ax[2].set_title("Run {}: joints against height".format(last["run"]), loc="left", fontsize=10, pad=28)
    for a_ in ax:
        a_.grid(alpha=0.25)
        for zf in FLOOR_Z:
            a_.axhline(zf, color="#a16207", ls="--", lw=0.8, alpha=0.7)
    ax[0].axhline(SENSOR_Z, color="grey", ls=":", lw=0.8)
    ax[0].text(0.02, FLOOR_Z[1] + 0.02, "upper floor", color="#a16207", fontsize=8, transform=ax[0].get_yaxis_transform())
    ax[0].text(0.02, FLOOR_Z[0] + 0.02, "lower floor", color="#a16207", fontsize=8, transform=ax[0].get_yaxis_transform())
    fig.suptitle("{}  --  Strategy F residual state (4 s after each record)".format(RUN_DIR), fontsize=10)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "deformed_shape_stratF.png")
    fig.savefig(out, dpi=150)
    print("-> " + out)


# ============================ MAIN ===================================
def find_runs():
    have = {}
    if os.path.isdir(RUN_DIR):
        for fn in os.listdir(RUN_DIR):
            m = re.match(r"stratC_run_(\d+)\.sav$", fn)
            if m:
                have[int(m.group(1))] = os.path.join(RUN_DIR, fn)
    sel = sorted(have) if RUNS is None else [r for r in RUNS if r in have]
    miss = [] if RUNS is None else [r for r in RUNS if r not in have]
    if miss:
        print("runs without a save (skipped): {}".format(miss))
    return [(r, have[r]) for r in sel]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    post_only = "--post" in sys.argv
    it = None
    if not post_only:
        try:
            import itasca as it
        except ImportError:
            print("not inside 3DEC: re-processing existing dumps (same as --post)"); post_only = True
    ref_log = os.path.join(OUT_DIR, "_dump_ref.log")
    logs = {}
    if not post_only:
        runs = find_runs()
        print("{} run saves in {}: {}".format(len(runs), RUN_DIR, [r for r, _ in runs]))
        if REDUMP or not os.path.isfile(ref_log):
            _restore(it, REF_SAV)
            props = probe_props(it)
            ref_log = run_dump(it, "ref", REF_SAV, props, restored=True)
        props = read_dump(ref_log)[2]
        for rn, sv in runs:
            tag = "run_{:02d}".format(rn)
            p = os.path.join(OUT_DIR, "_dump_{}.log".format(tag))
            logs[rn] = p if (os.path.isfile(p) and not REDUMP) else run_dump(it, tag, sv, props)
    else:
        for fn in sorted(os.listdir(OUT_DIR)):
            m = re.match(r"_dump_run_(\d+)\.log$", fn)
            if m and (RUNS is None or int(m.group(1)) in RUNS):
                logs[int(m.group(1))] = os.path.join(OUT_DIR, fn)
    ref_blocks, ref_subs, _, glob = read_dump(ref_log)
    bv = build_values(ref_subs, glob)
    lines = []
    P = lambda s="": (print(s), lines.append(s))
    P("reference {}: build tension {:.4g} Pa, cohesion {:.4g} Pa, kn {:.4g} Pa/m ({})".format(
        REF_SAV, bv["tension"], bv["cohesion"], bv["kn"], "most common value on the reference subcontacts"
        if (ref_subs and "tension" in ref_subs[0]) else "save FISH globals"))
    if glob:
        P("  save FISH globals: tension {:.4g}, cohesion {:.4g}, kn {:.4g}".format(glob["tension"], glob["cohesion"], glob["kn"]))
    if ref_subs and "tension" in ref_subs[0]:
        off = [s for s in ref_subs if s["g"] == "MasJ" and abs(s["tension"] - bv["tension"]) > 1e-3 * max(bv["tension"], 1)]
        P("reference check: {} of {} masonry-masonry subcontacts off the build tension".format(
            len(off), sum(1 for s in ref_subs if s["g"] == "MasJ")))
    global DMG
    DMG = DAMAGE_PROP
    if DMG == "auto":
        DMG = "stiffness-normal"
        if logs and bv["tension"] == bv["tension"] and bv["tension"] > 0:
            last = read_dump(logs[max(logs)])[1]
            if any(s["g"] == "MasJ" and 0 < s.get("tension", bv["tension"]) < (1 - 1e-6) * bv["tension"] for s in last):
                DMG = "tension"
    P("damage measure: '{}' below {:.0%} of its build value{}".format(
        DMG, DAMAGED_RATIO, "" if DAMAGE_PROP != "auto" else " (auto: the tension property of the last run is unchanged)"
        if DMG != "tension" else " (auto)"))
    zero_report(ref_subs, bv, P)
    results = []
    for rn in sorted(logs):
        blocks, subs, _, _ = read_dump(logs[rn])
        r = post("run_{:02d}".format(rn), rn, ref_blocks, bv, blocks, subs, P)
        if r:
            results.append(r)
    if not results:
        P("no run dumps processed"); return
    with open(os.path.join(OUT_DIR, "joints_by_run.csv"), "w", newline="") as f:
        rows = [j for r in results for j in r["joints"]]
        keys = sorted({k for j in rows for k in j}, key=lambda k: list(JKEYS).index(k) if k in JKEYS else 99)
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        [w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in j.items()}) for j in rows]
    onset = damage_onset(results, P)
    with open(os.path.join(OUT_DIR, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    plot(results, onset)


main()
