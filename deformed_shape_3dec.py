# -*- coding: ascii -*-
"""
deformed_shape_3dec.py -- deformed shape and hinge map from the pushover saves (3DEC 9).

WHAT IT DOES
  1. Restores the reference save (the pushover's RESTORE_SAV) and records every Masonry / T_B / Joist
     block centroid.
  2. Restores each pushover save in SAVES and records
       - every block centroid again  -> displacement = current - reference (large strain: positions update);
       - every subcontact of the joint groups the v5 driver defines in slot "vcheck"
         (MasJ = masonry-masonry, Pocket, Reaction_topTB, Reaction_base): position, normal displacement,
         normal force, area, and the contact normal (bed joints vs head joints).
  3. Post-processes:
       - course profile: mass-weighted lateral displacement dz of each course vs its original height,
         the rotation of each segment between courses and the kink (change of rotation) at each course;
       - bed-joint map: for each bed joint the largest opening, the open fraction of its area and the
         face it opens on (leading = the face the wall is pushed towards, trailing = the other one).
         A three-hinge mechanism shows: base and top hinges opening on the TRAILING face, the
         mid-span hinge opening on the LEADING face;
       - top joint, base and pocket openings; the displacement at the 2.06 m sensor height as a check
         against d_ctrl_mm in the pushover CSV.

OUTPUT  <OUT_DIR>/<tag>_blocks.csv    one row per block: id, group, original x/y/z, dx/dy/dz (mm)
        <OUT_DIR>/<tag>_courses.csv   course profile with rotations and kinks
        <OUT_DIR>/<tag>_joints.csv    bed-joint map (masonry), plus base / top joint rows
        <OUT_DIR>/<tag>_summary.txt   the printed tables
        <OUT_DIR>/deformed_shape.png  if matplotlib is available
        <OUT_DIR>/_dump_<tag>.log     raw dump (the post-processing can be re-run from these)

USAGE   inside 3DEC:   program call 'deformed_shape_3dec.py'    (or run it from the Python pane)
        outside 3DEC:  python deformed_shape_3dec.py --post     (re-processes the existing dump logs)

JOINT HEIGHTS (v2, 23 Oct)
  Joints are reported at their ORIGINAL height (y0). The dump now also writes the two block ids of every contact
  (block.contact.b1 / b2), so a joint sits exactly halfway between its two courses. Dumps written before this
  change have no ids; --post then assigns each subcontact geometrically (nearest course pair), which is exact away
  from the hinges but can share a wide-open joint's subcontacts with its neighbour (two adjacent rows with the same
  max opening). Re-dump inside 3DEC for exact rows. Pockets are split by floor (pocket_lower / pocket_upper), and
  the top-joint tilt and |N| sin(tilt) are printed: the pushover's F_top sums only the SHEAR z-components of the
  top joint, so once the top beam rotates F_tot over-reports the horizontal force by |N| sin(tilt).

ASSUMPTIONS TO CHECK ON THE FIRST RUN
  * block.subcontact.disp.norm > 0 is OPENING (the v5 driver's max_open_mm assumes the same).
  * block.contact.normal exists in your build; if it does not, the script falls back to a dump
    without normals and treats every masonry contact whose normal it cannot read as a bed joint
    candidate (head joints are then filtered by their small vertical extent -- see BED_ONLY).
  * Block ids are the same in the reference save and the pushover saves (same model lineage).
"""

import os, sys, math, csv

# ============================ CONFIG =================================
REF_SAV = "Part_I_MASON_v8_LargeStrain_springforce"
PO_DIR  = "pushover_results_v5_ft0_masonrypush"
SAVES = [   # (tag, save file, PUSH_DIR of that run)
    ("neg_end", PO_DIR + "/pushover_neg_end", -1),
    ("pos_end", PO_DIR + "/pushover_pos_end", +1),
]
OUT_DIR       = "deformed_shape"
T_WALL        = 0.21      # m
SENSOR_Y      = 2.06      # m, Top_Quarter_A/B height
OPEN_MM       = 0.01      # a subcontact counts as open above this
HINGE_MM      = 0.5       # a bed joint is reported as a hinge above this largest opening
COURSE_GAP    = 0.02      # m, centroid-height gap that separates two courses / two joints
BED_NY_MIN    = 0.8       # |normal_y| above this = bed joint
FLOOR_Y       = (0.38, 2.14)   # m, joist bearings (drawn on the figure)
FLOOR_SPLIT_Y = 1.26      # m, pockets below this = lower floor (0.38 m), above = upper floor (2.14 m)

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

FISH_CONTACTS = r"""
fish define dfm_cgroup(cx)
    local g = ''
    if block.contact.isgroup(cx, 'MasJ', 'vcheck') then
        g = 'MasJ'
    endif
    if block.contact.isgroup(cx, 'Pocket', 'vcheck') then
        g = 'Pocket'
    endif
    if block.contact.isgroup(cx, 'Reaction_topTB', 'vcheck') then
        g = 'TopTB'
    endif
    if block.contact.isgroup(cx, 'Reaction_base', 'vcheck') then
        g = 'Base'
    endif
    dfm_cgroup = g
end
fish define dfm_contacts
    loop foreach local cx block.contact.list
        local g = dfm_cgroup(cx)
        if g # '' then
            local ny = NORMAL_EXPR
            loop foreach local sc block.contact.subcontactlist(cx)
                local p = block.subcontact.pos(sc)
                io.out('@@C,' + g + ',' + string(block.contact.id(cx)) + ',' + string(ny) + ',' + string(comp.x(p)) + ',' + string(comp.y(p)) + ',' + string(comp.z(p)) + ',' + string(block.subcontact.disp.norm(sc)) + ',' + string(block.subcontact.force.norm(sc)) + ',' + string(block.subcontact.area(sc)) + BIDS_EXPR)
            endloop
        endif
    endloop
end
"""


def _cmdpath(p):
    return p.replace("\\", "/")


def run_dump(it, tag, save, contacts):
    log = os.path.join(OUT_DIR, "_dump_{}.log".format(tag))
    if os.path.exists(log):
        os.remove(log)
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(_cmdpath(save)))
    it.command("python-reset-state false")
    it.command(FISH_BLOCKS)
    have_normal = True
    if contacts:
        bids = "',' + string(block.id(block.contact.b1(cx))) + ',' + string(block.id(block.contact.b2(cx)))"
        done = False
        for nexpr, bexpr in (("comp.y(block.contact.normal(cx))", bids), ("comp.y(block.contact.normal(cx))", ""),
                             ("-9.0", bids), ("-9.0", "")):
            try:
                it.command(FISH_CONTACTS.replace("NORMAL_EXPR", nexpr).replace("BIDS_EXPR", bexpr or "''"))
                done = True
                break
            except Exception as e:
                print("contact dump variant failed ({}); trying a simpler one".format(e))
        if not done:
            raise RuntimeError("could not define the contact dump")
        have_normal = nexpr != "-9.0"
        print("contact dump: normals {}, block ids {}".format("yes" if have_normal else "no", "yes" if bexpr else "no"))
    it.command('program log-file "{}"'.format(_cmdpath(log)))
    it.command("program log on")
    it.command("[dfm_blocks]")
    if contacts:
        it.command("[dfm_contacts]")
    it.command("program log off")
    print("dumped {} -> {} (normals: {})".format(save, log, have_normal))
    return log


# ============================ PARSING ================================
def read_dump(log):
    blocks, subs = {}, []
    with open(log, errors="replace") as f:
        for line in f:
            i = line.find("@@B,")
            if i >= 0:
                s = line[i + 4:].strip().split(",")
                blocks[int(float(s[1]))] = (s[0], float(s[2]), float(s[3]), float(s[4]), float(s[5]))
                continue
            i = line.find("@@C,")
            if i >= 0:
                s = line[i + 4:].strip().split(",")
                d = dict(g=s[0], cid=int(float(s[1])), ny=float(s[2]), x=float(s[3]), y=float(s[4]),
                         z=float(s[5]), un=float(s[6]), fn=float(s[7]), area=float(s[8]), b1=None, b2=None)
                if len(s) >= 11 and s[9].strip() and s[10].strip():
                    d["b1"], d["b2"] = int(float(s[9])), int(float(s[10]))
                subs.append(d)
    return blocks, subs


def cluster(values, gap):
    """indices grouped by sorted value with a gap > gap between groups"""
    order = sorted(range(len(values)), key=lambda k: values[k])
    groups, cur = [], []
    for k in order:
        if cur and values[k] - values[cur[-1]] > gap:
            groups.append(cur); cur = []
        cur.append(k)
    if cur:
        groups.append(cur)
    return groups


# ============================ POST ===================================
def post(tag, push_dir, ref_blocks, blocks, subs, out_lines):
    P = lambda s="": (print(s), out_lines.append(s))
    lead = "+z" if push_dir > 0 else "-z"
    # ---- blocks ----
    rows = []
    for bid, (g, x, y, z, m) in blocks.items():
        if bid not in ref_blocks:
            continue
        _, x0, y0, z0, _ = ref_blocks[bid]
        rows.append(dict(id=bid, group=g, x0=x0, y0=y0, z0=z0, mass=m,
                         dx_mm=(x - x0) * 1e3, dy_mm=(y - y0) * 1e3, dz_mm=(z - z0) * 1e3, zc=z))
    with open(os.path.join(OUT_DIR, tag + "_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "group", "x0", "y0", "z0", "mass", "dx_mm", "dy_mm", "dz_mm"],
                           extrasaction="ignore")
        w.writeheader(); [w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}) for r in rows]
    mas = [r for r in rows if r["group"] == "Masonry"]
    if not mas:
        P("no Masonry blocks matched between the reference and {} -- block ids differ?".format(tag)); return None
    # ---- course profile ----
    crs = []
    for grp in cluster([r["y0"] for r in mas], COURSE_GAP):
        bl = [mas[k] for k in grp]; M = sum(b["mass"] for b in bl)
        wm = lambda key: sum(b[key] * b["mass"] for b in bl) / M
        crs.append(dict(y0=wm("y0"), n=len(bl), dz_mm=wm("dz_mm"), dz_min=min(b["dz_mm"] for b in bl),
                        dz_max=max(b["dz_mm"] for b in bl), dy_mm=wm("dy_mm"), zc=wm("zc")))
    for i, c in enumerate(crs):
        c["rot_mrad"] = float("nan") if i == 0 else (c["dz_mm"] - crs[i - 1]["dz_mm"]) / (c["y0"] - crs[i - 1]["y0"])
    for i, c in enumerate(crs):
        c["kink_mrad"] = (crs[i + 1]["rot_mrad"] - c["rot_mrad"]) if 0 < i < len(crs) - 1 else float("nan")
    with open(os.path.join(OUT_DIR, tag + "_courses.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["y0", "n", "dz_mm", "dz_min", "dz_max", "dy_mm", "rot_mrad", "kink_mrad"],
                           extrasaction="ignore")
        w.writeheader(); [w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in c.items()}) for c in crs]
    ys = [c["y0"] for c in crs]; dzs = [c["dz_mm"] for c in crs]
    d_sens = interp(SENSOR_Y, ys, dzs)
    P("=" * 78)
    P("{}: push {:+d} (leading face {}), {} masonry blocks in {} courses".format(tag, push_dir, lead, len(mas), len(crs)))
    P("lateral displacement at the sensor height {:.2f} m: {:.2f} mm  (compare d_ctrl_mm, last CSV row)".format(SENSOR_Y, d_sens))
    tb = [r for r in rows if r["group"] == "T_B"]
    if tb:
        P("top beam: dy {:+.2f} mm, dz {:+.2f} mm (should be ~0: velocity-z fixed)".format(tb[0]["dy_mm"], tb[0]["dz_mm"]))
    P("")
    P("  y0 (m)   dz (mm)  [min  max]     dy (mm)  rot (mrad)  kink (mrad)")
    for c in crs:
        P("  {:6.3f}  {:8.2f}  [{:6.2f} {:6.2f}]  {:7.2f}  {:9.2f}  {:9.2f}".format(
            c["y0"], c["dz_mm"], c["dz_min"], c["dz_max"], c["dy_mm"], c["rot_mrad"], c["kink_mrad"]))
    # ---- reference (undeformed) height of every subcontact ----
    # Joints are grouped by the height they had before the push, not by their current height: after a 10 deg
    # rotation the current heights of neighbouring joints overlap and would be merged into one row.
    yref_of = reference_heights(subs, mas, crs, ref_blocks, blocks)
    for s_, yr in zip(subs, yref_of):
        s_["yref"] = yr
    # ---- joints ----
    jrows = []
    if subs:
        def joint_stats(sub, label):
            A = sum(s["area"] for s in sub) or 1.0
            zc = sum(s["z"] * s["area"] for s in sub) / A
            yc = sum(s["yref"] * s["area"] for s in sub) / A
            smax = max(sub, key=lambda s: s["un"])
            op = [s for s in sub if s["un"] * 1e3 > OPEN_MM]
            Aop = sum(s["area"] for s in op)
            zw = (sum((s["z"] - zc) * s["un"] * s["area"] for s in op) / sum(s["un"] * s["area"] for s in op)) if op else 0.0
            face = "" if not op else ("+z" if zw > 0 else "-z")
            side = "" if not op else ("leading" if face == lead else "trailing")
            return dict(kind=label, y=yc, n_sub=len(sub), max_open_mm=smax["un"] * 1e3,
                        open_frac=Aop / A, face=face, side=side, z_maxopen_rel_mm=(smax["z"] - zc) * 1e3,
                        N_kN=-sum(s["fn"] for s in sub) / 1e3)
        mj = [s for s in subs if s["g"] == "MasJ"]
        have_n = any(s["ny"] > -2 for s in mj)
        bed = [s for s in mj if (abs(s["ny"]) >= BED_NY_MIN if have_n else True)]
        head = [s for s in mj if have_n and abs(s["ny"]) < BED_NY_MIN]
        for grp in cluster([s["yref"] for s in bed], COURSE_GAP):
            sub = [bed[k] for k in grp]
            if len(sub) < 4:
                continue
            jrows.append(joint_stats(sub, "bed"))
        for g, label, sel in (("Base", "base", None), ("TopTB", "top_joint", None),
                              ("Pocket", "pocket_lower", lambda s: s["yref"] < FLOOR_SPLIT_Y),
                              ("Pocket", "pocket_upper", lambda s: s["yref"] >= FLOOR_SPLIT_Y)):
            sub = [s for s in subs if s["g"] == g and (sel is None or sel(s))]
            if sub:
                jrows.append(joint_stats(sub, label))
        jrows.sort(key=lambda r: r["y"])
        with open(os.path.join(OUT_DIR, tag + "_joints.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(jrows[0])); w.writeheader()
            [w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()}) for r in jrows]
        P("")
        P("  joint       y0 (m)  max open (mm)  open frac  opens on        N (kN)   [normals: {}]".format("yes" if have_n else "NO -- head joints mixed in"))
        for r in jrows:
            flag = "  <- HINGE" if (r["kind"] in ("bed", "base") and r["max_open_mm"] > HINGE_MM) else ""
            P("  {:12s} {:6.3f} {:10.3f}     {:6.2f}    {:3s} {:9s}  {:7.2f}{}".format(
                r["kind"], r["y"], r["max_open_mm"], r["open_frac"], r["face"], r["side"], r["N_kN"], flag))
        if head:
            ho = [s for s in head if s["un"] * 1e3 > HINGE_MM]
            P("  head joints: {} subcontacts, {} open > {} mm (largest {:.2f} mm at y0 {:.2f} m)".format(
                len(head), len(ho), HINGE_MM, max(s["un"] for s in head) * 1e3,
                max(head, key=lambda s: s["un"])["yref"]))
        top = [s for s in subs if s["g"] == "TopTB"]
        if len(top) >= 3:
            th = plane_tilt(top)
            Ntop = -sum(s["fn"] for s in top) / 1e3
            P("  top joint tilt {:.1f} mrad (the top beam rotates with it if the joint stays closed); "
              "|N| sin(tilt) = {:.2f} kN".format(th * 1e3, abs(Ntop) * math.sin(abs(th))))
            P("    -> the pushover's F_top sums the subcontact SHEAR z-components only; the tilted normal force adds")
            P("       this much horizontal force that F_tot leaves out. On the descending branch compare |a| x m_tot.")
        hinges = [r for r in jrows if r["kind"] in ("bed", "base") and r["max_open_mm"] > HINGE_MM]
        P("")
        P("HINGES (bed joints opening > {} mm): ".format(HINGE_MM) + (", ".join(
            "{:.2f} m ({}, {:.1f} mm)".format(r["y"], r["side"], r["max_open_mm"]) for r in hinges) or "none"))
        P("three-hinge reading: base & top hinges open on the trailing face, the mid-span hinge on the leading face.")
    return dict(tag=tag, crs=crs, joints=jrows, push=push_dir)


def plane_tilt(sub):
    """tilt (rad) of a joint plane in the y-z plane from its subcontact positions (least squares y = c + k z)"""
    A = sum(s["area"] for s in sub)
    zm = sum(s["z"] * s["area"] for s in sub) / A
    ym = sum(s["y"] * s["area"] for s in sub) / A
    szz = sum(s["area"] * (s["z"] - zm) ** 2 for s in sub)
    szy = sum(s["area"] * (s["z"] - zm) * (s["y"] - ym) for s in sub)
    return math.atan(szy / szz) if szz > 0 else 0.0


def reference_heights(subs, mas, crs, ref_blocks, blocks):
    """undeformed height of each subcontact.
    With block ids in the dump: a masonry-masonry contact sits halfway between its two blocks' original centroids
    (a bed joint lands on the joint, a head joint on its course); a contact with one masonry block (pocket, base,
    top joint) is placed by the geometric rule below on that block.
    Without block ids (older dumps): a bed-joint subcontact takes the nearest masonry block A and the nearest block
    B in the course just above or below A (whichever is closer) and sits halfway between their original centroids;
    any other subcontact uses the geometric rule on the nearest masonry block -- its offset from that block's
    current centroid along the block's rotated vertical axis (course rotation) is added to the original height."""
    course_of, rot_of = {}, []
    for grp_i, grp in enumerate(cluster([r["y0"] for r in mas], COURSE_GAP)):
        for k in grp:
            course_of[mas[k]["id"]] = grp_i
    for i in range(len(crs)):
        v = [crs[j]["rot_mrad"] for j in (i, i + 1) if 0 <= j < len(crs) and crs[j]["rot_mrad"] == crs[j]["rot_mrad"]]
        rot_of.append(math.atan(sum(v) / len(v) / 1e3) if v else 0.0)
    mas_ids = [r["id"] for r in mas]
    cur = {b: blocks[b] for b in mas_ids if b in blocks}

    def local(bid, s):
        _, _, yc, zc, _ = cur[bid]
        th = rot_of[course_of[bid]]
        return ref_blocks[bid][2] + (s["z"] - zc) * math.sin(th) + (s["y"] - yc) * math.cos(th)

    d2 = lambda b, s: (cur[b][1] - s["x"]) ** 2 + (cur[b][2] - s["y"]) ** 2 + (cur[b][3] - s["z"]) ** 2
    by_course = {}
    for b in mas_ids:
        by_course.setdefault(course_of[b], []).append(b)

    def nearest(s, pool=None):
        return min(pool if pool is not None else mas_ids, key=lambda b: d2(b, s))
    out = []
    for s in subs:
        ids = [b for b in (s.get("b1"), s.get("b2")) if b is not None and b in cur]
        if len(ids) == 2:
            out.append(0.5 * (ref_blocks[ids[0]][2] + ref_blocks[ids[1]][2]))
        elif len(ids) == 1:
            out.append(local(ids[0], s))
        elif s["g"] == "MasJ" and not (-2 < s["ny"] and abs(s["ny"]) < BED_NY_MIN):
            a = nearest(s)
            ca = course_of[a]
            cand = [nearest(s, by_course[c]) for c in (ca - 1, ca + 1) if c in by_course]
            b = min(cand, key=lambda q: d2(q, s)) if cand else a
            out.append(0.5 * (ref_blocks[a][2] + ref_blocks[b][2]))
        else:
            out.append(local(nearest(s), s))
    return out


def interp(x, xs, ys):
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


def plot(results):
    try:
        import matplotlib
        matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print("matplotlib not available ({}); CSVs written, no figure".format(e)); return
    fig, ax = plt.subplots(1, 2, figsize=(11, 6.5), sharey=True)
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
    for k, r in enumerate(results):
        c = cols[k % len(cols)]
        ax[0].plot([q["dz_mm"] for q in r["crs"]], [q["y0"] for q in r["crs"]], "o-", ms=3, color=c, label=r["tag"])
        bed = [j for j in r["joints"] if j["kind"] in ("bed", "base")]
        if bed:
            ax[1].plot([j["max_open_mm"] for j in bed], [j["y"] for j in bed], "-", color=c, lw=1, alpha=0.6)
            for side, mk in (("leading", "^"), ("trailing", "v")):
                s = [j for j in bed if j["side"] == side]
                ax[1].plot([j["max_open_mm"] for j in s], [j["y"] for j in s], mk, ms=6, color=c, ls="none",
                           label="{} opens {}".format(r["tag"], side))
    ax[0].axhline(SENSOR_Y, color="grey", ls=":", lw=0.8)
    for a_ in ax:
        for yf in FLOOR_Y:
            a_.axhline(yf, color="#a16207", ls="--", lw=0.8, alpha=0.7)
    ax[1].text(0.98, FLOOR_Y[1] + 0.02, "upper floor", color="#a16207", fontsize=8,
               ha="right", transform=ax[1].get_yaxis_transform())
    ax[1].text(0.98, FLOOR_Y[0] + 0.02, "lower floor", color="#a16207", fontsize=8, ha="right",
               transform=ax[1].get_yaxis_transform())
    ax[0].set_xlabel("course lateral displacement dz (mm)"); ax[0].set_ylabel("original height (m)")
    ax[0].grid(alpha=0.25); ax[0].legend(fontsize=8); ax[0].set_title("Deformed shape", loc="left", fontsize=10)
    ax[1].axvline(HINGE_MM, color="grey", ls=":", lw=0.8)
    ax[1].set_xlabel("largest bed-joint opening (mm)"); ax[1].grid(alpha=0.25); ax[1].legend(fontsize=8)
    ax[1].set_title("Bed-joint openings: where the hinges are", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "deformed_shape.png"), dpi=150)
    print("-> " + os.path.join(OUT_DIR, "deformed_shape.png"))


# ============================ MAIN ===================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    post_only = "--post" in sys.argv
    if not post_only:
        try:
            import itasca as it
        except ImportError:
            print("not inside 3DEC: re-processing existing dumps (same as --post)"); post_only = True
    ref_log = os.path.join(OUT_DIR, "_dump_ref.log")
    logs = {}
    if not post_only:
        ref_log = run_dump(it, "ref", REF_SAV, contacts=False)
        for tag, save, _ in SAVES:
            if os.path.isfile(save) or os.path.isfile(save + ".sav"):
                logs[tag] = run_dump(it, tag, save, contacts=True)
            else:
                print("skip {}: {}(.sav) not found".format(tag, save))
    else:
        for tag, _, _ in SAVES:
            p = os.path.join(OUT_DIR, "_dump_{}.log".format(tag))
            if os.path.isfile(p):
                logs[tag] = p
    ref_blocks, _ = read_dump(ref_log)
    results = []
    for tag, _, push in SAVES:
        if tag not in logs:
            continue
        lines = []
        blocks, subs = read_dump(logs[tag])
        r = post(tag, push, ref_blocks, blocks, subs, lines)
        with open(os.path.join(OUT_DIR, tag + "_summary.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        if r:
            results.append(r)
    if results:
        plot(results)


main()
