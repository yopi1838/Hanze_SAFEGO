# -*- coding: ascii -*-
"""
pocket_def_jobs.py -- +/- pushover of the pK11 pocket build with DEFORMABLE masonry blocks in SMALL STRAIN (6 Nov 2026).
Not called directly: run_pK11def_push_pos.py / _neg.py set JOB and exec this file; each runs in its own 3DEC session.

BUILD: ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11_defSS.dat, built by the + session only (the - session waits for the save)
    -> Part_I_MASON_v8_SmallStrain_springforce_zup_pk2_ft0p05_pK11_def.sav
DRIVER: pushover_oop_3dec_v5_zup.py, read and patched in memory (not edited):
    * model large-strain OFF (the driver forces it on);
    * the push a*m is applied to the masonry GRIDPOINTS (block.gp.force.app, a * zone mass / 4 on each of a zone's
      4 gridpoints) instead of block.force.app on rigid blocks;
    * masonry mass, participating mass and largest velocity from the zones (mass = zone volume x block.mass/block.vol
      of its block; 3DEC has no block.zone.density) and their gridpoints (velocity);
    * damping: global 0 (the build's gravity/spring solves use global damping), then the driver's local 0.9;
    * RESTORE_SAV / OUT_DIR / PUSH_DIR / TOP_JOINT 'mason' / TENSION_MAX_PA 0.06e6.
    Everything else (stepping, servo, diagnostics, CSV) unchanged.
OUTPUT: pushover_results_v5zup_pk2_ft0p05_pK11_defSS_TopMason/pushover_pos.csv, pushover_neg.csv (+ _end saves)
STRAIN (6 Nov 2026): "SS" (default) as above; "LS" builds a copy of the DAT with large strain ON from 'model new'
(gravity and spring solves in large strain, as the rigid builds via mat_case.dat), save ..._LargeStrain_..._def, and
keeps the driver's own 'model large-strain on'; outputs ..._defLS_TopMason. Wrappers run_pK11defLS_push_pos/_neg.py.
POCKET (8 Nov 2026): "pK11" (default) as above; "glued" drops the pK11 pocket set from the generated DAT, so the
joist-masonry contacts keep the joint set every contact gets (f_t 0.05 MPa, cohesion 0.2 MPa, ks 2.5e9; the rigid
ft0p05 build's pockets); save ..._glued_def, outputs pushover_results_v5zup_pk2_ft0p05_glued_defLS_TopMason.
Wrappers run_gluedDefLS_push_pos/_neg.py.
NOTE: in small strain there is no geometric softening (no P-Delta, contact geometry not updated), so after the peak
the curve cannot descend as the large-strain runs do; read the elastic branch, the peak and the mechanism.
"""
import os, re, time

JOB = globals().get("JOB")
if JOB not in ("push+", "push-"):
    raise RuntimeError("pocket_def_jobs.py: JOB must be 'push+' or 'push-' (got {!r})".format(JOB))
VARIANT = globals().get("VARIANT", "")   # "" = baseline (joints as rigid build), "kser" = series-corrected joints
if VARIANT not in ("", "kser"):
    raise RuntimeError("pocket_def_jobs.py: VARIANT must be '' or 'kser' (got {!r})".format(VARIANT))
STRAIN = globals().get("STRAIN", "SS")   # "SS" small strain (build and push), "LS" large strain (build and push)
if STRAIN not in ("SS", "LS"):
    raise RuntimeError("pocket_def_jobs.py: STRAIN must be 'SS' or 'LS' (got {!r})".format(STRAIN))
SFX = "_" + VARIANT if VARIANT else ""
POCKET = globals().get("POCKET", "pK11")   # "pK11" (timber-mortar pocket set) or "glued" (pockets keep the joint set)
if POCKET not in ("pK11", "glued"):
    raise RuntimeError("pocket_def_jobs.py: POCKET must be 'pK11' or 'glued' (got {!r})".format(POCKET))
SAVE0 = "Part_I_MASON_v8_SmallStrain_springforce_zup_pk2_ft0p05_pK11_def"      # the save DAT0 writes
SAVE = "Part_I_MASON_v8_{}_springforce_zup_pk2_ft0p05_{}_def{}".format(
    "LargeStrain" if STRAIN == "LS" else "SmallStrain", POCKET, SFX)
DAT0 = "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11_defSS.dat"
GENERATED = bool(VARIANT) or STRAIN == "LS" or POCKET != "pK11"
DAT = ("ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_{}_def{}{}.dat".format(POCKET, STRAIN, SFX)
       if GENERATED else DAT0)   # written from DAT0
OUT = "pushover_results_v5zup_pk2_ft0p05_{}_def{}{}_TopMason".format(POCKET, STRAIN, SFX)
BASE = "pushover_oop_3dec_v5_zup.py"
MAX_WAIT_MIN = 360
TAG = "# [pocket_def_jobs {}]".format(JOB)
LOG = "pocket_def_jobs_{}_{}{}.log".format(POCKET, STRAIN, SFX)

PO_PUSH = '''
fish define po_push(a_push)
    ; every local declared ONCE, at the top
    local z = null
    local gp = null
    local i = 0
    local fz = 0.0
    loop foreach z block.zone.list
        if block.isgroup(block.zone.hostblock(z), 'Masonry') then
            loop i (1, 4)
                block.gp.force.app(block.zone.gp(z, i)) = vector(0.0, 0.0, 0.0)
            endloop
        endif
    endloop
    loop foreach z block.zone.list
        if block.isgroup(block.zone.hostblock(z), 'Masonry') then
            fz = 0.25 * a_push * block.zone.vol(z) * block.mass(block.zone.hostblock(z)) / block.vol(block.zone.hostblock(z))
            loop i (1, 4)
                gp = block.zone.gp(z, i)
                block.gp.force.app(gp) = block.gp.force.app(gp) + vector(0.0, fz, 0.0)
            endloop
        endif
    endloop
end
'''
WALL_MASS = '''
fish define wall_mass_part
    ; every local declared ONCE, at the top
    local z = null
    local i = 0
    local mz = 0.0
    local vy = 0.0
    global m_tot = 0.0
    global m_part = 0.0
    global wall_maxvel = 0.0
    loop foreach z block.zone.list
        if block.isgroup(block.zone.hostblock(z), 'Masonry') then
            mz = block.zone.vol(z) * block.mass(block.zone.hostblock(z)) / block.vol(block.zone.hostblock(z))
            vy = 0.0
            loop i (1, 4)
                vy = vy + 0.25 * block.gp.vel.y(block.zone.gp(z, i))
            endloop
            m_tot = m_tot + mz
            if (vy * push_sign) > 0.0 then
                m_part = m_part + mz
            endif
            if math.abs(vy) > wall_maxvel then
                wall_maxvel = math.abs(vy)
            endif
        endif
    endloop
end
'''


def log(s):
    line = "[{}] [{}def{}{} {}] {}".format(time.strftime("%H:%M:%S"), POCKET, STRAIN, SFX, JOB, s)
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def have_save():
    return os.path.isfile(SAVE + ".sav") or os.path.isfile(SAVE)


def ensure_save():
    """No lock file any more (a crashed build left one behind and the other session waited on it for ever).
    The + session is the only one that builds; the - session only waits until the save exists and its size has
    stopped changing for a minute."""
    if JOB == "push+":
        if have_save():
            log("save {} found".format(SAVE)); return
        if GENERATED:
            write_variant_dat()
        if not os.path.isfile(DAT):
            raise RuntimeError("{} not found".format(DAT))
        import itasca as it
        log("building the save: program call '{}'".format(DAT))
        it.command("python-reset-state false")
        it.command("program call '{}'".format(DAT))
        it.command("python-reset-state false")
        if not have_save():
            raise RuntimeError("{} did not write {}".format(DAT, SAVE))
        log("save built; check the DEFSS lines (zones, zone mass ~1308 kg, pocket stiffness-shear {})".format(
            "2.5e9" if POCKET == "glued" else "1.1e8"))
        return
    path = SAVE + ".sav" if os.path.isfile(SAVE + ".sav") or not os.path.isfile(SAVE) else SAVE
    log("waiting for the save from the + session ({})".format(path))
    t0, last, still = time.time(), -1, 0
    while still < 3:
        time.sleep(20)
        size = os.path.getsize(path) if os.path.isfile(path) else -1
        still = still + 1 if (size > 0 and size == last) else 0
        last = size
        if time.time() - t0 > MAX_WAIT_MIN * 60:
            raise RuntimeError("no save after {} min: start the + wrapper first (it builds the save)"
                               .format(MAX_WAIT_MIN))
    log("save ready after {:.1f} min".format((time.time() - t0) / 60.0))


PK11_BLOCK = ("; ---- pK11 pocket: timber-mortar set on the joist-masonry contacts (Arslan et al. 2024) ----\n"
              "block contact prop stiffness-shear 1.1e8 tension 0 cohesion 0.1e6 friction 31.0 ...\n"
              "    friction-residual 31.0 cohesion-residual 0 ...\n"
              "    range group-intersection 'Joist' 'Masonry'\n")


def write_variant_dat():
    """DAT0 with the variant's switches and save name -> DAT (rewritten at every build, so it follows DAT0).
    LS: both 'model large-strain off' command lines (after model new, before gravity) -> on.
    kser: KSERIES = 1."""
    with open(DAT0) as f:
        src = f.read()
    subs = [("\nmodel save '{}'\n".format(SAVE0), "\nmodel save '{}'\n".format(SAVE), 1)]
    if STRAIN == "LS":
        subs.append(("\nmodel large-strain off\n", "\nmodel large-strain on\n", 2))
    if VARIANT == "kser":
        subs.append(("[global KSERIES = 0]", "[global KSERIES = 1]", 1))
    if POCKET == "glued":   # the joist-masonry contacts keep the joint set given to every contact (f_t 0.05 MPa)
        subs.append((PK11_BLOCK, "; ---- glued pocket (POCKET 'glued'): the joist-masonry contacts keep the joint set above "
                                 "(f_t 0.05 MPa, cohesion 0.2 MPa, ks 2.5e9) ----\n", 1))
        subs.append(("(must be 1.1e8)", "(glued: must be 2.5e9, the joint set)", 1))
    for old, new, n in subs:
        if src.count(old) != n:
            raise RuntimeError("{}: expected {!r} {} time(s), found {}".format(DAT0, old, n, src.count(old)))
        src = src.replace(old, new)
    if STRAIN == "LS" and re.search(r"^model large-strain off", src, re.M):
        raise RuntimeError("{}: a 'model large-strain off' command is left".format(DAT))
    with open(DAT, "w") as f:
        f.write("; GENERATED by pocket_def_jobs.py (POCKET {}, STRAIN {}, VARIANT {!r}) from {} -- edit that file, not this one\n"
                .format(POCKET, STRAIN, VARIANT, DAT0) + src)
    log("wrote {} (STRAIN {}, KSERIES {}, save {})".format(DAT, STRAIN, int(VARIANT == "kser"), SAVE))


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def replace_fish(src, name, new):
    i = src.index("fish define {}".format(name))
    a = src.rindex('it.command("""', 0, i) + len('it.command("""')
    b = src.index('""")', i)
    return src[:a] + new + src[b:]


def patched(direction):
    with open(BASE) as f:
        src = f.read()
    if "Z-UP" not in src:
        raise RuntimeError("{} is not the z-up driver".format(BASE))
    old = 'it.command("model large-strain on")'
    if src.count(old) != 1:
        raise RuntimeError("{}: large-strain line not found once".format(BASE))
    if STRAIN == "SS":
        src = src.replace(old, 'it.command("model large-strain off")   ' + TAG)   # LS: the driver's own line stays
    # the build solves gravity and the spring force with GLOBAL damping: switch it off before the driver's local 0.9
    old = 'it.command("block mech damp local 0.9")'
    if src.count(old) != 1:
        raise RuntimeError("{}: damping line not found once".format(BASE))
    src = src.replace(old, 'it.command("block mech damp global 0.0"); ' + old + '   ' + TAG)
    src = replace_fish(src, "po_push(a_push)", PO_PUSH)
    src = replace_fish(src, "wall_mass_part", WALL_MASS)
    for k, v in (("RESTORE_SAV", SAVE), ("OUT_DIR", OUT), ("PUSH_DIR", direction), ("TOP_JOINT", "mason"),
                 ("TENSION_MAX_PA", 0.06e6)):
        src = set_last(src, k, v)
    return src


def run(direction):
    side = "pos" if direction > 0 else "neg"
    csvp = os.path.join(OUT, "pushover_{}.csv".format(side))
    endsav = os.path.join(OUT, "pushover_{}_end".format(side))
    if os.path.isfile(endsav + ".sav") or os.path.isfile(endsav):
        log("{} already finished -- nothing to do".format(csvp)); return
    if os.path.isfile(csvp):
        os.replace(csvp, csvp + ".partial"); log("{}: no end save, kept as .partial, pushing again".format(csvp))
    src = patched(direction)
    log("{} push, deformable masonry, {} strain: {} -> {}".format("+" if direction > 0 else "-", "large" if STRAIN == "LS" else "small", SAVE, csvp))
    t0 = time.time()
    exec(compile(src, "{} [{} def {} {}]".format(BASE, POCKET, STRAIN, side), "exec"),
         {"__name__": "__main__", "__file__": os.path.abspath(BASE)})
    log("done in {:.1f} min".format((time.time() - t0) / 60.0))


log("=" * 60)
ensure_save()
run(+1 if JOB == "push+" else -1)
