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
NOTE: in small strain there is no geometric softening (no P-Delta, contact geometry not updated), so after the peak
the curve cannot descend as the large-strain runs do; read the elastic branch, the peak and the mechanism.
"""
import os, re, time

JOB = globals().get("JOB")
if JOB not in ("push+", "push-"):
    raise RuntimeError("pocket_def_jobs.py: JOB must be 'push+' or 'push-' (got {!r})".format(JOB))
SAVE = "Part_I_MASON_v8_SmallStrain_springforce_zup_pk2_ft0p05_pK11_def"
DAT = "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_pK11_defSS.dat"
OUT = "pushover_results_v5zup_pk2_ft0p05_pK11_defSS_TopMason"
BASE = "pushover_oop_3dec_v5_zup.py"
MAX_WAIT_MIN = 360
TAG = "# [pocket_def_jobs {}]".format(JOB)
LOG = "pocket_def_jobs.log"

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
    line = "[{}] [pK11def {}] {}".format(time.strftime("%H:%M:%S"), JOB, s)
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
        if not os.path.isfile(DAT):
            raise RuntimeError("{} not found".format(DAT))
        import itasca as it
        log("building the save: program call '{}'".format(DAT))
        it.command("python-reset-state false")
        it.command("program call '{}'".format(DAT))
        it.command("python-reset-state false")
        if not have_save():
            raise RuntimeError("{} did not write {}".format(DAT, SAVE))
        log("save built; check the DEFSS lines (zones, zone mass ~1308 kg, pocket stiffness-shear 1.1e8)")
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
            raise RuntimeError("no save after {} min: start run_pK11def_push_pos.py first (it builds the save)"
                               .format(MAX_WAIT_MIN))
    log("save ready after {:.1f} min".format((time.time() - t0) / 60.0))


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
    src = src.replace(old, 'it.command("model large-strain off")   ' + TAG)
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
    log("{} push, deformable masonry, small strain: {} -> {}".format("+" if direction > 0 else "-", SAVE, csvp))
    t0 = time.time()
    exec(compile(src, "{} [pK11 def SS {}]".format(BASE, side), "exec"),
         {"__name__": "__main__", "__file__": os.path.abspath(BASE)})
    log("done in {:.1f} min".format((time.time() - t0) / 60.0))


log("=" * 60)
ensure_save()
run(+1 if JOB == "push+" else -1)
