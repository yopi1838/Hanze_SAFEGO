# -*- coding: ascii -*-
"""
pocket_jobs.py -- the three jobs on a POCKET variant of the pk2 f_t 0.05 MPa build, each in its OWN 3DEC session so they
run in parallel (5 Nov 2026). Not called directly: the wrappers set JOB (and POCKET) and exec this file.

    run_pK11_push_pos.py   JOB "push+"   pushover_oop_3dec_v5_zup.py, PUSH_DIR +1
    run_pK11_push_neg.py   JOB "push-"   pushover_oop_3dec_v5_zup.py, PUSH_DIR -1
    run_pK11_stratF.py     JOB "stratF"  run_stratF_zup.py, Strategy F runs 1-25, 1.5 % Maxwell

THE VARIANT (ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_<POCKET>.dat): the zup_pk2_ft0p05 save with the joist-pocket contacts
given the timber-mortar set of Arslan et al. (2024): ks 1.1e8 (pK11) / 3.6e8 (pK36) Pa/m, cohesion 0.1 MPa, friction
31 deg, tension 0; every other contact unchanged.

SAVE BUILT ONCE: the first session to find the save missing builds it (program call the .dat) under a lock file
<save>.building; the other sessions wait until the save exists and the lock is gone. A lock left by a crashed build
must be deleted by hand (the waiting session says so after LOCK_STALE_MIN).
OUTPUTS (no two sessions write the same file):
    pushover_results_v5zup_pk2_ft0p05_<POCKET>_TopMason/pushover_pos.csv, pushover_pos_end.sav  (push+)
                                                         pushover_neg.csv, pushover_neg_end.sav  (push-)
    stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_<POCKET>_tr/                                          (stratF)
A pushover counts as done only when its _end save exists: a CSV without it (a stopped run) is renamed *.partial and
the push starts again. Strategy F resumes from its own checkpoint.
"""
import os, re, time

POCKET = globals().get("POCKET", "pK11")
JOB = globals().get("JOB")
DAMP_PCT = globals().get("DAMP_PCT", 1.5)
if JOB not in ("push+", "push-", "stratF"):
    raise RuntimeError("pocket_jobs.py: JOB must be 'push+', 'push-' or 'stratF' (got {!r})".format(JOB))

SAVE = "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_" + POCKET
DAT = "ANALYSIS_PART_I_MASON_zup_pk2_ft0p05_{}.dat".format(POCKET)
LOCK = SAVE + ".building"
LOCK_STALE_MIN = 90
TAG = "# [pocket_jobs {} {}]".format(POCKET, JOB)
LOG = "pocket_jobs_{}.log".format(POCKET)


def log(s):
    line = "[{}] [{} {}] {}".format(time.strftime("%H:%M:%S"), POCKET, JOB, s)
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def have_save():
    return os.path.isfile(SAVE + ".sav") or os.path.isfile(SAVE)


def ensure_save():
    if have_save() and not os.path.isfile(LOCK):
        log("save {} found".format(SAVE))
        return
    if not os.path.isfile(DAT):
        raise RuntimeError("{} not found".format(DAT))
    try:
        fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, "built by {} at {}\n".format(JOB, time.ctime()).encode())
        os.close(fd)
        builder = True
    except FileExistsError:
        builder = False
    if builder:
        import itasca as it
        try:
            log("building the save: program call '{}'".format(DAT))
            it.command("python-reset-state false")
            it.command("program call '{}'".format(DAT))
            it.command("python-reset-state false")
        finally:
            os.remove(LOCK)
        if not have_save():
            raise RuntimeError("{} did not write {}".format(DAT, SAVE))
        log("save built; check the POCKET ... AFTER line above (stiffness-shear, cohesion 1e5)")
        return
    log("another session is building {}: waiting".format(SAVE))
    t0 = time.time()
    while os.path.isfile(LOCK) or not have_save():
        time.sleep(20)
        if os.path.isfile(LOCK) and time.time() - os.path.getmtime(LOCK) > LOCK_STALE_MIN * 60:
            raise RuntimeError("{} is older than {} min: the build probably crashed. Delete it and start again."
                               .format(LOCK, LOCK_STALE_MIN))
    log("save ready after {:.1f} min".format((time.time() - t0) / 60.0))


def set_last(src, name, val, base):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(base, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def run_push(direction):
    base = "pushover_oop_3dec_v5_zup.py"
    out = "pushover_results_v5zup_pk2_ft0p05_{}_TopMason".format(POCKET)
    side = "pos" if direction > 0 else "neg"
    csvp = os.path.join(out, "pushover_{}.csv".format(side))
    endsav = os.path.join(out, "pushover_{}_end".format(side))
    if os.path.isfile(endsav + ".sav") or os.path.isfile(endsav):
        log("{} already finished ({} exists) -- nothing to do".format(csvp, endsav))
        return
    if os.path.isfile(csvp):
        os.replace(csvp, csvp + ".partial")
        log("{} had no end save (stopped run): kept as .partial, pushing again".format(csvp))
    with open(base) as f:
        src = f.read()
    if "Z-UP" not in src:
        raise RuntimeError("{} is not the z-up driver".format(base))
    for k, v in (("RESTORE_SAV", SAVE), ("OUT_DIR", out), ("PUSH_DIR", direction), ("TOP_JOINT", "mason"),
                 ("TENSION_MAX_PA", 0.06e6)):
        src = set_last(src, k, v, base)
    log("{} push: {} -> {}".format("+" if direction > 0 else "-", SAVE, csvp))
    t0 = time.time()
    exec(compile(src, "{} [{} {}]".format(base, POCKET, side), "exec"),
         {"__name__": "__main__", "__file__": os.path.abspath(base)})
    log("done in {:.1f} min".format((time.time() - t0) / 60.0))


def run_stratF():
    base = "run_stratF_zup.py"
    key = "zup_pk2_ft0p05_" + POCKET
    entry = ('    "{}": ("{}.sav", "{}",\n'
             '                       (0.04e6, 0.06e6), (0.23455, 2.11091), (0.32250, 2.19886), 0.091),\n').format(key, SAVE, DAT)
    with open(base) as f:
        src = f.read()
    if '"{}"'.format(key) not in src:
        m = re.search(r"^_BUILDS[ \t]*=[ \t]*\{[^\n]*\n", src, re.M)
        if not m:
            raise RuntimeError("{}: '_BUILDS = {{' not found".format(base))
        src = src[:m.end()] + entry + src[m.end():]
    for k, v in (("BUILD", key), ("DAMP_PCT", DAMP_PCT), ("DRY_RUN", False)):
        src = set_last(src, k, v, base)
    log("Strategy F runs 1-25, {} % Maxwell, build {}".format(DAMP_PCT, key))
    exec(compile(src, "{} [{}]".format(base, key), "exec"), {"__name__": "__main__", "__file__": os.path.abspath(base)})


log("=" * 60)
ensure_save()
if JOB == "push+":
    run_push(+1)
elif JOB == "push-":
    run_push(-1)
else:
    run_stratF()
