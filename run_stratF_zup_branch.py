# -*- coding: ascii -*-
"""
run_stratF_zup_branch.py -- re-run Strategy F from a FINISHED run of an existing case with a different damping
(2 Nov 2026). Default: the pk2 f_t 0.05, 5 % case, restarted from its run-17 save, run 18 only, at 1.5 % Maxwell.
Same build, same state at the start of run 18, same record and tail: the damping is the only difference, so the
run-18 peaks (5 % in the source case vs 1.5 % here) say directly how much the damping holds the response down.

HOW: copies <SRC_DIR>/stratC_run_<FROM_RUN>.sav into a new folder, writes a checkpoint there that says "run FROM_RUN
done" for the new case/damping, and runs run_stratF_zup.py (NOT edited; read and patched in memory) with
BUILD / DAMP_PCT / OUT_DIR / CASE_ID / RUN_TO set here. The driver then resumes: restores the copied save, applies
the new Maxwell set (after zeroing local/global damping) and runs FROM_RUN+1 ... RUN_TO.
OUT_DIR: stratF_branch_r<FROM_RUN>_MAXWELL<pct x10>_<BUILD>_tr

USAGE (3DEC 9, Hanze_SAFEGO folder):  program call 'run_stratF_zup_branch.py'
      (DRY_RUN = True: patch + compile + show what would be copied/written, no 3DEC needed)
Compare:  python damping_vs_period_check.py --model "5 %=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr" \
              --model "1.5 % from run 17=stratF_branch_r17_MAXWELL15_zup_pk2_ft0p05_tr" --runs 18 --show 18
          and deformed_shape_stratF_zup.py with RUN_DIR set to the branch folder (damage at 2.29 m etc.)
"""
import os, re, json, shutil

# ============================ CONFIG =================================
DRY_RUN  = False
BUILD    = "zup_pk2_ft0p05"
SRC_DIR  = "stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr"   # the finished case to branch from
FROM_RUN = 17                                              # its save stratC_run_17.sav = state before run 18
RUN_TO   = 18                                              # last run of the branch (18 = one run; 21 = to the end of HU12)
DAMP_PCT = 1.5                                             # a key of run_stratF_zup.py's MAXWELL_SETS
BASE     = "run_stratF_zup.py"
# ======================================================================
TAG = "# [run_stratF_zup_branch]"
OUT_DIR = "stratF_branch_r{:02d}_MAXWELL{:02d}_{}_tr".format(FROM_RUN, int(round(DAMP_PCT * 10)), BUILD)
CASE_ID = "F_branch_r{:02d}_MAXWELL{:02d}_{}_tr".format(FROM_RUN, int(round(DAMP_PCT * 10)), BUILD)


def set_last(src, name, val):
    hits = list(re.finditer(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M))
    if not hits:
        raise RuntimeError("{}: no config line '{} = ...'".format(BASE, name))
    m = hits[-1]
    return src[:m.start()] + "{} = {!r}   {}".format(name, val, TAG) + src[m.end():]


def rep1(src, old, new):
    k = src.count(old)
    if k != 1:
        raise RuntimeError("{}: expected the line {!r} once, found {} -- update this runner".format(BASE, old.strip(), k))
    return src.replace(old, new)


def main():
    with open(BASE) as f:
        src = f.read()
    src = set_last(src, "BUILD", BUILD)
    src = set_last(src, "DAMP_PCT", DAMP_PCT)
    src = set_last(src, "DRY_RUN", DRY_RUN)
    src = rep1(src, 'OUT_DIR="stratF_full_LS_{}_{}_tr".format(_TAGR, BUILD),', "OUT_DIR={!r},   {}".format(OUT_DIR, TAG))
    src = rep1(src, 'CASE_ID="F_LS_{}_{}_tr".format(_TAGR, BUILD),', "CASE_ID={!r},   {}".format(CASE_ID, TAG))
    src = rep1(src, "RUN_TO=25,", "RUN_TO={},   {}".format(RUN_TO, TAG))
    # the patched config, to write a checkpoint the driver accepts (same case, damping, geometry)
    head = src.split('TAG = "# [run_stratF_zup]"')[0]
    ns = {}
    exec(compile(head, BASE + " [config]", "exec"), ns)
    ov = ns["OVERRIDES"]
    if ov["OUT_DIR"] != OUT_DIR or ov["RUN_TO"] != RUN_TO or ov.get("RUN_FROM", 1) != 1:
        raise RuntimeError("override check failed: {}".format({k: ov[k] for k in ("OUT_DIR", "RUN_FROM", "RUN_TO")}))
    sav_src = os.path.join(SRC_DIR, "stratC_run_{:02d}.sav".format(FROM_RUN))
    sav_dst = os.path.join(OUT_DIR, "stratC_run_{:02d}.sav".format(FROM_RUN))
    ck_src = os.path.join(SRC_DIR, "stratC_checkpoint.json")
    ck_dst = os.path.join(OUT_DIR, "stratC_checkpoint.json")
    print("[branch] {} run {:02d} -> {} : runs {:02d}-{:02d} at {} % ({})".format(
        SRC_DIR, FROM_RUN, OUT_DIR, FROM_RUN + 1, RUN_TO, DAMP_PCT, ov["MAXWELL_CMD"]))
    if not os.path.isfile(sav_src):
        raise RuntimeError("{} not found".format(sav_src))
    T_cur, summary = ns.get("T1_INIT", 0.09), []
    if os.path.isfile(ck_src):
        with open(ck_src) as f:
            st = json.load(f)
        T_cur = st.get("T_current", T_cur)
        summary = [s for s in st.get("summary", []) if int(s.get("run", 0)) <= FROM_RUN]
        if int(st.get("last_completed_run", 0)) < FROM_RUN:
            raise RuntimeError("{} has finished only up to run {} (< {})".format(SRC_DIR, st["last_completed_run"], FROM_RUN))
    state = {"last_completed_run": FROM_RUN, "T_current": T_cur, "summary": summary, "case": ov["CASE_ID"],
             "start_save": ov["START_SAVE"], "large_strain": ov["LARGE_STRAIN"], "maxwell": ov["MAXWELL_CMD"],
             "branched_from": sav_src}
    if DRY_RUN:
        print("[branch] DRY_RUN: would copy {} -> {} and write {} (last_completed_run {}, maxwell '{}')".format(
            sav_src, sav_dst, ck_dst, FROM_RUN, ov["MAXWELL_CMD"]))
    else:
        os.makedirs(OUT_DIR, exist_ok=True)
        if not os.path.isfile(sav_dst):
            shutil.copy2(sav_src, sav_dst)
            print("[branch] copied {} -> {}".format(sav_src, sav_dst))
        if not os.path.isfile(ck_dst):
            with open(ck_dst, "w", newline="\n") as f:
                json.dump(state, f, indent=2)
            print("[branch] checkpoint written: run {:02d} done, damping '{}'".format(FROM_RUN, ov["MAXWELL_CMD"]))
        else:
            print("[branch] checkpoint exists in {}: the driver resumes after its last finished run".format(OUT_DIR))
    code = compile(src, "{} [branch {}]".format(BASE, OUT_DIR), "exec")
    exec(code, {"__name__": "__main__", "__file__": os.path.abspath(BASE)})


main()
