# -*- coding: ascii -*-
"""
batch_ft0p15_pipeline.py -- one 3DEC session, six stages, run back to back on a single licence.

    0. pushover, POSITIVE, f_t = 0     (restores Part_I_MASON_v8_LargeStrain_springforce)
    1. pushover, NEGATIVE, f_t = 0     (same save)
    2. pushover, POSITIVE, f_t = 0.15  (restores Part_I_MASON_v8_LargeStrain_springforce_ft_0p15)
    3. pushover, NEGATIVE, f_t = 0.15  (same save)
    4. bilinear_idealisation.py on the two f_t = 0.15 branches (report; optionally the periods for stage 5)
    5. route2_Tm_ihsan_full25_3dec.py on the f_t = 0.15 save, both pushover branches, Maxwell 5 %, XI 0.05
       (a preflight-only pass first, logged, then the 25 runs)

HOW IT WORKS
    Each stage EXECUTES the existing script unchanged except for its CONFIG lines, which are rewritten in memory
    (NAME = value, matched at the start of a line; a name that is not found stops the batch, so a typo cannot
    silently run the default). No file on disk is edited. Every stage restores its own save, so the stages do not
    share model state.

    Pushovers use TOP_JOINT = "mason" (the build's joint). All four pushovers also write snapshot saves at SNAP_MM on
    the descending branch (for deformed_shape_3dec.py), inserted after the servo point is logged.
    Stage 0 writes into PO_DIR_FT0; if a pushover_pos.csv is already there it is renamed to
    pushover_pos_prev_<time>.csv first, so an earlier positive f_t = 0 run is kept.

RESTARTS
    <BATCH_DIR>/batch_state.json records the stages that finished. Run the file again and finished stages are
    skipped (set FORCE_STAGES to redo one). Stage 5 also resumes by itself from its own checkpoint in its OUT_DIR.

BEFORE STARTING
    * build Part_I_MASON_v8_LargeStrain_springforce_ft_0p15.sav (mat_case.dat: tension 0.15e6, cohesion ~0.3e6);
    * set DRY_RUN = True and run it once (inside 3DEC, or with plain python): it checks every file and every
      override and prints the patched config of each stage, without touching the model;
    * then DRY_RUN = False, and:   program call 'batch_ft0p15_pipeline.py'   (or run it from the Python pane).

The save checks at the start restore both saves once and read the joint tension stored by mat_case.dat (FISH
'tension'): ~0 for the f_t = 0 save, ~0.15 MPa for the new one. A wrong save stops the batch before any analysis.
"""
import os, re, sys, json, time, csv, math, traceback

# ============================ CONFIG =================================
DRY_RUN       = False
BATCH_DIR     = "batch_ft0p15"
FORCE_STAGES  = []                 # e.g. [4, 5] to redo stages that are marked done

SAVE_FT0      = "Part_I_MASON_v8_LargeStrain_springforce"          # f_t = 0 build (stage 1)
SAVE_FT015    = "Part_I_MASON_v8_LargeStrain_springforce_ft_0p15"  # f_t = 0.15 MPa build (stages 2, 3, 5)
FT0_MAX_PA    = 0.05e6             # check: tension in SAVE_FT0 below this
FT015_RANGE   = (0.10e6, 0.20e6)   # check: tension in SAVE_FT015 inside this range

PUSHOVER_SCRIPT = "pushover_oop_3dec_v4.py"
BILINvEAR_SCRIPT = "bilinear_idealisation.py"
R2TM_SCRIPT     = "r2m3_ihsan_full25.py"

PO_DIR_FT0    = "pushover_results_v5_ft0_masonrypush_TopMason"     # stages 0 and 1 (positive and negative f_t = 0)
PO_DIR_FT015  = "pushover_results_v5_ft0p15_masonrypush_TopMason"
SNAP_MM       = (5, 10, 20, 40)    # snapshot saves on every pushover's descending branch; () for none

BIL_DIR       = "bilinear_idealisation_ft0p15"
BIL_M_EFF     = 1200.0             # kg, the measured effective mass r2tm uses (the bilinear script's default is 1635)

R2TM_OUT_DIR  = "r2tm_TpTm_full25_LS_maxwell5_ft0p15"
MAXWELL_5PCT  = "block mech damp maxwell 0.0373 1.0005 0.0270 7.4315 0.0320 25.0004"   # 5 % flat 1-25 Hz
PERIODS_FROM  = "r2tm"             # "r2tm": r2tm computes T1/T_eff from the two CSVs (as in every previous run)
                                   # "bilinear": T1/T_eff of stage 4 (model average) passed as overrides
PUSHOVER_FORCE = "applied"        # force read from the pushover CSVs by stages 4 and 5:
                                   #   "applied":  |a| x m_tot (the push on the masonry at a settled servo point =
                                   #               the horizontal force the wall carries)
                                   #   "reaction": F_tot_kN as logged. F_tot sums only the SHEAR z-components at the
                                   #               top joint; once the top beam tilts, the top-joint normal force adds
                                   #               |N| sin(tilt) horizontally, which F_tot misses (+1.8 / +2.6 kN at
                                   #               80 mm on the 22 Oct elastic-top runs; deformed_shape_3dec.py, 23 Oct)
STOP_IF_FLAGGED = False            # True: stop before the 25 runs if the preflight flags any run

# per-stage overrides (the NAME = value lines rewritten in each script)
PO_COMMON = dict(TOP_JOINT="mason")
STAGES = {
    0: dict(kind="pushover", critical=False, title="pushover POS, f_t = 0",
            over=dict(PO_COMMON, RESTORE_SAV=SAVE_FT0, OUT_DIR=PO_DIR_FT0, PUSH_DIR=+1, TENSION_MAX_PA=FT0_MAX_PA),
            snaps=SNAP_MM),
    1: dict(kind="pushover", critical=False, title="pushover NEG, f_t = 0",
            over=dict(PO_COMMON, RESTORE_SAV=SAVE_FT0, OUT_DIR=PO_DIR_FT0, PUSH_DIR=-1, TENSION_MAX_PA=FT0_MAX_PA),
            snaps=SNAP_MM),
    2: dict(kind="pushover", critical=True, title="pushover POS, f_t = 0.15 MPa",
            over=dict(PO_COMMON, RESTORE_SAV=SAVE_FT015, OUT_DIR=PO_DIR_FT015, PUSH_DIR=+1, TENSION_MAX_PA=FT015_RANGE[1]),
            snaps=SNAP_MM),
    3: dict(kind="pushover", critical=True, title="pushover NEG, f_t = 0.15 MPa",
            over=dict(PO_COMMON, RESTORE_SAV=SAVE_FT015, OUT_DIR=PO_DIR_FT015, PUSH_DIR=-1, TENSION_MAX_PA=FT015_RANGE[1]),
            snaps=SNAP_MM),
    4: dict(kind="bilinear", critical=(PERIODS_FROM == "bilinear"), title="bilinear idealisation, f_t = 0.15 MPa"),
    5: dict(kind="r2tm", critical=True, title="Route 2 (T_p / T_m), f_t = 0.15 MPa, Maxwell 5 %"),
}
POS_CSV = os.path.join(PO_DIR_FT015, "pushover_pos.csv")
NEG_CSV = os.path.join(PO_DIR_FT015, "pushover_neg.csv")
R2TM_OVER = dict(OUT_DIR=R2TM_OUT_DIR, BASE_SAVE=SAVE_FT015 + ".sav", PUSHOVER_CSV=POS_CSV, PUSHOVER_CSV_NEG=NEG_CSV,
                 BACKBONE_COMBINE="mean", F_SOURCE="column", D_CUT_MM=None,
                 MAXWELL_CMD=MAXWELL_5PCT, XI=0.05, RESUME_AFTER=0, RESUME_SAVE="", STOP_AFTER=25,
                 VERIFY_BUILD=True, FREE_TB_ROTATIONS=False)
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None
os.makedirs(BATCH_DIR, exist_ok=True)
LOGF = os.path.join(BATCH_DIR, "batch_log.txt")
STATE = os.path.join(BATCH_DIR, "batch_state.json")


def log(msg=""):
    line = "[{}] {}".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line)
    with open(LOGF, "a") as f:
        f.write(line + "\n")


def load_state():
    if os.path.isfile(STATE):
        with open(STATE) as f:
            return json.load(f)
    return {"done": [], "results": {}}


def save_state(s):
    with open(STATE, "w") as f:
        json.dump(s, f, indent=1)


def save_exists(name):
    return os.path.isfile(name) or os.path.isfile(name + ".sav")


# ---------------- source patching ----------------
def patch(src, over, path):
    for name, val in over.items():
        pat = re.compile(r"^{}[ \t]*=.*$".format(re.escape(name)), re.M)
        m = pat.search(src)
        if m is None:
            raise RuntimeError("{}: config line '{} = ...' not found -- the script changed; update the batch".format(path, name))
        src = src[:m.start()] + "{} = {!r}   # [batch override]".format(name, val) + src[m.end():]
    return src


SNAP_ANCHOR = ('        rows.append(row("desc", a, d, Fb, Ft, Ftot, settled, m_part, a_pred, nch, dg))\n'
               '        flush()\n')
SNAP_CODE = ('        for _s in _BATCH_SNAP_MM:\n'
             '            if _s not in _batch_snapped and abs(d) >= _s:\n'
             '                it.command("model save \'{}\'".format(os.path.join(OUT_DIR, "pushover_{}_d{:02d}".format(\n'
             '                    "pos" if PUSH_DIR > 0 else "neg", int(_s))).replace("\\\\", "/")))\n'
             '                _batch_snapped.add(_s)\n'
             '                print("  snapshot save at d = {:.1f} mm (target {} mm)".format(d, _s))\n')


def add_snapshots(src, snaps, path):
    if not snaps:
        return src
    if src.count(SNAP_ANCHOR) != 1:
        log("  ! {}: snapshot anchor not found -- no snapshot saves (the run itself is unaffected)".format(path))
        return src
    src = src.replace(SNAP_ANCHOR, SNAP_ANCHOR + SNAP_CODE)
    return "_BATCH_SNAP_MM = {!r}\n_batch_snapped = set()\n".format(tuple(snaps)) + src


def run_script(path, over, snaps=(), expect_exit=False):
    with open(path) as f:
        src = f.read()
    src = add_snapshots(patch(src, over, path), snaps, path)
    for k, v in over.items():
        log("    {} = {!r}".format(k, v))
    if DRY_RUN:
        compile(src, path, "exec")
        log("    (dry run: patched source compiles; not executed)")
        return {}
    ns = {"__name__": "__main__", "__file__": os.path.abspath(path)}
    try:
        exec(compile(src, path, "exec"), ns)
    except SystemExit as e:
        if not expect_exit and e.code not in (0, None):
            raise RuntimeError("{} exited with {}".format(path, e.code))
    return ns


# ---------------- checks ----------------
def check_files():
    need = [PUSHOVER_SCRIPT, R2TM_SCRIPT, BILINEAR_SCRIPT]
    miss = [f for f in need if not os.path.isfile(f)]
    for s in (SAVE_FT0, SAVE_FT015):
        if not save_exists(s):
            miss.append(s + ".sav")
    if miss:
        raise RuntimeError("missing before start: " + ", ".join(miss))
    log("files: scripts and both saves found")


def check_save_tension(save, lo, hi):
    it.command("python-reset-state false")
    it.command("model restore '{}'".format(save))
    it.command("python-reset-state false")
    try:
        ft = float(it.fish.get("tension"))
    except Exception as e:
        raise RuntimeError("{}: FISH 'tension' not readable ({}) -- cannot verify the build".format(save, e))
    log("save {}: joint tension {:.3f} MPa (expected {:.2f}-{:.2f})".format(save, ft / 1e6, lo / 1e6, hi / 1e6))
    if not (lo <= ft <= hi):
        raise RuntimeError("{} has tension {:.3f} MPa, outside {:.2f}-{:.2f} MPa: wrong save".format(save, ft / 1e6, lo / 1e6, hi / 1e6))


def m_tot_of(path, default=1308.3):
    """masonry mass the v5 push normalises by, from the CSV's '# ... m_tot 1308.3 kg ...' line"""
    try:
        head = open(path).readline()
        return float(head.split("m_tot")[1].split("kg")[0])
    except Exception:
        return default


def strip_comment_csv(src, dst):
    """bilinear_idealisation.py reads the header from line 1; the v5 pushover CSV starts with a '#' line.
    With PUSHOVER_FORCE = "applied" the F_tot_kN column is replaced by sign(F_tot) |a| m_tot on converged rows
    (the logged value is kept as F_tot_logged_kN), so bilinear_idealisation.py reads the force the wall carries."""
    m_tot = m_tot_of(src)
    with open(src) as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines))
    fields = list(rows[0].keys()) + ["F_tot_logged_kN"] if rows else []
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            r["F_tot_logged_kN"] = r["F_tot_kN"]
            if PUSHOVER_FORCE == "applied" and r.get("converged", "1") == "1":
                Ft = float(r["F_tot_kN"])
                r["F_tot_kN"] = "{:.4f}".format(math.copysign(abs(float(r["a_mps2"])) * m_tot / 1e3, Ft if Ft else 1.0))
            w.writerow(r)
    log("  {} -> {} (force: {}{})".format(src, dst, PUSHOVER_FORCE,
                                         ", m_tot {:.1f} kg".format(m_tot) if PUSHOVER_FORCE == "applied" else ""))
    return dst


def csv_summary(path):
    rows = [r for r in csv.DictReader(l for l in open(path) if not l.lstrip().startswith("#"))]
    ok = [r for r in rows if r.get("converged", "1") == "1"]
    if not ok:
        return "no converged rows"
    F = [abs(float(r["F_tot_kN"])) for r in ok]; d = [abs(float(r["d_ctrl_mm"])) for r in ok]
    Fa = [abs(float(r["a_mps2"])) * m_tot_of(path) / 1e3 for r in ok]
    i, j = F.index(max(F)), Fa.index(max(Fa))
    return ("{} rows, peak |F_tot| {:.2f} kN at {:.1f} mm, peak |a| m_tot {:.2f} kN at {:.1f} mm; last point {:.1f} mm: "
            "|F_tot| {:.2f} / |a| m_tot {:.2f} kN").format(len(rows), F[i], d[i], Fa[j], d[j], d[-1], F[-1], Fa[-1])


# ---------------- stages ----------------
def stage_pushover(n, st):
    csv_old = os.path.join(st["over"]["OUT_DIR"], "pushover_{}.csv".format("pos" if st["over"]["PUSH_DIR"] > 0 else "neg"))
    if os.path.isfile(csv_old) and not DRY_RUN:
        keep = csv_old.replace(".csv", "_prev_{}.csv".format(time.strftime("%Y%m%d_%H%M%S")))
        os.rename(csv_old, keep)
        log("  existing {} kept as {}".format(csv_old, keep))
    run_script(PUSHOVER_SCRIPT, st["over"], st["snaps"])
    if DRY_RUN:
        return {}
    csvp = os.path.join(st["over"]["OUT_DIR"], "pushover_{}.csv".format("pos" if st["over"]["PUSH_DIR"] > 0 else "neg"))
    summ = csv_summary(csvp) if os.path.isfile(csvp) else "CSV missing"
    log("  {}: {}".format(csvp, summ))
    return {"csv": csvp, "summary": summ}


def stage_bilinear(n, st):
    os.makedirs(BIL_DIR, exist_ok=True)
    pos = strip_comment_csv(POS_CSV, os.path.join(BIL_DIR, "pushover_pos_clean.csv")) if os.path.isfile(POS_CSV) else None
    neg = strip_comment_csv(NEG_CSV, os.path.join(BIL_DIR, "pushover_neg_clean.csv")) if os.path.isfile(NEG_CSV) else None
    if pos is None and not DRY_RUN:
        raise RuntimeError("bilinear: positive f_t = 0.15 pushover CSV missing: " + POS_CSV)
    run_script(BILINEAR_SCRIPT, dict(MODEL_CSVS={"pos": pos, "neg": neg}, M_EFF=BIL_M_EFF, OUT_DIR=BIL_DIR))
    if DRY_RUN:
        return {}
    with open(os.path.join(BIL_DIR, "bilinear_periods.json")) as f:
        m = json.load(f)["model"]
    log("  bilinear (model, {}): K1 {:.2f} kN/mm, mu {:.2f}, T1 {:.4f} s, T_eff {:.4f} s (m_eff {:.0f} kg)".format(
        "+".join(m["branches_used"]), m["K1_kNmm"], m["mu"], m["T1_s"], m["Teff_s"], BIL_M_EFF))
    return {"T1_s": m["T1_s"], "Teff_s": m["Teff_s"], "K1_kNmm": m["K1_kNmm"], "mu": m["mu"]}


def stage_r2tm(n, st, state):
    over = dict(R2TM_OVER)
    if PUSHOVER_FORCE == "applied":
        over.update(F_SOURCE="applied", M_PUSH_KG=m_tot_of(POS_CSV), M_PUSH_KG_NEG=m_tot_of(NEG_CSV))
        log("  r2tm backbone force: |a| x m_tot ({:.1f} / {:.1f} kg)".format(over["M_PUSH_KG"], over["M_PUSH_KG_NEG"]))
    if PERIODS_FROM == "bilinear":
        b = state["results"].get("4", {})
        if "T1_s" not in b:
            raise RuntimeError("PERIODS_FROM = 'bilinear' but stage 4 has no periods")
        over.update(T1_S_OVERRIDE=b["T1_s"], TEFF_S_OVERRIDE=b["Teff_s"])
    # pass 1: preflight only (targets, amplitudes, flags) -- quick, logged before the long run
    log("  preflight pass")
    ns = run_script(R2TM_SCRIPT, dict(over, PREFLIGHT_ONLY=True), expect_exit=True)
    res = {}
    if not DRY_RUN:
        res = {"T1_s": ns.get("T1_S"), "Teff_s": ns.get("TEFF_S"), "periods": ns.get("PERIOD_NOTE"),
               "flagged": ns.get("FLAGGED", [])}
        log("  periods: {}".format(res["periods"]))
        log("  T1 {:.4f} s, T_eff {:.4f} s; flagged runs: {}".format(res["T1_s"], res["Teff_s"], res["flagged"] or "none"))
        if res["flagged"] and STOP_IF_FLAGGED:
            raise RuntimeError("preflight flagged runs {} and STOP_IF_FLAGGED is True".format(res["flagged"]))
    # pass 2: the 25 runs (resumes from its own checkpoint if this stage was interrupted)
    log("  dynamic pass (25 runs)")
    run_script(R2TM_SCRIPT, dict(over, PREFLIGHT_ONLY=False))
    return res


# ============================ MAIN ===================================
def main():
    log("=" * 78)
    log("batch_ft0p15_pipeline  {}".format("(DRY RUN)" if DRY_RUN else ""))
    state = load_state()
    for k in FORCE_STAGES:
        if k in state["done"]:
            state["done"].remove(k)
    check_files()
    if it is None and not DRY_RUN:
        raise RuntimeError("itasca not available: run inside 3DEC (or set DRY_RUN = True)")
    if it is not None:
        it.command("python-reset-state false")
        it.command("program automatic-model-save active off")
        if not DRY_RUN:
            check_save_tension(SAVE_FT0, -1.0, FT0_MAX_PA)
            check_save_tension(SAVE_FT015, *FT015_RANGE)
    for n in sorted(STAGES):
        st = STAGES[n]
        if n in state["done"]:
            log("stage {} ({}) already done -- skipped".format(n, st["title"])); continue
        log("-" * 78)
        log("stage {}: {}".format(n, st["title"]))
        t0 = time.time()
        try:
            if st["kind"] == "pushover":
                res = stage_pushover(n, st)
            elif st["kind"] == "bilinear":
                res = stage_bilinear(n, st)
            else:
                res = stage_r2tm(n, st, state)
        except Exception as e:
            log("stage {} FAILED after {:.1f} min: {}".format(n, (time.time() - t0) / 60, e))
            with open(LOGF, "a") as f:
                traceback.print_exc(file=f)
            if st["critical"]:
                log("stage {} is critical: batch stopped. Fix it and run the batch again (finished stages are skipped).".format(n))
                raise
            log("stage {} is not critical: continuing".format(n))
            continue
        log("stage {} done in {:.1f} min".format(n, (time.time() - t0) / 60))
        if not DRY_RUN:
            state["done"].append(n); state["results"][str(n)] = res; save_state(state)
    log("batch finished. Next: python postprocess_r2tm_full25.py --run-dir {} --control-dir <Strategy F folder>".format(R2TM_OUT_DIR))


main()
