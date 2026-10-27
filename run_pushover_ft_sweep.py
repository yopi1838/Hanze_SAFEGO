# -*- coding: ascii -*-
"""
run_pushover_ft_sweep.py -- positive pushovers on builds that differ ONLY in the joint tensile strength f_t.

Purpose: choose f_t from a curve instead of a single value. For each f_t in FT_LIST it
  1. writes ANALYSIS_PART_I_MASON_ft_<tag>.dat from the f_t = 0 build file (TEMPLATE): identical, except that right
     after "call 'mat_case.dat'" it sets   tension = f_t,  cohesion = COHESION,  Gt = 0.025 (2 f_t)^0.7 x 1000
     (mat_case.dat's own G_I expression), and it saves Part_I_MASON_v8_LargeStrain_springforce_ft_<tag>.sav;
     mat_case.dat itself is not edited;
  2. builds that save (program call) if it does not exist yet;
  3. runs pushover_oop_3dec_v5.py POSITIVE only from it (config lines rewritten in memory, as the batch does) into
     pushover_results_v5_ft_<tag>_TopMason/;
  4. adds a line to ft_sweep/ft_sweep_summary.csv: initial stiffness K_el (|a| m_tot against d over 0.1-1.5 mm, as
     the post-processing), T1 = 2 pi sqrt(L / K_el) with L = 1200 kg, peak |a| m_tot and its displacement,
     the force at 20 / 40 mm. The existing f_t = 0 and f_t = 0.15 positive pushovers are added for reference.
Everything except f_t is the f_t = 0 build: cohesion 0.2 MPa, same stiffness, joists, loads; the material table for
new contacts stays tension 0 / cohesion 0. A finished f_t (its pushover_pos.csv exists) is skipped; FORCE redoes it.

Targets (US-1, corrected base shear of 26 Oct): T1 about 0.091 s, peak about 10.6 kN total (8.7 kN wall part).
Expected: below about f_t = 0.06 MPa the cracking force stays under the rocking capacity, so the peak should stay
near the f_t = 0 value (10.9 kN) while T1 shortens; above it the cracking peak takes over (about 11 kN at 0.075,
13 kN at 0.1 -- Graziotti formula scaled to the DEM, to be checked here).

USAGE (3DEC 9, working folder = the one with ANALYSIS_PART_I_MASON.dat, GEO_WALL.sav, mat_case.dat and the
pushover script; needs the licence):
    DRY_RUN = True : writes the build .dat files, patches and compiles the pushover; no model work (plain python ok)
    DRY_RUN = False: program call 'run_pushover_ft_sweep.py'
"""
import os, re, csv, math, time

# ============================ CONFIG =================================
DRY_RUN       = False
FT_LIST       = (0.05e6, 0.075e6, 0.10e6)      # Pa, the builds to make and push
COHESION      = 0.2e6                          # Pa, kept at the f_t = 0 build's value: f_t is the only change
TEMPLATE      = "ANALYSIS_PART_I_MASON.dat"    # the f_t = 0 build (saves Part_I_MASON_v8_LargeStrain_springforce)
TEMPLATE_SAVE = "Part_I_MASON_v8_LargeStrain_springforce"
PUSHOVER_SCRIPT = "pushover_oop_3dec_v4.py"
PO_COMMON     = dict(TOP_JOINT="mason", PUSH_DIR=+1)
FORCE         = []                             # f_t values (Pa) to redo even if their CSV exists
REFERENCE     = [(0.0, "pushover_results_v5_ft0_masonrypush_TopMason/pushover_pos.csv"),
                 (0.15e6, "pushover_results_v5_ft0p15_masonrypush_TopMason/pushover_pos.csv")]
L_KG          = 1200.0                         # effective mass for T1 (measured mode shape)
K_EL_RANGE_MM = (0.1, 1.5)
OUT           = "ft_sweep"
# ======================================================================

try:
    import itasca as it
except ImportError:
    it = None


def tag(ft):
    return "{:.3f}".format(ft / 1e6).replace(".", "p")          # 0.05e6 -> '0p050'


def log(s=""):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), s)
    print(line)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "ft_sweep_log.txt"), "a") as f:
        f.write(line + "\n")


def build_dat(ft):
    with open(TEMPLATE) as f:
        src = f.read()
    calls = list(re.finditer(r"^[ \t]*call[ \t]+['\"]?mat_case\.dat['\"]?[^\n]*$", src, re.M | re.I))
    saves = list(re.finditer(r"^[ \t]*model[ \t]+save[ \t]+['\"]([^'\"]+)['\"][^\n]*$", src, re.M | re.I))
    if len(calls) != 1:
        raise RuntimeError("{}: expected one line \"call 'mat_case.dat'\", found {}".format(TEMPLATE, len(calls)))
    if not saves:
        raise RuntimeError("{}: no \"model save '...'\" line".format(TEMPLATE))
    sv = saves[-1]                                    # the final save of the build
    if sv.group(1) != TEMPLATE_SAVE:
        log("  note: {} saves '{}' (not '{}'); the sweep replaces that last save line".format(TEMPLATE, sv.group(1), TEMPLATE_SAVE))
    later = re.findall(r"^[ \t]*\[?[ \t]*(?:global[ \t]+)?(tension|cohesion)[ \t]*=", src[calls[0].end():sv.start()], re.M)
    if later:
        log("  ! {} assigns {} after the mat_case call: those lines would override the sweep -- check the template".format(
            TEMPLATE, sorted(set(later))))
    new_save = "{}_ft_{}".format(TEMPLATE_SAVE, tag(ft))
    over = ("\n; ---- run_pushover_ft_sweep: f_t = {:.3f} MPa, cohesion {:.3f} MPa (overrides mat_case.dat) ----\n"
            "[global tension  = {:g}]\n[global cohesion = {:g}]\n"
            "[global Gt = 0.025*((2*tension*1e-6)^0.7)*1000]\n"
            "[io.out('ft sweep build: tension ' + string(tension) + ' Pa, cohesion ' + string(cohesion) + ' Pa, G_I ' + string(Gt) + ' N/m')]"
            ).format(ft / 1e6, COHESION / 1e6, ft, COHESION)
    src = (src[:calls[0].end()] + over + src[calls[0].end():sv.start()] + "model save '{}'".format(new_save)
           + src[sv.end():])
    head = ("; ANALYSIS_PART_I_MASON_ft_{t}.dat -- written by run_pushover_ft_sweep.py from {tpl}.\n"
            "; Identical to it except f_t = {ft:.3f} MPa and cohesion {c:.3f} MPa for every joint, and the save name.\n"
            ).format(t=tag(ft), tpl=TEMPLATE, ft=ft / 1e6, c=COHESION / 1e6)
    dat = "ANALYSIS_PART_I_MASON_ft_{}.dat".format(tag(ft))
    with open(dat, "w") as f:
        f.write(head + src)
    return dat, new_save


def patch(src, over):
    for name, val in over.items():
        m = re.search(r"^{}[ \t]*=.*$".format(re.escape(name)), src, re.M)
        if m is None:
            raise RuntimeError("{}: config line '{} = ...' not found".format(PUSHOVER_SCRIPT, name))
        src = src[:m.start()] + "{} = {!r}   # [ft_sweep override]".format(name, val) + src[m.end():]
    return src


def pushover_src(save, out_dir, ft):
    with open(PUSHOVER_SCRIPT) as f:
        src = f.read()
    over = dict(PO_COMMON, RESTORE_SAV=save, OUT_DIR=out_dir, TENSION_MAX_PA=ft + 0.01e6)
    return patch(src, over), over


def summarise(ft, path):
    if not os.path.isfile(path):
        return None
    head = open(path).readline()
    try:
        m_tot = float(head.split("m_tot")[1].split("kg")[0])
    except Exception:
        m_tot = 1308.3
    rows = [r for r in csv.DictReader(l for l in open(path) if not l.lstrip().startswith("#"))]
    rows = [r for r in rows if r.get("converged", "1") == "1"]
    if not rows:
        return None
    d = [abs(float(r["d_ctrl_mm"])) for r in rows]
    F = [abs(float(r["a_mps2"])) * m_tot / 1e3 for r in rows]
    sel = [(x, y) for x, y in zip(d, F) if K_EL_RANGE_MM[0] <= x <= K_EL_RANGE_MM[1]]
    K = sum(x * y for x, y in sel) / sum(x * x for x, _ in sel) if len(sel) >= 2 else float("nan")   # kN/mm, through 0
    T1 = 2 * math.pi * math.sqrt(L_KG / (K * 1e6)) if K == K and K > 0 else float("nan")
    j = max(range(len(F)), key=lambda i: F[i])

    def at(x0):
        for i in range(1, len(d)):
            if d[i - 1] <= x0 <= d[i]:
                return F[i - 1] + (F[i] - F[i - 1]) * (x0 - d[i - 1]) / (d[i] - d[i - 1] or 1)
        return float("nan")
    return dict(ft_MPa=ft / 1e6, K_el_kN_per_mm=K, T1_s=T1, F_peak_kN=F[j], d_at_peak_mm=d[j],
                F_20mm_kN=at(20.0), F_40mm_kN=at(40.0), d_last_mm=d[-1], n_points=len(rows), csv=path)


def write_summary(rows):
    keys = ["ft_MPa", "K_el_kN_per_mm", "T1_s", "F_peak_kN", "d_at_peak_mm", "F_20mm_kN", "F_40mm_kN", "d_last_mm",
            "n_points", "csv"]
    p = os.path.join(OUT, "ft_sweep_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in sorted(rows, key=lambda r: r["ft_MPa"]):
            w.writerow({k: ("{:.4f}".format(v) if isinstance(v, float) else v) for k, v in r.items()})
    log("summary -> {}".format(p))
    for r in sorted(rows, key=lambda r: r["ft_MPa"]):
        log("  f_t {:.3f} MPa: K_el {:.2f} kN/mm, T1 {:.4f} s, peak {:.2f} kN at {:.1f} mm, F(20) {:.2f}, F(40) {:.2f} kN".format(
            r["ft_MPa"], r["K_el_kN_per_mm"], r["T1_s"], r["F_peak_kN"], r["d_at_peak_mm"], r["F_20mm_kN"], r["F_40mm_kN"]))


def main():
    log("=" * 70)
    log("run_pushover_ft_sweep  {}  f_t = {} MPa, cohesion {} MPa".format(
        "(DRY RUN)" if DRY_RUN else "", ", ".join("{:.3f}".format(x / 1e6) for x in FT_LIST), COHESION / 1e6))
    for f in (TEMPLATE, PUSHOVER_SCRIPT):
        if not os.path.isfile(f):
            raise RuntimeError("{} not found in {}".format(f, os.getcwd()))
    rows = []
    for ft in FT_LIST:
        log("-" * 70)
        dat, save = build_dat(ft)
        out_dir = "pushover_results_v5_ft_{}_TopMason".format(tag(ft))
        csvp = os.path.join(out_dir, "pushover_pos.csv")
        src, over = pushover_src(save, out_dir, ft)
        compile(src, PUSHOVER_SCRIPT, "exec")
        log("f_t {:.3f} MPa: {} -> {}.sav; pushover {}".format(ft / 1e6, dat, save,
                                                             ", ".join("{}={!r}".format(k, v) for k, v in over.items())))
        if DRY_RUN:
            continue
        if it is None:
            raise RuntimeError("itasca not available: run inside 3DEC (or DRY_RUN = True)")
        if os.path.isfile(csvp) and ft not in FORCE:
            log("  {} exists -- skipped (add {:g} to FORCE to redo)".format(csvp, ft))
        else:
            try:
                if not (os.path.isfile(save) or os.path.isfile(save + ".sav")):
                    t0 = time.time()
                    it.command("python-reset-state false")
                    it.command("program call '{}'".format(dat))
                    it.command("python-reset-state false")
                    log("  built {}.sav in {:.1f} min".format(save, (time.time() - t0) / 60))
                it.command("python-reset-state false")
                it.command("model restore '{}'".format(save))
                it.command("python-reset-state false")
                ftv = float(it.fish.get("tension"))
                if abs(ftv - ft) > 1.0:
                    raise RuntimeError("{}.sav holds tension {:.0f} Pa, expected {:.0f}".format(save, ftv, ft))
                if os.path.isfile(csvp):
                    os.rename(csvp, csvp.replace(".csv", "_prev_{}.csv".format(time.strftime("%Y%m%d_%H%M%S"))))
                t0 = time.time()
                try:
                    exec(compile(src, PUSHOVER_SCRIPT, "exec"), {"__name__": "__main__", "__file__": os.path.abspath(PUSHOVER_SCRIPT)})
                except SystemExit:
                    pass
                log("  pushover done in {:.1f} min".format((time.time() - t0) / 60))
            except Exception as e:
                log("  f_t {:.3f} FAILED: {} -- next one".format(ft / 1e6, e))
        r = summarise(ft, csvp)
        if r:
            rows.append(r)
    if DRY_RUN:
        log("dry run: build files written and the pushover compiles; set DRY_RUN = False to run")
        return
    for ft, p in REFERENCE:
        r = summarise(ft, p)
        if r:
            rows.append(r)
        else:
            log("  reference {} not found -- not in the summary".format(p))
    write_summary(rows)
    log("finished")


main()