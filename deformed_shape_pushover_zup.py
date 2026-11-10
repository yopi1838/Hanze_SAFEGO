# -*- coding: ascii -*-
"""
deformed_shape_pushover_zup.py -- deformed shape and joint damage of the z-up pushovers at their end saves
(80 mm: pushover_<side>_end.sav; continued to F = 0: pushover_<side>_end2.sav), 10 Nov 2026.

Reuses deformed_shape_stratF_zup.py unchanged: its functions are loaded from that file, and only the save list,
the reference save and the output folder are set here. The displacement is taken relative to the pocket variant's
own build save (the wall before the push), so it is the pushover deformation itself.

    code 1  pushover_pos_end   (+ push, 80 mm)          code 3  pushover_neg_end   (- push, 80 mm)
    code 2  pushover_pos_end2  (+ push, continued)      code 4  pushover_neg_end2  (- push, continued)
The code is the "run" number in the output files (run_01_courses.csv = + push at 80 mm, ...); saves.csv lists them.

USAGE  3DEC 9, Hanze_SAFEGO folder, no other model open:  program call 'deformed_shape_pushover_zup.py'
       outside 3DEC (re-process the dumps):                python deformed_shape_pushover_zup.py --post
OUTPUT <pushover folder>/deformed_shape/  (same files as the Strategy F version, plus saves.csv)
"""
import os, sys, csv

CASES = [   # pushover folder, build save the push started from
    ("pushover_results_v5zup_pk2_ft0p05_unclamp_TopMason", "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_unclamp"),
    ("pushover_results_v5zup_pk2_ft0p05_TopMason",         "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05"),
]
CODES = [(1, "pushover_pos_end"), (2, "pushover_pos_end2"), (3, "pushover_neg_end"), (4, "pushover_neg_end2")]

src = open("deformed_shape_stratF_zup.py").read()
cut = src.rstrip().rfind("\nmain()")
if cut < 0:
    raise RuntimeError("deformed_shape_stratF_zup.py: the closing 'main()' call was not found")
lib_src = src[:cut]
for a, b in (("Strategy F residual state (4 s after each record)", "pushover end state (1: + at 80 mm, 2: + continued, 3: - at 80 mm, 4: - continued)"),
             ("residual OOP displacement rel. to base (mm)", "OOP displacement rel. to base (mm)"),
             ("Residual deformed shape", "Deformed shape at the save"),
             ("Bed-joint damage by run", "Bed-joint damage by save code"),
             ("Run {}: joints against height", "Code {}: joints against height")):
    lib_src = lib_src.replace(a, b)      # figure labels only; nothing else in the functions changes

for folder, ref in CASES:
    if not os.path.isdir(folder):
        print("skip {}: folder not found".format(folder)); continue
    have = [(c, os.path.join(folder, n + ".sav")) for c, n in CODES if os.path.isfile(os.path.join(folder, n + ".sav"))]
    post_only = "--post" in sys.argv
    if not have and not post_only:
        print("skip {}: no pushover end save".format(folder)); continue
    ns = {"__name__": "dfm_lib", "__file__": os.path.abspath("deformed_shape_stratF_zup.py")}
    exec(compile(lib_src, "deformed_shape_stratF_zup.py [functions]", "exec"), ns)
    ns["RUN_DIR"] = folder
    ns["OUT_DIR"] = os.path.join(folder, "deformed_shape")
    ns["REF_SAV"] = ref
    ns["RUNS"] = None
    ns["find_runs"] = (lambda h=have: list(h))
    os.makedirs(ns["OUT_DIR"], exist_ok=True)
    with open(os.path.join(ns["OUT_DIR"], "saves.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["code", "save", "reference"]); [w.writerow([c, s, ref]) for c, s in have]
    print("=" * 80)
    print("{}: {} end saves {}, reference {}".format(folder, len(have), [os.path.basename(s) for _, s in have], ref))
    ns["main"]()
