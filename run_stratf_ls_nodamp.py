# -*- coding: ascii -*-
"""
run_stratF_LS_NODAMP.py
    Strategy F (FULL records), large strain, no damping.

Run inside 3DEC:
    python-reset-state false
    call 'run_stratF_LS_NODAMP.py'

Wrapper, not a copy: writes CONFIG to stratF_overrides.json (the file the
Strategy F driver reads) and then executes strategy_F_full_3dec.py. Its
companion differs only in MAXWELL_CMD, CASE_ID and OUT_DIR.
"""
import os, json, shutil

DRIVER   = "strategy_F_full_3dec_LS.py"
OVR_FILE = "stratF_overrides_g00.json"

CONFIG = {
    "CASE_ID":    "F_LS_NODAMP",
    "OUT_DIR":    "stratF_full_LS_NODAMP",
    "START_SAVE": "Part_I_MASON_v7.sav",
    "RUN_FROM": 1,
    "RUN_TO": 25,

    "LARGE_STRAIN": "on",
    # No damping. Contact dissipation only -- the single-variable reference
    # for the damped run.
    "MAXWELL_CMD": "",
    "DAMP_RATIO": 0.0,
    "DAMP_TYPE": "",

    # Moshfeghi et al. (2024) sec. 3.2. In Strategy F the record's velocity
    # table returns to zero at t = dur, so the table-stop instant is exact and
    # no detection rule is needed.
    "PERIOD_ID_METHOD": "experiment",
    "T1_init": 0.0948,
    "TAIL_SEC_RECORD": 2.5,
}

if not os.path.isfile(DRIVER):
    raise SystemExit("Cannot find {} in {}".format(DRIVER, os.getcwd()))

if os.path.exists(OVR_FILE):
    try:
        existing = json.load(open(OVR_FILE))
    except Exception:
        existing = None
    if existing != CONFIG:
        shutil.copyfile(OVR_FILE, OVR_FILE + ".bak")
        print("NOTE: existing {} differed and was backed up to {}.bak"
              .format(OVR_FILE, OVR_FILE))

with open(OVR_FILE, "w") as f:
    json.dump(CONFIG, f, indent=2)

print("=" * 70)
print("  {}".format(CONFIG["CASE_ID"]))
print("=" * 70)
print("  geometry   : large-strain {}".format(CONFIG["LARGE_STRAIN"]))
print("  damping    : {}".format(CONFIG["MAXWELL_CMD"] or "NONE"))
print("  period ID  : {}".format(CONFIG["PERIOD_ID_METHOD"]))
print("  output     : {}".format(CONFIG["OUT_DIR"]))
print("=" * 70)
print("  Check the driver's override banner below. Any key shown as")
print("  '! IGNORED' is NOT applied.")
print("=" * 70)

exec(compile(open(DRIVER).read(), DRIVER, "exec"), globals())