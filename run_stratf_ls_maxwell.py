# -*- coding: ascii -*-
"""
run_stratF_LS_MAXWELL.py
    Strategy F (FULL records), large strain, Maxwell damping 1.5% over 1-40 Hz.

Run inside 3DEC:
    python-reset-state false
    call 'run_stratF_LS_MAXWELL.py'

Wrapper, not a copy: writes CONFIG to stratF_overrides.json (the file the
Strategy F driver reads) and then executes strategy_F_full_3dec.py. Its
companion differs only in MAXWELL_CMD, CASE_ID and OUT_DIR.
"""
import os, json, shutil

DRIVER   = "strategy_F_full_3dec.py"
OVR_FILE = "stratF_overrides.json"

CONFIG = {
    "CASE_ID":    "F_LS_MAXWELL_1p5",
    "OUT_DIR":    "stratF_full_LS_MAXWELL_1p5",
    "START_SAVE": "Part_I_MASON_v7.sav",
    "RUN_FROM": 1,
    "RUN_TO": 25,

    "LARGE_STRAIN": "on",
    # 1.5% over 1-40 Hz, same coefficients as the pulse drivers.
    #
    # CAUTION specific to full records: JModelMason::getMaxShearStiffness()
    # returns a constant ks_ regardless of contact state, so if 3DEC applies
    # the Maxwell shear dashpot to open subcontacts a separated joint carries
    # viscous shear across the gap. A full record runs far longer than a
    # 3-cycle pulse, so a cracked band spends much more time open and the
    # artefact has more time to accumulate. Log energy-shear and compare
    # against the NODAMP companion before trusting this run.
    "MAXWELL_CMD": ("block mech damp maxwell "
                    "0.0120 1.0006 0.0093 9.0193 0.0104 40.0000"),
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