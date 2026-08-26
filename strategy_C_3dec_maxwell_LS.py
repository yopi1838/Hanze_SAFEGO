# -*- coding: ascii -*-
"""
run_stratC_LS_MAXWELL.py
    Strategy C, LARGE STRAIN, Maxwell damping 1.5% over 1-40 Hz,
    experiment-matched period identification.

Run inside 3DEC:
    python-reset-state false
    call 'run_stratC_LS_MAXWELL.py'

WHY THIS IS A WRAPPER AND NOT A COPY OF THE DRIVER
--------------------------------------------------
The physics lives in strategy_C_3dec_nodamp.py. Copying 1100 lines to change
three constants creates a second copy that drifts the first time either is
edited -- and the entire value of this pair of runs is that they differ in
EXACTLY ONE thing. Keeping one driver makes that guarantee structural rather
than a promise. The companion script run_stratC_LS_NODAMP.py is identical to
this file except for MAXWELL_CMD, CASE_ID and OUT_DIR; diff them and that is
all you will see.

WHAT IT DOES
    1. writes CONFIG below to stratC_overrides.json (the file the driver
       already reads, backing up any existing one), then
    2. executes strategy_C_3dec_nodamp.py.
"""
import os, json, shutil

DRIVER   = "strategy_C_3dec_nodamp.py"
OVR_FILE = "stratC_overrides_g00.json"

CONFIG = {
    "CASE_ID":  "LS_MAXWELL_1p5",
    "OUT_DIR":  "stratC_results_LS_MAXWELL_1p5",
    "BASE_SAVE": "Part_I_MASON_v7.sav",

    # --- the two variables under test -------------------------------------
    "LARGE_STRAIN": "on",
    # 1.5% of critical over 1-40 Hz, three (weight, frequency) Maxwell pairs.
    # Same coefficients as strategy_C_3dec_maxwell.py and sweep_cases.py, so
    # this is comparable to the existing damped set. Unlike mass-proportional
    # Rayleigh it stays near-constant across the band and does not over-damp
    # the slow rocking mode.
    "MAXWELL_CMD": ("block mech damp maxwell "
                    "0.0120 1.0006 0.0093 9.0193 0.0104 40.0000"),
    "DAMP_RATIO": 0.0,        # Rayleigh OFF; Maxwell is the damping here
    "DAMP_TYPE":  "",

    # --- period identification --------------------------------------------
    # Moshfeghi et al. (2024) sec. 3.2: three accelerometers, cropped at the
    # table stop, PSD -> singular values -> peak-pick, seeded once and tracked
    # forward. Validated on the experiment's own records to +0.2% (run 1) and
    # +1.3% (run 24) against Table 5.
    "PERIOD_ID_METHOD": "experiment",
    "PERIOD_SCHEME":    "adaptive",
    # Measured from run 1 of the Maxwell set by that method. NOT the historic
    # hard-coded 0.092, which was stale: it detunes the run-1 pulse and every
    # fallback after it.
    "T1_init": 0.0948,

    # --- held identical to the NODAMP companion ---------------------------
    "xi": 0.05,
    "n_cycles": 3.0,
    "USE_TWO_IM_CYCLES": False,
    "USE_RECORD_NCYCLES": False,
    "tail_sec": 2.5,
    "inter_run_gap": 0.5,
    "delta_t": 0.005,
    "TID_AGREE_TOL": 0.25,
    "TID_MAX_JUMP": 1.50,
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
print("  period ID  : {}  (T1_init = {} s)".format(
    CONFIG["PERIOD_ID_METHOD"], CONFIG["T1_init"]))
print("  output     : {}".format(CONFIG["OUT_DIR"]))
print("=" * 70)
print("  Check the driver's own override banner below. Any key reported as")
print("  '! IGNORED' is NOT being applied, and this is then not the case you")
print("  think it is -- stop rather than finding out after 25 runs.")
print("=" * 70)

exec(compile(open(DRIVER).read(), DRIVER, "exec"), globals())