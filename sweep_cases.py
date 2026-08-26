# -*- coding: ascii -*-
"""
Sensitivity case matrix for the Strategy C adaptive equivalent-pulse IDA.
=======================================================================
Design: one-at-a-time (OAT) screening around the corrected baseline, then a
2-level full factorial on whichever parameters the screening ranks highest.

Damage is CUMULATIVE across the protocol, so a case cannot be a subset of
runs -- every case is a full 1..25 sequence. That is what sets the cost, and
why the screening stage is OAT rather than factorial.

Baseline (case 00) is the corrected material of mat_case.dat plus the
two-IM cycle count and the guarded adaptive period.

MATERIAL keys are written into mat_case.dat and require a Part-I rebuild.
DRIVER keys go into stratC_overrides.json and do not.
"""

# v8 material, kept verbatim so the loading change can be tested against the
# v8 result without the material correction confounding it. tension/cohesion
# are the pre-correction values and Gc is the mis-transcribed expression.
V8_MATERIAL = dict(tension=0.200e6, cohesion=0.400e6, Gt=32.0,
                   Gc_expr="mis_transcribed")

BASELINE = dict(
    # --- material (Part-I rebuild required) ---
    tension=0.344e6, cohesion=0.750e6, Gt=55.0, Gc_expr="jafari_eq3",
    friction=36.9, fc_comp=20e6, large_strain="off",
    # --- driver ---
    USE_TWO_IM_CYCLES=True, n_cycles=3.0, PERIOD_SCHEME="adaptive",
    xi=0.05, DAMP_RATIO=0.0, DAMP_TYPE="", MAXWELL_CMD="",
)

MAXWELL = ("block mech damp maxwell 0.0120 1.0006 0.0093 9.0193 "
           "0.0104 40.0000")

# --- controlled loading test ---------------------------------------------
# Case 00 changes BOTH the material and the loading relative to NODAMP_v8, and
# the two push opposite ways: the corrected material (tension 0.200 -> 0.344,
# cohesion 0.400 -> 0.750 MPa) makes the wall stronger and will REDUCE
# displacement, while the two-IM cycle count raises the pulse acceleration
# 2.4-3.4x and should increase it. L00 holds the v8 material so the loading
# change is the single variable against v8, which is the clean test of the
# acceleration-threshold hypothesis.
#
# Read it on model EDP / Sd_target: v8 sat at 1.0-1.3 across the sequence
# (near-elastic). Above ~2 means the rocking mechanism has engaged.
CONTROL = [
    ("L00", "control", "two-IM cycles, v8 material", dict(V8_MATERIAL)),
    ("G00", "geometry", "large-strain ON, v8 material",
     dict(list(V8_MATERIAL.items()) + [("large_strain", "on")])),
    ("L01", "control", "v8 route + v8 material (v8 replica)",
     dict(list(V8_MATERIAL.items()) + [("USE_TWO_IM_CYCLES", False),
                                       ("n_cycles", 3.0)])),
]

# --- OAT screening: (case_id, group, what changed, {overrides}) ----------
SCREEN = [
    ("S01", "loading",  "fixed 3.0 cycles (v8 route)", dict(USE_TWO_IM_CYCLES=False, n_cycles=3.0)),
    ("S02", "loading",  "fixed 1.5 cycles",            dict(USE_TWO_IM_CYCLES=False, n_cycles=1.5)),
    ("S03", "loading",  "T_eq fixed period",           dict(PERIOD_SCHEME="fixed")),
    ("S04", "loading",  "xi_cal = 0.02",               dict(xi=0.02)),
    ("S05", "loading",  "xi_cal = 0.10",               dict(xi=0.10)),
    ("S06", "damping",  "Rayleigh 1%",                 dict(DAMP_RATIO=0.01)),
    ("S07", "damping",  "Rayleigh 3%",                 dict(DAMP_RATIO=0.03)),
    ("S08", "damping",  "Maxwell 1.5% (1-40 Hz)",      dict(MAXWELL_CMD=MAXWELL)),
    ("S09", "material", "tension 0.24 MPa (-30%)",     dict(tension=0.24e6)),
    ("S10", "material", "tension 0.44 MPa (+28%)",     dict(tension=0.44e6)),
    ("S11", "material", "cohesion 0.544 (Fig10b fit)", dict(cohesion=0.544e6)),
    ("S12", "material", "cohesion 1.03 (ratio 3.0)",   dict(cohesion=1.03e6)),
    ("S13", "material", "G_I 27.5 (x0.5)",             dict(Gt=27.5)),
    ("S14", "material", "G_I 110 (x2)",                dict(Gt=110.0)),
    ("S15", "material", "Gc 26.1 N/mm (Fig6a)",        dict(Gc_expr="jafari_fig6a")),
]

# --- factorial stage: filled after screening ----------------------------
# Pick the 3 highest-ranked parameters from the screening tornado and list
# their (low, high) levels here; expand() builds the 2^3 = 8 corner cases,
# of which the baseline corner is already run.
FACTORIAL_FACTORS = []      # e.g. [("tension",0.24e6,0.44e6), ...]


def expand_factorial(factors):
    import itertools
    if not factors:
        return []
    out = []
    for k, combo in enumerate(itertools.product(*[(0, 1)] * len(factors))):
        ov, tag = {}, []
        for (name, lo, hi), lvl in zip(factors, combo):
            ov[name] = hi if lvl else lo
            tag.append("{}{}".format(name[:4], "H" if lvl else "L"))
        out.append(("F{:02d}".format(k + 1), "factorial", "+".join(tag), ov))
    return out


def all_cases():
    cases = [("00", "baseline", "corrected baseline", {})]
    cases += CONTROL
    cases += SCREEN
    cases += expand_factorial(FACTORIAL_FACTORS)
    return cases


MATERIAL_KEYS = ("tension", "cohesion", "Gt", "Gc_expr", "friction",
                 "fc_comp", "large_strain")

if __name__ == "__main__":
    cs = all_cases()
    print("{} cases ({} need a Part-I rebuild)".format(
        len(cs), sum(1 for _, _, _, o in cs if any(k in MATERIAL_KEYS for k in o))))
    for cid, grp, desc, ov in cs:
        print("  {:<4} {:<9} {:<30} {}".format(cid, grp, desc, ov))