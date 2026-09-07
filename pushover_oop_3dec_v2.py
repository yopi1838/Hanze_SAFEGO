# -*- coding: ascii -*-
"""
pushover_oop_3dec_dispctrl.py -- DISPLACEMENT-CONTROLLED OOP pushover of
the US-1 wall (3DEC 9). Traces the FULL capacity curve, rising branch AND
post-peak softening, in one consistent loading condition.

WHY DISPLACEMENT CONTROL (and why the load-controlled v2.2 could not do it)
    v2.2 ramped a lateral body force (tilted gravity) and read the
    reaction. That is load control: past the limit point there is NO
    static equilibrium at any load >= the peak, so the solve diverges and
    the softening branch is unreachable. v2.2's ascending branch is still
    the validation anchor -- it peaked at F_tot = 29.0 kN, matching
    Table 4's 29.5 kN within 2%.
    Here displacement is the independent variable, so F is free to
    DECREASE as d grows -- that is what makes the descending branch
    traceable.

HOW
    The top beam T_B is driven at a small constant velocity in the push
    direction (an actuator platen). Gravity is VERTICAL ONLY -- no lateral
    body force any more; the lateral demand now comes entirely from the
    imposed top motion. Heavy local damping (0.9) + short dynamic chunks
    keep it quasi-static; a kinetic-energy monitor reports whether that
    held (KE_frac column; large -> slow VPUSH down).

    d_ctrl = mean(Top_Quarter_A/B) at 2.06 m -- the paper's Eq.(1) EDP,
             the same point the dynamic scheme matches Sd against.
    F_base = base + joist reaction (cstav), sign-flipped so it is POSITIVE
             with the drift -- this is the base shear, Fig 13's y-axis, and
             the column capacity_law consumes.
    F_appl = the T_B/Masonry joint force (topj) -- the actuator force. In
             equilibrium |F_appl| ~ |F_base|; the two tracking each other
             is the quasi-static check. They are NOT summed (that was the
             load-control quantity; here summing would double-count).

  ** LOAD-PATTERN CAVEAT -- state this in any write-up. **
    Driving T_B applies the lateral load at the TOP (2.68 m). v2.2's body
    force was mass-proportional (distributed). A top load has a longer
    lever about the base hinge, so this pushover may PEAK BELOW v2.2's
    29 kN -- that is the pattern difference, not an error. The check is
    built in: the run prints its own peak F_base next to the 29.0 kN
    load-control anchor. If they are close, the pattern effect is second
    order and the descending branch is trustworthy; if far apart, the
    capacity is pattern-sensitive and BOTH numbers must be reported. A
    mass-proportional displacement control (servo on a body-force
    multiplier) removes this caveat but needs an arc-length solver 3DEC
    does not provide out of the box -- flagged for Nicolo, not attempted
    here.

STARTING POINT / DAMPING / x2 CHECK / UNVERIFIED INTRINSICS: identical to
v2.2 -- see that file's header. Saves carry damp 0.9; never seed a
dynamic run from them.

OUTPUT  <OUT_DIR>/pushover_disp_<dir>.csv :
    d_ctrl_mm, F_base_kN, F_appl_kN, KE_frac, maxvel_mms
    F_base_kN is the capacity-law input (base shear, +ve with drift).
"""

import itasca as it
import os, csv, math

it.command("python-reset-state false")

# ============================ CONFIG =================================
RESTORE_SAV = "part_I_mason_LS"
OUT_DIR     = "pushover_results"
PUSH_DIR    = +1          # +1 / -1 : run both, separate executions
VPUSH       = 1.0e-4      # m/s imposed top velocity (quasi-static). If
                          # KE_frac stays high, halve it and rerun.
DT_STEP     = 0.5         # s of model time per recorded point
                          # -> d increment ~ VPUSH*DT_STEP = 0.05 mm/point
D_STOP_MM   = 40.0        # push comfortably past the 29.5 mm target
NPTS_MAX    = 1200        # safety cap
TOP_JOINT   = "elastic"
TOPJ_Y      = (2.57, 2.59)
MASONRY_RHO = 1885.0
KE_FRAC_WARN = 0.05       # KE / |push work| above this -> not quasi-static
F_ANCHOR_KN = 29.0        # v2.2 load-control peak, for the pattern check

os.makedirs(OUT_DIR, exist_ok=True)

# ============================ RESTORE ================================
if not (os.path.isfile(RESTORE_SAV) or os.path.isfile(RESTORE_SAV + ".sav")):
    raise RuntimeError("missing save file: {}(.sav)".format(RESTORE_SAV))
it.command("model restore '{}'".format(RESTORE_SAV))
it.command("python-reset-state false")
it.command("model large-strain on")
it.command("model dynamic active on")
it.command("block mech damp local 0.9")

# gravity stays VERTICAL -- the lateral demand is the imposed top motion
it.command("model gravity 0 -9.81 0")

# ---- top joint (same as v2.2) ---------------------------------------
if TOP_JOINT == "elastic":
    it.command("block contact jmodel assign elastic "
               "range group-intersection 'T_B' 'Masonry'")
    it.command("block contact property stiffness-normal @kn stiffness-shear @ks "
               "range group-intersection 'T_B' 'Masonry'")
else:
    it.command("block contact jmodel assign mohr "
               "range group-intersection 'T_B' 'Masonry'")
    it.command("block contact property stiffness-normal @kn stiffness-shear @ks "
               "tension 0 cohesion 0 friction 35 "
               "range group-intersection 'T_B' 'Masonry'")
it.command("block contact group 'TopJ' range pos-y {} {}".format(*TOPJ_Y))

# top-joint force (actuator) -- accumulate into a LOCAL, assign name once
it.command("""
fish define topj_shear
    local _tjs = 0.0
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'TopJ') then
            loop foreach local sc block.contact.subcontactlist(cx)
                _tjs = _tjs + block.subcontact.force.shear.z(sc)
            endloop
        endif
    endloop
    topj_shear = _tjs / 1000.0
end
""")

# quasi-static monitor: peak masonry block speed + kinetic energy
it.command("""
fish define wall_kinematics
    local _mv = 0.0
    local _ke = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry') then
            local vz = block.vel.z(b)
            local sp = math.abs(vz)
            if sp > _mv then
                _mv = sp
            endif
            _ke = _ke + 0.5 * {rho} * block.vol(b) * vz * vz
        endif
    endloop
    global wall_maxvel = _mv
    global wall_ke = _ke
end
""".format(rho=MASONRY_RHO))

# ---- drive T_B at constant velocity in the push direction -----------
# build already fixed rotation-y/z on T_B and set vel-x = 0; we override
# vel-z with the push. vel-y stays free so the beam can settle vertically.
it.command("block apply velocity-z {:.8f} range group 'T_B'"
           .format(VPUSH * PUSH_DIR))

# ============================ BASELINE ===============================
it.command("@Record_Disp")
d0 = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
            it.fish.get("Top_Quarter_B_Disp")) * 1000.0
F0b = it.fish.call_function("cstav")       # base+joist reaction, kN
F0t = it.fish.call_function("topj_shear")  # top-joint force, kN
sig_n = it.fish.get("ncstav")
print("baseline: d = {:.4f} mm, cstav = {:.4f} kN, topj = {:.4f} kN".format(
    d0, F0b, F0t))
print("avg base normal stress = {:.3f} MPa  "
      "(~0.10 correct; ~0.20 = x2 doubles it)".format(abs(sig_n) / 1e6))
print("driving T_B at {:.2e} m/s ({:+d}); vertical gravity only."
      .format(VPUSH, int(PUSH_DIR)))

# ============================ PUSH ===================================
rows = []
csv_path = os.path.join(OUT_DIR,
    "pushover_disp_{}.csv".format("pos" if PUSH_DIR > 0 else "neg"))
peak_F = 0.0
print("\ndisplacement-controlled push:")
for k in range(NPTS_MAX):
    it.command("model solve dynamic time {:.6f}".format(DT_STEP))
    it.command("@Record_Disp")
    it.command("@wall_kinematics")
    d = 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
               it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0
    # base shear as a positive resistance with the drift (Fig 13 sign):
    # cstav is the reaction (opposes motion), so resistance = -reaction,
    # then folded to the push direction.
    F_base = -(it.fish.call_function("cstav") - F0b) * PUSH_DIR
    F_appl = (it.fish.call_function("topj_shear") - F0t) * PUSH_DIR
    maxv = it.fish.get("wall_maxvel")
    ke = it.fish.get("wall_ke")
    push_work = abs(F_appl) * 1e3 * abs(d) / 1000.0 + 1e-9   # J, rough
    ke_frac = ke / push_work
    rows.append((round(d, 4), round(F_base, 4), round(F_appl, 4),
                 round(ke_frac, 5), round(maxv * 1000.0, 4)))
    if abs(F_base) > abs(peak_F):
        peak_F = F_base
    if k % 20 == 0 or abs(d) > D_STOP_MM:
        print("  d {:8.3f} mm  F_base {:8.3f} kN  F_appl {:8.3f} kN  "
              "KE_frac {:.4f}  maxv {:.3f} mm/s".format(
                  d, F_base, F_appl, ke_frac, maxv * 1000.0))
    if ke_frac > KE_FRAC_WARN and abs(d) > 1.0:
        print("  ** KE_frac {:.3f} > {} -- NOT quasi-static at d = {:.2f} mm. "
              "Halve VPUSH and rerun.".format(ke_frac, KE_FRAC_WARN, d))
    if abs(d) > D_STOP_MM:
        print("  stop: |d| > {} mm".format(D_STOP_MM))
        break

with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["d_ctrl_mm", "F_base_kN", "F_appl_kN", "KE_frac", "maxvel_mms"])
    w.writerows(rows)
it.command("model save '{}'".format(
    os.path.join(OUT_DIR, "pushover_disp_{}_end".format(
        "pos" if PUSH_DIR > 0 else "neg")).replace("\\", "/")))

print("-> {}".format(csv_path))
print("\nPEAK base shear (disp control) = {:.2f} kN".format(abs(peak_F)))
print("load-control anchor (v2.2)     = {:.2f} kN".format(F_ANCHOR_KN))
print("Table 4                        = 29.5 kN")
dev = abs(abs(peak_F) - F_ANCHOR_KN) / F_ANCHOR_KN * 100.0
if dev < 15.0:
    print("-> within {:.0f}% of the load-control peak: load pattern is a "
          "second-order effect, descending branch trustworthy.".format(dev))
else:
    print("-> {:.0f}% below the load-control peak: capacity IS pattern-"
          "sensitive. Report BOTH; the top-load number is a lower bound."
          .format(dev))
print("next: F_base_kN column feeds capacity_law directly (base shear, "
      "+ve with drift). K_sec(29.5 mm) now comes from the model, not the "
      "flag-extrapolated envelope.")