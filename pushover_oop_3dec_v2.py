# -*- coding: ascii -*-
"""
pushover_oop_3dec_v2.py -- OOP pushover of the US-1 wall (3DEC 9),
ascending AND descending branch.

WHAT CHANGED vs pushover_oop_3dec.py (v1), AND WHERE IT CAME FROM
    v1 was load-controlled: ramp lateral gravity, and the step where the
    static solve stops converging IS the peak -- but the curve ends there.
    K_sec beyond the peak was unreachable, which the capacity law needs
    for d up to ~30 mm if the model peaks earlier.

    This version adopts two ideas from Nicolo's pushover driver:
    1. QUASI-STATIC STEPPING WITH A BUDGET. Heavy local damping (0.9) and
       `model solve ratio ... cycles ...`: a step that cannot reach the
       ratio inside the budget does not kill the run -- it is LOGGED as
       unconverged (conv column) and the controller reacts. Equilibrium
       quality is a recorded quantity, not an assumption.
    2. THE PARTICIPATING-MASS DESCENDING CONTROLLER. Past the peak, load
       control cannot trace softening. Nicolo's controller measures the
       base shear V and the mass still moving WITH the push (m_part),
       then commands
           a_next = lambda * (V - a*(m_tot - m_part)) / m_part
       i.e. push only as hard as the still-resisting mass can carry,
       amplified by lambda (1.30) so the mechanism keeps failing. Exit
       when V ~ 0 or the control displacement passes D_STOP_MM.

    TRANSLATED, NOT COPIED: Nicolo's model is deformable (block.zone.*,
    gridpoint BCs, z-up, push in y, density 1835, half-model factor 0.5).
    Ours is rigid-block, y-up, push in z, Masonry density 1885, full
    model. Zone loops became block loops; his @shearforce is our @cstav
    (kN, from the build's own Reaction-group FISH); his manual restart at
    a fixed acceleration became an automatic switch on lost convergence.

STARTING POINT: model restore of part_I_mason_LS.sav -- built with the
equivalent-density file, gravity + spring load applied, solved to
equilibrium in large strain. Record_Disp / cstav / @kn / @ks all return
with the save; nothing is re-applied.

DAMPING WARNING: this driver sets `block mech damp local 0.9` for the
quasi-static push. Saves written by this run carry that damping --
NEVER seed a dynamic run from them.

THE x2 CHECK: the baseline printout reports ncstav from the restored
state. ~0.10 MPa -> precompression correct; ~0.20 MPa -> doubled (in the
dynamic runs too) -- resolve that before trusting this curve.

UNVERIFIED FISH INTRINSICS (rigid-block side): block.list,
block.isgroup, block.vol, block.vel.z. Your files verify the contact-
side family and block.gp.*; these are the standard names by the same
convention. If 3DEC rejects the @wall_mass_part define, send the error
text -- it will be a rename, not a logic change.

OUTPUT
    <OUT_DIR>/pushover_<dir>.csv :
        phase, a_mps2, d_ctrl_mm, F_base_kN, conv, m_part_kg
    phase asc = load-controlled ramp; phase desc = participating-mass
    controller. d_ctrl = mean(Top_Quarter_A/B) - baseline (Ch3/Ch4
    positions, paper Eq. 1). F_base = cstav - baseline (paper Eq. 2
    counterpart). RUN BOTH DIRECTIONS (PUSH_DIR = +1 / -1): Table 4
    reports them separately and the joists make the wall asymmetric.
"""

import itasca as it
import os, csv, math

it.command("python-reset-state false")

# ============================ CONFIG =================================
RESTORE_SAV = "part_I_mason_LS"   # .sav from the EQDENS build, solved in LS
OUT_DIR     = "pushover_results"
PUSH_DIR    = +1        # +1 / -1 : run both, separate executions
A_START     = 0.5       # m/s2, first ascending step
DA_ASC      = 0.5       # m/s2 per ascending step (expected peak ~22 m/s2
                        # if Table 4's 29.5 kN / 1318 kg transfers)
A_MAX       = 30.0      # m/s2 hard stop for the ascending ramp
RATIO       = 1e-5      # per-step convergence target
CYC_BUDGET  = 150000    # cycles allowed per step before it is declared
                        # unconverged (Nicolo's time cap, in cycles)
CONV_SWITCH = 0.10      # |a*m - V|/(a*m) above this after a full budget
                        # -> peak reached -> descending controller starts
LAMBDA_DESC = 1.30      # Nicolo's amplification on the descending branch
V_EXIT_KN   = 0.5       # descending exit: base shear essentially gone
D_STOP_MM   = 80.0      # or the control displacement is far past 30 mm
DESC_STEPS_MAX = 200    # safety cap on descending iterations
TOP_JOINT   = "elastic" # "elastic" or "mohr" (Fig 17: no cracking there)
TOPJ_Y      = (2.57, 2.59)
MASONRY_RHO = 1885.0    # kg/m3, from the build (line 24)

os.makedirs(OUT_DIR, exist_ok=True)

def fget(name):
    try:
        return it.fish.get(name)
    except Exception:
        return it.fish.get(name.lower())

# ============================ RESTORE ================================
if not (os.path.isfile(RESTORE_SAV) or os.path.isfile(RESTORE_SAV + ".sav")):
    raise RuntimeError("missing save file: {}(.sav)".format(RESTORE_SAV))
it.command("model restore '{}'".format(RESTORE_SAV))
it.command("python-reset-state false")
it.command("model large-strain on")     # guard; the save is already LS

# quasi-static: heavy local damping, Nicolo's setting. See DAMPING WARNING.
it.command("block mech damp local 0.9")

# ---- top joint (unchanged from v1) ----------------------------------
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

# ---- masonry mass + participating mass (Nicolo's zone loop, on blocks)
it.command("""
fish define wall_mass_part
    global m_tot = 0.0
    global m_part = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry')
            local mb = {rho} * block.vol(b)
            m_tot = m_tot + mb
            if math.sgn(block.vel.z(b)) == math.sgn(push_sign)
                m_part = m_part + mb
            endif
        endif
    endloop
end
""".format(rho=MASONRY_RHO))
it.command("[global push_sign = {:d}]".format(int(PUSH_DIR)))

# re-settle after the top-joint reassignment
it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CYC_BUDGET))

# ============================ BASELINE ===============================
it.command("@Record_Disp")
it.command("@cstav")
it.command("@wall_mass_part")
d0 = 0.5 * (fget("Top_Quarter_A_Disp") + fget("Top_Quarter_B_Disp")) * 1000.0
F0 = fget("cstav")
m_tot = fget("m_tot")
sig_n = fget("ncstav")
print("settled baseline: d = {:.4f} mm, F = {:.4f} kN, wall mass = {:.1f} kg "
      "(hand value 1318)".format(d0, F0, m_tot))
print("avg base normal stress = {:.3f} MPa  "
      "(~0.10 = precompression correct; ~0.20 = the x2 doubles it)"
      .format(abs(sig_n) / 1e6))
if abs(m_tot - 1318.0) / 1318.0 > 0.05:
    print("** WARNING: FISH wall mass differs from the hand value by >5% -- "
          "check block.vol/group before trusting the descending controller.")

# ============================ HELPERS ================================
def read_state():
    it.command("@Record_Disp")
    it.command("@cstav")
    it.command("@wall_mass_part")
    d = 0.5 * (fget("Top_Quarter_A_Disp") +
               fget("Top_Quarter_B_Disp")) * 1000.0 - d0
    F = fget("cstav") - F0            # kN, signed
    return d, F, fget("m_part")

def solve_step(a_mps2):
    it.command("model gravity 0 -9.81 {:.6f}".format(a_mps2))
    it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CYC_BUDGET))

rows = []
csv_path = os.path.join(OUT_DIR,
    "pushover_{}.csv".format("pos" if PUSH_DIR > 0 else "neg"))

# ============================ ASCENDING ==============================
print("\nascending: dir {:+d}, da {:.2f} m/s2, budget {} cycles/step"
      .format(int(PUSH_DIR), DA_ASC, CYC_BUDGET))
a = 0.0
a_peak = None
while a < A_MAX:
    a = A_START if a == 0.0 else a + DA_ASC
    solve_step(PUSH_DIR * a)
    d, F, m_part = read_state()
    # equilibrium quality: commanded inertial load vs measured base shear
    target_kN = a * m_tot / 1000.0
    conv = abs(target_kN - abs(F)) / target_kN if target_kN > 0 else 0.0
    rows.append(("asc", round(PUSH_DIR * a, 4), round(d, 4), round(F, 4),
                 round(conv, 5), round(m_part, 1)))
    print("  asc a {:6.2f}  d {:8.3f} mm  F {:8.3f} kN  conv {:.4f}"
          .format(PUSH_DIR * a, d, F, conv))
    if conv > CONV_SWITCH:
        a_peak = a
        print("  ** equilibrium lost at a = {:.2f} m/s2 (conv {:.3f} > {}) "
              "-- PEAK. Switching to the descending controller."
              .format(a, conv, CONV_SWITCH))
        break
    if abs(d) > D_STOP_MM:
        print("  stop: |d| > {} mm on the ascending branch".format(D_STOP_MM))
        break

# ============================ DESCENDING =============================
# Nicolo's participating-mass controller, on rigid blocks. The commanded
# acceleration follows what the still-resisting mass can carry:
#     a_next = LAMBDA * (V - a*(m_tot - m_part)) / m_part
# (all in N and m/s2 internally; cstav is kN, hence the 1e3.)
if a_peak is not None:
    a_i = a_peak
    for k in range(DESC_STEPS_MAX):
        solve_step(PUSH_DIR * a_i)
        d, F, m_part = read_state()
        rows.append(("desc", round(PUSH_DIR * a_i, 4), round(d, 4),
                     round(F, 4), None, round(m_part, 1)))
        print("  desc a {:6.2f}  d {:8.3f} mm  F {:8.3f} kN  m_part {:7.1f} kg"
              .format(PUSH_DIR * a_i, d, F, m_part))
        if abs(F) < V_EXIT_KN:
            print("  descending exit: |F| < {} kN".format(V_EXIT_KN))
            break
        if abs(d) > D_STOP_MM:
            print("  descending exit: |d| > {} mm".format(D_STOP_MM))
            break
        if m_part < 0.02 * m_tot:
            print("  descending exit: participating mass ~ 0")
            break
        V_N = abs(F) * 1e3
        a_new = (V_N - a_i * (m_tot - m_part)) / m_part
        a_i = max(LAMBDA_DESC * a_new, 0.05)
    it.command("model save '{}'".format(
        os.path.join(OUT_DIR, "pushover_{}_end".format(
            "pos" if PUSH_DIR > 0 else "neg")).replace("\\", "/")))

# ============================ OUTPUT =================================
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["phase", "a_mps2", "d_ctrl_mm", "F_base_kN", "conv", "m_part_kg"])
    w.writerows(rows)
print("-> {}".format(csv_path))
print("capacity_law reads this file directly once you rename/derive columns "
      "d_ctrl_mm and F_base_kN -- they are already named so; no edit needed.")
print("next: peak F vs Table 4 (+29.5 / -29.5 kN); K_sec(d) = F/d including "
      "the descending branch; T_sec = 2*pi*sqrt(m_eff/K_sec).")
