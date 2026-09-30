# -*- coding: ascii -*-
"""
pushover_oop_3dec.py  (v3) -- OOP pushover of the US-1 wall (3DEC 9),
ascending AND descending branch, ONE load pattern (mass-proportional).

BUILT ON NICOLO'S DRIVER (as translated in v2.2), WITH THE DESCENDING
BRANCH FIXED
    Kept from Nicolo / v2.2:
      * quasi-static stepping with a cycle BUDGET and heavy local damping;
        equilibrium quality is a logged quantity (conv, converged), never
        an assumption;
      * the participating mass m_part (masonry moving WITH the push) and
        his predictor  a_pred = lambda*(V - a*(m_tot - m_part))/m_part ;
      * the load: tilted gravity, i.e. a body force proportional to mass --
        the pattern the paper's Eq. 2 base shear corresponds to;
      * F_tot = cstav (base+joists) + topj_shear (top joint): the whole
        lateral reaction, comparable to Table 4 / Fig 13.

    What failed in v2.2, and why (rigid blocks, not a deformable model):
      1. The peak was declared only after the unconverged step had burnt
         its full 150k-cycle budget. In a rigid-block wall an unconverged
         step IS the mechanism running; by the time the budget was spent
         the wall was already far down the collapse path, so the
         controller started from a state that was not on the curve.
      2. Nicolo's predictor with lambda = 1.30 is OPEN-LOOP: it commands
         30% more acceleration than the still-resisting mass can hold.
         In a deformable wall that extra load is absorbed by crack growth
         and the step re-converges. In a rigid-block rocking wall there is
         no material softening left -- the descending branch is purely
         geometric (the toe migrates, the lever arm shrinks) -- and there
         is NO static equilibrium under any load above V(d)/m. The
         predicted step runs to D_STOP in one go. The same happens with
         m_part -> m_tot, which is what a rocking mechanism gives (the
         whole panel above the crack moves).
      The descending branch of a rocking wall is unstable under LOAD
      control by its nature; it can only be traced by controlling the
      RESPONSE. That is what v3 does.

    The fix (v3):
      * EARLY PEAK DETECTION: each ascending step solves in chunks; if the
        control displacement advances by more than D_RUNAWAY_MM within a
        step that has not reached RATIO, the mechanism is running -> peak,
        and the driver switches immediately, at a displacement still on
        the curve.
      * DISPLACEMENT-FOLLOWING SERVO ON THE SAME BODY FORCE: past the
        peak the independent variable becomes the control displacement
        d (mean of Top_Quarter_A/B at 2.06 m, paper Eq. 1). The lateral
        acceleration a is adjusted by feedback every CHUNK cycles,
            a <- a + KP*(d_target - d) - KV*(d - d_prev)   (PD, rate-limited,
                                                            clamped >= 0)
        until d sits at d_target and the wall is still; THEN (a, d, F_tot,
        m_part) is one point of the descending branch. Feedback is
        stable on the softening branch because the dynamic response to a
        is monotonic (more push -> forward, less push -> the precompression
        and self-weight bring it back) even though the static d(a) is
        not. The load pattern is unchanged from the ascending branch, so
        the two branches join without a pattern discontinuity (the T_B-
        driven displacement control of pushover_oop_3dec_dispctrl.py
        pushed the restrained SUPPORT and measured racking -- discarded).
      * Nicolo's predictor is still evaluated at the switch and logged
        (a_pred column) as a check on the servo's first descending point;
        it is no longer what drives the wall.

STARTING POINT: model restore of part_I_mason_LS.sav (EQDENS build,
gravity + spring load, solved to equilibrium in large strain). Record_Disp
/ cstav / ncstav / @kn / @ks return with the save; nothing is re-applied.

DAMPING WARNING: `block mech damp local 0.9` is set here. Saves written
by this run carry it -- NEVER seed a dynamic run from them.

THE x2 CHECK: the baseline printout reports ncstav: ~0.10 MPa = correct
precompression, ~0.20 MPa = doubled. Confirmed 0.103 MPa on your build.

TUNING (field, first run): KP, KV, CHUNK. Symptoms: d oscillates about
d_target -> halve KP, raise KV; the wall overshoots the target by more
than a few mm right after the peak -> raise KV and DA_SERVO_MAX (the D
term is what brakes the mechanism); the servo never settles (settle_chunks
hits SERVO_CHUNKS_MAX) -> raise D_TOL slightly or CHUNK; a saturates at 0
while d keeps growing -> past the instability point, branch ended, exit
is correct. Tested on a mock overdamped SDOF with a 29 kN peak and linear
geometric softening: traces the branch to F ~ 0 with |d - d_target| < D_TOL.

OUTPUT  <OUT_DIR>/pushover_<dir>.csv :
    phase, a_mps2, d_ctrl_mm, F_base_kN, F_top_kN, F_tot_kN, conv,
    converged, m_part_kg, a_pred_mps2, settle_chunks
  phase asc  = load-controlled ramp (Nicolo stepping)
  phase desc = servo-held points on the softening branch
  d_ctrl_mm  = mean(Top_Quarter_A/B) - baseline   [paper Eq. 1]
  F_base_kN  = cstav - baseline (base + joists);  F_top_kN = top joint
  F_tot_kN   = F_base + F_top  (capacity_law column f_tot_kn)
  conv       = |a*m_tot - |F_tot|| / (a*m_tot)   equilibrium diagnostic
RUN BOTH DIRECTIONS (PUSH_DIR = +1 / -1).
"""

import itasca as it
import os, csv, math

it.command("python-reset-state false")

# ============================ CONFIG =================================
RESTORE_SAV = "part_I_mason_LS_strong"
OUT_DIR     = "pushover_results_stronger"
PUSH_DIR    = -1          # +1 / -1 : run both, separate executions
MASONRY_RHO = 1885.0      # kg/m3, build line 24
TOP_JOINT   = "elastic"   # "elastic" or "mohr" (Fig 17: no cracking there)
TOPJ_Y      = (2.57, 2.59)

# ---- ascending (Nicolo stepping) ----
A_START     = 0.5         # m/s2
DA_ASC      = 0.5         # m/s2 per step
A_MAX       = 30.0        # m/s2 hard stop
RATIO       = 1e-5        # per-step convergence target
CYC_BUDGET  = 150000      # cycles per step before "unconverged"
CHUNK_ASC   = 5000        # cycles between runaway checks inside a step
D_RUNAWAY_MM = 3.0        # control-point advance within ONE unconverged
                          # step that means "the mechanism is running"
# ---- descending (servo) ----
DD_DESC     = 1.0         # mm per descending point
D_STOP_MM   = 80.0        # end of the traced branch
V_EXIT_KN   = 0.5         # or the resistance is gone
DESC_POINTS_MAX = 120
CHUNK       = 500         # cycles between servo updates
KP          = 0.30        # m/s2 per mm of displacement error   (P term)
KV          = 3.0         # m/s2 per (mm per chunk) of control velocity (D term:
                          # brakes the mechanism before the error has grown)
DA_SERVO_MAX = 0.50       # m/s2, max change of a per servo update
D_TOL       = 0.05        # mm, |d - d_target| for "on target"
N_STILL     = 4           # consecutive on-target, non-moving chunks = settled
SERVO_CHUNKS_MAX = 400    # per point (400 x 500 = 200k cycles)
LAMBDA_PRED = 1.30        # Nicolo's lambda -- for the LOGGED predictor only

os.makedirs(OUT_DIR, exist_ok=True)

# ============================ RESTORE ================================
if not (os.path.isfile(RESTORE_SAV) or os.path.isfile(RESTORE_SAV + ".sav")):
    raise RuntimeError("missing save file: {}(.sav)".format(RESTORE_SAV))
it.command("model restore '{}'".format(RESTORE_SAV))
it.command("python-reset-state false")
it.command("model large-strain on")       # the save is LS; guard only
it.command("block mech damp local 0.9")   # quasi-static (see warning)

# ---- top joint -------------------------------------------------------
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
it.command("[global push_sign = {:d}]".format(int(PUSH_DIR)))

# top-joint reaction (same pattern as the build's cstav; local accumulator,
# name assigned once -- avoids the nested-call recursion)
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

# masonry mass, participating mass (Nicolo's loop on blocks), peak speed
it.command("""
fish define wall_mass_part
    global m_tot = 0.0
    global m_part = 0.0
    global wall_maxvel = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry') then
            local mb = {rho} * block.vol(b)
            local vz = block.vel.z(b)
            m_tot = m_tot + mb
            if (vz * push_sign) > 0.0 then
                m_part = m_part + mb
            endif
            if math.abs(vz) > wall_maxvel then
                wall_maxvel = math.abs(vz)
            endif
        endif
    endloop
end
""".format(rho=MASONRY_RHO))

it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CYC_BUDGET))

# ============================ BASELINE ===============================
it.command("@Record_Disp")
it.command("@wall_mass_part")
d0  = 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1000.0
F0b = it.fish.call_function("cstav")
F0t = it.fish.call_function("topj_shear")
m_tot = it.fish.get("m_tot")
sig_n = it.fish.get("ncstav")
print("settled baseline: d = {:.4f} mm, F_base = {:.4f} kN, F_top = {:.4f} kN, "
      "wall mass = {:.1f} kg (hand value 1318)".format(d0, F0b, F0t, m_tot))
print("avg base normal stress = {:.3f} MPa (~0.10 correct; ~0.20 doubled)"
      .format(abs(sig_n) / 1e6))

# ============================ HELPERS ================================
def get_cycles():
    try:
        it.command("[global _cyc_now = mech.step]")
        return int(it.fish.get("_cyc_now"))
    except Exception:
        return int(it.cycle())

def set_a(a_mps2):
    it.command("model gravity 0 -9.81 {:.6f}".format(PUSH_DIR * a_mps2))

def ctrl_disp():
    """control displacement, mm, signed, relative to the settled baseline"""
    it.command("@Record_Disp")
    return 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
                  it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0

def read_state():
    d = ctrl_disp()
    it.command("@wall_mass_part")
    Fb = it.fish.call_function("cstav") - F0b
    Ft = it.fish.call_function("topj_shear") - F0t
    return d, Fb, Ft, it.fish.get("m_part"), it.fish.get("wall_maxvel")

def conv_of(a, Ftot):
    tgt = a * m_tot / 1000.0
    return abs(tgt - abs(Ftot)) / tgt if tgt > 0 else 0.0

def solve_step(a_mps2):
    """Nicolo's budgeted quasi-static step with EARLY runaway detection.
    Returns (converged, ran_away)."""
    set_a(a_mps2)
    c0 = get_cycles()
    d_start = ctrl_disp()
    while get_cycles() - c0 < CYC_BUDGET:
        c1 = get_cycles()
        it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CHUNK_ASC))
        if get_cycles() - c1 < CHUNK_ASC:      # ratio reached inside the chunk
            return True, False
        if (ctrl_disp() - d_start) * PUSH_DIR > D_RUNAWAY_MM:
            return False, True                  # mechanism running -> peak
    return False, False                         # budget spent, no runaway

def nicolo_predictor(a, Ftot, m_part):
    V_N = abs(Ftot) * 1e3
    m_p = max(m_part, 0.05 * m_tot)
    return max(LAMBDA_PRED * (V_N - a * (m_tot - m_p)) / m_p, 0.0)

def servo_to(d_target, a):
    """Hold the wall at d_target by feedback on the body-force multiplier.
    Returns (settled, a, chunks_used)."""
    still = 0
    d_prev = ctrl_disp()
    for j in range(SERVO_CHUNKS_MAX):
        set_a(a)
        it.command("model cycle {:d}".format(CHUNK))
        d = ctrl_disp()
        err = (d_target - d) * PUSH_DIR           # >0: needs more push
        vel = (d - d_prev) * PUSH_DIR             # mm per chunk, >0: advancing
        da = max(-DA_SERVO_MAX, min(DA_SERVO_MAX, KP * err - KV * vel))
        a = max(0.0, min(A_MAX, a + da))
        moving = abs(d - d_prev) > 0.5 * D_TOL
        d_prev = d
        still = still + 1 if (abs(err) < D_TOL and not moving) else 0
        if still >= N_STILL:
            return True, a, j + 1
    return False, a, SERVO_CHUNKS_MAX

rows = []
csv_path = os.path.join(OUT_DIR, "pushover_{}.csv".format("pos" if PUSH_DIR > 0 else "neg"))

# ============================ ASCENDING ==============================
print("\nascending: dir {:+d}, da {:.2f} m/s2, budget {} cycles/step, "
      "runaway check every {} cycles (> {} mm)".format(
          int(PUSH_DIR), DA_ASC, CYC_BUDGET, CHUNK_ASC, D_RUNAWAY_MM))
a = 0.0
a_last_ok = 0.0
peak = False
while a < A_MAX:
    a = A_START if a == 0.0 else a + DA_ASC
    ok, ran = solve_step(a)
    d, Fb, Ft, m_part, vmax = read_state()
    Ftot = Fb + Ft
    a_pred = nicolo_predictor(a, Ftot, m_part)
    rows.append(("asc", round(PUSH_DIR * a, 4), round(d, 4), round(Fb, 4), round(Ft, 4),
                 round(Ftot, 4), round(conv_of(a, Ftot), 5), 1 if ok else 0,
                 round(m_part, 1), round(a_pred, 4), ""))
    print("  asc a {:6.2f}  d {:8.3f} mm  Fb {:7.3f}  Ft {:7.3f}  Ftot {:8.3f} kN  "
          "conv {:.4f}  m_part {:6.0f}  {}".format(
              PUSH_DIR * a, d, Fb, Ft, Ftot, conv_of(a, Ftot), m_part,
              "ok" if ok else ("RUNAWAY" if ran else "UNCONVERGED")))
    if ok:
        a_last_ok = a
        if abs(d) > D_STOP_MM:
            print("  stop: |d| > {} mm on the ascending branch".format(D_STOP_MM)); break
        continue
    peak = True
    print("  ** PEAK: a = {:.2f} m/s2 ({}). Last converged a = {:.2f}. "
          "Nicolo predictor a_pred = {:.2f} (logged, not applied). "
          "Switching to the displacement servo at d = {:.2f} mm.".format(
              a, "runaway" if ran else "budget", a_last_ok, a_pred, d))
    break

# ============================ DESCENDING (servo) =====================
if peak:
    d_sw = ctrl_disp()
    a = a_last_ok                       # arrest: below the peak load
    d_target = d_sw                     # first: hold where we are
    for k in range(DESC_POINTS_MAX):
        settled, a, nch = servo_to(d_target, a)
        d, Fb, Ft, m_part, vmax = read_state()
        Ftot = Fb + Ft
        a_pred = nicolo_predictor(a, Ftot, m_part)
        rows.append(("desc", round(PUSH_DIR * a, 4), round(d, 4), round(Fb, 4), round(Ft, 4),
                     round(Ftot, 4), round(conv_of(a, Ftot), 5), 1 if settled else 0,
                     round(m_part, 1), round(a_pred, 4), nch))
        print("  desc d_t {:7.2f}  d {:8.3f} mm  a {:6.3f}  Ftot {:8.3f} kN  "
              "m_part {:6.0f}  a_pred {:5.2f}  {} ({} chunks)".format(
                  d_target, d, PUSH_DIR * a, Ftot, m_part, a_pred,
                  "settled" if settled else "NOT settled", nch))
        if abs(Ftot) < V_EXIT_KN and settled:
            print("  descending exit: |F_tot| < {} kN".format(V_EXIT_KN)); break
        if abs(d) > D_STOP_MM:
            print("  descending exit: |d| > {} mm".format(D_STOP_MM)); break
        if a <= 0.0 and (d - d_target) * PUSH_DIR > 5.0 * D_TOL:
            print("  descending exit: a = 0 and the wall still advances -- "
                  "past the instability point (resistance exhausted)"); break
        d_target = d_target + PUSH_DIR * DD_DESC
    it.command("model save '{}'".format(
        os.path.join(OUT_DIR, "pushover_{}_end".format(
            "pos" if PUSH_DIR > 0 else "neg")).replace("\\", "/")))

# ============================ OUTPUT =================================
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["phase", "a_mps2", "d_ctrl_mm", "F_base_kN", "F_top_kN", "F_tot_kN",
                "conv", "converged", "m_part_kg", "a_pred_mps2", "settle_chunks"])
    w.writerows(rows)
print("-> {}".format(csv_path))
print("capacity_law reads d_ctrl_mm / F_tot_kN directly. Peak F_tot vs Table 4 "
      "(+29.5 / -29.5 kN); K_sec(d) = F/d along both branches.")