# -*- coding: ascii -*-
"""
pushover_oop_3dec_v5.py -- OOP pushover of the US-1 wall (3DEC 9), v3 servo driver unchanged,
mass-proportional push applied to the MASONRY ONLY (v5), with diagnostics.

v5 change vs v4:
  The lateral load is no longer a tilted gravity. Gravity stays vertical and the push is an
  applied force a*m_i on every masonry block (block.force.app, z component). Reason: with
  tilted gravity the joists also receive a*m along their own axis; on the near-frictionless
  roller nothing holds them in the pocket but ~0.14 kN of friction, and at a ~ 5.7 m/s2 the
  0.58 kN drive walked a joist out of its pocket at 17-24 mm in every run (pocket_open_mm
  jumping to ~280 mm). The joists now carry only their weight and are dragged by the pocket,
  as in the test. The spring force on T_B (block.force.app, -10350 N in y) is untouched:
  po_push only sets the z component of masonry blocks.
  Consequences:
    * F_tot (base + top reaction) keeps its meaning and is the comparison quantity;
    * 'a' now normalises by the masonry mass only (plateau ~10 kN -> a ~ 7.7 m/s2); for the
      regime boundaries against table PTA keep F / (m_wall + m_joists), i.e. ~0.58 g;
    * a pocket exit (pocket_open_mm > POCKET_EXIT_MM) now aborts the run with a message
      instead of silently continuing on a different structure.

v4 diagnostics kept: joist_N_kN / joist_Vz_kN (joist-support contacts), F_base_only_kN,
max_open_mm (masonry-masonry joints), pocket_open_mm (joist-masonry pocket contacts); vcheck
slot by elimination; TOP_JOINT "elastic" | "mohr" | "mason". Run both PUSH_DIR.
"""

import itasca as it
import os, csv, math

it.command("python-reset-state false")

# ============================ CONFIG =================================
RESTORE_SAV = "Part_I_MASON_v8_LargeStrain_springforce"   # tension = 0, c = 0.2, Gs = 100 build
OUT_DIR     = "pushover_results_v5_ft0_masonrypush_TopMason"
PUSH_DIR    = +1          # +1 / -1 : run both, separate executions
TOP_JOINT   = "mason"   # "elastic" (v3 default), "mohr", or "mason" (leave the build's joint)
TOPJ_Y      = (2.57, 2.59)
POCKET_EXIT_MM = 50.0     # a pocket contact opening beyond this = a joist has left its pocket -> abort
TENSION_MAX_PA = 0.05e6   # abort if the save's joint tension is above this (bonded build by mistake)

# ---- ascending (Nicolo stepping) ----
A_START     = 0.5
DA_ASC      = 0.5
A_MAX       = 30.0
RATIO       = 1e-5
CYC_BUDGET  = 150000
CHUNK_ASC   = 5000
D_RUNAWAY_MM = 3.0
# ---- descending (servo) ----
DD_DESC     = 1.0
D_STOP_MM   = 80.0
V_EXIT_KN   = 0.5
DESC_POINTS_MAX = 120
CHUNK       = 500
KP          = 0.30
KV          = 3.0
DA_SERVO_MAX = 0.50
D_TOL       = 0.05
N_STILL     = 4
SERVO_CHUNKS_MAX = 400
LAMBDA_PRED = 1.30

os.makedirs(OUT_DIR, exist_ok=True)

# ============================ RESTORE ================================
if not (os.path.isfile(RESTORE_SAV) or os.path.isfile(RESTORE_SAV + ".sav")):
    raise RuntimeError("missing save file: {}(.sav)".format(RESTORE_SAV))
it.command("model restore '{}'".format(RESTORE_SAV))
it.command("python-reset-state false")
it.command("model large-strain on")       # the save is already LS; guard only
it.command("model gravity 0 -9.81 0")     # v5: gravity stays vertical throughout
it.command("block mech damp local 0.9")   # quasi-static (never seed a dynamic run from these saves)

# build guard
try:
    _ft = float(it.fish.get("tension"))
    print("joint tension in this save: {:.3f} MPa".format(_ft / 1e6))
    if _ft > TENSION_MAX_PA:
        raise RuntimeError("save has joint tension {:.3f} MPa: bonded build, not the f_t = 0 build".format(_ft / 1e6))
except RuntimeError:
    raise
except Exception as _e:
    print("WARNING: cannot read FISH 'tension' ({}); build not verified".format(_e))

# ---- top joint -------------------------------------------------------
if TOP_JOINT == "elastic":
    it.command("block contact jmodel assign elastic "
               "range group-intersection 'T_B' 'Masonry'")
    it.command("block contact property stiffness-normal @kn stiffness-shear @ks "
               "range group-intersection 'T_B' 'Masonry'")
elif TOP_JOINT == "mohr":
    it.command("block contact jmodel assign mohr "
               "range group-intersection 'T_B' 'Masonry'")
    it.command("block contact property stiffness-normal @kn stiffness-shear @ks "
               "tension 0 cohesion 0 friction 35 "
               "range group-intersection 'T_B' 'Masonry'")
else:
    pass                                   # "mason": the build's joint stays
it.command("block contact group 'TopJ' range pos-y {} {}".format(*TOPJ_Y))

# vcheck slot by elimination
it.command("""
block contact group "MasJ" slot "vcheck"
block contact group "Reaction_base"  slot "vcheck" range pos-y -0.01 0.01
block contact group "Reaction_topTB" slot "vcheck" range group-intersection 'T_B' 'Masonry'
block contact group "Pocket"         slot "vcheck" range group-intersection 'Joist' 'Masonry'
block contact group "Reaction_joist" slot "vcheck" range group-intersection 'Joist' 'Joist_S1'
block contact group "Reaction_joist" slot "vcheck" range group-intersection 'Joist' 'Joist_S2'
""")
it.command("[global push_sign = {:d}]".format(int(PUSH_DIR)))

# top-joint reaction (unchanged from v3)
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

# diagnostics: joist-support forces, base-only shear, largest joint and pocket openings
it.command("""
fish define po_diag
    local jn = 0.0
    local jv = 0.0
    local fb = 0.0
    local uo = 0.0
    local up = 0.0
    local sc = null
    local un = 0.0
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, 'Reaction_joist', 'vcheck') then
            loop foreach sc block.contact.subcontactlist(cx)
                jn = jn + block.subcontact.force.norm(sc)
                jv = jv + block.subcontact.force.shear.z(sc)
            endloop
        endif
        if block.contact.isgroup(cx, 'Reaction_base', 'vcheck') then
            loop foreach sc block.contact.subcontactlist(cx)
                fb = fb + block.subcontact.force.shear.z(sc)
            endloop
        endif
        if block.contact.isgroup(cx, 'MasJ', 'vcheck') then
            loop foreach sc block.contact.subcontactlist(cx)
                un = block.subcontact.disp.norm(sc)
                if un > uo then
                    uo = un
                endif
            endloop
        endif
        if block.contact.isgroup(cx, 'Pocket', 'vcheck') then
            loop foreach sc block.contact.subcontactlist(cx)
                un = block.subcontact.disp.norm(sc)
                if un > up then
                    up = un
                endif
            endloop
        endif
    endloop
    global joist_N_kN    = jn / 1000.0
    global joist_Vz_kN   = jv / 1000.0
    global base_only_kN  = fb / 1000.0
    global max_open_mm   = uo * 1000.0
    global pocket_open_mm = up * 1000.0
    po_diag = jn
end
""")

# v5: the push -- applied force a*m on masonry blocks only (z component); gravity untouched.
# block.force.app is SET, so the masonry blocks' applied force is exactly (0, 0, a*m).
# T_B keeps its (0, -10350, 0) spring force because it is not in group 'Masonry'.
it.command("""
fish define po_push(a_push)
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry') then
            block.force.app(b) = vector(0.0, 0.0, a_push * block.mass(b))
        endif
    endloop
end
""")

# masonry mass, participating mass (mass from the blocks themselves)
it.command("""
fish define wall_mass_part
    global m_tot = 0.0
    global m_part = 0.0
    global wall_maxvel = 0.0
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry') then
            local mb = block.mass(b)
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
""")

it.command("[po_push(0.0)]")
it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CYC_BUDGET))

# ============================ BASELINE ===============================
it.command("@Record_Disp")
it.command("@wall_mass_part")
it.command("@po_diag")
d0  = 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1000.0
F0b = it.fish.call_function("cstav")
F0t = it.fish.call_function("topj_shear")
F0bo = it.fish.get("base_only_kN")
m_tot = it.fish.get("m_tot")
print("settled baseline: d = {:.4f} mm, F_base = {:.4f} kN, F_top = {:.4f} kN, "
      "masonry mass = {:.1f} kg (push normalisation)".format(d0, F0b, F0t, m_tot))
print("joist supports: N = {:.3f} kN (static 3.22 expected), max joint opening = {:.4f} mm, "
      "pocket opening = {:.4f} mm".format(it.fish.get("joist_N_kN"), it.fish.get("max_open_mm"),
                                          it.fish.get("pocket_open_mm")))
print("push mode: applied force a*m on MASONRY only; gravity vertical; joists carry weight only")

# ============================ HELPERS ================================
def get_cycles():
    try:
        it.command("[global _cyc_now = mech.step]")
        return int(it.fish.get("_cyc_now"))
    except Exception:
        return int(it.cycle())

def set_a(a_mps2):
    it.command("[po_push({:.6f})]".format(PUSH_DIR * a_mps2))

def ctrl_disp():
    it.command("@Record_Disp")
    return 0.5 * (it.fish.get("Top_Quarter_A_Disp") +
                  it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0

def read_state():
    d = ctrl_disp()
    it.command("@wall_mass_part")
    it.command("@po_diag")
    Fb = it.fish.call_function("cstav") - F0b
    Ft = it.fish.call_function("topj_shear") - F0t
    diag = (it.fish.get("joist_N_kN"), it.fish.get("joist_Vz_kN"),
            it.fish.get("base_only_kN") - F0bo, it.fish.get("max_open_mm"),
            it.fish.get("pocket_open_mm"))
    return d, Fb, Ft, it.fish.get("m_part"), it.fish.get("wall_maxvel"), diag

def conv_of(a, Ftot):
    tgt = a * m_tot / 1000.0
    return abs(tgt - abs(Ftot)) / tgt if tgt > 0 else 0.0

def pocket_exit(dg):
    return dg[4] > POCKET_EXIT_MM

def solve_step(a_mps2):
    set_a(a_mps2)
    c0 = get_cycles()
    d_start = ctrl_disp()
    while get_cycles() - c0 < CYC_BUDGET:
        c1 = get_cycles()
        it.command("model solve ratio {:g} cycles {:d}".format(RATIO, CHUNK_ASC))
        if get_cycles() - c1 < CHUNK_ASC:
            return True, False
        if (ctrl_disp() - d_start) * PUSH_DIR > D_RUNAWAY_MM:
            return False, True
    return False, False

def nicolo_predictor(a, Ftot, m_part):
    V_N = abs(Ftot) * 1e3
    m_p = max(m_part, 0.05 * m_tot)
    return max(LAMBDA_PRED * (V_N - a * (m_tot - m_p)) / m_p, 0.0)

def servo_to(d_target, a):
    still = 0
    d_prev = ctrl_disp()
    for j in range(SERVO_CHUNKS_MAX):
        set_a(a)
        it.command("model cycle {:d}".format(CHUNK))
        d = ctrl_disp()
        err = (d_target - d) * PUSH_DIR
        vel = (d - d_prev) * PUSH_DIR
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
HEADER = ["phase", "a_mps2", "d_ctrl_mm", "F_base_kN", "F_top_kN", "F_tot_kN",
          "conv", "converged", "m_part_kg", "a_pred_mps2", "settle_chunks",
          "joist_N_kN", "joist_Vz_kN", "F_base_only_kN", "max_open_mm", "pocket_open_mm"]

def flush():
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# push: applied force a*m on masonry only (block.force.app z); gravity vertical; "
                    "m_tot {:.1f} kg; save {}; top joint {}".format(m_tot, RESTORE_SAV, TOP_JOINT)])
        w.writerow(HEADER); w.writerows(rows)

def row(phase, a, d, Fb, Ft, Ftot, ok, m_part, a_pred, nch, dg):
    return (phase, round(PUSH_DIR * a, 4), round(d, 4), round(Fb, 4), round(Ft, 4),
            round(Ftot, 4), round(conv_of(a, Ftot), 5), 1 if ok else 0,
            round(m_part, 1), round(a_pred, 4), nch,
            round(dg[0], 4), round(dg[1], 4), round(dg[2], 4), round(dg[3], 4), round(dg[4], 4))

aborted = False

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
    d, Fb, Ft, m_part, vmax, dg = read_state()
    Ftot = Fb + Ft
    a_pred = nicolo_predictor(a, Ftot, m_part)
    rows.append(row("asc", a, d, Fb, Ft, Ftot, ok, m_part, a_pred, "", dg))
    flush()
    print("  asc a {:6.2f}  d {:8.3f} mm  Fb {:7.3f}  Ft {:7.3f}  Ftot {:8.3f} kN  conv {:.4f}  "
          "joistN {:6.3f}  open {:6.3f}  pocket {:6.3f} mm  {}".format(
              PUSH_DIR * a, d, Fb, Ft, Ftot, conv_of(a, Ftot), dg[0], dg[3], dg[4],
              "ok" if ok else ("RUNAWAY" if ran else "UNCONVERGED")))
    if pocket_exit(dg):
        print("  ** ABORT: pocket opening {:.1f} mm > {} mm -- a joist has left its pocket; "
              "points from here on are not the wall's curve.".format(dg[4], POCKET_EXIT_MM))
        aborted = True
        break
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
if peak and not aborted:
    d_sw = ctrl_disp()
    a = a_last_ok
    d_target = d_sw
    for k in range(DESC_POINTS_MAX):
        settled, a, nch = servo_to(d_target, a)
        d, Fb, Ft, m_part, vmax, dg = read_state()
        Ftot = Fb + Ft
        a_pred = nicolo_predictor(a, Ftot, m_part)
        rows.append(row("desc", a, d, Fb, Ft, Ftot, settled, m_part, a_pred, nch, dg))
        flush()
        print("  desc d_t {:7.2f}  d {:8.3f} mm  a {:6.3f}  Ftot {:8.3f} kN  joistN {:6.3f}  "
              "open {:6.3f}  pocket {:6.3f} mm  {} ({} chunks)".format(
                  d_target, d, PUSH_DIR * a, Ftot, dg[0], dg[3], dg[4],
                  "settled" if settled else "NOT settled", nch))
        if pocket_exit(dg):
            print("  ** ABORT: pocket opening {:.1f} mm > {} mm -- a joist has left its pocket; "
                  "curve valid up to the previous point.".format(dg[4], POCKET_EXIT_MM))
            aborted = True
            break
        if abs(Ftot) < V_EXIT_KN and settled:
            print("  descending exit: |F_tot| < {} kN".format(V_EXIT_KN)); break
        if abs(d) > D_STOP_MM:
            print("  descending exit: |d| > {} mm".format(D_STOP_MM)); break
        if a <= 0.0 and (d - d_target) * PUSH_DIR > 5.0 * D_TOL:
            print("  descending exit: a = 0 and the wall still advances -- past the instability point"); break
        d_target = d_target + PUSH_DIR * DD_DESC
    it.command("model save '{}'".format(
        os.path.join(OUT_DIR, "pushover_{}_end".format(
            "pos" if PUSH_DIR > 0 else "neg")).replace("\\", "/")))

flush()
print("-> {}".format(csv_path))
if aborted:
    print("Run ended on a pocket exit. If it happened on the masonry-only push, the joist was dragged "
          "out by the pocket itself, which is a model question, not a loading artefact.")