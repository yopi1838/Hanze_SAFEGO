# -*- coding: ascii -*-
"""
pushover_continue_zup.py -- continue a finished z-up pushover (pushover_oop_3dec_v5_zup.py) past its 80 mm stop, from
its end save, down to the instability point (F = 0) or a new D_STOP_MM (9 Nov 2026).

Not called directly: a wrapper sets CASE and PUSH_DIR and execs this file (run_continue_<case>_<pos|neg>.py).
3DEC 9, Hanze_SAFEGO folder:  program call 'run_continue_unclamp_pos.py'

WHAT IT RESTORES
    <OUT_DIR>/pushover_<side>_end.sav   written by the driver right after its last read_state (no cycles in between),
                                        so the restored state IS the last row of <OUT_DIR>/pushover_<side>.csv.
    The save holds the model AND all FISH functions/variables (po_push, wall_mass_part, po_diag, topj_shear, cstav,
    Record_Disp, push_sign, the contact groups and the applied force a*m on the masonry). It does NOT hold the driver's
    Python variables, i.e. the baselines d0, F0b, F0t, F0bo taken before the first push. They are rebuilt here from the
    CSV's last row:  baseline = (value read now) - (value in the last row).  The first state read after the restore
    then reproduces the last row exactly; the script prints both so the match can be checked.

WHAT IT DOES
    The driver's descending servo, unchanged (same gains, same 1 mm steps, same settle test), from the last row's
    displacement and push, until
        |F_tot| < V_EXIT_KN with the wall settled  (instability reached), or
        a = 0 and the wall still advances          (past the instability point), or
        |d| > D_STOP_MM, or
        pocket opening > POCKET_EXIT_MM            (a joist has left its pocket: points after it are not the wall).
    POCKET_EXIT_MM: the driver used 50 mm. The positive pushes reach ~40 mm pocket opening at 80 mm and grow ~0.5 mm/mm,
    so 50 mm would stop them near 100 mm. The joists sit about 105 mm deep in the 210 mm wall (joist from Y 0.105 m),
    so 90 mm is used here; check the pocket_open_mm column against that depth before reading the tail.
    Damping: local 0.9 as in the driver, global 0 (quasi-static servo; equilibrium points do not depend on it).

OUTPUT (the original CSV and end save are not touched)
    <OUT_DIR>/pushover_<side>_cont.csv       the continuation rows only (phase "cont")
    <OUT_DIR>/pushover_<side>_full.csv       the original rows + the continuation, same header as the driver
    <OUT_DIR>/pushover_<side>_end2.sav       end save of the continuation
"""
import itasca as it
import os, csv

it.command("python-reset-state false")

# ============================ CONFIG =================================
CASES = {   # CASE -> pushover folder of the finished run
    "unclamp": "pushover_results_v5zup_pk2_ft0p05_unclamp_TopMason",
    "glued":   "pushover_results_v5zup_pk2_ft0p05_TopMason",
    "pK11":    "pushover_results_v5zup_pk2_ft0p05_pK11_TopMason",
}
CASE     = globals().get("CASE", "unclamp")
PUSH_DIR = int(globals().get("PUSH_DIR", +1))     # +1 towards the joists, -1 away
OUT_DIR  = CASES[CASE]
SIDE     = "pos" if PUSH_DIR > 0 else "neg"

D_STOP_MM       = 200.0    # past the mechanism's instability displacement (140-166 mm)
POCKET_EXIT_MM  = 90.0     # see the header (driver: 50 mm)
V_EXIT_KN       = 0.5
DD_DESC         = 1.0
DESC_POINTS_MAX = 160
# servo, as in pushover_oop_3dec_v5_zup.py
A_MAX = 30.0; CHUNK = 500; KP = 0.30; KV = 3.0; DA_SERVO_MAX = 0.50; D_TOL = 0.05; N_STILL = 4
SERVO_CHUNKS_MAX = 400; LAMBDA_PRED = 1.30

END_SAV  = os.path.join(OUT_DIR, "pushover_{}_end".format(SIDE)).replace("\\", "/")
CSV_IN   = os.path.join(OUT_DIR, "pushover_{}.csv".format(SIDE))
CSV_CONT = os.path.join(OUT_DIR, "pushover_{}_cont.csv".format(SIDE))
CSV_FULL = os.path.join(OUT_DIR, "pushover_{}_full.csv".format(SIDE))
END2_SAV = os.path.join(OUT_DIR, "pushover_{}_end2".format(SIDE)).replace("\\", "/")

# ============================ LAST ROW OF THE FINISHED RUN ===========
if not (os.path.isfile(END_SAV + ".sav") or os.path.isfile(END_SAV)):
    raise RuntimeError("missing end save {}.sav -- the run did not finish, or the save stayed on another machine".format(END_SAV))
if not os.path.isfile(CSV_IN):
    raise RuntimeError("missing {}".format(CSV_IN))
with open(CSV_IN) as f:
    lines = f.read().splitlines()
HEAD_LINE = lines[0]
data = list(csv.DictReader(l for l in lines[1:] if l.strip()))
HEADER = list(data[0].keys())
last = data[-1]
if last["phase"] != "desc":
    raise RuntimeError("{}: last row is phase '{}', not 'desc' -- nothing to continue".format(CSV_IN, last["phase"]))
d_last = float(last["d_ctrl_mm"]); a_last = abs(float(last["a_mps2"]))
Fb_last, Ft_last = float(last["F_base_kN"]), float(last["F_top_kN"]); Fbo_last = float(last["F_base_only_kN"])
print("continue {} {}: last row d {:.3f} mm, a {:.4f} m/s2, F_tot {:.3f} kN, pocket {:.2f} mm".format(
    CASE, SIDE, d_last, a_last, float(last["F_tot_kN"]), float(last["pocket_open_mm"])))

# ============================ RESTORE ================================
it.command("model restore '{}'".format(END_SAV))
it.command("python-reset-state false")
it.command("model large-strain on")
it.command("model gravity 0 0 -9.81")
it.command("block mech damp global 0.0")
it.command("block mech damp local 0.9")
it.command("[global push_sign = {:d}]".format(PUSH_DIR))

# ============================ BASELINES FROM THE LAST ROW ============
it.command("@Record_Disp")
it.command("@wall_mass_part")
it.command("@po_diag")
d_abs = 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1000.0
d0 = d_abs - d_last
F0b = it.fish.call_function("cstav") - Fb_last
F0t = it.fish.call_function("topj_shear") - Ft_last
F0bo = it.fish.get("base_only_kN") - Fbo_last
m_tot = it.fish.get("m_tot")
print("baselines rebuilt: d0 {:.4f} mm, F0b {:.4f}, F0t {:.4f}, F0bo {:.4f} kN; masonry mass {:.1f} kg".format(
    d0, F0b, F0t, F0bo, m_tot))


# ============================ HELPERS (as the driver) ================
def set_a(a_mps2):
    it.command("[po_push({:.6f})]".format(PUSH_DIR * a_mps2))


def ctrl_disp():
    it.command("@Record_Disp")
    return 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1000.0 - d0


def read_state():
    d = ctrl_disp()
    it.command("@wall_mass_part")
    it.command("@po_diag")
    Fb = it.fish.call_function("cstav") - F0b
    Ft = it.fish.call_function("topj_shear") - F0t
    diag = (it.fish.get("joist_N_kN"), it.fish.get("joist_Vz_kN"), it.fish.get("base_only_kN") - F0bo,
            it.fish.get("max_open_mm"), it.fish.get("pocket_open_mm"))
    return d, Fb, Ft, it.fish.get("m_part"), it.fish.get("wall_maxvel"), diag


def conv_of(a, Ftot):
    tgt = a * m_tot / 1000.0
    return abs(tgt - abs(Ftot)) / tgt if tgt > 0 else 0.0


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


def row(phase, a, d, Fb, Ft, Ftot, ok, m_part, a_pred, nch, dg):
    return (phase, round(PUSH_DIR * a, 4), round(d, 4), round(Fb, 4), round(Ft, 4), round(Ftot, 4),
            round(conv_of(a, Ftot), 5), 1 if ok else 0, round(m_part, 1), round(a_pred, 4), nch,
            round(dg[0], 4), round(dg[1], 4), round(dg[2], 4), round(dg[3], 4), round(dg[4], 4))


rows = []


def flush():
    note = "{} | continued from {} to D_STOP {} mm (pushover_continue_zup.py)".format(HEAD_LINE.strip(), END_SAV, D_STOP_MM)
    with open(CSV_CONT, "w", newline="") as f:
        w = csv.writer(f); w.writerow([note]); w.writerow(HEADER); w.writerows(rows)
    with open(CSV_FULL, "w", newline="") as f:
        w = csv.writer(f); w.writerow([note]); w.writerow(HEADER)
        w.writerows([[r[k] for k in HEADER] for r in data]); w.writerows(rows)


# ============================ CHECK: RESTORED STATE = LAST ROW =======
d, Fb, Ft, m_part, vmax, dg = read_state()
print("restored state: d {:.4f} mm (row {:.4f}), F_tot {:.4f} kN (row {:.4f}), joist N {:.3f} kN (row {}), "
      "pocket {:.3f} mm (row {})".format(d, d_last, Fb + Ft, float(last["F_tot_kN"]), dg[0], last["joist_N_kN"],
                                         dg[4], last["pocket_open_mm"]))

# ============================ DESCENDING, CONTINUED ==================
a = a_last
d_target = d_last + PUSH_DIR * DD_DESC
why = "point limit {}".format(DESC_POINTS_MAX)
for k in range(DESC_POINTS_MAX):
    settled, a, nch = servo_to(d_target, a)
    d, Fb, Ft, m_part, vmax, dg = read_state()
    Ftot = Fb + Ft
    rows.append(row("cont", a, d, Fb, Ft, Ftot, settled, m_part, nicolo_predictor(a, Ftot, m_part), nch, dg))
    flush()
    print("  cont d_t {:7.2f}  d {:8.3f} mm  a {:6.3f}  Ftot {:8.3f} kN  joistN {:6.3f}  pocket {:6.3f} mm  {} ({} chunks)".format(
        d_target, d, PUSH_DIR * a, Ftot, dg[0], dg[4], "settled" if settled else "NOT settled", nch))
    if dg[4] > POCKET_EXIT_MM:
        why = "pocket opening {:.1f} mm > {} mm (a joist has left its pocket; curve valid to the previous point)".format(dg[4], POCKET_EXIT_MM)
        break
    if abs(Ftot) < V_EXIT_KN and settled:
        why = "|F_tot| < {} kN: instability reached".format(V_EXIT_KN); break
    if abs(d) > D_STOP_MM:
        why = "|d| > {} mm".format(D_STOP_MM); break
    if a <= 0.0 and (d - d_target) * PUSH_DIR > 5.0 * D_TOL:
        why = "a = 0 and the wall still advances: past the instability point"; break
    d_target = d_target + PUSH_DIR * DD_DESC
print("  exit: " + why)
it.command("model save '{}'".format(END2_SAV))
flush()
print("-> {}  and  {}".format(CSV_CONT, CSV_FULL))
