# -*- coding: ascii -*-
"""
check_frame_zup.py -- does the Z-UP build reproduce the old build? (3DEC 9; restores the two saves, small push only)

The Z-up frame (Geo_Prep_zup.dat, ANALYSIS_PART_I_MASON_zup.dat) is the same wall rotated:
    X = 1.292 - x (in-plane),  Y = z (out-of-plane, +towards the joists),  Z = y (vertical).
Physics does not care about the frame, so every quantity below must agree between the two saves to solver tolerance
(a few per mille). A mismatch points at a direction that was not converted.

For each save, after 'model restore':
  A. settled state: number of blocks / contacts / subcontacts; masonry, joist, top-beam mass; vertical reactions at
     the base, the joist supports (per floor) and the top joint; normal force in the lower and upper pockets;
     positions of the five monitoring points (mapped back to the old frame for the comparison).
  B. a small elastic out-of-plane push: a = PUSH_A on every masonry block (force a*m, OOP, +towards the joists, the
     v5 pushover's load), solved to equilibrium; then the OOP displacement at the 2.06 m sensors (Record_Disp's
     Top_Quarter_A/B, same material points in both builds), K = a m_tot / d, T1 = 2 pi sqrt(1200 kg / K), the base
     shear (cstav), and the pocket and support forces again.
The saves are only read; nothing is written over them. Output: CHECK_CSV and a side-by-side table.

USAGE (3DEC 9, folder with both saves):   program call 'check_frame_zup.py'
"""
import os, csv, math
import itasca as it

CASES = [   # (label, save, vertical axis, OOP axis)
    ("old", "Part_I_MASON_v8_LargeStrain_springforce", "y", "z"),
    ("zup", "Part_I_MASON_v8_LargeStrain_springforce_zup", "z", "y"),
]
PUSH_A = 1.0                 # m/s2, elastic range of the f_t = 0 build (rocking needs ~7.7)
SOLVE = "model solve ratio 1e-6 cycles 200000"
M_EFF = 1200.0
CHECK_CSV = "check_frame_zup.csv"


def fish_defs(v, o):
    """FISH helpers for vertical axis v and OOP axis o"""
    it.command("""
fish define fchk_sum(gname, gslot, lo, hi)
    local f = 0.0
    loop foreach local cx block.contact.list
        if block.contact.isgroup(cx, gname, gslot)
            loop foreach local sc block.contact.subcontactlist(cx)
                local yv = comp.{v}(block.subcontact.pos(sc))
                if yv >= lo
                    if yv <= hi
                        f = f + block.subcontact.force.norm(sc)
                    endif
                endif
            endloop
        endif
    endloop
    fchk_sum = f / 1000.0
end
fish define fchk_mass(gname)
    local m = 0.0
    loop foreach local b block.list
        if block.isgroup(b, gname)
            m = m + block.mass(b)
        endif
    endloop
    fchk_mass = m
end
fish define fchk_counts
    local nb = 0
    loop foreach local b block.list
        nb = nb + 1
    endloop
    global fchk_nb = nb
    local nc = 0
    local ns = 0
    loop foreach local cx block.contact.list
        nc = nc + 1
        loop foreach local sc block.contact.subcontactlist(cx)
            ns = ns + 1
        endloop
    endloop
    global fchk_nc = nc
    global fchk_ns = ns
end
fish define fchk_push(a)
    loop foreach local b block.list
        if block.isgroup(b, 'Masonry')
            block.force.app(b) = vector({pushvec})
        endif
    endloop
end
""".format(v=v, pushvec=", ".join("a * block.mass(b)" if ax == o else "0.0" for ax in "xyz")))


def fval(expr):
    """value of a FISH expression (via a global; no Python-side argument passing)"""
    it.command("[global _fchk_v = {}]".format(expr))
    return float(it.fish.get("_fchk_v"))


def state(v, o):
    it.command("block contact group 'fchk_pocket' slot 'fchk' range group-intersection 'Joist' 'Masonry'")
    it.command("block contact group 'fchk_s1' slot 'fchk' range group-intersection 'Joist' 'Joist_S1'")
    it.command("block contact group 'fchk_s2' slot 'fchk' range group-intersection 'Joist' 'Joist_S2'")
    it.command("block contact group 'fchk_top' slot 'fchk' range group-intersection 'T_B' 'Masonry'")
    r = {}
    r["base_N_kN"] = fval("fchk_sum('Reaction_base', 'vcheck', -0.01, 0.01)")
    r["support_low_N_kN"] = fval("fchk_sum('fchk_s1', 'fchk', -10.0, 10.0)")
    r["support_up_N_kN"] = fval("fchk_sum('fchk_s2', 'fchk', -10.0, 10.0)")
    r["top_joint_N_kN"] = fval("fchk_sum('fchk_top', 'fchk', -10.0, 10.0)")
    r["pocket_low_N_kN"] = fval("fchk_sum('fchk_pocket', 'fchk', -10.0, 1.3)")
    r["pocket_up_N_kN"] = fval("fchk_sum('fchk_pocket', 'fchk', 1.3, 10.0)")
    return r


def sensors():
    it.command("@Record_Disp")
    return dict(dA_mm=float(it.fish.get("Top_Quarter_A_Disp")) * 1e3, dB_mm=float(it.fish.get("Top_Quarter_B_Disp")) * 1e3,
                dmid_mm=float(it.fish.get("Mid_Disp")) * 1e3, dbot_mm=float(it.fish.get("Bot_Quarter_Disp")) * 1e3)


def to_old(p, v, o):
    """a point of this frame in old coordinates (x, y_vertical, z_oop)"""
    if v == "y":
        return p
    return (1.292 - p[0], p[2], p[1])


rows = []
for label, sv, v, o in CASES:
    if not (os.path.isfile(sv) or os.path.isfile(sv + ".sav")):
        print("[check_frame_zup] ! {} missing -- skipped".format(sv)); continue
    it.command("model restore '{}'".format(sv))
    it.command("python-reset-state false")
    fish_defs(v, o)
    it.command("[fchk_counts]")
    r = dict(label=label, save=sv, n_blocks=int(it.fish.get("fchk_nb")), n_contacts=int(it.fish.get("fchk_nc")),
             n_subcontacts=int(it.fish.get("fchk_ns")))
    for g in ("Masonry", "Joist", "T_B"):
        r["mass_{}_kg".format(g)] = fval("fchk_mass('{}')".format(g))
    r.update(state(v, o))
    it.command("@Record_Disp")
    for nm in ("monitoring_point_3a", "monitoring_point_3b", "monitoring_point_2"):
        it.command("[global _fq = block.gp.pos({})]".format(nm))
        q = it.fish.get("_fq")
        try:
            q = (float(q[0]), float(q[1]), float(q[2]))
        except Exception:
            q = (float(q.x()), float(q.y()), float(q.z()))
        po = to_old(q, v, o)
        r["{}_old_xyz".format(nm)] = "({:.4f} {:.4f} {:.4f})".format(*po)
    s0 = sensors()
    F0 = float(it.fish.call_function("cstav"))
    # elastic OOP push
    it.command("[fchk_push({:.6f})]".format(PUSH_A))
    it.command(SOLVE)
    s1 = sensors()
    d = 0.5 * ((s1["dA_mm"] - s0["dA_mm"]) + (s1["dB_mm"] - s0["dB_mm"]))
    m_tot = r["mass_Masonry_kg"]
    K = PUSH_A * m_tot / 1e3 / d if d else float("nan")             # kN/mm
    r.update(push_d_mm=d, push_dmid_mm=s1["dmid_mm"] - s0["dmid_mm"], push_dbot_mm=s1["dbot_mm"] - s0["dbot_mm"],
             K_kN_per_mm=K, T1_s=2 * math.pi * math.sqrt(M_EFF / (K * 1e6)) if K == K and K > 0 else float("nan"),
             push_base_shear_kN=float(it.fish.call_function("cstav")) - F0)
    for k, val in state(v, o).items():
        r["push_" + k] = val
    rows.append(r)
    print("[check_frame_zup] {} done".format(label))

if rows:
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(CHECK_CSV, "w", newline="\n") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rows:
            w.writerow([("{:.5g}".format(r[k]) if isinstance(r.get(k), float) else r.get(k, "")) for k in keys])
    print("\n{:30s}".format("quantity") + "".join("{:>22s}".format(r["label"]) for r in rows) + "   rel. diff")
    for k in keys:
        if k in ("label", "save"):
            continue
        vals = [r.get(k) for r in rows]
        line = "{:30s}".format(k) + "".join("{:>22s}".format("{:.5g}".format(x) if isinstance(x, float) else str(x)) for x in vals)
        if len(vals) == 2 and all(isinstance(x, (int, float)) for x in vals) and abs(vals[0]) > 1e-9:
            rd = (vals[1] - vals[0]) / abs(vals[0])
            line += "   {:+.2%}{}".format(rd, "   <-- CHECK" if abs(rd) > 0.01 else "")
        print(line)
    print("-> " + CHECK_CSV)
