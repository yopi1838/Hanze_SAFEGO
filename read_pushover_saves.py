# -*- coding: ascii -*-
"""
read_pushover_saves.py -- per-floor joist lever state read from pushover saves (3DEC 9; restores only, no cycling).

Why (29 Oct): the slip-work pushover (f_t = 0.05) logged the far-support load of the joists summed over both floors
(3.2 kN at rest -> 4.9 kN positive push, 0.9 kN negative). The lever model needs it PER FLOOR, with the joists'
rotations, to tell which floor is pressed and which lifts off. This reads it from the saves the run left behind.

For each save in SAVES, after 'model restore':
    support_low_N_kN / support_up_N_kN        far-support load per floor (Joist - Joist_S1 / Joist_S2, sum of N)
    pocket_<low|up>_<top|bot|back|side>_N_kN  normal force on each pocket face (ceiling, floor, back = joist end, sides)
    pocket_<low|up>_back_open_max_mm          largest normal opening of the back-face contacts (joist end gap)
    rot_<joist|mas>_<low|up>_mrad, rot_rel_*  rotations as in slip_work_tracker.py (gridpoint displacement since the
                                              build, so compare an end save with its start save: the change is the push)
Classification and rotations come from slip_work_tracker.py (same groups, same faces, same gridpoints).
Output: SAVES_CSV (one row per save) and a printed table with end - start per direction.

USAGE (3DEC 9, working folder = the one with the saves and slip_work_tracker.py):
    program call 'read_pushover_saves.py'
"""
import os, sys, csv, importlib

# ============================ CONFIG =================================
RUN_DIR  = "pushover_results_v5_ft_0p050_TopMason_slipwork"
SAVES = [   # (label, save, direction)
    ("pos_start", "Part_I_MASON_v8_LargeStrain_springforce_ft_0p050", "pos"),
    ("pos_end",   os.path.join(RUN_DIR, "pushover_pos_end"), "pos"),
    ("neg_start", "Part_I_MASON_v8_LargeStrain_springforce_ft_0p050_n", "neg"),
    ("neg_end",   os.path.join(RUN_DIR, "pushover_neg_end"), "neg"),
]
SAVES_CSV = os.path.join(RUN_DIR, "lever_state_from_saves.csv")
# ======================================================================

import itasca as it
sys.path.insert(0, os.getcwd())
import slip_work_tracker as swt
swt = importlib.reload(swt)

KEYS = ["support_low", "support_up"] + ["pocket_{}_{}".format(p, f) for p in ("low", "up") for f in ("top", "bot", "back", "side")]


def read_state():
    swt._set_groups()
    N = dict((k, 0.0) for k in KEYS)
    nsc = dict((k, 0) for k in KEYS)
    opn = {"pocket_low_back": 0.0, "pocket_up_back": 0.0}
    for sc in it.block.subcontact.list():
        r = swt._region(sc, installed=False)
        if r not in ("support", "pocket_low", "pocket_up"):
            continue
        k = swt._sub(r, sc)
        if k not in N:
            continue
        N[k] += float(sc.force_norm()) / 1e3
        nsc[k] += 1
        if k in opn:
            try:
                opn[k] = max(opn[k], float(sc.disp_norm()) * 1e3)
            except Exception:
                pass
    swt._rot_ok["ok"] = None                 # gridpoints are looked up again in every save
    rot = swt.rotations()
    row = {}
    for k in KEYS:
        row[k + "_N_kN"] = N[k]
        row[k + "_n"] = nsc[k]
    for k, v in opn.items():
        row[k + "_open_max_mm"] = v
    for k, v in rot.items():
        row["rot_{}_mrad".format(k)] = v
    try:
        it.command("@Record_Disp")
        row["d_sensor_mm"] = 0.5 * (it.fish.get("Top_Quarter_A_Disp") + it.fish.get("Top_Quarter_B_Disp")) * 1e3
    except Exception:
        row["d_sensor_mm"] = float("nan")
    return row


rows = []
for label, sv, br in SAVES:
    if not (os.path.isfile(sv) or os.path.isfile(sv + ".sav")):
        print("[read_pushover_saves] ! {} missing -- skipped".format(sv)); continue
    it.command("model restore '{}'".format(sv.replace("\\", "/")))
    it.command("python-reset-state false")
    r = read_state(); r["label"] = label; r["save"] = sv; r["dir"] = br
    rows.append(r)
    print("[read_pushover_saves] {}: supports low/up {:.2f} / {:.2f} kN; lower pocket ceiling/floor/back/side "
          "{:.1f} / {:.1f} / {:.1f} / {:.1f} kN; rel. rotation low/up {:.2f} / {:.2f} mrad".format(
              label, r["support_low_N_kN"], r["support_up_N_kN"], r["pocket_low_top_N_kN"], r["pocket_low_bot_N_kN"],
              r["pocket_low_back_N_kN"], r["pocket_low_side_N_kN"], r["rot_rel_low_mrad"], r["rot_rel_up_mrad"]))

if rows:
    keys = ["label", "dir", "save", "d_sensor_mm"] + [k for k in rows[0] if k not in ("label", "dir", "save", "d_sensor_mm")]
    os.makedirs(os.path.dirname(SAVES_CSV) or ".", exist_ok=True)
    with open(SAVES_CSV, "w", newline="\n") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rows:
            w.writerow([("{:.4f}".format(r[k]) if isinstance(r.get(k), float) else r.get(k, "")) for k in keys])
    print("[read_pushover_saves] -> {}".format(SAVES_CSV))
    by = dict((r["label"], r) for r in rows)
    for br in ("pos", "neg"):
        a, b = by.get(br + "_start"), by.get(br + "_end")
        if a and b:
            print("  {} push, end - start: support low {:+.2f} kN, up {:+.2f} kN | joist rotation low {:+.1f}, up {:+.1f} mrad "
                  "| relative to pocket low {:+.1f}, up {:+.1f} mrad".format(
                      br, b["support_low_N_kN"] - a["support_low_N_kN"], b["support_up_N_kN"] - a["support_up_N_kN"],
                      b["rot_joist_low_mrad"] - a["rot_joist_low_mrad"], b["rot_joist_up_mrad"] - a["rot_joist_up_mrad"],
                      b["rot_rel_low_mrad"] - a["rot_rel_low_mrad"], b["rot_rel_up_mrad"] - a["rot_rel_up_mrad"]))
