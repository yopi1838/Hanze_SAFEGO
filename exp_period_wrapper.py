# -*- coding: ascii -*-
"""Staged wrapper around exp_period_evolution.analyze_run with JSON cache.
Usage: python exp_period_wrapper.py TEST_NO [EXP_DIR] LO HI
       python exp_period_wrapper.py finalize   (applies monotonicity, CSVs)

EXP_DIR defaults to ./EXP_DATA (override: SAFEGO_EXP_DATA).
State is cached in ./.cache/period_state.json.
"""
import sys, json, csv
from pathlib import Path
import safego_paths as sp

STATE = sp.state_file("period")
HERE = sp.ROOT
OUT = sp.postproc_dir()

def load():
    return json.loads(STATE.read_text()) if STATE.exists() else {"9": {}, "12": {}}

if sys.argv[1] == "finalize":
    st = load()
    for test in ("9", "12"):
        rows = [st[test][k] for k in sorted(st[test], key=int)]
        # monotonicity constraint (same as exp_period_evolution.run_all)
        T_prev, ncorr = None, 0
        for r in rows:
            T_raw = r.get("T1_best")
            r["T1_raw"] = T_raw
            if T_raw is None:
                r["T1_best"] = T_prev; r["corrected"] = True; ncorr += 1
                continue
            if T_prev is not None and T_raw < T_prev * 0.97:
                r["T1_best"] = T_prev; r["corrected"] = True; ncorr += 1
            else:
                T_prev = T_raw; r["corrected"] = False
        print("test %s: %d corrected" % (test, ncorr))
        fout = OUT / ("exp_Test%s_period_psd.csv" % test)
        with open(str(fout), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["run", "T1_raw_s", "T1_best_s", "corrected"])
            for r in rows:
                w.writerow([r["run"],
                            round(r["T1_raw"], 5) if r.get("T1_raw") else "",
                            round(r["T1_best"], 5) if r.get("T1_best") else "",
                            "Y" if r["corrected"] else ""])
        print("->", fout)
    sys.exit()

import exp_period_evolution as pe
# EXP_DIR is optional: "TEST_NO LO HI" or "TEST_NO EXP_DIR LO HI"
if len(sys.argv) >= 5:
    test_no, exp_dir, lo, hi = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
else:
    test_no, exp_dir, lo, hi = sys.argv[1], str(sp.exp_data_dir()), int(sys.argv[2]), int(sys.argv[3])
st = load()
for rn in range(lo, hi + 1):
    if str(rn) in st[test_no]:
        continue
    r = pe.analyze_run(rn, exp_dir, test_no=int(test_no))
    if r is None:
        print("run %d: missing" % rn); continue
    keep = {"run": rn, "T1_best": r.get("T1_best"),
            "f1_best": r.get("f1_best"), "f1_source": r.get("f1_source"),
            "t_fv": r.get("t_free_vib_start")}
    st[test_no][str(rn)] = keep
    print("test %s run %2d: T1=%s  (src=%s)" % (
        test_no, rn,
        "%.4f s" % r["T1_best"] if r.get("T1_best") else "--",
        r.get("f1_source")))
    STATE.write_text(json.dumps(st))
