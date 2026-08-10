# -*- coding: ascii -*-
"""
Pre-generate 3DEC velocity tables from the full measured shake-table records.
=============================================================================
This is the OUTSIDE-3DEC half of the full-record strategy (Strategy F). It
turns each Test{T}Run{N}.xlsx into a conditioned velocity table that
strategy_F_3dec.py imports and applies. The split exists because reading xlsx
needs openpyxl, which 3DEC's embedded interpreter does not ship, and the raw
records are hundreds of MB. Generate the tables once here, then the driver only
touches small text files.

The conditioning is exactly exp_table_to_vel.py (differentiate channel 5, window
to the significant duration, remove the mean velocity, cosine-taper both ends),
imported rather than re-implemented so the two cannot drift.

OUTPUT (into --out-dir, default the repo root)
    vel_exp_T{T}_run{NN}.txt          one 3DEC velocity table per run
    record_tables_manifest_US{n}.csv  run, table_file, dur_s, n_samples, dt,
                                       pgd_amp_mm, pgv_mps, fidelity

USAGE
    python make_record_tables.py                    # US-1 (Test 9), runs 1-24
    python make_record_tables.py --test 12          # US-2 (Test 12), 1-25
    python make_record_tables.py --runs 22 24       # just the FR76 runs
    python make_record_tables.py --out-dir stratF_results_US1

Style note: runs outside 3DEC, but kept .format()/ascii for consistency.
"""
import argparse
import csv
from pathlib import Path

import exp_table_to_vel as etv

try:
    import safego_paths as sp
except ImportError:
    sp = None


class _Args(object):
    """Minimal stand-in for the argparse namespace exp_table_to_vel.build wants."""
    def __init__(self, taper, lowpass, from_acc):
        self.t0 = None
        self.t1 = None
        self.taper = taper
        self.lowpass = lowpass
        self.from_acc = from_acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", type=int, default=9, choices=(9, 12))
    ap.add_argument("--exp-dir", default=None)
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--runs", type=int, nargs=2, default=None,
                    metavar=("LO", "HI"))
    ap.add_argument("--taper", type=float, default=0.25)
    ap.add_argument("--lowpass", type=float, default=None)
    ap.add_argument("--from-acc", action="store_true")
    args = ap.parse_args()

    if args.exp_dir:
        exp_dir = Path(args.exp_dir)
    elif sp is not None:
        exp_dir = sp.exp_data_dir()
    else:
        exp_dir = Path("EXP_DATA")

    out_dir = Path(args.out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    us = 1 if args.test == 9 else 2
    lo, hi = args.runs if args.runs else (1, 25)
    build_args = _Args(args.taper, args.lowpass, args.from_acc)

    manifest = []
    for rn in range(lo, hi + 1):
        src = exp_dir / "Test{}Run{}.xlsx".format(args.test, rn)
        if not src.is_file():
            print("  missing, skipped: {}".format(src.name))
            continue
        t, d, a = etv.read_run(src)
        # Window on the UNION of acceleration-energy and displacement-energy
        # significant durations. Acceleration alone misses a large low-
        # frequency table swing (US-1 run 18: an 18 mm excursion at t~11 s that
        # the acc window drops, reporting 2 mm). For a full-record replay the
        # displacement swing is the part that matters most to a rocking wall,
        # so it must be inside the window.
        # Significant duration = the standard Arias 5-95% window on the table
        # acceleration, then extended if necessary so the single largest
        # displacement swing is guaranteed to be inside it (with a 1 s margin).
        # Arias alone is the accepted significant-duration measure but is high-
        # frequency biased and drops US-1 run 18's low-frequency 18 mm swing;
        # the guard restores it without the drift-inflation that a displacement-
        # energy window suffers (which stretched run 24 to ~88 s).
        import numpy as _np
        i0, i1 = etv.arias_window(a, 0.05, 0.95)
        _dt = float(_np.median(_np.diff(t)))
        _dd = d - _np.polyval(_np.polyfit(t, d, 1), t)   # detrended
        _ipk = int(_np.argmax(_np.abs(_dd - _dd[0])))
        _m = int(round(1.0 / _dt))
        i0 = max(0, min(i0, _ipk - _m))
        i1 = min(len(t) - 1, max(i1, _ipk + _m))
        build_args.t0 = float(t[i0])
        build_args.t1 = float(t[i1])
        tt, v, info = etv.build(t, d, a, build_args)

        name = "exp_T{}_run{:02d}".format(args.test, rn)
        table_file = out_dir / "vel_exp_T{}_run{:02d}.txt".format(args.test, rn)
        etv.write_table(table_file, name, tt, v)

        manifest.append({
            "run": rn,
            "table_file": table_file.name,
            "dur_s": round(info["dur"], 4),
            "n_samples": info["n"],
            "dt": round(info["dt"], 6),
            "pgd_amp_mm": round(info["d_win_pp_mm"] / 2.0, 4),
            "pgv_mps": round(info["v_peak"], 5),
            "fidelity": round(info["fidelity"], 4),
        })
        print("  run {:2d}: {:.1f} s, {} pts  ->  {}  "
              "(PGD {:.2f} mm, fidelity {:.3f})".format(
                  rn, info["dur"], info["n"], table_file.name,
                  info["d_win_pp_mm"] / 2.0, info["fidelity"]))

    if not manifest:
        raise SystemExit("no records found under {}".format(exp_dir))

    man_path = out_dir / "record_tables_manifest_US{}.csv".format(us)
    cols = ["run", "table_file", "dur_s", "n_samples", "dt",
            "pgd_amp_mm", "pgv_mps", "fidelity"]
    with open(str(man_path), "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for m in manifest:
            w.writerow([m[c] for c in cols])

    total = sum(m["dur_s"] for m in manifest)
    print("\n-> {} tables + {}".format(len(manifest), man_path.name))
    print("   total model time to run them all: {:.0f} s "
          "({:.0f} runs x {:.1f} s mean)".format(
              total, len(manifest), total / len(manifest)))


if __name__ == "__main__":
    main()
