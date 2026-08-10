# -*- coding: ascii -*-
"""
Strategy F (FULL) : significant-duration-truncated record tables, runs 1-25.
=============================================================================
The OUTSIDE-3DEC half of the full-sequence fidelity strategy. It turns every
Test{T}Run{N}.xlsx measured shake-table motion into a conditioned velocity
table that strategy_F_full_3dec.py replays in order, run 1 through 25. The
motion applied is the earthquake the specimen actually saw (HU12 / EC40 / FR76,
already scaled and reproduced by the table), so nothing is synthesised.

TRUNCATION (this is what makes it "full" rather than the earlier strategy F)
    The window is NOT the textbook Arias 5-95% band. Two deliberate changes:

      1. START AT t = 0, not at the 5% Arias time. Keeping the record from the
         first sample preserves the low-frequency build-up that a 5% lower cut
         discards. US-1 run 18 is the case that forced this: its single largest
         table swing (17.97 mm) sits early and a 5% cut reported it as 2.13 mm.
         From t = 0 the window returns the full 17.97 mm.

      2. END AT THE VELOCITY ZERO-CROSSING NEAREST THE 95% ARIAS TIME. The 95%
         Arias time marks the end of the strong motion, but cutting exactly
         there generally lands on a non-zero table velocity, which would hit the
         model with a velocity step. Snapping the cut to the closest zero of the
         table velocity makes the applied motion start and end at rest, with no
         step and no artificial jerk. Endpoints are then clamped to exactly 0.

    So the window is [0, t_zero], t_zero = argmin_over_zero_crossings |t - t95|.

DRIFT (default on, --no-drift-correct to disable)
    Cutting at the 95% Arias zero-crossing leaves the table at whatever position
    it had reached (US-1 run 24: -16.9 mm) -- that offset is real, but replaying
    24-25 records in sequence would let the base WALK by the sum of these
    offsets. A zero-mean half-sine baseline, sin(pi t / T_win), is subtracted
    from the velocity to cancel the net displacement. It vanishes at both ends,
    so the zero-velocity endpoints are preserved exactly, and its amplitude is
    <0.4% of peak velocity on every US-1 run, so the motion is not distorted.
    Turn it off to apply the raw truncated motion including the base offset.

VELOCITY SOURCE
    Channel 5, 'Disp - Shake table' (m), central-difference differentiated.
    Robust: differentiating a displacement only amplifies high-frequency noise
    (--lowpass handles it) and cannot integrate an accelerometer offset into a
    drift. Arias intensity for the 95% time is taken from channel 12 (the table
    acceleration), which is what "significant duration" is defined on.

OUTPUT (into --out-dir, default the folder record_tables_full/)
    vel_expfull_T{T}_run{NN}.txt          one 3DEC velocity table per run
    record_tables_full_manifest_US{n}.csv run, table_file, dur_s, n_samples,
                                          dt, pgd_amp_mm, pgv_mps, net_drift_mm,
                                          t95_s, t_cut_s
    Schema is a superset of make_record_tables.py's, so
    strategy_F_full_3dec.py (and postprocess_stratC.py) read it unchanged.

USAGE
    python make_record_tables_full.py                 # US-1 (Test 9), 1-25
    python make_record_tables_full.py --test 12       # US-2 (Test 12), 1-25
    python make_record_tables_full.py --runs 1 25
    python make_record_tables_full.py --no-drift-correct
    python make_record_tables_full.py --lowpass 25    # Hz, if dv/dt is noisy

Style note: Python-2/3 compatible (.format(), no f-strings), ascii, LF newlines,
to match the rest of the codebase.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np

import exp_table_to_vel as etv

try:
    import safego_paths as sp
except ImportError:
    sp = None


def arias_cumfrac(a):
    """Normalised cumulative Arias intensity (0..1) of an acceleration trace.
    Mean removed on the leading samples so an instrument offset does not tilt
    the energy curve."""
    a = a - np.mean(a[:200])
    e = np.cumsum(a ** 2)
    if e[-1] <= 0:
        return np.zeros_like(e)
    return e / e[-1]


def trapz(y, dx):
    """Trapezoidal integral, written out so it does not depend on np.trapz /
    np.trapezoid (the name moved between numpy versions)."""
    if len(y) < 2:
        return 0.0
    return float(np.sum(0.5 * (y[1:] + y[:-1])) * dx)


def window_and_condition(t, d, a, lowpass_hz=None, drift_correct=True):
    """Return (tt, v, info) for the significant-duration-truncated table motion.

    tt starts at 0, v starts and ends at exactly 0. See module docstring."""
    dt = float(np.median(np.diff(t)))

    # velocity from the measured table displacement
    v_full = np.gradient(d, dt)
    if lowpass_hz:
        v_full = etv.lowpass(v_full, dt, lowpass_hz)

    # 95% Arias time on the table acceleration
    ec = arias_cumfrac(a)
    i95 = int(np.searchsorted(ec, 0.95))
    i95 = min(max(i95, 1), len(t) - 1)

    # nearest velocity zero-crossing to the 95% Arias index
    s = np.sign(v_full)
    s[s == 0] = 1.0
    zc = np.where(np.diff(s) != 0)[0]     # v[i], v[i+1] straddle zero
    if len(zc) == 0:
        i_cut = i95
    else:
        i_cut = int(zc[np.argmin(np.abs(zc - i95))])
    i_cut = max(i_cut, 2)

    tt = t[:i_cut + 1] - t[0]
    v = v_full[:i_cut + 1].copy()
    d_win = d[:i_cut + 1]

    # clamp endpoints to rest (they already sit on a zero-crossing / at t=0)
    v[0] = 0.0
    v[-1] = 0.0

    T_win = float(tt[-1]) if len(tt) > 1 else dt
    drift_raw = trapz(v, dt)              # net table displacement of raw window

    if drift_correct and T_win > 0:
        # zero-endpoint half-sine baseline that integrates to drift_raw:
        #   integral_0^T sin(pi t / T) dt = 2 T / pi  ->  k = drift_raw*pi/(2T)
        k = drift_raw * math.pi / (2.0 * T_win)
        v = v - k * np.sin(math.pi * tt / T_win)
        v[0] = 0.0
        v[-1] = 0.0
        bump_frac = abs(k) / (np.max(np.abs(v)) + 1e-30)
    else:
        bump_frac = 0.0
    drift_final = trapz(v, dt)

    d_recon = np.concatenate(([0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)))
    pp = lambda x: float(np.max(x) - np.min(x))

    info = {
        "dt": dt,
        "n": len(v),
        "dur": T_win,
        "i95": i95, "t95": float(t[i95] - t[0]),
        "i_cut": i_cut, "t_cut": float(t[i_cut] - t[0]),
        "v_peak": float(np.max(np.abs(v))),
        "pgd_amp_mm": pp(d_win) * 1000.0 / 2.0,
        "d_win_pp_mm": pp(d_win) * 1000.0,
        "d_full_pp_mm": pp(d) * 1000.0,
        "d_recon_pp_mm": pp(d_recon) * 1000.0,
        "fidelity": pp(d_recon) / max(pp(d_win), 1e-12),
        "drift_raw_mm": drift_raw * 1000.0,
        "drift_final_mm": drift_final * 1000.0,
        "bump_frac": bump_frac,
    }
    return tt, v, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", type=int, default=9, choices=(9, 12))
    ap.add_argument("--exp-dir", default=None)
    ap.add_argument("--out-dir", default="record_tables_full",
                    help="folder for the tables + manifest "
                         "(default: record_tables_full)")
    ap.add_argument("--runs", type=int, nargs=2, default=None,
                    metavar=("LO", "HI"))
    ap.add_argument("--lowpass", type=float, default=None,
                    help="low-pass corner in Hz applied to the velocity")
    ap.add_argument("--no-drift-correct", action="store_true",
                    help="apply the raw truncated motion, keeping the base "
                         "offset at the 95%% Arias cut")
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
    drift_correct = not args.no_drift_correct

    manifest = []
    for rn in range(lo, hi + 1):
        src = exp_dir / "Test{}Run{}.xlsx".format(args.test, rn)
        if not src.is_file():
            print("  missing, skipped: {}".format(src.name))
            continue
        t, d, a = etv.read_run(src)
        tt, v, info = window_and_condition(
            t, d, a, lowpass_hz=args.lowpass, drift_correct=drift_correct)

        name = "expfull_T{}_run{:02d}".format(args.test, rn)
        table_file = out_dir / "vel_expfull_T{}_run{:02d}.txt".format(
            args.test, rn)
        etv.write_table(table_file, name, tt, v)

        manifest.append({
            "run": rn,
            "table_file": table_file.name,
            "dur_s": round(info["dur"], 4),
            "n_samples": info["n"],
            "dt": round(info["dt"], 6),
            "pgd_amp_mm": round(info["pgd_amp_mm"], 4),
            "pgv_mps": round(info["v_peak"], 5),
            "net_drift_mm": round(info["drift_final_mm"], 4),
            "t95_s": round(info["t95"], 3),
            "t_cut_s": round(info["t_cut"], 3),
        })
        print("  run {:2d}: 0 -> {:.2f} s ({} pts)  PGD {:.2f} mm  PGV {:.3f} m/s"
              "  drift {:+.2f}->{:+.3f} mm  ->  {}".format(
                  rn, info["dur"], info["n"], info["pgd_amp_mm"],
                  info["v_peak"], info["drift_raw_mm"], info["drift_final_mm"],
                  table_file.name))

    if not manifest:
        raise SystemExit("no records found under {}".format(exp_dir))

    man_path = out_dir / "record_tables_full_manifest_US{}.csv".format(us)
    cols = ["run", "table_file", "dur_s", "n_samples", "dt",
            "pgd_amp_mm", "pgv_mps", "net_drift_mm", "t95_s", "t_cut_s"]
    with open(str(man_path), "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for m in manifest:
            w.writerow([m[c] for c in cols])

    total = sum(m["dur_s"] for m in manifest)
    print("\n-> {} tables + {}".format(len(manifest), man_path.name))
    print("   drift correction: {}".format(
        "off (raw offsets kept)" if args.no_drift_correct else "on (half-sine)"))
    print("   total model time to run them all: {:.0f} s "
          "({:.0f} runs x {:.1f} s mean)".format(
              total, len(manifest), total / len(manifest)))


if __name__ == "__main__":
    main()