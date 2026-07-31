# -*- coding: ascii -*-
"""
Recompute model tilt using the EXPERIMENTAL definitions, from channel CSVs
that have already been exported. Works on existing run sets -- no 3DEC
re-run needed.

Why: the model and the experiment were reporting two different quantities
under the same name "tilt".

  chord tilt      atan(d_top / H).  What the displacement transducers give
                  (exp_Test9_tilt.csv), and what the model's tilt_full_wall
                  already is. These two ARE comparable.

  inclinometer    exp_US1_tiltmeter.csv. A gravity-referenced MEMS device
                  bonded to one masonry unit. Its settled value is that
                  unit's rotation -- NOT the wall's chord tilt. Its peak
                  value is not a rotation at all: during shaking the device
                  reads horizontal acceleration as apparent tilt.

This script produces both, so each can be compared against the right
experimental series:

  1. chord tilt, experiment definition
     d_top = mean(Channel_3, Channel_4)   (the experiment averages the two
     top quarter points; the model's FISH used the centreline instead)
     d_table subtracted, matching the frame-mounted transducers.

  2. emulated inclinometer
     theta_apparent = atan( a_z/g + tan(theta_true) )
     a_z is differentiated from the model's velocity channel at the sensor
     height. theta_true needs a genuine local-rotation channel from the
     model (tilt_local_incl in instrument_tilt_v2.dat); without it, only
     the acceleration artefact is reproduced, which is still enough to show
     that the measured PEAKS are artefact-dominated.

Usage
    python tilt_experiment_match.py [SIM_DIR] [--test {9,12}] [--out CSV]

SIM_DIR defaults to safego_paths.sim_dir(). Output goes to
<SIM>/postproc/sim_tilt_expmatch.csv unless --out is given.
"""
import argparse, csv, glob, math, os
import numpy as np

try:
    import safego_paths as sp
    _HAVE_SP = True
except ImportError:
    _HAVE_SP = False

G = 9.80665
TAIL_FRAC = 0.05                 # matches postprocess_stratC.py
RUN_RE_PREFIX = "Run"

# Sensor elevations differ by specimen (see exp_tilt_from_raw.py).
TEST_CFG = {9:  {"bot_z": 0.66, "label": "US-1 (Test 9)"},
            12: {"bot_z": 0.60, "label": "US-2 (Test 12)"}}

H_FULL = 2.06                    # table -> top sensor chord


def read_hist(path):
    """3DEC history export: 2 header lines, then whitespace-separated
    time / value columns."""
    d = np.genfromtxt(path, skip_header=2)
    if d.ndim != 2 or d.shape[1] < 2:
        raise ValueError("unexpected shape in {}".format(path))
    m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    return d[m, 0], d[m, 1]


def find_channel(folder, pattern):
    hits = glob.glob(os.path.join(folder, "*" + pattern + "*"))
    return hits[0] if hits else None


def metrics(v, ref=None):
    """peak = max |v - v(t0)| within the run; resid = mean of the last
    TAIL_FRAC, relative to ref (or to this run's start if ref is None).
    Same convention as postprocess_stratC.py."""
    if len(v) == 0:
        return float("nan"), float("nan")
    base = v[0] if ref is None else ref
    peak = float(np.max(np.abs(v - v[0])))
    resid = float(np.mean(v[int((1.0 - TAIL_FRAC) * len(v)):]) - base)
    return peak, resid


def _lowpass(x, dt, fc):
    """Zero-phase low-pass. Mimics the finite bandwidth of a real MEMS
    inclinometer, without which the differentiated velocity is dominated by
    rocking-impact spikes that no physical device would register."""
    if fc is None or fc <= 0 or dt <= 0:
        return x
    fs = 1.0 / dt
    if fc >= 0.5 * fs:                           # already below Nyquist
        return x
    try:
        from scipy.signal import butter, filtfilt
        b, a = butter(3, fc / (0.5 * fs), btype="low")
        return filtfilt(b, a, x)
    except Exception:
        # boxcar fallback if scipy is unavailable
        w = max(1, int(round(fs / fc)))
        k = np.ones(w) / w
        return np.convolve(x, k, mode="same")


def emulate_inclinometer(t, vel_z, theta_true_deg=None, fc=25.0):
    """Reproduce what a gravity-referenced MEMS inclinometer would report.

    The device senses the direction of apparent gravity. Under horizontal
    acceleration a_z it reads

        theta_apparent = atan( (a_z + g*sin(theta)) / (g*cos(theta)) )
                       = atan( a_z/g + tan(theta) )

    so during shaking the reading is dominated by a_z, and only settles to
    the true rotation once the table stops. This is why a measured peak of
    9.9 deg at run 1 coexists with 0.33 mm of wall movement, and why run 3
    reads 26.8 deg while the wall moves 0.018 mm.

    CAVEAT: the peak value depends strongly on the assumed device bandwidth
    `fc`, which is not documented for this instrument. Treat the emulated
    peak as an order-of-magnitude demonstration that the reading is
    acceleration-dominated, not as a calibrated prediction. The settled
    value is much less sensitive to fc and is the one worth comparing.
    """
    if len(t) < 3:
        return np.zeros_like(vel_z)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    a_z = np.gradient(vel_z, t)                  # m/s^2
    a_z = _lowpass(a_z, dt, fc)
    if theta_true_deg is None:
        tan_true = 0.0
    else:
        tan_true = np.tan(np.radians(theta_true_deg))
    return np.degrees(np.arctan(a_z / G + tan_true))


def process_run(folder, bot_z, incl_fc=25.0):
    """Return a dict of experiment-matched metrics for one run folder."""
    def ch(pat):
        p = find_channel(folder, pat)
        return read_hist(p) if p else (None, None)

    t3, c3 = ch("Channel_3_")
    t4, c4 = ch("Channel_4_")
    t5, c5 = ch("Channel_5_DispTable")
    t1, c1 = ch("Channel_1_DispBot")
    t2, c2 = ch("Channel_2_DispMid")
    t19, c19 = ch("Channel_19")
    t17, v17 = ch("Channel_17_AccTop")           # velocity, despite the name
    tloc, cloc = ch("tilt_local_incl")           # only if instrument_tilt_v2 was used

    if c3 is None or c4 is None or c5 is None:
        return None

    n = min(len(x) for x in (c3, c4, c5) if x is not None)
    c3, c4, c5 = c3[:n], c4[:n], c5[:n]
    t = t3[:n]

    d_top_exp = 0.5 * (c3 + c4)                  # experiment definition
    rel_exp = (d_top_exp - c5) * 1000.0          # mm
    tilt_full_exp = np.degrees(np.arctan((d_top_exp - c5) / H_FULL))

    out = {}
    out["peak_rel_disp_expdef_mm"], out["resid_rel_disp_expdef_mm"] = metrics(rel_exp)
    out["peak_tilt_full_expdef_deg"], out["resid_tilt_full_expdef_deg"] = metrics(tilt_full_exp)

    # the as-built centreline definition, for the difference it makes
    if c19 is not None:
        m = min(n, len(c19))
        tilt_full_centre = np.degrees(np.arctan((c19[:m] - c5[:m]) / H_FULL))
        out["peak_tilt_full_centreline_deg"], _ = metrics(tilt_full_centre)
        out["centreline_bulge_max_mm"] = float(np.max(np.abs(c19[:m] - d_top_exp[:m])) * 1000.0)
    # twist / in-plane asymmetry between the two top quarter points
    out["quarterpoint_spread_max_mm"] = float(np.max(np.abs(c3 - c4)) * 1000.0)

    # segment tilts with the specimen's own sensor elevation
    if c1 is not None and c2 is not None:
        m = min(n, len(c1), len(c2))
        tb = np.degrees(np.arctan((c1[:m] - c5[:m]) / bot_z))
        tl = np.degrees(np.arctan((c2[:m] - c1[:m]) / (1.26 - bot_z)))
        out["peak_tilt_bot_seg_deg"], out["resid_tilt_bot_seg_deg"] = metrics(tb)
        out["peak_tilt_low_seg_deg"], out["resid_tilt_low_seg_deg"] = metrics(tl)

    # emulated inclinometer
    if v17 is not None:
        m = min(len(t17), len(v17))
        theta_true = None
        if cloc is not None:
            k = min(m, len(cloc))
            theta_true = cloc[:k]
            m = k
        incl = emulate_inclinometer(t17[:m], v17[:m], theta_true, fc=incl_fc)
        out["incl_emulated_peak_deg"] = float(np.max(np.abs(incl)))
        out["incl_emulated_resid_deg"] = float(np.mean(incl[int((1 - TAIL_FRAC) * len(incl)):]))
        out["incl_true_rotation_available"] = "Y" if cloc is not None else "N"
    return out


def main():
    ap = argparse.ArgumentParser()
    default_sim = str(sp.sim_dir()) if _HAVE_SP else "stratC_results_NODAMP_v6_NEW"
    ap.add_argument("sim_dir", nargs="?", default=default_sim)
    ap.add_argument("--test", type=int, default=9, choices=(9, 12),
                    help="which specimen's sensor layout to match (default 9)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--incl-fc", dest="incl_fc", type=float, default=25.0,
                    help="assumed inclinometer bandwidth in Hz for the emulation "
                         "(default 25; set 0 to disable filtering)")
    args = ap.parse_args()

    sim = args.sim_dir
    bot_z = TEST_CFG[args.test]["bot_z"]
    out_csv = args.out
    if out_csv is None:
        pp = os.path.join(sim, "postproc")
        if not os.path.isdir(pp):
            os.makedirs(pp)
        out_csv = os.path.join(pp, "sim_tilt_expmatch.csv")

    folders = sorted(f for f in glob.glob(os.path.join(sim, RUN_RE_PREFIX + "*"))
                     if os.path.isdir(f))
    if not folders:
        raise SystemExit("no Run* folders under {}".format(sim))

    print("matching sensor layout of {} (bottom sensor at {:.2f} m)".format(
        TEST_CFG[args.test]["label"], bot_z))

    rows, fields = [], ["run"]
    for f in folders:
        try:
            run_no = int(os.path.basename(f)[3:5])
        except ValueError:
            continue
        r = process_run(f, bot_z, args.incl_fc)
        if r is None:
            print("  run {:02d}: channels missing, skipped".format(run_no))
            continue
        r["run"] = run_no
        rows.append(r)
        for k in r:
            if k not in fields:
                fields.append(k)
        print("  run {:02d}: peak_rel {:.3f} mm | tilt_full {:.4f} deg | bulge {:.3f} mm".format(
            run_no, r.get("peak_rel_disp_expdef_mm", float('nan')),
            r.get("peak_tilt_full_expdef_deg", float('nan')),
            r.get("centreline_bulge_max_mm", float('nan'))))

    with open(out_csv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["run"]):
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v)
                        for k, v in r.items()})
    print("\n-> {} ({} runs)".format(out_csv, len(rows)))
    print("\nCompare:")
    print("  *_expdef_*            vs  exp_Test{9,12}_tilt.csv      (chord tilt)")
    print("  incl_emulated_*       vs  exp_US{1,2}_tiltmeter.csv    (inclinometer)")
    print("  never compare chord tilt against the inclinometer columns.")


if __name__ == "__main__":
    main()
