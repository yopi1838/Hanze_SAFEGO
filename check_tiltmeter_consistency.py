# -*- coding: ascii -*-
"""
Cross-check the tiltmeter's tilt channel against its own accelerometers.

The logger file (US1_Tilt_values.csv / US2_Tilt_values.csv) records, every 5 s,
a two-axis tilt reading AND two independent 3-axis accelerometer packs. At rest
the accelerometer vector IS the gravity direction, so

    theta_static = atan2(a_y, a_z)

is a direct, independent measurement of the same rotation the tilt channel
reports. This script compares them.

It exists because on both specimens the two disagree by a factor of 20-30 --
see MODELLING_DISCREPANCIES.md section 1. Since the permanent-rotation curve
derived from the tilt channel is the headline experimental evidence for
ratcheting, that disagreement has to be resolved before the model comparison
means anything.

Usage
    python check_tiltmeter_consistency.py [TILT_CSV ...]

With no arguments it checks EXP_DATA/US1_Tilt_values.csv and
EXP_DATA/US2_Tilt_values.csv.
"""
import csv, math, sys, os
import datetime as dt
import statistics as st

try:
    import safego_paths as sp
    _HAVE_SP = True
except ImportError:
    _HAVE_SP = False

ACC_THR_FACTOR = 5.0        # same run-detection rule as process_tiltmeter.py
RUN_GAP_S = 30.0
SETTLE_N = 6                # logger windows used for the post-run settled value


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter=';'):
            d = {"time": r["time"]}
            for k, v in r.items():
                if k == "time":
                    continue
                try:
                    d[k] = float(v.replace(",", "."))
                except (ValueError, AttributeError):
                    d[k] = float("nan")
            rows.append(d)
    if not rows:
        raise SystemExit("no rows in {}".format(path))
    t0 = dt.datetime.strptime(rows[0]["time"], "%d-%m-%Y_%H:%M:%S")
    for r in rows:
        r["s"] = (dt.datetime.strptime(r["time"], "%d-%m-%Y_%H:%M:%S") - t0).total_seconds()
    return rows


def swing(r):
    return max(r["sensor1 amaxy1"] - r["sensor1 aminy1"],
               r["sensor1 amaxx1"] - r["sensor1 aminx1"])


def acc_angle(r, pack):
    """Static tilt implied by the gravity vector, in degrees."""
    return math.degrees(math.atan2(r["sensor1 avaly" + pack],
                                   r["sensor1 avalz" + pack]))


def detect_runs(rows):
    sw = [swing(r) for r in rows]
    base = st.median(sw)
    act = [i for i, v in enumerate(sw) if v > ACC_THR_FACTOR * base]
    if not act:
        return [], base
    runs, start, prev = [], act[0], act[0]
    for i in act[1:]:
        if rows[i]["s"] - rows[prev]["s"] > RUN_GAP_S:
            runs.append((start, prev))
            start = i
        prev = i
    runs.append((start, prev))
    return runs, base


def check(path):
    rows = load(path)
    runs, base = detect_runs(rows)
    print("=" * 78)
    print(os.path.basename(path))
    print("=" * 78)
    print("  {} logger windows over {:.0f} min, {} run windows detected".format(
        len(rows), rows[-1]["s"] / 60.0, len(runs)))

    # sanity: is the accelerometer really seeing gravity at rest?
    pre = rows[:runs[0][0]] if runs else rows[:50]
    norms = [math.sqrt(r["sensor1 avalx1"] ** 2 + r["sensor1 avaly1"] ** 2
                       + r["sensor1 avalz1"] ** 2) for r in pre]
    scale = st.mean(norms)
    print("  resting |a| = {:.0f} counts (sd {:.0f})  -> {:.0f} counts/g".format(
        scale, st.pstdev(norms), scale))
    print("  clipping: accel rails at 32760 counts ({:.2f} g); "
          "tilt channel rails near +/-25.75 deg".format(32760.0 / scale))

    b_t = st.median([r["sensor1 tvaly"] for r in pre])
    b1 = st.median([acc_angle(r, "1") for r in pre])
    b2 = st.median([acc_angle(r, "2") for r in pre])

    print("\n  {:<4} {:>14} {:>13} {:>13} {:>9}".format(
        "run", "tilt ch (deg)", "accel1 (deg)", "accel2 (deg)", "ratio"))
    for k, (a, b) in enumerate(runs, 1):
        seg = rows[b + 1: b + 1 + SETTLE_N]
        if not seg:
            continue
        t = st.median([r["sensor1 tvaly"] for r in seg]) / 1000.0 - b_t / 1000.0
        a1 = st.median([acc_angle(r, "1") for r in seg]) - b1
        a2 = st.median([acc_angle(r, "2") for r in seg]) - b2
        ratio = t / a1 if abs(a1) > 0.005 else float("nan")
        print("  {:<4d} {:>14.3f} {:>13.3f} {:>13.3f} {:>9.1f}".format(k, t, a1, a2, ratio))

    # final value, measured well after all shaking has stopped
    if runs:
        tail = [r for r in rows if r["s"] > rows[runs[-1][1]]["s"] + 60] or rows[-SETTLE_N:]
    else:
        tail = rows[-SETTLE_N:]
    t = st.median([r["sensor1 tvaly"] for r in tail]) / 1000.0 - b_t / 1000.0
    a1 = st.median([acc_angle(r, "1") for r in tail]) - b1
    a2 = st.median([acc_angle(r, "2") for r in tail]) - b2
    sd_t = st.pstdev([r["sensor1 tvaly"] for r in tail])
    print("\n  FINAL permanent rotation, at full rest ({} windows, tilt sd {:.0f} mdeg):".format(
        len(tail), sd_t))
    print("     tilt channel        : {:+8.3f} deg   <- the published curve".format(t))
    print("     accelerometer pack 1: {:+8.3f} deg".format(a1))
    print("     accelerometer pack 2: {:+8.3f} deg".format(a2))
    if abs(a1) > 1e-4:
        print("     DISAGREEMENT        : {:8.1f} x".format(t / a1))
    print("\n  The two accelerometer packs are independent chips. Where they agree with")
    print("  each other but not with the tilt channel, the tilt channel is the outlier.")


def main():
    args = sys.argv[1:]
    if not args:
        if _HAVE_SP:
            d = sp.exp_data_dir()
            args = [str(d / "US1_Tilt_values.csv"), str(d / "US2_Tilt_values.csv")]
        else:
            args = ["EXP_DATA/US1_Tilt_values.csv", "EXP_DATA/US2_Tilt_values.csv"]
    for p in args:
        if os.path.isfile(p):
            check(p)
            print()
        else:
            print("skipped (not found): {}".format(p))


if __name__ == "__main__":
    main()
