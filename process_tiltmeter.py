# -*- coding: ascii -*-
"""
Extract per-run tiltmeter values (sensor 18) for US-1 from
US1_Tilt_values.csv (5 s logger; semicolon-separated, decimal comma;
tilt columns in millidegrees). Run windows are detected from the
accelerometer swing channels, exactly as in compare_tilt_robust.py.

Metrics per run:
  peak_abs_deg : max |tilt - local pre-run baseline| within the run window
  residual_deg : settled (median) tilt in the quiet window after the run,
                 relative to the global initial baseline (cumulative lean)

Output: <SIM>/postproc/exp_US1_tiltmeter.csv (name follows the input stem)
Usage:  python process_tiltmeter.py [TILT_CSV]
        TILT_CSV defaults to ./EXP_DATA/US1_Tilt_values.csv; pass
        ./EXP_DATA/US2_Tilt_values.csv for the US-2 specimen.
"""
import sys, csv
import numpy as np
import pandas as pd
from pathlib import Path

import safego_paths as sp

HERE = sp.ROOT
DEF_CSV = str(sp.tilt_csv("US1"))   # US2: pass EXP_DATA/US2_Tilt_values.csv
# output name is derived from the input filename (US1_... -> exp_US1_tiltmeter.csv)

# NOTE: the OOP axis of sensor 18 is the tiltmeter's Y axis. The settled
# cumulative tvaly reproduces the published US-1 tilt-angle curve
# (2.77 deg at Run 21, 6.03 deg at Run 24); tvalx gives values ~6x smaller
# (in-plane). compare_tilt_robust.py's use of tvalx is incorrect.
TILT_VAL, TILT_MIN, TILT_MAX = ('sensor1 tvaly', 'sensor1 tminy', 'sensor1 tmaxy')
ACC_COLS = [('sensor1 amaxx1', 'sensor1 aminx1'),
            ('sensor1 amaxy1', 'sensor1 aminy1')]
ACC_THR_FACTOR = 5.0
RUN_GAP_S = 30.0
PRE_BASELINE_N = 3


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEF_CSV
    stem = Path(csv_path).stem.split("_")[0]  # "US1" / "US2"
    out_csv = sp.postproc_dir() / ("exp_{}_tiltmeter.csv".format(stem))
    df = pd.read_csv(csv_path, sep=';', decimal=',')
    df['datetime'] = pd.to_datetime(df['time'], format='%d-%m-%Y_%H:%M:%S')
    sec = (df['datetime'] - df['datetime'].iloc[0]).dt.total_seconds().values

    sw = None
    for cmax, cmin in ACC_COLS:
        s = (df[cmax] - df[cmin]).values.astype(float)
        sw = s if sw is None else np.maximum(sw, s)
    thr = float(np.median(sw)) * ACC_THR_FACTOR
    active = np.where(sw > thr)[0]
    runs, s0, p = [], active[0], active[0]
    for i in active[1:]:
        if (sec[i] - sec[p]) > RUN_GAP_S:
            runs.append((s0, p)); s0 = i
        p = i
    runs.append((s0, p))
    print("detected {} run windows (threshold {:.0f})".format(len(runs), thr))

    tval = df[TILT_VAL].values.astype(float)
    tmin = df[TILT_MIN].values.astype(float)
    tmax = df[TILT_MAX].values.astype(float)
    g0 = float(np.median(tval[:5]))

    rows = []
    for k, (s, e) in enumerate(runs):
        pre = tval[max(0, s - PRE_BASELINE_N):s]
        lb = float(np.median(pre)) if len(pre) else g0
        pk_pos = (float(np.max(tmax[s:e + 1])) - lb) / 1000.0
        pk_neg = (float(np.min(tmin[s:e + 1])) - lb) / 1000.0
        nxt = runs[k + 1][0] if k + 1 < len(runs) else len(tval)
        settle = tval[e + 1:nxt]
        res = ((float(np.median(settle)) if len(settle)
                else float(tval[min(e, len(tval) - 1)])) - g0) / 1000.0
        rows.append({"run": k + 1,
                     "t_start_s": round(float(sec[s]), 1),
                     "dur_s": round(float(sec[e] - sec[s]), 1),
                     "peak_pos_deg": round(pk_pos, 4),
                     "peak_neg_deg": round(pk_neg, 4),
                     "peak_abs_deg": round(max(abs(pk_pos), abs(pk_neg)), 4),
                     "residual_deg": round(res, 4)})
        print("run %2d: peak_abs=%.3f deg  residual(cum)=%.3f deg" % (
            k + 1, rows[-1]["peak_abs_deg"], rows[-1]["residual_deg"]))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(str(out_csv), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("-> {}".format(out_csv))


if __name__ == "__main__":
    main()
