# check_overload.py -- the tests' spring precompression (channels 6 'Overload left', 7 'Overload right', kN) per run:
# value before shaking, peak / minimum during the run, value after, and its relation to the wall displacement.
import numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
rows = []
for test, lab, nmax in ((9, "US1", 24), (12, "US2", 25)):
    for n in range(1, nmax + 1):
        df = pd.read_excel("EXP_DATA/Test{}Run{}.xlsx".format(test, n)); df.columns = df.columns.map(str).str.strip()
        t = df["Time"].values.astype(float); fs = 1 / np.median(np.diff(t))
        L, R_ = df["6"].values.astype(float), df["7"].values.astype(float)
        S = L + R_
        b, a = butter(4, 10.0 / (0.5 * fs)); Sf = filtfilt(b, a, S)
        d = 0.5 * (df["3"].values + df["4"].values).astype(float) * 1e3; d = filtfilt(b, a, d - d[0])
        n0 = int(1.0 * fs)
        s0, s1 = float(np.median(S[:n0])), float(np.median(S[-n0:]))
        i = int(np.argmax(np.abs(d)))
        k = np.polyfit(np.abs(d), Sf, 1)[0] if np.ptp(d) > 0.5 else np.nan
        rows.append((lab, n, float(np.median(L[:n0])), float(np.median(R_[:n0])), s0, float(Sf.max()), float(Sf.min()), s1,
                     float(Sf[i]), float(np.abs(d).max()), k))
print("test run | left0 right0 | sum: before  max   min   after | sum at peak |d| | peak|d| mm | dSum/d|d| kN/mm")
for r in rows:
    print("{} {:2d} | {:6.2f} {:6.2f} | {:6.2f} {:6.2f} {:6.2f} {:6.2f} | {:6.2f} | {:6.2f} | {:6.3f}".format(*r))
