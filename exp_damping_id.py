"""exp_damping_id.py -- small-amplitude damping ratio of US-1 / US-2 from the raw test channels (read only).

Two independent estimates, both on the RELATIVE acceleration of the top-quarter sensor (a17 - a12, absolute channels):
  A. ring-down after the table stops (t_table_true_stop_s from exp_envelope_peaks/ringdown_summary.csv):
     band-pass 4-30 Hz, positive peaks above 4x the tail noise, exponential fit ln A(t) -> zeta = -slope T / (2 pi)
     (per-cycle log decrements are also written; they scatter strongly because the free motion is micro-metres)
  B. transfer function table -> wall (H1 = S_xy / S_xx, 8 s Welch segments, 0.125 Hz resolution), first mode peak
     in 5-20 Hz, half-power bandwidth zeta = (f2 - f1) / (2 fn); coherence at the peak is reported.
     Trusted rows: peak at 8.5-12.5 Hz (the wall mode) and coherence >= 0.85.
Caveats: both describe the closed-joint (elastic) wall at very small amplitude; neither measures the energy lost per
rocking impact at large amplitude, which the DEM models through its contacts, not through Maxwell damping.
The 0.125 Hz resolution limits B to zeta >~ 0.6 %.

USAGE   python exp_damping_id.py <raw_dir> <ringdown_summary.csv> <out_dir>
        python exp_damping_id.py EXP_DATA exp_envelope_peaks/ringdown_summary.csv exp_damping
"""
import sys, os, glob, re
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks, csd, welch, coherence

TESTS = (("US1", "Test9"), ("US2", "Test12"))


def load(p):
    df = pd.read_excel(p); df.columns = df.columns.map(str).str.strip()
    t = pd.to_numeric(df["Time"], errors="coerce").values
    a12 = pd.to_numeric(df["12"], errors="coerce").values; a17 = pd.to_numeric(df["17"], errors="coerce").values
    m = np.isfinite(t) & np.isfinite(a12) & np.isfinite(a17)
    return t[m], a12[m] * 9.81, a17[m] * 9.81


def ringdown(t, x, t0):
    fs = 1 / np.median(np.diff(t))
    b, a = butter(4, [4 / (fs / 2), 30 / (fs / 2)], "band")
    y = filtfilt(b, a, x - np.mean(x)); w = t >= t0 + 0.05; ys, ts = y[w], t[w]
    if len(ys) < fs:
        return None, []
    noise = np.std(ys[-int(min(1.0, len(ys) / fs / 4) * fs):])
    pk, _ = find_peaks(ys, distance=int(fs / 25)); pk = pk[ys[pk] > 4 * noise]
    if len(pk) < 4:
        return None, []
    A, T = ys[pk], np.diff(ts[pk]); cyc = []
    k = len(pk)
    for i in range(len(pk) - 1):
        if T[i] > 0.2:
            k = i + 1; break
        d = np.log(A[i] / A[i + 1])
        if d > 0:
            cyc.append(dict(zeta=d / np.sqrt(4 * np.pi ** 2 + d ** 2), amp_mm=0.5 * (A[i] + A[i + 1]) / (2 * np.pi / T[i]) ** 2 * 1e3, T=T[i]))
    if k < 4:
        return None, cyc
    sl = np.polyfit(ts[pk[:k]] - ts[pk[0]], np.log(A[:k]), 1)[0]; Tm = float(np.median(T[:k - 1]))
    return dict(n_peaks=k, T_ring=Tm, zeta_ring_fit=-sl * Tm / (2 * np.pi), amp0_mm=float(A[0] / (2 * np.pi / Tm) ** 2 * 1e3)), cyc


def frf(t, x, y):
    fs = 1 / np.median(np.diff(t)); nper = int(fs * 8)
    x = x - np.mean(x); y = y - np.mean(y)
    f, Pxy = csd(x, y, fs, nperseg=nper); _, Pxx = welch(x, fs, nperseg=nper); _, C = coherence(x, y, fs, nperseg=nper)
    H = np.abs(Pxy / Pxx); band = (f > 5) & (f < 20); i = int(np.argmax(np.where(band, H, 0))); fn, Hp = f[i], H[i]
    half = Hp / np.sqrt(2); j, k = i, i
    while j > 0 and H[j] > half: j -= 1
    while k < len(H) - 1 and H[k] > half: k += 1
    f1 = np.interp(half, [H[j], H[j + 1]], [f[j], f[j + 1]]); f2 = np.interp(half, [H[k], H[k - 1]], [f[k], f[k - 1]])
    return dict(f_peak_Hz=fn, T_frf=1 / fn, zeta_halfpower=(f2 - f1) / (2 * fn), coherence=C[i])


def main():
    raw, rsum, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    rs = pd.read_csv(rsum); rows, cycles = [], []
    for lab, test in TESTS:
        files = sorted(glob.glob(os.path.join(raw, test + "Run*.xlsx")), key=lambda p: int(re.search(r"Run(\d+)", p).group(1)))
        for p in files:
            if "processed" in p.lower():
                continue
            n = int(re.search(r"Run(\d+)", p).group(1)); t, a12, a17 = load(p); rel = a17 - a12
            r = dict(test=lab, run=n); r.update(frf(t, a12, rel))
            q = rs[(rs.specimen == lab) & (rs.run == n)]
            if not q.empty:
                rd, cyc = ringdown(t, rel, float(q.t_table_true_stop_s.iloc[0]))
                if rd:
                    r.update(rd)
                cycles += [dict(test=lab, run=n, **c) for c in cyc]
            rows.append(r); print("{} run {:2d}: FRF {:.2f} Hz zeta {:.3f} (coh {:.2f}) | ring-down fit {}".format(
                lab, n, r["f_peak_Hz"], r["zeta_halfpower"], r["coherence"],
                "{:.3f}".format(r["zeta_ring_fit"]) if "zeta_ring_fit" in r else "-"))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(out, "damping_runs.csv"), index=False)
    pd.DataFrame(cycles).to_csv(os.path.join(out, "ringdown_cycles.csv"), index=False)
    ok = df[(df.f_peak_Hz.between(8.5, 12.5)) & (df.coherence >= 0.85)]
    print("\ntrusted FRF rows (8.5-12.5 Hz, coherence >= 0.85):")
    for lab in ("US1", "US2"):
        s = ok[ok.test == lab].zeta_halfpower
        rf = df[df.test == lab].zeta_ring_fit.dropna()
        print("  {}: half-power zeta median {:.4f} (range {:.4f}-{:.4f}, {} runs) | ring-down fit median {:.4f} ({} runs)".format(
            lab, s.median(), s.min(), s.max(), len(s), rf.median(), len(rf)))
    print("->", out)


if __name__ == "__main__":
    main()
