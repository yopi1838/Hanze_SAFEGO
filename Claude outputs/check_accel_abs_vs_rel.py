"""check_accel_abs_vs_rel.py -- are the wall accelerometers (ch 15/16/17) ABSOLUTE or RELATIVE to the table (ch 12)?

The processed-file script assumed 'relative' and added ch12. An accelerometer on the wall measures absolute
acceleration unless the channel was derived. The test: below ~3 Hz the wall moves with the table as a rigid
body, so in a 0.5-3 Hz band
    absolute channel  ->  ch_j ~= ch12        (ratio ~1, correlation ~+1)
    relative channel  ->  ch_j ~= 0           (ratio << 1)
Then the run's peak base shear is recomputed with the group's masses under BOTH assumptions, and with the
old calibrated masses for reference, so the number that went into the loops can be traced.

Usage:  python check_accel_abs_vs_rel.py Test9Run24.xlsx [more raw run files ...]
        python check_accel_abs_vs_rel.py EXPRAW/Test9Run24.xlsx EXPRAW/Test12Run24.xlsx
Raw file layout (as make_processed_globalzero.py reads it): row 1 = [None, 'Time', 1, 2, ..., 17]; data below.
Outputs: per file, a table on screen and <stem>_abs_vs_rel.png next to the file.
"""
import sys, os
import numpy as np
from scipy.signal import butter, filtfilt, coherence
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G0 = 9.81
BAND = (0.5, 3.0)                 # rigid-body band, Hz
M_GROUP = {15: 455.03, 16: 670.07, 17: 606.47}   # group's sheet: wall 1/4-1/2-1/4 at 18.8 kN/m3 + slabs + beam, no factor
M_OLD   = {15: 946.679, 16: 743.092, 17: 936.500} # processed files: midpoint split at 1619 x 2.317 (calibrated to 29.5 kN)

def load_raw(path):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True); hdr = next(it)
    col = {}
    for j, h in enumerate(hdr):
        if h is None: continue
        if isinstance(h, str) and h.strip().lower() == "time": col["time"] = j
        else:
            try: col[int(float(h))] = j
            except (TypeError, ValueError): pass
    rows = [r for r in it if r[col["time"]] is not None]
    A = np.array([[r[j] if r[j] is not None else np.nan for j in range(len(hdr))] for r in rows], float)
    return A[:, col["time"]], {k: A[:, j] for k, j in col.items() if k != "time"}

def bandpass(x, fs, lo, hi):
    b, a = butter(3, [lo / (fs / 2), hi / (fs / 2)], "band"); return filtfilt(b, a, x)

def analyse(path):
    t, ch = load_raw(path); fs = 1.0 / float(np.median(np.diff(t)))
    a12 = ch[12]
    # strong-motion window: where the table acceleration exceeds 5% of its peak
    m = np.abs(a12) > 0.05 * np.abs(a12).max(); i0, i1 = np.where(m)[0][[0, -1]]; s = slice(i0, i1)
    print(f"\n{os.path.basename(path)}: fs {fs:.0f} Hz, table PGA {np.abs(a12).max():.3f} g, window {t[i0]:.1f}-{t[i1]:.1f} s")
    print(f"  rigid-body band {BAND[0]}-{BAND[1]} Hz:   {'ch':>4} {'rms_j/rms_12':>12} {'slope j~12':>10} {'corr':>6} {'coh':>6}   verdict")
    x12 = bandpass(a12, fs, *BAND)[s]
    verdicts = {}
    for j in (15, 16, 17):
        xj = bandpass(ch[j], fs, *BAND)[s]
        ratio = np.sqrt(np.mean(xj**2)) / np.sqrt(np.mean(x12**2))
        slope = float(np.polyfit(x12, xj, 1)[0]); corr = float(np.corrcoef(x12, xj)[0, 1])
        f, C = coherence(a12[s], ch[j][s], fs=fs, nperseg=min(2048, (i1 - i0) // 4)); coh = float(np.mean(C[(f >= BAND[0]) & (f <= BAND[1])]))
        v = "ABSOLUTE (moves with the table)" if (ratio > 0.6 and corr > 0.7) else ("RELATIVE (table removed)" if ratio < 0.35 else "ambiguous")
        verdicts[j] = v
        print(f"  {'':34s}{j:4d} {ratio:12.2f} {slope:10.2f} {corr:6.2f} {coh:6.2f}   {v}")
    # peak base shear under both assumptions
    print(f"  peak base shear (kN), this run:")
    print(f"    {'masses':>34s} {'ch15-17 absolute as recorded':>30s} {'ch15-17 relative (+ch12)':>26s}")
    for lab, M in (("group's sheet (1732 kg, no factor)", M_GROUP), ("old processed files (2626 kg)", M_OLD)):
        V_abs = sum(M[j] * ch[j] for j in (15, 16, 17)) * G0
        V_rel = sum(M[j] * (ch[j] + a12) for j in (15, 16, 17)) * G0
        print(f"    {lab:>34s} {np.nanmax(np.abs(V_abs))/1e3:30.2f} {np.nanmax(np.abs(V_rel))/1e3:26.2f}")
    # figure: low-pass traces and scatter
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[2, 1]))
    tt = t[s]; ax[0].plot(tt, x12, color="#18212b", lw=1.4, label="ch12 table")
    for j, c in ((15, "#4a7c59"), (16, "#2b5f8a"), (17, "#b4491a")):
        ax[0].plot(tt, bandpass(ch[j], fs, *BAND)[s], color=c, lw=0.9, alpha=0.9, label=f"ch{j} ({verdicts[j].split(' ')[0].lower()})")
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel(f"acceleration, {BAND[0]}-{BAND[1]} Hz band (g)"); ax[0].legend(fontsize=8, ncol=4, frameon=False)
    ax[0].set_title(f"{os.path.basename(path)}: do the wall channels follow the table at rigid-body frequencies?", loc="left", fontsize=10, fontweight="bold")
    x17 = bandpass(ch[17], fs, *BAND)[s]
    ax[1].plot(x12, x17, ".", ms=2, color="#b4491a", alpha=0.5); lim = np.abs(x12).max()
    ax[1].plot([-lim, lim], [-lim, lim], "k:", lw=1, label="1:1 (absolute)"); ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel("ch12 table (g)"); ax[1].set_ylabel("ch17 top quarter (g)"); ax[1].legend(fontsize=8, frameon=False); ax[1].set_aspect("equal")
    for a in ax:
        a.grid(alpha=0.25)
        for sp in ("top", "right"): a.spines[sp].set_visible(False)
    fig.tight_layout(); out = os.path.splitext(path)[0] + "_abs_vs_rel.png"; fig.savefig(out, dpi=150); plt.close(fig); print(f"  -> {out}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    for p in sys.argv[1:]: analyse(p)
