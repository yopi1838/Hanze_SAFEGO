# -*- coding: utf-8 -*-
"""compare_stratF_pocket.py -- Strategy F at 1.5 % Maxwell: glued pocket vs pK11 pocket, against US-1 and US-2.

Run from the Hanze_SAFEGO root:  python compare_stratF_pocket.py
Reads
  stratF_postproc_pk2_ft0p05_5vs1p5/stratF_metrics.csv   (glued pocket, variant "pk2 f t 0.05, 1.5 %")
  stratF_postproc_pK11_1p5/stratF_metrics.csv            (pK11 pocket; postprocess_stratF_full25.py, same settings)
  stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_pK11_tr/RunNN_*/ Channel_8..11 (joist-wall slip, mm; FISH ch8..ch11_slip)
  EXP_DATA/Test9RunNN.xlsx, Test12RunNN.xlsx channels 8..11 (m)
Slip, per run: peak |s - s_start| and residual s_end - s_start, 10 Hz zero-phase low-pass on both sides, magnitudes
(the sign conventions of the test LVDTs and the model channels are not matched). Slab 1 = channels 8, 9; slab 2 = 10, 11
(largest of the pair).
Output: stratF_pocket_compare/fig_stratF_pocket.png, stratF_pocket_compare.csv. Linear axes.
"""
import csv, glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

sys.path.insert(0, ".")
import postprocess_r2tm_full25 as R

OUT = "stratF_pocket_compare"
GLUED = ("stratF_postproc_pk2_ft0p05_5vs1p5/stratF_metrics.csv", "pk2 f t 0.05, 1.5 %")
PK11 = ("stratF_postproc_pK11_1p5/stratF_metrics.csv", "pK11 pocket, 1.5 %")
GLUED5 = ("stratF_postproc_pk2_ft0p05_5vs1p5/stratF_metrics.csv", "pk2 f t 0.05, 5 %")
PK11_DIR = "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_pK11_tr"
SLIP = ["Channel_8_LeftJoistBelowSlab1", "Channel_9_LeftJoistAboveSlab1",
        "Channel_10_LeftJoistBelowSlab2", "Channel_11_LeftJoistBelowSlab2"]
LP = 10.0
RUNS = range(1, 26)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def metrics(path, variant):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["variant"] == variant:
                rows[int(r["run"])] = r
    return rows


def lowpass(t, v):
    fs = 1.0 / np.median(np.diff(t))
    if fs <= 2.5 * LP:
        return v
    b, a = butter(4, LP / (0.5 * fs))
    return filtfilt(b, a, v)


def slip_stats(t, v):
    v = lowpass(np.asarray(t, float), np.asarray(v, float))
    n0 = max(3, len(v) // 200)
    s = v - np.mean(v[:n0])
    return float(np.max(np.abs(s))), float(abs(np.mean(s[-n0:])))


def model_slip(n):
    rf = R.run_folder(PK11_DIR, n)
    out = []
    for c in SLIP:
        t, v = R.find_hist(rf, c) if rf else (None, None)
        out.append(slip_stats(t, v) if t is not None else (np.nan, np.nan))
    return out


def test_slip(test, n):
    import pandas as pd
    p = os.path.join("EXP_DATA", "{}Run{}.xlsx".format(test, n))
    if not os.path.isfile(p):
        return [(np.nan, np.nan)] * 4
    df = pd.read_excel(p)
    df.columns = [str(c).strip() for c in df.columns]
    t = df["Time"].values.astype(float)
    return [slip_stats(t, df[str(k)].values.astype(float) * 1e3) for k in (8, 9, 10, 11)]


def pair_max(s, which):          # which: 0 peak, 1 residual; slab1 = ch 8, 9; slab2 = ch 10, 11
    return (np.nanmax([s[0][which], s[1][which]]), np.nanmax([s[2][which], s[3][which]]))


def main():
    os.makedirs(OUT, exist_ok=True)
    G, P, G5 = metrics(*GLUED), metrics(*PK11), metrics(*GLUED5)
    rows = []
    for n in RUNS:
        g, p, g5 = G.get(n, {}), P.get(n, {}), G5.get(n, {})
        ms, t1, t2 = model_slip(n), test_slip("Test9", n), test_slip("Test12", n)
        r = dict(run=n, record=p.get("record", g.get("record", "")), scale=p.get("scale", ""),
                 d_glued5=fnum(g5.get("peak_incr_mm")), F_glued5=fnum(g5.get("peak_F_kN")),
                 d_glued=fnum(g.get("peak_incr_mm")), d_pK11=fnum(p.get("peak_incr_mm")),
                 d_US1=fnum(p.get("US1_peak_mm")), d_US2=fnum(p.get("US2_peak_mm")),
                 F_glued=fnum(g.get("peak_F_kN")), F_pK11=fnum(p.get("peak_F_kN")),
                 F_US1=fnum(p.get("US1_peak_F_kN")), F_US2=fnum(p.get("US2_peak_F_kN")))
        for lab, s in (("pK11", ms), ("US1", t1), ("US2", t2)):
            (r["slip1_pk_" + lab], r["slip2_pk_" + lab]) = pair_max(s, 0)
            (r["slip1_res_" + lab], r["slip2_res_" + lab]) = pair_max(s, 1)
        rows.append(r)
    with open(os.path.join(OUT, "stratF_pocket_compare.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})

    x = np.array([r["run"] for r in rows])
    A = lambda k: np.array([r[k] for r in rows], float)
    st = dict(glued5=dict(color="#0F766E", lw=1.4, alpha=0.6, ls="--", marker="s", mfc="white", ms=4,
                          label="DEM glued pocket, 5 %"),
              glued=dict(color="#0F766E", lw=2.2, alpha=0.6, marker="s", ms=4, label="DEM glued pocket, 1.5 %"),
              pK11=dict(color="#B91C1C", lw=2.2, alpha=0.6, marker="D", ms=4, label="DEM pK11 pocket, 1.5 %"),
              US1=dict(color="#111827", lw=1.0, marker="o", mfc="white", ms=5, label="US-1"),
              US2=dict(color="#1D4ED8", lw=1.0, ls="--", marker="^", mfc="white", ms=5, label="US-2"))
    fig, ax = plt.subplots(2, 2, figsize=(15, 9.5))
    for a, key, yl in ((ax[0, 0], "d", "peak |d| at 2.06 m, per run (mm)"), (ax[0, 1], "F", "peak |F|, group formula (kN)")):
        for s in ("glued5", "glued", "pK11", "US1", "US2"):
            a.plot(x, A(key + "_" + s), **st[s])
        a.set_ylabel(yl)
    for a, kind, yl in ((ax[1, 0], "pk", "peak joist-wall slip in the run (mm)"),
                        (ax[1, 1], "res", "residual joist-wall slip left by the run (mm)")):
        for s in ("pK11", "US1", "US2"):
            sty = dict(st[s]); lab = sty.pop("label")
            a.plot(x, A("slip2_{}_{}".format(kind, s)), label=lab + ", slab 2 (ch 10-11)", **sty)
            sty2 = dict(sty); sty2.update(alpha=0.35 if s == "pK11" else 0.45, ls=":")
            a.plot(x, A("slip1_{}_{}".format(kind, s)), label=lab + ", slab 1 (ch 8-9)", **sty2)
        a.set_ylabel(yl)
    for a in ax.flat:
        a.set_xticks(list(RUNS)); a.set_xlim(0.5, 25.5); a.set_ylim(bottom=0)
        a.grid(color="#E8E8E8", lw=0.7); a.set_xlabel("run")
        for xv in (21.5,):
            a.axvline(xv, color="#9CA3AF", lw=0.8)
    ax[0, 0].legend(fontsize=8, loc="upper left"); ax[1, 0].legend(fontsize=7.5, loc="upper left", ncol=1)
    fig.suptitle("Strategy F, pk2 f_t 0.05: glued pocket (5 % and 1.5 % Maxwell) vs pK11 pocket (1.5 %) against US-1 / US-2 "
                 "(runs 1-21 HU12/EC40, 22-25 FR76)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_stratF_pocket.png"), dpi=150)

    def gm(a, b, sel):
        a, b = a[sel], b[sel]; k = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        return float(np.exp(np.mean(np.log(a[k] / b[k])))) if k.any() else np.nan
    for lab, sel in (("all runs", x > 0), ("runs 1-21", x <= 21), ("runs 22-25", x >= 22)):
        print("{:10s} disp model/US-1: glued 5% {:.2f} glued 1.5% {:.2f} pK11 {:.2f} | force model/US-1: glued 5% {:.2f} glued 1.5% {:.2f} pK11 {:.2f}".format(
            lab, gm(A("d_glued5"), A("d_US1"), sel), gm(A("d_glued"), A("d_US1"), sel), gm(A("d_pK11"), A("d_US1"), sel),
            gm(A("F_glued5"), A("F_US1"), sel), gm(A("F_glued"), A("F_US1"), sel), gm(A("F_pK11"), A("F_US1"), sel)))
    for r in rows:
        if r["run"] >= 15:
            print("run {:2d}  d glued {:6.2f} pK11 {:6.2f} US1 {:6.2f} US2 {:6.2f} | F glued {:5.2f} pK11 {:5.2f} US1 {:5.2f} US2 {:5.2f}"
                  " | slip2 pk/res pK11 {:5.2f}/{:5.2f} US1 {:5.2f}/{:5.2f} US2 {:5.2f}/{:5.2f}".format(
                      r["run"], r["d_glued"], r["d_pK11"], r["d_US1"], r["d_US2"], r["F_glued"], r["F_pK11"], r["F_US1"],
                      r["F_US2"], r["slip2_pk_pK11"], r["slip2_res_pK11"], r["slip2_pk_US1"], r["slip2_res_US1"],
                      r["slip2_pk_US2"], r["slip2_res_US2"]))
    print("-> " + OUT)


if __name__ == "__main__":
    main()
