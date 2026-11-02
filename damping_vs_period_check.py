# -*- coding: ascii -*-
"""
damping_vs_period_check.py -- is the Strategy F model's under-response a DAMPING or a PERIOD (stiffness) problem?
(2 Nov 2026; plain python, no 3DEC)

For every run it compares, at the 2.06 m sensor height:
  * the MEASURED relative displacement (0.5 (ch3 + ch4)) and relative acceleration (a17 - a12) of the test;
  * the DEM (Strategy F run folder: rel_disp_top_exp_mm, vz_acc17, Channel_12_AccTable), time-aligned on the
    table acceleration by cross-correlation;
  * linear oscillators driven by the MEASURED table acceleration (0.3-30 Hz):
        T_MODEL  (the model's T1) at 5 % and at 1.5 %   -> what damping alone changes
        T_TEST1  (the test's run-1 period) at 1.5 %     -> what the initial period alone changes
        T_eq, xi_eq of that test run (exp_damping/equivalent_damping_runs_acc.csv, if present)
                                                         -> the best linear stand-in for the test itself
  * the model's free decay in the tail after each record (TAIL_S seconds): xi from the 4-30 Hz relative-acceleration
    peaks (log-decrement fit), next to the test's ring-down after its table stops (exp_damping/damping_runs.csv).
READING: if the 5 % and 1.5 % oscillators give nearly the same peaks, damping is not what limits the response;
if the T_TEST1 / T_eq oscillators reach the test's peaks and the T_MODEL ones do not, the gap is period (stiffness).

USAGE
  python damping_vs_period_check.py --model "pk2 f_t 0.05=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr"
  python damping_vs_period_check.py --model "f_t 0.05=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05_tr" \
                                    --model "f_t 0.05 fix=stratF_full_LS_MAXWELL50_zup_pk2_ft0p05fix_tr" --runs 1-21
  options: --test Test9 | Test12 (US-1 default, US-2)   --exp EXP_DATA   --out damping_period_check
           --runs 1-25   --show 2,9,13,17,21 (runs drawn as time histories)   --tail 4.0   --t-model 0.090
           --t-test1 (default: the test's run-1 T_eq from the csv, else 0.106 for Test9 / 0.087 for Test12)
OUTPUT  <out>/damping_period_<test>.csv            one row per (model, run)
        <out>/fig_timehist_<test>_<model>.png       the --show runs: test, DEM, oscillators
        <out>/fig_summary_<test>.png                peaks per run (test, DEMs, oscillators) and the free-decay xi
"""
import os, re, glob, argparse
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks, correlate, lsim
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FS = 200.0           # common sampling rate (Hz): the test files' rate
BAND = (0.3, 30.0)


def bp(x, lo=BAND[0], hi=BAND[1], fs=FS):
    b, a = butter(4, [lo / (0.5 * fs), min(hi, 0.45 * fs) / (0.5 * fs)], "band")
    return filtfilt(b, a, x - np.mean(x))


def read_hist(path):
    """3DEC history export: two header lines, then time / value"""
    d = np.loadtxt(path, skiprows=2)
    return d[:, 0], d[:, 1]


def find_file(folder, *needles):
    for n in needles:
        hit = sorted(glob.glob(os.path.join(folder, "*{}*.csv".format(n))))
        if hit:
            return hit[0]
    return None


def load_model(folder, run):
    rd = sorted(glob.glob(os.path.join(folder, "Run{:02d}_*".format(run))))
    rd = [p for p in rd if os.path.isdir(p)]
    if not rd:
        return None
    f = rd[0]
    p12 = find_file(f, "Channel_12_AccTable"); p17 = find_file(f, "vz_acc17", "acc17"); pd_ = find_file(f, "rel_disp_top_exp_mm")
    if not (p12 and p17 and pd_):
        print("  {}: missing channel files (Channel_12_AccTable / vz_acc17 / rel_disp_top_exp_mm)".format(f)); return None
    t, v12 = read_hist(p12); _, v17 = read_hist(p17); _, dr = read_hist(pd_)
    tu = np.arange(t[0], t[-1], 1.0 / FS)
    v12 = np.interp(tu, t, v12); v17 = np.interp(tu, t, v17); dr = np.interp(tu, t, dr)
    a12 = np.gradient(v12, 1.0 / FS); a17 = np.gradient(v17, 1.0 / FS)       # exported as velocities
    return dict(t=tu, a12=a12, arel=a17 - a12, d=dr - dr[0], tag=os.path.basename(f))


def load_test(exp_dir, test, run):
    p = os.path.join(exp_dir, "{}Run{}.xlsx".format(test, run))
    if not os.path.isfile(p):
        return None
    df = pd.read_excel(p); df.columns = df.columns.map(str).str.strip()
    t = pd.to_numeric(df["Time"], errors="coerce").values
    g = lambda c: pd.to_numeric(df[c], errors="coerce").values
    a12 = g("12") * 9.81; a17 = g("17") * 9.81; d = 0.5 * (g("3") + g("4")) * 1e3        # m -> mm
    m = np.isfinite(t) & np.isfinite(a12) & np.isfinite(a17) & np.isfinite(d)
    t, a12, a17, d = t[m], a12[m], a17[m], d[m]
    tu = np.arange(t[0], t[-1], 1.0 / FS)
    a12, a17, d = (np.interp(tu, t, x) for x in (a12, a17, d))
    return dict(t=tu, a12=a12, arel=a17 - a12, d=d - np.mean(d[:int(FS)]))


def oscillator(t, ag, T, xi):
    """relative displacement (mm) of a linear oscillator under base acceleration ag (m/s2)"""
    w = 2 * np.pi / T
    _, y, _ = lsim(([[0.0, 1.0], [-w * w, -2 * xi * w]], [[0.0], [-1.0]], [[1.0, 0.0]], [[0.0]]), ag, t)
    return y * 1e3


def free_decay_xi(t, y, t0, T_guess, dur=2.0):
    """xi from the free decay after t0 (4-30 Hz peaks of |y|, log-linear fit); nan if too few peaks"""
    yb = bp(y, 4.0, 30.0); w = (t >= t0 + 0.05) & (t <= t0 + dur)
    ys, ts = yb[w], t[w]
    if len(ys) < FS * 0.5:
        return np.nan, np.nan
    pk, _ = find_peaks(np.abs(ys), distance=max(1, int(0.3 * T_guess * FS)))
    if len(pk) < 6:
        return np.nan, np.nan
    A = np.abs(ys[pk]); keep = A > 0.03 * A[0]; tp, A = ts[pk][keep], A[keep]
    if len(A) < 6:
        return np.nan, np.nan
    Tm = 2 * np.median(np.diff(tp))
    return -np.polyfit(tp, np.log(A), 1)[0] * Tm / (2 * np.pi), Tm


def parse_runs(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        elif part.strip():
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", action="append", required=True, help='"label=folder" (repeatable)')
    ap.add_argument("--test", default="Test9")
    ap.add_argument("--exp", default="EXP_DATA")
    ap.add_argument("--out", default="damping_period_check")
    ap.add_argument("--runs", default="1-25")
    ap.add_argument("--show", default="2,9,13,17,21")
    ap.add_argument("--tail", type=float, default=4.0, help="free tail after each record in the model (s)")
    ap.add_argument("--t-model", type=float, default=0.090)
    ap.add_argument("--t-test1", type=float, default=None)
    ap.add_argument("--eq-csv", default=os.path.join("exp_damping", "equivalent_damping_runs_acc.csv"))
    ap.add_argument("--ring-csv", default=os.path.join("exp_damping", "damping_runs.csv"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    models = []
    for m in a.model:
        lab, _, fol = m.partition("=")
        models.append((lab.strip(), (fol or lab).strip()))
    tname = {"Test9": "US-1", "Test12": "US-2"}.get(a.test, a.test)
    eq = {}
    if os.path.isfile(a.eq_csv):
        e = pd.read_csv(a.eq_csv)
        for _, r in e[e["test"] == tname].iterrows():
            eq[int(r["run"])] = (float(r["T_eq_s"]), float(r["xi_eq"]))
    ring = {}
    if os.path.isfile(a.ring_csv):
        e = pd.read_csv(a.ring_csv)
        for _, r in e[e["test"] == tname.replace("-", "")].iterrows():
            ring[int(r["run"])] = float(r["zeta_ring_fit"]) if pd.notna(r["zeta_ring_fit"]) else np.nan
    T1 = a.t_test1 or (eq.get(1, (None,))[0] or {"Test9": 0.106, "Test12": 0.087}.get(a.test, 0.10))
    print("{}: model T1 {:.3f} s, test run-1 period {:.3f} s; T_eq per run {}; ring-down per run {}".format(
        tname, a.t_model, T1, "found" if eq else "NOT found ({})".format(a.eq_csv),
        "found" if ring else "NOT found ({})".format(a.ring_csv)))
    runs = parse_runs(a.runs); show = parse_runs(a.show)
    rows, store = [], {}
    for run in runs:
        te = load_test(a.exp, a.test, run)
        if te is None:
            print("run {}: no {}Run{}.xlsx -- skipped".format(run, a.test, run)); continue
        ag = bp(te["a12"]); de = bp(te["d"]); ae = bp(te["arel"])
        osc = {"T{:.3f} 5%".format(a.t_model): oscillator(te["t"], ag, a.t_model, 0.05),
               "T{:.3f} 1.5%".format(a.t_model): oscillator(te["t"], ag, a.t_model, 0.015),
               "T{:.3f} 1.5%".format(T1): oscillator(te["t"], ag, T1, 0.015)}
        if run in eq:
            osc["T_eq {:.3f} xi_eq {:.1%}".format(*eq[run])] = oscillator(te["t"], ag, eq[run][0], eq[run][1])
        base = dict(run=run, test_peak_d_mm=np.max(np.abs(de)), test_peak_arel=np.max(np.abs(ae)),
                    test_T_eq=eq.get(run, (np.nan, np.nan))[0], test_xi_eq=eq.get(run, (np.nan, np.nan))[1],
                    test_xi_ringdown=ring.get(run, np.nan))
        for k, y in osc.items():
            base["osc_" + k.replace(" ", "_") + "_peak_mm"] = np.max(np.abs(y))
        for lab, fol in models:
            mo = load_model(fol, run)
            if mo is None:
                continue
            x = bp(te["a12"], 0.5, 20.0); y = bp(mo["a12"], 0.5, 20.0)
            c = correlate(x, y, mode="full"); lag = (np.argmax(c) - (len(y) - 1)) / FS
            tm = mo["t"] + lag; dm = bp(mo["d"]); am = bp(mo["arel"])
            xi_m, T_m = free_decay_xi(tm, mo["arel"], tm[-1] - a.tail, a.t_model)
            r = dict(base, model=lab, model_run=mo["tag"], lag_s=lag, model_peak_d_mm=np.max(np.abs(dm)),
                     model_peak_arel=np.max(np.abs(am)), model_xi_freedecay=xi_m, model_T_freedecay=T_m)
            r["ratio_d_model_test"] = r["model_peak_d_mm"] / r["test_peak_d_mm"]
            rows.append(r)
            if run in show:
                store.setdefault(lab, {})[run] = (te["t"], de, np.interp(te["t"], tm, dm, left=np.nan, right=np.nan), osc)
            print("run {:2d} {:14s}: peak d  test {:6.2f}  model {:6.2f} mm  | osc T{:.3f}: 5% {:5.2f}, 1.5% {:5.2f}  "
                  "T{:.3f} 1.5%: {:5.2f}{} | model free-decay xi {:.1%}".format(
                      run, lab[:14], r["test_peak_d_mm"], r["model_peak_d_mm"], a.t_model,
                      base["osc_T{:.3f}_5%_peak_mm".format(a.t_model)], base["osc_T{:.3f}_1.5%_peak_mm".format(a.t_model)],
                      T1, base["osc_T{:.3f}_1.5%_peak_mm".format(T1)],
                      ("  T_eq: {:5.2f}".format(np.max(np.abs(osc[[k for k in osc if k.startswith('T_eq')][0]])))
                       if run in eq else ""), xi_m))
    if not rows:
        print("nothing to compare"); return
    R = pd.DataFrame(rows)
    csvp = os.path.join(a.out, "damping_period_{}.csv".format(a.test)); R.to_csv(csvp, index=False)
    print("-> " + csvp)

    # ---- time histories of the --show runs, one figure per model ----
    cols = {"T{:.3f} 5%".format(a.t_model): ("#2a78d6", "--"), "T{:.3f} 1.5%".format(a.t_model): ("#2a78d6", ":"),
            "T{:.3f} 1.5%".format(T1): ("#1baf7a", "-")}
    for lab, d in store.items():
        rs = sorted(d)
        fig, ax = plt.subplots(len(rs), 1, figsize=(12, 2.9 * len(rs)), squeeze=False)
        for k, run in enumerate(rs):
            t, de, dm, osc = d[run]; ip = int(np.nanargmax(np.abs(de))); w = (t >= t[ip] - 1.5) & (t <= t[ip] + 3.0)
            axk = ax[k, 0]
            axk.plot(t[w], de[w], color="#111111", lw=1.6, label="{} measured".format(tname))
            axk.plot(t[w], dm[w], color="#c2410c", lw=1.2, label="DEM: " + lab)
            for kk, y in osc.items():
                c, ls = cols.get(kk, ("#7c3aed", "-."))
                axk.plot(t[w], y[w], color=c, ls=ls, lw=0.9, label="oscillator " + kk)
            axk.set_title("Run {}: peak {} {:.2f} mm, DEM {:.2f} mm, ".format(run, tname, np.nanmax(np.abs(de)), np.nanmax(np.abs(dm)))
                          + ", ".join("{} {:.2f}".format(kk, np.max(np.abs(y))) for kk, y in osc.items()), loc="left", fontsize=9)
            axk.set_ylabel("rel. disp. 2.06 m (mm)"); axk.grid(alpha=0.25)
            if k == 0:
                axk.legend(fontsize=7.5, ncol=3, loc="lower right")
        ax[-1, 0].set_xlabel("test time (s)")
        fig.tight_layout()
        fp = os.path.join(a.out, "fig_timehist_{}_{}.png".format(a.test, re.sub(r"[^A-Za-z0-9]+", "_", lab).strip("_")))
        fig.savefig(fp, dpi=140); plt.close(fig); print("-> " + fp)

    # ---- summary: peaks per run and free-decay damping ----
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    first = R.drop_duplicates("run").sort_values("run")
    ax[0].plot(first["run"], first["test_peak_d_mm"], "o-", color="#111111", lw=1.8, label=tname + " measured")
    pal = ["#c2410c", "#7c3aed", "#0891b2", "#a16207"]
    for i, (lab, _) in enumerate(models):
        s = R[R["model"] == lab].sort_values("run")
        if len(s):
            ax[0].plot(s["run"], s["model_peak_d_mm"], "s-", color=pal[i % len(pal)], lw=1.4, label="DEM: " + lab)
    for kk, (c, ls) in cols.items():
        col = "osc_" + kk.replace(" ", "_") + "_peak_mm"
        if col in first:
            ax[0].plot(first["run"], first[col], color=c, ls=ls, lw=1.1, label="oscillator " + kk)
    eqc = [c for c in first.columns if c.startswith("osc_T_eq")]
    if eqc:
        v = first[eqc].bfill(axis=1).iloc[:, 0]
        ax[0].plot(first["run"], v, color="#7c3aed", ls="-.", lw=1.1, label="oscillator at the test's own T_eq, xi_eq")
    ax[0].set_yscale("log"); ax[0].set_xlabel("run"); ax[0].set_ylabel("peak relative displacement at 2.06 m (mm)")
    ax[0].set_title("Peaks: damping (5 % vs 1.5 %) against period (model T1 vs test)", loc="left", fontsize=10)
    ax[0].grid(alpha=0.25, which="both"); ax[0].legend(fontsize=7.5)
    for i, (lab, _) in enumerate(models):
        s = R[R["model"] == lab].sort_values("run")
        if len(s):
            ax[1].plot(s["run"], 100 * s["model_xi_freedecay"], "s-", color=pal[i % len(pal)], label="DEM free decay: " + lab)
    ax[1].plot(first["run"], 100 * first["test_xi_ringdown"], "o", color="#111111", label=tname + " ring-down after table stop")
    ax[1].plot(first["run"], 100 * first["test_xi_eq"], "o--", color="#6b7280", ms=4, lw=1,
               label=tname + " equivalent damping during shaking (forced fit)")
    ax[1].set_xlabel("run"); ax[1].set_ylabel("damping ratio (%)"); ax[1].grid(alpha=0.25); ax[1].legend(fontsize=7.5)
    ax[1].set_title("Damping: model free decay vs test", loc="left", fontsize=10)
    fig.tight_layout(); fp = os.path.join(a.out, "fig_summary_{}.png".format(a.test)); fig.savefig(fp, dpi=140)
    print("-> " + fp)


if __name__ == "__main__":
    main()
