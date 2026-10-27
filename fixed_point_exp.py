"""fixed_point_exp.py -- the fixed point d* = Sd(T_sec(d*)) run on the TEST's own loop backbone, fed with the
TABLE motion the test actually received (ch12), and compared with the test's own peak displacements.

This checks the method independently of the DEM model: if the fixed point cannot predict the test from
the test's own backbone and the test's own table motion, it cannot be expected to do better on the model.

INPUT  <raw_dir>/Test9Run<N>.xlsx            ch12 = table acceleration [g]
       <bb_dir>/loop_backbone_curve.csv      from loop_backbone.py (sign, d_mm, F_p90_kN, ...)
       <bb_dir>/loop_backbone_points.csv     measured peaks per run
OPTIONS --m <kg>        effective mass (default 1582 = median pp-secant m_eff of runs 1, 2, 6)
        --kel <kN/mm>   elastic stiffness below the first envelope bin (default 7.53, run 1 pp-secant)
        --xi <0.05>     spectral damping
METHOD T_sec(d) = 2 pi sqrt(m d / F(d)), F(d) = p90 branch envelope made non-decreasing (cumulative max);
       for d below the first bin F = K_el d. Sd(T) of ch12 by the exact piecewise-linear recurrence
       (Nigam & Jennings 1969) on a 400-point log grid 0.04-3 s. g(d) = Sd(T_sec(d)) - d scanned on
       6000 points 0.02-120 mm; first + -> - zero = stable d*, the next - -> + zero = tipping.
       Three backbones: '+' side, '-' side, mean of the two.
OUTPUT <out>/fixed_point_exp.csv, fig_fixed_point_exp.png
"""
import sys, os, re, glob, csv
import numpy as np
from scipy.signal import ss2tf, lfilter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G0 = 9.81


def arg(flag, default):
    return type(default)(sys.argv[sys.argv.index(flag) + 1]) if flag in sys.argv else default


def load_ch12(path):
    import openpyxl
    it = openpyxl.load_workbook(path, read_only=True).worksheets[0].iter_rows(values_only=True)
    h = list(next(it))
    jt = [j for j, x in enumerate(h) if isinstance(x, str) and x.strip().lower() == "time"][0]
    j12 = [j for j, x in enumerate(h) if x is not None and not isinstance(x, str) and int(float(x)) == 12][0]
    A = np.array([[r[jt], np.nan if r[j12] is None else r[j12]] for r in it if r[jt] is not None], float)
    a = np.nan_to_num(A[:, 1]); a -= np.median(a[: int(2.0 / np.median(np.diff(A[:, 0])))])
    return np.median(np.diff(A[:, 0])), a * G0


def sd_spectrum(ag, dt, T, xi):
    """relative-displacement spectrum, exact recurrence for piecewise-linear ground acceleration."""
    out = np.empty(len(T))
    for k, Tn in enumerate(T):
        w = 2 * np.pi / Tn; wd = w * np.sqrt(1 - xi ** 2); e = np.exp(-xi * w * dt)
        s, c = np.sin(wd * dt), np.cos(wd * dt)
        A11 = e * (xi / np.sqrt(1 - xi ** 2) * s + c); A12 = e / wd * s
        A21 = -w / np.sqrt(1 - xi ** 2) * e * s; A22 = e * (c - xi / np.sqrt(1 - xi ** 2) * s)
        t1 = (2 * xi ** 2 - 1) / (w ** 2 * dt); t2 = 2 * xi / (w ** 3 * dt)
        B11 = e * ((t1 + xi / w) * s / wd + (t2 + 1 / w ** 2) * c) - t2
        B12 = -e * (t1 * s / wd + t2 * c) - 1 / w ** 2 + t2
        B21 = e * ((t1 + xi / w) * (c - xi / np.sqrt(1 - xi ** 2) * s) - (t2 + 1 / w ** 2) * (wd * s + xi * w * c)) + 1 / (w ** 2 * dt)
        B22 = -e * (t1 * (c - xi / np.sqrt(1 - xi ** 2) * s) - t2 * (wd * s + xi * w * c)) - 1 / (w ** 2 * dt)
        # x_{k+1} = A x_k + B1 p_k + B2 p_{k+1}; with y_k = x_k - B2 p_k this is a standard (A, A B2 + B1, C, C B2) system
        A = np.array([[A11, A12], [A21, A22]]); B1 = np.array([B11, B21]); B2 = np.array([B12, B22])
        b, a = ss2tf(A, (A @ B2 + B1).reshape(2, 1), np.array([[1.0, 0.0]]), np.array([[B2[0]]]))
        mx = np.abs(lfilter(b[0], a, -ag)).max()
        out[k] = mx
    return out


def backbone(rows, kel):
    """p90 envelope -> non-decreasing (cumulative max), never stiffer than K_el (below ~0.15 mm the
    envelope carries force with no displacement: potentiometer stiction), flat beyond the last bin."""
    d, F = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
    o = np.argsort(d); d, F = d[o], np.maximum.accumulate(np.maximum(F[o], 0.0))
    d = np.concatenate([[1e-4], d]); F = np.concatenate([[kel * 1e-4], F])
    def f(x):
        x = np.asarray(x, float)
        return np.minimum(np.interp(x, d, F), kel * x)
    return f


def fixed_point(dgrid, G):
    g = G - dgrid; s = np.sign(g)
    i = np.where((s[:-1] > 0) & (s[1:] <= 0))[0]
    if len(i) == 0: return np.nan, np.nan
    k = i[0]; dstar = dgrid[k] - g[k] * (dgrid[k + 1] - dgrid[k]) / (g[k + 1] - g[k])
    j = np.where((s[k + 1:-1] < 0) & (s[k + 2:] >= 0))[0]
    tip = dgrid[k + 1 + j[0]] if len(j) else np.nan
    return dstar, tip


def main():
    raw_dir, bb_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    m = arg("--m", 1582.0); kel = arg("--kel", 7.53); xi = arg("--xi", 0.05)
    cur = {"pos": [], "neg": []}
    for r in csv.DictReader(open(os.path.join(bb_dir, "loop_backbone_curve.csv"))):
        cur[r["sign"]].append((float(r["d_mm"]), float(r["F_p90_kN"])))
    meas = {int(r["run"]): r for r in csv.DictReader(open(os.path.join(bb_dir, "loop_backbone_points.csv")))}
    fb = {"+": backbone(cur["pos"], kel), "-": backbone(cur["neg"], kel)}
    fb["mean"] = lambda x: 0.5 * (fb["+"](x) + fb["-"](x))
    dgrid = np.linspace(0.02, 120, 6000)
    Tsec = {k: 2 * np.pi * np.sqrt(m * dgrid * 1e-3 / (f(dgrid) * 1e3)) for k, f in fb.items()}
    Tg = np.geomspace(0.04, 3.0, 400)
    rows, SD = [], {}
    files = sorted(glob.glob(os.path.join(raw_dir, "Test9Run*.xlsx")), key=lambda p: int(re.search(r"Run(\d+)\.xlsx", p).group(1)))
    for p in files:
        n = int(re.search(r"Run(\d+)\.xlsx", p).group(1))
        dt, ag = load_ch12(p); sd = sd_spectrum(ag, dt, Tg, xi) * 1e3; SD[n] = sd
        mr = meas[n]; dp, dn = float(mr["d_pk_pos"]), -float(mr["d_pk_neg"])
        row = [n, np.abs(ag).max() / G0, np.interp(0.0911, Tg, sd), dp, dn, max(dp, dn)]
        for k in ("+", "-", "mean"):
            G = np.interp(Tsec[k], Tg, sd)
            ds, tip = fixed_point(dgrid, G)
            row += [ds, tip, float(np.interp(ds, dgrid, Tsec[k])) if np.isfinite(ds) else np.nan]
        rows.append(row)
        print("run {:2d} PGA {:.3f} g  Sd(T1) {:6.2f} mm | meas +{:5.2f} -{:5.2f} | d* + {:6.2f}  - {:6.2f}  mean {:6.2f} (T* {:.3f} s, tip {})".format(
            n, row[1], row[2], dp, dn, row[6], row[9], row[12], row[14], "{:.1f}".format(row[13]) if np.isfinite(row[13]) else "-"))
    hdr = ["run", "PGA_g", "Sd_T1_mm", "meas_pos_mm", "meas_neg_mm", "meas_max_mm"]
    for k in ("pos", "neg", "mean"): hdr += ["dstar_" + k, "tip_" + k, "Tstar_" + k]
    with open(os.path.join(out, "fixed_point_exp.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(hdr)
        for r in rows: w.writerow([r[0]] + ["{:.4f}".format(x) for x in r[1:]])
    R = np.array(rows, float)
    # figure
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    a = ax[0]
    for k, c in (("+", "#1f6f8b"), ("-", "#b4491a"), ("mean", "k")):
        a.semilogx(dgrid, Tsec[k], color=c, lw=1.6, label="T_sec, backbone " + k)
    a.axhline(0.0911, color="grey", ls=":", lw=0.8); a.text(0.03, 0.095, "T1 = 0.091 s (Table 5)", fontsize=7, color="grey")
    a.set_xlabel("d (mm)"); a.set_ylabel("T_sec (s)"); a.grid(alpha=0.25, which="both"); a.legend(fontsize=8)
    a.set_title("Secant period on the loop backbone, m_eff = {:.0f} kg".format(m), loc="left", fontsize=10)
    a = ax[1]
    for n, c in ((19, "#2b5f8a"), (21, "#4a8f5a"), (23, "#c08a20"), (24, "#b4491a")):
        G = np.interp(Tsec["mean"], Tg, SD[n]); a.loglog(dgrid, G, color=c, lw=1.4, label="run {}".format(n))
        a.loglog(meas_max := R[R[:, 0] == n, 5][0], meas_max, "o", color=c, mec="k")
    a.loglog(dgrid, dgrid, "k--", lw=0.8, label="45 deg")
    a.set_xlim(0.1, 120); a.set_ylim(0.1, 300); a.grid(alpha=0.25, which="both"); a.legend(fontsize=8)
    a.set_xlabel("d (mm)"); a.set_ylabel("Sd(T_sec(d)) (mm)")
    a.set_title("Demand along the backbone (mean); o = measured max |d| on the 45 deg line", loc="left", fontsize=10)
    a = ax[2]
    a.loglog(R[:, 5], R[:, 12], "ks", ms=5, label="d* mean backbone vs measured max")
    a.loglog(R[:, 3], R[:, 6], "o", color="#1f6f8b", ms=4, label="d* (+ backbone) vs measured +")
    a.loglog(R[:, 4], R[:, 9], "o", color="#b4491a", ms=4, label="d* (- backbone) vs measured -")
    for r in R: a.annotate(str(int(r[0])), (r[5], r[12]), fontsize=6, xytext=(3, 2), textcoords="offset points")
    lim = [0.05, 200]; a.plot(lim, lim, "k--", lw=0.8)
    for f_ in (0.5, 2): a.plot(lim, [x * f_ for x in lim], ":", color="grey", lw=0.7)
    a.set_xlim(lim); a.set_ylim(lim); a.grid(alpha=0.25, which="both"); a.legend(fontsize=7)
    a.set_xlabel("measured peak d (mm)"); a.set_ylabel("fixed point d* (mm)")
    a.set_title("Test backbone + achieved table motion -> test peaks", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_fixed_point_exp.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
