"""t1_teff_envelopes.py -- T1 and T_eff from the experimental envelope peaks of US-1 / US-2, compared with
the backbones digitised from Fig. 13 of Moshfeghi et al. (2024). All forces without the 2.317348 factor.

SOURCES (per specimen)
  fig13     Fig. 13 digitised. US-1: us1_fig13_digitised.csv (drift_pct, F_kN, d_mm), F / 2.317348.
            US-2: digitised here from the published figure (dashed grey), drift x 20.6 mm/% (d at 2.06 m), F / 2.317348.
  peaks     per-sign peaks of each run, NOT simultaneous (d_max_pos with F_max_pos ...) = the Fig. 13 convention,
            from envelope_peaks_US{1,2}.csv (exp_envelope_peaks.py, processed_globalzero / 2.317348)
  simult    F at the instant of d_max, per sign (simultaneous pairs), same CSVs
  both      the two simultaneous pairs of every run and sign together: (d at F_max, F_max) and (d_max, F at d_max);
            F_y, d_y then come from the force-peak pairs and d_u from the displacement-peak pairs
  Each peak set twice: d at the GLOBAL zero as stored, and d relative to the run start (d - d_start).
  Runs 1-24 for both specimens (what Fig. 13 shows); US-2 also with run 25, which the figure leaves out.
RECIPE (the one used for the 20 Oct handover, so the numbers are comparable)
  K_el   least-squares slope through the origin of the points with 0.1 <= |d| <= 1.5 mm, both signs pooled
         (points with d and F of opposite sign are dropped; if a sign has none in range, its first point is used)
  F_y,d_y  peak |F| per sign and the d there;  d_u = largest |d| per sign (end of the envelope);  mean of the two signs
  T1_m   = 2 pi sqrt(m / K_el), m = 1133.3 kg = sum of the file masses once the factor is removed
           (408.52 + 320.66 + 404.13); m_impl = K_el (T1_Table5 / 2 pi)^2 is the mass the measured period implies
  T_y    = T1 sqrt(K_el / K_y), K_y = F_y / d_y;   T_eff = T_y sqrt(d_u / d_y) = T1 sqrt(K_el d_u / F_y)
         computed with T1 = Table 5 (US-1 0.091 s, US-2 0.087 s), i.e. anchored on the measured period
OUTPUT <out>/t1_teff_envelopes.csv, fig_t1_teff_envelopes.png, us2_fig13_digitised.csv
USAGE  python t1_teff_envelopes.py <envelope_dir> <us1_fig13_csv> <out_dir>
"""
import sys, os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FAC = 2.317348
H_MM = 20.6                       # mm per % drift (top-quarter sensor at 2.06 m)
M_FILE = 408.52 + 320.66 + 404.13
T1_TAB5 = {"US1": 0.091, "US2": 0.087}
# US-2, Fig. 13 dashed grey, read from the published image (drift %, F kN at paper scale).
# Axis calibration from the gridlines: 242 px per 1 % drift, 8.43 px per kN. Between 0 and -22 kN the
# negative branch is hidden under the red ST curves; the positive branch below 8 kN likewise.
US2_FIG13 = [(-0.390, -27.0), (-0.340, -28.1), (-0.265, -27.51), (-0.108, -22.44), (0.0, 0.0),
             (0.033, 8.66), (0.045, 12.15), (0.062, 15.61), (0.088, 19.03), (0.130, 22.30),
             (0.190, 25.31), (0.266, 28.06), (0.368, 30.30), (0.465, 27.30)]


def recipe(dp, Fp, dn, Fn, T1):
    """dp, Fp: positive branch points (d > 0) in order of |d|; dn, Fn likewise for the negative branch."""
    out = {}
    sel_d, sel_F = [], []
    for d, F in ((dp, Fp), (dn, Fn)):
        ok = (d * F > 0) & (np.abs(d) >= 0.1) & (np.abs(d) <= 1.5)
        if ok.sum() == 0:
            j = np.argmax(d * F > 0); ok = np.zeros_like(d, bool); ok[j] = True
        sel_d.append(d[ok]); sel_F.append(F[ok])
    dd, FF = np.concatenate(sel_d), np.concatenate(sel_F)
    K = float(np.sum(dd * FF) / np.sum(dd * dd))
    Kp = float(np.sum(sel_d[0] * sel_F[0]) / np.sum(sel_d[0] ** 2)); Kn = float(np.sum(sel_d[1] * sel_F[1]) / np.sum(sel_d[1] ** 2))
    side = {}
    for lab, d, F in (("pos", dp, Fp), ("neg", dn, Fn)):
        i = int(np.argmax(np.abs(F))); side[lab] = (abs(F[i]), abs(d[i]), float(np.max(np.abs(d))))
    Fy = np.mean([side[s][0] for s in side]); dy = np.mean([side[s][1] for s in side]); du = np.mean([side[s][2] for s in side])
    Ky = Fy / dy
    out.update(K_el=K, K_pos=Kp, K_neg=Kn, Fy=Fy, dy=dy, du=du, Ky=Ky, mu=du / dy,
               T1_m=2 * np.pi * np.sqrt(M_FILE / (K * 1e6)), m_impl=K * 1e6 * (T1 / 2 / np.pi) ** 2,
               T1=T1, Ty=T1 * np.sqrt(K / Ky), Teff=T1 * np.sqrt(K * du / Fy))
    for s in side:
        Fys, dys, dus = side[s]
        out["Teff_" + s] = T1 * np.sqrt(K * dus / Fys); out["Fy_" + s] = Fys; out["dy_" + s] = dys; out["du_" + s] = dus
    return out


def branches_from_curve(d, F):
    o = np.argsort(d); d, F = d[o], F[o]
    p = d > 0; n = d < 0
    return d[p], F[p], d[n][::-1], F[n][::-1]


def branches_from_peaks(rows, kind, rel):
    dp, Fp, dn, Fn = [], [], [], []
    for r in rows:
        off = float(r["d_start_mm"]) if rel else 0.0
        if kind == "both":          # both simultaneous pairs: (d at F_max, F_max) and (d_max, F at d_max)
            for key_d, key_F in (("d_at_Fmax_pos_mm", "F_max_pos_kN"), ("d_max_pos_mm", "F_at_dmax_pos_kN")):
                dp.append(float(r[key_d]) - off); Fp.append(float(r[key_F]))
            for key_d, key_F in (("d_at_Fmax_neg_mm", "F_max_neg_kN"), ("d_max_neg_mm", "F_at_dmax_neg_kN")):
                dn.append(float(r[key_d]) - off); Fn.append(float(r[key_F]))
            continue
        if kind == "peaks":
            dp.append(float(r["d_max_pos_mm"]) - off); Fp.append(float(r["F_max_pos_kN"]))
            dn.append(float(r["d_max_neg_mm"]) - off); Fn.append(float(r["F_max_neg_kN"]))
        else:
            dp.append(float(r["d_max_pos_mm"]) - off); Fp.append(float(r["F_at_dmax_pos_kN"]))
            dn.append(float(r["d_max_neg_mm"]) - off); Fn.append(float(r["F_at_dmax_neg_kN"]))
    # run order = loading order; the envelope is read in run order, the last run gives d_u
    return np.array(dp), np.array(Fp), np.array(dn), np.array(Fn)


def main():
    env, f13, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    g = np.genfromtxt(f13, delimiter=",", names=True)
    fig13 = {"US1": (g["d_mm"], g["F_kN"] / FAC),
             "US2": (np.array([p[0] for p in US2_FIG13]) * H_MM, np.array([p[1] for p in US2_FIG13]) / FAC)}
    with open(os.path.join(out, "us2_fig13_digitised.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["drift_pct", "F_kN_paper_scale", "d_mm", "F_kN_div_2317"])
        for dr, F in US2_FIG13: w.writerow([dr, F, "{:.3f}".format(dr * H_MM), "{:.3f}".format(F / FAC)])
    res = []; curves = {}
    for spec in ("US1", "US2"):
        T1 = T1_TAB5[spec]
        d, F = fig13[spec]; br = branches_from_curve(d, F)
        res.append((spec, "Fig. 13 digitised", recipe(*br, T1))); curves[(spec, "fig13")] = br
        allrows = list(csv.DictReader(open(os.path.join(env, "envelope_peaks_{}.csv".format(spec)))))
        nmax = int(allrows[-1]["run"])
        # Fig. 13 stops at run 24 for both specimens (US-2 run 25 is in the files but not in the figure)
        for last in sorted({24, nmax}):
            rows = [r for r in allrows if int(r["run"]) <= last]
            for kind, lab in (("peaks", "peaks per sign (non-simult.)"), ("simult", "F at d_max (simultaneous)"), ("both", "both simultaneous pairs (F_max,d) + (d_max,F)")):
                for rel in (False, True):
                    br = branches_from_peaks(rows, kind, rel)
                    name = "{}, {}, runs 1-{}".format(lab, "run-rel. d" if rel else "global-zero d", last)
                    res.append((spec, name, recipe(*br, T1)))
                    if last == 24: curves[(spec, kind, rel)] = br
    keys = ["K_el", "K_pos", "K_neg", "Fy", "dy", "du", "mu", "T1_m", "m_impl", "T1", "Ty", "Teff", "Teff_pos", "Teff_neg"]
    with open(os.path.join(out, "t1_teff_envelopes.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["specimen", "source"] + keys)
        for spec, name, r in res: w.writerow([spec, name] + ["{:.4f}".format(r[k]) for k in keys])
    print("{:4s} {:58s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>5s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s}".format(
        "spec", "source", "K_el", "K+", "K-", "F_y", "d_y", "d_u", "mu", "T1_m", "m_impl", "T_y", "T_eff", "Teff+", "Teff-"))
    for spec, name, r in res:
        print("{:4s} {:58s} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:5.2f} {:6.3f} {:6.0f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}".format(
            spec, name, r["K_el"], r["K_pos"], r["K_neg"], r["Fy"], r["dy"], r["du"], r["mu"], r["T1_m"], r["m_impl"], r["Ty"], r["Teff"], r["Teff_pos"], r["Teff_neg"]))
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
    for a, spec in zip(ax, ("US1", "US2")):
        dp, Fp, dn, Fn = curves[(spec, "fig13")]
        a.plot(np.r_[dn[::-1], 0, dp], np.r_[Fn[::-1], 0, Fp], "-", color="#b4491a", lw=2,
               label="Fig. 13 digitised / 2.317" + (" (points)" if spec == "US2" else ""), marker="o" if spec == "US2" else None, ms=3)
        for kind, c, mk, lab in (("peaks", "#2b5f8a", "s", "peaks per sign, global zero"), ("simult", "#4a8f5a", "o", "F at d_max, global zero")):
            dp, Fp, dn, Fn = curves[(spec, kind, False)]
            a.plot(dp, Fp, mk, color=c, ms=4, mec="k", mew=0.3, label=lab); a.plot(dn, Fn, mk, color=c, ms=4, mec="k", mew=0.3)
        r = [x for x in res if x[0] == spec and x[1] == "Fig. 13 digitised"][0][2]
        dd = np.linspace(-3, 3, 10); a.plot(dd, r["K_el"] * dd, ":", color="#b4491a", lw=1)
        a.axhline(0, color="k", lw=0.5); a.axvline(0, color="k", lw=0.5); a.grid(alpha=0.25)
        a.set_title("{}: runs 1-24, T1 (Table 5) = {:.3f} s".format(spec, T1_TAB5[spec]), loc="left", fontsize=10)
        a.set_xlabel("d at 2.06 m (mm)")
    ax[0].set_ylabel("F (kN), factor 2.317 removed"); ax[0].legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_t1_teff_envelopes.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
