"""periods_on_spectra.py -- the 25 scaled target spectra of the protocol with the periods that matter marked:
the previous Route 2 model's bilinear T1 / T_eff, the test's own T1 / T_eff (effective mass L = sum m phi,
both simultaneous pairs), and the Table 5 periods.

SPECTRA  spectrum_<REC>.csv (T_s, Sd_m), 5 %-damped record spectra, scaled by the protocol factors
         (handover of 20 Oct: HU12 runs 1-21, EC40 runs 3/6/8/11, FR76 runs 22-25).
PERIODS  Route 2 model (bilinear_periods.json, m_eff 1635 kg, old bonded wall): T1 0.112, T_eff 0.241 s
         experiment, this session (t1_teff_with_Leff.csv, L = 1200 kg, "both simultaneous pairs", runs 1-24)
OUTPUT   <out>/fig_periods_on_spectra.png, periods_compared.csv
USAGE    python periods_on_spectra.py <repo_dir> <t1_teff_with_Leff.csv> <out_dir>
"""
import sys, os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROTOCOL = [(1, "HU12", 0.50), (2, "HU12", 0.75), (3, "EC40", 0.20), (4, "HU12", 1.00), (5, "HU12", 1.25),
            (6, "EC40", 0.30), (7, "HU12", 1.50), (8, "EC40", 0.40), (9, "HU12", 1.75), (10, "HU12", 2.00),
            (11, "EC40", 0.50), (12, "HU12", 2.25), (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00),
            (16, "HU12", 3.50), (17, "HU12", 4.00), (18, "HU12", 4.50), (19, "HU12", 5.00), (20, "HU12", 5.50),
            (21, "HU12", 6.00), (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75), (25, "FR76", 2.00)]
REC_COL = {"HU12": "#2b5f8a", "EC40": "#e08a1e", "FR76": "#b4491a"}


def main():
    repo, lcsv, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    spec = {r: np.genfromtxt(os.path.join(repo, "spectrum_{}.csv".format(r)), delimiter=",", names=True) for r in REC_COL}
    model = json.load(open(os.path.join(repo, "Claude outputs", "bilinear_periods.json")))["model"]
    rows = list(csv.DictReader(open(lcsv)))
    exp = {}
    for r in rows:
        if r["source"].startswith("both simultaneous") and "global-zero d, runs 1-24" in r["source"] and r["mass_set"].startswith("v2"):
            exp[r["specimen"]] = (float(r["T1_L_s"]), float(r["Teff_L_s"]), float(r["Ty_L_s"]))
    periods = [("Route 2 model (bonded wall, m 1635 kg)", model["T1_s"], model["Teff_s"], "k"),
               ("ft0 / c0.2 model v5 (L 1200 kg, d_u 76 mm)", 0.0964, 0.5645, "#6a4c93"),
               ("ft0 / c0.2 model v5 (d_u = test's 27.4 mm)", 0.0964, 0.3397, "#b08fd8"),
               ("US-1 test, L = 1200 kg", exp["US1"][0], exp["US1"][1], "#1f6f8b"),
               ("US-2 test, L = 1200 kg", exp["US2"][0], exp["US2"][1], "#b4491a")]
    tab5 = {"US-1": 0.091, "US-2": 0.087}
    # table of Sd at each period, every run
    with open(os.path.join(out, "periods_compared.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "T1_s", "Teff_s", "Ty_s"])
        w.writerow(["Route 2 model", "{:.4f}".format(model["T1_s"]), "{:.4f}".format(model["Teff_s"]), ""])
        w.writerow(["ft0 model v5, L=1200, d_u 76 mm", "0.0964", "0.5645", "0.2247"])
        w.writerow(["ft0 model v5, L=1200, d_u 27.4 mm", "0.0964", "0.3397", "0.2247"])
        for k in ("US1", "US2"): w.writerow([k + " test (both pairs, L=1200)", "{:.4f}".format(exp[k][0]), "{:.4f}".format(exp[k][1]), "{:.4f}".format(exp[k][2])])
        w.writerow([]); w.writerow(["run", "record", "scale"] + ["Sd({}, {}={:.3f})_mm".format(n.split(" ")[0] + n.split(" ")[1], lab, T) for n, T1, Te, c in periods for lab, T in (("T1", T1), ("Teff", Te))])
        for n, rec, s in PROTOCOL:
            w.writerow([n, rec, s] + ["{:.2f}".format(s * np.interp(T, spec[rec]["T_s"], spec[rec]["Sd_m"]) * 1e3) for _, T1, Te, _ in periods for T in (T1, Te)])
    print("{:40s} {:>7s} {:>7s}".format("case", "T1", "T_eff"))
    for n, T1, Te, c in periods: print("{:40s} {:7.3f} {:7.3f}".format(n, T1, Te))
    print("\nSd at T1/Teff (mm):  run     " + "   ".join(n.split(" (")[0][:14] for n, *_ in periods))
    for n, rec, s in PROTOCOL:
        if n in (1, 5, 10, 15, 19, 21, 22, 23, 24, 25):
            v = [s * np.interp(T, spec[rec]["T_s"], spec[rec]["Sd_m"]) * 1e3 for _, T1, Te, _ in periods for T in (T1, Te)]
            print("  {:2d} {} x{:.2f}  ".format(n, rec, s) + "   ".join("{:5.1f}/{:5.1f}".format(v[2 * k], v[2 * k + 1]) for k in range(len(periods))))
    # figure
    fig, ax = plt.subplots(figsize=(17, 6.5))
    for n, rec, s in PROTOCOL:
        T = spec[rec]["T_s"]; Sd = s * spec[rec]["Sd_m"] * 1e3; m = T <= 0.6
        ax.plot(T[m], Sd[m], "-", color=REC_COL[rec], lw=0.9, alpha=0.35 + 0.6 * n / 25)
        ax.text(0.602, min(Sd[np.searchsorted(T, 0.6)], 49.5), str(n), fontsize=6, color=REC_COL[rec], va="center", clip_on=True)
    ls = {"T1": "-", "Teff": "--"}
    for name, T1, Te, c in periods:
        ax.axvline(T1, color=c, lw=1.6, ls="-"); ax.axvline(Te, color=c, lw=1.6, ls="--")
        ax.text(T1 + 0.002, 48, "T1 {:.3f}".format(T1), color=c, fontsize=8, rotation=90, va="top")
        ax.text(Te + 0.002, 48, "T_eff {:.3f}".format(Te), color=c, fontsize=8, rotation=90, va="top")
    for k, T in tab5.items():
        ax.axvline(T, color="grey", lw=0.8, ls=":"); ax.text(T - 0.0045, 2, "Table 5 " + k, color="grey", fontsize=7, rotation=90, va="bottom")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, lw=1.6, label=n) for n, _, _, c in periods]
    h += [Line2D([], [], color="k", lw=1.2, ls="-", label="T1"), Line2D([], [], color="k", lw=1.2, ls="--", label="T_eff")]
    h += [Line2D([], [], color=REC_COL[r], lw=1.5, label=r) for r in REC_COL]
    ax.legend(handles=h, fontsize=7, loc="upper left", ncol=2)
    ax.set_xlim(0, 0.62); ax.set_ylim(0, 50); ax.set_clip_on(True); ax.grid(alpha=0.25)
    ax.set_xlabel("oscillator period T (s)"); ax.set_ylabel("scaled spectral displacement Sd (mm), 5 % damping")
    ax.set_title("Protocol spectra (25 runs, run number at the right edge) with the model's and the tests' T1 / T_eff", loc="left", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_periods_on_spectra.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
