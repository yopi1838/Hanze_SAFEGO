"""exp_effective_mass.py -- effective mass of US-1 / US-2 from the measured mode shape, and T1 / T_eff with it.

WHY  The envelope plots base shear F = sum_j m_j a_j against the displacement d at 2.06 m. In the first mode
     a_j = phi_j * a_top (phi normalised to 1 at the 2.06 m accelerometer, ch17, same height as the displacement
     sensors), so F = L * d_dotdot with  L = sum_j m_j phi_j.  Hence K_el = L w^2 and
         T1 = 2 pi sqrt(L / K_el).
     L is NOT the modal mass sum m phi^2 (that pairs with the modal force, not the base shear).
MODE SHAPE  from the raw files (EXP_DATA/Test<N>Run<r>.xlsx), channels 15/16/17 (0.66 / 1.26 / 2.06 m), in the
     wall's free vibration after the strong motion (table term below 5 % of its peak, ch12). Band-pass 5-20 Hz,
     phi_j = <a_j a_17> / <a_17 a_17> (least squares, keeps the sign). Also during the strong motion for comparison.
MASSES  two sets: the group's sheet (v2: 455.03 / 670.07 / 606.47 kg, slabs + beam lumped on the wall sensors) and
     the old files' masses without the factor (408.52 / 320.66 / 404.13 kg, wall only).
T1, T_eff  T1 = 2 pi sqrt(L / K_el) for each envelope source in t1_teff_envelopes.csv; T_eff = T1 sqrt(K_el d_u / F_y)
     with the same source's K_el, d_u, F_y; compared with Table 5.
OUTPUT <out>/mode_shape_runs.csv, effective_mass_summary.csv, t1_teff_with_Leff.csv, fig_effective_mass.png
USAGE  python exp_effective_mass.py <raw_dir> <t1_teff_envelopes.csv> <out_dir>
"""
import sys, os, re, glob, csv
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MASSES = {"v2 sheet (1732 kg)": {15: 455.03, 16: 670.07, 17: 606.47}, "old files / 2.317 (1133 kg)": {15: 408.52, 16: 320.66, 17: 404.13}}
TESTS = (("US1", "Test9", 0.091), ("US2", "Test12", 0.087))
HEIGHTS = {15: 0.66, 16: 1.26, 17: 2.06}


def load_raw(path):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True); hdr = list(next(it))
    col = {}
    for j, h in enumerate(hdr):
        if isinstance(h, str) and h.strip().lower() == "time": col["t"] = j
        elif h is not None and not isinstance(h, str): col[int(float(h))] = j
    rows = [r for r in it if r[col["t"]] is not None]
    A = np.array([[np.nan if (j >= len(r) or r[j] is None) else r[j] for j in (col["t"], col[12], col[15], col[16], col[17])] for r in rows], float)
    A = np.nan_to_num(A)
    return A[:, 0], A[:, 1], {15: A[:, 2], 16: A[:, 3], 17: A[:, 4]}


def bp(x, fs, lo, hi):
    b, a = butter(4, [lo / (0.5 * fs), hi / (0.5 * fs)], "band"); return filtfilt(b, a, x - np.mean(x))


def mode_shape(path):
    t, a12, ch = load_raw(path); fs = 1.0 / np.median(np.diff(t))
    b, a = butter(4, 4.0 / (0.5 * fs)); env = np.abs(hilbert(filtfilt(b, a, a12 - np.median(a12[: int(2 * fs)]))))
    strong = env > 0.05 * env.max(); i_on, i_off = np.where(strong)[0][[0, -1]]
    free = np.zeros(len(t), bool); free[i_off + int(0.1 * fs): i_off + int(3.0 * fs)] = True   # first 3 s of free vibration
    forced = np.zeros(len(t), bool); forced[i_on:i_off] = True
    f = {k: bp(v, fs, 5.0, 20.0) for k, v in ch.items()}
    out = {}
    for lab, m in (("free", free), ("forced", forced)):
        ref = f[17][m]
        out[lab] = {k: float(np.dot(f[k][m], ref) / np.dot(ref, ref)) for k in (15, 16, 17)}
        out[lab]["rms17_g"] = float(np.std(ref))
    return out


def main():
    raw, tt_csv, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    shapes = {}
    with open(os.path.join(out, "mode_shape_runs.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["specimen", "run", "phi15_free", "phi16_free", "rms17_free_g", "phi15_forced", "phi16_forced", "rms17_forced_g"])
        for spec, test, T1 in TESTS:
            files = sorted(glob.glob(os.path.join(raw, test + "Run*.xlsx")), key=lambda p: int(re.search(r"Run(\d+)\.xlsx", p).group(1)))
            print("== {}".format(spec))
            for p in files:
                n = int(re.search(r"Run(\d+)\.xlsx", p).group(1)); s = mode_shape(p); shapes[(spec, n)] = s
                w.writerow([spec, n] + ["{:.4f}".format(x) for x in (s["free"][15], s["free"][16], s["free"]["rms17_g"], s["forced"][15], s["forced"][16], s["forced"]["rms17_g"])])
                print("  run {:2d}  free: phi = [{:.3f} {:.3f} 1] (rms17 {:.4f} g)   forced: phi = [{:.3f} {:.3f} 1]".format(
                    n, s["free"][15], s["free"][16], s["free"]["rms17_g"], s["forced"][15], s["forced"][16]))
    # effective mass per specimen: median mode shape over the runs with a usable ring-down (rms17 > 2e-4 g), free segment
    Leff = {}
    with open(os.path.join(out, "effective_mass_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["specimen", "segment", "phi15_median", "phi16_median", "n_runs", "mass_set", "L_eff_kg", "modal_mass_kg", "sum_m_kg"])
        for spec, test, T1 in TESTS:
            for seg in ("free", "forced"):
                S = [shapes[k][seg] for k in shapes if k[0] == spec and shapes[k][seg]["rms17_g"] > 2e-4]
                p15, p16 = np.median([s[15] for s in S]), np.median([s[16] for s in S])
                for ms, M in MASSES.items():
                    L = M[15] * p15 + M[16] * p16 + M[17]; mm = M[15] * p15 ** 2 + M[16] * p16 ** 2 + M[17]
                    w.writerow([spec, seg, "{:.4f}".format(p15), "{:.4f}".format(p16), len(S), ms, "{:.1f}".format(L), "{:.1f}".format(mm), "{:.1f}".format(sum(M.values()))])
                    if seg == "free": Leff[(spec, ms)] = L
                    print("{} {:6s}: phi = [{:.3f} {:.3f} 1] ({} runs)  {:28s} L = sum m phi = {:6.0f} kg   (modal mass {:5.0f}, sum m {:5.0f})".format(
                        spec, seg, p15, p16, len(S), ms, L, mm, sum(M.values())))
    # T1 and T_eff with L
    rows = list(csv.DictReader(open(tt_csv)))
    with open(os.path.join(out, "t1_teff_with_Leff.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["specimen", "source", "mass_set", "L_eff_kg", "K_el_kN_mm", "T1_L_s", "T1_table5_s", "T1_ratio", "Fy_kN", "dy_mm", "du_mm", "Ty_L_s", "Teff_L_s", "Teff_table5_s"])
        print("\n{:4s} {:58s} {:>6s} {:>6s} {:>7s} {:>6s} {:>7s} {:>7s} {:>7s}".format("spec", "source", "L kg", "K_el", "T1(L)", "T1/T5", "T_y(L)", "Teff(L)", "Teff(T5)"))
        for r in rows:
            spec = r["specimen"]; T5 = {"US1": 0.091, "US2": 0.087}[spec]; K = float(r["K_el"]); Fy, dy, du = float(r["Fy"]), float(r["dy"]), float(r["du"])
            for ms in MASSES:
                L = Leff[(spec, ms)]; T1 = 2 * np.pi * np.sqrt(L / (K * 1e6)); Ty = T1 * np.sqrt(K * dy / Fy); Te = T1 * np.sqrt(K * du / Fy)
                w.writerow([spec, r["source"], ms, "{:.1f}".format(L), "{:.3f}".format(K), "{:.4f}".format(T1), T5, "{:.3f}".format(T1 / T5), Fy, dy, du, "{:.4f}".format(Ty), "{:.4f}".format(Te), r["Teff"]])
                if ms.startswith("v2"):
                    print("{:4s} {:58s} {:6.0f} {:6.2f} {:7.4f} {:6.2f} {:7.4f} {:7.4f} {:7.4f}".format(spec, r["source"], L, K, T1, T1 / T5, Ty, Te, float(r["Teff"])))
    # figure: mode shapes per run + L per run
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for spec, test, T1 in TESTS:
        runs = sorted(n for s_, n in shapes if s_ == spec); c = "#1f6f8b" if spec == "US1" else "#b4491a"
        for n in runs:
            s = shapes[(spec, n)]["free"]
            if s["rms17_g"] > 2e-4:
                ax[0].plot([0, s[15], s[16], 1.0, 0.0], [0, 0.66, 1.26, 2.06, 2.58], "-", color=c, alpha=0.25, lw=0.8)
        M = MASSES["v2 sheet (1732 kg)"]
        Ls = [M[15] * shapes[(spec, n)]["free"][15] + M[16] * shapes[(spec, n)]["free"][16] + M[17] for n in runs]
        ok = [shapes[(spec, n)]["free"]["rms17_g"] > 2e-4 for n in runs]
        ax[1].plot(np.array(runs)[ok], np.array(Ls)[ok], "o-", color=c, ms=4, label="{} , L = sum m phi, v2 masses".format(spec))
    ax[0].set_xlabel("phi (free vibration, normalised at 2.06 m)"); ax[0].set_ylabel("height (m)"); ax[0].grid(alpha=0.25)
    ax[0].set_title("Measured first-mode shape, every run (top beam fixed at 2.58 m, base at 0)", loc="left", fontsize=10)
    ax[0].axvline(0, color="k", lw=0.5)
    for lab, M in MASSES.items():
        ax[1].axhline(sum(M.values()), color="grey", ls=":", lw=0.8); ax[1].text(1, sum(M.values()) + 15, "sum m, " + lab, fontsize=7, color="grey")
    ax[1].set_xlabel("run"); ax[1].set_ylabel("L_eff (kg)"); ax[1].grid(alpha=0.25); ax[1].legend(fontsize=8)
    ax[1].set_title("Effective mass for base shear vs d(2.06 m)", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_effective_mass.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
