"""route2_Tp_preflight.py -- Route 2 with the record's predominant period as the carrier of BOTH waves.

    wave A : one sine cycle at T_p(record), amplitude such that Sd_A(T1)    = s * Sd_rec(T1)
    wave B : one sine cycle at T_p(record), amplitude such that Sd_B(T_eff) = s * Sd_rec(T_eff)
    signal : A then B, same phase, DELTA_T 0.005 s, written as 3DEC velocity tables.

Differs from the original Route 2 (route2_full25_LS_3dec.py), where wave A ran AT T1 and wave B AT T_eff, and from
the fixed-point driver (carrier = max(T_p, T*)). T1 / T_eff here come from the ft0 / c 0.2 MPa pushover with the
test's effective mass L = 1200 kg (see exp_envelope_peaks/periods_compared.csv).

T_p = period of the peak 5 %-damped pseudo-acceleration of each record (vel_*.txt): HU12 0.121, EC40 0.251, FR76 0.261 s.

OUTPUT  <out>/route2_Tp_preflight.csv  (per run: targets, amplitudes, PGAs, the combined signal's Sd at T1 and T_eff)
        <out>/RunNN_REC_sX_vel.txt      (3DEC tables, same format as the Route 2 drivers)
        <out>/periods_Tp.json           (T1, T_eff, T_p per record -- read by route2_Tp_3dec.py)
        <out>/fig_route2_Tp_preflight.png
USAGE   python route2_Tp_preflight.py <repo_dir> <out_dir> [--T1 0.1014] [--Teff 0.3572] [--cycles 1]
"""
import sys, os, csv, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DELTA_T, XI = 0.005, 0.05
TP = {"HU12": 0.121, "EC40": 0.251, "FR76": 0.261}
PGA_REC = {"HU12": 0.130, "EC40": 0.375, "FR76": 0.352}
PROTOCOL = [(1, "HU12", 0.50), (2, "HU12", 0.75), (3, "EC40", 0.20), (4, "HU12", 1.00), (5, "HU12", 1.25),
            (6, "EC40", 0.30), (7, "HU12", 1.50), (8, "EC40", 0.40), (9, "HU12", 1.75), (10, "HU12", 2.00),
            (11, "EC40", 0.50), (12, "HU12", 2.25), (13, "HU12", 2.50), (14, "HU12", 2.75), (15, "HU12", 3.00),
            (16, "HU12", 3.50), (17, "HU12", 4.00), (18, "HU12", 4.50), (19, "HU12", 5.00), (20, "HU12", 5.50),
            (21, "HU12", 6.00), (22, "FR76", 1.00), (23, "FR76", 1.50), (24, "FR76", 1.75), (25, "FR76", 2.00)]


def arg(flag, default):
    return type(default)(sys.argv[sys.argv.index(flag) + 1]) if flag in sys.argv else default


def sd_spectrum(ag, dt, Ts, xi=XI):
    """Newmark average-acceleration SDOF, as in the Route 2 drivers."""
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float))
    w = 2.0 * np.pi / Ts; c = 2.0 * xi * w; w2 = w * w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0] - c * v - w2 * u; umax = np.abs(u)
    for i in range(1, len(ag)):
        u_new = u + dt * v + 0.25 * dt * dt * a; v_new = v + 0.5 * dt * a
        a_new = (-ag[i] - c * v_new - w2 * u_new) / (1.0 + 0.5 * dt * c + 0.25 * dt * dt * w2)
        u = u_new + 0.25 * dt * dt * a_new; v = v_new + 0.5 * dt * a_new; a = a_new
        umax = np.maximum(umax, np.abs(u))
    return umax


def grid_period(T):
    spc = max(1, int(round(T / DELTA_T))); return spc * DELTA_T, spc


def cycles(Vamp, T, n):
    Tg, spc = grid_period(T); tt = np.arange(n * spc + 1) * DELTA_T
    return tt, Vamp * np.sin(2.0 * math.pi * tt / Tg), Tg


def calibrate(T_carrier, T_eval, Sd_target, n):
    """velocity amplitude of n cycles at T_carrier whose Sd at T_eval (with a 2 s free tail) equals Sd_target."""
    _, v, Tg = cycles(1.0, T_carrier, n); vt = np.r_[v, np.zeros(int(2.0 / DELTA_T))]
    sd0 = float(sd_spectrum(np.gradient(vt, DELTA_T), DELTA_T, [T_eval])[0])
    return Sd_target / sd0, Tg


def main():
    repo, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
    T1, TE, NC = arg("--T1", 0.1014), arg("--Teff", 0.3572), arg("--cycles", 1)
    SPEC = {}
    for r in TP:
        d = np.genfromtxt(os.path.join(repo, "spectrum_{}.csv".format(r)), delimiter=",", skip_header=1); SPEC[r] = (d[:, 0], d[:, 1])
    sd_rec = lambda r, T: float(np.interp(T, SPEC[r][0], SPEC[r][1]))
    json.dump(dict(T1_s=T1, Teff_s=TE, Tp_s=TP, cycles=NC, xi=XI, delta_t=DELTA_T,
                   source="ft0/c0.2 pushover (pushover_pos_new.csv), K_el 4.60 kN/mm, F_y 10.17 kN, d_u 27.4 mm, L 1200 kg"),
              open(os.path.join(out, "periods_Tp.json"), "w"), indent=1)
    rows = []; sigs = {}
    print("T1 {:.4f}  T_eff {:.4f}  cycles {}  T_p {}".format(T1, TE, NC, TP))
    print("{:>3s} {:4s} {:>5s} | {:>7s} {:>7s} | {:>7s} {:>7s} {:>6s} {:>6s} | {:>8s} {:>8s} | {:>6s} {:>6s}".format(
        "run", "rec", "s", "SdA_tgt", "SdB_tgt", "VA", "VB", "PGA_A", "PGA_B", "Sd(T1)", "Sd(Teff)", "PGArec", "amp"))
    for n, rec, s in PROTOCOL:
        Tp = TP[rec]; SdA, SdB = s * sd_rec(rec, T1), s * sd_rec(rec, TE)
        VA, Tg = calibrate(Tp, T1, SdA, NC); VB, _ = calibrate(Tp, TE, SdB, NC)
        _, vA, _ = cycles(VA, Tp, NC); _, vB, _ = cycles(VB, Tp, NC)
        v = np.r_[vA, vB[1:]]; t = np.arange(len(v)) * DELTA_T
        vt = np.r_[v, np.zeros(int(2.5 / DELTA_T))]
        tot = sd_spectrum(np.gradient(vt, DELTA_T), DELTA_T, [T1, TE])
        pgaA, pgaB = np.abs(np.gradient(vA, DELTA_T)).max() / 9.81, np.abs(np.gradient(vB, DELTA_T)).max() / 9.81
        lbl = "Run{:02d}_{}_s{}".format(n, rec, "{:.2f}".format(s).replace(".", "p"))
        with open(os.path.join(out, lbl + "_vel.txt"), "w", newline="\n") as f:
            f.write("{}\n{}\t0\n".format(lbl, len(t)))
            for a_, b_ in zip(t, v): f.write("{:.6f}\t{:.9e}\n".format(a_, b_))
        row = dict(run=n, record=rec, scale=s, Tp_s=Tg, T1_s=T1, Teff_s=TE, SdA_target_mm=SdA * 1e3, SdB_target_mm=SdB * 1e3,
                   VA_mps=VA, VB_mps=VB, PGA_A_g=pgaA, PGA_B_g=pgaB, PGA_rec_scaled_g=s * PGA_REC[rec],
                   SdT1_total_mm=tot[0] * 1e3, SdTeff_total_mm=tot[1] * 1e3, table_dur_s=t[-1],
                   ratio_SdT1=tot[0] / SdA, ratio_SdTeff=tot[1] / SdB, amp_B_over_rec=pgaB / (s * PGA_REC[rec]))
        rows.append(row); sigs[n] = (t, v)
        print("{:3d} {:4s} {:5.2f} | {:7.2f} {:7.2f} | {:7.4f} {:7.4f} {:6.2f} {:6.2f} | {:8.2f} {:8.2f} | {:6.2f} {:6.1f}x".format(
            n, rec, s, SdA * 1e3, SdB * 1e3, VA, VB, pgaA, pgaB, tot[0] * 1e3, tot[1] * 1e3, s * PGA_REC[rec], row["amp_B_over_rec"]))
    with open(os.path.join(out, "route2_Tp_preflight.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); [w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}) for r in rows]
    # figure
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    for n in (5, 15, 21, 23, 25):
        t, v = sigs[n]; ax[0].plot(t, v, lw=1, label="run {}".format(n))
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("table velocity (m/s)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)
    ax[0].set_title("A + B at T_p (one cycle each)", loc="left", fontsize=10)
    rn = [r["run"] for r in rows]
    ax[1].plot(rn, [r["PGA_A_g"] for r in rows], "o-", ms=4, label="PGA wave A"); ax[1].plot(rn, [r["PGA_B_g"] for r in rows], "s-", ms=4, label="PGA wave B")
    ax[1].plot(rn, [r["PGA_rec_scaled_g"] for r in rows], "k--", lw=1, label="record PGA x s"); ax[1].axhline(0.9, color="grey", ls=":", lw=0.8)
    ax[1].set_xlabel("run"); ax[1].set_ylabel("g"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25); ax[1].set_title("Wave PGAs vs the record's", loc="left", fontsize=10)
    ax[2].plot(rn, [r["ratio_SdT1"] for r in rows], "o-", ms=4, label="Sd(T1) of A+B / target"); ax[2].plot(rn, [r["ratio_SdTeff"] for r in rows], "s-", ms=4, label="Sd(T_eff) of A+B / target")
    ax[2].axhline(1, color="k", lw=0.8); ax[2].set_xlabel("run"); ax[2].legend(fontsize=8); ax[2].grid(alpha=0.25)
    ax[2].set_title("Interaction of the two waves (calibrated separately)", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_route2_Tp_preflight.png"), dpi=150)
    print("->", out)


if __name__ == "__main__":
    main()
