"""Prepare the controlled 2 x 2 Route 2 experiment, OUTSIDE 3DEC.

python route2_comparison/prepare_signals.py
Only wave B's carrier and spectral target location vary. Wave A stays at T1.
The adaptive target is a secant-spectrum hypothesis, not standard N2 or a
dynamic stability criterion. No DEM response is used to fit these signals.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import lfilter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "prepared"
DT = 0.001                 # velocity table spacing, not the DEM integration step
TAIL = 2.5                 # included in spectral calibration
XI = 0.05
CRACK_T_RATIO = 1.30        # inherited hypothesis; not newly fitted
TP_RANGE = (0.05, 1.0)     # fixed for every record; peak of pseudo-acceleration
FIT_TOL = 0.05             # BOTH target ordinates must be within 5%
METHODS = {
    "M1_fixed_wall": ("fixed", "wall"),
    "M2_fixed_record": ("fixed", "record"),
    "M3_adaptive_wall": ("adaptive", "wall"),
    "M4_adaptive_record": ("adaptive", "record"),
}
COLORS = ["#536878", "#ce7626", "#277d76", "#8756a3"]
PROTOCOL = [(1,"HU12",.5),(2,"HU12",.75),(3,"EC40",.2),(4,"HU12",1),
    (5,"HU12",1.25),(6,"EC40",.3),(7,"HU12",1.5),(8,"EC40",.4),
    (9,"HU12",1.75),(10,"HU12",2),(11,"EC40",.5),(12,"HU12",2.25),
    (13,"HU12",2.5),(14,"HU12",2.75),(15,"HU12",3),(16,"HU12",3.5),
    (17,"HU12",4),(18,"HU12",4.5),(19,"HU12",5),(20,"HU12",5.5),
    (21,"HU12",6),(22,"FR76",1),(23,"FR76",1.5),(24,"FR76",1.75),(25,"FR76",2)]


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def oscillator_history(acc_intervals, dt, period, xi=XI):
    """Exact state transition for constant acceleration in each time interval.

    This is the derivative of the linearly interpolated velocity TABLE that is
    applied to 3DEC. Returned displacements are at interval boundaries; halving
    dt checks the error from locating a peak between those boundaries.
    """
    w = 2 * np.pi / period
    A = np.array([[0., 1.], [-w*w, -2*xi*w]])
    aug = np.zeros((3, 3))
    aug[:2, :2] = A
    aug[1, 2] = -1.
    E = expm(aug * dt)
    Ad, Bd = E[:2, :2], E[:2, 2]
    b = [Bd[0], Ad[0, 1]*Bd[1] - Ad[1, 1]*Bd[0]]
    a = [1., -np.trace(Ad), np.linalg.det(Ad)]
    return np.r_[0., lfilter(b, a, acc_intervals)]


def response_histories(v, dt, periods, tail=TAIL):
    a = np.r_[np.diff(v)/dt, np.zeros(int(round(tail/dt)))]
    return np.array([oscillator_history(a, dt, T) for T in periods])


def spectrum(v, dt, periods, tail=TAIL):
    return np.max(np.abs(response_histories(v, dt, periods, tail)), axis=1)


def basis(carrier_a, carrier_b, dt=DT):
    """One positive-first velocity cycle each, sequential, continuous velocity."""
    na, nb = [max(4, int(round(T/dt))) for T in (carrier_a, carrier_b)]
    va, vb = np.zeros(na+nb+1), np.zeros(na+nb+1)
    va[:na+1] = np.sin(2*np.pi*np.arange(na+1)/na)
    vb[na:] = np.sin(2*np.pi*np.arange(nb+1)/nb)
    va[[0, na, -1]] = 0.
    vb[[0, na, -1]] = 0.
    return np.arange(len(va))*dt, va, vb, na*dt, nb*dt


def fit_amplitudes(va, vb, periods, targets):
    """Global grid plus bounded refinements in mixing ratio; scale is analytic.

    Linear oscillator HISTORIES superpose; spectral PEAKS do not. Search all
    mixing ratios, including a vanished wave. An infeasible two-target fit is
    retained as a diagnostic and explicitly marked, never described as matched.
    """
    ha = response_histories(va, DT, periods)
    hb = response_histories(vb, DT, periods)
    unit_a = targets[0] / np.max(np.abs(ha[0]))
    unit_b = targets[1] / np.max(np.abs(hb[1]))
    ha, hb = ha*unit_a, hb*unit_b

    def evaluate(r):
        s = np.max(np.abs(r*ha + (1-r)*hb), axis=1)
        gain = np.exp(np.mean(np.log(targets/s)))
        return float(np.mean(np.log(gain*s/targets)**2)), gain

    grid = np.linspace(0., 1., 1001)
    cost = np.array([evaluate(r)[0] for r in grid])
    candidates = [(cost[0], 0.), (cost[-1], 1.)]
    for i in range(1, len(grid)-1):
        if cost[i] <= cost[i-1] and cost[i] <= cost[i+1]:
            opt = minimize_scalar(lambda r: evaluate(r)[0],
                bounds=(grid[i-1], grid[i+1]), method="bounded",
                options={"xatol": 1e-12})
            candidates.append((float(opt.fun), float(opt.x)))
    _, r = min(candidates)
    _, gain = evaluate(r)
    amps = np.array([gain*r*unit_a, gain*(1-r)*unit_b])
    achieved = np.max(np.abs(amps[0]*ha/unit_a + amps[1]*hb/unit_b), axis=1)
    return amps, achieved


def inputs():
    bj = json.loads((ROOT/"bilinear_idealisation/bilinear_periods.json").read_text())
    bb = np.genfromtxt(ROOT/"pushover_results/pushover_pos.csv", delimiter=",", names=True)
    x, F = bb["d_ctrl_mm"], np.abs(bb["F_tot_kN"])
    keep = np.isfinite(x) & np.isfinite(F) & (x > 0)
    x, F = x[keep], F[keep]
    order = np.argsort(x)
    x, F = np.r_[0., x[order]], np.r_[0., F[order]]
    if np.any(np.diff(x) <= 0) or x[-1] < 70:
        raise ValueError("Backbone must have unique displacements covering 0..70 mm")
    specs = {r: np.loadtxt(ROOT/("spectrum_"+r+".csv"), delimiter=",", skiprows=1)
             for r in ("HU12", "EC40", "FR76")}
    return bj, x, F, specs


def sd(spec, periods):
    if np.any(np.asarray(periods) < spec[0, 0]) or np.any(np.asarray(periods) > spec[-1, 0]):
        raise ValueError("Target period outside supplied spectrum; no extrapolation")
    return np.interp(periods, spec[:, 0], spec[:, 1])


def adaptive_targets(bj, x, F, specs):
    """Freeze one common proxy-state schedule before running any of the methods."""
    mass = bj["rule"]["m_eff"]
    dy = bj["model"]["per_branch"]["pos"]["dy_mm"]
    k = np.interp(.5, x, F)/.5
    cracked = False
    rows = []
    for run, rec, scale in PROTOCOL:
        def period(d):
            force = np.interp(d, x, F)
            if cracked:
                force = np.minimum(force, k*np.asarray(d)/CRACK_T_RATIO**2)
            return 2*np.pi*np.sqrt(mass/(force*1e6/np.asarray(d)))

        def g(d):
            return scale*sd(specs[rec], period(d))*1000-d

        grid = np.geomspace(.02, 70., 1500)
        gg = g(grid)
        down, up = [], []
        for i in range(len(grid)-1):
            if gg[i] > 0 >= gg[i+1]:
                down.append(brentq(g, grid[i], grid[i+1]))
            elif gg[i] <= 0 < gg[i+1]:
                up.append(brentq(g, grid[i], grid[i+1]))
        if not down:
            raise ValueError("No first downward crossing for run %d; do not invent a target" % run)
        d = down[0]
        rows.append(dict(run=run, record=rec, scale=scale,
            proxy_backbone="cracked" if cracked else "virgin", d_star_mm=d,
            T_star_s=float(period(d)), upper_crossing_mm=up[0] if up else None,
            demand_above_at_70mm=bool(gg[-1] > 0), crack_period_ratio=CRACK_T_RATIO))
        cracked = cracked or d >= dy
    return rows


def record_properties(specs):
    rows = []
    for rec, filename in [("HU12", "vel_HU.txt"),("EC40", "vel_EC.txt"),("FR76", "vel_FR.txt")]:
        spec = specs[rec]
        mask = (spec[:, 0] >= TP_RANGE[0]) & (spec[:, 0] <= TP_RANGE[1])
        b = spec[mask]
        tp = b[np.argmax(b[:, 1]*(2*np.pi/b[:, 0])**2), 0]
        raw = np.loadtxt(ROOT/filename, skiprows=2)
        dt = np.median(np.diff(raw[:, 0]))
        acc = np.diff(raw[:, 1])/np.diff(raw[:, 0])
        f = np.fft.rfftfreq(len(acc), dt)
        power = np.abs(np.fft.rfft(acc))**2
        band = (f > .25) & (f < 20.)
        tm = np.sum(power[band]/f[band])/np.sum(power[band])
        rows.append(dict(record=rec, Tp_s=float(tp), Tm_s=float(tm),
            PGA_g=float(np.max(np.abs(acc))/9.81), PGV_mps=float(np.max(np.abs(raw[:, 1]))),
            duration_s=float(raw[-1, 0]-raw[0, 0]), table=filename))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Never change inputs underneath an already started DEM comparison.
    if list((OUT.parent/"results").glob("*/checkpoint.json")):
        raise RuntimeError("DEM checkpoints exist; use a new comparison directory to regenerate")
    bj, x, F, specs = inputs()
    targets = adaptive_targets(bj, x, F, specs)
    props = record_properties(specs)
    properties = {p["record"]: p for p in props}
    t1, teff = bj["model"]["T1_s"], bj["model"]["Teff_s"]
    write_csv(OUT/"adaptive_targets.csv", targets)
    write_csv(OUT/"record_properties.csv", props)
    rows, waves = [], {}
    for method, (target_rule, carrier_rule) in METHODS.items():
        (OUT/method).mkdir(exist_ok=True)
        for tr in targets[21:]:
            run, rec, scale = tr["run"], tr["record"], tr["scale"]
            tb = teff if target_rule == "fixed" else tr["T_star_s"]
            cb = tb if carrier_rule == "wall" else properties[rec]["Tp_s"]
            periods = np.array([t1, tb])
            demand = scale*sd(specs[rec], periods)
            t, va, vb, ca_actual, cb_actual = basis(t1, cb)
            amps, got = fit_amplitudes(va, vb, periods, demand)
            v = amps[0]*va + amps[1]*vb
            error = np.abs(got/demand-1)
            # Halve response integration spacing without changing the applied table.
            tfine = np.arange(2*(len(t)-1)+1)*(DT/2)
            vf = np.interp(tfine, t, v)
            refined = spectrum(vf, DT/2, periods)
            numerical_error = float(np.max(np.abs(refined/got-1)))
            if numerical_error > .005:
                raise RuntimeError("SDOF peak sampling not converged to 0.5%")
            label = "Run%02d_%s_s%s" % (run, rec, ("%.2f" % scale).replace(".", "p"))
            table = OUT/method/(label+"_vel.txt")
            with table.open("w", newline="\n") as f:
                f.write(label+"\n%d\t0\n" % len(t))
                np.savetxt(f, np.c_[t, v], fmt=["%.6f", "%.12e"], delimiter="\t")
            pga = float(np.max(np.abs(np.diff(v)/DT))/9.81)
            pgv = float(np.max(np.abs(v)))
            rows.append(dict(method=method, run=run, record=rec, scale=scale,
                target_rule=target_rule, carrier_rule=carrier_rule,
                T_A_target_s=t1, T_B_target_s=tb, T_A_carrier_s=ca_actual, T_B_carrier_s=cb_actual,
                target_A_mm=float(demand[0]*1000), target_B_mm=float(demand[1]*1000),
                achieved_A_mm=float(got[0]*1000), achieved_B_mm=float(got[1]*1000),
                amplitude_A_mps=float(amps[0]), amplitude_B_mps=float(amps[1]),
                max_target_error_pct=float(max(error)*100), fit_ok=bool(max(error) <= FIT_TOL),
                numerical_error_pct=numerical_error*100,
                PGA_g=pga, PGV_mps=pgv, PGA_ratio_to_record=pga/(scale*properties[rec]["PGA_g"]),
                PGV_ratio_to_record=pgv/(scale*properties[rec]["PGV_mps"]),
                pulse_duration_s=float(t[-1]), tail_s=TAIL,
                net_table_displacement_mm=float(np.trapz(v, t)*1000),
                table=table.relative_to(OUT.parent).as_posix(), sha256=sha(table)))
            waves[(method, run)] = (t, v, periods, demand)
    write_csv(OUT/"signal_comparison.csv", rows)
    source_files = ["bilinear_idealisation/bilinear_periods.json", "pushover_results/pushover_pos.csv"]
    source_files += ["spectrum_"+r+".csv" for r in specs] + [p["table"] for p in props]
    manifest = dict(schema=1, runs=[22,23,24,25], methods=list(METHODS), dt_s=DT,
        xi=XI, tail_s=TAIL, n_cycles_A=1, n_cycles_B=1, phase="same, positive velocity first",
        record_carrier="peak pseudo-acceleration period of supplied Sd spectrum",
        predominant_period_range_s=list(TP_RANGE), fit_tolerance_fraction=FIT_TOL,
        adaptive_memory="common precomputed proxy schedule; independent of DEM results",
        crack_period_ratio=CRACK_T_RATIO, automatic_tipping_switch=False,
        source_hashes={f:sha(ROOT/f) for f in source_files}, signals=rows)
    write_json(OUT/"manifest.json", manifest)

    # Representative run: target/carrier structure is identical for FR76 22..25,
    # with amplitudes proportional to the intensity for this frozen proxy state.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    dense = np.geomspace(.05, .6, 180)
    for ax, (method, color) in zip(axes.flat, zip(METHODS, COLORS)):
        t, v, periods, demand = waves[(method, 22)]
        syn = spectrum(v, DT, dense)*1000
        ax.plot(dense, sd(specs["FR76"], dense)*1000, color="black", lw=1.5, label="FR76 record spectrum")
        ax.plot(dense, syn, color=color, lw=2, label="Combined two-wave spectrum")
        ax.scatter(periods, demand*1000, color="black", marker="x", s=70, label="Two targets")
        ax.set(title=method.replace("_", " "), xlabel="Oscillator period (s)", ylabel="5% spectral displacement (mm)")
        ax.set_xlim(.05, .6); ax.set_ylim(0, 50); ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Route 2 input comparison | FR76 x1.00 | spectra, not DEM predictions")
    fig.savefig(OUT/"spectra_run22.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(4, 1, figsize=(11, 8), sharex=True, layout="constrained")
    for ax, method, color in zip(axes, METHODS, COLORS):
        t, v, _, _ = waves[(method, 22)]
        row = next(r for r in rows if r["method"]==method and r["run"]==22)
        ax.plot(t, v, color=color, lw=1.8)
        ax.axhline(0, color="gray", lw=.5)
        ax.set_ylabel("Velocity (m/s)")
        ax.set_title("%s | B carrier %.3f s | max target error %.1f%%" %
            (method.replace("_", " "), row["T_B_carrier_s"], row["max_target_error_pct"]), loc="left", fontsize=10)
        ax.grid(alpha=.2)
    axes[-1].set_xlabel("Time (s); the 2.5 s zero-input calibration tail is omitted here")
    fig.suptitle("Prepared inputs | FR76 x1.00 | wave A remains at T1")
    fig.savefig(OUT/"waveforms_run22.png", dpi=180)
    plt.close(fig)
    print("Prepared 16 signals. No 3DEC simulations have been run.")
    for r in rows:
        if r["run"] == 22:
            print("%s: B target %.3f s / %.2f mm; carrier %.3f s; fit error %.2f%%; PGA ratio %.2f" %
                (r["method"], r["T_B_target_s"], r["target_B_mm"], r["T_B_carrier_s"],
                 r["max_target_error_pct"], r["PGA_ratio_to_record"]))


if __name__ == "__main__":
    main()
