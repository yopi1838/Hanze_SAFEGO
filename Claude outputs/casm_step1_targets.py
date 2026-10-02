"""casm_step1_targets.py -- capacity-adaptive spectral matching, STEP 1 (no 3DEC).

For every run of the 25-run protocol: solve the capacity-spectrum fixed point
    d* = scale * Sd_record( T_sec(d*) ),   T_sec(d) = 2 pi sqrt( m_eff / (F(d)/d) )
with F(d) the pushover backbone (pushover_pos.csv, |F_tot_kN| vs d_ctrl_mm) and
Sd_record the 5%-damped displacement spectrum of the record (spectrum_*.csv).
Reports d*, T*, and the current Route 2 targets (Sd at T1 and at T_eff) for comparison,
plus the model/record/test peaks we know, so the targets can be judged before any run.

Also: record mean periods T_m (Rathje), PGA, PGV of HU12 / FR76 for the carrier/bounds.
Outputs: casm_step1_targets.csv, fig_casm_step1.png
"""
import json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROTOCOL = [(1,"HU12",0.50),(2,"HU12",0.75),(3,"EC40",0.20),(4,"HU12",1.00),(5,"HU12",1.25),(6,"EC40",0.30),
            (7,"HU12",1.50),(8,"EC40",0.40),(9,"HU12",1.75),(10,"HU12",2.00),(11,"EC40",0.50),(12,"HU12",2.25),
            (13,"HU12",2.50),(14,"HU12",2.75),(15,"HU12",3.00),(16,"HU12",3.50),(17,"HU12",4.00),(18,"HU12",4.50),
            (19,"HU12",5.00),(20,"HU12",5.50),(21,"HU12",6.00),(22,"FR76",1.00),(23,"FR76",1.50),(24,"FR76",1.75),(25,"FR76",2.00)]
# what we know (mm, peak at 2.06 m, within-run).  Exact where the postprocess printed them; others read off the plot (~).
KNOWN = {
    "record on model":   {21: 17.0, 22: 3.7, 23: 5.9, 24: 8.3, 25: 41.6},
    "Route 2 (as run)":  {21: 21.5, 22: 8.1, 23: 28.6, 24: 32.1, 25: 38.6},
    "US-1 (Test 9)":     {21: 17.0, 22: 4.0, 23: 13.15, 24: 29.48},
    "US-2 (Test 12)":    {21: 10.3, 22: 4.0, 23: 5.9, 24: 7.7, 25: 31.0},
}
APPROX_RUNS = {"record on model": [], "Route 2 (as run)": [], "US-1 (Test 9)": [21, 22], "US-2 (Test 12)": [22, 23, 24, 25]}

bj = json.load(open("bilinear_periods.json")); M_EFF = bj["rule"]["m_eff"]; M_ALT = bj["model"]["m_alt"]
T1, TEFF = bj["model"]["T1_s"], bj["model"]["Teff_s"]; DY, DU = bj["model"]["per_branch"]["pos"]["dy_mm"], bj["model"]["per_branch"]["pos"]["du_mm"]

# ---- backbone: |F_tot| vs d, monotone in d (asc then desc branch), linear interpolation, through the origin
po = np.genfromtxt("pushover_pos.csv", delimiter=",", skip_header=1, usecols=(2, 5))
d_bb = np.r_[0.0, po[:, 0]]; F_bb = np.r_[0.0, np.abs(po[:, 1])]
o = np.argsort(d_bb); d_bb, F_bb = d_bb[o], F_bb[o]
def F_of(d):  return float(np.interp(d, d_bb, F_bb))
def T_sec(d, m=M_EFF):
    k = F_of(d) * 1e3 / (d * 1e-3)          # N/m
    return 2 * math.pi * math.sqrt(m / k)
print(f"backbone: K1 at 0.5 mm {F_of(0.5)/0.5:.2f} kN/mm -> T_sec {T_sec(0.5):.4f} s (bilinear T1 {T1:.4f});  "
      f"T_sec at dy {T_sec(DY):.3f}, at 10 mm {T_sec(10):.3f}, at du {T_sec(DU):.3f} (bilinear T_eff {TEFF:.3f}), at 30 mm {T_sec(30):.3f}")

SPEC = {k: np.genfromtxt(f"spectrum_{k}.csv", delimiter=",", skip_header=1) for k in ("HU12", "EC40", "FR76")}
def Sd(rec, T): s = SPEC[rec]; return float(np.interp(T, s[:, 0], s[:, 1])) * 1e3     # mm

def fixed_point(rec, scale, m=M_EFF, dmax=70.0):
    """all crossings of g(d) = scale*Sd(T_sec(d)) - d on a log grid; return first (stable) and the list."""
    dd = np.logspace(math.log10(0.02), math.log10(dmax), 1200)
    g = np.array([scale * Sd(rec, T_sec(d, m)) - d for d in dd])
    xs = []
    for i in range(len(dd) - 1):
        if g[i] > 0 >= g[i + 1]:                   # demand curve crosses the 1:1 line from above -> stable
            xs.append(float(np.interp(0.0, [g[i + 1], g[i]], [dd[i + 1], dd[i]])))
    if not xs:                                     # demand above d everywhere up to dmax: runaway (no linearised solution)
        return None, []
    return xs[0], xs


# ---- damage memory (property-based): once the method's own d* has exceeded dy, the initial branch is
# softened by (T1/T_end)^2 with T_end/T1 = CRACK_T_RATIO (the ring-down elongation Strategy C measured
# after cracking, 1.30). The rocking plateau and descent are geometry + overburden and stay as they are.
CRACK_T_RATIO = 1.30
K1 = F_of(0.5) / 0.5
def F_state(d, cracked):
    return min(K1 * d / CRACK_T_RATIO ** 2, F_of(d)) if cracked else F_of(d)
def T_state(d, cracked, m=M_EFF):
    return 2 * math.pi * math.sqrt(m / (F_state(d, cracked) * 1e3 / (d * 1e-3)))
def crossings(rec, scale, cracked, m=M_EFF, dmax=70.0):
    """stable solutions (demand crosses the 1:1 line from above), tipping points (from below),
    runaway = demand still above d at dmax."""
    dd = np.logspace(math.log10(0.02), math.log10(dmax), 1500)
    g = np.array([scale * Sd(rec, T_state(d, cracked, m)) - d for d in dd])
    st = [float(np.interp(0, [g[i+1], g[i]], [dd[i+1], dd[i]])) for i in range(len(dd)-1) if g[i] > 0 >= g[i+1]]
    un = [float(dd[i]) for i in range(len(dd)-1) if g[i] <= 0 < g[i+1]]
    return st, un, bool(g[-1] > 0)

rows = []; cracked = False
for run, rec, sc in PROTOCOL:
    st, un, ra = crossings(rec, sc, cracked); st_alt, _, _ = crossings(rec, sc, cracked, M_ALT)
    ds = st[0] if st else None; Ts = T_state(ds, cracked) if ds else float("nan")
    rows.append(dict(run=run, record=rec, scale=sc, backbone="cracked" if cracked else "virgin",
                     d_star_mm=ds, T_star_s=Ts, tipping_mm=(un[0] if un else None), runaway_beyond_tipping=ra,
                     d_star_m1318_mm=(st_alt[0] if st_alt else None),
                     Sd_at_T1_mm=sc * Sd(rec, T1), Sd_at_Teff_mm=sc * Sd(rec, TEFF),
                     state="elastic" if (ds or 1e9) < DY else ("cracked" if (ds or 1e9) < DU else "beyond du")))
    if ds and ds >= DY: cracked = True
print(f"\n{'run':>3} {'rec':>5} {'sc':>5} {'backbone':>8} | {'d*':>6} {'T*':>6} {'tipping':>7} {'runaway':>7} | {'R2 Sd@T1':>8} {'R2 Sd@Teff':>10} | {'rec/model':>9} {'US-1':>6} {'US-2':>6}")
for r in rows:
    k = r["run"]; tp = f"{r['tipping_mm']:.1f}" if r["tipping_mm"] else "--"
    print(f"{k:3d} {r['record']:>5} {r['scale']:5.2f} {r['backbone']:>8} | {r['d_star_mm'] if r['d_star_mm'] else float('nan'):6.2f} {r['T_star_s']:6.3f} {tp:>7} {str(r['runaway_beyond_tipping']):>7} | "
          f"{r['Sd_at_T1_mm']:8.2f} {r['Sd_at_Teff_mm']:10.2f} | {KNOWN['record on model'].get(k, float('nan')):9.1f} {KNOWN['US-1 (Test 9)'].get(k, float('nan')):6.1f} {KNOWN['US-2 (Test 12)'].get(k, float('nan')):6.1f}")
# ---- record properties for the carrier / bounds
def rec_props(fn):
    d = np.loadtxt(fn, skiprows=2); t, v = d[:, 0], d[:, 1]; dt = np.median(np.diff(t)); a = np.gradient(v, dt)
    f = np.fft.rfftfreq(len(a), dt); A = np.abs(np.fft.rfft(a)); m = (f > 0.25) & (f < 20)
    Tm = np.sum(A[m] ** 2 / f[m]) / np.sum(A[m] ** 2)
    return dict(PGV=abs(v).max(), PGA_g=abs(a).max() / 9.81, Tm=Tm, dur=t[-1])
for k, fn in (("HU12", "vel_HU.txt"), ("FR76", "vel_FR.txt")):
    p = rec_props(fn); print(f"{k}: PGV {p['PGV']:.3f} m/s  PGA {p['PGA_g']:.3f} g  T_m {p['Tm']:.3f} s  dur {p['dur']:.1f} s")

with open("casm_step1_targets.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---- figure: targets vs known peaks, and the backbone / secant-period map
F = 1.25
fig, ax = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw=dict(width_ratios=[1.6, 1]))
a = ax[0]; runs = [r["run"] for r in rows]
a.plot(runs, [r["Sd_at_Teff_mm"] for r in rows], "-", color="#b4491a", lw=1.2, alpha=0.6, label="Route 2 target: Sd at T_eff (cracked state, every run)")
a.plot(runs, [r["Sd_at_T1_mm"] for r in rows], "--", color="#b4491a", lw=1.0, alpha=0.6, label="Route 2 target: Sd at T1")
a.plot(runs, [r["d_star_mm"] if r["d_star_mm"] else np.nan for r in rows], "o-", color="#18212b", lw=2, ms=6, label="capacity-adaptive target d* (fixed point)")
cols = {"record on model": "#4a7c59", "Route 2 (as run)": "#b4491a", "US-1 (Test 9)": "#2b5f8a", "US-2 (Test 12)": "#f0a020"}
mk = {"record on model": "s", "Route 2 (as run)": "^", "US-1 (Test 9)": "o", "US-2 (Test 12)": "D"}
for lab, vals in KNOWN.items():
    a.plot(list(vals), list(vals.values()), mk[lab], color=cols[lab], ms=8, mfc="none" if lab.startswith("US") else cols[lab], label=lab + " (peak, known runs)")
tip = [(r["run"], r["tipping_mm"]) for r in rows if r["tipping_mm"]]
a.plot([t[0] for t in tip], [t[1] for t in tip], "_", color="#8a8f97", ms=14, mew=2, label="tipping displacement (demand outruns capacity above it)")
for r in rows:
    if r["runaway_beyond_tipping"]:
        a.annotate("runaway above", (r["run"], r["tipping_mm"]), xytext=(0, 26 if r["run"] % 2 else 12), textcoords="offset points", ha="center", fontsize=7.5 * F, color="#8a8f97")
a.axhline(DY, color="#8a8f97", ls=":", lw=1); a.text(12.5, DY * 1.12, f"dy {DY:.1f}", fontsize=9 * F, color="#8a8f97")
a.axhline(DU, color="#8a8f97", ls=":", lw=1); a.text(12.5, DU * 1.12, f"du {DU:.1f}", fontsize=9 * F, color="#8a8f97")
a.set_yscale("log"); a.set_ylim(0.05, 80); a.set_xticks(runs); a.set_xlabel("run", fontsize=11 * F); a.set_ylabel("peak OOP displacement at 2.06 m (mm)", fontsize=11 * F)
a.set_title("(a) what each run is asked to reach", loc="left", fontweight="bold", fontsize=12 * F); a.legend(fontsize=8 * F, loc="lower right", frameon=True, framealpha=0.9)
a.axvspan(21.5, 25.5, color="#f2efe8", zorder=0)
b = ax[1]
dd = np.logspace(-1.3, math.log10(70), 300)
b.plot(dd, [T_state(x, False) for x in dd], color="#18212b", lw=2, label=f"virgin backbone, m_eff {M_EFF:.0f} kg")
b.plot(dd, [T_state(x, True) for x in dd], color="#4a7c59", lw=2, label=f"cracked: initial stiffness /{CRACK_T_RATIO}^2, plateau unchanged")
b.plot(dd, [T_state(x, False, M_ALT) for x in dd], color="#8a8f97", lw=1.2, ls="--", label=f"virgin, m = {M_ALT:.0f} kg")
b.axhline(T1, color="#b4491a", ls=":", lw=1); b.text(0.06, T1 * 1.03, "T1 (bilinear)", fontsize=9 * F, color="#b4491a")
b.axhline(TEFF, color="#b4491a", ls=":", lw=1); b.text(0.06, TEFF * 1.03, "T_eff (bilinear)", fontsize=9 * F, color="#b4491a")
for r in rows:
    if r["d_star_mm"]: b.plot(r["d_star_mm"], r["T_star_s"], "o", color="#4a7c59" if r["record"] != "FR76" else "#b4491a", ms=5)
b.set_xscale("log"); b.set_xlabel("displacement d (mm)", fontsize=11 * F); b.set_ylabel("secant period (s)", fontsize=11 * F); b.set_ylim(0.08, 0.5)
b.set_title("(b) the wall's period follows its amplitude", loc="left", fontweight="bold", fontsize=12 * F); b.legend(fontsize=8.5 * F, frameon=False, loc="upper left")
for x in ax:
    x.grid(alpha=0.25, which="both"); x.tick_params(labelsize=9.5 * F)
    for sp in ("top", "right"): x.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("fig_casm_step1.png", dpi=150); print("-> fig_casm_step1.png, casm_step1_targets.csv")
