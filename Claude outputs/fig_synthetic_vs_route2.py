"""fig_synthetic_vs_route2.py -- a spectrum-compatible synthetic for FR76 x1.00 (option 3), next to
the record and the Route 2 signal. No structure period enters the synthetic: the target is the record's
5%-damped Sd over 0.05-1.0 s, the envelope is the record's Arias build-up, the phases are random."""
import numpy as np, math, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DT = 0.005; XI = 0.05; SEED = 1
N_SINES = 14; T_LO, T_HI = 0.05, 1.0
T1, TEFF = 0.112, 0.241

def sd_spectrum(ag, dt, Ts, xi=XI):
    ag = np.asarray(ag, float); Ts = np.atleast_1d(np.asarray(Ts, float)); w = 2*np.pi/Ts; c = 2*xi*w; w2 = w*w
    u = np.zeros_like(w); v = np.zeros_like(w); a = -ag[0]-c*v-w2*u; umax = np.abs(u)
    for i in range(1, len(ag)):
        un = u+dt*v+0.25*dt*dt*a; vn = v+0.5*dt*a; an = (-ag[i]-c*vn-w2*un)/(1+0.5*dt*c+0.25*dt*dt*w2)
        u = un+0.25*dt*dt*an; v = vn+0.5*dt*an; a = an; umax = np.maximum(umax, np.abs(u))
    return umax

# ---- the record
d = np.loadtxt("vel_FR.txt", skiprows=2); t_r, v_r = d[:, 0], d[:, 1]; dt_r = float(np.median(np.diff(t_r))); a_r = np.gradient(v_r, dt_r)
Ts = np.r_[np.arange(0.05, 0.30, 0.005), np.arange(0.30, 1.01, 0.01)]
Sd_r = sd_spectrum(a_r, dt_r, Ts) * 1e3
# Arias build-up -> envelope and significant duration
IA = np.cumsum(a_r**2) * dt_r; IA /= IA[-1]
t5, t95 = t_r[np.searchsorted(IA, 0.05)], t_r[np.searchsorted(IA, 0.95)]
rate = np.gradient(IA, dt_r); k = int(0.5 / dt_r); rate_s = np.convolve(rate, np.ones(k)/k, mode="same")
print(f"FR76: PGV {abs(v_r).max():.3f} m/s  PGA {abs(a_r).max()/9.81:.3f} g  D5-95 {t5:.2f}-{t95:.2f} s = {t95-t5:.2f} s")

# ---- the synthetic: velocity = envelope * sum of sines, amplitudes iterated to the record's Sd at the N periods
t_lo, t_hi = max(0, t5 - 0.5), min(t_r[-1], t95 + 0.5)
tt = np.arange(0, t_hi - t_lo, DT); env = np.interp(tt + t_lo, t_r, np.sqrt(np.clip(rate_s, 0, None)))
env /= env.max(); env[:int(0.2/DT)] *= np.linspace(0, 1, int(0.2/DT)); env[-int(0.3/DT):] *= np.linspace(1, 0, int(0.3/DT))
Tn = np.logspace(math.log10(T_LO), math.log10(T_HI), N_SINES); wn = 2*np.pi/Tn
rng = np.random.default_rng(SEED); ph = rng.uniform(0, 2*np.pi, N_SINES)
Sd_tgt = np.interp(Tn, Ts, Sd_r)
Vn = np.full(N_SINES, 0.02)
def build(Vn):
    v = env * np.sum([Vn[i]*np.sin(wn[i]*tt + ph[i]) for i in range(N_SINES)], axis=0)
    v -= np.linspace(0, 1, len(v)) * v[-1]            # end at zero velocity (no residual table drift)
    return v
for it in range(120):
    v_s = build(Vn); S = sd_spectrum(np.gradient(v_s, DT), DT, Tn) * 1e3
    Vn *= (Sd_tgt / S) ** 0.4
v_s = build(Vn); a_s = np.gradient(v_s, DT); Sd_s = sd_spectrum(a_s, DT, Ts) * 1e3
err = sd_spectrum(a_s, DT, Tn) * 1e3 / Sd_tgt
print(f"synthetic: {N_SINES} sines over {T_LO}-{T_HI} s, {tt[-1]:.2f} s long; Sd match at the {N_SINES} periods: {err.min():.2f}-{err.max():.2f} of target")
print(f"           PGV {abs(v_s).max():.3f} m/s ({abs(v_s).max()/abs(v_r).max():.2f} x record)  PGA {abs(a_s).max()/9.81:.3f} g ({abs(a_s).max()/abs(a_r).max():.2f} x record)")

# ---- Route 2 as run (run 22 from the log): one cycle at T1 (Vamp 0.04879) + one at T_eff (0.16973)
def cyc(V, T):
    spc = int(round(T/DT)); x = np.arange(spc+1)*DT; return V*np.sin(2*np.pi*x/(spc*DT))
v_2 = np.r_[cyc(0.04879, 0.11), cyc(0.16973, 0.24)[1:]]; a_2 = np.gradient(v_2, DT); Sd_2 = sd_spectrum(a_2, DT, Ts) * 1e3
print(f"Route 2 run 22: {len(v_2)*DT:.2f} s, PGV {abs(v_2).max():.3f}, PGA {abs(a_2).max()/9.81:.2f} g")
# Sd at the wall's periods
for lab, S in (("record", Sd_r), ("synthetic", Sd_s), ("Route 2", Sd_2)):
    print(f"  {lab:10s} Sd@T1 {np.interp(T1, Ts, S):5.2f}  @T_eff {np.interp(TEFF, Ts, S):5.2f}  @0.4 {np.interp(0.4, Ts, S):5.1f}  @0.6 {np.interp(0.6, Ts, S):5.1f} mm")

with open("synthetic_FR76_x1p00_vel.txt", "w", newline="\n") as f:
    f.write("synthetic_FR76_x1p00\n{}\t0\n".format(len(tt)))
    for a_, b_ in zip(tt, v_s): f.write("{:.6f}\t{:.9e}\n".format(a_, b_))

# ---- figure
F = 1.25
fig = plt.figure(figsize=(17, 9.5)); gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1], height_ratios=[1, 1, 1])
C = {"record": "#18212b", "synthetic": "#2b5f8a", "route2": "#b4491a"}
ax0 = fig.add_subplot(gs[0, 0]); ax0.plot(t_r, v_r, color=C["record"], lw=0.9); ax0.axvspan(t5, t95, color="#f2efe8", zorder=0)
ax0.set_title(f"(a) FR76 x1.00 record  --  17.3 s, PGV {abs(v_r).max():.2f} m/s, PGA {abs(a_r).max()/9.81:.2f} g; shaded D5-95 = {t95-t5:.1f} s", loc="left", fontweight="bold", fontsize=11*F)
ax1 = fig.add_subplot(gs[1, 0], sharex=ax0); ax1.plot(tt + t_lo, v_s, color=C["synthetic"], lw=0.9)
ax1.set_title(f"(b) spectrum-compatible synthetic  --  {tt[-1]:.1f} s, {N_SINES} sines, PGV {abs(v_s).max():.2f}, PGA {abs(a_s).max()/9.81:.2f} g", loc="left", fontweight="bold", fontsize=11*F)
ax2 = fig.add_subplot(gs[2, 0], sharex=ax0); ax2.plot(np.arange(len(v_2))*DT + t5, v_2, color=C["route2"], lw=1.4)
ax2.set_title(f"(c) Route 2 as run  --  0.35 s, 1 cycle at T1 + 1 at T_eff, PGV {abs(v_2).max():.2f}, PGA {abs(a_2).max()/9.81:.2f} g", loc="left", fontweight="bold", fontsize=11*F)
for ax in (ax0, ax1, ax2):
    ax.set_ylim(-0.3, 0.3); ax.set_ylabel("table velocity (m/s)", fontsize=10*F); ax.grid(alpha=0.25)
ax2.set_xlabel("time (s)", fontsize=10*F); ax2.set_xlim(0, t_r[-1])
axs = fig.add_subplot(gs[:, 1])
axs.plot(Ts, Sd_r, color=C["record"], lw=2.4, label="record")
axs.plot(Ts, Sd_s, color=C["synthetic"], lw=2, label="synthetic (matched at the dots)"); axs.plot(Tn, Sd_tgt, "o", color=C["synthetic"], ms=5)
axs.plot(Ts, Sd_2, color=C["route2"], lw=2, label="Route 2 as run (matched at T1, T_eff)")
for T, lab in ((T1, "T1"), (TEFF, "T_eff")):
    axs.axvline(T, color="#8a8f97", ls=":", lw=1); axs.text(T, 66, lab, ha="center", fontsize=9*F, color="#8a8f97")
axs.set_xscale("log"); axs.set_xlim(0.05, 1.0); axs.set_ylim(0, 70); axs.set_xlabel("period (s)", fontsize=10*F); axs.set_ylabel("5%-damped Sd (mm)", fontsize=10*F)
axs.set_title("(d) what each input demands of a linear SDOF", loc="left", fontweight="bold", fontsize=11*F); axs.legend(fontsize=9*F, frameon=False); axs.grid(alpha=0.25, which="both")
for ax in (ax0, ax1, ax2, axs):
    ax.tick_params(labelsize=9*F)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("fig_synthetic_vs_route2.png", dpi=150); print("-> fig_synthetic_vs_route2.png, synthetic_FR76_x1p00_vel.txt")
