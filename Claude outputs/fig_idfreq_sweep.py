"""Frequency identification of the EUCENTRE SIN model: as built vs kn/ks x 0.31.

Reads the 3DEC ACC_MID exports (2 header lines, whitespace), identifies each
ring-down (zero-crossing period + zero-padded FFT), then fits the
two-compliance model  1/f^2 = C_b + C_j/s  through the two points and solves
it for the two measured frequencies.
"""
import numpy as np, math, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UP = "/root/.claude/uploads/08af91e6-b1c7-5c94-b6cd-980805357292/"
RUNS = [("as built (kn 91 GPa/m)", 1.00, UP + "33c11cf3-ACC_MID.csv", "#b4491a"),
        ("kn, ks x 0.31",          0.31, UP + "35869cbe-ACC_MID.csv", "#2b5f8a")]
MEAS = [("SIN-01-00, 0.1 MPa (target)", 14.27, "#4a7c59"), ("SIN-03-00, 0.3 MPa virgin", 18.75, "#8a8f97")]
PULSE_SEC, BAND = 0.03, (5.0, 35.0)
F = 1.5

def read(p):
    d = np.genfromtxt(p, skip_header=2); m = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1]); return d[m, 0], d[m, 1]
def ring(t, y, t0=PULSE_SEC + 0.02):
    m = t >= t0; return t[m], y[m] - y[m].mean()
def f_zc(t, y):
    s = np.sign(y); i = np.where((s[1:] > 0) & (s[:-1] <= 0))[0]
    tc = t[i] - y[i] * (t[i + 1] - t[i]) / (y[i + 1] - y[i]); per = np.diff(tc)
    per = per[per > 0.25 * np.median(per)]          # drop double crossings from noise at zero
    return 1 / np.mean(per), len(per)
def spec(t, y):
    dt = np.median(np.diff(t)); n = len(y); npad = max(8 * n, 1 << 16)
    f = np.fft.rfftfreq(npad, dt); A = np.abs(np.fft.rfft(y * np.hanning(n), npad)); return f, A / A[(f >= BAND[0]) & (f <= BAND[1])].max()

res = []
for lab, s, p, c in RUNS:
    t, y = read(p); tr, yr = ring(t, y); fz, n = f_zc(tr, yr); f, A = spec(tr, yr)
    m = (f >= BAND[0]) & (f <= BAND[1]); ff = f[m][A[m].argmax()]
    res.append(dict(label=lab, s=s, f_zc=fz, f_fft=ff, n=n, t=tr, y=yr, f=f, A=A, c=c))
    print(f"{lab:28s} s={s:.2f}  zero-cross {fz:.2f} Hz  FFT {ff:.2f} Hz  ({n} cycles)")

# two-compliance fit through the two points: 1/f^2 = Cb + Cj/s  (normalised so Cb+Cj = 1 at s=1)
f0, f1 = res[0]["f_zc"], res[1]["f_zc"]; s1 = res[1]["s"]
r = (f0 / f1) ** 2                       # = Cb + Cj/s1 with Cb + Cj = 1
Cj = (r - 1) / (1 / s1 - 1); Cb = 1 - Cj
fit = lambda s: f0 / np.sqrt(Cb + Cj / s)
solve = lambda ft: Cj / ((f0 / ft) ** 2 - Cb)
print(f"fit: C_blocks={Cb:.3f}  C_joints={Cj:.3f}   exponent f~s^{math.log(f1/f0)/math.log(s1):.2f}")
for lab, fm, _ in MEAS:
    print(f"  {lab:30s} {fm:.2f} Hz -> kn_scale {solve(fm):.3f}")

fig, ax = plt.subplots(1, 3, figsize=(18.0, 6.4), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))
a = ax[0]
for R in res:
    a.plot(R["t"], R["y"] / np.abs(R["y"]).max(), color=R["c"], lw=1.3, label=f'{R["label"]}: {R["f_zc"]:.1f} Hz')
a.set_xlim(0.05, 0.45); a.set_ylim(-1.05, 1.55); a.set_xlabel("time after pulse (s)", fontsize=12 * F); a.set_ylabel("mid-height velocity (normalised)", fontsize=12 * F)
a.set_title("(a) ring-down, 30 ms 0.05 g pulse", fontsize=13 * F, fontweight="bold", loc="left")
a.legend(fontsize=11 * F, loc="upper right", frameon=False)
a = ax[1]
for R in res:
    a.plot(R["f"], R["A"], color=R["c"], lw=1.6)
    a.axvline(R["f_zc"], color=R["c"], ls=":", lw=1)
for lab, fm, c in MEAS:
    a.axvline(fm, color=c, ls="--", lw=1.8); a.text(fm, 1.04, f"{fm:.2f}", color=c, ha="center", fontsize=11 * F)
a.set_xlim(5, 30); a.set_ylim(0, 1.12); a.set_xlabel("frequency (Hz)", fontsize=12 * F); a.set_ylabel("amplitude spectrum", fontsize=12 * F)
a.set_title("(b) spectrum; measured dashed", fontsize=13 * F, fontweight="bold", loc="left")
a = ax[2]
ss = np.linspace(0.15, 1.05, 200); a.plot(ss, fit(ss), color="#18212b", lw=1.6, label=r"fit $1/f^2 = C_b + C_j/s$")
for R in res:
    a.plot(R["s"], R["f_zc"], "o", color=R["c"], ms=9, zorder=5)
for lab, fm, c in MEAS:
    sm = solve(fm); a.axhline(fm, color=c, ls="--", lw=1.4); a.plot(sm, fm, "s", color=c, ms=8, zorder=5)
    a.annotate(f"s = {sm:.2f}", (sm, fm), xytext=(8 * F, -16 * F), textcoords="offset points", color=c, fontsize=11 * F, fontweight="bold")
a.set_xlabel("kn, ks scale s", fontsize=12 * F); a.set_ylabel("first frequency (Hz)", fontsize=12 * F); a.set_xlim(0.15, 1.05); a.set_ylim(10, 22)
a.set_title(f"(c) blocks: {100*Cb:.0f}% of compliance", fontsize=13 * F, fontweight="bold", loc="left")
a.legend(fontsize=11 * F, loc="lower right", frameon=False)
for x in ax:
    x.tick_params(labelsize=10.5 * F); x.grid(alpha=0.25)
    for sp in ("top", "right"): x.spines[sp].set_visible(False)
fig.tight_layout()
out = "/tmp/work/SL/fig_idfreq_kn0p31.png"; fig.savefig(out, dpi=150); print("->", out)
with open("/tmp/work/SL/idfreq_fit.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["item", "value"])
    for R in res: w.writerow([f"f_zc s={R['s']}", f"{R['f_zc']:.3f}"]); w.writerow([f"f_fft s={R['s']}", f"{R['f_fft']:.3f}"])
    w.writerow(["C_blocks", f"{Cb:.4f}"]); w.writerow(["C_joints", f"{Cj:.4f}"])
    for lab, fm, _ in MEAS: w.writerow([f"kn_scale for {fm} Hz", f"{solve(fm):.4f}"])
