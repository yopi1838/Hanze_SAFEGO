import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0, '/home/claude/duct'); sys.path.insert(0, '/home/claude/mech')
from vw import curve, config
from rocking_model import doherty_trilinear
W, O = 13.16, 10.35 + 1.485 + 1.18
P = '/home/claude/v14/pushover_results_v5_ft0_masonrypush/'
pos = [r for r in csv.DictReader(open(P + 'pushover_pos.csv')) if r['converged'] != '0']
neg = [r for r in csv.DictReader(open(P + 'pushover_neg.csv')) if r['converged'] != '0']
dp = np.array([abs(float(r['d_ctrl_mm'])) for r in pos]); Fp = np.array([abs(float(r['F_tot_kN'])) for r in pos])
ap = np.array([abs(float(r['a_mps2'])) for r in pos]); jN = np.array([float(r['joist_N_kN']) for r in pos])
dn = np.array([abs(float(r['d_ctrl_mm'])) for r in neg]); Fn = np.array([abs(float(r['F_tot_kN'])) for r in neg])
v5 = np.genfromtxt('/home/claude/env/v5_pushover_digitised.csv', delimiter=',', names=True)

# joist rotational stiffness back-calculated from the logged joist support reaction (pos run, hinge assumed at 2.0 m)
H1, L = 2.0, 1.0
th = np.linspace(1e-5, 0.03, 300); dd = []
for a in th:
    lo, up, _ = config(a, H1); dd.append((up(np.array([0.21, 2.06]))[0] - 0.21) * 1e3)
dd = np.array(dd); th_at = np.interp(dp, dd, th)
dN = jN - jN[0]; m = dp > 3
k_r = np.sum(dN[m] * L * 1e3 * th_at[m]) / np.sum(th_at[m] ** 2)        # Nm/rad, LS through the origin
print('k_r back-calculated %.3g Nm/rad' % k_r)

fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw=dict(width_ratios=[2.1, 1]))
a0 = ax[0]
a0.plot(v5['d_mm'], v5['F_kN'], '-', color='#9aa3ad', lw=1.6, label='DEM v5 + (digitised from your figure, to 80 mm)')
a0.plot(dp, Fp, 'o-', color='#2a78d6', ms=3.5, lw=2, label='DEM ft0 +, gravity tilt (F_tot; joist mass loaded)')
a0.plot(dn, Fn, 's-', color='#eb6834', ms=3.5, lw=2, label='DEM ft0 −, masonry-only (F_tot)')
for h1, ls, c in ((1.508, ':', '#6250d6'), (2.0, '--', '#6250d6')):
    d, F = curve(h1, W, O, thmax=0.11)
    a0.plot(d, F, ls, color=c, lw=1.6, label='rigid two-block, hinge at %.2f m: F0 %.1f kN, slope %.3f kN/mm' % (h1, F[0], np.polyfit(d[(d > 40) & (d < 80)], F[(d > 40) & (d < 80)], 1)[0]))
    if h1 == 2.0:
        uz = np.interp(0, -F, d); pts = np.array(doherty_trilinear(F[0], uz, 'new'))
        a0.plot(pts[:, 0], pts[:, 1], '-', color='#1baf7a', lw=1.8, label='Doherty trilinear on the 2.0 m line (new): plateau %.1f kN to %.0f mm' % (pts[1, 1], pts[2, 0]))
a0.set_xlim(0, 110); a0.set_ylim(0, 18); a0.grid(alpha=0.25)
a0.set_xlabel('|displacement| at 2.06 m, mm'); a0.set_ylabel('horizontal force, kN')
a0.set_title('The DEM curve is the rocking line of a hinge at ~2.0 m, rounded as Doherty describes', loc='left', fontsize=10.5)
a0.legend(fontsize=7.6, loc='upper right', framealpha=0.95)
a1 = ax[1]
a1.plot(dp, dN, 'o-', color='#2a78d6', ms=3.5, lw=2, label='joist support reaction − baseline (pos run)')
a1.plot(dp, k_r * th_at / L / 1e3, '--', color='#c2410c', lw=1.5, label='k_r θ / L,  k_r = %.1e Nm/rad, L = 1 m' % k_r)
a1.set_xlabel('|displacement| at 2.06 m, mm'); a1.set_ylabel('ΔN on the joist supports, kN'); a1.grid(alpha=0.25)
a1.set_title('The joists press harder on their supports\nas the wall rotates (pos run)', loc='left', fontsize=10.5)
a1.legend(fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig('/home/claude/mech/fig_pushover_ductility.png', dpi=160)
