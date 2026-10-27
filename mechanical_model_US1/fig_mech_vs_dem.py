import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rocking_model import graziotti, bilinear_elastic
t,b,h,zs=210.0,1292.0,2580.0,2060.0; W=13.16; O=(10.35+1.485+1.18)*1000/b; Wl=W*1000/b
K=4.6   # DEM K_el over 0.1-1.5 mm (both branches)
fig,ax=plt.subplots(figsize=(11,6))
v5=np.genfromtxt('/home/claude/env/v5_pushover_digitised.csv',delimiter=',',names=True)
ax.plot(v5['d_mm'],v5['F_kN'],'-',color='#2a78d6',lw=2.2,label='DEM pushover +, masonry-only (F_tot)')
neg=[r for r in csv.DictReader(open('/home/claude/v14/pushover_results_v5_ft0_masonrypush/pushover_neg.csv')) if r.get('converged','1')!='0']
dn=np.array([abs(float(r['d_ctrl_mm'])) for r in neg]); Fa=np.array([abs(float(r['a_mps2']))*1.308 for r in neg])
ax.plot(dn,Fa,'-',color='#eb6834',lw=2.2,label='DEM pushover −, masonry-only (applied a·1308 kg)')
for e,ls,name in ((0.0,'-','e = 0'),(t/2,'--','e = t/2')):
    r=graziotti(t,h,Wl,O,0.0,1,1,e); F0=r['F0']*b/1000; us=r['uins']*(h-zs)/(h-r['h1'])
    uy,Fy,_=bilinear_elastic(F0,us,K)
    ax.plot([0,uy,us],[0,Fy,0],ls,color='#6250d6',lw=1.8,label='three-hinge strip (Graziotti) + K %.1f, %s: F0 %.1f kN'%(K,name,F0))
# base rocking + top-course sliding under the fixed beam, top friction from the DEM (R_top 3.7 kN)
Wd=1.308*9.81; Od=10.35+0.213*9.81; Rt=3.7
u=np.linspace(0,110,300); F=(Wd*(0.21-2.58*u/2060)+Od*0.21)/2.58+2*Rt
uy=(F[0])/K; Fm=np.minimum(K*u,F)
ax.plot(u,Fm,'-',color='#1baf7a',lw=1.8,label='base rocking + top-course sliding, R_top 3.7 kN: F0 %.1f kN'%F[0])
ax.plot(u,np.minimum(K*u,F-2*Rt),':',color='#1baf7a',lw=1.6,label='   same mechanism without top friction: %.1f kN'%(F[0]-2*Rt))
ax.set_xlim(0,110); ax.set_ylim(0,14); ax.grid(alpha=0.25)
ax.set_xlabel('|displacement| at 2.06 m, mm'); ax.set_ylabel('horizontal force, kN')
ax.set_title('Why the DEM pushover is so ductile: the mechanism it forms, not the load',loc='left',fontsize=11)
ax.legend(fontsize=8,loc='upper right',framealpha=0.95)
fig.tight_layout(); fig.savefig('/home/claude/mech/fig_mech_vs_dem_mechanisms.png',dpi=160)
