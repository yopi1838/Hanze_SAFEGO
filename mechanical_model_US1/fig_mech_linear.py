import csv, math, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rocking_model import graziotti, bilinear_elastic, doherty_trilinear
t,b,h,zs=210.0,1292.0,2580.0,2060.0; W=13.16; O=(10.35+1.485+1.18)*1000/b; Wl=W*1000/b
K_T1=1200*(2*math.pi/0.091)**2/1e6
fig,axs=plt.subplots(1,2,figsize=(14,5.6),sharey=True)
v5=np.genfromtxt('/home/claude/env/v5_pushover_digitised.csv',delimiter=',',names=True)
neg=[r for r in csv.DictReader(open('/home/claude/v14/pushover_results_v5_ft0_masonrypush/pushover_neg.csv')) if r.get('converged','1')!='0']
dn=np.array([abs(float(r['d_ctrl_mm'])) for r in neg]); Fa=np.array([abs(float(r['a_mps2']))*1.308 for r in neg])
cur=list(csv.DictReader(open('/home/claude/lb/out_v2/loop_backbone_curve.csv')))
for ax,xmax in zip(axs,(110,30)):
    ax.plot(v5['d_mm'],v5['F_kN'],'-',color='#2a78d6',lw=2,label='DEM pushover +, masonry-only')
    ax.plot(dn,Fa,'-',color='#eb6834',lw=2,label='DEM pushover −, applied a·1308 kg')
    for s,ls in (('pos','-'),('neg',':')):
        rr=[r for r in cur if r['sign']==s]; ax.plot([float(r['d_mm']) for r in rr],[float(r['F_p90_kN']) for r in rr],ls,color='#1baf7a',lw=1.6,label='US-1 test envelope, '+('+' if s=='pos' else '−'))
    for e,ls,name in ((0.0,'-','e = 0'),(t/2,'--','e = t/2')):
        r=graziotti(t,h,Wl,O,0.0,1,1,e); F0=r['F0']*b/1000; us=r['uins']*(h-zs)/(h-r['h1'])
        uy,Fy,_=bilinear_elastic(F0,us,K_T1)
        ax.plot([0,uy,us],[0,Fy,0],ls,color='#6250d6',lw=1.8,label='rigid body + elastic K %.1f kN/mm, %s'%(K_T1,name))
        p=doherty_trilinear(F0,us,'new'); ax.plot(*zip(*p),ls,color='#e87ba4',lw=1.8,label='Doherty trilinear (new), %s'%name)
    ax.axhspan(8.3,11.6,color='#999',alpha=0.1,lw=0)
    ax.set_xlim(0,xmax); ax.set_ylim(0,14); ax.grid(alpha=0.25); ax.set_xlabel('|displacement| at 2.06 m, mm')
axs[0].set_ylabel('horizontal force, kN'); axs[0].legend(fontsize=7.5,loc='upper right')
axs[0].set_title('Rigid-body strip with a linear branch added (US-1, f_t = 0)',loc='left',fontsize=10.5)
axs[1].set_title('zoom, 0–30 mm',loc='left',fontsize=10.5)
fig.tight_layout(); fig.savefig('/home/claude/mech/fig_mech_model_linear_branch.png',dpi=160)
