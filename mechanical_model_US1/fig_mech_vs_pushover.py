import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rocking_model import graziotti
t,b,h,zs=210.0,1292.0,2580.0,2060.0; W=13.16; O_top=10.35+1.485; Fj=1.18
Wl=W*1000/b; O=(O_top+Fj)*1000/b
C=dict(dem="#2a78d6",neg="#eb6834",test="#1baf7a",gz="#6250d6",a2="#eda100",onur="#e34948")
fig,ax=plt.subplots(figsize=(11,6.2))
# DEM pushovers (masonry-only push, d at 2.06 m)
v5=np.genfromtxt('/home/claude/env/v5_pushover_digitised.csv',delimiter=',',names=True)
ax.plot(v5['d_mm'],v5['F_kN'],'-',color=C['dem'],lw=2,label='DEM pushover +, masonry-only (v5, digitised)')
neg=list(csv.DictReader(open('/home/claude/v14/pushover_results_v5_ft0_masonrypush/pushover_neg.csv')))
neg=[r for r in neg if r.get('converged','1')!='0']
dn=np.array([abs(float(r['d_ctrl_mm'])) for r in neg]); Fa=np.array([abs(float(r['a_mps2']))*1.308 for r in neg]); Ft=np.array([abs(float(r['F_tot_kN'])) for r in neg])
ax.plot(dn,Fa,'-',color=C['neg'],lw=2,label='DEM pushover −, masonry-only: applied a·1308 kg')
ax.plot(dn,Ft,'--',color=C['neg'],lw=1.4,label='DEM pushover −: F_tot (reactions)')
# test envelope (v2 force, simultaneous branches, p90), both signs folded
cur=list(csv.DictReader(open('/home/claude/lb/out_v2/loop_backbone_curve.csv')))
for s,ls,lab in (('pos','-','US-1 test envelope, + (v2 force, loop branches p90)'),('neg',':','US-1 test envelope, −')):
    rr=[r for r in cur if r['sign']==s]; ax.plot([float(r['d_mm']) for r in rr],[float(r['F_p90_kN']) for r in rr],ls,color=C['test'],lw=1.8,label=lab)
# Graziotti single strip (A1), rigid-body lines F0 (1 - u/u_ins), u read at 2.06 m
for e,ls,lab in ((0.0,'-','Graziotti strip, e = 0 (RB0)'),(t/2,'--','Graziotti strip, e = t/2 (RBΔ)')):
    r=graziotti(t,h,Wl,O,0.0,1,1,e); F0=r['F0']*b/1000; us=r['uins']*(h-zs)/(h-r['h1'])
    ax.plot([0,us],[F0,0],ls,color=C['gz'],lw=2,label='%s: %.1f kN, h1 %.2f m'%(lab,F0,r['h1']/1000))
# A2 hinge at the top joist
for e,ls in ((0.0,'-'),(t/2,'--')):
    h1=2140.0; F0=(2/h1*(Wl+O)*t+O/(h-h1)*(t+2*e))*b/1000
    ui=(2/h1*(Wl+O)*t+O*(t+2*e)/(h-h1))/(2/h1*(Wl+O)+2*O/(h-h1)); us=ui*zs/h1
    ax.plot([0,us],[F0,0],ls,color=C['a2'],lw=1.6,label='strip, hinge at top joist (2.14 m), e = %s: %.1f kN'%('0' if e==0 else 't/2',F0))
# Onur split (B): capacity only
for F,lab in ((16.76,'e = 0'),(21.88,'e = t/2')):
    ax.plot([0],[F],'s',color=C['onur'],ms=8,mec='white',mew=1.5,clip_on=False)
    ax.annotate("Onur's cavity split applied as is, %s: %.1f kN"%(lab,F),(0,F),xytext=(8,0),textcoords='offset points',va='center',fontsize=8.5,color='#333')
ax.axhspan(8.3,11.6,color='#999999',alpha=0.10,lw=0); ax.text(42,11.8,'test plateau 8.3–11.6 kN (runs 16–24)',ha='left',fontsize=8,color='#555')
ax.set_xlim(0,110); ax.set_ylim(0,23); ax.grid(alpha=0.25)
ax.set_xlabel('|displacement| at 2.06 m (top-quarter sensor), mm'); ax.set_ylabel('horizontal force, kN')
ax.set_title('US-1 solid wall: rigid-body mechanical model (Graziotti et al. 2016) vs DEM pushover and test',loc='left',fontsize=11)
ax.legend(fontsize=7.6,loc='upper right',ncol=1,framealpha=0.95)
fig.tight_layout(); fig.savefig('/home/claude/mech/fig_mech_model_vs_pushover.png',dpi=160)
