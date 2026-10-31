import csv, io, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M=1.3083
def rd(p):
    L=open(p).read().splitlines()
    rows=list(csv.DictReader(io.StringIO("\n".join(l for l in L if l.strip() and not l.lstrip().lstrip('"').startswith("#")))))
    d=np.array([abs(float(r["d_ctrl_mm"])) for r in rows]); a=np.array([abs(float(r["a_mps2"])) for r in rows])
    ok=np.array([int(float(r["converged"]))==1 for r in rows])
    o=np.argsort(d[ok]); return d[ok][o], (a[ok]*M)[o]
def peaks(p):
    pos,neg=[],[]
    for r in csv.DictReader(open(p)):
        d0=float(r["d_start_mm"])
        pos.append((abs(float(r["d_at_Fmax_pos_mm"])-d0),abs(float(r["F_max_pos_kN"]))))
        neg.append((abs(float(r["d_at_Fmax_neg_mm"])-d0),abs(float(r["F_max_neg_kN"]))))
    return np.array(pos),np.array(neg)
INK,MUT,GRID,SURF="#0b0b0b","#52514e","#e7e6e2","#fcfcfb"
C1,C2,C3,TST="#2a78d6","#eb6834","#1baf7a","#8a8984"
cur={"pos":[("old frame, old pockets","pushover_results_v5_ft0_masonrypush_TopMason/pushover_pos.csv",C1,"-"),
            (None,"pushover_results_v5_ft0_masonrypush/pushover_pos.csv",C1,(0,(4,3))),
            ("z-up, old pockets","pushover_results_v5zup_ft0_TopMason/pushover_pos.csv",C2,"-"),
            ("z-up, pockets at test heights (pk2)","pushover_results_v5zup_pk2_ft0_TopMason/pushover_pos.csv",C3,"-")],
     "neg":[("old frame, old pockets","pushover_results_v5_ft0_masonrypush_TopMason/pushover_neg.csv",C1,"-"),
            ("z-up, pockets at test heights (pk2)","pushover_results_v5zup_pk2_ft0_TopMason/pushover_neg.csv",C3,"-")]}
pp,pn=peaks("exp_envelope_peaks_v3/envelope_peaks_US1.csv")
fig,axs=plt.subplots(1,2,figsize=(13,5.2),sharey=True,facecolor=SURF)
for ax,br,title,pk in ((axs[0],"pos","Positive push (towards the joists)",pp),(axs[1],"neg","Negative push",pn)):
    ax.set_facecolor(SURF)
    ax.scatter(pk[:,0],pk[:,1],s=26,facecolor="none",edgecolor=TST,lw=1.2,zorder=2,label="US-1 test: peak force per run")
    for lab,p,c,ls in cur[br]:
        d,F=rd(p); ax.plot(d,F,color=c,lw=2,ls=ls,zorder=3,label=lab)
    ax.set_title(title,loc="left",color=INK,fontsize=12,fontweight="bold")
    ax.set_xlim(0,82); ax.set_ylim(0,13); ax.set_xlabel("displacement at the 2.06 m sensor (mm)",color=MUT)
    ax.grid(color=GRID,lw=0.8); ax.set_axisbelow(True)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for s in ("left","bottom"): ax.spines[s].set_color("#b5b4ae")
    ax.tick_params(colors=MUT)
axs[0].set_ylabel("applied push |a| x m_masonry (kN)",color=MUT)
axs[0].text(44,6.5,"old frame, elastic-top run\n(same as mason top to 27 mm)",fontsize=8,color=MUT)
h,l=axs[0].get_legend_handles_labels(); axs[0].legend(h,l,fontsize=8.5,frameon=False,loc="lower left")
h,l=axs[1].get_legend_handles_labels(); axs[1].legend(h,l,fontsize=8.5,frameon=False,loc="lower left")
fig.suptitle("f_t = 0 pushovers: frame change (old -> z-up) and pocket heights (old -> test, pk2)",x=0.01,ha="left",color=INK,fontsize=13)
fig.tight_layout(); fig.savefig("fig_pushover_zup_pk2.png",dpi=150,facecolor=SURF)
