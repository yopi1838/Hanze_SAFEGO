"""Two-wave Route 2 variants: same Sd targets (record Sd at T1 and at T_eff), different wave periods.
   S  : structure periods (current)          A@T1,   B@T_eff
   SJ : structure periods, jointly tuned     (2x2 solve on the combined signal)
   P  : both waves at the record's T_p       A@T_p,  B@T_p
   PM : record short + long periods          A@T_p,  B@T_m (Rathje)
"""
import numpy as np, csv
DT=0.005; XI=0.05; T1,TEFF=0.112,0.241
def sd_spectrum(ag, dt, Ts, xi=XI):
    ag=np.asarray(ag,float); Ts=np.atleast_1d(np.asarray(Ts,float)); w=2*np.pi/Ts; c=2*xi*w; w2=w*w
    u=np.zeros_like(w); v=np.zeros_like(w); a=-ag[0]-c*v-w2*u; umax=np.abs(u)
    for i in range(1,len(ag)):
        un=u+dt*v+0.25*dt*dt*a; vn=v+0.5*dt*a; an=(-ag[i]-c*vn-w2*un)/(1+0.5*dt*c+0.25*dt*dt*w2)
        u=un+0.25*dt*dt*an; v=vn+0.5*dt*an; a=an; umax=np.maximum(umax,np.abs(u))
    return umax
def cyc(V,T):
    spc=int(round(T/DT)); tt=np.arange(spc+1)*DT; return V*np.sin(2*np.pi*tt/(spc*DT))
def join(vA,vB): return np.r_[vA,vB[1:]]
def sd(v,Ts): return sd_spectrum(np.gradient(v,DT),DT,Ts)*1e3
def tune_sep(TA,TB,SdA,SdB):
    a=SdA/sd(cyc(1,TA),[T1])[0]; b=SdB/sd(cyc(1,TB),[TEFF])[0]; return join(cyc(a,TA),cyc(b,TB))
def tune_joint(TA,TB,SdA,SdB):
    # peak-response is not linear in (a,b) because peaks may occur in different cycles; iterate a 2x2 secant
    a,b=SdA/sd(cyc(1,TA),[T1])[0], SdB/sd(cyc(1,TB),[TEFF])[0]
    for _ in range(40):
        v=join(cyc(a,TA),cyc(b,TB)); s=sd(v,[T1,TEFF]); a*=SdA/s[0]; b*=SdB/s[1]
        if abs(s[0]-SdA)<1e-3 and abs(s[1]-SdB)<1e-3: break
    return join(cyc(a,TA),cyc(b,TB))
chk=[T1,TEFF,0.30,0.40,0.50,0.60]
def report(name,v,rec):
    S=sd(v,chk); a=np.gradient(v,DT)
    print(f"  {name:38s}"+"".join(f"{x:8.2f}" for x in S)+f"  PGV {abs(v).max():.3f}  PGA {abs(a).max()/9.81:.2f} g  dur {len(v)*DT:.2f} s")
    return [name]+[f"{x:.3f}" for x in S]+[f"{abs(v).max():.4f}",f"{abs(a).max()/9.81:.3f}"]
rows=[]
for rec,Tp,Tm,SdT1,SdTe,Sdrec in [("FR76 x1.00",0.26,0.374,1.73,13.39,{T1:1.73,TEFF:13.39,0.30:18.0,0.40:33.5,0.50:44.4,0.60:38.5}),
                                   ("HU12 x1.00",0.12,None,0.97,1.85,None)]:
    print(f"\n{rec}: targets Sd(T1)={SdT1} Sd(T_eff)={SdTe}   T_p={Tp}  T_m={Tm}")
    print("  "+" "*38+"".join(f"{'Sd@%.2f'%T:>8s}" for T in chk))
    if Sdrec: print("  "+f"{'record':38s}"+"".join(f"{Sdrec[T]:8.2f}" for T in chk))
    rows.append(report("S  : A@T1  B@T_eff (current)",tune_sep(T1,TEFF,SdT1,SdTe),rec))
    rows.append(report("SJ : A@T1  B@T_eff, joint tune",tune_joint(T1,TEFF,SdT1,SdTe),rec))
    rows.append(report(f"P  : A@T_p B@T_p  (T_p {Tp})",tune_sep(Tp,Tp,SdT1,SdTe),rec))
    rows.append(report(f"PJ : A@T_p B@T_p, joint tune",tune_joint(Tp,Tp,SdT1,SdTe),rec))
    if Tm: rows.append(report(f"PM : A@T_p B@T_m ({Tm})",tune_sep(Tp,Tm,SdT1,SdTe),rec))
    if Tm: rows.append(report(f"PMJ: A@T_p B@T_m, joint tune",tune_joint(Tp,Tm,SdT1,SdTe),rec))
with open("route2_twowave_variants.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["variant"]+[f"Sd_{T:.3f}" for T in chk]+["PGV","PGA_g"]); w.writerows(rows)
