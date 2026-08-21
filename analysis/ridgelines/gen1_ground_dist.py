#!/usr/bin/env python3
"""Statistical characterization of gen1's (2008 leaf-off) ground-class return distribution,
FOREST vs OPEN, and their comparison (slope-normal d, 1 cm, ground_fine_pooled.npz).

Physical model: a ground CORE Gaussian(mu,sig) smeared DOWN by a penetration exponential
(left-EMG = Gaussian MINUS Exp(tau)).  mu = ground-surface level, sig = measurement
precision, tau = penetration depth scale.  Forest should carry a real tau (leaf-off pulses
reaching soil through bare canopy); open should be ~symmetric measurement noise (tau small).

    ./lidar-icp/bin/python analysis/ridgelines/gen1_ground_dist.py
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import erfc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]

def wq(c, q):                                   # weighted quantile from a histogram
    cdf = np.cumsum(c)/c.sum(); return np.interp(q, cdf, fc)
def gauss(d,mu,sig): return np.exp(-0.5*((d-mu)/sig)**2)/(sig*np.sqrt(2*np.pi))
def remg(d,m,sig,tau):
    z=(d-m)/sig - sig/tau
    with np.errstate(over="ignore", invalid="ignore"):
        return (0.5/tau)*np.exp(np.clip(0.5*(sig/tau)**2-(d-m)/tau, -700, 700))*erfc(-z/np.sqrt(2))
def lemg(d,mu,sig,tau): return remg(-d,-mu,sig,tau)   # Gaussian - Exp: downward penetration tail

def analyze(key, lbl):
    c = Z[key].astype(float); N = c.sum(); p = c/N
    mean=(p*fc).sum(); std=np.sqrt((p*(fc-mean)**2).sum())
    skew=(p*((fc-mean)/std)**3).sum(); kurt=(p*((fc-mean)/std)**4).sum()-3
    med=wq(c,.5); mode=fc[np.argmax(c)]
    # NMAD: weighted median of |d-med|, x1.4826
    order=np.argsort(np.abs(fc-med)); ad=np.abs(fc-med)[order]; wcdf=np.cumsum(c[order])/N
    nmad=1.4826*np.interp(.5,wcdf,ad)
    bowley=(wq(c,.9)+wq(c,.1)-2*med)/(wq(c,.9)-wq(c,.1))
    def nll(f): ff=np.clip(f,1e-12,None)*FW; return -(c*np.log(ff)).sum()
    opt=dict(method="Nelder-Mead",options=dict(maxiter=30000,xatol=1e-7,fatol=1e-3))
    r1=minimize(lambda x:nll(gauss(fc,x[0],np.exp(x[1]))),[mode,np.log(0.1)],**opt); a1=4+2*r1.fun
    # left-EMG (ground core + DOWNWARD penetration tail) — tests the one-sided hypothesis
    rE=minimize(lambda x:nll(lemg(fc,x[0],np.exp(x[1]),np.exp(x[2]))),[mode,np.log(0.05),np.log(0.08)],**opt); aE=6+2*rE.fun
    tauE=np.exp(rE.x[2])
    # core + broad halo (both means free): the model that fits a symmetric heavy tail
    def m2(x):
        mu,ls,mu2,lsh,lw=x; sg=np.exp(ls); sh=np.exp(lsh); w=1/(1+np.exp(-lw))
        return nll(w*gauss(fc,mu,sg)+(1-w)*gauss(fc,mu2,sh))
    r2=minimize(m2,[mode,np.log(0.06),mode,np.log(0.25),1.0],**opt); a2=10+2*r2.fun
    mu,sg,mu2,sh,w=r2.x[0],np.exp(r2.x[1]),r2.x[2],np.exp(r2.x[3]),1/(1+np.exp(-r2.x[4]))
    if w<0.5: mu,sg,mu2,sh,w=mu2,sh,mu,sg,1-w      # comp1 = narrow core
    print(f"\n=== gen1 {lbl} (n={N:,.0f}) ===")
    print(f"  mode {mode:+.3f}  median {med:+.3f}  mean {mean:+.3f} m")
    print(f"  std {std*1000:.0f}  NMAD {nmad*1000:.0f}  IQR {(wq(c,.75)-wq(c,.25))*1000:.0f} mm   "
          f"skew {skew:+.2f}  Bowley {bowley:+.3f}  exkurt {kurt:+.1f}")
    print(f"  pctl: p05 {wq(c,.05):+.3f} p10 {wq(c,.10):+.3f} p25 {wq(c,.25):+.3f} p50 {med:+.3f} "
          f"p75 {wq(c,.75):+.3f} p90 {wq(c,.90):+.3f} p95 {wq(c,.95):+.3f}")
    print(f"  tail balance: low (<mode-0.2m) {p[fc<mode-0.2].sum()*100:.1f}%  vs high (>mode+0.2m) {p[fc>mode+0.2].sum()*100:.1f}%")
    print(f"  1-Gaussian sig {np.exp(r1.x[1])*1000:.0f} mm (AIC {a1:,.0f})")
    print(f"  penetration-tail (left-EMG) tau {tauE*1000:.0f} mm  dAIC {aE-a1:+,.0f} "
          f"{'REJECTED (not one-sided)' if aE>=a1 else 'supported'}")
    print(f"  core+halo: CORE mu {mu:+.4f} sig {sg*1000:.0f} mm w {w:.2f} | HALO mu {mu2:+.4f} sig {sh*1000:.0f} mm w {1-w:.2f}"
          f"  dAIC {a2-a1:+,.0f}")
    return dict(c=c,N=N,mode=mode,med=med,mean=mean,std=std,nmad=nmad,skew=skew,bowley=bowley,
                kurt=kurt,mu=mu,sig=sg,w=w,mu2=mu2,sh=sh,tauE=tauE,sig1=np.exp(r1.x[1]))

F=analyze("gen1_1","FOREST"); O=analyze("gen1_2","OPEN")

print("\n=== COMPARISON (forest vs open) ===")
print(f"  ground mode:       forest {F['mode']:+.3f}  open {O['mode']:+.3f} m  -> {'SAME' if abs(F['mode']-O['mode'])<0.011 else 'differ'}")
print(f"  core mu (level):   forest {F['mu']:+.3f}  open {O['mu']:+.3f} m  -> forest {(F['mu']-O['mu'])*1000:+.0f} mm")
print(f"  core sig:          forest {F['sig']*1000:.0f}  open {O['sig']*1000:.0f} mm  -> forest {F['sig']/O['sig']:.1f}x wider")
print(f"  halo sig / weight: forest {F['sh']*1000:.0f} mm @ {(1-F['w'])*100:.0f}%  open {O['sh']*1000:.0f} mm @ {(1-O['w'])*100:.0f}%")
print(f"  total std:         forest {F['std']*1000:.0f}  open {O['std']*1000:.0f} mm  ({F['std']/O['std']:.1f}x)")
print(f"  NMAD (robust):     forest {F['nmad']*1000:.0f}  open {O['nmad']*1000:.0f} mm  ({F['nmad']/O['nmad']:.1f}x)")
print(f"  skew:              forest {F['skew']:+.2f} (heavy tails BOTH sides)  open {O['skew']:+.2f} (upper-tail skew)")
print(f"  penetration-tail:  REJECTED both -> forest breadth is broad ~symmetric scatter, not a one-sided tail")

# ---- figure: forest vs open, linear + log, with core+halo fits ----------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 11), sharey='row')
for col,(R,lbl) in enumerate([(F,"FOREST"),(O,"OPEN")]):
    d=R["c"]/R["N"]/FW
    core=R["w"]*gauss(fc,R["mu"],R["sig"]); halo=(1-R["w"])*gauss(fc,R["mu2"],R["sh"]); mix=core+halo
    for row,scale in enumerate(["linear","log"]):
        a=axes[row,col]; plot=a.semilogx if scale=="log" else a.plot
        yv=lambda z:(np.where(z>0,z,np.nan) if scale=="log" else z)
        plot(yv(d),fc,"C0",lw=1.6,label=f"gen1 {lbl.lower()} data")
        plot(yv(mix),fc,"k--",lw=1.3,label="core+halo fit")
        plot(yv(core),fc,"C2",lw=1.0,alpha=.8,label=f"core sig {R['sig']*1000:.0f} mm ({R['w']*100:.0f}%)")
        plot(yv(halo),fc,"C1",lw=1.0,alpha=.8,label=f"halo sig {R['sh']*1000:.0f} mm ({(1-R['w'])*100:.0f}%)")
        a.axhline(R["mode"],color="C0",ls=":",lw=1); a.axhline(0,color="k",lw=.5); a.set_ylim(-0.6,0.5)
        if scale=="log": a.set_xlim(1e-2,2e1); a.grid(alpha=.3,which="both")
        else: a.grid(alpha=.3)
        a.set_xlabel(f"density (1/m{', log' if scale=='log' else ''})")
        a.set_title(f"gen1 {lbl} — {scale}  (NMAD {R['nmad']*1000:.0f} mm, skew {R['skew']:+.1f})")
        a.legend(fontsize=8)
axes[0,0].set_ylabel("slope-normal d (m)"); axes[1,0].set_ylabel("slope-normal d (m)")
fig.suptitle("gen1 (2008 leaf-off) ground return: FOREST broad+heavy-tailed (canopy scatter); "
             "OPEN tighter, upper-tail skew — SAME mode", y=0.995)
fig.savefig("figures/refdatum/gen1_forest_vs_open.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_forest_vs_open.png")
