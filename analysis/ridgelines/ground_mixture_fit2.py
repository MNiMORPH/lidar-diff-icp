#!/usr/bin/env python3
"""Fit gen2 (and gen1) ground-class fine (1 cm) distributions; test whether a second
'plant' component is resolvable, and extract the ground level.  Reuses the pooled 1 cm
histograms from ground_mixture_fit.py (no re-streaming).

Frame-INVARIANT diagnostics (independent of the reference plane, which is gen2 bare earth):
  - Bowley skew (p90+p10-2p50)/(p90-p10): a one-sided plant shoulder -> positive skew.
  - per-percentile rise gen2-gen1: an added UPPER shoulder makes the rise grow toward high
    percentiles; deeper gen1 PENETRATION makes it grow toward low percentiles.
  - 1-Gaussian vs 2-Gaussian AIC: is a plant component even identifiable?

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/ground_mixture_fit2.py
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]
def g(gen,s): return Z[f"{gen}_{s}"]

def q(c,p):
    cdf=np.cumsum(c)/c.sum(); return np.interp(p,cdf,fc)
def bowley(c):
    p10,p50,p90=q(c,.1),q(c,.5),q(c,.9); return (p90+p10-2*p50)/(p90-p10)
def gauss(d,mu,sig): return np.exp(-0.5*((d-mu)/sig)**2)/(sig*np.sqrt(2*np.pi))

def nll1(p,c):
    mu,ls=p; f=np.clip(gauss(fc,mu,np.exp(ls)),1e-12,None)*FW; return -(c*np.log(f)).sum()
def nll2(p,c):
    mu,ls,off,lsp,lw=p; sg=np.exp(ls); sp=np.exp(lsp); w=1/(1+np.exp(-lw))
    f=w*gauss(fc,mu,sg)+(1-w)*gauss(fc,mu+np.exp(off),sp)
    f=np.clip(f,1e-12,None)*FW; return -(c*np.log(f)).sum()

def fit(gen,s,lbl):
    c=g(gen,s); mode=fc[np.argmax(c)]; med=q(c,.5)
    r1=minimize(nll1,[mode,np.log(0.06)],args=(c,),method="Nelder-Mead",
                options=dict(maxiter=8000,xatol=1e-6,fatol=1e-3))
    mu1,ls1=r1.x; aic1=2*2+2*r1.fun
    r2=minimize(nll2,[mode,np.log(0.05),np.log(0.15),np.log(0.15),0.6],args=(c,),
                method="Nelder-Mead",options=dict(maxiter=20000,xatol=1e-6,fatol=1e-3))
    mu,ls,off,lsp,lw=r2.x; sg=np.exp(ls); sp=np.exp(lsp); w=1/(1+np.exp(-lw)); mup=mu+np.exp(off)
    aic2=2*5+2*r2.fun
    print(f"  {lbl:12s}: mode {mode:+.3f} med {med:+.3f}  skew {bowley(c):+.3f}  "
          f"| 1G mu {mu1:+.4f} sig {np.exp(ls1):.3f} AIC {aic1:.0f}  "
          f"| 2G ground mu {mu:+.4f} sig {sg:.3f} w {w:.2f}; plant +{np.exp(off):.3f}m sig {sp:.3f}  "
          f"dAIC {aic2-aic1:+.0f}")
    return dict(mu1=mu1,sig1=np.exp(ls1),mode=mode,med=med,mu=mu,sg=sg,w=w,mup=mup,sp=sp)

print("Fine (1 cm) ground-class fits (dAIC<0 => 2-Gaussian justified):")
R={}
for gen in ["gen2","gen1"]:
    for s,lbl in [(1,f"{gen} forest"),(2,f"{gen} open")]:
        R[(gen,s)]=fit(gen,s,lbl)

print("\nPer-percentile rise gen2-gen1 (mm) — trend reveals mechanism:")
for s,lbl in [(1,"forest"),(2,"open")]:
    a=g("gen1",s); b=g("gen2",s)
    row=" ".join(f"p{int(p*100):02d} {(q(b,p)-q(a,p))*1000:+4.0f}" for p in [.1,.25,.5,.75,.9])
    print(f"  {lbl:6s}: {row}   (grows toward LOW p => gen1 penetration; toward HIGH p => gen2 shoulder)")

# ---- clean overlay figure, non-circular presentation --------------------------------
fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,s,lbl in [(axes[0],1,"FOREST"),(axes[1],2,"OPEN")]:
    for gen,col in [("gen1","C0"),("gen2","C3")]:
        c=g(gen,s); ax.plot(c/c.sum()/FW,fc,col,lw=1.7,label=f"{gen} ground-class")
        ax.axhline(R[(gen,s)]["med"],color=col,ls=":",lw=1)
    ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.5,0.6)
    ax.set_xlabel("density (1/m)"); ax.set_title(f"{lbl}: fine 1 cm ground-class")
    ax.legend(fontsize=9)
axes[0].set_ylabel("slope-normal d (m)  [plane = gen2 bare earth]")
fig.suptitle("Fine ground-class: gen1 deeper low tail (penetration); upper tails ~coincide")
fig.savefig("figures/refdatum/ground_fine_overlay.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/ground_fine_overlay.png")
