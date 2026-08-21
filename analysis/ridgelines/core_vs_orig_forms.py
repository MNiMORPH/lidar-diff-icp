#!/usr/bin/env python3
"""Compare ground-return distribution FORMS on the robust CORE strata vs the original
penetration strata: does removing marginal/contaminated cells sharpen the forms?
CORE from ground_fine_core.npz, ORIGINAL from ground_fine_csf_all.npz (same window/transform).
Reports moments + robust stats + log-derivative tail form; overlays the distributions.

    ./lidar-icp/bin/python analysis/ridgelines/core_vs_orig_forms.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CORE=np.load("data/derived/elba_fulldensity/ground_fine_core.npz")
ORIG=np.load("data/derived/elba_fulldensity/ground_fine_csf_all.npz")
fc=CORE["fc"]; FW=fc[1]-fc[0]
def wq(c,q): cdf=np.cumsum(c)/c.sum(); return np.interp(q,cdf,fc)
def logderiv(c,h=0.03):
    N=c.sum(); f=c/N/FW; ok=c>0; x=fc[ok]; y=np.log(f[ok]); w0=c[ok]; s=np.full(fc.shape,np.nan)
    for i,d0 in enumerate(fc):
        k=np.exp(-0.5*((x-d0)/h)**2)*w0
        if k.sum()<50: continue
        X=x-d0; Sw=k.sum();Sx=(k*X).sum();Sxx=(k*X*X).sum();Sy=(k*y).sum();Sxy=(k*X*y).sum()
        den=Sw*Sxx-Sx*Sx
        if abs(den)>1e-12: s[i]=(Sw*Sxy-Sx*Sy)/den
    return s
def form(c,mode,rng):
    s=logderiv(c); m=(fc>=mode+rng[0])&(fc<=mode+rng[1])&np.isfinite(s)
    if m.sum()<3: return "n/a"
    sm=s[m];dm=fc[m]-mode
    return "EXP" if np.std(sm)/abs(np.mean(sm))<np.std(sm*dm)/abs(np.mean(sm*dm)) else "POW"
def describe(c,lbl):
    N=c.sum(); p=c/N; mean=(p*fc).sum(); std=np.sqrt((p*(fc-mean)**2).sum())
    skew=(p*((fc-mean)/std)**3).sum(); med=wq(c,.5); mode=fc[np.argmax(c)]
    o=np.argsort(np.abs(fc-med)); nmad=1.4826*np.interp(.5,np.cumsum(c[o])/N,np.abs(fc-med)[o])
    print(f"    {lbl:16s} n={N:>9,.0f}  mode {mode:+.3f} med {med:+.3f}  sig {std*1000:3.0f} NMAD {nmad*1000:3.0f}mm "
          f"skew {skew:+.2f}  tails {form(c,mode,(-0.30,-0.12))}/{form(c,mode,(0.12,0.30))}")

print("Ground-return forms: CORE vs ORIGINAL strata (gen1 CSF cloth, gen2 internal)")
for key,lab in [("g1csf","gen1 CSF"),("g2gnd","gen2")]:
    for s,strat in [(1,"FOREST"),(2,"OPEN")]:
        print(f"  {lab} {strat}:")
        describe(ORIG[f"{key}_{s}"].astype(float), "original")
        describe(CORE[f"{key}_{s}"].astype(float), "CORE")

# ---- figure: core (bold) vs original (faint), log density, ground zoom --------------
fig,axes=plt.subplots(2,2,figsize=(13,11),sharey='row')
for jr,(key,lab) in enumerate([("g1csf","gen1 CSF cloth ground"),("g2gnd","gen2 internal ground")]):
    for jc,(s,strat) in enumerate([(1,"FOREST"),(2,"OPEN")]):
        ax=axes[jr,jc]
        for src,style,al,tag in [(ORIG,":",0.5,"original"),(CORE,"-",1.0,"CORE")]:
            c=src[f"{key}_{s}"].astype(float); d=c/c.sum()/FW
            col="C0" if s==1 else "C1"
            ax.semilogx(np.where(d>0,d,np.nan),fc,col,ls=style,lw=1.6,alpha=al,label=tag)
        ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.5,0.5); ax.set_xlim(1e-2,3e1)
        ax.set_xlabel("density (1/m, log)"); ax.set_title(f"{lab} — {strat}"); ax.legend(fontsize=8); ax.grid(alpha=.3,which="both")
axes[0,0].set_ylabel("slope-normal d (m)"); axes[1,0].set_ylabel("slope-normal d (m)")
fig.suptitle("Ground-return distributions: CORE (solid) vs ORIGINAL (dotted) strata",y=0.995)
fig.savefig("figures/refdatum/core_vs_orig_forms.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/core_vs_orig_forms.png")
