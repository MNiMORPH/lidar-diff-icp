#!/usr/bin/env python3
"""Plot & characterize the CSF gen1 ground return distribution and the ALL-returns
distributions (both epochs), from ground_fine_csf_all.npz.

Fig 1: gen1 GROUND by classifier — CSF cloth vs 2008 vendor vs gen2 internal — forest+open, lin+log.
Fig 2: ALL returns gen1 vs gen2 — forest+open, log density (ground peak + low veg + understory).
Also prints moments + log-derivative tail form for the CSF ground.

    ./lidar-icp/bin/python analysis/ridgelines/ground_csf_all_analyze.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_csf_all.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]
def g(n,s): return Z[f"{n}_{s}"].astype(float)
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

def stats(c,lbl):
    N=c.sum(); p=c/N; mean=(p*fc).sum(); std=np.sqrt((p*(fc-mean)**2).sum())
    skew=(p*((fc-mean)/std)**3).sum(); kurt=(p*((fc-mean)/std)**4).sum()-3
    med=wq(c,.5); mode=fc[np.argmax(c)]
    o=np.argsort(np.abs(fc-med)); nmad=1.4826*np.interp(.5,np.cumsum(c[o])/N,np.abs(fc-med)[o])
    s=logderiv(c)
    def tf(rng):
        m=(fc>=mode+rng[0])&(fc<=mode+rng[1])&np.isfinite(s)
        if m.sum()<3: return "n/a"
        sm=s[m];dm=fc[m]-mode
        return "EXP" if np.std(sm)/abs(np.mean(sm))<np.std(sm*dm)/abs(np.mean(sm*dm)) else "POW"
    print(f"  {lbl:16s} n={N:>10,.0f}  mode {mode:+.3f} med {med:+.3f}  sig {std*1000:.0f} NMAD {nmad*1000:.0f}mm "
          f"skew {skew:+.2f} exk {kurt:+.1f}  tails lo/hi {tf((-0.30,-0.12))}/{tf((0.12,0.30))}")

print("gen1 GROUND-return distributions by classifier (1 cm, slope-normal):")
for s,lbl in [(1,"FOREST"),(2,"OPEN")]:
    print(f" {lbl}:")
    stats(g("g1csf",s),  "gen1 CSF cloth")
    stats(g("g1vend",s), "gen1 vendor cl2")
    stats(g("g2gnd",s),  "gen2 internal")

# ---- Fig 1: gen1 ground by classifier + gen2, forest/open, lin+log -------------------
fig,axes=plt.subplots(2,2,figsize=(13,11),sharey='row')
series=[("g1csf","C0","gen1 CSF cloth"),("g1vend","C4","gen1 vendor cl2"),("g2gnd","C3","gen2 internal")]
for jcol,(s,strat) in enumerate([(1,"FOREST"),(2,"OPEN")]):
    for irow,scale in enumerate(["linear","log"]):
        ax=axes[irow,jcol]; plot=ax.semilogx if scale=="log" else ax.plot
        for n,col,lab in series:
            c=g(n,s); d=c/c.sum()/FW; yv=(np.where(d>0,d,np.nan) if scale=="log" else d)
            plot(yv,fc,col,lw=1.5,label=lab); ax.axhline(wq(c,.5),color=col,ls=":",lw=.8)
        ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.5,0.5)
        if scale=="log": ax.set_xlim(1e-2,3e1); ax.grid(alpha=.3,which="both")
        else: ax.grid(alpha=.3)
        ax.set_xlabel(f"density (1/m{', log' if scale=='log' else ''})")
        ax.set_title(f"{strat} ground return — {scale}"); ax.legend(fontsize=8)
axes[0,0].set_ylabel("slope-normal d (m)"); axes[1,0].set_ylabel("slope-normal d (m)")
fig.suptitle("gen1 GROUND: CSF cloth vs 2008 vendor vs gen2 internal (plane = gen2 bare earth)",y=0.995)
fig.savefig("figures/refdatum/ground_csf_vs_vendor.png",dpi=130,bbox_inches="tight"); plt.close(fig)

# ---- Fig 2: ALL returns gen1 vs gen2, forest/open, log (ground + veg) ----------------
fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,(s,strat) in zip(axes,[(1,"FOREST"),(2,"OPEN")]):
    for n,col,lab in [("g1all","C0","gen1 all returns"),("g2all","C3","gen2 all returns")]:
        c=g(n,s); d=c/c.sum()/FW; ax.semilogx(np.where(d>0,d,np.nan),fc,col,lw=1.5,label=lab)
    for n,col in [("g1csf","C0"),("g2gnd","C3")]:
        ax.axhline(wq(g(n,s),.5),color=col,ls=":",lw=.9)
    ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.6,3.5); ax.set_xlim(1e-3,3e1)
    ax.set_xlabel("density (1/m, log)"); ax.set_title(f"{strat}: ALL returns (dotted=ground median)")
    ax.grid(alpha=.3,which="both"); ax.legend(fontsize=9)
axes[0].set_ylabel("slope-normal d (m)")
fig.suptitle("ALL returns near ground (ground spike + low veg + understory), gen1 vs gen2",y=1.0)
fig.savefig("figures/refdatum/all_returns_fine.png",dpi=130,bbox_inches="tight"); plt.close(fig)

print("\nALL-returns near-ground fraction (mass in the fine window that is within |d|<0.15 m):")
for s,lbl in [(1,"FOREST"),(2,"OPEN")]:
    for n,e in [("g1all","gen1"),("g2all","gen2")]:
        c=g(n,s); near=c[np.abs(fc)<0.15].sum()/c.sum()
        print(f"  {lbl} {e} all: {near*100:.1f}% of windowed returns within 0.15 m of ground plane")
print("\nwrote figures/refdatum/ground_csf_vs_vendor.png, all_returns_fine.png")
