#!/usr/bin/env python3
"""NEGATIVE RESULT: a mixture-derived transfer curve pct0(span) scores WORSE than the logistic.

Andy, 2026-10-18: "Having now seen the data, can you come up with a more appropriate and
perhaps elegant fitting function?"

THE IDEA. Section D measured gen1's floodplain PDF as a symmetric core plus an exponential
upper tail of nearly constant e-folding length. If that is the whole structure, then span
and pct0 are not two variables but two readouts of ONE latent quantity -- the tail
fraction f -- and the transfer function is the parametric curve

    d ~ (1-f) Normal(mu, sigma) + f ExpNormal(mu, sigma, L)
    span(f) = Q90 - Q10,   pct0(f) = P(d <= 0)

with three physical parameters (mu, sigma, L) instead of the logistic's four empirical ones.

THE RESULT (block-CV, identical folds and population to floodplain_transfer_function.py):

    window   mixture CV R2   logistic CV R2
    1x1          0.114           0.130
    3x3          0.105           0.147
    5x5          0.063           0.124

WHY IT LOSES, and this is the finding: a FIXED core width cannot reach the data's range.
With sigma ~70 mm the core alone spans 180 mm, but bare cells have spans of 100; and the
mixture saturates at span ~436 while vegetated cells reach 660. The core WIDENS with
vegetation as well -- per cell at 3x3, p50-p10 runs 60 -> 180 mm from the least- to the
most-vegetated quartile (section D). So span is set by TWO latent quantities, core width and
tail fraction, and no one-parameter family through (span, pct0) can follow it. That is also
why span explains only ~14% of per-cell variance in pct0 under any form.

Kept as a record; the logistic stays. Scratch-derived, so it reads the table directly.
"""
import numpy as np, pyarrow.parquet as pq
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares
from scipy.stats import norm, exponnorm
S="data/derived/elbaext_delong_star_nosurf"; LO,HI,DZ=-1500.,1500.,20.
edges=np.arange(LO,HI+DZ,DZ); ctr=0.5*(edges[:-1]+edges[1:]); NB=ctr.size
sl=np.load(f"{S}/slope.npy"); fp=np.load(f"{S}/floodplain_mask.npy").astype(bool); sp=np.load(f"{S}/gen1_spread_mm.npy"); ny,nx=sl.shape
t=pq.read_table(f"{S}/beam_offset_table.parquet",columns=["cell","d_mm_corr","in_grid"]); g=t["in_grid"].to_numpy().astype(bool)
cell=t["cell"].to_numpy()[g]; d=t["d_mm_corr"].to_numpy()[g].astype(float); ok=np.isfinite(d)&(d>LO)&(d<HI); cell,d=cell[ok],d[ok]
H=np.zeros((ny*nx,NB),np.int32); np.add.at(H,(cell,np.clip(((d-LO)/DZ).astype(int),0,NB-1)),1); H=H.reshape(ny,nx,NB)
flat=fp&(sl<=2)&np.isfinite(sp)&(H.sum(2)>=3)
iy,ix=np.nonzero(flat); block=(iy//100)*10000+ix//100; ub,binv=np.unique(block,return_inverse=True)
fold=np.random.default_rng(0).permutation(np.arange(ub.size)%5)[binv]
# fine grid for the model curve (10 mm) so span is resolved better than the data's 20 mm bins
fc=np.arange(-1500,1500,2.0)
def curve(mu,sig,L,fs=np.linspace(0,0.999,400)):
    g=norm.pdf(fc,mu,sig); e=exponnorm.pdf(fc,L/sig,loc=mu,scale=sig)
    out=[]
    for f in fs:
        pdf=(1-f)*g+f*e; P=pdf/pdf.sum(); C=np.cumsum(P)
        out.append((fc[np.searchsorted(C,.9)]-fc[np.searchsorted(C,.1)], 100*P[fc<=0].sum()))
    out=np.array(out); return out[:,0], out[:,1]
def G(p,x):
    s,pc=curve(*p); o=np.argsort(s); s,pc=s[o],pc[o]
    return np.interp(x,s,pc)   # beyond max span the model saturates at its floor
def fitG(x,y):
    r=least_squares(lambda p:G(p,x)-y,x0=[-15,60,180],bounds=([-300,10,20],[300,400,1000]),loss="soft_l1",f_scale=10)
    return r.x
def logistic(p,x): return p[1]+(p[0]-p[1])/(1+np.exp(np.clip((x-p[2])/p[3],-50,50)))
def fitL(x,y):
    return least_squares(lambda p:logistic(p,x)-y,x0=[70,20,200,80],loss="soft_l1",f_scale=10,bounds=([0,0,0,1],[100,100,2000,1000])).x
def cv(x,y,fitf,modf):
    pred=np.empty_like(y)
    for f in range(5):
        tr=fold!=f; pred[~tr]=modf(fitf(x[tr],y[tr]),x[~tr])
    res=y-pred; return 1-res.var()/y.var(), 1.4826*np.median(np.abs(res-np.median(res)))
for k in (1,3,5):
    Hp=H.astype(np.float32) if k==1 else np.rint(uniform_filter(H.astype(np.float32),size=(k,k,1),mode="nearest")*k*k)
    F=Hp.reshape(-1,NB)[flat.ravel()]; tot=F.sum(1); C=np.cumsum(F,1); q=lambda p: ctr[np.argmax(C>=p*tot[:,None],axis=1)]
    x=q(.9)-q(.1); y=100*F[:,ctr<=0].sum(1)/tot
    # subsample for the CV fits (the mixture curve is ~100x costlier than the logistic); score on all held-out cells
    pm=fitG(x,y); pl=fitL(x,y)
    r2g,ng=cv(x,y,fitG,G); r2l,nl=cv(x,y,fitL,logistic)
    print(f"{k}x{k}: MIXTURE mu {pm[0]:.1f} sigma {pm[1]:.1f} L {pm[2]:.1f}  CV R2 {r2g:.3f} NMAD {ng:.1f}   |   LOGISTIC CV R2 {r2l:.3f} NMAD {nl:.1f}", flush=True)
    if k==3:
        bi=np.minimum((np.argsort(np.argsort(x))*20)//x.size,19)
        print(f"   {'obs span':>9}{'obs pct0':>9}{'mixture':>9}{'logistic':>9}")
        for b in range(20):
            xs=np.median(x[bi==b]); print(f"   {xs:>9.0f}{np.median(y[bi==b]):>9.1f}{G(pm,np.array([xs]))[0]:>9.1f}{logistic(pl,xs):>9.1f}")
        s,pc=curve(*pm,fs=np.linspace(0,1,11)); print("   model f->(span,pct0):", [(f"{a:.0f}",f"{b:.0f}") for a,b in zip(s,pc)])
    del Hp,F
