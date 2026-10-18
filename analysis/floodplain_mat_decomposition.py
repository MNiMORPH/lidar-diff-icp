#!/usr/bin/env python3
"""Is the plateau the height of a low vegetation component? Decompose the dense-cell core.

Andy, 2026-10-18: "Is the plateau telling us the height of the low component of the
vegetation? We are not hitting bare earth, but rather, vegetation at some height above the
ground?"

TEST. Fix a GROUND component's shape from the barest span decile (Gaussian, mu -17.6,
sigma 49.5 mm). In each span decile fit: ground (fixed shape, free weight) + a MAT Gaussian
(free centre, width, weight) + an exponential tail above the mat. If dense cells hold a
sharp ground surface under a mat carrying >~20% of returns, it must appear as a shoulder
with a constrained weight.

RESULT: NOT IDENTIFIABLE. From the second decile up the ground weight collapses to 0.07-0.09
and is unconstrained; a one-component model (single broad core + tail) fits every decile as
well (R2 0.999 vs 1.000). The dense-cell core is ONE hump, width ~130 mm against 50 on bare
ground, with no resolvable narrow ground peak beneath it. Either ground returns are few
(< ~10%) or the surface gen1 sees has no sharp base. Kept as the record that this was tried.
"""
import numpy as np, pyarrow.parquet as pq
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares
from scipy.stats import norm, exponnorm
S="data/derived/elbaext_delong_star_nosurf"; LO,HI,DZ=-2000.,3000.,20.
edges=np.arange(LO,HI+DZ,DZ); ctr=0.5*(edges[:-1]+edges[1:]); NB=ctr.size
sl=np.load(f"{S}/slope.npy"); fp=np.load(f"{S}/floodplain_mask.npy").astype(bool); sp=np.load(f"{S}/gen1_spread_mm.npy"); ny,nx=sl.shape
t=pq.read_table(f"{S}/beam_offset_table.parquet",columns=["cell","d_mm_corr","in_grid"]); g=t["in_grid"].to_numpy().astype(bool)
cell=t["cell"].to_numpy()[g]; d=t["d_mm_corr"].to_numpy()[g].astype(float); ok=np.isfinite(d)&(d>LO)&(d<HI); cell,d=cell[ok],d[ok]
H=np.zeros((ny*nx,NB),np.int32); np.add.at(H,(cell,np.clip(((d-LO)/DZ).astype(int),0,NB-1)),1); H=H.reshape(ny,nx,NB)
flat=fp&(sl<=2)&np.isfinite(sp)&(H.sum(2)>=3)
Hp=np.rint(uniform_filter(H.astype(np.float32),size=(3,3,1),mode="nearest")*9); F=Hp.reshape(-1,NB)[flat.ravel()]; tot=F.sum(1)
C=np.cumsum(F,1); q=lambda p: ctr[np.argmax(C>=p*tot[:,None],axis=1)]
span=q(.9)-q(.1); y=100*F[:,ctr<=0].sum(1)/tot
NBIN=10; bi=np.minimum((np.argsort(np.argsort(span))*NBIN)//span.size,NBIN-1)
P=F/tot[:,None]; pdf=np.array([P[bi==b].mean(0) for b in range(NBIN)])/DZ
# 1. GROUND component shape from the barest decile: one Gaussian + exponential tail
def one(p,x): mu,s,L,f=p; return (1-f)*norm.pdf(x,mu,s)+f*exponnorm.pdf(x,L/s,loc=mu,scale=s)
r0=least_squares(lambda p:(one(p,ctr)-pdf[0])*1e3,[-10,50,150,0.2],bounds=([-300,5,10,0],[300,300,2000,1]))
mu_g,s_g=r0.x[0],r0.x[1]; print(f"GROUND component, fixed from the barest span decile: mu_g {mu_g:.1f} mm, sigma_g {s_g:.1f} mm (fit R2 {1-((one(r0.x,ctr)-pdf[0])**2).sum()/((pdf[0]-pdf[0].mean())**2).sum():.3f})")
# 2. per decile: ground (fixed shape, free weight) + MAT Gaussian (free mu, sigma, weight) + exponential tail above the mat
def two(p,x):
    wg,mu_m,s_m,L,wt=p; wm=max(1-wg-wt,0)
    return wg*norm.pdf(x,mu_g,s_g)+wm*norm.pdf(x,mu_m,s_m)+wt*exponnorm.pdf(x,L/s_m,loc=mu_m,scale=s_m)
def onefree(p,x): return one(p,x)
print(f"\n{'decile':>6}{'span':>6}{'obs pct0':>9} | {'w_ground':>9}{'mat mu':>8}{'mat sig':>8}{'w_mat':>7}{'w_tail':>7}{'L':>6}{'R2 2-comp':>10}{'R2 1-comp':>10} | {'pred pct0':>10}{'ground-only below0':>20}")
for b in range(NBIN):
    r=least_squares(lambda p:(two(p,ctr)-pdf[b])*1e3,[0.5,40,60,150,0.2],bounds=([0,-100,5,10,0],[1,600,300,2000,1]))
    r1=least_squares(lambda p:(one(p,ctr)-pdf[b])*1e3,[-10,60,150,0.3],bounds=([-300,5,10,0],[300,300,2000,1]))
    wg,mu_m,s_m,L,wt=r.x; wm=max(1-wg-wt,0)
    ss=((pdf[b]-pdf[b].mean())**2).sum(); R2=1-((two(r.x,ctr)-pdf[b])**2).sum()/ss; R1=1-((one(r1.x,ctr)-pdf[b])**2).sum()/ss
    m=two(r.x,ctr); pm=m*DZ; pm/=pm.sum(); pred=100*pm[ctr<=0].sum()
    print(f"{b:>6}{np.median(span[bi==b]):>6.0f}{np.median(y[bi==b]):>9.1f} | {wg:>9.2f}{mu_m:>8.0f}{s_m:>8.0f}{wm:>7.2f}{wt:>7.2f}{L:>6.0f}{R2:>10.3f}{R1:>10.3f} | {pred:>10.1f}{100*wg*norm.cdf(0,mu_g,s_g):>20.1f}")
print("\nground-only below0 = the share of ALL returns that are ground-component returns below the gen2 surface (w_ground x 0.63-ish);")
print("the difference from pred pct0 is the mat's lower flank below gen2.")
