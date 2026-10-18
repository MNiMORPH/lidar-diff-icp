#!/usr/bin/env python3
"""Is the lower plateau of pct0(span) real, and can a no-plateau, global-parameter form
replace the logistic? Three tests, recorded together (2026-10-18). NEGATIVE on the second.

Andy: "Denser and denser vegetation should continue to prevent the laser beams from
hitting the surface instead of simply plateauing to a certain number hitting the surface.
This could lead us to a functional form that manages the fit with global parameters only."

TEST 1 (window artifact?): rebuild the histograms with the window widened from +/-1.5 m to
+/-5 m. Only 0.74% of flat-floodplain returns had been clipped; the plateau is unchanged.
Not the window.

TEST 2 (no-plateau forms in span): threshold-exponential, stretched exponential and
power-law, each with global parameters only, under the same 500 m block CV as the
logistic. All score WORSE (0.115-0.135 vs 0.153); exp-to-floor ties the logistic exactly.
Within the observed range the data flatten between spans of 440 and 580 mm, and the top
4% of cells (span > 600) are a different population whose pct0 RISES to 30-42 -- not
denser sedge. The reason span cannot keep registering density is algebraic: for a core
plus an exponential tail, p90 - p10 saturates near core + L ln(10) ~ 570 mm as the tail
fraction -> 1. The proxy saturates; the physics may not.

TEST 3 (a proxy that does not saturate): Beer-Lambert p_hi exp(-x/h0) in the p90 HEIGHT
reaches CV R2 0.656 and keeps falling at the top (12.9 -> 8.1) -- but the ALGEBRAIC NULL
(every cell the same bare-floodplain shape, shifted only to reproduce its own pct0) gives
R2 1.000 for p90 height and 0.99 for the tail integral measured from gen2. Those proxies
are measured relative to the gen2 surface and so READ THE LEVEL: they are pct0 restated.
Every translation-invariant candidate (p90 - mode; tail length as mean excess above the
mode; tail mass above the mode; their product) scores <= 0.13, below span, and none keeps
falling in its top decile.

CONCLUSION. gen1's own returns carry no vegetation-density measure that is independent of
the gen2 level AND keeps registering past where span saturates. Beer-Lambert P = exp(-k tau)
is the right physics, but tau must come from somewhere other than gen1's cloud, and the
cloud's own P is what pct0 already measures (pct0 ~ 0.63 P). The logistic in span stays as
the description; its floor is the floor of the OBSERVED range, not a claim about denser
vegetation than this floodplain holds.

Three sections follow, as run; each reads the beam table directly and is meant to be
executed one section at a time (copy a section out, or run the file and read all three).
"""

# ============================== TEST 1 and 2: window, no-plateau forms in span
import numpy as np, pyarrow.parquet as pq
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares
S="data/derived/elbaext_delong_star_nosurf"
sl=np.load(f"{S}/slope.npy"); fp=np.load(f"{S}/floodplain_mask.npy").astype(bool); sp=np.load(f"{S}/gen1_spread_mm.npy"); ny,nx=sl.shape
t=pq.read_table(f"{S}/beam_offset_table.parquet",columns=["cell","d_mm_corr","in_grid"]); g=t["in_grid"].to_numpy().astype(bool)
cell=t["cell"].to_numpy()[g]; d=t["d_mm_corr"].to_numpy()[g].astype(float); fin=np.isfinite(d)
print(f"returns with |d| > 1500 mm on flat floodplain cells (clipped by the old window): ", end="")
def build(LO,HI,DZ):
    edges=np.arange(LO,HI+DZ,DZ); ctr=0.5*(edges[:-1]+edges[1:]); NB=ctr.size
    ok=fin&(d>LO)&(d<HI); H=np.zeros((ny*nx,NB),np.int32); np.add.at(H,(cell[ok],np.clip(((d[ok]-LO)/DZ).astype(int),0,NB-1)),1)
    return H.reshape(ny,nx,NB),ctr,NB
H,ctr,NB=build(-5000,5000,20)
flat=fp&(sl<=2)&np.isfinite(sp)&(H.sum(2)>=3)
fcell=np.zeros(ny*nx,bool); fcell[np.flatnonzero(flat.ravel())]=True
onfp=fcell[cell]&fin; print(f"{np.mean(np.abs(d[onfp])>1500)*100:.2f}% of {onfp.sum():,}")
iy,ix=np.nonzero(flat); block=(iy//100)*10000+ix//100; ub,binv=np.unique(block,return_inverse=True)
fold=np.random.default_rng(0).permutation(np.arange(ub.size)%5)[binv]
Hp=np.rint(uniform_filter(H.astype(np.float32),size=(3,3,1),mode="nearest")*9); F=Hp.reshape(-1,NB)[flat.ravel()]; tot=F.sum(1)
C=np.cumsum(F,1); q=lambda p: ctr[np.argmax(C>=p*tot[:,None],axis=1)]
x=q(.9)-q(.1); y=100*F[:,ctr<=0].sum(1)/tot; p90=q(.9); p10=q(.1)
bi=np.minimum((np.argsort(np.argsort(x))*20)//x.size,19)
print("\n3x3, window widened to +/-5 m. Per bin: median span, pct0, p10, p90, and share of cells whose p90 exceeds +1.5 m")
print(f"  {'span':>6}{'pct0':>7}{'p10':>6}{'p90':>6}{'p90>1.5m':>10}{'n':>6}")
for b in range(20):
    m=bi==b; print(f"  {np.median(x[m]):>6.0f}{np.median(y[m]):>7.1f}{np.median(p10[m]):>6.0f}{np.median(p90[m]):>6.0f}{100*np.mean(p90[m]>1500):>9.1f}%{m.sum():>6}")
# finer look at the top: 10 bins within the top decile of span
top=x>=np.percentile(x,90); xt=x[top]; yt=y[top]; bt=np.minimum((np.argsort(np.argsort(xt))*10)//xt.size,9)
print("\n  top decile of span, 10 sub-bins:")
for b in range(10):
    m=bt==b; print(f"  {np.median(xt[m]):>6.0f}{np.median(yt[m]):>7.1f}{'':>12}{'':>10}{m.sum():>6}")
# candidate global-parameter forms
forms={
 "logistic (4)":      (lambda p,x: p[1]+(p[0]-p[1])/(1+np.exp(np.clip((x-p[2])/p[3],-50,50))), [70,20,250,60], ([0,0,0,1],[100,100,3000,2000])),
 "thresh-exp (3)":    (lambda p,x: p[0]*np.exp(-np.maximum(x-p[1],0)/p[2]), [60,150,250], ([0,0,10],[100,2000,5000])),
 "stretched-exp (3)": (lambda p,x: p[0]*np.exp(-(x/p[1])**p[2]), [60,400,2], ([0,10,0.2],[100,5000,10])),
 "exp-to-floor (4)":  (lambda p,x: p[1]+(p[0]-p[1])*np.exp(-np.maximum(x-p[2],0)/p[3]), [60,20,150,200], ([0,0,0,10],[100,100,2000,5000])),
 "power (3)":         (lambda p,x: p[0]*(1+x/p[1])**(-p[2]), [70,200,2], ([0,1,0.1],[100,5000,20])),
}
def fit(f,x0,bd,x,y): return least_squares(lambda p:f(p,x)-y,x0=x0,bounds=bd,loss="soft_l1",f_scale=10).x
print(f"\n  {'form':<20}{'CV R2':>7}{'NMAD':>6}  params (fit on all)            pred at span 400 / 600 / 900")
for nm,(f,x0,bd) in forms.items():
    pred=np.empty_like(y)
    for k in range(5):
        tr=fold!=k; pred[~tr]=f(fit(f,x0,bd,x[tr],y[tr]),x[~tr])
    res=y-pred; r2=1-res.var()/y.var(); nm_=1.4826*np.median(np.abs(res-np.median(res)))
    p=fit(f,x0,bd,x,y)
    print(f"  {nm:<20}{r2:>7.3f}{nm_:>6.1f}  {np.array2string(p,precision=1,max_line_width=200):<30} {f(p,400):>5.1f} / {f(p,600):>5.1f} / {f(p,900):>5.1f}")

# ============================== TEST 3a: level-dependent proxies (Beer-Lambert in p90 etc.)
import numpy as np, pyarrow.parquet as pq
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares
S="data/derived/elbaext_delong_star_nosurf"; LO,HI,DZ=-5000.,5000.,20.
edges=np.arange(LO,HI+DZ,DZ); ctr=0.5*(edges[:-1]+edges[1:]); NB=ctr.size
sl=np.load(f"{S}/slope.npy"); fp=np.load(f"{S}/floodplain_mask.npy").astype(bool); sp=np.load(f"{S}/gen1_spread_mm.npy"); ny,nx=sl.shape
t=pq.read_table(f"{S}/beam_offset_table.parquet",columns=["cell","d_mm_corr","in_grid"]); g=t["in_grid"].to_numpy().astype(bool)
cell=t["cell"].to_numpy()[g]; d=t["d_mm_corr"].to_numpy()[g].astype(float); ok=np.isfinite(d)&(d>LO)&(d<HI); cell,d=cell[ok],d[ok]
H=np.zeros((ny*nx,NB),np.int32); np.add.at(H,(cell,np.clip(((d-LO)/DZ).astype(int),0,NB-1)),1); H=H.reshape(ny,nx,NB)
flat=fp&(sl<=2)&np.isfinite(sp)&(H.sum(2)>=3)
iy,ix=np.nonzero(flat); block=(iy//100)*10000+ix//100; ub,binv=np.unique(block,return_inverse=True)
fold=np.random.default_rng(0).permutation(np.arange(ub.size)%5)[binv]
Hp=np.rint(uniform_filter(H.astype(np.float32),size=(3,3,1),mode="nearest")*9); F=Hp.reshape(-1,NB)[flat.ravel()]; tot=F.sum(1)
C=np.cumsum(F,1); q=lambda p: ctr[np.argmax(C>=p*tot[:,None],axis=1)]
span=q(.9)-q(.1); p90=q(.9); y=100*F[:,ctr<=0].sum(1)/tot
above=ctr>100.0; na=F[:,above].sum(1); hbar=np.where(na>0,(F[:,above]*ctr[above]).sum(1)/np.maximum(na,1),0.0)   # mean height of returns >+100 mm
fabove=na/tot                                                   # tail mass above +100 mm (NOT independent of y)
hmass=hbar*fabove                                               # height x mass: a column-amount proxy (tail integral above +100)
proxies={"span (p90-p10)":span,"p90 height":np.maximum(p90,0),"mean height of returns >100 mm":hbar,"tail integral: mean height x mass above 100 mm":hmass}
forms={
 "logistic (4)": (lambda p,x: p[1]+(p[0]-p[1])/(1+np.exp(np.clip((x-p[2])/p[3],-50,50))), [65,20,250,60], ([0,0,0,1],[100,100,5000,5000])),
 "Beer-Lambert p*exp(-x/h0) (2)": (lambda p,x: p[0]*np.exp(-x/p[1]), [65,300], ([0,1],[100,20000])),
 "BL with threshold (3)": (lambda p,x: p[0]*np.exp(-np.maximum(x-p[1],0)/p[2]), [65,50,300], ([0,0,1],[100,5000,20000])),
}
def fit(f,x0,bd,x,y): return least_squares(lambda p:f(p,x)-y,x0=x0,bounds=bd,loss="soft_l1",f_scale=10).x
def cvscore(f,x0,bd,x,y,sel):
    pred=np.full_like(y,np.nan)
    for k in range(5):
        tr=sel&(fold!=k); te=sel&(fold==k); pred[te]=f(fit(f,x0,bd,x[tr],y[tr]),x[te])
    res=(y-pred)[sel]; return 1-res.var()/y[sel].var(), 1.4826*np.median(np.abs(res-np.median(res)))
mixed=span>600
print(f"3x3, window +/-5 m, {flat.sum():,} cells; 'mixed' = span > 600 mm: {mixed.sum():,} cells ({100*mixed.mean():.1f}%), median pct0 there {np.median(y[mixed]):.1f}")
for sel_name,sel in (("ALL cells",np.ones_like(mixed)),("EXCLUDING span>600 (sensitivity only)",~mixed)):
    print(f"\n== {sel_name} ==")
    print(f"  {'proxy x':<48}{'form':<30}{'CV R2':>7}{'NMAD':>6}  params")
    for pn,x in proxies.items():
        for fn,(f,x0,bd) in forms.items():
            r2,nm=cvscore(f,x0,bd,x,y,sel); p=fit(f,x0,bd,x[sel],y[sel])
            print(f"  {pn:<48}{fn:<30}{r2:>7.3f}{nm:>6.1f}  {np.array2string(p,precision=1,max_line_width=200)}", flush=True)
# does the proxy keep registering density where span saturates? binned medians of pct0 vs each proxy, top decile only
print("\nTop decile of each proxy, 5 sub-bins: does pct0 keep falling?")
for pn,x in proxies.items():
    top=x>=np.percentile(x,90); xt=x[top]; yt=y[top]; bt=np.minimum((np.argsort(np.argsort(xt))*5)//xt.size,4)
    print(f"  {pn:<48}"+"  ".join(f"x {np.median(xt[bt==b]):.0f}: {np.median(yt[bt==b]):.1f}" for b in range(5)))

# ============================== TEST 3b: translation-invariant proxies and the ALGEBRAIC NULL
import numpy as np, pyarrow.parquet as pq, sys
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares
sys.path.insert(0,"analysis"); from hsm import hsm_from_hist
S="data/derived/elbaext_delong_star_nosurf"; LO,HI,DZ=-5000.,5000.,20.
edges=np.arange(LO,HI+DZ,DZ); ctr=0.5*(edges[:-1]+edges[1:]); NB=ctr.size
sl=np.load(f"{S}/slope.npy"); fp=np.load(f"{S}/floodplain_mask.npy").astype(bool); sp=np.load(f"{S}/gen1_spread_mm.npy"); ny,nx=sl.shape
t=pq.read_table(f"{S}/beam_offset_table.parquet",columns=["cell","d_mm_corr","in_grid"]); g=t["in_grid"].to_numpy().astype(bool)
cell=t["cell"].to_numpy()[g]; d=t["d_mm_corr"].to_numpy()[g].astype(float); ok=np.isfinite(d)&(d>LO)&(d<HI); cell,d=cell[ok],d[ok]
H=np.zeros((ny*nx,NB),np.int32); np.add.at(H,(cell,np.clip(((d-LO)/DZ).astype(int),0,NB-1)),1); H=H.reshape(ny,nx,NB)
flat=fp&(sl<=2)&np.isfinite(sp)&(H.sum(2)>=3)
iy,ix=np.nonzero(flat); block=(iy//100)*10000+ix//100; ub,binv=np.unique(block,return_inverse=True)
fold=np.random.default_rng(0).permutation(np.arange(ub.size)%5)[binv]
Hp=np.rint(uniform_filter(H.astype(np.float32),size=(3,3,1),mode="nearest")*9); F=Hp.reshape(-1,NB)[flat.ravel()]; tot=F.sum(1)
C=np.cumsum(F,1); q=lambda p: ctr[np.argmax(C>=p*tot[:,None],axis=1)]
span=q(.9)-q(.1); p90=q(.9); p10=q(.1); y=100*F[:,ctr<=0].sum(1)/tot
mode=np.array([hsm_from_hist(F[i],ctr) for i in range(F.shape[0])])
# translation-invariant, non-saturating candidates
thr=mode+100.0
mask=ctr[None,:]>thr[:,None]; na=(F*mask).sum(1)
Lx=np.where(na>0,((F*mask)*(ctr[None,:]-thr[:,None])).sum(1)/np.maximum(na,1),0.0)    # mean excess above mode+100: exponential tail's L (MLE)
p90m=p90-mode                                                                          # p90 measured from the cell's own mode
fab=na/tot                                                                             # share of returns above mode+100 (tail mass, from the mode)
Lmass=Lx*fab                                                                           # tail integral measured from the mode
forms={"logistic (4)":(lambda p,x: p[1]+(p[0]-p[1])/(1+np.exp(np.clip((x-p[2])/p[3],-50,50))),[65,20,250,60],([0,0,0,1],[100,100,5000,5000])),
       "Beer-Lambert (2)":(lambda p,x: p[0]*np.exp(-x/p[1]),[65,300],([0,1],[100,20000]))}
def fit(f,x0,bd,x,y): return least_squares(lambda p:f(p,x)-y,x0=x0,bounds=bd,loss="soft_l1",f_scale=10).x
def cvs(f,x0,bd,x,y):
    pred=np.empty_like(y)
    for k in range(5):
        tr=fold!=k; pred[~tr]=f(fit(f,x0,bd,x[tr],y[tr]),x[~tr])
    r=y-pred; return 1-r.var()/y.var(),1.4826*np.median(np.abs(r-np.median(r)))
print("== TRANSLATION-INVARIANT proxies (measured from the cell's own mode, not from gen2) ==")
print(f"  {'proxy':<44}{'form':<18}{'CV R2':>7}{'NMAD':>6}  params")
for pn,x in (("span p90-p10",span),("p90 - mode",p90m),("tail length L: mean excess above mode+100",Lx),("tail mass above mode+100 (share)",fab),("tail integral from mode: L x mass",Lmass)):
    for fn,(f,x0,bd) in forms.items():
        r2,nm=cvs(f,x0,bd,x,y); p=fit(f,x0,bd,x,y)
        print(f"  {pn:<44}{fn:<18}{r2:>7.3f}{nm:>6.1f}  {np.array2string(p,precision=2,max_line_width=200)}",flush=True)
    top=x>=np.percentile(x,90); xt=x[top]; yt=y[top]; bt=np.minimum((np.argsort(np.argsort(xt))*5)//xt.size,4)
    print(f"  {'':<44}top decile: "+"  ".join(f"x {np.median(xt[bt==b]):.0f}: {np.median(yt[bt==b]):.1f}" for b in range(5)))
# ALGEBRAIC NULL for the level-dependent proxies: one fixed shape, shifted per cell so that pct0 matches the observed value
print("\n== ALGEBRAIC NULL: every cell = the bare-floodplain average shape, SHIFTED to reproduce its observed pct0. ==")
print("   If a proxy's R2 here is as high as in the data, the proxy is reading the level, not the vegetation.")
bare=span<np.percentile(span,25); shape=(F[bare]/tot[bare,None]).mean(0); shape/=shape.sum(); Cs=np.cumsum(shape)
fine=np.arange(-5000,5000,1.0); Cf=np.interp(fine,ctr,Cs)
u=np.array([fine[np.searchsorted(Cf,max(min(v/100,Cf[-1]-1e-9),Cf[0]+1e-9))] for v in y])   # shift u: CDF(0 - u) = pct0  -> surface sits at -u in the shape
# with the shape shifted by s = -u... define per-cell cloud = shape translated so that P(d<=0)=pct0: translate by tau where CDF_shape(-tau)=pct0 -> tau=-u
tau=-u
null={"p90 height (null)":ctr[np.argmax(Cs>=.9)]+tau,"mean height >100 mm from gen2 (null)":None,"span (null)":np.full_like(tau,ctr[np.argmax(Cs>=.9)]-ctr[np.argmax(Cs>=.1)])}
# mean height of returns above +100 mm (relative to gen2) under the null, by direct integration of the shifted shape
h=np.array([((shape*((ctr+tt)>100))*(ctr+tt)).sum()/max((shape*((ctr+tt)>100)).sum(),1e-12) for tt in tau]); null["mean height >100 mm from gen2 (null)"]=h
mass=np.array([(shape*((ctr+tt)>100)).sum() for tt in tau]); null["tail integral from gen2 (null)"]=h*mass
for pn,x in null.items():
    if np.std(x)==0: print(f"  {pn:<44}constant under the null (translation-invariant): R2 undefined"); continue
    for fn,(f,x0,bd) in forms.items():
        r2,nm=cvs(f,x0,bd,x,y); print(f"  {pn:<44}{fn:<18}{r2:>7.3f}{nm:>6.1f}",flush=True)
