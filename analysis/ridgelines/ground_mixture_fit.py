#!/usr/bin/env python3
"""Reconstruct the gen2 ground-class return distribution as a two-component mixture:
    (1) a GROUND-surface Gaussian  N(mu_g, sigma_g)   -- the true bare-earth return
    (2) a PLANT component above ground (leaf-on groundcover/litter that leaked into the
        ground class), tried as (a) a second Gaussian and (b) an upward-decaying
        exponentially-modified Gaussian (EMG) sharing the measurement sigma.

Goal (Andy): recover the LIKELY GROUND LEVEL and its statistic for extraction -- i.e. the
ground-Gaussian mean mu_g, which sits BELOW the raw ground-class median (the median is
pulled up by the plant shoulder).  Report which raw percentile mu_g corresponds to, so the
pipeline can pick a matched ground estimator.

Resolution: the saved histograms are 0.25 m (too coarse for sigma_g~0.1 m).  This script
re-streams the class-2 clouds and accumulates POOLED 1 cm histograms per land-cover
stratum, using the SAME slope-normal transform as slope_normal_returns.py (copied verbatim).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/ground_mixture_fit.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
from scipy.special import erfc, logsumexp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ---- grid + transform (verbatim from slope_normal_returns.py) -----------------------
NY, NX = 700, 508; X0, Y0 = 577492.8, 4882737.6; RES = 5.0
GEN2 = "data/after/3dep2021_fulldensity.laz"; GEN1 = "data/before/4342-29-64.laz"
CHUNK = 5_000_000
Zg = np.load("data/derived/elba_fulldensity/z_after.npy")
Zg_filled = Zg.copy(); m = ~np.isfinite(Zg_filled)
if m.any():
    Zg_filled = Zg_filled[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zg_filled, RES)
cos_slope = 1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zg_flat = Zg_filled.ravel(); gx_flat = gx.ravel(); gy_flat = gy.ravel(); cos_flat = cos_slope.ravel()

# ---- per-cell stratum: 0 none, 1 forest, 2 open (same defs as decomposition) --------
pen = np.load("data/derived/elba_fulldensity/penetration.npy")
fld = np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
strat = np.zeros(NY*NX, np.int8)
sf = (pen < 0.25) & ~fld & np.isfinite(pen)
so = (pen >= 0.45) & ~fld & np.isfinite(pen)
strat[sf.ravel()] = 1; strat[so.ravel()] = 2

# ---- fine 1 cm pooled histograms, class 2 only, per stratum -------------------------
FLO, FHI, FW = -0.8, 1.6, 0.01
fedges = np.arange(FLO, FHI + 0.5*FW, FW); fc = 0.5*(fedges[:-1]+fedges[1:]); NF = fc.size
H = {("gen1",1): np.zeros(NF), ("gen1",2): np.zeros(NF),
     ("gen2",1): np.zeros(NF), ("gen2",2): np.zeros(NF)}

def accumulate(path, gen):
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            cl = np.asarray(pts.classification)
            g = cl == 2
            if not g.any(): continue
            x = np.asarray(pts.x, np.float64)[g]; y = np.asarray(pts.y, np.float64)[g]
            z = np.asarray(pts.z, np.float64)[g]
            ix = ((x-X0)/RES).astype(np.int64); iy = ((y-Y0)/RES).astype(np.int64)
            keep = (ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)
            ix=ix[keep]; iy=iy[keep]; x=x[keep]; y=y[keep]; z=z[keep]
            cell = iy*NX+ix
            xc = X0+(ix+0.5)*RES; yc = Y0+(iy+0.5)*RES
            Zp = Zg_flat[cell] + gx_flat[cell]*(x-xc) + gy_flat[cell]*(y-yc)
            d = (z - Zp)*cos_flat[cell]
            st = strat[cell]
            fb = np.searchsorted(fedges, d, side="right")-1
            ok = (fb>=0)&(fb<NF)
            for s in (1,2):
                sel = ok & (st==s)
                if sel.any():
                    np.add.at(H[(gen,s)], fb[sel], 1)
    print(f"  {gen}: forest={H[(gen,1)].sum():,.0f}  open={H[(gen,2)].sum():,.0f} class-2 returns pooled")

print("streaming class-2 clouds (1 cm pooled)...")
accumulate(GEN2, "gen2"); accumulate(GEN1, "gen1")
np.savez_compressed("data/derived/elba_fulldensity/ground_fine_pooled.npz",
                    fedges=fedges, fc=fc, **{f"{g}_{s}": H[(g,s)] for (g,s) in H})

# ---- mixture models ------------------------------------------------------------------
def gauss(d, mu, sig): return np.exp(-0.5*((d-mu)/sig)**2)/(sig*np.sqrt(2*np.pi))
def emg(d, mu, sig, tau):   # right-tailed exponentially-modified Gaussian, mean mu+tau
    z = (d-mu)/sig - sig/tau
    return (0.5/tau)*np.exp(0.5*(sig/tau)**2 - (d-mu)/tau)*erfc(-z/np.sqrt(2))
def sp(x): return np.log1p(np.exp(-abs(x))) + max(x,0)  # softplus, stable

def nll_2gauss(p, d, c):
    mu, ls, off, lsp, lw = p
    sg=np.exp(ls); sp2=np.exp(lsp); w=1/(1+np.exp(-lw))
    f = w*gauss(d,mu,sg) + (1-w)*gauss(d, mu+np.exp(off), sp2)
    f = np.clip(f,1e-12,None)*FW
    return -(c*np.log(f)).sum()
def nll_emg(p, d, c):
    mu, ls, lt, lw = p
    sg=np.exp(ls); tau=np.exp(lt); w=1/(1+np.exp(-lw))
    f = w*gauss(d,mu,sg) + (1-w)*emg(d,mu,sg,tau)
    f = np.clip(f,1e-12,None)*FW
    return -(c*np.log(f)).sum()

def fit(counts, tag):
    c = counts.copy(); d = fc
    mode = fc[np.argmax(c)]
    # 2-Gaussian
    r1 = minimize(nll_2gauss, [mode, np.log(0.08), np.log(0.15), np.log(0.15), 0.5],
                  args=(d,c), method="Nelder-Mead", options=dict(maxiter=20000,xatol=1e-5,fatol=1e-3))
    mu,ls,off,lsp,lw = r1.x; sg=np.exp(ls); sp2=np.exp(lsp); w=1/(1+np.exp(-lw)); mup=mu+np.exp(off)
    k1=5; aic1=2*k1+2*r1.fun
    # EMG-plant
    r2 = minimize(nll_emg, [mode, np.log(0.08), np.log(0.15), 0.5],
                  args=(d,c), method="Nelder-Mead", options=dict(maxiter=20000,xatol=1e-5,fatol=1e-3))
    mu2,ls2,lt2,lw2 = r2.x; sg2=np.exp(ls2); tau=np.exp(lt2); w2=1/(1+np.exp(-lw2))
    k2=4; aic2=2*k2+2*r2.fun
    # raw stats
    cdf = np.cumsum(c)/c.sum()
    raw_med = np.interp(0.5, cdf, d)
    # percentile of the fitted ground mean in the raw distribution
    pct_mu = np.interp(mu, d, cdf)*100
    print(f"\n=== {tag} (n={c.sum():,.0f}) ===")
    print(f"  raw: mode {mode:+.3f}  median {raw_med:+.3f} m")
    print(f"  [2-Gaussian] ground mu_g={mu:+.4f}  sig_g={sg:.4f}  w_ground={w:.2f}  | "
          f"plant mu_p={mup:+.3f} sig_p={sp2:.3f}   AIC={aic1:.0f}")
    print(f"  [Gauss+EMG ] ground mu_g={mu2:+.4f}  sig_g={sg2:.4f}  w_ground={w2:.2f}  | "
          f"plant tau={tau:.3f} (mean +{tau:.2f} m)   AIC={aic2:.0f}   {'<-better' if aic2<aic1 else ''}")
    print(f"  EXTRACTION: ground level mu_g = {mu:+.4f} m (2G) / {mu2:+.4f} m (EMG); "
          f"= raw p{pct_mu:.0f}.  Raw median is +{(raw_med-mu)*1000:.0f} mm above true ground.")
    best = ("2G",mu,sg,w,mup,sp2,None) if aic1<=aic2 else ("EMG",mu2,sg2,w2,None,None,tau)
    return dict(mode=mode, raw_med=raw_med, mu=mu, sg=sg, w=w, mup=mup, sp=sp2,
                mu2=mu2, sg2=sg2, w2=w2, tau=tau, aic1=aic1, aic2=aic2, pct_mu=pct_mu, best=best)

res = {}
for s,lbl in [(1,"gen2 FOREST"),(2,"gen2 OPEN")]:
    res[s] = fit(H[("gen2",s)], lbl)

# ---- figure: fit overlay, gen2 forest + open ----------------------------------------
fig, axes = plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,s,lbl in [(axes[0],1,"gen2 FOREST"),(axes[1],2,"gen2 OPEN")]:
    c = H[("gen2",s)]; dens = c/c.sum()/FW; R=res[s]
    ax.plot(dens, fc, "0.4", lw=1.2, label="gen2 ground-class (1 cm)")
    g = R["w"]*gauss(fc,R["mu"],R["sg"]); pl = (1-R["w"])*gauss(fc,R["mup"],R["sp"])
    ax.plot(g, fc, "C0", lw=1.8, label=f"ground N({R['mu']:+.2f},{R['sg']:.2f})")
    ax.plot(pl, fc, "C2", lw=1.5, label="plant (2nd Gaussian)")
    ax.plot(g+pl, fc, "C3--", lw=1.3, label="mixture")
    ax.axhline(R["mu"], color="C0", ls=":", lw=1)
    ax.axhline(R["raw_med"], color="k", ls=":", lw=1, label=f"raw median {R['raw_med']:+.2f}")
    ax.set_ylim(-0.6,1.2); ax.set_xlabel("density (1/m)"); ax.set_title(lbl); ax.legend(fontsize=8)
axes[0].set_ylabel("slope-normal d (m)")
fig.suptitle("gen2 ground-class = ground Gaussian + plant component; mu_g = extracted ground")
fig.savefig("figures/refdatum/ground_mixture_fit.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/ground_mixture_fit.png and data/.../ground_fine_pooled.npz")
