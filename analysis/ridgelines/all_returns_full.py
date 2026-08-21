#!/usr/bin/env python3
"""ALL returns, NO vertical window — full slope-normal profile ground->canopy-top, both epochs,
per land-cover stratum.  Same transform as before (plane = gen2 bare earth).  1 cm bins over
-3..60 m so nothing is truncated.  Saves all_returns_full.npz; plots the full column (log density).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/all_returns_full.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; GEN2="data/after/3dep2021_fulldensity.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zf=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
pen=np.load("data/derived/elba_fulldensity/penetration.npy"); fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
strat=np.zeros(NY*NX,np.int8)
strat[((pen<0.25)&~fld&np.isfinite(pen)).ravel()]=1; strat[((pen>=0.45)&~fld&np.isfinite(pen)).ravel()]=2

FLO,FHI,FW=-3.0,60.0,0.01; fedges=np.arange(FLO,FHI+0.5*FW,FW); fc=0.5*(fedges[:-1]+fedges[1:]); NF=fc.size
H={(g,s):np.zeros(NF) for g in("g1","g2") for s in(1,2)}
def acc(path,g):
    n=0
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
            cl=np.asarray(pts.classification)
            ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
            keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
            ix=ix[keep];iy=iy[keep];x=x[keep];y=y[keep];z=z[keep]; n+=keep.sum()
            cell=iy*NX+ix; xc=X0+(ix+0.5)*RES; yc=Y0+(iy+0.5)*RES
            d=(z-(Zf[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*cosf[cell]
            st=strat[cell]; fb=np.searchsorted(fedges,d,side="right")-1; ok=(fb>=0)&(fb<NF)
            for s in(1,2):
                ss=ok&(st==s)
                if ss.any(): np.add.at(H[(g,s)],fb[ss],1)
    print(f"  {g}: kept {n:,}")
print("streaming gen1 (all)..."); acc(GEN1,"g1")
print("streaming gen2 (all)..."); acc(GEN2,"g2")
for g in("g1","g2"):
    for s in(1,2): print(f"  {g} strat{s}: {H[(g,s)].sum():,.0f} returns, max d with returns = {fc[np.max(np.nonzero(H[(g,s)]))]:.1f} m")
np.savez_compressed("data/derived/elba_fulldensity/all_returns_full.npz",fedges=fedges,fc=fc,
                    **{f"{g}_{s}":H[(g,s)] for g in("g1","g2") for s in(1,2)})

# ---- plot: full column, log density, forest + open ----------------------------------
fig,axes=plt.subplots(1,2,figsize=(13,8),sharey=True)
for ax,(s,lbl) in zip(axes,[(1,"FOREST"),(2,"OPEN")]):
    for g,col,lab in [("g1","C0","gen1 (2008 leaf-off) all"),("g2","C3","gen2 (2021 leaf-on) all")]:
        c=H[(g,s)]; d=c/c.sum()/FW; ax.semilogx(np.where(d>0,d,np.nan),fc,col,lw=1.3,label=lab)
    ax.axhline(0,color="k",lw=.5); ax.set_ylim(-1,40); ax.set_xlim(1e-5,3e1)
    ax.set_xlabel("density (1/m, log)"); ax.set_title(f"{lbl}: ALL returns, full column"); ax.legend(fontsize=9); ax.grid(alpha=.3,which="both")
axes[0].set_ylabel("slope-normal d (m)  [plane = gen2 bare earth]")
fig.suptitle("ALL returns, no vertical window — full ground-to-canopy profile",y=0.98)
fig.savefig("figures/refdatum/all_returns_full.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("wrote figures/refdatum/all_returns_full.png")
