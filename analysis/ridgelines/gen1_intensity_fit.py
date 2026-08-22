#!/usr/bin/env python3
"""Model gen1 per-cell ground-return INTENSITY as a linear function of scan angle and
ground-return fraction (penetration), to see how deterministic it is and to separate the
geometry part (angle) from the canopy part (fraction-at-fixed-angle).

Caveat: angle and fraction are ~-0.84 correlated (the source of the 0.6 split), so joint
coefficients are entangled; we report nested models to see what each ADDS.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_intensity_fit.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy")
N=NY*NX; tot=np.zeros(N); gnd=np.zeros(N); isum=np.zeros(N); angsum=np.zeros(N)
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); cl=np.asarray(pts.classification)
        inten=np.asarray(pts.intensity).astype(float); ang=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7); cell=iy[keep]*NX+ix[keep]
        clk=cl[keep]; ik=inten[keep]; ak=ang[keep]; g=clk==2
        np.add.at(tot,cell,1); np.add.at(gnd,cell,g.astype(float)); np.add.at(isum,cell,np.where(g,ik,0.0))
        np.add.at(angsum,cell,ak)
pen=np.where(tot>0,gnd/np.maximum(tot,1),np.nan)
gint=np.where(gnd>0,isum/np.maximum(gnd,1),np.nan)
ang=np.where(tot>0,angsum/np.maximum(tot,1),np.nan)

def fit(y,cols,names):
    X=np.column_stack([np.ones_like(y)]+cols); beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    yh=X@beta; r2=1-np.sum((y-yh)**2)/np.sum((y-np.mean(y))**2); rmse=np.sqrt(np.mean((y-yh)**2))
    return beta,r2,rmse,yh

for msk,lbl in [((~fld.ravel())&np.isfinite(pen)&(gnd>=5),"ALL cells"),
                ((g2pen.ravel()<0.25)&(~fld.ravel())&np.isfinite(pen)&(gnd>=5),"FOREST cells")]:
    p=pen[msk]; a=ang[msk]; I=gint[msk]
    print(f"\n=== gen1 ground intensity fit — {lbl} (n={msk.sum()}) ===")
    print(f"  ranges: intensity {I.min():.0f}-{I.max():.0f} (mean {I.mean():.1f}); angle {a.min():.0f}-{a.max():.0f}; pen {p.min():.2f}-{p.max():.2f}")
    print(f"  corr(angle,pen)={np.corrcoef(a,p)[0,1]:+.2f}  corr(I,angle)={np.corrcoef(I,a)[0,1]:+.2f}  corr(I,pen)={np.corrcoef(I,p)[0,1]:+.2f}")
    for cols,names in [([p],["pen"]),([a],["ang"]),([a,p],["ang","pen"]),
                       ([a,p,a*p],["ang","pen","ang*pen"])]:
        b,r2,rmse,_=fit(I,cols,names)
        terms=" ".join(f"{names[i]}:{b[i+1]:+.2f}" for i in range(len(names)))
        print(f"  I ~ {'+'.join(names):18s}: R2={r2:.3f} RMSE={rmse:.1f}   b0={b[0]:+.1f} {terms}")

# figure for FOREST: predicted vs actual for the ang+pen model
msk=(g2pen.ravel()<0.25)&(~fld.ravel())&np.isfinite(pen)&(gnd>=5)
p=pen[msk]; a=ang[msk]; I=gint[msk]
b,r2,rmse,yh=fit(I,[a,p],["ang","pen"])
fig,ax=plt.subplots(1,2,figsize=(13,6))
hb=ax[0].hexbin(yh,I,gridsize=45,bins="log",cmap="viridis",mincnt=1); lims=[I.min(),np.percentile(I,99)]
ax[0].plot(lims,lims,"r--"); ax[0].set_xlim(lims); ax[0].set_ylim(lims)
ax[0].set_xlabel("predicted intensity (ang+pen)"); ax[0].set_ylabel("actual intensity"); ax[0].set_title(f"forest: I ~ ang+pen  (R2={r2:.2f}, RMSE={rmse:.1f})"); fig.colorbar(hb,ax=ax[0])
# intensity vs angle, colored by pen
sc=ax[1].hexbin(a,I,C=p,gridsize=40,cmap="plasma",mincnt=1); ax[1].set_xlabel("|scan angle| (deg)"); ax[1].set_ylabel("ground intensity")
ax[1].set_title("intensity vs angle, colored by penetration"); fig.colorbar(sc,ax=ax[1],label="penetration")
fig.suptitle("gen1 ground intensity vs scan angle + ground-return fraction (forest)",y=1.0)
fig.savefig("figures/refdatum/gen1_intensity_fit.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_intensity_fit.png")
