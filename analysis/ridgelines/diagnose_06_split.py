#!/usr/bin/env python3
"""Diagnose the split at ~0.6 gen1 ground-return fraction in the penetration-vs-intensity cloud.
Per cell, collect the candidates that could create two populations:
  - RETURN STRUCTURE: single(nr=1) vs multi(nr>1) pulse mix -> sets hard bounds on the fraction
    (all-single -> frac 1.0; all-double canopy+ground -> frac 0.5; triples push below 0.5).
  - point density (tot returns / cell), flight-line overlap (# distinct point_source_id), scan angle.
Then: histogram of penetration (is it bimodal at 0.6?) and what differs across the 0.6 boundary.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/diagnose_06_split.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)

N=NY*NX
tot=np.zeros(N); gnd=np.zeros(N); isum=np.zeros(N)
gnd_single=np.zeros(N); gnd_multi=np.zeros(N)   # ground returns from single vs multi-return pulses
allret_multi=np.zeros(N); angsum=np.zeros(N)
psid=[dict() for _ in range(0)]                  # skip per-cell psid dict (memory); track globally below
psid_sets=None
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); cl=np.asarray(pts.classification)
        inten=np.asarray(pts.intensity).astype(float); nr=np.asarray(pts.number_of_returns)
        ang=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
        cell=iy[keep]*NX+ix[keep]; clk=cl[keep]; ik=inten[keep]; nrk=nr[keep]; ak=ang[keep]
        g=clk==2
        np.add.at(tot,cell,1); np.add.at(gnd,cell,g.astype(float)); np.add.at(isum,cell,np.where(g,ik,0.0))
        np.add.at(gnd_single,cell,(g&(nrk==1)).astype(float)); np.add.at(gnd_multi,cell,(g&(nrk>1)).astype(float))
        np.add.at(allret_multi,cell,(nrk>1).astype(float)); np.add.at(angsum,cell,ak)
pen=np.where(tot>0,gnd/np.maximum(tot,1),np.nan)
gint=np.where(gnd>0,isum/np.maximum(gnd,1),np.nan)
multifrac=np.where(tot>0,allret_multi/np.maximum(tot,1),np.nan)   # fraction of ALL returns that are multi-return
gsingfrac=np.where(gnd>0,gnd_single/np.maximum(gnd,1),np.nan)     # fraction of GROUND returns that are single-return
dens=tot/(RES*RES); meang=np.where(tot>0,angsum/np.maximum(tot,1),np.nan)

m=(~fld.ravel())&np.isfinite(pen)&(gnd>=5)
pe=pen[m]; gi=gint[m]; mf=multifrac[m]; gsf=gsingfrac[m]; dn=dens[m]; ag=meang[m]
print(f"cells: {m.sum()}")

# histogram of penetration
h,edges=np.histogram(pe,bins=50,range=(0,1))
print("\npenetration histogram (is it bimodal at ~0.6?):")
for i in range(0,50,2):
    bar="#"*int(60*h[i]/h.max())
    print(f"  {edges[i]:.2f}-{edges[i+2]:.2f}: {bar} {h[i]+h[i+1]}")

print("\nacross the 0.6 boundary (pen<0.6 vs >0.6):")
lo=pe<0.6; hi=pe>=0.6
for name,arr in [("multi-return frac (all)",mf),("ground-single frac",gsf),("density (pts/m2)",dn),("|scan angle|",ag),("ground intensity",gi)]:
    print(f"  {name:22s}: pen<0.6 {np.median(arr[lo]):.2f}   pen>=0.6 {np.median(arr[hi]):.2f}")

# is pen essentially 1 - (multi structure)? check pen vs multifrac
print(f"\n  corr(pen, multi-return frac) = {np.corrcoef(pe,mf)[0,1]:+.3f}  (strong negative => fraction is set by return structure)")

print(f"  corr(pen, |scan angle|) = {np.corrcoef(pe,ag)[0,1]:+.3f} ; corr(pen, density) = {np.corrcoef(pe,dn)[0,1]:+.3f}")
fig,ax=plt.subplots(1,3,figsize=(19,5.5))
ax[0].hist(pe,bins=60,range=(0,1),color="steelblue"); ax[0].axvline(0.6,color="r",ls="--")
ax[0].set_xlabel("gen1 ground-return fraction"); ax[0].set_ylabel("cells"); ax[0].set_title("penetration histogram (bimodal)")
hb=ax[1].hexbin(pe,ag,gridsize=45,bins="log",cmap="magma",mincnt=1); ax[1].axvline(0.6,color="c",ls="--")
ax[1].set_xlabel("ground-return fraction"); ax[1].set_ylabel("|scan angle| (deg)")
ax[1].set_title(f"fraction vs SCAN ANGLE (r={np.corrcoef(pe,ag)[0,1]:+.2f})"); fig.colorbar(hb,ax=ax[1])
sc=ax[2].hexbin(pe,gi,C=ag,gridsize=45,cmap="coolwarm",mincnt=1,vmin=0,vmax=14); ax[2].axvline(0.6,color="k",ls="--")
ax[2].set_xlabel("ground-return fraction"); ax[2].set_ylabel("ground intensity")
ax[2].set_title("pen-vs-intensity colored by |scan angle|"); fig.colorbar(sc,ax=ax[2],label="|scan angle| (deg)")
fig.suptitle("The ~0.6 split in gen1 ground-return fraction is SCAN GEOMETRY (swath edge/overlap vs nadir)",y=1.02)
fig.savefig("figures/refdatum/diagnose_06_split.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/diagnose_06_split.png")
