#!/usr/bin/env python3
"""Relative ground elevation (per-cell gen1 ground median d) vs three gen1-internal covariates:
scan angle, ground-return fraction (penetration), and ground-return intensity. 3 panels.
Per-cell d from slope_normal_returns.npz (gen1 class-2 median); angle/pen/intensity streamed
from the gen1 cloud. Forest cells.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_elev_vs_3.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy")
dmed=np.load("data/derived/elba_fulldensity/slope_normal_returns.npz")["gen1_ground_median_d"]*1000  # mm
N=NY*NX; tot=np.zeros(N); gnd=np.zeros(N); isum=np.zeros(N); angsum=np.zeros(N)
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); cl=np.asarray(pts.classification)
        inten=np.asarray(pts.intensity).astype(float); ang=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7); cell=iy[keep]*NX+ix[keep]
        clk=cl[keep]; ik=inten[keep]; ak=ang[keep]; g=clk==2
        np.add.at(tot,cell,1); np.add.at(gnd,cell,g.astype(float)); np.add.at(isum,cell,np.where(g,ik,0.0)); np.add.at(angsum,cell,ak)
pen=np.where(tot>0,gnd/np.maximum(tot,1),np.nan)
gint=np.where(gnd>0,isum/np.maximum(gnd,1),np.nan)
ang=np.where(tot>0,angsum/np.maximum(tot,1),np.nan)
d=dmed.ravel()

msk=(g2pen.ravel()<0.25)&(~fld.ravel())&np.isfinite(pen)&np.isfinite(d)&(gnd>=8)
D=d[msk]; A=ang[msk]; P=pen[msk]; I=gint[msk]
print(f"forest cells: {msk.sum()}")
print(f"corr(d, angle)={np.corrcoef(D,A)[0,1]:+.3f}  corr(d, pen)={np.corrcoef(D,P)[0,1]:+.3f}  corr(d, intensity)={np.corrcoef(D,I)[0,1]:+.3f}")

def med_line(x,y,nb=11):
    e=np.quantile(x,np.linspace(0,1,nb)); mx=[];my=[]
    for i in range(len(e)-1):
        b=(x>=e[i])&(x<e[i+1] if i<len(e)-2 else x<=e[i+1])
        if b.sum()<50: continue
        mx.append(np.median(x[b])); my.append(np.median(y[b]))
    return mx,my

fig,ax=plt.subplots(1,3,figsize=(19,6))
for a,(xv,xl,r) in zip(ax,[(A,"|scan angle| (deg)",np.corrcoef(D,A)[0,1]),
                           (P,"ground-return fraction (penetration)",np.corrcoef(D,P)[0,1]),
                           (I,"ground-return intensity",np.corrcoef(D,I)[0,1])]):
    hb=a.hexbin(xv,D,gridsize=45,bins="log",cmap="viridis",mincnt=1)
    mx,my=med_line(xv,D); a.plot(mx,my,"r.-",lw=2,label="median")
    a.set_xlabel(xl); a.set_ylabel("relative ground elevation d (mm)"); a.set_ylim(np.percentile(D,1),np.percentile(D,99))
    a.set_title(f"d vs {xl.split('(')[0].strip()}  (r={r:+.2f})"); a.legend(fontsize=8); fig.colorbar(hb,ax=a)
fig.suptitle("gen1 relative ground elevation vs scan angle, ground-fraction, intensity (forest, per cell)",y=1.0)
fig.savefig("figures/refdatum/gen1_elev_vs_3.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("wrote figures/refdatum/gen1_elev_vs_3.png")
