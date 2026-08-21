#!/usr/bin/env python3
"""Per-cell gen1 ground-return FRACTION (penetration = gen1 class-2 / gen1 non-noise) vs gen1
ground-return INTENSITY (mean intensity of gen1 class-2 returns in the cell). Both purely
gen1-internal. Tests whether gen1's own canopy density and its ground-return strength are
linked (attenuation: denser canopy -> weaker ground returns -> lower fraction).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_pen_vs_intensity.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy")   # only to pick forest cells

tot=np.zeros(NY*NX); gnd=np.zeros(NY*NX); isum=np.zeros(NY*NX)
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); cl=np.asarray(pts.classification)
        inten=np.asarray(pts.intensity).astype(float)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
        cell=(iy[keep]*NX+ix[keep]); clk=cl[keep]; ik=inten[keep]
        np.add.at(tot,cell,1); g=clk==2
        np.add.at(gnd,cell,g.astype(float)); np.add.at(isum,cell,np.where(g,ik,0.0))
pen=np.where(tot>0,gnd/np.maximum(tot,1),np.nan)
gint=np.where(gnd>0,isum/np.maximum(gnd,1),np.nan)
pen=pen.reshape(NY,NX); gint=gint.reshape(NY,NX)

forest=(g2pen<0.25)&~fld&np.isfinite(g2pen)&np.isfinite(pen)&np.isfinite(gint)&(gnd.reshape(NY,NX)>=5)
allc=~fld&np.isfinite(pen)&np.isfinite(gint)&(gnd.reshape(NY,NX)>=5)
for msk,lbl in [(allc,"all cells"),(forest,"forest cells")]:
    p=pen[msk]; g=gint[msk]
    print(f"{lbl}: n={msk.sum()}  corr(gen1_pen, gen1_ground_intensity) = {np.corrcoef(p,g)[0,1]:+.3f}")
    e=np.quantile(p,np.linspace(0,1,9))
    for i in range(len(e)-1):
        b=(p>=e[i])&(p<e[i+1] if i<len(e)-2 else p<=e[i+1])
        if b.sum()<50: continue
        print(f"    pen {e[i]:.3f}-{e[i+1]:.3f} (med {np.median(p[b]):.3f}): ground intensity med {np.median(g[b]):.1f}")

fig,ax=plt.subplots(1,2,figsize=(14,6))
for a,(msk,lbl) in zip(ax,[(allc,"all cells"),(forest,"forest cells")]):
    p=pen[msk]; g=gint[msk]
    hb=a.hexbin(p,g,gridsize=45,bins="log",cmap="viridis",mincnt=1)
    e=np.quantile(p,np.linspace(0,1,11)); mx=[];my=[]
    for i in range(len(e)-1):
        b=(p>=e[i])&(p<e[i+1] if i<len(e)-2 else p<=e[i+1])
        if b.sum()<30: continue
        mx.append(np.median(p[b])); my.append(np.median(g[b]))
    a.plot(mx,my,"r.-",lw=2,label="median")
    a.set_xlabel("gen1 ground-return fraction (penetration)"); a.set_ylabel("gen1 ground-return intensity (mean)")
    a.set_title(f"{lbl}  (r={np.corrcoef(p,g)[0,1]:+.2f})"); a.legend(fontsize=9)
    fig.colorbar(hb,ax=a,label="log10 cell count")
fig.suptitle("gen1 (leaf-off): ground-return fraction vs ground-return intensity — both internal to gen1",y=1.0)
fig.savefig("figures/refdatum/gen1_pen_vs_intensity.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_pen_vs_intensity.png")
