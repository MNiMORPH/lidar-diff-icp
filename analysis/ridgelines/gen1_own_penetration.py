#!/usr/bin/env python3
"""Redo the gen1 ground-sink test against gen1's OWN (leaf-off) penetration, not gen2's
leaf-on penetration.  gen1 was leaf-off (Nov 2008): its pulses passed through BARE branches.
Binning gen1 on gen2's leaf-on penetration conflates two different canopy states.

Compute gen1_pen = (gen1 class-2 returns)/(gen1 non-noise returns) per cell, then bin gen1
ground d by gen1_pen AND (for contrast) gen2_pen. If the sink follows gen1_pen -> it's the
bare-branch canopy gen1 actually traversed; if only gen2_pen -> it's tied to leaf-on
structure gen1 never saw (i.e. not a gen1-traversal effect).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_own_penetration.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zff=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
gen2pen=np.load("data/derived/elba_fulldensity/penetration.npy")  # gen2 leaf-on
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)

# --- pass: gen1 per-cell totals (for gen1 penetration) + per-return ground d ---
tot=np.zeros(NY*NX); gnd=np.zeros(NY*NX)
D=[];CELL=[]
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
        cl=np.asarray(pts.classification)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
        ix=ix[keep];iy=iy[keep];x=x[keep];y=y[keep];z=z[keep];cl=cl[keep]; cell=iy*NX+ix
        np.add.at(tot,cell,1); g=cl==2; np.add.at(gnd,cell,g.astype(float))
        # ground d
        cg=cell[g]; xg=x[g]; yg=y[g]; zg=z[g]
        xc=X0+((cg%NX)+0.5)*RES; yc=Y0+((cg//NX)+0.5)*RES
        d=(zg-(Zff[cg]+gxf[cg]*(xg-xc)+gyf[cg]*(yg-yc)))*cosf[cg]
        D.append(d); CELL.append(cg)
gen1pen=np.where(tot>0,gnd/np.maximum(tot,1),np.nan).reshape(NY,NX)
d=np.concatenate(D)*1000; cell=np.concatenate(CELL)
g1p=gen1pen.ravel()[cell]; g2p=gen2pen.ravel()[cell]
forest=((gen2pen<0.25)&~fld&np.isfinite(gen2pen)).ravel()[cell]   # keep same forest def for comparability

print(f"gen1 forest ground returns: {forest.sum():,}")
print(f"gen1 penetration (leaf-off) over forest: median {np.nanmedian(g1p[forest]):.2f}  "
      f"[gen2 leaf-on median {np.nanmedian(g2p[forest]):.2f}]")
print(f"corr(gen1_pen, gen2_pen) over forest cells = {np.corrcoef(g1p[forest],g2p[forest])[0,1]:+.3f}")

def bin_med(axis,dd,lbl):
    e=np.quantile(axis,np.linspace(0,1,7)); print(f"\n gen1 ground d binned by {lbl}:")
    mx=[];my=[]
    for i in range(len(e)-1):
        b=(axis>=e[i])&(axis<e[i+1] if i<len(e)-2 else axis<=e[i+1])
        if b.sum()<500: continue
        mx.append(np.median(axis[b])); my.append(np.median(dd[b]))
        print(f"   {lbl} {e[i]:.3f}-{e[i+1]:.3f} (med {np.median(axis[b]):.3f}): d {np.median(dd[b]):+6.1f} mm  n={b.sum()}")
    return np.array(mx),np.array(my)
fm=forest&np.isfinite(g1p)&np.isfinite(g2p)
x1,y1=bin_med(g1p[fm],d[fm],"gen1_pen (leaf-off)")
x2,y2=bin_med(g2p[fm],d[fm],"gen2_pen (leaf-on)")
print(f"\n sink range: vs gen1_pen {y1.max()-y1.min():.0f} mm ; vs gen2_pen {y2.max()-y2.min():.0f} mm")
print(f" corr(d,gen1_pen)={np.corrcoef(d[fm],g1p[fm])[0,1]:+.3f}  corr(d,gen2_pen)={np.corrcoef(d[fm],g2p[fm])[0,1]:+.3f}")

fig,ax=plt.subplots(figsize=(9,6))
ax.plot(x1,y1,"C0o-",label="vs gen1 own penetration (leaf-off, bare branches)")
ax.plot(x2,y2,"C3s--",label="vs gen2 penetration (leaf-on) — the misleading axis")
ax.set_xlabel("penetration = ground-return fraction  (LEFT = denser)"); ax.invert_xaxis()
ax.set_ylabel("gen1 ground median d (mm)")
ax.set_title("gen1 ground sink: against its OWN leaf-off canopy vs against gen2's leaf-on canopy")
ax.legend(fontsize=9); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/gen1_own_penetration.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_own_penetration.png")
