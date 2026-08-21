#!/usr/bin/env python3
"""Plot gen1 ground-return distributions for the two classifications — CSF cloth vs internal
(2008 vendor class-2) — in ALL forest vs CORE forest, to see the all->core shift per classifier.
gen1 clouds only (CSF from elba.las; vendor class-2 from the raw 2008 tile). Same slope-normal
transform (plane = gen2 bare earth); 1 cm bins.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_csf_vs_internal.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
CSF="data/csf_cache/elba.las"; GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zf=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()

pen=np.load("data/derived/elba_fulldensity/penetration.npy"); fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
cf=np.load("data/derived/elba_fulldensity/core_forest.npy")
is_forest=((pen<0.25)&~fld&np.isfinite(pen)).ravel()      # all forest
is_core=cf.ravel()                                        # core forest (subset)

FLO,FHI,FW=-0.8,1.6,0.01; fe=np.arange(FLO,FHI+0.5*FW,FW); fc=0.5*(fe[:-1]+fe[1:]); NF=fc.size
H={k:np.zeros(NF) for k in ["csf_all","csf_core","vend_all","vend_core"]}
def acc(path,gnd_filter,pre):
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            x=np.asarray(pts.x,np.float64);y=np.asarray(pts.y,np.float64);z=np.asarray(pts.z,np.float64)
            cl=np.asarray(pts.classification)
            g = np.ones(len(cl),bool) if gnd_filter is None else (cl==gnd_filter)
            if not g.any(): continue
            x=x[g];y=y[g];z=z[g]
            ix=((x-X0)/RES).astype(np.int64);iy=((y-Y0)/RES).astype(np.int64)
            keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); ix=ix[keep];iy=iy[keep];x=x[keep];y=y[keep];z=z[keep]
            cell=iy*NX+ix; xc=X0+(ix+0.5)*RES; yc=Y0+(iy+0.5)*RES
            d=(z-(Zf[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*cosf[cell]
            fb=np.searchsorted(fe,d,side="right")-1; ok=(fb>=0)&(fb<NF)
            fa=ok&is_forest[cell]; co=ok&is_core[cell]
            if fa.any(): np.add.at(H[pre+"_all"],fb[fa],1)
            if co.any(): np.add.at(H[pre+"_core"],fb[co],1)
print("CSF cloth..."); acc(CSF,None,"csf")
print("vendor class-2..."); acc(GEN1,2,"vend")

def med(c): c=c.astype(float); cdf=np.cumsum(c)/c.sum(); return np.interp(.5,cdf,fc)*1000
print("\ngen1 ground median (mm), forest all -> core:")
print(f"  CSF cloth : {med(H['csf_all']):+.1f} -> {med(H['csf_core']):+.1f}   (shift {med(H['csf_core'])-med(H['csf_all']):+.1f})")
print(f"  internal  : {med(H['vend_all']):+.1f} -> {med(H['vend_core']):+.1f}   (shift {med(H['vend_core'])-med(H['vend_all']):+.1f})")

fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,(suf,lbl) in zip(axes,[("all","ALL forest"),("core","CORE forest")]):
    for pre,col,tag in [("csf","C0","CSF cloth ground"),("vend","C3","internal (2008 class-2)")]:
        c=H[f"{pre}_{suf}"]; d=c/c.sum()/FW
        ax.semilogx(np.where(d>0,d,np.nan),fc,col,lw=1.7,label=f"{tag}  (med {med(c):+.0f} mm)")
        ax.axhline(med(c)/1000,color=col,ls=":",lw=1)
    ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.5,0.4); ax.set_xlim(1e-2,3e1)
    ax.set_xlabel("density (1/m, log)"); ax.set_title(f"gen1 ground — {lbl}"); ax.legend(fontsize=9); ax.grid(alpha=.3,which="both")
axes[0].set_ylabel("slope-normal d (m)  [plane = gen2 bare earth]")
fig.suptitle("gen1 ground return: CSF cloth vs internal (vendor class-2), all vs core forest",y=1.0)
fig.savefig("figures/refdatum/gen1_csf_vs_internal.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_csf_vs_internal.png")
