#!/usr/bin/env python3
"""NEAR-NADIR CSF cloth ground returns only (|scan angle|<=2 deg). Reconstruct forest-floor
elevation histograms for ALL forest vs CORE forest. Tests whether the all->core sink survives
when the oblique under-penetration is removed (near-nadir only).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_csf_nadir_all_vs_core.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
CSF="data/csf_cache/elba.las"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zff=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy").ravel()
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool).ravel()
allf=(g2pen<0.25)&(~fld)&np.isfinite(g2pen); core=np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()

Da=[];Dc=[]; Aa=[];Ac=[]   # near-nadir (D*) and all-angle (A*) for all-forest / core
with laspy.open(CSF) as f:
    for pts in f.chunk_iterator(CHUNK):
        aa=np.abs(np.asarray(pts.scan_angle).astype(float)*0.006)
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        ing=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); cell=np.where(ing,iy*NX+ix,0)
        nad=ing&(aa<=2.0)
        def dof(sel):
            c=cell[sel]; xc=X0+((c%NX)+0.5)*RES; yc=Y0+((c//NX)+0.5)*RES
            return (z[sel]-(Zff[c]+gxf[c]*(x[sel]-xc)+gyf[c]*(y[sel]-yc)))*cosf[c]*1000
        for sel,DD in [(nad&allf[cell],Da),(nad&core[cell],Dc),(ing&allf[cell],Aa),(ing&core[cell],Ac)]:
            if sel.any(): DD.append(dof(sel))
da=np.concatenate(Da); dc=np.concatenate(Dc); aa_=np.concatenate(Aa); ac=np.concatenate(Ac)
print("CSF cloth ground forest-floor elevation, same method:")
print(f"  ALL-ANGLE : all forest median {np.median(aa_):+.1f} (n={aa_.size:,})  core {np.median(ac):+.1f} (n={ac.size:,})  shift {np.median(ac)-np.median(aa_):+.1f} mm")
print(f"  NEAR-NADIR: all forest median {np.median(da):+.1f} (n={da.size:,})  core {np.median(dc):+.1f} (n={dc.size:,})  shift {np.median(dc)-np.median(da):+.1f} mm")
print(f"  => near-nadir shrinks the all->core sink from {np.median(ac)-np.median(aa_):+.1f} to {np.median(dc)-np.median(da):+.1f} mm")

fig,ax=plt.subplots(figsize=(10,6))
b=np.arange(-350,251,10)
for d,col,l in [(da,"C2","ALL forest"),(dc,"C1","CORE forest")]:
    ax.hist(d,bins=b,density=True,histtype="step",lw=2,color=col,
            label=f"{l}  (n={d.size:,}, med {np.median(d):+.0f} mm)")
    ax.axvline(np.median(d),color=col,ls=":",lw=1.5)
ax.set_xlabel("forest-floor elevation d (mm)  [slope-normal, vs gen2 bare earth]"); ax.set_ylabel("density")
ax.set_title("NEAR-NADIR CSF cloth ground: forest-floor elevation, all forest vs core forest")
ax.legend(); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/gen1_csf_nadir_all_vs_core.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_csf_nadir_all_vs_core.png")
