#!/usr/bin/env python3
"""CSF CLOTH ground returns (data/csf_cache/elba.las), pooled into NEAR-NADIR (|angle|<=2 deg)
and OBLIQUE (|angle|>=8 deg). In CORE forest, show histograms of forest-floor elevation
(slope-normal d vs the gen2 bare-earth plane) for the two return populations.

CSF cache is LAS point-format 7 -> uses 'scan_angle' (scaled 0.006 deg/unit), not scan_angle_rank.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_csf_nadir_oblique_hist.py
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
core=np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()

Dn=[];Do=[]; amax=0
with laspy.open(CSF) as f:
    for pts in f.chunk_iterator(CHUNK):
        sa=np.asarray(pts.scan_angle).astype(float)*0.006   # -> degrees
        aa=np.abs(sa); amax=max(amax,aa.max())
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        ok=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); cell=np.where(ok,iy*NX+ix,0); cm=ok&core[cell]
        for sel,DD in [(cm&(aa<=2.0),Dn),(cm&(aa>=8.0),Do)]:
            if not sel.any(): continue
            c=cell[sel]; xc=X0+((c%NX)+0.5)*RES; yc=Y0+((c//NX)+0.5)*RES
            d=(z[sel]-(Zff[c]+gxf[c]*(x[sel]-xc)+gyf[c]*(y[sel]-yc)))*cosf[c]
            DD.append(d*1000)
dn=np.concatenate(Dn); do=np.concatenate(Do)
print(f"CSF |scan_angle| max = {amax:.1f} deg")
def stat(d,l): print(f"  {l}: n={d.size:,}  mode~{np.round(np.median(d),0):+.0f} median {np.median(d):+.1f}  mean {np.mean(d):+.1f} mm")
print("core-forest CSF forest-floor elevation:"); stat(dn,"near-nadir (<=2)"); stat(do,"oblique (>=8)")
print(f"  oblique - nadir (median): {np.median(do)-np.median(dn):+.1f} mm")

fig,ax=plt.subplots(figsize=(10,6))
b=np.arange(-350,251,10)
for d,col,l in [(dn,"C0","near-nadir (|angle|<=2°)"),(do,"C3","oblique (|angle|>=8°)")]:
    ax.hist(d,bins=b,density=True,histtype="step",lw=2,color=col,
            label=f"{l}  (n={d.size:,}, med {np.median(d):+.0f} mm)")
    ax.axvline(np.median(d),color=col,ls=":",lw=1.5)
ax.set_xlabel("forest-floor elevation d (mm)  [slope-normal, vs gen2 bare earth]"); ax.set_ylabel("density")
ax.set_title("CSF cloth ground in CORE forest: forest-floor elevation, near-nadir vs oblique returns")
ax.legend(); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/gen1_csf_nadir_oblique_hist.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_csf_nadir_oblique_hist.png")
