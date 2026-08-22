#!/usr/bin/env python3
"""Near-nadir only: relative ground elevation vs return intensity (gen1, forest).
Restricting to |scan angle| <= 2 deg removes the scan-angle/range/obliquity confound, so any
intensity->elevation (range-walk) relationship should show cleaner. Per-return.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_nadir_elev_intensity.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0; NADIR=2.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zff=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
is_forest=np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()   # CORE forest only

D=[];I=[]
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        cl=np.asarray(pts.classification); ang=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        sel=(cl==2)&(ang<=NADIR)
        if not sel.any(): continue
        x=np.asarray(pts.x,np.float64)[sel]; y=np.asarray(pts.y,np.float64)[sel]; z=np.asarray(pts.z,np.float64)[sel]
        inten=np.asarray(pts.intensity).astype(float)[sel]
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); cell=iy[keep]*NX+ix[keep]
        fm=is_forest[cell]
        if not fm.any(): continue
        cell=cell[fm]; x=x[keep][fm]; y=y[keep][fm]; z=z[keep][fm]
        xc=X0+((cell%NX)+0.5)*RES; yc=Y0+((cell//NX)+0.5)*RES
        d=(z-(Zff[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*cosf[cell]
        D.append(d*1000); I.append(inten[keep][fm])
d=np.concatenate(D); I=np.concatenate(I)
print(f"near-nadir (|angle|<={NADIR}) CORE forest ground returns: {d.size:,}")
print(f"corr(d, intensity) = {np.corrcoef(d,I)[0,1]:+.3f}")
print("\nmedian d by intensity bin (near-nadir):")
e=np.quantile(I,np.linspace(0,1,11))
mx=[];my=[]
for i in range(len(e)-1):
    b=(I>=e[i])&(I<e[i+1] if i<len(e)-2 else I<=e[i+1])
    if b.sum()<300: continue
    mx.append(np.median(I[b])); my.append(np.median(d[b]))
    print(f"  intensity {e[i]:5.0f}-{e[i+1]:5.0f} (med {np.median(I[b]):4.0f}): d {np.median(d[b]):+7.1f} mm  n={b.sum()}")

fig,ax=plt.subplots(figsize=(9,6))
hb=ax.hexbin(I,d,gridsize=50,bins="log",cmap="viridis",mincnt=1,extent=(0,60,-400,300))
ax.plot(mx,my,"r.-",lw=2.5,ms=10,label="median d")
ax.set_xlim(0,60); ax.set_ylim(-350,250)
ax.set_xlabel("ground-return intensity"); ax.set_ylabel("relative ground elevation d (mm)")
ax.set_title(f"gen1 near-nadir (|angle|<={NADIR:.0f}deg) CORE forest: elevation vs intensity  (r={np.corrcoef(d,I)[0,1]:+.2f}, n={d.size:,})")
ax.legend(); fig.colorbar(hb,label="log10 count")
fig.savefig("figures/refdatum/gen1_nadir_elev_intensity_core.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_nadir_elev_intensity_core.png")
