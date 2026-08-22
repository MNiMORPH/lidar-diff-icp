#!/usr/bin/env python3
"""OBLIQUE CSF cloth ground returns only (|scan angle|>=8 deg). Forest-floor elevation
histograms for FARMLAND vs FOREST (all) vs CORE FOREST — the oblique counterpart to the
near-nadir all/core forest plots.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_csf_oblique_strata_hist.py
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
strata={"farmland":(g2pen>=0.45)&(~fld)&np.isfinite(g2pen),
        "forest (all)":(g2pen<0.25)&(~fld)&np.isfinite(g2pen),
        "core forest":np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()}
D={s:[] for s in strata}
with laspy.open(CSF) as f:
    for pts in f.chunk_iterator(CHUNK):
        aa=np.abs(np.asarray(pts.scan_angle).astype(float)*0.006); obl=aa>=8.0
        if not obl.any(): continue
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        ok=obl&(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); cell=np.where(ok,iy*NX+ix,0)
        for s,mask in strata.items():
            sel=ok&mask[cell]
            if not sel.any(): continue
            c=cell[sel]; xc=X0+((c%NX)+0.5)*RES; yc=Y0+((c//NX)+0.5)*RES
            d=(z[sel]-(Zff[c]+gxf[c]*(x[sel]-xc)+gyf[c]*(y[sel]-yc)))*cosf[c]
            D[s].append(d*1000)
D={s:np.concatenate(v) for s,v in D.items()}
print("OBLIQUE CSF cloth ground, forest-floor elevation:")
for s in strata: print(f"  {s:14s}: n={D[s].size:,}  median {np.median(D[s]):+.1f}  mean {np.mean(D[s]):+.1f} mm")

fig,ax=plt.subplots(figsize=(10,6)); b=np.arange(-350,251,10)
cols={"farmland":"C2","forest (all)":"C0","core forest":"C1"}
for s in strata:
    ax.hist(D[s],bins=b,density=True,histtype="step",lw=2,color=cols[s],
            label=f"{s}  (n={D[s].size:,}, med {np.median(D[s]):+.0f} mm)")
    ax.axvline(np.median(D[s]),color=cols[s],ls=":",lw=1.5)
ax.set_xlabel("forest-floor elevation d (mm)  [slope-normal, vs gen2 bare earth]"); ax.set_ylabel("density")
ax.set_title("OBLIQUE (|angle|>=8°) CSF cloth ground: floor elevation — farmland vs forest vs core forest")
ax.legend(); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/gen1_csf_oblique_strata_hist.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_csf_oblique_strata_hist.png")
