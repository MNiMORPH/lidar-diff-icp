#!/usr/bin/env python3
"""Rebuild the fine (1 cm) slope-normal pooled histograms per land-cover stratum for:
  g1csf  = gen1 GROUND as classified by the pipeline's CSF cloth (data/csf_cache/elba.las)
  g1vend = gen1 ground by the 2008 VENDOR classification (class 2 in 4342-29-64.laz)  [reference]
  g1all  = gen1 ALL returns (4342-29-64.laz, every point)
  g2gnd  = gen2 ground by the INTERNAL 3DEP classification (class 2 in the fulldensity cloud)
  g2all  = gen2 ALL returns (fulldensity cloud, every point)
Same slope-normal transform as slope_normal_returns.py (plane = gen2 bare earth), copied verbatim.
Saves data/derived/elba_fulldensity/ground_fine_csf_all.npz.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/ground_fine_csf_all.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

NY, NX = 700, 508; X0, Y0 = 577492.8, 4882737.6; RES = 5.0
CSF   = "data/csf_cache/elba.las"                 # gen1 CSF ground (all class 2)
GEN1  = "data/before/4342-29-64.laz"              # gen1 raw (all returns + vendor class)
GEN2  = "data/after/3dep2021_fulldensity.laz"     # gen2 full density (all returns + internal class)
CHUNK = 5_000_000
Zg = np.load("data/derived/elba_fulldensity/z_after.npy")
Zg_filled = Zg.copy(); m = ~np.isfinite(Zg_filled)
if m.any(): Zg_filled = Zg_filled[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zg_filled, RES); cos_slope = 1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zg_flat=Zg_filled.ravel(); gx_flat=gx.ravel(); gy_flat=gy.ravel(); cos_flat=cos_slope.ravel()

pen = np.load("data/derived/elba_fulldensity/penetration.npy")
fld = np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
strat = np.zeros(NY*NX, np.int8)
strat[((pen<0.25)&~fld&np.isfinite(pen)).ravel()] = 1     # forest
strat[((pen>=0.45)&~fld&np.isfinite(pen)).ravel()] = 2    # open

FLO, FHI, FW = -0.8, 4.0, 0.01
fedges = np.arange(FLO, FHI+0.5*FW, FW); fc = 0.5*(fedges[:-1]+fedges[1:]); NF = fc.size
names = ["g1csf","g1vend","g1all","g2gnd","g2all"]
H = {(n,s): np.zeros(NF) for n in names for s in (1,2)}

def accumulate(path, targets):
    """targets: list of (name, filt) where filt(cl)->bool mask (None = all)."""
    tot=0
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
            cl=np.asarray(pts.classification)
            ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
            keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
            ix=ix[keep];iy=iy[keep];x=x[keep];y=y[keep];z=z[keep];cl=cl[keep]; tot+=keep.sum()
            cell=iy*NX+ix; xc=X0+(ix+0.5)*RES; yc=Y0+(iy+0.5)*RES
            d=(z-(Zg_flat[cell]+gx_flat[cell]*(x-xc)+gy_flat[cell]*(y-yc)))*cos_flat[cell]
            st=strat[cell]; fb=np.searchsorted(fedges,d,side="right")-1; ok=(fb>=0)&(fb<NF)
            for name,filt in targets:
                sel = ok if filt is None else (ok & filt(cl))
                for s in (1,2):
                    ss = sel & (st==s)
                    if ss.any(): np.add.at(H[(name,s)], fb[ss], 1)
    print(f"  {path}: kept {tot:,}")

print("streaming CSF cache (gen1 cloth ground)...");  accumulate(CSF,  [("g1csf", None)])
print("streaming gen1 raw (all + vendor ground)..."); accumulate(GEN1, [("g1all", None), ("g1vend", lambda cl: cl==2)])
print("streaming gen2 full density (all + ground)..."); accumulate(GEN2, [("g2all", None), ("g2gnd", lambda cl: cl==2)])

for n in names:
    print(f"  {n}: forest {H[(n,1)].sum():,.0f}  open {H[(n,2)].sum():,.0f}")
np.savez_compressed("data/derived/elba_fulldensity/ground_fine_csf_all.npz",
                    fedges=fedges, fc=fc, **{f"{n}_{s}": H[(n,s)] for n in names for s in (1,2)})
print("saved data/derived/elba_fulldensity/ground_fine_csf_all.npz")
