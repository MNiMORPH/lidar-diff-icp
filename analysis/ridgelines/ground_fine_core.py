#!/usr/bin/env python3
"""Re-accumulate the fine (1 cm) slope-normal ground-return histograms stratified by the
ROBUST CORE strata (core_forest.npy / core_open.npy) instead of the raw penetration masks,
to see how much the cleaner (uncontaminated) classification sharpens the distribution forms.
gen1 ground = CSF cloth (elba.las); gen2 ground = internal class 2. Also keeps all-returns.
Same slope-normal transform (plane = gen2 bare earth). Saves ground_fine_core.npz.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/ground_fine_core.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
CSF="data/csf_cache/elba.las"; GEN1="data/before/4342-29-64.laz"; GEN2="data/after/3dep2021_fulldensity.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zf=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()

# CORE strata as the per-cell label (1 forest, 2 open, 0 neither)
cf=np.load("data/derived/elba_fulldensity/core_forest.npy")
co=np.load("data/derived/elba_fulldensity/core_open.npy")
strat=np.zeros(NY*NX,np.int8); strat[cf.ravel()]=1; strat[co.ravel()]=2

FLO,FHI,FW=-0.8,4.0,0.01; fedges=np.arange(FLO,FHI+0.5*FW,FW); fc=0.5*(fedges[:-1]+fedges[1:]); NF=fc.size
names=["g1csf","g1all","g2gnd","g2all"]; H={(n,s):np.zeros(NF) for n in names for s in(1,2)}
def acc(path,targets):
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            x=np.asarray(pts.x,np.float64);y=np.asarray(pts.y,np.float64);z=np.asarray(pts.z,np.float64)
            cl=np.asarray(pts.classification)
            ix=((x-X0)/RES).astype(np.int64);iy=((y-Y0)/RES).astype(np.int64)
            keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7)
            ix=ix[keep];iy=iy[keep];x=x[keep];y=y[keep];z=z[keep];cl=cl[keep]
            cell=iy*NX+ix; xc=X0+(ix+0.5)*RES; yc=Y0+(iy+0.5)*RES
            d=(z-(Zf[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*cosf[cell]
            st=strat[cell]; fb=np.searchsorted(fedges,d,side="right")-1; ok=(fb>=0)&(fb<NF)
            for name,filt in targets:
                sel=ok if filt is None else (ok&filt(cl))
                for s in(1,2):
                    ss=sel&(st==s)
                    if ss.any(): np.add.at(H[(name,s)],fb[ss],1)
print("CSF (gen1 cloth ground)...");  acc(CSF,  [("g1csf",None)])
print("gen1 raw (all)...");           acc(GEN1, [("g1all",None)])
print("gen2 full density (all+gnd)..."); acc(GEN2,[("g2all",None),("g2gnd",lambda cl:cl==2)])
for n in names: print(f"  {n}: forest {H[(n,1)].sum():,.0f}  open {H[(n,2)].sum():,.0f}")
np.savez_compressed("data/derived/elba_fulldensity/ground_fine_core.npz",fedges=fedges,fc=fc,
                    **{f"{n}_{s}":H[(n,s)] for n in names for s in(1,2)})
print("saved data/derived/elba_fulldensity/ground_fine_core.npz")
