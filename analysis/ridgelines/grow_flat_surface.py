#!/usr/bin/env python3
"""Seed-grown flat reference surface for vertical tie points.
From a seed cell, flood-fill outward over cells that are (a) low-slope and (b) co-planar with
a plane continuously refit to the growing region (|z - plane| < tol). This follows the actual
feature (parking lot, roof, road) and stops at curbs / edges / vegetation, adapting the
footprint instead of a fixed window. Then computes the gen1-gen2 vertical tie (median gen1
d_mm) over the grown surface.

    ./lidar-icp/bin/python analysis/ridgelines/grow_flat_surface.py
"""
import numpy as np
from collections import deque

NX,NY=508,700; X0,Y0=577492.8,4882737.6; RES=5.0
Z=np.load("data/derived/elba_fulldensity/z_after.npy")      # gen2 bare earth (m)
SL=np.load("data/derived/elba_fulldensity/slope.npy")       # deg
A=np.load("data/derived/elba_fulldensity/gen1_csf_angles.npz")
cell=A["cell"]; dmm=A["d_mm"]; psid=A["point_source_id"]
# per-cell gen1 offset (median d_mm) and count
off=np.full(NX*NY,np.nan); cnt=np.zeros(NX*NY,int)
o=np.argsort(cell); cs=cell[o]; ds=dmm[o]; uq,st=np.unique(cs,return_index=True)
for k,c in enumerate(uq):
    s=st[k]; e=st[k+1] if k+1<len(st) else len(cs); off[c]=np.median(ds[s:e]); cnt[c]=e-s
off=off.reshape(NY,NX); cnt=cnt.reshape(NY,NX)

def grow(six,siy,max_slope=2.5,plane_tol=0.04,max_cells=800,refit=15):
    def ok(x,y): return 0<=x<NX and 0<=y<NY and np.isfinite(Z[y,x])
    if not ok(six,siy): return None
    reg={(six+dx,siy+dy) for dx in(-1,0,1) for dy in(-1,0,1) if ok(six+dx,siy+dy) and SL[siy+dy,six+dx]<max_slope}
    reg.add((six,siy))
    def fit(R):
        P=np.array(list(R)); A_=np.c_[np.ones(len(P)),P[:,0],P[:,1]]; zz=Z[P[:,1],P[:,0]]
        c,*_=np.linalg.lstsq(A_,zz,rcond=None); return c
    coef=fit(reg); q=deque(reg); k=0
    while q and len(reg)<max_cells:
        cx,cy=q.popleft()
        for dx,dy in((1,0),(-1,0),(0,1),(0,-1)):
            nx,ny=cx+dx,cy+dy
            if (nx,ny) in reg or not ok(nx,ny) or SL[ny,nx]>=max_slope: continue
            if abs(Z[ny,nx]-(coef[0]+coef[1]*nx+coef[2]*ny))<plane_tol:
                reg.add((nx,ny)); q.append((nx,ny)); k+=1
                if k%refit==0: coef=fit(reg)
    return reg

def tie(name,six,siy):
    reg=grow(six,siy)
    if reg is None: print(f"  {name}: seed off-grid"); return
    P=np.array(list(reg)); ix=P[:,0]; iy=P[:,1]
    o=off[iy,ix]; o=o[np.isfinite(o)]
    cells=set((iy*NX+ix).tolist()); m=np.isin(cell,list(cells))
    ext_e=(ix.max()-ix.min()+1)*RES; ext_n=(iy.max()-iy.min()+1)*RES
    z=Z[iy,ix]; A_=np.c_[np.ones(len(z)),ix,iy]; c,*_=np.linalg.lstsq(A_,z,rcond=None); rough=np.std(z-A_@c)
    print(f"  {name}: grown {len(reg)} cells (~{ext_e:.0f}x{ext_n:.0f} m), planar roughness {rough*1000:.0f} mm")
    if o.size>=3:
        print(f"       tie offset median {np.median(o):+.1f} mm  NMAD {1.4826*np.median(np.abs(o-np.median(o))):.1f}  "
              f"(n_cells {o.size}, gen1 returns {m.sum()}, swaths {np.unique(psid[m]).tolist()})  -> lift +{-np.median(o):.1f} mm")

# in-grid seeds
print("Seed-grown flat-surface ties (max_slope 2.5deg, plane_tol 40mm):")
tie("point 1        (E579671 N4885761)", 435,604)
tie("small lot      (E579144 N4883394)", 330,131)
