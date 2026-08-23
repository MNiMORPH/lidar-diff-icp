#!/usr/bin/env python3
"""(a) Residual vertical offset at Andy's stable tie surfaces on the ELBAEXT grid.

For each tie point, SEARCH THE NEIGHBOURHOOD for the flattest local patch (the recorded
coordinate may sit on a curb/edge/slope), seed the flood-fill there, grow a co-planar
low-slope surface on the gen2 z_after, and take the median elbaext DoD (gen2-gen1,
reference_plane tie) over it = the RESIDUAL offset AFTER the elbaext datum tie.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/elbaext_tie_offsets.py
"""
import json, numpy as np
from collections import deque
D="data/derived/elbaext/"
meta=json.load(open(D+"meta.json"))
NX,NY=int(meta["nx"]),int(meta["ny"]); RES=float(meta["res"])
X0,Y0=float(meta["bounds"][0]),float(meta["bounds"][1])
Z=np.load(D+"z_after.npy"); SL=np.load(D+"slope.npy"); DOD=np.load(D+"dod.npy")*1000.0        # reference_plane (mm)
DODG=np.load(D+"dod_geoid.npy")*1000.0  # geoid (mm)
print(f"elbaext grid {NX}x{NY} res {RES}  origin ({X0:.0f},{Y0:.0f})")

tp=np.load("data/derived/elba_fulldensity/andy_tie_points.npz")
E,N=tp["easting"],tp["northing"]; cert=tp["certain"]; lat,lon=tp["lat"],tp["lon"]
ix0=np.floor((E-X0)/RES).astype(int); iy0=np.floor((N-Y0)/RES).astype(int)
ing=(ix0>=0)&(ix0<NX)&(iy0>=0)&(iy0<NY)
print(f"tie points: {len(E)} total, {ing.sum()} inside elbaext grid\n")

def ok(x,y): return 0<=x<NX and 0<=y<NY and np.isfinite(Z[y,x])

def find_flat_seed(ix,iy,max_move_m=10.0,max_slope=2.5,win=2):
    R=int(np.ceil(max_move_m/RES))
    """Search +-R cells for the flattest local patch: min planar-residual over a
    (2*win+1) window among low-slope finite cells; tie-break toward the point."""
    best=None; bestscore=1e9
    for dy in range(-R,R+1):
        for dx in range(-R,R+1):
            if RES*np.hypot(dx,dy)>max_move_m: continue
            x,y=ix+dx,iy+dy
            if not ok(x,y) or SL[y,x]>=max_slope: continue
            P=[(x+a,y+b) for a in range(-win,win+1) for b in range(-win,win+1) if ok(x+a,y+b)]
            if len(P)<(2*win+1)**2*0.7: continue
            P=np.array(P); A=np.c_[np.ones(len(P)),P[:,0],P[:,1]]
            c,*_=np.linalg.lstsq(A,Z[P[:,1],P[:,0]],rcond=None); resid=np.std(Z[P[:,1],P[:,0]]-A@c)*1000
            score=resid+0.15*RES*np.hypot(dx,dy)   # mm + gentle distance penalty (mm/ m)
            if score<bestscore: bestscore=score; best=(x,y,resid)
    return best

def grow(sx,sy,max_slope=2.5,plane_tol=0.04,max_cells=1200,refit=15):
    if not ok(sx,sy): return None
    reg={(sx+dx,sy+dy) for dx in(-1,0,1) for dy in(-1,0,1) if ok(sx+dx,sy+dy) and SL[sy+dy,sx+dx]<max_slope}
    reg.add((sx,sy))
    def fit(R):
        P=np.array(list(R)); A=np.c_[np.ones(len(P)),P[:,0],P[:,1]]; c,*_=np.linalg.lstsq(A,Z[P[:,1],P[:,0]],rcond=None); return c
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

print(f"{'#':>2} {'lat,lon':>21} {'cert':>5} | {'seed move':>9} | {'grown surface':>20} {'rough':>6} | ref_plane vs geoid (mm)")
for i in range(len(E)):
    if not ing[i]:
        print(f"{i+1:>2} {lat[i]:8.4f},{lon[i]:8.4f} {str(bool(cert[i])):>5} | OUTSIDE grid"); continue
    fs=find_flat_seed(ix0[i],iy0[i])
    if fs is None:
        print(f"{i+1:>2} {lat[i]:8.4f},{lon[i]:8.4f} {str(bool(cert[i])):>5} | no flat patch within 10 m"); continue
    sx,sy,sres=fs; movem=RES*np.hypot(sx-ix0[i],sy-iy0[i])
    reg=grow(sx,sy)
    if reg is None or len(reg)<3:
        print(f"{i+1:>2} {lat[i]:8.4f},{lon[i]:8.4f} {str(bool(cert[i])):>5} | seed@{movem:.0f}m but no growth"); continue
    P=np.array(list(reg)); jx,jy=P[:,0],P[:,1]
    ext_e=(jx.max()-jx.min()+1)*RES; ext_n=(jy.max()-jy.min()+1)*RES
    z=Z[jy,jx]; A=np.c_[np.ones(len(z)),jx,jy]; c,*_=np.linalg.lstsq(A,z,rcond=None); rough=np.std(z-A@c)*1000
    dr=DOD[jy,jx]; dg=DODG[jy,jx]; fin=np.isfinite(dr)&np.isfinite(dg); dr=dr[fin]; dg=dg[fin]
    mr=np.median(dr); mg=np.median(dg)
    ser=1.2533*1.4826*np.median(np.abs(dr-mr))/np.sqrt(dr.size)
    seg=1.2533*1.4826*np.median(np.abs(dg-mg))/np.sqrt(dg.size)
    print(f"{i+1:>2} {lat[i]:8.4f},{lon[i]:8.4f} {str(bool(cert[i])):>5} | {movem:5.0f}m | "
          f"ref {mr:+6.1f}±{ser:3.0f} | geoid {mg:+6.1f}±{seg:3.0f} | shift {mg-mr:+5.1f} (n={dr.size})")
