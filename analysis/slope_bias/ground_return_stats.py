#!/usr/bin/env python3
"""GROUND-CLASSIFIED-RETURN STATISTICS in 5x5 m cells (Task 4, veg-lift hypothesis).

Question: do the class-2 ("ground") returns show a low-vegetation column, and is it worse
in gen2 (2021-05-01 green-up) than gen1 (2008-11 dormant) on the SAME forested cells? That
is the signature of a seasonal ground-cover lift that biases the gridded ground high.

Slope removed at the source: per 5 m cell we fit a plane to the class-2 points and look at
residuals r = z - plane (r has zero mean per cell by LSQ). A veg column shows as an UPWARD
tail (positive p95/p99, fraction of returns well above the plane), stronger in gen2-forest
than gen1-forest and than open cells. Same AOI grid as dod.npy / slope.npy / bedrock.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/ground_return_stats.py
"""
import laspy, numpy as np

RES = 5.0
X0, Y0 = 577492.8, 4882737.6
slope = np.load("data/derived/elba_fulldensity/slope.npy")
ny, nx = slope.shape; NC = ny * nx
GEN1 = "data/before/4342-29-64.laz"
GEN2 = "data/after/3dep2021_fulldensity.laz"
HED = np.arange(-0.60, 1.201, 0.01)                 # residual histogram edges (m)
HC = 0.5 * (HED[:-1] + HED[1:])

def cu(px, py):
    ix = ((px - X0) / RES).astype(np.int64); iy = ((py - Y0) / RES).astype(np.int64)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    u = px - (X0 + ix * RES); v = py - (Y0 + iy * RES)
    return iy * nx + ix, u, v, ok

def read_class2(fn, chunk=20_000_000):
    """yield (cell, u, v, z) for class-2 returns, chunked."""
    with laspy.open(fn) as fh:
        for pts in fh.chunk_iterator(chunk):
            cl = np.asarray(pts.classification); m = cl == 2
            if not m.any(): continue
            px = np.asarray(pts.x)[m]; py = np.asarray(pts.y)[m]; z = np.asarray(pts.z)[m].astype(np.float64)
            c, u, v, ok = cu(px, py)
            yield c[ok], u[ok], v[ok], z[ok]

def totals(fn, chunk=20_000_000):
    t = np.zeros(NC, np.int64)
    with laspy.open(fn) as fh:
        for pts in fh.chunk_iterator(chunk):
            c, _, _, ok = cu(np.asarray(pts.x), np.asarray(pts.y))
            t += np.bincount(c[ok], minlength=NC)
    return t

def plane_fit(fn):
    S = {k: np.zeros(NC) for k in ("n","u","v","uu","uv","vv","z","uz","vz")}
    for c,u,v,z in read_class2(fn):
        S["n"]+=np.bincount(c,minlength=NC); S["u"]+=np.bincount(c,u,NC); S["v"]+=np.bincount(c,v,NC)
        S["uu"]+=np.bincount(c,u*u,NC); S["uv"]+=np.bincount(c,u*v,NC); S["vv"]+=np.bincount(c,v*v,NC)
        S["z"]+=np.bincount(c,z,NC); S["uz"]+=np.bincount(c,u*z,NC); S["vz"]+=np.bincount(c,v*z,NC)
    n=S["n"]; ok=n>=8
    A=np.zeros((NC,3,3)); b=np.zeros((NC,3))
    A[:,0,0]=n;A[:,0,1]=S["u"];A[:,0,2]=S["v"];A[:,1,0]=S["u"];A[:,1,1]=S["uu"];A[:,1,2]=S["uv"]
    A[:,2,0]=S["v"];A[:,2,1]=S["uv"];A[:,2,2]=S["vv"];b[:,0]=S["z"];b[:,1]=S["uz"];b[:,2]=S["vz"]
    coef=np.zeros((NC,3)); idx=np.where(ok)[0]; coef[idx]=np.linalg.solve(A[idx],b[idx])
    return coef, n

def residual_stats(fn, coef, zones):
    """per-cell moments + per-zone residual histogram."""
    M = {k: np.zeros(NC) for k in ("n","r2","r3")}
    H = {z: np.zeros(len(HC)) for z in zones}
    zid = np.full(NC, -1, np.int8)
    for i,(zk,zm) in enumerate(zones.items()): zid[zm]=i
    zkeys=list(zones)
    for c,u,v,z in read_class2(fn):
        r = z - (coef[c,0]+coef[c,1]*u+coef[c,2]*v)
        M["n"]+=np.bincount(c,minlength=NC); M["r2"]+=np.bincount(c,r*r,NC); M["r3"]+=np.bincount(c,r**3,NC)
        zi=zid[c]
        for i,zk in enumerate(zkeys):
            sel=zi==i
            if sel.any(): H[zk]+=np.histogram(r[sel],bins=HED)[0]
    return M, H

def pct(h, q):
    cx=np.cumsum(h); tot=cx[-1]
    if tot==0: return np.nan
    return HC[np.searchsorted(cx, q*tot)]
def fabove(h, thr):
    tot=h.sum();  return np.nan if tot==0 else h[HC>thr].sum()/tot

print("slope + gen2 penetration ...", flush=True)
t2 = totals(GEN2)
coef2, n2 = plane_fit(GEN2)
pen = np.divide(n2, np.maximum(t2,1))                       # gen2 ground fraction per cell
sl = slope.ravel()
zones = {
 "forestSteep": (pen<0.25)&(sl>=15)&(sl<40)&(n2>=8),
 "forestFlat":  (pen<0.25)&(sl<15)&(n2>=8),
 "openSteep":   (pen>=0.40)&(sl>=15)&(sl<40)&(n2>=8),
 "openFlat":    (pen>=0.40)&(sl<15)&(n2>=8),
}
for zk,zm in zones.items(): print(f"  zone {zk}: {int(zm.sum())} cells")

print("gen2 residual stats ...", flush=True)
M2,H2 = residual_stats(GEN2, coef2, zones)
print("gen1 plane fit + residual stats (same zones/cells) ...", flush=True)
coef1,n1 = plane_fit(GEN1)
M1,H1 = residual_stats(GEN1, coef1, zones)

def cellstd(M,zm):
    g=zm&(M["n"]>=8); v=np.divide(M["r2"][g],M["n"][g]); return np.median(np.sqrt(np.clip(v,0,None)))
def cellskew(M,zm):
    g=zm&(M["n"]>=12); v=np.divide(M["r2"][g],M["n"][g]); m3=np.divide(M["r3"][g],M["n"][g])
    return np.median(np.divide(m3,np.maximum(v**1.5,1e-9)))

print("\n=== class-2 return statistics per 5 m cell, slope removed (plane residual) ===")
print("  veg column => bigger cell-std, positive skew, high p95/p99, big frac>10/25 cm")
hdr=f"{'zone':12s} {'epoch':5s} {'ret/cell':>8} {'cellStd':>7} {'skew':>6} "\
    f"{'p50':>6} {'p90':>6} {'p95':>6} {'p99':>6} {'>10cm':>6} {'>25cm':>6}"
print(hdr)
for zk,zm in zones.items():
    for ep,(M,H,nn) in (("gen1",(M1,H1[zk],n1)),("gen2",(M2,H2[zk],n2))):
        medret=np.median(nn[zm&(nn>=8)]) if (zm&(nn>=8)).any() else np.nan
        print(f"{zk:12s} {ep:5s} {medret:8.0f} {cellstd(M,zm):7.3f} {cellskew(M,zm):6.2f} "
              f"{pct(H,.50):6.3f} {pct(H,.90):6.3f} {pct(H,.95):6.3f} {pct(H,.99):6.3f} "
              f"{100*fabove(H,.10):5.1f}% {100*fabove(H,.25):5.1f}%")
    print()
print("KEY: compare gen2-forest vs gen1-forest AND vs gen2-open. A seasonal ground-cover "
      "lift = gen2 forest has a taller/heavier upward tail than gen1 forest and than open.")
