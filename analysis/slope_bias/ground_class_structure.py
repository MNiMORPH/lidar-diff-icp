#!/usr/bin/env python3
"""GROUND-CLASS STRUCTURE PROBE (Task 4) — is there low vegetation inside the class-2
ground returns that biases the gridded ground HIGH?

Competing/complementary hypothesis to pure sparsity: the +tan(slope) bias may be low veg
(grass, brush, leaf litter, understory) MISCLASSIFIED as ground. That biases the ground
surface up and MORE points won't fix it (same contamination) -- but a robust LOW estimator
would. Signature: within a cell, class-2 residuals to the local ground PLANE have a
positive-skewed tail (a core near true ground + points sitting 0.15-1 m above it), and
it is stronger on forested cells than open ones at the same slope.

SLOPE-AWARE BY CONSTRUCTION: we fit a plane to the class-2 points IN EACH 5 m CELL and
characterize residuals to that plane, so within-cell tilt is removed at the source (not
differenced out afterward). Run on the FULL-density cloud (real ground structure).

  pass 1: per-cell normal equations -> local ground plane (a + b*u + c*v)
  pass 2: residual r = z - plane; per-cell tail fractions, skew, RMS, max

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/ground_class_structure.py
"""
import laspy, numpy as np
import warnings; warnings.filterwarnings("ignore")

FN = "data/after/3dep2021_fulldensity.laz"
RES = 5.0; STEEP = 15.0; LOWPEN = 0.25; MINN = 12; CHUNK = 20_000_000

h = laspy.open(FN).header
x0, y0 = h.mins[0], h.mins[1]
nx = int((h.maxs[0] - x0) / RES) + 1
ny = int((h.maxs[1] - y0) / RES) + 1
NC = nx * ny

def cu(px, py):
    ix = np.clip(((px - x0) / RES).astype(np.int64), 0, nx - 1)
    iy = np.clip(((py - y0) / RES).astype(np.int64), 0, ny - 1)
    u = px - (x0 + ix * RES); v = py - (y0 + iy * RES)      # local coords in cell (0..RES)
    return iy * nx + ix, u.astype(np.float64), v.astype(np.float64)

# ---- PASS 1: per-cell plane normal equations (class-2 only) + total returns -------------
S = {k: np.zeros(NC) for k in ("n","u","v","uu","uv","vv","z","uz","vz")}
tot = np.zeros(NC, np.int64)
print("PASS 1: ground-plane fit + totals ...", flush=True)
with laspy.open(FN) as fh:
    for pts in fh.chunk_iterator(CHUNK):
        cl = np.asarray(pts.classification); keep = cl != 7
        px = np.asarray(pts.x)[keep]; py = np.asarray(pts.y)[keep]; cl = cl[keep]
        c_all, _, _ = cu(px, py)
        tot += np.bincount(c_all, minlength=NC)
        g2 = cl == 2
        c, u, v = cu(px[g2], py[g2]); z = np.asarray(pts.z)[keep][g2].astype(np.float64)
        S["n"]  += np.bincount(c, minlength=NC)
        S["u"]  += np.bincount(c, u, NC);  S["v"]  += np.bincount(c, v, NC)
        S["uu"] += np.bincount(c, u*u, NC); S["uv"] += np.bincount(c, u*v, NC)
        S["vv"] += np.bincount(c, v*v, NC); S["z"]  += np.bincount(c, z, NC)
        S["uz"] += np.bincount(c, u*z, NC); S["vz"] += np.bincount(c, v*z, NC)
        print("  ...chunk", flush=True)

# solve 3x3 per cell for [a,b,c] where enough points
n = S["n"]; ok = n >= MINN
A = np.zeros((NC, 3, 3)); b = np.zeros((NC, 3))
A[:,0,0]=n; A[:,0,1]=S["u"]; A[:,0,2]=S["v"]
A[:,1,0]=S["u"]; A[:,1,1]=S["uu"]; A[:,1,2]=S["uv"]
A[:,2,0]=S["v"]; A[:,2,1]=S["uv"]; A[:,2,2]=S["vv"]
b[:,0]=S["z"]; b[:,1]=S["uz"]; b[:,2]=S["vz"]
coef = np.zeros((NC, 3))
idx = np.where(ok)[0]
coef[idx] = np.linalg.solve(A[idx], b[idx])                 # a,b,c per cell
print(f"  fitted planes on {ok.sum():,} cells (n>={MINN})", flush=True)

# ---- PASS 2: residuals to the local plane -> tails, skew, RMS ---------------------------
R = {k: np.zeros(NC) for k in ("cnt","g15","g30","g50","sr","sr2","sr3","rmax")}
R["rmax"][:] = -1e9
print("PASS 2: residual structure ...", flush=True)
with laspy.open(FN) as fh:
    for pts in fh.chunk_iterator(CHUNK):
        cl = np.asarray(pts.classification); g2 = (cl == 2)
        if not g2.any(): continue
        px = np.asarray(pts.x)[g2]; py = np.asarray(pts.y)[g2]
        z = np.asarray(pts.z)[g2].astype(np.float64)
        c, u, v = cu(px, py)
        m = ok[c]
        c=c[m]; u=u[m]; v=v[m]; z=z[m]
        a=coef[c,0]; bb=coef[c,1]; cc=coef[c,2]
        r = z - (a + bb*u + cc*v)
        R["cnt"] += np.bincount(c, minlength=NC)
        R["g15"] += np.bincount(c, (r>0.15).astype(float), NC)
        R["g30"] += np.bincount(c, (r>0.30).astype(float), NC)
        R["g50"] += np.bincount(c, (r>0.50).astype(float), NC)
        R["sr"]  += np.bincount(c, r, NC); R["sr2"] += np.bincount(c, r*r, NC)
        R["sr3"] += np.bincount(c, r**3, NC)
        np.maximum.at(R["rmax"], c, r)
        print("  ...chunk", flush=True)

# per-cell stats
cnt = R["cnt"]; good = cnt >= MINN
mean = np.divide(R["sr"], np.maximum(cnt,1))
var  = np.divide(R["sr2"], np.maximum(cnt,1)) - mean**2
std  = np.sqrt(np.clip(var, 0, None))
m3   = np.divide(R["sr3"], np.maximum(cnt,1)) - 3*mean*var - mean**3
skew = np.divide(m3, np.maximum(std**3, 1e-9))
frac15 = np.divide(R["g15"], np.maximum(cnt,1))
frac30 = np.divide(R["g30"], np.maximum(cnt,1))

# zones (slope from cached 2 m surface, block-mean to 5 m)
p1 = np.load("data/derived/elba/recov_pass1.npz"); s2 = p1["slope"]; f = int(RES/2)
ny2,nx2 = s2.shape; Ny,Nx = ny*f, nx*f
pad = np.full((Ny,Nx), np.nan, np.float32); pad[:min(Ny,ny2),:min(Nx,nx2)] = s2[:min(Ny,ny2),:min(Nx,nx2)]
slope = np.nanmean(pad.reshape(ny,f,nx,f), axis=(1,3)).ravel()
pen = np.divide(S["n"], np.maximum(tot,1))                   # ground fraction per cell

zones = {
 "steepForest": good & (slope>=STEEP) & (pen<LOWPEN),
 "steepOpen":   good & (slope>=STEEP) & (pen>=0.40),
 "flatForest":  good & (slope<STEEP)  & (pen<LOWPEN),
 "flatOpen":    good & (slope<STEEP)  & (pen>=0.40),
}
print("\n=== class-2 residual structure to the LOCAL GROUND PLANE (slope removed) ===")
print("  positive skew + tail fraction = low veg sitting IN the ground class (biases high)")
print(f"{'zone':13s} {'cells':>8} {'medRMS':>7} {'medSkew':>8} "
      f"{'med%>15cm':>9} {'med%>30cm':>9} {'medMax':>7}")
for zk, zm in zones.items():
    if zm.sum()==0: print(f"{zk:13s} none"); continue
    print(f"{zk:13s} {zm.sum():8,} {np.median(std[zm]):7.3f} {np.median(skew[zm]):8.2f} "
          f"{100*np.median(frac15[zm]):8.1f}% {100*np.median(frac30[zm]):8.1f}% "
          f"{np.median(R['rmax'][zm]):7.2f}")
print("\nRead: if steepForest (and flatForest) show markedly higher RMS/skew/tail than the "
      "OPEN controls at the SAME slope, low veg is contaminating the ground class -> a "
      "robust LOW estimator (or reclassify) fixes it; density alone does not.")
