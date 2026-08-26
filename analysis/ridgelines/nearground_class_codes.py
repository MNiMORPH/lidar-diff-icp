#!/usr/bin/env python3
"""What did each epoch's ground classifier actually DECIDE, above and below the surface?

The rank arithmetic says a quarter of gen2's non-ground near-ground returns sit BELOW the
surface the classifier drew. This asks what those returns are called, and does the
matching test for gen1 -- where WE do the classification (PDAL ELM + CSF, ground.py), so
the assumptions are ours to state: ``filters.elm`` marks isolated LOW points as noise
(class 7) and they are dropped BEFORE the cloth runs, and CSF then drapes a cloth from
above, which can only come to rest ON the point cloud. Both steps encode the same
assumption -- the lowest spatially COHERENT surface is the ground, and an isolated low
point is a blunder, not terrain.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/nearground_class_codes.py
"""
import argparse, json, os
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--gen2", default="data/after/3dep2021_fulldensity.laz")
ap.add_argument("--gen1-csf", default="data/csf_cache/elba.las")
ap.add_argument("--chunk", type=int, default=3_000_000)
A = ap.parse_args()

NAMES = {0: "never classified", 1: "unassigned", 2: "ground", 3: "low veg", 4: "med veg",
         5: "high veg", 6: "building", 7: "noise", 9: "water", 10: "rail", 11: "road",
         12: "overlap", 13: "wire guard", 14: "wire cond", 17: "bridge", 18: "high noise"}


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r, int(j.get("ny") or round((b[3]-b[1])/r)), \
                   int(j.get("nx") or round((b[2]-b[0])/r))
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES, NY, NX = grid(A.tile)
cube = np.load(f"{A.tile}/nearground_cells_sn.npz")
cells = cube["cells"]; zlo = float(cube["zlo"]); zhi = float(cube["zhi"])
_zf = np.load(f"{A.tile}/z_after.npy").copy()
_m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
keep = np.zeros(NY*NX, bool); keep[cells] = True


def heights(pts, drop_noise):
    cl = np.asarray(pts.classification)
    sel = (cl != 7) if drop_noise else np.ones(cl.size, bool)
    x = np.asarray(pts.x)[sel]; y = np.asarray(pts.y)[sel]; z = np.asarray(pts.z)[sel]
    ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
    ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
    cc = np.where(ing, iy*NX + ix, 0)
    ok = ing & keep[cc]
    cc = cc[ok]
    xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
    h = (z[ok] - (zflat[cc] + gxf[cc]*(x[ok]-xc) + gyf[cc]*(y[ok]-yc))) / nnorm[cc]
    w = (h >= zlo) & (h < zhi)
    return h[w], cl[sel][ok][w]


below = {}; above = {}
with laspy.open(A.gen2) as f:
    for pts in f.chunk_iterator(A.chunk):
        h, cl = heights(pts, drop_noise=False)
        for c in np.unique(cl):
            m = cl == c
            below[c] = below.get(c, 0) + int((h[m] < 0).sum())
            above[c] = above.get(c, 0) + int((h[m] >= 0).sum())
tot = sum(below.values()) + sum(above.values())
print(f"\ngen2 near-ground window, {tot:,} returns, by CLASSIFICATION and side of the surface")
print(f"  {'class':22s} {'below':>12s} {'above':>12s} {'% of window':>11s} {'% of class below':>17s}")
for c in sorted(set(below) | set(above), key=lambda k: -(below.get(k, 0)+above.get(k, 0))):
    b = below.get(c, 0); a = above.get(c, 0)
    if b + a < tot/1000:
        continue
    print(f"  {c:2d} {NAMES.get(c,'?'):18s} {b:12,d} {a:12,d} {100*(a+b)/tot:10.1f}% "
          f"{100*b/max(a+b,1):16.1f}%")
nong_b = sum(v for k, v in below.items() if k not in (2, 7))
nong_a = sum(v for k, v in above.items() if k not in (2, 7))
print(f"  non-ground (excl. noise): {100*nong_b/max(nong_b+nong_a,1):.1f}% of it lies BELOW "
      f"the surface the classifier drew")

if os.path.exists(A.gen1_csf):
    g = laspy.read(A.gen1_csf)
    h, _ = heights(g, drop_noise=False)
    hc = np.load(f"{A.tile}/nearground_cells_sn.npz")["H1"]
    dz = float(cube["dz"]); k0 = int(round((0.0 - zlo)/dz))
    allb = hc[:, :k0].sum(); alla = hc[:, k0:].sum()
    print(f"\ngen1: our CSF ground (post-ELM) vs ALL returns in the same window")
    print(f"  CSF ground kept {h.size:,} returns; the all-return column holds "
          f"{int(hc.sum()):,}  -> w_g = {h.size/max(int(hc.sum()),1):.2f}")
    print(f"  of the CSF GROUND returns, {100*(h < 0).mean():.1f}% lie below the gen2 surface; "
          f"of ALL gen1 returns, {100*allb/max(allb+alla,1):.1f}% do")
