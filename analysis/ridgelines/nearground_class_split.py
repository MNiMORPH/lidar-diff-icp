#!/usr/bin/env python3
"""Why gen2's ground sits BELOW the median of gen2's own near-ground column.

`z_after` is the per-cell MEDIAN of the *ground-classified* returns (`ground_q = 0.50`,
`after_ground = "class2"`, verified in pipeline.py). The near-ground column, by contrast,
holds ALL returns. Everything the classifier excludes lies above the surface, so if a
fraction ``w_g`` of the window is ground and the ground returns are symmetric about the
surface,

    rank of the surface in the FULL column = 0.5 * w_g

This script tests that prediction directly by rebuilding gen2's near-ground column split
into class-2 and non-class-2, and checking three things:
  1. the surface sits at rank ~0.50 of the CLASS-2 column (the anchor is what we think);
  2. the non-class-2 mass is almost entirely ABOVE the surface (one-sided contamination);
  3. 0.5 * w_g reproduces the measured full-column rank.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/nearground_class_split.py --tile data/derived/elba_fulldensity \\
        --gen2 data/after/3dep2021_fulldensity.laz
"""
import argparse, json, os
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--gen2", default="data/after/3dep2021_fulldensity.laz")
ap.add_argument("--chunk", type=int, default=3_000_000)
A = ap.parse_args()

STRATA = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
          ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01)]


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
cells = cube["cells"]; edges = cube["edges"]
zlo = float(cube["zlo"]); dz = float(cube["dz"]); NZ = edges.size - 1

_zf = np.load(f"{A.tile}/z_after.npy").copy()
_m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
lut = np.full(NY*NX, -1, np.int32); lut[cells] = np.arange(cells.size, dtype=np.int32)

Hg = np.zeros((cells.size, NZ), np.int32)      # class 2
Hn = np.zeros((cells.size, NZ), np.int32)      # everything else (noise class 7 dropped)
with laspy.open(A.gen2) as f:
    for pts in f.chunk_iterator(A.chunk):
        cl = np.asarray(pts.classification)
        good = cl != 7
        x = np.asarray(pts.x)[good]; y = np.asarray(pts.y)[good]; z = np.asarray(pts.z)[good]
        g2 = cl[good] == 2
        ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
        ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        cid = lut[iy[ing]*NX + ix[ing]]
        sub = cid >= 0
        if not sub.any():
            continue
        cc = (iy[ing]*NX + ix[ing])[sub]
        xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
        h = (z[ing][sub] - (zflat[cc] + gxf[cc]*(x[ing][sub]-xc) + gyf[cc]*(y[ing][sub]-yc))) \
            / nnorm[cc]
        zi = np.floor((h - zlo)/dz).astype(np.int64)
        m = (zi >= 0) & (zi < NZ)
        gg = g2[ing][sub][m]
        np.add.at(Hg, (cid[sub][m][gg], zi[m][gg]), 1)
        np.add.at(Hn, (cid[sub][m][~gg], zi[m][~gg]), 1)
print(f"gen2 window: {Hg.sum():,} class-2 + {Hn.sum():,} other returns", flush=True)

cover = np.load(f"{A.tile}/canopy_cover_pfs.npy").ravel()[cells]
stable, _ = reference_cells(A.tile, cells=cells, slope_max=90.0)
Cg = np.cumsum(Hg, 1).astype(float); Cn = np.cumsum(Hn, 1).astype(float)
ng = Cg[:, -1]; nn = Cn[:, -1]; nt = ng + nn
k0 = int(round((0.0 - zlo)/dz))
ok = stable & (nt > 0) & (ng > 0) & np.isfinite(cover)

f_full = (Cg[:, k0-1] + Cn[:, k0-1]) / np.maximum(nt, 1)     # rank of surface, ALL returns
f_g = Cg[:, k0-1] / np.maximum(ng, 1)                        # rank of surface, class-2 only
w_g = ng / np.maximum(nt, 1)                                 # ground fraction of the window
n_above = 1 - Cn[:, k0-1] / np.maximum(nn, 1)                # of non-ground, fraction above

print(f"\n{os.path.basename(A.tile)}: is the surface's rank just half the ground fraction?")
print(f"{'stratum':14s} {'cells':>7s} | {'rank(class2)':>12s} {'w_g':>5s} {'0.5*w_g':>8s} "
      f"{'rank(all)':>9s} | {'non-grd above':>13s}")
for nm, lo, hi in STRATA + [("ALL", -9, 9)]:
    m = ok & (cover > lo) & (cover <= hi) if nm != "ALL" else ok
    if m.sum() < 100:
        continue
    mm = m & (nn > 0)
    print(f"{nm:14s} {m.sum():7,d} | {np.median(f_g[m]):12.2f} {np.median(w_g[m]):5.2f} "
          f"{0.5*np.median(w_g[m]):8.2f} {np.median(f_full[m]):9.2f} | "
          f"{np.median(n_above[mm]):13.2f}")
