#!/usr/bin/env python3
"""Per-cell near-ground return HISTOGRAMS for both epochs -- the object, not a summary.

p05 and p50 were a guess at which statistics matter. This stores the whole near-ground
distribution per cell, per epoch, at 0.02 m, so any statistic (mode, skew, width, the
difference of the two CDFs, the height where the epochs diverge) can be tried afterwards
without re-reading a multi-GB cloud. One streaming pass per epoch; the cube is small.

Cells are the reference population -- divides with |curv_laplacian| <= --curv-max -- since
that is where the offset is measured. ALL returns in those cells are kept: no
classification, no ground/vegetation decision, no height threshold beyond the window.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/nearground_cells.py --tile data/derived/elba_fulldensity \
        --gen1 data/before/4342-29-64.laz --gen2 data/after/3dep2021_fulldensity.laz
"""
import argparse, json, os
import numpy as np, laspy

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--gen1", required=True); ap.add_argument("--gen2", required=True)
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--zlo", type=float, default=-1.0)
ap.add_argument("--zhi", type=float, default=2.0)
ap.add_argument("--dz", type=float, default=0.02)
ap.add_argument("--chunk", type=int, default=3_000_000)
ap.add_argument("--out", default=None)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
# --out names a PRODUCT OF THIS TILE, so a bare filename resolves inside the tile
# directory -- the same place the default goes. It used to be handed to np.savez_compressed
# verbatim, so `--out nearground_cells_sn.npz` (exactly what lidar-diff-workflow --plan
# prints) wrote 5 MB into the current working directory, and every consumer then failed
# with FileNotFoundError looking in the tile. The command succeeded and the product went
# somewhere nothing reads: worse than an error. Pass a path with a directory component to
# put it anywhere else.
OUT = A.out or "nearground_cells.npz"
if not os.path.dirname(OUT):
    OUT = os.path.join(A.tile, OUT)


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r, int(j.get("ny") or round((b[3]-b[1])/r)), \
                   int(j.get("nx") or round((b[2]-b[0])/r))
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES, NY, NX = grid(A.tile)
# SLOPE-NORMAL residual, the same quantity as d_mm: subtract the local PLANE (not just the
# cell-centre elevation) and divide by |n|. Without the plane term, intra-cell terrain relief
# (+-3.54 m x tan(slope) = +-2 m at 30 deg) swamps the near-ground PDF on any slope. Low
# vegetation is only compressed by cos(slope) in this frame -- a few cm -- and staying in one
# frame keeps the PDF commensurate with the offset it is meant to predict.
_zg = np.load(f"{A.tile}/z_after.npy")
_zf = _zg.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    from scipy.ndimage import distance_transform_edt
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
zflat = _zf.ravel()
ridge = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
curv = np.abs(np.load(f"{A.tile}/curv_laplacian.npy").ravel())
keep_cell = ridge & np.isfinite(curv) & (curv <= A.curv_max) & np.isfinite(zflat)
cells = np.flatnonzero(keep_cell)
lut = np.full(NY*NX, -1, np.int32); lut[cells] = np.arange(cells.size, dtype=np.int32)
edges = np.arange(A.zlo, A.zhi + 0.5*A.dz, A.dz); NZ = edges.size - 1
print(f"{TILE}: {cells.size:,} reference cells, {NZ} height bins of {A.dz:g} m "
      f"({A.zlo:g} to {A.zhi:g} m)")


def cube(path, label):
    H = np.zeros((cells.size, NZ), np.int32)
    n_in = 0
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(A.chunk):
            cl = np.asarray(pts.classification)
            good = cl != 7
            x = np.asarray(pts.x)[good]; y = np.asarray(pts.y)[good]; z = np.asarray(pts.z)[good]
            ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
            ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
            cid = lut[iy[ing]*NX + ix[ing]]
            sub = cid >= 0
            if not sub.any():
                continue
            cc = (iy[ing]*NX + ix[ing])[sub]
            xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
            xs = x[ing][sub]; ys = y[ing][sub]
            h = (z[ing][sub] - (zflat[cc] + gxf[cc]*(xs - xc) + gyf[cc]*(ys - yc))) / nnorm[cc]
            zi = np.floor((h - A.zlo)/A.dz).astype(np.int64)
            m = (zi >= 0) & (zi < NZ)
            np.add.at(H, (cid[sub][m], zi[m]), 1)
            n_in += int(m.sum())
    print(f"  {label}: {n_in:,} returns landed in the window", flush=True)
    return H


H1 = cube(A.gen1, "gen1")
H2 = cube(A.gen2, "gen2")
np.savez_compressed(OUT, cells=cells, edges=edges, H1=H1, H2=H2,
                    zlo=A.zlo, zhi=A.zhi, dz=A.dz, curv_max=A.curv_max)
print(f"wrote {OUT}  ({H1.sum():,} gen1 + {H2.sum():,} gen2 returns over {cells.size:,} cells)")
