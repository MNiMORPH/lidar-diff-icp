#!/usr/bin/env python3
"""Does classifying BOTH epochs the same way close the gen1-gen2 gap?

Compares, on the reference divide cells of the CSF patch:
  * `z_after`  -- gen2 ground as the pipeline has it (median of the VENDOR class-2)
  * `z_csf`    -- gen2 ground from OUR CSF (identical filter and parameters to gen1)
and asks whether the difference tracks canopy cover, i.e. whether the cover-dependent
offset we have been trying to correct is partly a difference in classification RULE.

Everything is in the cube's slope-normal frame, so `h = 0` is `z_after` by construction
and the per-cell median of the CSF ground's `h` IS `z_csf - z_after`.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/gen2_csf_compare.py
"""
import argparse, json, os
import numpy as np, laspy, pyarrow.parquet as pq
from scipy.ndimage import distance_transform_edt

from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--csf", default="data/csf_cache/elba_gen2_patch.las")
ap.add_argument("--min-pts", type=int, default=5)
A = ap.parse_args()

j = json.load(open(f"{A.tile}/corrections.json"))
b = j["bounds"]; RES = float(j["res_m"])
X0, Y0 = b[0], b[1]
zf = np.load(f"{A.tile}/z_after.npy"); NY, NX = zf.shape
_zf = zf.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)

g = laspy.read(A.csf)
x = np.asarray(g.x); y = np.asarray(g.y); z = np.asarray(g.z)
ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
cc = (iy[ing]*NX + ix[ing])
xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
h = (z[ing] - (zflat[cc] + gxf[cc]*(x[ing]-xc) + gyf[cc]*(y[ing]-yc)))/nnorm[cc] * 1000.0

order = np.argsort(cc, kind="stable")
cs = cc[order]; hs = h[order]
uniq, start = np.unique(cs, return_index=True)
stop = np.r_[start[1:], cs.size]
dz_csf = np.full(NY*NX, np.nan); npts = np.zeros(NY*NX, np.int64)
for u, a, bb in zip(uniq, start, stop):
    if bb - a >= A.min_pts:
        dz_csf[u] = np.median(hs[a:bb])
    npts[u] = bb - a

cover = np.load(f"{A.tile}/canopy_cover_pfs.npy").ravel()
dod = np.load(f"{A.tile}/dod.npy").ravel() * 1000.0          # gen2 - gen1, + = rose
stable, _ = reference_cells(A.tile, slope_max=90.0)
ok = stable & np.isfinite(dz_csf) & np.isfinite(cover) & np.isfinite(dod)
print(f"patch reference cells with >= {A.min_pts} CSF ground pts: {ok.sum():,} "
      f"(median {np.median(npts[ok]):.0f} pts/cell)")
print("\ngen2 ground: OUR CSF minus the VENDOR class-2 median, and what it does to the DoD")
print(f"  {'stratum':14s} {'cells':>6s} | {'z_csf - z_after':>15s} | {'DoD now':>8s} "
      f"{'DoD if CSF':>11s} | {'change':>7s}")
rows = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
        ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01)]
for nm, lo, hi in rows + [("ALL", -9, 9)]:
    m = ok & (cover > lo) & (cover <= hi) if nm != "ALL" else ok
    if m.sum() < 20:
        continue
    d = np.median(dz_csf[m]); D = np.median(dod[m])
    print(f"  {nm:14s} {m.sum():6,d} | {d:+15.1f} | {D:+8.1f} {D + d:+11.1f} | {d:+7.1f}")
print("\n(DoD = gen2 - gen1, + = elevation rose. Moving gen2's ground down by |z_csf-z_after|")
print(" moves the DoD down by the same amount. A cover-dependent shift means part of the")
print(" apparent hillslope gain is a CLASSIFIER difference, not terrain.)")
