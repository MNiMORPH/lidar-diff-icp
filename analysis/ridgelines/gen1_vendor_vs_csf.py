#!/usr/bin/env python3
"""The same classifier question on the epoch we are actually correcting: gen1.

gen1 ships with a vendor TerraScan classification (its LAS header still says so) and we
re-classify it ourselves with PDAL ELM + CSF. If the two disagree in a way that tracks
canopy cover, that is a classification artefact masquerading as ground change.

Both grounds are reduced identically -- per-cell MEDIAN of the slope-normal residual to
the same gen2 plane -- so the difference is the classifier and nothing else.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/gen1_vendor_vs_csf.py
"""
import argparse, json
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--raw", default="data/before/4342-29-64.laz")
ap.add_argument("--csf", default="data/csf_cache/elba.las")
ap.add_argument("--min-pts", type=int, default=3)
A = ap.parse_args()

j = json.load(open(f"{A.tile}/corrections.json"))
b = j["bounds"]; RES = float(j["res_m"]); X0, Y0 = b[0], b[1]
zf = np.load(f"{A.tile}/z_after.npy"); NY, NX = zf.shape
_zf = zf.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)


def cell_median(path, class2_only):
    f = laspy.read(path)
    cl = np.asarray(f.classification)
    sel = (cl == 2) if class2_only else np.ones(cl.size, bool)
    x = np.asarray(f.x)[sel]; y = np.asarray(f.y)[sel]; z = np.asarray(f.z)[sel]
    ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
    ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
    cc = iy[ing]*NX + ix[ing]
    xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
    h = (z[ing] - (zflat[cc] + gxf[cc]*(x[ing]-xc) + gyf[cc]*(y[ing]-yc)))/nnorm[cc]*1000.0
    o = np.argsort(cc, kind="stable"); cs = cc[o]; hs = h[o]
    u, st = np.unique(cs, return_index=True); sp = np.r_[st[1:], cs.size]
    med = np.full(NY*NX, np.nan); n = np.zeros(NY*NX, np.int64)
    for uu, a, bb in zip(u, st, sp):
        n[uu] = bb - a
        if bb - a >= A.min_pts:
            med[uu] = np.median(hs[a:bb])
    return med, n


v_med, v_n = cell_median(A.raw, class2_only=True)     # vendor TerraScan ground
c_med, c_n = cell_median(A.csf, class2_only=False)    # our ELM + CSF ground (file is ground only)
cover = np.load(f"{A.tile}/canopy_cover_pfs.npy").ravel()
stable, _ = reference_cells(A.tile, slope_max=90.0)
ok = stable & np.isfinite(v_med) & np.isfinite(c_med) & np.isfinite(cover)
print(f"cells with both grounds (>= {A.min_pts} pts each): {ok.sum():,}; "
      f"median pts/cell vendor {np.median(v_n[ok]):.0f}, CSF {np.median(c_n[ok]):.0f}")
print("\ngen1 ground: OUR CSF minus the VENDOR TerraScan class-2, per cover stratum (mm)")
print(f"  {'stratum':14s} {'cells':>7s} | {'csf - vendor':>12s} | {'vendor n':>8s} {'csf n':>6s}")
for nm, lo, hi in [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
                   ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01), ("ALL", -9, 9)]:
    m = ok & (cover > lo) & (cover <= hi) if nm != "ALL" else ok
    if m.sum() < 50:
        continue
    print(f"  {nm:14s} {m.sum():7,d} | {np.median(c_med[m]-v_med[m]):+12.1f} | "
          f"{np.median(v_n[m]):8.0f} {np.median(c_n[m]):6.0f}")
