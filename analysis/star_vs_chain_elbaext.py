#!/usr/bin/env python3
"""Run the STAR frame on elbaext and compare it with the CHAIN that align_swaths solves.

Reports only. Touches no product and no pipeline code.

Both are fitted on the SAME raw gen1 cloud (the CSF cache, before any registration) and
the SAME gen2 reference, so the two solutions are directly comparable -- which is the
whole point, and the thing this session kept getting wrong by comparing numbers from
different estimators on different populations.

    ./lidar-icp/bin/python analysis/star_vs_chain_elbaext.py \
        --tile elbaext_delong --csf data/csf_cache/elbaext.las --valley-top-m 230.0
"""
import argparse
import json

import numpy as np
from scipy.ndimage import gaussian_filter

from lidar_diff_icp import groundest, terrain
from lidar_diff_icp.starframe import fit_star_shifts

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True)
ap.add_argument("--csf", required=True)
ap.add_argument("--valley-top-m", type=float, required=True,
                help="CITED per tile, never chosen here; elbaext = 230.0 (refcells.py)")
ap.add_argument("--ground-q", type=float, default=0.50)
ap.add_argument("--iterations", type=int, default=2)
A = ap.parse_args()

import laspy

D = f"data/derived/{A.tile}"
c = json.load(open(f"{D}/corrections.json"))
b = c["bounds"]; res = float(c["res_m"]); X0, Y0 = b[0], b[1]
Zref = np.load(f"{D}/z_after_differenced.npy"); ny, nx = Zref.shape
grid = (X0, Y0, res, nx, ny)
stable = np.load(f"{D}/stable.npy").astype(bool)

Zf = terrain.terrain_masks(Zref, res, valley_top_m=A.valley_top_m)["filled"]
Zreg = gaussian_filter(Zf, 1.2)
plane = (Zreg.ravel(), np.gradient(Zreg, res, axis=1).ravel(),
         np.gradient(Zreg, res, axis=0).ravel())


def ground_of(xx, yy, zz):
    return groundest.estimate_ground(xx, yy, zz, grid, A.ground_q, "slope_normal",
                                     plane=plane)


f = laspy.read(A.csf)
x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
sid = np.asarray(f.point_source_id)
ground = np.asarray(f.classification) == 2
print(f"{A.tile}: {ground.sum():,} CSF ground points, "
      f"{len(np.unique(sid[ground]))} swaths, {stable.sum():,} stable cells\n")

corr, rows = fit_star_shifts(x, y, z, sid, ground, Zref, ground_of, grid, stable,
                             iterations=A.iterations, verbose=True)

CH = c["per_swath_internal_alignment_dxdydz_m"]
print("\nSTAR vs CHAIN, same cloud, same reference")
print(f"  {'swath':>7}{'star dx':>10}{'chain dx':>10}{'star dy':>10}{'chain dy':>10}"
      f"{'star dz':>10}{'chain dz':>10}   (mm)")
sx, cx, sz_, cz_ = [], [], [], []
for s in sorted(corr):
    st_ = corr[s]; ch = CH.get(str(s))
    if ch is None:
        continue
    sx.append(1000*st_[0]); cx.append(1000*ch[0])
    sz_.append(1000*st_[2]); cz_.append(1000*ch[2])
    print(f"  {s:>7}{1000*st_[0]:>10.0f}{1000*ch[0]:>10.0f}{1000*st_[1]:>10.0f}"
          f"{1000*ch[1]:>10.0f}{1000*st_[2]:>10.1f}{1000*ch[2]:>10.1f}")
sx, cx = np.array(sx), np.array(cx)
sz_, cz_ = np.array(sz_), np.array(cz_)
print(f"\n  dx spread:  star {sx.max()-sx.min():>8.0f} mm    chain {cx.max()-cx.min():>8.0f} mm")
print(f"  dz spread:  star {sz_.max()-sz_.min():>8.1f} mm    chain {cz_.max()-cz_.min():>8.1f} mm")
print(f"  corr(star dx, chain dx) = {np.corrcoef(sx,cx)[0,1]:+.3f}")
print(f"  corr(star dz, chain dz) = {np.corrcoef(sz_,cz_)[0,1]:+.3f}")
print("\n  NOTE the gauge: only DIFFERENCES between swaths are meaningful in either")
print("  solution. Compare the SPREAD and the SHAPE, never the level.")
