#!/usr/bin/env python3
"""Pedagogical M3C2 vs. gridded-vertical comparison on bare sloped ground.

M3C2 (Lague et al. 2013) measures change between two point clouds NOT as a
vertical z(x,y) difference but along the LOCAL SURFACE NORMAL:

  1. at each "core point", estimate the surface normal from a neighborhood
     (radius = normal scale D) of the reference cloud;
  2. project a cylinder (radius d) along that normal and take the mean position
     of each cloud's points inside it;
  3. the M3C2 distance is the gap between those two means, measured ALONG the
     normal (not vertically);
  4. a level of detection (LoD) is formed per point from the roughness (point
     scatter) in each cylinder + a registration term.

We run it on bare, open, sloped cells (where the method should help most) and
compare its noise to the vertical gridded difference at the same places.
"""
import numpy as np
import laspy
import py4dgeo

from lidar_diff_icp import io, coreg
from lidar_diff_icp.swathdiff import _median_grid


def _nmad(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def main():
    res = 2.0
    X0, Y0, X1, Y1 = 577492.8, 4882737.6, 580035.0, 4886238.3
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    # --- 2008: internally align to swath 135, then rigid-tie to 2021 ---
    f = laspy.read("data/before/4342-29-64.laz")
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); cl8 = np.asarray(f.classification)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    terr = (nr8 == 1) & ~np.isin(cl8, [5, 6, 9])
    pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8),
                       np.zeros_like(ps8), io.MN_2008_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz

    g = laspy.read("data/after/3dep2021_fulltile.laz")
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    gm = (np.asarray(g.classification) == 2) & (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    Z21 = _median_grid(x2[gm], y2[gm], z2[gm], res, X0, Y0, nx, ny)
    inb = terr & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
    Z08 = _median_grid(xc[inb], yc[inb], zc[inb], res, X0, Y0, nx, ny)
    c = coreg.nuth_kaab(Z21, Z08, res)          # rigid tie
    xc += c.dx; yc += c.dy; zc += c.dz

    # --- point clouds for M3C2 ---
    p2008 = np.column_stack([xc[inb], yc[inb], zc[inb]]).astype(np.float64)
    p2021 = np.column_stack([x2[gm], y2[gm], z2[gm]]).astype(np.float64)
    print(f"2008 terrain pts: {len(p2008):,}   2021 ground pts: {len(p2021):,}")

    # core points: subsample the 2008 cloud onto a 4 m grid
    key = (np.floor(p2008[:, 0] / 4) * 100000 + np.floor(p2008[:, 1] / 4)).astype(np.int64)
    _, idx = np.unique(key, return_index=True)
    core = p2008[idx]
    print(f"core points (4 m): {len(core):,}")

    e2008 = py4dgeo.Epoch(p2008)
    e2021 = py4dgeo.Epoch(p2021)
    m3c2 = py4dgeo.M3C2(epochs=(e2008, e2021), corepoints=core,
                        normal_radii=(4.0,), cyl_radii=(2.0,), max_distance=10.0)
    dist, unc = m3c2.run()
    lod = unc["lodetection"]
    normals = m3c2.directions()                 # per-core normal vectors
    tilt = np.degrees(np.arccos(np.clip(np.abs(normals[:, 2]), 0, 1)))  # normal vs vertical
    print(f"mean normal tilt from vertical: {np.nanmean(tilt):.1f} deg "
          f"(M3C2 follows the surface, not z)")

    # --- slope + gridded-vertical diff at each core point (for the comparison) ---
    gy, gx = np.gradient(Z21, res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    vdiff = Z21 - _median_grid(xc[inb], yc[inb], zc[inb], res, X0, Y0, nx, ny)
    col = np.clip(((core[:, 0] - X0) / res).astype(int), 0, nx - 1)
    row = np.clip(((core[:, 1] - Y0) / res).astype(int), 0, ny - 1)
    sl_c = slope[row, col]; vd_c = vdiff[row, col]

    print("\n  stable-ground NOISE by slope band: M3C2 (normal) vs gridded-vertical")
    print(f"  {'slope band':>12} {'n':>7} {'M3C2 NMAD':>10} {'vert NMAD':>10} {'M3C2 LoD':>9}")
    for lo, hi in [(0, 3), (3, 8), (8, 15), (15, 30)]:
        s = np.isfinite(dist) & np.isfinite(vd_c) & (sl_c >= lo) & (sl_c < hi) \
            & (np.abs(dist) < 1.0)              # exclude gross change/outliers
        if s.sum() < 200:
            continue
        print(f"  {lo:>4}-{hi:<3} deg {s.sum():>7,} {_nmad(dist[s]):>10.4f} "
              f"{_nmad(vd_c[s]):>10.4f} {np.nanmedian(lod[s]):>9.4f}")


if __name__ == "__main__":
    main()
