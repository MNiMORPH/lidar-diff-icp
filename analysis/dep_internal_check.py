#!/usr/bin/env python3
"""Check a survey's own internal flight-line consistency (before using it as a
reference). Nuth & Kaeaeb between overlapping flight lines on GROUND returns.

3DEP classifies most points as 1 (unclassified) rather than into veg classes, so
the terrain surface is taken from ground class 2 (not the single-return-terrain
filter used for the 2008 data). Reuses the validated core (coreg, variogram).

Example:
    python analysis/dep_internal_check.py data/after/3dep2021_subpatch.laz --res 3
"""
import argparse
from itertools import combinations
import numpy as np
import laspy

from lidar_diff_icp.swathdiff import _median_grid
from lidar_diff_icp import coreg, variogram as vg


def _nmad(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("laz")
    ap.add_argument("--res", type=float, default=3.0)
    ap.add_argument("--ground-class", type=int, default=2)
    ap.add_argument("--min-overlap-cells", type=int, default=2000)
    args = ap.parse_args()

    f = laspy.read(args.laz)
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    ps = np.asarray(f.point_source_id); cls = np.asarray(f.classification)
    ground = cls == args.ground_class
    lines = [int(v) for v in np.unique(ps[ground]) if (ground & (ps == v)).sum() > 5000]
    print(f"{args.laz}: flight lines with ground {lines}  (res {args.res} m)")
    print(f"  {'pair':>13} {'overlap_cells':>13} {'dz0':>7} {'nmad0':>7} "
          f"{'dx':>7} {'dy':>7} {'dz':>7} {'nmad1':>7} {'range':>6}")

    for a, b in combinations(lines, 2):
        ma = ground & (ps == a); mb = ground & (ps == b)
        x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
        y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
        if x1 <= x0 or y1 <= y0:
            continue
        nx = int(np.ceil((x1 - x0) / args.res)); ny = int(np.ceil((y1 - y0) / args.res))
        zr = _median_grid(x[ma], y[ma], z[ma], args.res, x0, y0, nx, ny)
        zs = _median_grid(x[mb], y[mb], z[mb], args.res, x0, y0, nx, ny)
        d0 = zr - zs; m0 = np.isfinite(d0)
        if m0.sum() < args.min_overlap_cells:
            continue
        c = coreg.nuth_kaab(zr, zs, args.res)
        r = zr - (coreg._shift_grid(zs, c.dx, c.dy, args.res) + c.dz)
        m = np.isfinite(r)
        gy, gx = np.mgrid[0:ny, 0:nx]
        xs = x0 + (gx[m] + 0.5) * args.res; ys = y0 + (gy[m] + 0.5) * args.res
        mod = vg.fit_spherical(*vg.empirical_variogram(xs, ys, r[m], 400))
        print(f"  {a:>6}-{b:<6} {m0.sum():>13,} {np.nanmedian(d0):>7.3f} "
              f"{_nmad(d0[m0]):>7.3f} {c.dx:>7.3f} {c.dy:>7.3f} {c.dz:>7.3f} "
              f"{_nmad(r[m]):>7.3f} {mod.range_:>6.0f}")
    print("  (dx,dy,dz ~ 0 => no detectable inter-flightline offset; nmad1 is the")
    print("   surface roughness on the overlap, not a registration error.)")


if __name__ == "__main__":
    main()
