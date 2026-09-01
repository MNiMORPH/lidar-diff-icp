#!/usr/bin/env python3
"""Write a tile's slope.npy: surface slope in degrees from the gap-filled gen2 grid.

slope.npy was declared a base input of the derived-product chain, as if difference_dem
produced it. It does not. The only producer was analysis/slope_bias/fulldensity_regrid.py,
hardcoded to elba_fulldensity -- which is why carlton and cook have dod.npy, lod.npy,
z_after.npy and corrections.json but no slope, and why every downstream step there reports
its inputs absent.

The quantity needs nothing but z_after and the cell size, and the repo already agrees on the
definition: fulldensity_regrid.py:54 and gen1_save_angles_slope.py:52 compute it the same
way, and this recipe reproduces the shipped elba slope.npy BIT-IDENTICALLY (0 of 355,600
cells differ). Gaps in z_after are filled by nearest neighbour first, as both do, so the
gradient is not taken across a hole.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/make_slope.py \
        --tile elba_fulldensity --check
"""
import argparse, json, os

import numpy as np
from scipy.ndimage import distance_transform_edt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True, help="tile name under data/derived/")
ap.add_argument("--check", action="store_true",
                help="compare against the existing slope.npy and write nothing")
ap.add_argument("--out", default=None, help="default <tile>/slope.npy")
A = ap.parse_args()

D = f"data/derived/{A.tile}"
OUT = A.out or f"{D}/slope.npy"


def res_of(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p))
            return float(j.get("res") or j.get("res_m"))
    raise SystemExit(f"no grid meta for {tile}; cannot know the cell size")


RES = res_of(A.tile)
z = np.load(f"{D}/z_after.npy")
zf = z.copy()
gap = ~np.isfinite(zf)
if gap.any():
    zf = zf[tuple(distance_transform_edt(gap, return_distances=False, return_indices=True))]
gy, gx = np.gradient(zf, RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))

print(f"{A.tile}: {z.shape} at {RES} m; {int(gap.sum()):,} gap cells filled by nearest "
      f"neighbour before differencing")
print(f"  slope deg: median {np.median(slope):.2f}  p90 {np.percentile(slope, 90):.2f}  "
      f"max {slope.max():.2f}")

if A.check:
    if not os.path.exists(OUT):
        raise SystemExit(f"--check but {OUT} does not exist; nothing to compare against")
    old = np.load(OUT)
    print(f"\nagainst the existing {OUT}:")
    print(f"  bit-identical: {np.array_equal(slope, old, equal_nan=True)}")
    d = slope - old
    f = np.isfinite(d)
    print(f"  cells differing: {int((d[f] != 0).sum()):,} of {int(f.sum()):,}")
    print("\nnothing written (--check)")
else:
    np.save(OUT, slope)
    print(f"\nwrote {OUT}")
