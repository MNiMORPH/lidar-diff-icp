#!/usr/bin/env python3
"""Write a tile's penetration.npy -- the producer that never existed.

`penetration.npy` had no producer in the source: `canopy.ground_penetration()` computed it
and `run_all_sites.py` called that function but persisted only the derived leaf-on flag, so
the layer existed for one tile as an interactive by-product. Reproducing the shipped elba
file from that function confirmed the origin exactly -- 354,923 cells finite in both, ZERO
differing -- with one deviation: 677 cells that have NO gen2 returns hold 0.0 in the shipped
file where the function returns NaN. Under the cut every consumer uses (forest = pen < 0.25)
those read as maximally closed canopy. Writing the layer from here fixes that by
construction, because NaN is what the function returns.

--check compares against the existing file and writes NOTHING, which is how the origin above
was established; use it before overwriting a layer that 34 tracked files consume.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/make_penetration.py \
        --tile elba_fulldensity --after data/after/3dep2021_fulldensity.laz --check

The AFTER cloud must be the FULL-RETURN gen2 cloud, not a class-2 extract: penetration is
ground returns over all non-noise returns, so a class-2 file gives 1.0 everywhere. The script
refuses if the cloud it is given contains only ground.
"""
import argparse, json, os
import numpy as np

from lidar_diff_icp.canopy import ground_penetration

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True, help="tile name under data/derived/")
ap.add_argument("--after", required=True,
                help="FULL-RETURN gen2 cloud. Not derivable from the tile directory, so it "
                     "is required rather than guessed.")
ap.add_argument("--check", action="store_true",
                help="compare against the existing penetration.npy and write nothing")
ap.add_argument("--out", default=None, help="output path (default <tile>/penetration.npy)")
A = ap.parse_args()

D = f"data/derived/{A.tile}"
OUT = A.out or f"{D}/penetration.npy"


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]
            r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], int(round((b[2]-b[0])/r)), int(round((b[3]-b[1])/r)), r
    raise SystemExit(f"no grid meta for {tile}")


X0, Y0, NX, NY, RES = grid(A.tile)
z = np.load(f"{D}/z_after.npy")
if z.shape != (NY, NX):
    raise SystemExit(f"grid meta says ({NY}, {NX}) but z_after.npy is {z.shape}")
print(f"{A.tile}: grid ({NY}, {NX}) at {RES} m from ({X0}, {Y0})")

frac = ground_penetration(A.after, (X0, Y0, X0 + NX*RES, Y0 + NY*RES), RES, NX, NY)
fin = np.isfinite(frac)
if fin.any() and np.nanmin(frac) == 1.0:
    raise SystemExit(f"{A.after} gives penetration 1.0 everywhere, so it holds only ground "
                     f"returns. Penetration needs the FULL-RETURN cloud.")
print(f"  measured {int(fin.sum()):,} cells; {int((~fin).sum()):,} have NO returns and are "
      f"NaN (not 0.0 -- absent is not a measurement)")
if fin.any():
    print(f"  median {np.median(frac[fin]):.4f}   min {frac[fin].min():.4f}   "
          f"max {frac[fin].max():.4f}")

if A.check:
    if not os.path.exists(OUT):
        raise SystemExit(f"--check but {OUT} does not exist; nothing to compare against")
    old = np.load(OUT)
    both = fin & np.isfinite(old)
    print(f"\nagainst the existing {OUT}:")
    print(f"  finite in BOTH               {int(both.sum()):,}")
    print(f"  of those, cells that DIFFER  {int((frac[both] != old[both]).sum()):,}")
    nan_now = ~fin & np.isfinite(old)
    print(f"  now NaN, previously finite   {int(nan_now.sum()):,}"
          + (f"  (previously {np.unique(old[nan_now])})" if nan_now.any() else ""))
    print(f"  now finite, previously NaN   {int((fin & ~np.isfinite(old)).sum()):,}")
    print("\nnothing written (--check)")
else:
    np.save(OUT, frac)
    print(f"\nwrote {OUT}")
