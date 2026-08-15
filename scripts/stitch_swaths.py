#!/usr/bin/env python3
"""Stitch a tile's flight-line passes into one internally consistent frame.

Solves the free-network Nuth & Kaeaeb adjustment (coreg.align_swaths), applies
the per-swath shift, verifies that the swath-overlap offsets collapse, and
writes a single merged bare-earth DEM. The whole group's absolute datum (offset
from the 2021 3DEP) is left free -- tie it separately.

Example:
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
      python scripts/stitch_swaths.py data/before/4342-29-64.laz \
      --out data/derived/4342-29-64_stitched_dem.tif
"""
import argparse
from itertools import combinations
from pathlib import Path
import numpy as np

from lidar_diff_icp import io, coreg
from lidar_diff_icp.swathdiff import _median_grid

RES = 2.0
EXCLUDE = (5, 6, 9)


def _pair_offset(x, y, z, terr, psid, a, b, res):
    ma = terr & (psid == a); mb = terr & (psid == b)
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        return None
    nx = int(np.ceil((x1 - x0) / res)); ny = int(np.ceil((y1 - y0) / res))
    d = (_median_grid(x[ma], y[ma], z[ma], res, x0, y0, nx, ny)
         - _median_grid(x[mb], y[mb], z[mb], res, x0, y0, nx, ny))
    return float(np.nanmedian(d))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("laz")
    ap.add_argument("--res", type=float, default=RES)
    ap.add_argument("--ref", type=int, default=None,
                    help="swath id to pin as local reference (default: lowest id)")
    ap.add_argument("--out", help="output GeoTIFF for the merged DEM")
    args = ap.parse_args()

    pc = io.read_tile(args.laz)
    ref = args.ref if args.ref is not None else int(pc.swaths.min())
    corrections, edges, mis = coreg.align_swaths(pc, res=args.res, ref=ref)
    print(f"{Path(args.laz).name}: swaths {pc.swaths.tolist()}")
    print(f"per-swath correction into the common frame (reference swath {ref} pinned to 0):")
    print(f"  {'swath':>6} {'Dx_m':>8} {'Dy_m':>8} {'Dz_m':>8}")
    for s, (dx, dy, dz) in corrections.items():
        print(f"  {s:>6} {dx:>8.3f} {dy:>8.3f} {dz:>8.3f}")
    print(f"network misclosure (per-edge residual): "
          f"max |dx,dy,dz| = {np.abs(mis).max(axis=0).round(4).tolist()} m"
          f"  ({len(edges)} edges, {len(pc.swaths)} nodes)")

    # verify overlaps collapse: vertical offset before vs after
    terr = ~np.isin(pc.classification, EXCLUDE)
    xc, yc, zc = coreg.apply_alignment(pc, corrections)
    print("\nswath-overlap vertical offset (median dz), before -> after:")
    for a, b in combinations(pc.swaths.tolist(), 2):
        before = _pair_offset(pc.x, pc.y, pc.z, terr, pc.point_source_id, a, b, args.res)
        if before is None:
            continue
        after = _pair_offset(xc, yc, zc, terr, pc.point_source_id, a, b, args.res)
        print(f"  {a}-{b}:  {before:+.3f} m -> {after:+.3f} m")

    if args.out:
        _write_dem(xc, yc, zc, terr, args.res, args.out)


def _write_dem(x, y, z, terr, res, out):
    import rasterio
    from rasterio.transform import from_origin
    xg, yg, zg = x[terr], y[terr], z[terr]
    x0, y0 = xg.min(), yg.min()
    nx = int(np.ceil((xg.max() - x0) / res)); ny = int(np.ceil((yg.max() - y0) / res))
    dem = _median_grid(xg, yg, zg, res, x0, y0, nx, ny)
    dem = np.flipud(dem)  # rasterio row 0 = north
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out, "w", driver="GTiff", height=ny, width=nx, count=1,
                       dtype="float32", crs="EPSG:26915",
                       transform=from_origin(x0, y0 + ny * res, res, res),
                       nodata=np.nan) as dst:
        dst.write(dem.astype("float32"), 1)
    print(f"\nwrote merged internally-consistent DEM: {out}  ({nx}x{ny} @ {res} m)")


if __name__ == "__main__":
    main()
