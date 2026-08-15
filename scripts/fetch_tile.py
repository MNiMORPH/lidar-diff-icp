#!/usr/bin/env python3
"""Fetch the MN 2008 lidar tile covering a coordinate (or a named tile).

Examples
--------
    python scripts/fetch_tile.py --lon -92.004137 --lat 44.101944 --out data/before
    python scripts/fetch_tile.py --easting 579705.72 --northing 4883677.71
    python scripts/fetch_tile.py --tile 4342-29-64
"""
import argparse
from pathlib import Path

from lidar_diff_icp import MN_2008_CRS, tiles


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lon", type=float)
    p.add_argument("--lat", type=float)
    p.add_argument("--easting", type=float)
    p.add_argument("--northing", type=float)
    p.add_argument("--tile", help="explicit AAAA-BB-CC tile name")
    p.add_argument("--out", default="data/before", help="output directory")
    p.add_argument("--cache", default="data/tile_index_cache.json",
                   help="tile bbox index cache (JSON)")
    args = p.parse_args()

    if args.tile:
        name = args.tile
    else:
        if args.easting is not None and args.northing is not None:
            e, n = args.easting, args.northing
        elif args.lon is not None and args.lat is not None:
            from pyproj import Transformer
            e, n = Transformer.from_crs(
                "EPSG:4326", MN_2008_CRS, always_xy=True
            ).transform(args.lon, args.lat)
        else:
            p.error("give --tile, or --easting/--northing, or --lon/--lat")
        name = tiles.find_tile(e, n, cache=args.cache)
        print(f"coordinate falls in tile {name}")

    dest = tiles.download_tile(name, args.out)
    print(f"downloaded {dest} ({Path(dest).stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
