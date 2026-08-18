#!/usr/bin/env python3
"""Fetch the MN 2008 lidar tile covering a coordinate (or a named tile).

The county is picked from the coordinate automatically (FCC Census Area lookup,
verified against the live MnGeo county listing); override with --county for a
point on a county line, or when giving --easting/--northing (which can't be
geocoded without a lon/lat).

Examples
--------
    python scripts/fetch_tile.py --lon -92.004137 --lat 44.101944 --out data/before
    python scripts/fetch_tile.py --lon -90.17978 --lat 48.04963   # NE MN (Cook Co.)
    python scripts/fetch_tile.py --easting 579705.72 --northing 4883677.71 --county winona
    python scripts/fetch_tile.py --tile 4342-29-64 --county winona
"""
import argparse
from pathlib import Path

from lidar_diff_icp import MN_GEN1_CRS, tiles


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lon", type=float)
    p.add_argument("--lat", type=float)
    p.add_argument("--easting", type=float)
    p.add_argument("--northing", type=float)
    p.add_argument("--tile", help="explicit AAAA-BB-CC tile name")
    p.add_argument("--county", help="MnGeo county dir for the download path "
                   "(default: resolved from the coordinate); needed only with --tile)")
    p.add_argument("--out", default="data/before", help="output directory")
    p.add_argument("--cache", default=tiles.DEFAULT_TILE_INDEX_CACHE,
                   help="statewide tile-centroid index cache (CSV)")
    args = p.parse_args()

    from pyproj import Transformer
    county = args.county
    if args.tile:
        name = args.tile
        if county is None:
            p.error("give --county with an explicit --tile")
    else:
        if args.lon is not None and args.lat is not None:
            lon, lat = args.lon, args.lat
            e, n = Transformer.from_crs("EPSG:4326", MN_GEN1_CRS,
                                        always_xy=True).transform(lon, lat)
        elif args.easting is not None and args.northing is not None:
            e, n = args.easting, args.northing
            lon, lat = Transformer.from_crs(MN_GEN1_CRS, "EPSG:4326",
                                            always_xy=True).transform(e, n)
        else:
            p.error("give --tile (+--county), or --easting/--northing, or --lon/--lat")
        name = tiles.find_tile(e, n, cache=args.cache)     # statewide centroid index
        if county is None:
            county = tiles.county_for_lonlat(lon, lat)     # county only for the download path
        print(f"coordinate -> tile {name} (county {county})")

    dest = tiles.download_tile(name, args.out, county=county)
    print(f"downloaded {dest} ({Path(dest).stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
