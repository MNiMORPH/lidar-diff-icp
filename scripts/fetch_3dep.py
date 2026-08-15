#!/usr/bin/env python3
"""Fetch a patch of USGS 3DEP lidar from the public EPT archive (AWS), by bbox.

Reads only the octree nodes overlapping the bounds (bandwidth-efficient),
decimates via --resolution, and reprojects EPSG:3857 (EPT storage) -> EPSG:26915
so the patch is in true metres and matches the 2008 tile. Preserves PointSourceId
and GpsTime (LAS 1.4) so flightline-consistency analysis is possible.

Requires PDAL -> run in the conda `lidar-icp` env:
    source ~/anaconda3/etc/profile.d/conda.sh && conda activate lidar-icp
    python scripts/fetch_3dep.py \
      --url https://s3-us-west-2.amazonaws.com/usgs-lidar-public/MN_SEDriftless_2_2021/ept.json \
      --bounds 577492.8 4882737.6 580035.0 4886238.3 --resolution 1.0 \
      --out data/after/3dep2021_4342-29-64.laz
"""
import argparse, json, os
os.environ.update(GDAL_HTTP_MAX_RETRY="10", GDAL_HTTP_RETRY_DELAY="2",
                  GDAL_HTTP_TIMEOUT="60", GDAL_HTTP_CONNECTTIMEOUT="15")
import pdal
from pyproj import Transformer


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", required=True, help="EPT ept.json URL")
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"), help="EPSG:26915")
    ap.add_argument("--resolution", type=float, default=1.0,
                    help="min point spacing to read (m); larger = fewer points")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    minx, miny, maxx, maxy = a.bounds

    # EPT bounds must be in the EPT SRS (EPSG:3857); take the envelope of corners
    t = Transformer.from_crs("EPSG:26915", "EPSG:3857", always_xy=True)
    xs, ys = [], []
    for cx in (minx, maxx):
        for cy in (miny, maxy):
            X, Y = t.transform(cx, cy); xs.append(X); ys.append(Y)
    bounds = f"([{min(xs)},{max(xs)}],[{min(ys)},{max(ys)}])"

    pipe = [
        {"type": "readers.ept", "filename": a.url, "bounds": bounds,
         "resolution": a.resolution},
        {"type": "filters.reprojection", "in_srs": "EPSG:3857",
         "out_srs": "EPSG:26915"},
        {"type": "writers.las", "filename": a.out, "minor_version": 4,
         "dataformat_id": 6, "compression": True},
    ]
    print("bounds (3857):", bounds, flush=True)
    n = None
    for attempt in range(1, 4):                     # retry flaky S3 node reads
        try:
            n = pdal.Pipeline(json.dumps(pipe)).execute()
            break
        except RuntimeError as exc:
            print(f"attempt {attempt} failed: {exc}", flush=True)
    if n is None:
        raise SystemExit("EPT read failed after retries")
    print(f"wrote {a.out}: {n:,} points", flush=True)


if __name__ == "__main__":
    main()
