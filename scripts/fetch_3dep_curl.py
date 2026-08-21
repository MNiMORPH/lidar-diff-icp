#!/usr/bin/env python3
"""Robust 3DEP EPT patch fetch by direct curl of octree node tiles.

PDAL's readers.ept over S3 is unreliable in some environments (intermittent node
read failures), while plain HTTPS GETs are solid. This walks the EPT hierarchy,
downloads the node .laz tiles overlapping a bbox (depths up to --max-depth),
decodes with laspy, reprojects EPSG:3857 -> EPSG:26915, and writes one LAS,
**streaming per tile** so peak RAM is ~one node tile rather than the whole cloud
(the dense 3DEP can be hundreds of millions of points). Preserves PointSourceId,
classification, and return numbers.

Give --base explicitly, or --auto to resolve the covering gen2 project from the
bbox (with a mandatory coverage check). Run with the PROJ fix if a conda base
leaks it:
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal python scripts/fetch_3dep_curl.py \
      --auto --bounds 578014 4883738 579514 4885238 --max-depth 9 \
      --out data/after/3dep2021_subpatch.laz
    # or pin the project:
    #   --base https://s3-us-west-2.amazonaws.com/usgs-lidar-public/MN_SEDriftless_2_2021
"""
import argparse, json, subprocess, urllib.request
from pathlib import Path
import numpy as np
import laspy
from pyproj import Transformer


def _get_json(url):
    return json.load(urllib.request.urlopen(url, timeout=60))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", help="EPT base URL (dir with ept.json). Omit with "
                    "--auto to resolve the covering gen2 3DEP project from --bounds.")
    ap.add_argument("--auto", action="store_true",
                    help="auto-resolve the gen2 3DEP EPT project covering --bounds "
                         "(most recent non-mosaic project) and REQUIRE it to fully "
                         "cover the bbox before downloading")
    ap.add_argument("--ept-cache", default=None,
                    help="cache file for the EPT project-boundary index (--auto)")
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"), help="EPSG:26915")
    ap.add_argument("--max-depth", type=int, default=9)
    ap.add_argument("--max-tiles", type=int, default=None,
                    help="cap on node tiles fetched (default: NO cap -- fetch every "
                         "overlapping tile for complete coverage). A cap that "
                         "truncates keeps the lexically-first tiles and prints a "
                         "loud incomplete-coverage warning (this silently caused a "
                         "half-tile coverage gap before).")
    ap.add_argument("--workers", type=int, default=16, help="parallel download workers")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not a.base and not a.auto:
        ap.error("give --base <EPT url>, or --auto to resolve it from --bounds")
    if a.auto:
        from lidar_diff_icp import threedep
        # bbox (26915) corners -> lon/lat; resolve the covering gen2 project and
        # require full coverage before we download anything.
        to_ll = Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
        pts = [to_ll.transform(cx, cy) for cx in a.bounds[0::2] for cy in a.bounds[1::2]]
        lons = [p[0] for p in pts]; lats = [p[1] for p in pts]
        bbox_ll = (min(lons), min(lats), max(lons), max(lats))
        clon, clat = (min(lons) + max(lons)) / 2, (min(lats) + max(lats)) / 2
        ref = threedep.resolve_reference(clon, clat, bbox_ll, cache=a.ept_cache)
        a.base = ref["url"].rsplit("/ept.json", 1)[0]
        print(f"auto-resolved gen2 reference: {ref['name']} (year {ref['latest']}); "
              f"boundary fully covers the tile bbox", flush=True)

    ept = _get_json(f"{a.base}/ept.json")
    b = ept["bounds"]                                   # [x0,y0,z0,x1,y1,z1] in 3857
    side = b[3] - b[0]
    # bbox (26915) -> 3857 envelope
    t = Transformer.from_crs("EPSG:26915", "EPSG:3857", always_xy=True)
    xs, ys = zip(*[t.transform(cx, cy) for cx in a.bounds[0::2] for cy in a.bounds[1::2]])
    bx0, bx1, by0, by1 = min(xs), max(xs), min(ys), max(ys)

    def nb(d, x, y):                                    # node xy extent in 3857
        s = side / (2 ** d)
        return b[0] + x * s, b[1] + y * s, b[0] + (x + 1) * s, b[1] + (y + 1) * s
    def overlap(d, x, y):
        n = nb(d, x, y)
        return not (n[2] < bx0 or n[0] > bx1 or n[3] < by0 or n[1] > by1)

    # collect overlapping existing nodes, fetching sub-hierarchies as needed
    keep = []
    seen_hier = set()
    frontier = [_get_json(f"{a.base}/ept-hierarchy/0-0-0-0.json")]
    while frontier:
        H = frontier.pop()
        for key, cnt in H.items():
            d, x, y, z = map(int, key.split("-"))
            if d > a.max_depth or not overlap(d, x, y):
                continue
            keep.append(key)
            # descend into a stored sub-hierarchy at the file boundary
            if cnt == -1 and key not in seen_hier and d <= a.max_depth:
                seen_hier.add(key)
                try:
                    frontier.append(_get_json(f"{a.base}/ept-hierarchy/{key}.json"))
                except Exception as e:
                    print(f"  sub-hierarchy {key} failed: {e}")
    keep = sorted(set(keep))
    n_overlap = len(keep)
    if a.max_tiles is not None and n_overlap > a.max_tiles:
        import sys
        print(f"WARNING: {n_overlap} node tiles overlap the bbox but --max-tiles="
              f"{a.max_tiles} keeps only the lexically-first {a.max_tiles} -> "
              f"{n_overlap - a.max_tiles} DROPPED, INCOMPLETE COVERAGE. Remove "
              f"--max-tiles (default) or raise it to cover the whole area.",
              file=sys.stderr, flush=True)
        keep = keep[: a.max_tiles]
    print(f"{len(keep)} overlapping node tiles (depth<= {a.max_depth})", flush=True)

    tmp = Path(a.out).parent / "_ept_tiles"; tmp.mkdir(parents=True, exist_ok=True)

    def dl(key, dst):
        for _ in range(4):                           # resilient to transient DNS/network
            r = subprocess.run(
                ["curl", "-sS", "--retry", "5", "--retry-delay", "2",
                 "--retry-all-errors", "--retry-connrefused", "--max-time", "120",
                 "-o", str(dst), f"{a.base}/ept-data/{key}.laz"])
            if r.returncode == 0 and dst.exists() and dst.stat().st_size > 0:
                return
            dst.unlink(missing_ok=True)
        raise RuntimeError(f"failed to download {key}")

    from concurrent.futures import ThreadPoolExecutor
    need = [k for k in keep if not (tmp / f"{k}.laz").exists()]
    if need:
        print(f"downloading {len(need)} tiles ({a.workers} workers)...", flush=True)
        with ThreadPoolExecutor(max_workers=a.workers) as ex:
            list(ex.map(lambda k: dl(k, tmp / f"{k}.laz"), need))
        print("  downloads complete", flush=True)

    # Stream: read each tile, reproject 3857 -> 26915, clip to the bbox, and write
    # incrementally. Peak RAM is one tile, not the whole cloud -- the dense 3DEP
    # can be hundreds of millions of points (461M on the MN River valley tile), so
    # accumulating everything then concatenating OOMs. Fixed header offsets (the
    # bbox min, which every clipped point is >=) avoid needing all data up front.
    tr = Transformer.from_crs("EPSG:3857", "EPSG:26915", always_xy=True)
    from laspy.header import GpsTimeType
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.global_encoding.gps_time_type = GpsTimeType.STANDARD   # 3DEP gps_time = Adjusted Standard GPS Time
    hdr.offsets = [a.bounds[0], a.bounds[1], 0.0]; hdr.scales = [0.01, 0.01, 0.01]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    total = 0; psids = set()
    with laspy.open(a.out, mode="w", header=hdr) as writer:
        for i, key in enumerate(keep):
            dst = tmp / f"{key}.laz"
            if not dst.exists():
                dl(key, dst)
            for _ in range(2):                           # heal truncated/partial tiles
                try:
                    f = laspy.read(str(dst)); break
                except Exception as e:
                    print(f"  re-downloading {key} (bad read: {e})", flush=True)
                    dst.unlink(missing_ok=True); dl(key, dst)
            else:
                raise RuntimeError(f"{key} unreadable after re-download")
            E, N = tr.transform(np.asarray(f.x), np.asarray(f.y))
            m = (E >= a.bounds[0]) & (E <= a.bounds[2]) & (N >= a.bounds[1]) & (N <= a.bounds[3])
            if not m.any():
                continue
            ch = laspy.LasData(hdr)
            ch.x = E[m]; ch.y = N[m]; ch.z = np.asarray(f.z)[m]
            ps = np.asarray(f.point_source_id)[m]
            ch.point_source_id = ps.astype(np.uint16)
            ch.classification = np.asarray(f.classification)[m].astype(np.uint8)
            ch.return_number = np.asarray(f.return_number)[m].astype(np.uint8)
            ch.number_of_returns = np.asarray(f.number_of_returns)[m].astype(np.uint8)
            # PRESERVE acquisition time + scan geometry + intensity (a merge that drops
            # these is why the flight date and scan angle were lost -- read from source).
            ch.gps_time = np.asarray(f.gps_time)[m]
            ch.scan_angle_rank = np.asarray(f.scan_angle_rank)[m].astype(np.int8)
            ch.intensity = np.asarray(f.intensity)[m].astype(np.uint16)
            writer.write_points(ch.points)
            total += int(m.sum()); psids.update(np.unique(ps).tolist())
            if (i + 1) % 20 == 0:
                print(f"  streamed {i+1}/{len(keep)} tiles, {total:,} pts in bbox", flush=True)
    print(f"{total:,} points in bbox; PointSourceId uniques: {len(psids)}", flush=True)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
