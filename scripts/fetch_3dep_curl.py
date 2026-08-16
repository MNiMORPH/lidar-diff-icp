#!/usr/bin/env python3
"""Robust 3DEP EPT patch fetch by direct curl of octree node tiles.

PDAL's readers.ept over S3 is unreliable in some environments (intermittent node
read failures), while plain HTTPS GETs are solid. This walks the EPT hierarchy,
downloads the node .laz tiles overlapping a bbox (depths up to --max-depth),
decodes with laspy, reprojects EPSG:3857 -> EPSG:26915, and writes one LAZ.
Preserves PointSourceId and GPS time.

Run with the PROJ fix if a conda base leaks it:
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal python scripts/fetch_3dep_curl.py \
      --base https://s3-us-west-2.amazonaws.com/usgs-lidar-public/MN_SEDriftless_2_2021 \
      --bounds 578014 4883738 579514 4885238 --max-depth 9 \
      --out data/after/3dep2021_subpatch.laz
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
    ap.add_argument("--base", required=True, help="EPT base URL (dir with ept.json)")
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"), help="EPSG:26915")
    ap.add_argument("--max-depth", type=int, default=9)
    ap.add_argument("--max-tiles", type=int, default=160)
    ap.add_argument("--workers", type=int, default=16, help="parallel download workers")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

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
    keep = sorted(set(keep))[: a.max_tiles]
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

    xs_, ys_, zs_, psid, cls, rn, nr = ([] for _ in range(7))
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
        xs_.append(np.asarray(f.x)); ys_.append(np.asarray(f.y)); zs_.append(np.asarray(f.z))
        psid.append(np.asarray(f.point_source_id)); cls.append(np.asarray(f.classification))
        rn.append(np.asarray(f.return_number)); nr.append(np.asarray(f.number_of_returns))
        if (i + 1) % 20 == 0:
            print(f"  read {i+1}/{len(keep)}", flush=True)
    X = np.concatenate(xs_); Y = np.concatenate(ys_); Z = np.concatenate(zs_)
    PS = np.concatenate(psid); CL = np.concatenate(cls)
    RN = np.concatenate(rn); NR = np.concatenate(nr)

    # reproject 3857 -> 26915 and clip to the exact bbox
    tr = Transformer.from_crs("EPSG:3857", "EPSG:26915", always_xy=True)
    E, N = tr.transform(X, Y)
    m = (E >= a.bounds[0]) & (E <= a.bounds[2]) & (N >= a.bounds[1]) & (N <= a.bounds[3])
    E, N, Z, PS, CL, RN, NR = E[m], N[m], Z[m], PS[m], CL[m], RN[m], NR[m]
    print(f"{E.size:,} points in bbox; PointSourceId uniques: {np.unique(PS).size}", flush=True)

    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [E.min(), N.min(), Z.min()]; hdr.scales = [0.01, 0.01, 0.01]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = E, N, Z
    las.point_source_id = PS.astype(np.uint16); las.classification = CL.astype(np.uint8)
    las.return_number = RN.astype(np.uint8); las.number_of_returns = NR.astype(np.uint8)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    las.write(a.out)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
