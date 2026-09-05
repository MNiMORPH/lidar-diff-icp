#!/usr/bin/env python3
"""Did we download ALL the 3DEP points over each tile, or only some?

Whole-file density cannot answer this: whitewater's truncated file averaged 11.39 returns
per m2 across a tile that was 15.45 west and 5.52 east, and would have passed any
tile-average test. This asks the source instead. The EPT hierarchy JSON carries a POINT
COUNT PER NODE, so the number of points available over a bbox can be summed WITHOUT
downloading any of them, and compared with what our file actually holds.

Reports per site: points available in the bbox by octree depth, points in our file, and the
ratio. A ratio near 1 means we have it all; well below 1 means the fetch was truncated --
by --max-tiles (which keeps the lexically-first tiles and drops the rest) or by a
a depth cap shallower than the data.

    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
        ./lidar-icp/bin/python analysis/ept_coverage_check.py
"""
import argparse, json, sys, urllib.request
import laspy
from pyproj import Transformer

sys.path.insert(0, "scripts")
from run_all_sites import header_bounds
from lidar_diff_icp import threedep
from lidar_diff_icp.sites import SITES, site
from lidar_diff_icp import completeness

ap = argparse.ArgumentParser()
ap.add_argument("--only", nargs="*", help="site names; default all")
ap.add_argument("--cache", default="data/ept_index_cache.json")
ap.add_argument("--write", action="store_true",
                help="record the measurement in each tile's data_completeness.json, so "
                     "the gate has something to read. Without this the run only prints.")
A = ap.parse_args()


def _get(url):
    return json.load(urllib.request.urlopen(url, timeout=60))


def available(bounds):
    to_ll = Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
    pts = [to_ll.transform(cx, cy) for cx in bounds[0::2] for cy in bounds[1::2]]
    lons = [p[0] for p in pts]; lats = [p[1] for p in pts]
    ref = threedep.resolve_reference((min(lons)+max(lons))/2, (min(lats)+max(lats))/2,
                                     (min(lons), min(lats), max(lons), max(lats)),
                                     cache=A.cache)
    base = ref["url"].rsplit("/ept.json", 1)[0]
    ept = _get(f"{base}/ept.json"); b = ept["bounds"]; side = b[3] - b[0]
    t = Transformer.from_crs("EPSG:26915", "EPSG:3857", always_xy=True)
    xs, ys = zip(*[t.transform(cx, cy) for cx in bounds[0::2] for cy in bounds[1::2]])
    bx0, bx1, by0, by1 = min(xs), max(xs), min(ys), max(ys)

    def ov(d, x, y):
        s = side / (2 ** d)
        n = (b[0]+x*s, b[1]+y*s, b[0]+(x+1)*s, b[1]+(y+1)*s)
        return not (n[2] < bx0 or n[0] > bx1 or n[3] < by0 or n[1] > by1)

    per_depth = {}; nodes = 0; seen = set(); inb = 0.0

    def frac_in(d, x, y):
        s_ = side / (2 ** d)
        n = (b[0]+x*s_, b[1]+y*s_, b[0]+(x+1)*s_, b[1]+(y+1)*s_)
        a = (max(0.0, min(n[2], bx1)-max(n[0], bx0))
             * max(0.0, min(n[3], by1)-max(n[1], by0)))
        return a / (s_*s_)
    frontier = [_get(f"{base}/ept-hierarchy/0-0-0-0.json")]
    while frontier:
        H = frontier.pop()
        for key, cnt in H.items():
            d, x, y, z = map(int, key.split("-"))
            if not ov(d, x, y):
                continue
            nodes += 1
            if cnt == -1:
                if key not in seen:
                    seen.add(key)
                    try:
                        frontier.append(_get(f"{base}/ept-hierarchy/{key}.json"))
                    except Exception as e:
                        print(f"    sub-hierarchy {key} failed: {e}")
            else:
                per_depth[d] = per_depth.get(d, 0) + cnt
                inb += cnt * frac_in(d, x, y)
    return ref["name"], nodes, per_depth, inb


names = A.only or list(SITES)
print("ratio = our file / the AREA-WEIGHTED share of the node set inside the bbox.")
print("Nodes clip the bbox, so raw 'available' overstates what a complete fetch would hold;")
print("this cost two sites being read as short when they were complete.")
print(f"\n{'site':13s} {'project':26s} {'node pts':>13s} {'in-bbox est':>13s} "
      f"{'our file':>13s} {'ratio':>7s}  deepest")
for nm in names:
    S = site(nm)
    before, after, bnds = S.gen1, S.gen2, S.bounds
    if bnds is None:
        bnds = header_bounds(before, 5.0)
    try:
        proj, nodes, per_depth, inb = available(bnds)
    except Exception as e:
        print(f"{nm:13s} FAILED: {type(e).__name__}: {str(e)[:70]}"); continue
    avail = sum(per_depth.values())
    with laspy.open(after) as f:
        have = f.header.point_count
    deep = max(per_depth) if per_depth else -1
    print(f"{nm:13s} {proj[:26]:26s} {avail:13,} {inb:13,.0f} {have:13,} "
          f"{have/max(inb,1):7.2f}  d{deep}")
    if A.write:
        # gen2 only: this asks the 3DEP EPT source, which has no gen1 counterpart. gen1's
        # completeness is a different question (a delivered county tile, not a fetch) and
        # is deliberately left unrecorded rather than guessed at.
        completeness.write(
            S.tile_dir, epoch="gen2", cloud=after, points_in_file=have,
            points_available=inb,
            measured_by=f"analysis/ept_coverage_check.py against {proj}",
            note="area-weighted share of the EPT node set inside the bbox; nodes clip the "
                 "box, so a raw node-count sum overstates a complete fetch")
        print(f"              -> wrote {completeness.record_path(S.tile_dir)}")
