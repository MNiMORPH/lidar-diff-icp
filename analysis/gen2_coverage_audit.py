#!/usr/bin/env python3
"""Per-cell gen2 return density, and what changes across a coverage boundary.

Written after whitewater's class-2 spread grid showed a hard vertical seam: west of easting
586362 the median cell holds ~115 class-2 returns, east of it ~22, while the class-2 FRACTION
only falls 0.34 -> 0.20. A 5x drop in count against a 1.7x drop in fraction means fewer
returns arriving, not a different classifier -- an acquisition boundary. This reports, per
side, the return density and the flight lines and acquisition times that produced it, so the
boundary can be attributed rather than guessed at.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/gen2_coverage_audit.py \
        --tile data/derived/whitewater --gen2 data/after/3dep_4358_fulltile.laz --split-x 586362
"""
import argparse, json
import numpy as np, laspy

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True)
ap.add_argument("--gen2", required=True)
ap.add_argument("--split-x", type=float, required=True,
                help="easting of the boundary to test. MINE unless stated: read off the "
                     "class-2 count grid, not fitted.")
ap.add_argument("--chunk", type=int, default=5_000_000)
A = ap.parse_args()

j = json.load(open(f"{A.tile}/corrections.json")); b = j["bounds"]; RES = float(j["res_m"])
tot = {"west": 0, "east": 0}; cls2 = {"west": 0, "east": 0}
psid = {"west": {}, "east": {}}
gps = {"west": [np.inf, -np.inf], "east": [np.inf, -np.inf]}
with laspy.open(A.gen2) as f:
    for pts in f.chunk_iterator(A.chunk):
        x = np.asarray(pts.x); y = np.asarray(pts.y)
        ing = (x >= b[0]) & (x < b[2]) & (y >= b[1]) & (y < b[3])
        if not ing.any():
            continue
        x = x[ing]; cl = np.asarray(pts.classification)[ing]
        ps = np.asarray(pts.point_source_id)[ing]; gt = np.asarray(pts.gps_time)[ing]
        for side, sel in (("west", x < A.split_x), ("east", x >= A.split_x)):
            if not sel.any():
                continue
            tot[side] += int(sel.sum()); cls2[side] += int((cl[sel] == 2).sum())
            u, c = np.unique(ps[sel], return_counts=True)
            for k, v in zip(u.tolist(), c.tolist()):
                psid[side][k] = psid[side].get(k, 0) + v
            gps[side][0] = min(gps[side][0], float(gt[sel].min()))
            gps[side][1] = max(gps[side][1], float(gt[sel].max()))

area = {"west": (A.split_x - b[0]) * (b[3] - b[1]),
        "east": (b[2] - A.split_x) * (b[3] - b[1])}
print(f"{A.tile}, split at easting {A.split_x:.0f}")
for side in ("west", "east"):
    n, n2, ar = tot[side], cls2[side], area[side]
    print(f"\n  {side.upper():5s}  area {ar/1e6:.2f} km2")
    print(f"    all returns  {n:12,}   {n/ar:7.2f} per m2")
    print(f"    class 2      {n2:12,}   {n2/ar:7.2f} per m2   fraction {n2/max(n,1):.3f}")
    print(f"    gps_time     {gps[side][0]:.1f} .. {gps[side][1]:.1f}  "
          f"(span {gps[side][1]-gps[side][0]:.1f} s)")
    tp = sum(psid[side].values())
    print(f"    flight lines ({len(psid[side])}):")
    for k in sorted(psid[side], key=lambda k: -psid[side][k]):
        print(f"        psid {k:6d}  {psid[side][k]:12,}  {100*psid[side][k]/tp:5.1f}%")
print(f"\n  DENSITY RATIO west/east: all {(tot['west']/area['west'])/(tot['east']/area['east']):.2f}x"
      f"   class2 {(cls2['west']/area['west'])/(cls2['east']/area['east']):.2f}x")
