#!/usr/bin/env python3
"""Merge lidar tiles into one cloud, preserving EVERY point dimension.

The elbaext gen1 merge was built by copying a hand-picked field subset, which silently
dropped scan_angle_rank (-> incidence collapses to slope). This does it generically: read
each tile, crop to bounds+buffer, and concatenate ALL dimensions of the point format, so no
attribute is ever lost. Tiles must share a point format (checked).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/merge_gen1_tiles.py \
        OUT.laz  buffer_m  bminx bminy bmaxx bmaxy  tile1.laz tile2.laz ...
"""
import sys, numpy as np, laspy

OUT = sys.argv[1]; BUF = float(sys.argv[2])
bx0, by0, bx1, by1 = map(float, sys.argv[3:7]); TILES = sys.argv[7:]
X0, Y0, X1, Y1 = bx0 - BUF, by0 - BUF, bx1 + BUF, by1 + BUF

hdr = None; chunks = {}          # dim name -> list of cropped arrays
n_total = 0
for t in TILES:
    f = laspy.read(t)
    if hdr is None:
        hdr = laspy.LasHeader(point_format=f.header.point_format.id, version=str(f.header.version))
        hdr.scales = f.header.scales; hdr.offsets = f.header.offsets
        dims = list(f.point_format.dimension_names)
        chunks = {d: [] for d in dims}
    assert f.point_format.id == hdr.point_format.id, f"{t}: point format mismatch"
    x = np.asarray(f.x); y = np.asarray(f.y)
    m = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)
    for d in chunks:
        chunks[d].append(np.asarray(getattr(f, d))[m])
    n_total += int(m.sum())
    print(f"  {t}: kept {int(m.sum()):,} / {len(x):,}")

out = laspy.LasData(hdr)
for d, parts in chunks.items():
    setattr(out, d, np.concatenate(parts))
out.write(OUT)
sa = np.asarray(out.scan_angle_rank) if "scan_angle_rank" in chunks else np.asarray(out.scan_angle)
print(f"wrote {OUT}  ({n_total:,} pts, {len(chunks)} dims, scan_angle nonzero%={100*(sa!=0).mean():.0f})")
