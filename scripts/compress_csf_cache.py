#!/usr/bin/env python3
"""Compress the CSF ground cache in place: every data/csf_cache/<name>.las becomes <name>.laz.

Andy, 2026-11-13 ("Run the entire process automatedly"), after the disk reached 5 GB free: the
cache was 10.1 GB of uncompressed LAS and one tile compressed 11x (216 -> 19 MB in 2 s). Every
reader goes through laspy, which reads LAZ, and the path literals were renamed in the same commit
series, so the names are honest.

Per file: read the LAS, write the LAZ, read the LAZ back, and compare point count, point format and
the scaled integer X Y Z, classification, point_source_id, gps_time, return_number,
number_of_returns, intensity and scan angle arrays element by element. Only when every array is
identical is the LAS deleted. A file whose round trip fails is left as it was and reported.
Idempotent: an existing verified .laz is skipped. Prints one line per file and a total.
"""
import glob
import os
import sys
import time

import laspy
import numpy as np

CACHE = sys.argv[1] if len(sys.argv) > 1 else "data/csf_cache"
FIELDS = ["X", "Y", "Z", "classification", "point_source_id", "gps_time", "return_number", "number_of_returns", "intensity"]


def same(a, b):
    if len(a) != len(b) or a.header.point_format.id != b.header.point_format.id: return False
    for f in FIELDS + ["scan_angle", "scan_angle_rank"]:
        try: xa, xb = np.asarray(a[f]), np.asarray(b[f])
        except Exception: continue   # dimension absent in this point format
        if not np.array_equal(xa, xb): return False
    return True


before = after = 0; failed = []
for src in sorted(glob.glob(os.path.join(CACHE, "*.las"))):
    dst = src[:-4] + ".laz"; t = time.time(); s0 = os.path.getsize(src)
    a = laspy.read(src)
    if not os.path.exists(dst): a.write(dst)
    b = laspy.read(dst)
    ok = same(a, b); del b
    if ok:
        os.remove(src); before += s0; after += os.path.getsize(dst)
        print(f"{os.path.basename(src):36s} {s0/1e6:7.0f} MB -> {os.path.getsize(dst)/1e6:5.0f} MB  {len(a):>11,d} pts verified, LAS removed  ({time.time()-t:.0f} s)", flush=True)
    else:
        failed.append(src); print(f"{os.path.basename(src):36s} ROUND TRIP DIFFERS; LAS kept, LAZ removed", flush=True); os.remove(dst)
    del a
print(f"TOTAL {before/1e9:.2f} GB -> {after/1e9:.2f} GB; failed: {failed if failed else 'none'}")
