#!/usr/bin/env python3
"""Extract the gen2 (3DEP) ASPRS class-2 GROUND to its own file for elbaext.

Why: with stream=True the after cloud is never in memory, so
references.flat_hard_cells cannot find the pavement/pad reference surfaces, and the
pipeline (correctly) REFUSES to fall back to the deactivated parabola tie
(pipeline.py PARABOLA guardrail). The fix the error message itself prescribes is to
run NON-streaming with the gen2 ground in memory. The full 415M-point 3DEP cloud is
too large to load, but the class-2 GROUND subset (132M pts, 31.8%) as its own file is
tractable and is exactly what read_after_ground(mode="class2") consumes.

Streams the source in chunks (O(one chunk) RAM), keeps Classification==2 in the
elbaext buffer, preserves point_source_id / gps_time / return numbers / class.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/elbaext_extract_gen2_ground.py
"""
import numpy as np, laspy

SRC = "data/after/elbaext_3dep_fulldensity.laz"
OUT = "data/after/elbaext_3dep_fd_class2.laz"

with laspy.open(SRC) as fh:
    hdr_in = fh.header
    out_hdr = laspy.LasHeader(point_format=hdr_in.point_format.id,
                              version=str(hdr_in.version))
    out_hdr.scales = hdr_in.scales; out_hdr.offsets = hdr_in.offsets
    total = 0
    with laspy.open(OUT, mode="w", header=out_hdr) as w:
        for pts in fh.chunk_iterator(20_000_000):
            m = np.asarray(pts.classification) == 2
            if not m.any():
                continue
            w.write_points(pts[m])
            total += int(m.sum())
print(f"wrote {OUT}  N={total:,} class-2 ground points")
