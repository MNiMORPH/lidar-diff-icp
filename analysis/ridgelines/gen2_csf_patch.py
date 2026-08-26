#!/usr/bin/env python3
"""Re-classify gen2 with the SAME ground filter we use on gen1 (PDAL ELM + CSF), on one
400 m Elba patch, and measure how much of the gen1-gen2 gap that closes.

Motivation: the two epochs' "ground" surfaces are not the same kind of object. gen2's is
the median of the vendor's class-2 (which discards 52% of the near-ground window, mostly
from above); gen1's is the median of our CSF ground (which retains 97%). Part of the
offset we have been trying to correct may be that difference in RULE, not in the terrain.
The way to find out is to apply ONE rule to both epochs.

Patch chosen for cover spread: 1,716 reference divide cells, 1,154 open and 381 dense.
Read straight from the COPC per sub-tile so no multi-GB intermediate is ever materialised;
200 m cores with a 100 m halo keep each CSF near ~3M points (shared-laptop guardrail).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/gen2_csf_patch.py
"""
import argparse, json, os, shutil, subprocess, tempfile
import numpy as np, laspy

from lidar_diff_icp.ground import classify_ground_csf, find_pdal

ap = argparse.ArgumentParser()
ap.add_argument("--copc", default="data/after/3dep2021_fulldensity.copc.laz")
ap.add_argument("--out", default="data/csf_cache/elba_gen2_patch.las")
ap.add_argument("--x0", type=float, default=577592.8)
ap.add_argument("--y0", type=float, default=4884537.6)
ap.add_argument("--size", type=float, default=400.0)
ap.add_argument("--core", type=float, default=200.0)
ap.add_argument("--halo", type=float, default=100.0)
A = ap.parse_args()

pdal = find_pdal()
X1 = A.x0 + A.size; Y1 = A.y0 + A.size
tmp = tempfile.mkdtemp(prefix="gen2csf_")
w = None
total = 0
try:
    ncore = int(round(A.size / A.core))
    for j in range(ncore):
        for i in range(ncore):
            cx0 = A.x0 + i*A.core; cx1 = cx0 + A.core
            cy0 = A.y0 + j*A.core; cy1 = cy0 + A.core
            crop = os.path.join(tmp, f"t{i}{j}.las")
            pipe = {"pipeline": [
                {"type": "readers.copc", "filename": A.copc,
                 "bounds": f"([{cx0-A.halo},{cx1+A.halo}],[{cy0-A.halo},{cy1+A.halo}])"},
                {"type": "writers.las", "filename": crop}]}
            pj = os.path.join(tmp, "crop.json")
            json.dump(pipe, open(pj, "w"))
            subprocess.run([pdal, "pipeline", pj], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            n_in = laspy.open(crop).header.point_count
            print(f"tile[{i},{j}] read {n_in:,} pts -> CSF ...", flush=True)
            gp = classify_ground_csf(crop)          # SAME defaults as gen1: ELM + cloth
            g = laspy.read(gp)
            gx = np.asarray(g.x); gy = np.asarray(g.y)
            keep = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            if w is None:
                hdr = laspy.LasHeader(point_format=g.header.point_format.id,
                                      version=str(g.header.version))
                hdr.scales = g.header.scales; hdr.offsets = g.header.offsets
                w = laspy.open(A.out, mode="w", header=hdr)
            out = laspy.LasData(hdr)
            out.x = gx[keep]; out.y = gy[keep]; out.z = np.asarray(g.z)[keep]
            for d in g.point_format.dimension_names:
                if d in ("X", "Y", "Z"):
                    continue
                setattr(out, d, np.asarray(getattr(g, d))[keep])
            w.write_points(out.points)
            total += int(keep.sum())
            print(f"tile[{i},{j}] ground-in-core {keep.sum():,}  (cum {total:,})", flush=True)
            os.remove(crop)
            shutil.rmtree(os.path.dirname(gp), ignore_errors=True)
finally:
    if w is not None:
        w.close()
    shutil.rmtree(tmp, ignore_errors=True)
print(f"wrote {A.out}  N={total:,} CSF ground points over the {A.size:g} m patch")
