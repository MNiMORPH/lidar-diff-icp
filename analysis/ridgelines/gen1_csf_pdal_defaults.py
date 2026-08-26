#!/usr/bin/env python3
"""Re-classify elba gen1 ground with PDAL's OWN CSF defaults, for comparison.

Our `classify_ground_csf` defaults depart from PDAL's on three knobs, all in the
"keep more, follow lower" direction: rigidness 1 (vs 3), threshold 1.5 (vs 0.5),
hdiff 0.5 (vs 0.3). This runs the stock configuration so the effect of that choice
can be measured rather than assumed -- on the ground surface AND on the steep-cell
coverage the loosening was meant to protect.

Tiled 2x2 with a 150 m halo (CSF is spatially local; ~2M points per tile keeps peak
RAM modest on a shared machine).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/gen1_csf_pdal_defaults.py
"""
import argparse, os, shutil
import numpy as np, laspy

from lidar_diff_icp.ground import classify_ground_csf

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="data/before/4342-29-64.laz")
ap.add_argument("--out", default="data/csf_cache/elba_gen1_pdaldefaults.las")
ap.add_argument("--nx", type=int, default=2)
ap.add_argument("--ny", type=int, default=2)
ap.add_argument("--overlap", type=float, default=150.0)
# PDAL filters.csf stock defaults, verified with `pdal --options filters.csf`
ap.add_argument("--rigidness", type=int, default=3)
ap.add_argument("--threshold", type=float, default=None)
ap.add_argument("--hdiff", type=float, default=None)
ap.add_argument("--resolution", type=float, default=None)
ap.add_argument("--outlier", action="store_true", help="add PDAL filters.outlier before CSF")
A = ap.parse_args()

os.makedirs("data/csf_cache", exist_ok=True)
tmpd = "data/derived/_csf_pdaldef_tmp"
os.makedirs(tmpd, exist_ok=True)

f = laspy.read(A.src)
x = np.asarray(f.x); y = np.asarray(f.y)
X0, X1 = x.min(), x.max(); Y0, Y1 = y.min(), y.max()
dx = (X1 - X0) / A.nx; dy = (Y1 - Y0) / A.ny
print(f"src N={len(x):,}  rigidness={A.rigidness} threshold={A.threshold} "
      f"hdiff={A.hdiff} resolution={A.resolution}", flush=True)

hdr = laspy.LasHeader(point_format=f.header.point_format.id, version=str(f.header.version))
hdr.scales = f.header.scales; hdr.offsets = f.header.offsets


def copy_all_dims(src, out_hdr, mask):
    out = laspy.LasData(out_hdr)
    out.x = np.asarray(src.x)[mask]; out.y = np.asarray(src.y)[mask]
    out.z = np.asarray(src.z)[mask]
    for d in src.point_format.dimension_names:
        if d in ("X", "Y", "Z"):
            continue
        setattr(out, d, np.asarray(getattr(src, d))[mask])
    return out


w = None; total = 0
try:
    for j in range(A.ny):
        for i in range(A.nx):
            cx0, cx1 = X0 + i*dx, X0 + (i+1)*dx
            cy0, cy1 = Y0 + j*dy, Y0 + (j+1)*dy
            hm = ((x >= cx0-A.overlap) & (x <= cx1+A.overlap) &
                  (y >= cy0-A.overlap) & (y <= cy1+A.overlap))
            tin = f"{tmpd}/tile_{i}_{j}.laz"
            copy_all_dims(f, hdr, hm).write(tin)
            print(f"tile[{i},{j}] halo N={hm.sum():,} -> CSF ...", flush=True)
            gp = classify_ground_csf(tin, resolution=A.resolution, rigidness=A.rigidness,
                                     threshold=A.threshold, hdiff=A.hdiff,
                                     outlier=A.outlier)
            g = laspy.read(gp)
            gx = np.asarray(g.x); gy = np.asarray(g.y)
            core = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            if i == 0:        core |= (gx < cx0)
            if i == A.nx-1:   core |= (gx >= cx1)
            if j == 0:        core |= (gy < cy0)
            if j == A.ny-1:   core |= (gy >= cy1)
            if w is None:
                ghdr = laspy.LasHeader(point_format=g.header.point_format.id,
                                       version=str(g.header.version))
                ghdr.scales = g.header.scales; ghdr.offsets = g.header.offsets
                w = laspy.open(A.out, mode="w", header=ghdr)
            w.write_points(copy_all_dims(g, ghdr, core).points)
            total += int(core.sum())
            print(f"tile[{i},{j}] ground-in-core {core.sum():,} (cum {total:,})", flush=True)
            os.remove(tin); shutil.rmtree(os.path.dirname(gp), ignore_errors=True)
finally:
    if w is not None:
        w.close()
    shutil.rmtree(tmpd, ignore_errors=True)
print(f"wrote {A.out}  N={total:,} ground points "
      f"({100*total/len(x):.1f}% of the cloud; our tuned run keeps ~88%)", flush=True)
