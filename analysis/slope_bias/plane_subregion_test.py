#!/usr/bin/env python3
"""SUBREGION PLANE-ESTIMATOR TEST (Task 4, decisive estimator check).

Does the +tan(S) bias come from the slope_normal MEDIAN "landing at a tilt-correlated spot
on a steep cell" (pipeline.py:514)? The tilt-unbiased estimators (plane/poly2) were only
ever tried on the DECIMATED cloud, where 6 gen2 pts/cell fall below the plane-fit threshold
and the code FALLS BACK to the median (pipeline.py:533,568) -- so they were never actually
exercised. On the FULL cloud (120 pts/cell) they are well-constrained.

Test on a forested-slope SUBREGION (plane/poly2 need stream=False -> must fit in RAM, so we
clip first). Same before/after, same corrections; only the ground ESTIMATOR changes:
  slope_normal (median)  vs  plane (per-cell LSQ plane, read at cell centre, tilt-unbiased).
If plane flattens the tan(S) trend -> estimator artifact. If it survives -> estimator out.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/plane_subregion_test.py
"""
import os, time, numpy as np, laspy
from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.canopy import ground_penetration

RES = 5.0
SUB = (577992.8, 4883237.6, 579492.8, 4884737.6)     # 1500x1500 m, grid-aligned
SCR = "/tmp/claude-1000/-home-awickert-projects-lidar-diff-icp/0804ccd6-abab-448a-9abe-ef4acc47d513/scratchpad"
G1, G2 = "data/before/4342-29-64.laz", "data/after/3dep2021_fulldensity.laz"
G1C, G2C = f"{SCR}/g1_sub.laz", f"{SCR}/g2_sub.laz"

def clip(src, dst, bounds, chunk=20_000_000):
    x0, y0, x1, y1 = bounds
    with laspy.open(src) as f:
        hdr = f.header
        oh = laspy.LasHeader(point_format=hdr.point_format, version=hdr.version)
        oh.scales = hdr.scales; oh.offsets = hdr.offsets
        w = laspy.open(dst, mode="w", header=oh)
        kept = 0
        for pts in f.chunk_iterator(chunk):
            x = np.asarray(pts.x); y = np.asarray(pts.y)
            m = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
            if m.any():
                w.write_points(pts[m]); kept += int(m.sum())
        w.close()
    return kept

for src, dst in ((G1, G1C), (G2, G2C)):
    if not os.path.exists(dst):
        t = time.time(); n = clip(src, dst, SUB)
        print(f"clipped {src.split('/')[-1]} -> {n:,} pts in {time.time()-t:.0f}s", flush=True)

COMMON = dict(res=RES, ground_source="csf", after_ground="class2", stream=False,
              robust_stable=True, csf_cache=f"{SCR}/g1_sub_csf.las")
runs = {}
for tag, est in (("slope_normal(median)", "slope_normal"), ("plane", "plane")):
    t = time.time()
    print(f"[{tag}] difference_dem ...", flush=True)
    runs[tag] = difference_dem(G1C, G2C, SUB, ground=est, **COMMON)
    r = runs[tag]; ex = np.isfinite(r["dod"])
    print(f"[{tag}] {time.time()-t:.0f}s  medDoD={np.nanmedian(r['dod'][ex])*1000:+.1f} mm  "
          f"sigma={r['stable_sigma']:.3f}", flush=True)

# slope + forest from the full-density subregion
base = runs["plane"]; nx, ny = base["nx"], base["ny"]; bnds = base["bounds"]
Z = base["z_after"].copy(); nm = ~np.isfinite(Z)
if nm.any():
    from scipy.ndimage import distance_transform_edt as edt
    Z = Z[tuple(edt(nm, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Z, RES); slope = np.degrees(np.arctan(np.hypot(gx, gy)))
pen = ground_penetration(G2C, bnds, RES, nx, ny); forest = pen < 0.25

BINS = [0, 10, 15, 20, 25, 30, 90]
print(f"\n=== median DoD (mm) vs slope on FORESTED cells ({int((forest&np.isfinite(base['dod'])).sum())} cells) ===")
print(f"{'slope':>8} {'n':>6} | " + " | ".join(f"{t:>18}" for t in runs))
for lo, hi in zip(BINS[:-1], BINS[1:]):
    row = f"{lo:>3}-{hi:<3} "
    n0 = None; cells = []
    for tag in runs:
        d = runs[tag]["dod"]; m = np.isfinite(d) & forest & (slope >= lo) & (slope < hi)
        cells.append((m.sum(), np.median(d[m])*1000 if m.any() else np.nan))
    row = f"{lo:>3}-{hi:<3} {cells[0][0]:>6} | " + " | ".join(f"{c[1]:>+18.1f}" for c in cells)
    print(row)
print("\nVERDICT: if 'plane' is ~flat vs slope while slope_normal(median) rises, the +tan(S) "
      "bias is the median's tilt-spot sensitivity under differential epoch sampling -> "
      "estimator fix. If plane ALSO rises, the estimator is exonerated.")
