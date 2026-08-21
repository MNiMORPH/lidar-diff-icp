#!/usr/bin/env python3
"""STEP 1 (v2) — trace ridgelines with the Scherler & Schwanghart (2020) divide network,
via Andy's verified reimplementation (rivernetworkx.dreich.drainage_divides).

Divides are basin BOUNDARIES of the channel network (the dual of the drainage net), so
sub-channelization-threshold farm-furrow micro-topography cannot create spurious ridges --
the failure mode of the naive inverted-flow-accumulation approach (which traced tillage
furrows in the fields). QC each threshold by furrow contamination (open+flat field cells).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/trace_ridgelines.py
"""
import sys, os, numpy as np
sys.path.insert(0, "/home/awickert/dataanalysis/r.fluvial")
from rivernetworkx import dreich as D
from scipy.ndimage import uniform_filter, distance_transform_edt, gaussian_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lidar_diff_icp.viz import hillshade

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load("data/derived/elba_fulldensity/z_after.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
pen = np.load("data/derived/elba_fulldensity/penetration.npy")
ny, nx = z.shape
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]

# smooth ONLY for tracing (kill sub-~15 m tillage furrows so D8 doesn't route along them);
# convexity in Step 2 is measured on the ORIGINAL DEM.
ztrace = gaussian_filter(zf.astype(np.float64), sigma=2.5)   # ~12 m
filled = D.fill(ztrace, nodata=-9999.0, cellsize=RES)
fi = D.build_flowinfo(filled, nodata=-9999.0, cellsize=RES)
tpi_large = zf - uniform_filter(zf, size=61)
field = (pen >= 0.45) & (slope < 5)                # open + flat = cultivated field

print(f"{'threshold':>10} {'ridge cells':>11} {'%on highs':>10} {'%furrow(field)':>14} {'%forest':>8}")
results = {}
for T in [50, 100, 200, 500, 1000]:
    divide, lab = D.drainage_divides(fi, threshold=T)
    results[T] = divide
    if divide.sum() == 0:
        print(f"{T:>10}  (none)"); continue
    print(f"{T:>10} {int(divide.sum()):>11} {100*np.mean(tpi_large[divide]>0):>9.0f}% "
          f"{100*np.mean(field[divide]):>13.0f}% {100*np.mean(pen[divide]<0.25):>7.0f}%")

T = 200
ridge = results[T]
print(f"\nCHOSEN threshold={T}: {int(ridge.sum())} ridge cells, "
      f"{100*np.mean(tpi_large[ridge]>0):.0f}% on highs, {100*np.mean(field[ridge]):.0f}% in fields")
np.save("data/derived/elba_fulldensity/ridge_mask.npy", ridge)

hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES)
fig, ax = plt.subplots(figsize=(11, 13))
ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ov = np.zeros((ny, nx, 4)); ov[ridge] = (0.9, 0.1, 0.1, 1.0)
ax.imshow(ov, extent=ext, origin="lower")
ax.set_title(f"Step 1 (S&S divides): ridgelines (red), threshold={T}\n"
             f"{int(ridge.sum())} cells, {100*np.mean(tpi_large[ridge]>0):.0f}% on highs, "
             f"{100*np.mean(field[ridge]):.0f}% in fields")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
fig.savefig("figures/ridgelines_step1.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote figures/ridgelines_step1.png")
