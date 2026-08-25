#!/usr/bin/env python3
"""Render a saved DoD in the project's standard figure form.

Mirrors the figure that scripts/gridded_ground_dod.py produces at the end of a pipeline
run -- same layout, colormap, limits and titles -- but reads the saved dod.npy/lod.npy so an
existing product can be plotted without re-running difference_dem. The convention it follows
is the one fixed in the README: DoD is after - before, RED = EROSION, BLUE = DEPOSITION, on
a standard NW (315/45) hillshade.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/plot_dod.py --tile data/derived/elba_fulldensity
"""
import argparse, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--dod", default="dod.npy")
ap.add_argument("--lod", default="lod.npy")
ap.add_argument("--v", type=float, default=0.3, help="DoD colour limit (m), as in the CLI figure")
ap.add_argument("--out", default=None)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))

for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
    p = f"{A.tile}/{fn}"
    if os.path.exists(p):
        j = json.load(open(p)); b = j["bounds"]; res = float(j.get("res") or j.get("res_m"))
        X0, Y0 = b[0], b[1]; break
else:
    raise SystemExit(f"no grid meta in {A.tile}")

dod = np.load(f"{A.tile}/{A.dod}")
lod = np.load(f"{A.tile}/{A.lod}")
Z21 = np.load(f"{A.tile}/z_after.npy")
ny, nx = dod.shape
ok = np.isfinite(dod)
print(f"{TILE}: {ok.sum():,} cells   median {1000*np.median(dod[ok]):+.1f} mm   "
      f"NMAD {1000*1.4826*np.median(np.abs(dod[ok]-np.median(dod[ok]))):.1f} mm")
print(f"   |DoD| > LoD on {100*np.mean(np.abs(dod[ok]) > lod[ok]):.1f}% of cells")

hs = hillshade(Z21, res, X0, Y0)                # GDAL NW shaded relief, oriented by the geotransform
ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = A.v
fig, ax = plt.subplots(1, 2, figsize=(15, 9))
ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
ax[0].set_title("DEM of Difference (gridded ground): gen2 − gen1 (m)\nred = erosion, blue = deposition")
fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
ax[1].set_title("level of detection (m)")
fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
out = A.out or f"figures/dod_{TILE}.png"
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"wrote {out}")
