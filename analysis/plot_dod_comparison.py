#!/usr/bin/env python3
"""Corrected vs uncorrected DoD, side by side, on one colour scale.

Exists because the standard site figure pairs a DoD with its LoD, and a tile can have a
corrected DoD but NO corrected LoD -- xdem's error model is a Delaunay interpolation that
goes degenerate when the correction shrinks the stable set (whitewater: 60,568 -> 48,498
cells, "initial simplex is flat"). Pairing the corrected DoD with the uncorrected LoD to
fill the panel would be a mislabel, so this figure compares the two DoDs instead.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/plot_dod_comparison.py \
        --tile data/derived/whitewater
"""
import argparse, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lidar_diff_icp.figures import grid_of
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True)
ap.add_argument("--a", default="dod_gen2_median.npy", help="left panel (uncorrected)")
ap.add_argument("--b", default="dod_cover_q2.npy", help="right panel (corrected)")
ap.add_argument("--v", type=float, default=0.3, help="colour limit, m. MINE, matches figures.py")
ap.add_argument("--figdir", default="figures/sites")
A = ap.parse_args()

name = os.path.basename(A.tile.rstrip("/"))
X0, Y0, res, nx, ny = grid_of(A.tile)
Z21 = np.load(f"{A.tile}/z_after.npy")
hs = hillshade(Z21, res, X0, Y0, fill_gaps=False)
ext = (X0, X0 + nx * res, Y0, Y0 + ny * res)

fig, ax = plt.subplots(1, 2, figsize=(15, 9), sharex=True, sharey=True)
for a, fn, ttl in ((ax[0], A.a, "gen2 at its MEDIAN (uncorrected)"),
                   (ax[1], A.b, "gen2 at q(SD), piecewise (corrected)")):
    d = np.load(f"{A.tile}/{fn}")
    f = np.isfinite(d)
    nm = 1.4826 * np.median(np.abs(d[f] - np.median(d[f]))) * 1000
    a.imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im = a.imshow(d, extent=ext, origin="lower", cmap="RdBu", vmin=-A.v, vmax=A.v)
    a.set_title(f"{ttl}\n[{fn}]  {int(f.sum()):,} cells   NMAD {nm:.1f} mm   "
                f"median {np.median(d[f])*1000:+.1f} mm")
    a.set_xlabel("Easting (m)")
    fig.colorbar(im, ax=a, shrink=0.55, extend="both")
ax[0].set_ylabel("Northing (m)")
fig.suptitle(f"{name}: DEM of Difference, gen2 - gen1 (m).  red = erosion, blue = deposition",
             y=0.94)
os.makedirs(A.figdir, exist_ok=True)
out = f"{A.figdir}/{name}_dod_corrected_vs_median.png"
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote", out)
