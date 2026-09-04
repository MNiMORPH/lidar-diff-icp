#!/usr/bin/env python3
"""Map gen2 return density over a tile, to see a truncated download as a PICTURE.

A truncated fetch is invisible in a tile average -- whitewater's shipped file read a healthy
11.39 returns/m2 while being 15.45 west of one easting and 5.52 east. It is obvious in a map.
Compares any number of clouds on one colour scale over the same grid.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/gen2_density_map.py \
        --bounds 709531 5323589 711986 5327144 --name cook \
        --clouds data/after_ne/ne_3dep.laz data/after_ne/ne_3dep_fulldensity.laz \
        --labels "before re-fetch" "full density"
"""
import argparse, os
import numpy as np, laspy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--bounds", nargs=4, type=float, required=True)
ap.add_argument("--clouds", nargs="+", required=True)
ap.add_argument("--labels", nargs="+", default=None)
ap.add_argument("--name", required=True)
ap.add_argument("--cell", type=float, default=50.0,
                help="density-map cell, m. MINE: coarse enough that the count is stable, "
                     "fine enough to resolve an octree-node edge.")
ap.add_argument("--figdir", default="figures/sites")
A = ap.parse_args()
B = A.bounds
nx = int(np.ceil((B[2]-B[0])/A.cell)); ny = int(np.ceil((B[3]-B[1])/A.cell))
labels = A.labels or [os.path.basename(c) for c in A.clouds]

grids = []
for c in A.clouds:
    g = np.zeros((ny, nx), np.int64)
    with laspy.open(c) as f:
        for pts in f.chunk_iterator(5_000_000):
            x = np.asarray(pts.x); y = np.asarray(pts.y)
            m = (x >= B[0]) & (x < B[2]) & (y >= B[1]) & (y < B[3])
            if not m.any():
                continue
            np.add.at(g, (((y[m]-B[1])/A.cell).astype(int),
                          ((x[m]-B[0])/A.cell).astype(int)), 1)
    d = g / (A.cell**2)
    grids.append(d)
    v = d[d > 0]
    print(f"{os.path.basename(c):34s} total {g.sum():13,}  per m2: p10 {np.percentile(v,10):6.2f} "
          f"median {np.median(v):6.2f}  p90 {np.percentile(v,90):6.2f}")

vmax = float(np.percentile(np.concatenate([g[g > 0] for g in grids]), 99))
fig, ax = plt.subplots(1, len(grids), figsize=(7.4*len(grids), 7.6), squeeze=False)
ext = (B[0], B[2], B[1], B[3])
for a, d, lb in zip(ax[0], grids, labels):
    im = a.imshow(np.where(d > 0, d, np.nan), extent=ext, origin="lower",
                  cmap="viridis", vmin=0, vmax=vmax)
    a.set_title(f"{lb}\n{d.sum()*A.cell**2:,.0f} pts   median {np.median(d[d>0]):.1f} /m2")
    a.set_xlabel("Easting (m)")
    fig.colorbar(im, ax=a, shrink=0.6, extend="max", label="returns per m2")
ax[0][0].set_ylabel("Northing (m)")
fig.suptitle(f"{A.name}: gen2 return density, {A.cell:g} m cells, common colour scale", y=0.97)
fig.tight_layout(rect=(0, 0, 1, 0.94))
os.makedirs(A.figdir, exist_ok=True)
out = f"{A.figdir}/{A.name}_gen2_density.png"
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote", out)
