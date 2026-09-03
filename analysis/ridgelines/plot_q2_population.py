#!/usr/bin/env python3
"""Map WHERE a q2 fit was calculated: the cells behind it, on the tile's own hillshade.

The population is not re-specified here. Every selection setting is read back out of the
fit's own JSON -- cover layer, valley treatment, per-cell minima -- and the cuts are then
re-applied with the same reference_cells() call the fit used. If this map and the fit ever
disagreed, that would be a defect rather than a plotting choice, so there is no way to pass
a different population on the command line.

Cells are coloured by the covariate that was fitted, so the map answers two questions at
once: which ground the relation rests on, and whether the cover range is spread over it or
concentrated in one corner.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/plot_q2_population.py --tile data/derived/whitewater \
        --fit q2_cover_fit_lowveg.json
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyarrow.parquet as pq

from lidar_diff_icp.refcells import reference_cells
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True)
ap.add_argument("--fit", default="q2_cover_fit_lowveg.json",
                help="the fit JSON whose population is to be drawn")
ap.add_argument("--out", default=None)
A = ap.parse_args()

D = A.tile.rstrip("/")
SITE = os.path.basename(D)
cal = json.load(open(os.path.join(D, A.fit)))
st = cal["settings"]
cover_layer = cal.get("cover_layer", "canopy_cover_pfs.npy")

z = np.load(f"{D}/z_after.npy")
NY, NX = z.shape
# Grid origin, so the hillshade is illuminated geographically NW rather than in array space.
_meta = next(json.load(open(f"{D}/{f}")) for f in ("meta.json", "corrections_geoid.json",
                                                   "corrections.json")
             if os.path.exists(f"{D}/{f}"))
_b = _meta["bounds"]; RES = float(_meta.get("res") or _meta.get("res_m"))
X0, Y0 = float(_b[0]), float(_b[1])
cube = np.load(f"{D}/nearground_cells_sn.npz")
cells = cube["cells"]
cover = np.load(f"{D}/{cover_layer}").ravel()[cells]

vt = st["valley_top_m"]
vt = vt if vt == "antimode" else float(vt)
stable, rep = reference_cells(D, cells=cells, slope_max=90.0,
                              exclude_valley=not st["include_valley"], valley_top_m=vt)

# The two per-cell minima, applied exactly as the fit applies them.
t = pq.read_table(f"{D}/beam_offset_table.parquet", columns=["cell", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
n1 = np.bincount(t["cell"].to_numpy()[g], minlength=z.size)
ng = np.load(f"{D}/nearground_gen2_class_split.npz")["Hg"].sum(1)
ok = (stable & (n1[cells] >= max(1, st["min_gen1_returns"]))
      & (ng >= max(1, st["min_gen2_returns"])) & np.isfinite(cover))

n_fit = int(ok.sum())
if n_fit != cal["population"]["cells_fitted"]:
    print(f"WARNING: this map has {n_fit:,} cells but the fit recorded "
          f"{cal['population']['cells_fitted']:,}. The tile's products have changed since "
          f"the fit was made; the map is of the CURRENT products, not of that fit.")

grid = np.full(z.size, np.nan)
grid[cells[ok]] = cover[ok]
grid = grid.reshape(NY, NX)

zf = z.copy()
bad = ~np.isfinite(zf)
if bad.any():
    from scipy.ndimage import distance_transform_edt as edt
    zf = zf[tuple(edt(bad, return_distances=False, return_indices=True))]

fig, ax = plt.subplots(figsize=(8.4, 8.4 * NY / NX))
ax.imshow(hillshade(zf, RES, X0, Y0), cmap="gray", origin="lower", vmin=0, vmax=1)
im = ax.imshow(np.ma.masked_invalid(grid), origin="lower", cmap="viridis",
               vmin=0.0, vmax=float(np.nanpercentile(cover[ok], 98)))
cb = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
cb.set_label(os.path.splitext(cover_layer)[0], fontsize=9)

cuts = ", ".join(f"{k} −{v:,}" for k, v in rep.items()
                 if k not in ("start", "kept") and v)
ax.set_title(f"{SITE}: the {n_fit:,} cells the q2({os.path.splitext(cover_layer)[0]}) fit "
             f"rests on, of {rep['start']:,} in the near-ground cube\n"
             f"cuts applied: {cuts}", fontsize=8.5, loc="left")
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
out = A.out or f"analysis/ridgelines/q2_population_{SITE}_{os.path.splitext(cover_layer)[0]}.png"
fig.savefig(out, dpi=140)
print(f"cells drawn: {n_fit:,} of {rep['start']:,} in the cube")
print(f"cuts: {cuts}")
print(f"wrote {out}")
