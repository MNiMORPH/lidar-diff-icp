#!/usr/bin/env python3
"""The standard site figures, plus the vegetation-LoD landscape masking.

Andy, 2026-09-16: "could you produce the standard set of figures, which includes the new
error-related landscape masking?"

WHAT THE MASKING IS. gen1 cannot characterise its own vegetation -- 16 returns per 5 m cell
in the -1..+2 m band against gen2's 224, 97.8% single-return -- so on vegetated flat
floodplain the 2008 ground is NOT MEASURED and a ~78 mm LoD over it is false precision.
`vegetation_lod.inflate_lod` adds gen1's own per-cell lift (median - p05 of its slope-normal
offsets) in QUADRATURE, on FLAT FLOODPLAIN ONLY, so bank erosion keeps the LoD it had.

Writes `lod_veg.npy` / `change_veg.npy` BESIDE the base products, never over them, and draws
the two standard figures from them plus a third showing where the inflation applied and what
it silenced. The base `dod.npy` is UNCHANGED -- this is an error statement, not a correction.

    ./lidar-icp/bin/python scripts/figures_with_vegetation_lod.py \
        --tile data/derived/elbaext_independent --slope-max-deg 2
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lidar_diff_icp import figures, vegetation_lod as V
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True)
ap.add_argument("--lift", required=True, help="per-cell gen1 lift raster, mm")
ap.add_argument("--slope-max-deg", type=float, required=True,
                help="MINE, and stated: where the inflation applies. It decides how much "
                     "flat floodplain stops being measurable. No default.")
ap.add_argument("--figdir", default="figures/sites")
A = ap.parse_args()

D = A.tile.rstrip("/")
name = os.path.basename(D)
dod = np.load(f"{D}/dod.npy"); lod = np.load(f"{D}/lod.npy")
stable = np.load(f"{D}/stable.npy") if os.path.exists(f"{D}/stable.npy") else None
fp = np.load(f"{D}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{D}/slope.npy") if os.path.exists(f"{D}/slope.npy") else \
     np.load("data/derived/elbaext/slope.npy")
lift = np.load(A.lift)
res = float(json.load(open(f"{D}/corrections.json"))["res_m"])

lod_v, applied, info = V.inflate_lod(lod, lift, fp, sl, slope_max_deg=A.slope_max_deg)
np.save(f"{D}/lod_veg.npy", lod_v)

det = detect_change_standard(dod, lod_v, stable if stable is not None else np.zeros_like(dod, bool), res)
change_v = det["change"]
np.save(f"{D}/change_veg.npy", change_v)
base = np.load(f"{D}/change.npy").astype(bool)
silenced = base & ~change_v.astype(bool)
gained = change_v.astype(bool) & ~base
with open(f"{D}/lod_veg.json", "w") as fh:
    json.dump({**info, "regions_before": None, "regions_after": len(det["regions"]),
               "cells_detected_before": int(base.sum()),
               "cells_detected_after": int(change_v.sum()),
               "cells_silenced": int(silenced.sum()), "cells_gained": int(gained.sum()),
               "lift_raster": A.lift,
               "note": "an ERROR statement, not a correction: dod.npy is unchanged"},
              fh, indent=1)
    fh.write("\n")

fa = figures.dod_lod_figure(D, A.figdir, name, lod_name="lod_veg.npy", suffix="_veg")
fb = figures.change_figure(D, A.figdir, name, change_name="change_veg.npy", suffix="_veg")

# third panel: WHERE the masking acts, and what it removed
X0, Y0, r, nx, ny = figures.grid_of(D)
Z21 = np.load(f"{D}/z_after.npy")
hs = hillshade(Z21, r, X0, Y0, fill_gaps=True)
ext = (X0, X0 + nx * r, Y0, Y0 + ny * r)
fig, ax = plt.subplots(1, 2, figsize=(19, 8.5))
ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray")
im = ax[0].imshow(np.where(applied, 1000 * (lod_v - lod), np.nan), extent=ext,
                  origin="lower", cmap="magma", vmin=0, vmax=250)
ax[0].set_title(f"{name}: vegetation LoD inflation (mm)\n"
                f"flat floodplain, slope <= {A.slope_max_deg:g} deg; "
                f"{info['cells_inflated']:,} cells, median lift {info['median_lift_mm']:.0f} mm")
fig.colorbar(im, ax=ax[0], shrink=0.6, extend="max", label="LoD added (mm)")
ax[1].imshow(hs, extent=ext, origin="lower", cmap="gray")
ax[1].imshow(np.where(silenced, 1.0, np.nan), extent=ext, origin="lower",
             cmap="autumn", vmin=0, vmax=1, alpha=0.85)
ax[1].set_title(f"detections SILENCED by the vegetation LoD\n"
                f"{int(silenced.sum()):,} cells removed, {int(gained.sum()):,} gained; "
                f"{int(base.sum()):,} -> {int(change_v.sum()):,} detected")
for a_ in ax:
    a_.set_xlabel("Easting (m)"); a_.set_ylabel("Northing (m)")
out = f"{A.figdir}/{name}_vegetation_lod.png"
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"  {int(base.sum()):,} -> {int(change_v.sum()):,} detected cells "
      f"({int(silenced.sum()):,} silenced, {int(gained.sum()):,} gained)")
print(f"  wrote {fa}\n        {fb}\n        {out}")
