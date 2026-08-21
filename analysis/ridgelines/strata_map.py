#!/usr/bin/env python3
"""Map the land-cover strata used throughout the return-distribution analysis:
  FOREST = penetration < 0.25 ; OPEN = penetration >= 0.45 ; intermediate (0.25-0.45) and
  floodplain are excluded.  Penetration = per-cell ground-return fraction.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/strata_map.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from lidar_diff_icp.viz import hillshade

D = "data/derived/elba_fulldensity/"
RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load(D + "z_after.npy"); ny, nx = z.shape
pen = np.load(D + "penetration.npy"); fld = np.load(D + "floodplain_mask.npy").astype(bool)

forest = (pen < 0.25) & ~fld & np.isfinite(pen)
open_  = (pen >= 0.45) & ~fld & np.isfinite(pen)
inter  = (pen >= 0.25) & (pen < 0.45) & ~fld & np.isfinite(pen)

ext = (X0, X0+nx*RES, Y0, Y0+ny*RES)
hs = hillshade(z, RES, X0, Y0, fill_gaps=True)
ov = np.zeros((ny, nx, 4))
ov[forest] = (0.13, 0.55, 0.13, 0.55)   # forest green
ov[open_]  = (0.93, 0.79, 0.30, 0.55)   # open tan
ov[inter]  = (0.55, 0.55, 0.55, 0.45)   # intermediate grey
ov[fld]    = (0.20, 0.45, 0.85, 0.40)   # floodplain blue (excluded)

fig, ax = plt.subplots(figsize=(11, 13))
ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ax.imshow(ov, extent=ext, origin="lower")
ax.set_xlabel("Easting (m, UTM 15N)"); ax.set_ylabel("Northing (m)")
n_all = forest.sum()+open_.sum()+inter.sum()+fld.sum()
ax.set_title(f"Analysis strata (penetration = ground-return fraction)\n"
             f"forest {forest.sum():,} | open {open_.sum():,} | intermediate {inter.sum():,} | floodplain {fld.sum():,} cells")
ax.legend(handles=[
    Patch(facecolor=(0.13,0.55,0.13,.7), label="FOREST  (pen < 0.25)"),
    Patch(facecolor=(0.93,0.79,0.30,.7), label="OPEN  (pen ≥ 0.45)"),
    Patch(facecolor=(0.55,0.55,0.55,.6), label="intermediate (excluded)"),
    Patch(facecolor=(0.20,0.45,0.85,.6), label="floodplain (excluded)"),
], loc="upper right", fontsize=9)
fig.savefig("figures/refdatum/strata_map.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote figures/refdatum/strata_map.png")
