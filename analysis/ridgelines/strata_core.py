#!/usr/bin/env python3
"""Robust ("core") land-cover strata: only the PURE INTERIOR of significant homogeneous
clusters, to minimize cross-cover contamination for the return-distribution analysis.

Keeps the original per-cell masks (forest pen<0.25, open pen>=0.45) and ADDS core masks:
  1. Reclassify contamination: open cells carrying above-ground vegetation (veg_frac > VEG_MAX)
     are trees/structures embedded in fields -> demote to 'other' so they seed buffers too.
  2. Erode each class from ANY non-same cover (distance-transform >= R_CELLS): a cell survives
     only if everything within R is the same class -> not marginal to another cover.
  3. Keep only connected components >= A_MIN cells -> significant clusters only.
Saves core_forest.npy / core_open.npy; maps core vs original; validates purity.

Grid geometry comes from the tile's own corrections.json. Neither input is optional:
`penetration.npy` DEFINES the two classes (forest pen<0.25, open pen>=0.45) and
`canopy_struct.npz` DEFINES the tree-in-open demotion, so a run missing either would
produce masks that are not the masks this script names. Both refuse rather than default.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/strata_core.py \
        --tile elba_fulldensity
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import distance_transform_edt, label as cclabel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from lidar_diff_icp.viz import hillshade

_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default="elba_fulldensity",
                 help="tile name under data/derived/; grid geometry is read from its "
                      "corrections.json, so no origin is hardcoded")
ARGS = _ap.parse_args()
TILE = ARGS.tile
D = f"data/derived/{TILE}/"


def _grid(tile):                                    # (X0,Y0,RES) from the tile's own meta
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]
            return b[0], b[1], float(j.get("res") or j.get("res_m"))
    raise SystemExit(f"no grid meta for {tile}: none of meta.json, corrections_geoid.json, "
                     f"corrections.json exists under data/derived/{tile}")


def _req(name, why):
    """Both inputs DEFINE the strata, so neither has a run-without form: a substitute would
    silently rename the product rather than reproduce it."""
    p = D + name
    if not os.path.exists(p):
        raise SystemExit(f"{p} is missing. {why} There is no --without for it, because the "
                         f"core masks are DEFINED by it; produce it for {TILE} first.")
    return p


X0, Y0, RES = _grid(TILE)
R_CELLS = 3          # 15 m purity radius: all cells within R must be same cover
A_MIN   = 50         # >= 50 cells (1250 m^2) connected -> a significant cluster
VEG_MAX = 0.02       # open cells with veg_frac above this = tree/structure -> demote

z = np.load(D + "z_after.npy"); ny, nx = z.shape
pen = np.load(_req("penetration.npy", "It sets the forest/open cut (pen<0.25 / pen>=0.45)."))
fld = np.load(_req("floodplain_mask.npy",
                   "It removes valley-bottom cells from both classes.")).astype(bool)
veg = np.load(_req("canopy_struct.npz",
                   "Its veg_frac demotes trees embedded in fields."))["veg_frac"]
fin = np.isfinite(pen)

forest0 = (pen < 0.25) & ~fld & fin                       # original masks (kept)
open0   = (pen >= 0.45) & ~fld & fin
tree_in_open = open0 & (veg > VEG_MAX)                     # embedded trees/structures

# class label: 1 forest, 2 clean-open, 0 other (intermediate, floodplain, nodata, tree-in-open)
cls = np.zeros((ny, nx), np.int8)
cls[forest0] = 1
cls[open0 & ~tree_in_open] = 2

def core(mask, r, amin):
    dist = distance_transform_edt(mask)                   # cells: distance to nearest non-mask
    interior = mask & (dist >= r)                         # pure interior (erosion by r)
    lab, n = cclabel(interior)
    sizes = np.bincount(lab.ravel())
    big = sizes >= amin; big[0] = False
    return interior & big[lab]

core_forest = core(cls == 1, R_CELLS, A_MIN)
core_open   = core(cls == 2, R_CELLS, A_MIN)
np.save(D + "core_forest.npy", core_forest); np.save(D + "core_open.npy", core_open)

print(f"original : forest {forest0.sum():>7,}  open {open0.sum():>7,}")
print(f"core     : forest {core_forest.sum():>7,}  open {core_open.sum():>7,}  "
      f"(kept {core_forest.sum()/forest0.sum()*100:.0f}% / {core_open.sum()/open0.sum()*100:.0f}%)")
print(f"tree-in-open cells demoted: {tree_in_open.sum():,}")
# purity validation: above-ground content
print("\nPURITY (mean veg_frac, fraction of cells with any tall veg p95>3m):")
p95 = np.load(D + "canopy_struct.npz")["canopy_height_p95"]
for m, lbl in [(open0,"open orig"),(core_open,"open CORE"),(forest0,"forest orig"),(core_forest,"forest CORE")]:
    print(f"  {lbl:12s}: mean veg_frac {veg[m].mean():.4f}   cells p95>3m {(p95[m]>3).mean()*100:4.1f}%   "
          f"median pen {np.median(pen[m]):.2f}")

# ---- map: original (faint) with core (bold) overlaid --------------------------------
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES); hs = hillshade(z, RES, X0, Y0, fill_gaps=True)
fig, axes = plt.subplots(1, 2, figsize=(20, 12), sharex=True, sharey=True)
for ax, mode in zip(axes, ["original", "core"]):
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    ov = np.zeros((ny, nx, 4))
    if mode == "original":
        ov[forest0] = (0.13,0.55,0.13,0.55); ov[open0] = (0.93,0.79,0.30,0.55)
        ttl = f"ORIGINAL  (forest {forest0.sum():,} / open {open0.sum():,})"
    else:
        ov[forest0] = (0.13,0.55,0.13,0.15); ov[open0] = (0.93,0.79,0.30,0.15)   # faint context
        ov[core_forest] = (0.05,0.45,0.05,0.85); ov[core_open] = (0.95,0.70,0.10,0.85)
        ttl = f"CORE  (forest {core_forest.sum():,} / open {core_open.sum():,})  R={R_CELLS*RES:.0f}m, cluster>={A_MIN}"
    ax.imshow(ov, extent=ext, origin="lower")
    ax.set_title(ttl); ax.set_xlabel("Easting (m)")
axes[0].set_ylabel("Northing (m)")
axes[1].legend(handles=[Patch(facecolor=(0.05,0.45,0.05,.85),label="core FOREST"),
                        Patch(facecolor=(0.95,0.70,0.10,.85),label="core OPEN"),
                        Patch(facecolor=(0.6,0.6,0.6,.4),label="original (trimmed away)")],
               loc="upper right", fontsize=10)
fig.suptitle("Robust core strata: pure interior of significant homogeneous clusters", y=0.99)
_fig = ("figures/refdatum/strata_core.png" if TILE == "elba_fulldensity"
        else f"figures/refdatum/strata_core_{TILE}.png")
fig.savefig(_fig, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {_fig} ; saved {D}core_forest.npy, {D}core_open.npy")
