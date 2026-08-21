#!/usr/bin/env python3
"""RECOVERABILITY TEST (Task 4, gate) — is there ground under class-2 to recover?

The +tan(slope) DoD bias is leaf-on canopy ground-sparsity: on forested slopes the
2021 (green-up) class-2 ground surface reads HIGH because true ground is under-sampled.
Before building any robust under-canopy filter (Kraus-Pfeifer), answer the cheap
prerequisite question:

    Among gen2 pulses that PENETRATED the canopy (number_of_returns > 1), does the
    LAST echo sit BELOW the class-2 ground surface?

  - If penetrating last echoes (esp. the NOT-ground ones, class != 2) cluster a modest
    depth BELOW the class-2 surface on forested slopes -> there is lower coherent ground
    a robust filter could pull to => recovery is POSSIBLE.
  - If they sit AT or ABOVE class-2 (class-2 already found the floor) -> nothing lower to
    recover => the honest answer is COARSEN + LoD, not re-filtering.

Non-circular: the class-2 surface is gridded from vendor class-2 points; we then ask
whether a DIFFERENT set (class-1 last echoes) falls below it. Control: open/flat cells
(ground well sampled) should show NO sub-ground last echoes.

Memory-safe: two chunked passes over the 183 M-point cloud (never holds it whole).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/recoverability_test.py
"""
import laspy, numpy as np, pandas as pd
from scipy.ndimage import distance_transform_edt, uniform_filter

FN = "data/after/3dep2021_fulldensity.laz"
RES = 2.0            # ground-reference grid cell (m)
GQ = 0.20           # robust-low percentile for the class-2 ground surface
STEEP = 15.0        # deg; slope split
LOWPEN = 0.25       # ground/total return fraction below which a cell is "canopy-starved"
CHUNK = 20_000_000

h = laspy.open(FN).header
x0, y0 = h.mins[0], h.mins[1]
nx = int((h.maxs[0] - x0) / RES) + 1
ny = int((h.maxs[1] - y0) / RES) + 1
NC = nx * ny
print(f"grid {nx} x {ny} = {NC:,} cells at {RES} m", flush=True)

def cellidx(px, py):
    ix = np.clip(((px - x0) / RES).astype(np.int64), 0, nx - 1)
    iy = np.clip(((py - y0) / RES).astype(np.int64), 0, ny - 1)
    return iy * nx + ix

# ---- PASS 1: per-cell class-2 z (for the ground surface) + per-cell return totals -------
# For p20 without storing everything: keep only class-2 (cell,z) pairs (~51 M) streamed
# into lists, plus running total/ground counts for penetration.
tot = np.zeros(NC, np.int64)
grd = np.zeros(NC, np.int64)
c2_cell_parts, c2_z_parts = [], []
print("PASS 1: class-2 surface + penetration ...", flush=True)
with laspy.open(FN) as fh:
    for pts in fh.chunk_iterator(CHUNK):
        cl = np.asarray(pts.classification)
        keep = cl != 7
        cx = np.asarray(pts.x)[keep]; cy = np.asarray(pts.y)[keep]
        cz = np.asarray(pts.z)[keep].astype(np.float32); cl = cl[keep]
        c = cellidx(cx, cy)
        tot += np.bincount(c, minlength=NC)
        g2 = cl == 2
        grd += np.bincount(c[g2], minlength=NC)
        c2_cell_parts.append(c[g2].astype(np.int32))
        c2_z_parts.append(cz[g2])
        print(f"  ...{tot.sum():,} pts", flush=True)

c2_cell = np.concatenate(c2_cell_parts); c2_z = np.concatenate(c2_z_parts)
del c2_cell_parts, c2_z_parts
zg = pd.Series(c2_z).groupby(c2_cell).quantile(GQ)
del c2_cell, c2_z
grid = np.full(NC, np.nan, np.float32); grid[zg.index.values] = zg.values
grid = grid.reshape(ny, nx)
nanmask = ~np.isfinite(grid)
filled = grid.copy()
if nanmask.any():
    filled = filled[tuple(distance_transform_edt(nanmask, return_distances=False,
                                                 return_indices=True))]
sm = uniform_filter(filled, size=3)
gy, gx = np.gradient(sm, RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
with np.errstate(invalid="ignore", divide="ignore"):
    pen = (grd / np.maximum(tot, 1)).reshape(ny, nx)
zref_flat = filled.ravel(); slope_flat = slope.ravel(); pen_flat = pen.ravel()
print(f"  ground surface built ({np.isfinite(grid).sum():,} cells populated)", flush=True)
np.savez("data/derived/elba/recov_pass1.npz", zref=filled, slope=slope, pen=pen,
         nx=nx, ny=ny)

# ---- PASS 2: residual of penetrating last echoes vs the class-2 surface -----------------
# 2x2 zones ISOLATE slope from canopy: steep x {forest,open}, flat x {forest,open}.
# If steep+OPEN also shows deep sub-ground last echoes -> it's within-cell slope geometry,
# NOT canopy float -> recovery claim weakens. If only steep+FOREST is deep -> canopy float.
STEEP_ = STEEP; FOREST = LOWPEN; OPEN = 0.40
def zones(c):
    st = slope_flat[c] >= STEEP_; fl = ~st
    fo = pen_flat[c] < FOREST;    op = pen_flat[c] >= OPEN
    return {"steepForest": st & fo, "steepOpen": st & op,
            "flatForest": fl & fo,  "flatOpen": fl & op}

ZK = ["steepForest", "steepOpen", "flatForest", "flatOpen"]
res_lists = {f"{z}_{g}": [] for z in ZK for g in ("ng", "g")}
cellmin = {z: [] for z in ZK}     # (cell,resid) stacks for per-cell deepest
print("PASS 2: penetrating-last-echo residuals ...", flush=True)
with laspy.open(FN) as fh:
    for pts in fh.chunk_iterator(CHUNK):
        cl = np.asarray(pts.classification)
        rn = np.asarray(pts.return_number); nr = np.asarray(pts.number_of_returns)
        sel = (cl != 7) & (nr > 1) & (rn == nr)           # penetrating last echo
        if not sel.any():
            continue
        cz = np.asarray(pts.z)[sel].astype(np.float32)
        c = cellidx(np.asarray(pts.x)[sel], np.asarray(pts.y)[sel])
        cl = cl[sel]; resid = cz - zref_flat[c]
        zm = zones(c)
        for z in ZK:
            for gk, gm in (("ng", cl != 2), ("g", cl == 2)):
                m = zm[z] & gm
                res_lists[f"{z}_{gk}"].append(resid[m])
            mm = zm[z]
            cellmin[z].append(np.stack([c[mm].astype(np.int64), resid[mm]]))
        print("  ...chunk done", flush=True)

def report(tag, arr):
    if arr.size == 0:
        print(f"  {tag:36s} n=0"); return
    below = arr < -0.05
    print(f"  {tag:36s} n={arr.size:>10,}  median={np.median(arr):+.3f}  "
          f"p10={np.percentile(arr,10):+.3f}  frac<-5cm={below.mean():5.1%}  "
          f"med|below={np.median(arr[below]) if below.any() else float('nan'):+.3f}")

print("\n=== residual of penetrating last echo vs class-2 p20 ground surface (m) ===")
print("    (+ above class-2 ground; - BELOW it = lower ground to recover)")
print("    2x2 isolates SLOPE from CANOPY. Compare steepForest vs steepOpen.\n")
for z in ZK:
    print(f"[{z}]")
    report("last echo NON-ground (class!=2)", np.concatenate(res_lists[f"{z}_ng"]))
    report("last echo ground (class==2)",     np.concatenate(res_lists[f"{z}_g"]))
    print()

print("=== per-cell: deepest penetrating-last-echo vs class-2 surface ===")
for z in ZK:
    cr = np.concatenate(cellmin[z], axis=1)
    lo = pd.Series(cr[1]).groupby(cr[0].astype(np.int64)).min().values
    print(f"  {z:12s}: n={lo.size:>8,} cells  "
          f"frac deepest >10cm below class-2 = {np.mean(lo < -0.10):5.1%}  "
          f"median deepest = {np.median(lo):+.3f} m")

print("\nVERDICT GUIDE: canopy-float recovery is REAL only if steepForest shows sub-class-2 "
      "last echoes that steepOpen does NOT (same slope, differ only in canopy). If steepOpen "
      "is just as deep, the sub-ground signal is within-cell slope geometry, not recoverable "
      "ground -> answer = coarsen+LoD.")
