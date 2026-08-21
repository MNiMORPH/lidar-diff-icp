#!/usr/bin/env python3
"""COARSEN RESOLUTION SIZING (Task 4) — at what cell size do closed-canopy steep cells
hold enough ground returns that a few pulses got through the veg?

Andy's question: pick the coarsening resolution from the DATA. The binding constraint is
under-canopy ground sparsity on the forested slopes where the 2021 (green-up) leaf-on bias
lives. So measure, as a function of cell size, the number of class-2 GROUND returns per
cell on the HARDEST cells (closed canopy AND steep), and find where "at least a few" get
through for essentially all of them.

  ground return per cell = a pulse that penetrated to ground => "got through the veg".

For each resolution: bin ground (class 2) and total returns; penetration = ground/total;
slope from the cached 2 m surface, block-averaged. Hard cells = slope>=15 deg AND
penetration<0.25. Report the ground-count distribution there.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/coarsen_resolution.py
"""
import laspy, numpy as np
from scipy.ndimage import uniform_filter

FN = "data/after/3dep2021_fulldensity.laz"
RES_LIST = [5, 10, 12, 15, 20, 25, 30]
STEEP = 15.0
LOWPEN = 0.25
CHUNK = 20_000_000

h = laspy.open(FN).header
x0, y0 = h.mins[0], h.mins[1]
xmax, ymax = h.maxs[0], h.maxs[1]

# slope from cached 2 m ground surface
p1 = np.load("data/derived/elba/recov_pass1.npz")
slope2 = p1["slope"]                      # (ny2, nx2) at 2 m
RES2 = 2.0

def dims(res):
    nx = int((xmax - x0) / res) + 1
    ny = int((ymax - y0) / res) + 1
    return nx, ny

# accumulate ground + total bincounts per resolution in one chunked pass
acc = {res: {"g": np.zeros(dims(res)[0] * dims(res)[1], np.int64),
             "t": np.zeros(dims(res)[0] * dims(res)[1], np.int64)} for res in RES_LIST}
print("counting ground/total per cell across resolutions ...", flush=True)
with laspy.open(FN) as fh:
    for pts in fh.chunk_iterator(CHUNK):
        cl = np.asarray(pts.classification)
        keep = cl != 7
        px = np.asarray(pts.x)[keep]; py = np.asarray(pts.y)[keep]; cl = cl[keep]
        g2 = cl == 2
        for res in RES_LIST:
            nx, ny = dims(res)
            ix = np.clip(((px - x0) / res).astype(np.int64), 0, nx - 1)
            iy = np.clip(((py - y0) / res).astype(np.int64), 0, ny - 1)
            c = iy * nx + ix
            acc[res]["t"] += np.bincount(c, minlength=nx * ny)
            acc[res]["g"] += np.bincount(c[g2], minlength=nx * ny)
        print(f"  ...chunk", flush=True)

def slope_coarse(res):
    """Block-mean the 2 m slope grid onto the res grid (approximate, for masking)."""
    nx, ny = dims(res)
    f = int(round(res / RES2))
    ny2, nx2 = slope2.shape
    Ny, Nx = ny * f, nx * f
    pad = np.full((Ny, Nx), np.nan, np.float32)
    pad[:min(Ny, ny2), :min(Nx, nx2)] = slope2[:min(Ny, ny2), :min(Nx, nx2)]
    return np.nanmean(pad.reshape(ny, f, nx, f), axis=(1, 3))

print("\n=== ground returns per cell on CLOSED-CANOPY STEEP cells "
      "(slope>=15 deg, penetration<0.25) ===")
print(f"{'res(m)':>6} {'hardcells':>10} {'%>=1':>6} {'%>=3':>6} {'%>=6':>6} "
      f"{'p10':>5} {'p25':>5} {'median':>7}")
for res in RES_LIST:
    nx, ny = dims(res)
    g = acc[res]["g"].reshape(ny, nx).astype(float)
    t = acc[res]["t"].reshape(ny, nx).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        pen = g / np.maximum(t, 1)
    sc = slope_coarse(res)
    hard = (sc >= STEEP) & (pen < LOWPEN) & (t > 0)
    gc = g[hard]
    if gc.size == 0:
        print(f"{res:6d}  no hard cells"); continue
    print(f"{res:6d} {gc.size:10.0f} {100*np.mean(gc>=1):5.0f}% {100*np.mean(gc>=3):5.0f}% "
          f"{100*np.mean(gc>=6):5.0f}% {np.percentile(gc,10):5.0f} {np.percentile(gc,25):5.0f} "
          f"{np.median(gc):7.0f}")

print("\nRead: pick the smallest res where nearly all hard cells clear the 'few pulses "
      "through the veg' bar (e.g. p10>=3-6, %>=6 near 100).")
