"""Canopy/understory STRUCTURE metrics from the full unclassified 2021 cloud,
and a test of which best predicts the forest-floor measurement offset (dz/dt)
in the lidar DoD.

Streams the 183 M-point cloud in chunks, accumulating per-cell counts and a
per-cell height histogram (for the p95). Excludes class-7 noise. Height above
ground = z - z_after[row,col], with z_after gap-filled by nearest neighbour.

Run:
  cd /home/awickert/projects/lidar-diff-icp
  env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/canopy_struct.py
"""
import numpy as np
import laspy
from scipy import ndimage

LAZ = "data/after/3dep2021_fulldensity.laz"
DERIVED = "data/derived/elba_fulldensity"
OUT = f"{DERIVED}/canopy_struct.npz"

NY, NX = 700, 508
X0, Y0, RES = 577492.8, 4882737.6, 5.0
CELL_AREA = RES * RES  # 25 m^2

# Height histogram bins for p95: 0..40 m in 0.5 m bins, plus overflow.
HBIN = 0.5
HMAX = 40.0
NBINS = int(HMAX / HBIN)  # 80 bins covering [0,40); returns >=40 m go to last bin
NCELLS = NY * NX

# ---- gap-fill ground surface -------------------------------------------------
z_after = np.load(f"{DERIVED}/z_after.npy")
nan = ~np.isfinite(z_after)
if nan.any():
    idx = ndimage.distance_transform_edt(nan, return_distances=False,
                                         return_indices=True)
    z_after = z_after[tuple(idx)]
z_after_flat = z_after.ravel()

# ---- accumulators (float32 counts) -------------------------------------------
total = np.zeros(NCELLS, dtype=np.float64)          # all returns (class != 7)
ground = np.zeros(NCELLS, dtype=np.float64)         # class-2 returns
n_under = np.zeros(NCELLS, dtype=np.float64)        # h in (0.5, 2.0]
n_mid = np.zeros(NCELLS, dtype=np.float64)          # h in (2, 5]
n_veg = np.zeros(NCELLS, dtype=np.float64)          # h > 0.5
first_total = np.zeros(NCELLS, dtype=np.float64)    # first returns
first_lowgap = np.zeros(NCELLS, dtype=np.float64)   # first returns with h < 0.5
hist = np.zeros((NCELLS, NBINS), dtype=np.float32)  # per-cell height histogram

f = laspy.open(LAZ)
processed = 0
for pts in f.chunk_iterator(20_000_000):
    cls = np.asarray(pts.classification)
    keep = cls != 7
    x = np.asarray(pts.x)[keep]
    y = np.asarray(pts.y)[keep]
    z = np.asarray(pts.z)[keep]
    cls = cls[keep]
    rn = np.asarray(pts.return_number)[keep]

    ix = ((x - X0) / RES).astype(np.int64)
    iy = ((y - Y0) / RES).astype(np.int64)
    inb = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
    ix = ix[inb]; iy = iy[inb]; z = z[inb]; cls = cls[inb]; rn = rn[inb]
    lin = iy * NX + ix

    h = z - z_after_flat[lin]

    np.add.at(total, lin, 1)
    np.add.at(ground, lin[cls == 2], 1)

    is_under = (h > 0.5) & (h <= 2.0)
    is_mid = (h > 2.0) & (h <= 5.0)
    is_veg = h > 0.5
    np.add.at(n_under, lin[is_under], 1)
    np.add.at(n_mid, lin[is_mid], 1)
    np.add.at(n_veg, lin[is_veg], 1)

    is_first = rn == 1
    np.add.at(first_total, lin[is_first], 1)
    np.add.at(first_lowgap, lin[is_first & (h < 0.5)], 1)

    # height histogram: clip to [0, HMAX), bin index
    hb = np.clip(h, 0.0, HMAX - 1e-6)
    bi = (hb / HBIN).astype(np.int64)
    bi = np.clip(bi, 0, NBINS - 1)
    # only bin non-negative heights (below-ground returns -> bin 0, but they are
    # ground/noise; include them at bin 0 so p95 is over ALL returns as specified)
    np.add.at(hist, (lin, bi), 1)

    processed += len(z)
    print(f"  processed {processed:,} points", flush=True)

# ---- derive metrics ----------------------------------------------------------
def safe_div(a, b):
    out = np.full_like(a, np.nan, dtype=np.float64)
    m = b > 0
    out[m] = a[m] / b[m]
    return out

understory_frac = safe_div(n_under, total).reshape(NY, NX)
midstory_frac = safe_div(n_mid, total).reshape(NY, NX)
veg_frac = safe_div(n_veg, total).reshape(NY, NX)
ground_return_density = (ground / CELL_AREA).reshape(NY, NX)  # per m^2
low_gap = safe_div(first_lowgap, first_total).reshape(NY, NX)

# p95 of height per cell from the histogram
cum = np.cumsum(hist, axis=1)
tot_h = cum[:, -1]
target = 0.95 * tot_h
# first bin whose cumulative count >= target
p95_bin = np.full(NCELLS, np.nan)
has = tot_h > 0
# searchsorted per row via broadcasting
idxb = np.argmax(cum >= target[:, None], axis=1)
# bin centre in m
p95 = (idxb + 0.5) * HBIN
p95[~has] = np.nan
canopy_height_p95 = p95.reshape(NY, NX)

np.savez(
    OUT,
    canopy_height_p95=canopy_height_p95,
    understory_frac=understory_frac,
    midstory_frac=midstory_frac,
    veg_frac=veg_frac,
    ground_return_density=ground_return_density,
    low_gap=low_gap,
)
print("saved", OUT, flush=True)
