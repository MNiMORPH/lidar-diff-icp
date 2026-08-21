"""Per-return slope-normal distance to the local ground surface, binned per cell.

For every lidar return in both epochs (gen2 2021 3DEP, gen1 2008 MN), compute the
perpendicular distance to the gen2 bare-earth ground plane (a common reference frame),
and accumulate per-cell histograms of that distance. The histograms are the reusable
basis; compact per-cell summaries are stored alongside for quick use.

Run:
  cd /home/awickert/projects/lidar-diff-icp && \
  env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/slope_normal_returns.py
"""
import numpy as np
import laspy
from scipy.ndimage import distance_transform_edt

# ---- grid ----
NY, NX = 700, 508
X0, Y0 = 577492.8, 4882737.6
RES = 5.0

# ---- histogram bins (the saved reusable basis) ----
BIN_LO, BIN_HI, BIN_W = -1.0, 40.0, 0.25
edges = np.arange(BIN_LO, BIN_HI + 0.5 * BIN_W, BIN_W)  # inclusive of 40.0
NBINS = edges.size - 1
print(f"bins: {NBINS} from {edges[0]} to {edges[-1]} step {BIN_W}")

# ---- fine internal ground-band bins (NOT saved) ----
# The coarse 0.25 m basis is too granular to resolve the per-cell ground median/p10
# offset (which is a few-cm signal, and the differenceable quantity between epochs).
# We accumulate ground returns in a narrow, fine (1 cm) histogram to recover the
# ground-return position summaries at ~1 cm precision.
FLO, FHI, FW = -1.5, 2.5, 0.01
fedges = np.arange(FLO, FHI + 0.5 * FW, FW)
NFBINS = fedges.size - 1
fcentres = 0.5 * (fedges[:-1] + fedges[1:])
print(f"fine ground bins: {NFBINS} from {fedges[0]} to {fedges[-1]} step {FW}")

CHUNK = 20_000_000

GEN2 = "data/after/3dep2021_fulldensity.laz"
GEN1 = "data/before/4342-29-64.laz"
OUT = "data/derived/elba_fulldensity/slope_normal_returns.npz"

# ---- reference ground plane (gen2 bare earth) ----
Zg = np.load("data/derived/elba_fulldensity/z_after.npy")  # (NY,NX) float64
assert Zg.shape == (NY, NX), Zg.shape
nanmask = np.isnan(Zg)
# fill NaNs by nearest valid neighbour
_, (iy_src, ix_src) = distance_transform_edt(nanmask, return_indices=True)
Zg_filled = Zg.copy()
Zg_filled[nanmask] = Zg[iy_src[nanmask], ix_src[nanmask]]

# per-cell slope components (m/m): gx=d/dEast(cols), gy=d/dNorth(rows)
gy, gx = np.gradient(Zg_filled, RES)
cos_slope = 1.0 / np.sqrt(1.0 + gx * gx + gy * gy)  # (NY,NX)

# cell-centre coordinates
ix_col = np.arange(NX)
iy_row = np.arange(NY)
xc_of_col = X0 + (ix_col + 0.5) * RES  # (NX,)
yc_of_row = Y0 + (iy_row + 0.5) * RES  # (NY,)

# flat-index lookups of per-cell fields for fast gather
Zg_flat = Zg_filled.ravel()
gx_flat = gx.ravel()
gy_flat = gy.ravel()
cos_flat = cos_slope.ravel()


def make_hist_array():
    return np.zeros((NY, NX, NBINS), dtype=np.uint32)


def accumulate(path, hist_all, hist_ground, fine_ground):
    """Stream a cloud, compute slope-normal d per return, bin into per-cell hists.

    Fills in place: hist_all, hist_ground (both (NY,NX,NBINS) uint32, coarse saved
    basis) and fine_ground ((NY,NX,NFBINS) uint32, fine internal ground histogram for
    precise summaries). class 7 (noise) is excluded. Ground = classification==2.
    """
    n_total = 0
    n_kept = 0
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(CHUNK):
            x = np.asarray(pts.x, dtype=np.float64)
            y = np.asarray(pts.y, dtype=np.float64)
            z = np.asarray(pts.z, dtype=np.float64)
            cl = np.asarray(pts.classification)
            n_total += x.size

            # cell indices
            ix = ((x - X0) / RES).astype(np.int64)
            iy = ((y - Y0) / RES).astype(np.int64)

            # keep in-grid and non-noise
            keep = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY) & (cl != 7)
            if not keep.any():
                continue
            ix = ix[keep]; iy = iy[keep]
            x = x[keep]; y = y[keep]; z = z[keep]; cl = cl[keep]
            n_kept += ix.size

            cell = iy * NX + ix  # flat cell index

            # predicted ground elevation: tilted plane extrapolated from cell centre
            xc = X0 + (ix + 0.5) * RES
            yc = Y0 + (iy + 0.5) * RES
            Zp = Zg_flat[cell] + gx_flat[cell] * (x - xc) + gy_flat[cell] * (y - yc)
            r = z - Zp
            d = r * cos_flat[cell]  # slope-normal distance (m)

            # bin index; clip so <BIN_LO and >=BIN_HI drop out (searchsorted -> 0..NBINS)
            b = np.searchsorted(edges, d, side="right") - 1
            valid = (b >= 0) & (b < NBINS)
            cellv = cell[valid]
            bv = b[valid]
            gv = (cl[valid] == 2)

            # flat 3D index into (NY,NX,NBINS): cell*NBINS + bin
            flat3 = cellv * NBINS + bv
            np.add.at(hist_all.reshape(-1), flat3, 1)
            if gv.any():
                np.add.at(hist_ground.reshape(-1), flat3[gv], 1)

            # fine ground-band histogram (class 2 only) for precise summaries
            gmask = (cl == 2)
            if gmask.any():
                dg = d[gmask]
                cellg = cell[gmask]
                fb = np.searchsorted(fedges, dg, side="right") - 1
                fv = (fb >= 0) & (fb < NFBINS)
                fflat = cellg[fv] * NFBINS + fb[fv]
                np.add.at(fine_ground.reshape(-1), fflat, 1)

    print(f"  {path}: total={n_total:,} kept(in-grid,non-noise)={n_kept:,}")


# ---- accumulate both epochs ----
def make_fine_array():
    return np.zeros((NY, NX, NFBINS), dtype=np.uint32)


print("gen2 (streaming 182.9M returns)...")
gen2_all = make_hist_array()
gen2_ground = make_hist_array()
gen2_fine = make_fine_array()
accumulate(GEN2, gen2_all, gen2_ground, gen2_fine)

print("gen1 (streaming)...")
gen1_all = make_hist_array()
gen1_ground = make_hist_array()
gen1_fine = make_fine_array()
accumulate(GEN1, gen1_all, gen1_ground, gen1_fine)


# ---- compact per-cell summaries from the histograms ----
centres = 0.5 * (edges[:-1] + edges[1:])  # (NBINS,)


def per_cell_quantile(hist, q, bin_centres):
    """Per-cell q-quantile of d from a (NY,NX,nb) count histogram.

    Returns the centre of the first bin whose cumulative count reaches q*n; NaN where
    empty. Bin-centre precision equals the histogram's bin width, so pass a FINE
    histogram (fedges/fcentres) where a precise value is needed.
    """
    n = hist.sum(axis=2)  # (NY,NX)
    out = np.full((NY, NX), np.nan, dtype=np.float32)
    cum = np.cumsum(hist, axis=2).astype(np.float64)  # (NY,NX,nb)
    tgt = q * n  # (NY,NX)
    nz = n > 0
    ge = cum >= tgt[:, :, None]
    idx = np.argmax(ge, axis=2)  # (NY,NX); 0 if none (but nz guarantees some)
    out[nz] = bin_centres[idx[nz]].astype(np.float32)
    return out


print("computing per-cell summaries...")
# ground median/p10 from the FINE (1 cm) histogram for a differenceable few-cm signal
gen2_ground_median_d = per_cell_quantile(gen2_fine, 0.50, fcentres)
gen2_ground_p10_d = per_cell_quantile(gen2_fine, 0.10, fcentres)
gen1_ground_median_d = per_cell_quantile(gen1_fine, 0.50, fcentres)
gen1_ground_p10_d = per_cell_quantile(gen1_fine, 0.10, fcentres)

# counts
gen2_n_all = gen2_all.sum(axis=2).astype(np.uint32)
gen2_n_ground = gen2_ground.sum(axis=2).astype(np.uint32)
gen1_n_all = gen1_all.sum(axis=2).astype(np.uint32)
gen1_n_ground = gen1_ground.sum(axis=2).astype(np.uint32)

# gen2 veg / canopy structure from gen2_all histogram
# d>0.5 : "above ground" fraction
b_half = int(np.searchsorted(edges, 0.5, side="right") - 1)   # bin containing 0.5 boundary
# use edges to identify counts strictly above 0.5 and within (0.5,2]
# fraction with d>0.5: sum of bins with centre>0.5 ... but bins straddle; use edge threshold.
# We bin by edges; count of returns with d>0.5 == returns in bins whose lower edge >= 0.5,
# plus partial of the straddling bin. For a fraction summary, use bin left-edge >= 0.5.
left_edges = edges[:-1]
mask_gt_half = left_edges >= 0.5              # bins entirely above 0.5
mask_under = (left_edges >= 0.5) & (edges[1:] <= 2.0)  # (0.5,2] understory bins

n2 = gen2_n_all.astype(np.float64)
with np.errstate(invalid="ignore", divide="ignore"):
    gen2_vegfrac = (gen2_all[:, :, mask_gt_half].sum(axis=2).astype(np.float64) / n2).astype(np.float32)
    gen2_understory_frac = (gen2_all[:, :, mask_under].sum(axis=2).astype(np.float64) / n2).astype(np.float32)
gen2_vegfrac[gen2_n_all == 0] = np.nan
gen2_understory_frac[gen2_n_all == 0] = np.nan

gen2_canopy_p95_d = per_cell_quantile(gen2_all, 0.95, centres)

# ---- save ----
print(f"saving {OUT} ...")
np.savez_compressed(
    OUT,
    edges=edges,
    gen2_all=gen2_all,
    gen2_ground=gen2_ground,
    gen1_all=gen1_all,
    gen1_ground=gen1_ground,
    gen2_ground_median_d=gen2_ground_median_d,
    gen2_ground_p10_d=gen2_ground_p10_d,
    gen1_ground_median_d=gen1_ground_median_d,
    gen1_ground_p10_d=gen1_ground_p10_d,
    gen2_vegfrac=gen2_vegfrac,
    gen2_understory_frac=gen2_understory_frac,
    gen2_canopy_p95_d=gen2_canopy_p95_d,
    gen2_n_all=gen2_n_all,
    gen2_n_ground=gen2_n_ground,
    gen1_n_all=gen1_n_all,
    gen1_n_ground=gen1_n_ground,
)
print("saved.")

# ---- sanity report ----
print("\n===== SANITY REPORT =====")
s = np.hypot(gx, gy)  # slope magnitude (m/m)
flat = s < 0.05  # <~2.9 deg
# "open" cells: gen2 veg fraction low -> few returns above ground
open_cell = (gen2_vegfrac < 0.05) & (gen2_n_ground > 5)
forest_cell = (gen2_vegfrac > 0.5) & (gen2_n_ground > 5)

openflat = open_cell & flat
print(f"open&flat cells: {openflat.sum()}")
print(f"  gen2 ground median-d on open&flat: {np.nanmedian(gen2_ground_median_d[openflat]):+.4f} m")
print(f"  gen1 ground median-d on open&flat: {np.nanmedian(gen1_ground_median_d[openflat & (gen1_n_ground > 5)]):+.4f} m")

# gen2-minus-gen1 ground median-d difference, forested vs open
both = (gen2_n_ground > 5) & (gen1_n_ground > 5)
diff = gen2_ground_median_d - gen1_ground_median_d
op = both & open_cell
fo = both & forest_cell
print(f"\ngen2 - gen1 ground median-d (slope-normal offset):")
print(f"  open cells   (n={op.sum():6d}): median {np.nanmedian(diff[op]):+.4f} m")
print(f"  forest cells (n={fo.sum():6d}): median {np.nanmedian(diff[fo]):+.4f} m")
print("=========================")
