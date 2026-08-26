#!/usr/bin/env python3
"""Continuous (float) percentile fit for matching gen1 and gen2 ground, both directions.

Two ways to make the DoD unbiased, fitted as REAL numbers rather than on a 0.01 grid:
  A. raise gen1  -- find q1 such that median_cells[ Q1(q1) - z_after ] = 0
  B. lower gen2  -- hold gen1 at q = 0.50 and find q2 such that
                    median_cells[ Q1(0.50) - Q2(q2) ] = 0

Q1 is a per-cell linear-interpolated quantile of the REGISTERED per-return offsets
(`d_mm_corr`, gen1 CSF ground). Q2 is the same on gen2's class-2 near-ground histogram
(2 cm bins, interpolated within the bin so Q2 is continuous in q). Both are monotone in q,
so the pooled median residual is monotone and Brent's method finds the root exactly.

    ./lidar-icp/bin/python analysis/ridgelines/percentile_float_fit.py
"""
import argparse
import numpy as np, pyarrow.parquet as pq
from scipy.optimize import brentq

from lidar_diff_icp.binstats import block_ids
from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
# Per-cell sample-size requirements, DEFINITIONAL at 1 and discretionary above it: Q1 is a
# per-cell quantile of the gen1 offsets and needs >=1 return; Q2 is a quantile of the gen2
# near-ground histogram and needs >=1 point. Higher values re-select toward open ground,
# because the cells they drop are the canopy-enriched ones. Kept identical to
# q2_cover_fit.py's, which reads the same cube: two different cuts on one cube make the two
# results incomparable.
ap.add_argument("--min-gen1", type=int, default=1,
                help="gen1 returns needed per cell (1 = definitional: a quantile needs a point)")
ap.add_argument("--min-gen2", type=int, default=1,
                help="gen2 class-2 near-ground points needed per cell (1 = definitional)")
A = ap.parse_args()

# Fine cover bins, INCLUDING the sparse high-cover tail. Pooling ">0.35" hides it:
# 92% of that bin is 0.35-0.50, so the cells with the largest effect vanish into the
# average. Rarity is a fact about the landscape, not about the measurement -- n enters
# once, through the (cluster-robust) standard error.
STRATA = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
          ("mid   .20-.35", 0.20, 0.35), ("0.35-0.50", 0.35, 0.50),
          ("0.50-0.65", 0.50, 0.65), (">0.65", 0.65, 1.01)]


def ragged_sorted(cell, val, ncell):
    """Flat array of per-cell sorted values, with per-cell offsets and counts."""
    o = np.lexsort((val, cell))
    cs = cell[o]; vs = val[o]
    n = np.bincount(cs, minlength=ncell)
    off = np.r_[0, np.cumsum(n)[:-1]]
    return vs, off, n


def ragged_quantile(vs, off, n, q, sel):
    """Linear-interpolated quantile (numpy 'linear' convention) for the selected cells."""
    nn = n[sel].astype(float)
    pos = q * (nn - 1.0)
    lo = np.floor(pos).astype(np.int64); hi = np.minimum(lo + 1, nn.astype(np.int64) - 1)
    f = pos - lo
    b = off[sel]
    return vs[b + lo] * (1 - f) + vs[b + hi] * f


def hist_quantile(C, ntot, q, zlo, dz):
    """Quantile of a histogram row, interpolated WITHIN the bin, in mm."""
    r = q * ntot
    k = (C >= r[:, None]).argmax(1)
    below = np.where(k > 0, C[np.arange(C.shape[0]), np.maximum(k - 1, 0)], 0.0)
    inbin = C[np.arange(C.shape[0]), k] - below
    f = np.where(inbin > 0, (r - below) / np.maximum(inbin, 1e-9), 0.0)
    return (zlo + (k + np.clip(f, 0, 1)) * dz) * 1000.0


D = A.tile
cube = np.load(f"{D}/nearground_cells_sn.npz")
cells = cube["cells"]; dz = float(cube["dz"]); zlo = float(cube["zlo"])
zf = np.load(f"{D}/z_after.npy"); N = zf.size
cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()[cells]

t = pq.read_table(f"{D}/beam_offset_table.parquet",
                  columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
ce = t["cell"].to_numpy()[g]; dc = t["d_mm_corr"].to_numpy()[g].astype(float)
vs, off, n1 = ragged_sorted(ce, dc, N)

split = np.load(f"{D}/nearground_gen2_class_split.npz")
assert np.array_equal(split["cells"], cells), "class-split cube is on different cells"
Hg = split["Hg"]; Cg = np.cumsum(Hg, 1).astype(float); ng = Cg[:, -1]

stable, _ = reference_cells(D, cells=cells, slope_max=90.0)
ok = stable & (n1[cells] >= A.min_gen1) & (ng >= A.min_gen2) & np.isfinite(cover)
sel = cells[ok]
print(f"{D}: {ok.sum():,} stable cells "
      f"(gen1 >= {A.min_gen1} returns, gen2 class-2 >= {A.min_gen2})")

Cs = Cg[ok]; ns = ng[ok]


def bias_gen1(q):                      # A: raise gen1, target = z_after (h = 0)
    return float(np.median(ragged_quantile(vs, off, n1, q, sel)))


g1_med = ragged_quantile(vs, off, n1, 0.50, sel)


def bias_gen2(q):                      # B: lower gen2 to meet gen1's median
    return float(np.median(g1_med - hist_quantile(Cs, ns, q, zlo, dz)))


q1 = brentq(bias_gen1, 0.35, 0.75, xtol=1e-6)
q2 = brentq(bias_gen2, 0.20, 0.70, xtol=1e-6)
print(f"\nA. raise gen1 : q1* = {q1:.4f}   (pooled median residual "
      f"{bias_gen1(q1):+.3f} mm)")
print(f"B. lower gen2 : q2* = {q2:.4f}   (pooled median residual "
      f"{bias_gen2(q2):+.3f} mm)")

rA = ragged_quantile(vs, off, n1, q1, sel)
rB = g1_med - hist_quantile(Cs, ns, q2, zlo, dz)
blk = block_ids(sel, zf.shape[1], 5.0, 50.0)


def med_se(v, b):
    """Cluster-robust SE of a median: spread of 50 m block medians."""
    ub, inv = np.unique(b, return_inverse=True)
    bm = np.array([np.median(v[inv == i]) for i in range(ub.size)])
    return float(np.std(bm, ddof=1) / np.sqrt(ub.size)), ub.size


print(f"\nresidual by cover at the float optima (mm, 0 = perfect; +- cluster-robust SE)")
print(f"  {'stratum':14s} {'cells':>7s} {'blocks':>6s} | {'A: gen1 q1*':>16s} "
      f"{'B: gen2 q2*':>16s} | {'A q1*':>7s} {'B q2*':>7s}")
cv = cover[ok]
for nm, lo, hi in STRATA + [("ALL", -9, 9)]:
    m = (cv > lo) & (cv <= hi) if nm != "ALL" else np.ones(cv.size, bool)
    if not m.sum():
        continue
    qa = brentq(lambda q: float(np.median(ragged_quantile(vs, off, n1, q, sel[m]))),
                0.05, 0.95, xtol=1e-6)
    qb = brentq(lambda q: float(np.median(g1_med[m] - hist_quantile(Cs[m], ns[m], q, zlo, dz))),
                0.02, 0.95, xtol=1e-6)
    sa, nb = med_se(rA[m], blk[m]); sb, _ = med_se(rB[m], blk[m])
    print(f"  {nm:14s} {m.sum():7,d} {nb:6,d} | {np.median(rA[m]):8.1f} +- {sa:5.1f} "
          f"{np.median(rB[m]):8.1f} +- {sb:5.1f} | {qa:7.4f} {qb:7.4f}")
np.savez(f"{D}/percentile_float_fit.npz", q1=q1, q2=q2, cells=cells, ok=ok)
print(f"\nsaved {D}/percentile_float_fit.npz")
