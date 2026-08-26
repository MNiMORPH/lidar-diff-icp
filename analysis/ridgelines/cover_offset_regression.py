#!/usr/bin/env python3
"""Fit offset vs percent forest cover DETERMINISTICALLY: no bins, no thresholds, no seed.

Every earlier version of this answer went through binned medians, and each binning choice
(edges, minimum counts, a cover cut, a return-count filter) turned out to move the result --
the segmented-vs-smooth contest and the apparent high-cover acceleration were both decided
by choices rather than by data. This fits the returns directly.

Two lines, both closed-form or deterministically iterated, no tuning:

    OLS   minimises squared error -- the plainest possible fit
    LAD   minimises absolute error (the median-regression analogue of every binned median
          used so far), by IRLS from the OLS start, fixed tolerance and iteration cap

Uncertainty is a DELETE-ONE-BLOCK JACKKNIFE over square spatial blocks: refit with each
block removed, and take the spread of the refits. That accounts for returns inside one
woodlot not being independent, and unlike a bootstrap it uses no random numbers, so the
error bars are reproducible to the last digit. For OLS the refits are exact and O(1) per
block from sufficient statistics.

The only remaining choices are the physical ones -- which ground is declared not to have
changed -- and they are printed, not buried.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/cover_offset_regression.py --tile data/derived/elbaext
"""
import argparse, json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

from lidar_diff_icp import binstats as bs

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--offset", default="corr", choices=("raw", "corr"))
ap.add_argument("--inc-max", type=float, default=5.0)
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--block-m", type=float, default=50.0)
ap.add_argument("--per-cell", action="store_true",
                help="fit per-cell medians instead of raw returns (robustness check)")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = ("" if TILE == "elba_fulldensity" else f"_{TILE}") + ("_percell" if A.per_cell else "")
DCOL = "d_mm_corr" if A.offset == "corr" else "d_mm"

meta = next(json.load(open(f"{A.tile}/{fn}")) for fn in
            ("meta.json", "corrections_geoid.json", "corrections.json")
            if os.path.exists(f"{A.tile}/{fn}"))
res = float(meta.get("res") or meta.get("res_m"))
nx = int(meta.get("nx") or round((meta["bounds"][2]-meta["bounds"][0])/res))

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", DCOL, "incidence", "canopy_cover", "curv_laplacian", "in_grid"])
df = df[df.in_grid.values]
sel = (np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()[df.cell.to_numpy()]
       & (df.curv_laplacian.abs().to_numpy() <= A.curv_max)
       & (df.incidence.to_numpy() < A.inc_max)
       & np.isfinite(df[DCOL].to_numpy()) & np.isfinite(df.canopy_cover.to_numpy()))
s = df[sel]
cell = s.cell.to_numpy(); y = s[DCOL].to_numpy(float); x = s.canopy_cover.to_numpy(float)
blk = bs.block_ids(cell, nx=nx, res=res, block_m=A.block_m)
unit = "returns"
if A.per_cell:
    g = pd.DataFrame({"cell": cell, "y": y, "x": x, "b": blk}).groupby("cell")
    y = g.y.median().to_numpy(); x = g.x.first().to_numpy(); blk = g.b.first().to_numpy()
    unit = "cell medians"

print("=" * 90)
print(f"OFFSET vs COVER, fitted directly  [{TILE}, {DCOL}, {unit}]")
print(f"ground declared unchanged: divides AND |Laplacian|<={A.curv_max:g} AND incidence<{A.inc_max:g}°")
print(f"{len(y):,} {unit}; cover {x.min():.3f}–{x.max():.3f}; no bins, no cover cut, no seed")
print("=" * 90)


def ols(x_, y_):
    n = x_.size; sx = x_.sum(); sy = y_.sum(); sxx = (x_*x_).sum(); sxy = (x_*y_).sum()
    den = n*sxx - sx*sx
    b = (n*sxy - sx*sy) / den
    return (sy - b*sx) / n, b


def lad(x_, y_, iters=60, tol=1e-9):
    """Least absolute deviations by IRLS from the OLS start. Deterministic."""
    a, b = ols(x_, y_)
    for _ in range(iters):
        r = y_ - (a + b*x_)
        w = 1.0 / np.maximum(np.abs(r), 1e-6)
        sw = w.sum(); swx = (w*x_).sum(); swy = (w*y_).sum()
        swxx = (w*x_*x_).sum(); swxy = (w*x_*y_).sum()
        den = sw*swxx - swx*swx
        nb = (sw*swxy - swx*swy) / den
        na = (swy - nb*swx) / sw
        if abs(nb-b) < tol and abs(na-a) < tol:
            a, b = na, nb; break
        a, b = na, nb
    return a, b


ub, inv = np.unique(blk, return_inverse=True)
nblk = ub.size
# exact delete-one-block OLS from sufficient statistics
n_t = y.size; sx_t = x.sum(); sy_t = y.sum(); sxx_t = (x*x).sum(); sxy_t = (x*y).sum()
n_b = np.bincount(inv, minlength=nblk).astype(float)
sx_b = np.bincount(inv, weights=x, minlength=nblk)
sy_b = np.bincount(inv, weights=y, minlength=nblk)
sxx_b = np.bincount(inv, weights=x*x, minlength=nblk)
sxy_b = np.bincount(inv, weights=x*y, minlength=nblk)
n_j = n_t - n_b; sx_j = sx_t - sx_b; sy_j = sy_t - sy_b
sxx_j = sxx_t - sxx_b; sxy_j = sxy_t - sxy_b
den_j = n_j*sxx_j - sx_j**2
good = den_j > 0
b_j = (n_j[good]*sxy_j[good] - sx_j[good]*sy_j[good]) / den_j[good]
a_j = (sy_j[good] - b_j*sx_j[good]) / n_j[good]
k = b_j.size
jack = lambda t, full: np.sqrt((k-1)/k * np.sum((t - t.mean())**2))
a0, b0 = ols(x, y)
print(f"\nOLS   d = {a0:+.2f} {b0:+.2f} * cover   "
      f"(jackknife SE: intercept ±{jack(a_j,a0):.2f}, slope ±{jack(b_j,b0):.2f}; {k:,} blocks)")
al, bl = lad(x, y)
print(f"LAD   d = {al:+.2f} {bl:+.2f} * cover   (median regression; the binned-median analogue)")
print(f"\n  slope in mm per 10% cover:  OLS {b0/10:+.2f}   LAD {bl/10:+.2f}")
print(f"  predicted offset:  " + "  ".join(f"{c:.0%}: OLS {a0+b0*c:+.0f} / LAD {al+bl*c:+.0f} mm"
                                           for c in (0.1, 0.3, 0.5)))

fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=130)
# min_n=1 is binned_stats' definitional floor. These check bins are equipopulated by
# construction (quantile edges), so a minimum count removes nothing here today -- which is
# exactly why carrying one is dead weight that would start deleting the sparse tail the
# moment the tile, the bin count or the selection changed.
ck = bs.binned_stats(x, y, bs.quantile_edges(x, 12, first_edge=0.02), block=blk, min_n=1)
ax.errorbar(ck.x, ck.y, yerr=ck.se, fmt="o", ms=5, capsize=3, color="0.35", zorder=4,
            label="binned medians ± block SE (check only, not fitted)")
xs = np.linspace(0, x.max(), 200)
ax.plot(xs, a0 + b0*xs, "-", lw=2.0, color="C0", label=f"OLS  {a0:+.1f} {b0:+.1f}·cover")
ax.plot(xs, al + bl*xs, "-", lw=2.0, color="C3", label=f"LAD  {al:+.1f} {bl:+.1f}·cover")
ax.axhline(0, color="k", lw=.7)
ax.set_xlabel("canopy cover (fraction)")
ax.set_ylabel("offset d (mm)   [gen1 − gen2; + = lower in 2021]")
ax.set_title(f"offset vs cover, fitted to {unit} directly — {TILE}\n"
             f"divides, |Laplacian|≤{A.curv_max:g}, incidence<{A.inc_max:g}°; "
             f"delete-one-block jackknife ({A.block_m:g} m)", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=.3)
out = f"figures/refdatum/cover_offset_regression{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
