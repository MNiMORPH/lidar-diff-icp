#!/usr/bin/env python3
"""The offset-vs-cover relationship, shape-free: isotonic regression.

Every parametric attempt at this curve forced a decision the data should have made -- a
break point, a bin edge, a cover cut, a functional form -- and each decision moved the
answer. Isotonic regression removes all of them. It has no bins, no bandwidth, no
thresholds, no random seed, and exactly one assumption: that increasing canopy cover does
not make the offset less negative. That is physically expected (more canopy, more of the
2021 ground surface sitting in vegetation) and is what the data show.

It is also the right tool for the concern a straight line handles badly: dense-forest cells
are rare, so an abundance-weighted line is pulled to the sparse-forest regime and
under-predicts the dense end. Isotonic fits where the data are, at whatever level they sit,
without a rare regime having to out-vote a common one to be represented.

Fitted to PER-CELL medians so that within-cell return noise is already reduced and no cell
counts more than another for being over-flown more often. Uncertainty is a delete-one-block
jackknife over 50 m blocks -- deterministic, and it accounts for cells in one woodlot not
being independent.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/cover_offset_isotonic.py --tile data/derived/elbaext
"""
import argparse, json, os
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

from lidar_diff_icp import binstats as bs

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--offset", default="corr", choices=("raw", "corr"))
ap.add_argument("--inc-max", type=float, default=5.0)
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--block-m", type=float, default=50.0)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = "" if TILE == "elba_fulldensity" else f"_{TILE}"
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
g = pd.DataFrame({"cell": s.cell.to_numpy(), "y": s[DCOL].to_numpy(float),
                  "x": s.canopy_cover.to_numpy(float)}).groupby("cell")
y = g.y.median().to_numpy(); x = g.x.first().to_numpy(); cells = np.array(list(g.groups.keys()))
blk = bs.block_ids(cells, nx=nx, res=res, block_m=A.block_m)
print("=" * 88)
print(f"ISOTONIC OFFSET vs COVER  [{TILE}, {DCOL}]")
print(f"ground declared unchanged: divides AND |Laplacian|<={A.curv_max:g} AND incidence<{A.inc_max:g}°")
print(f"{len(y):,} cell medians ({sel.sum():,} returns); cover {x.min():.3f}–{x.max():.3f}")
print("no bins, no bandwidth, no cover cut, no seed; only monotonicity assumed")
print("=" * 88)

fit = lambda xx, yy: IsotonicRegression(increasing=False, out_of_bounds="clip").fit(xx, yy)
iso = fit(x, y)
GRID = np.array([0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80])
curve = iso.predict(GRID)

ub = np.unique(blk)
jk = np.empty((ub.size, GRID.size))
for i, u in enumerate(ub):
    k = blk != u
    jk[i] = fit(x[k], y[k]).predict(GRID)
kk = ub.size
se = np.sqrt((kk - 1) / kk * np.sum((jk - jk.mean(axis=0)) ** 2, axis=0))

print(f"\n{'cover':>7s} {'offset (mm)':>13s} {'jackknife SE':>14s} {'cells at or above':>19s}")
for c, v, e in zip(GRID, curve, se):
    print(f"{c:>7.2f} {v:>13.1f} {e:>14.1f} {int((x >= c).sum()):>19,}")

# what a straight line would have said, for contrast
def ols(xx, yy):
    n = xx.size; sx = xx.sum(); sy = yy.sum(); sxx = (xx*xx).sum(); sxy = (xx*yy).sum()
    b = (n*sxy - sx*sy) / (n*sxx - sx*sx); return (sy - b*sx)/n, b
a0, b0 = ols(x, y)
print(f"\n  a straight line through the same cells: d = {a0:+.1f} {b0:+.1f}*cover")
print(f"  {'cover':>7s} {'isotonic':>10s} {'line':>8s} {'line error':>12s}")
for c in (0.10, 0.30, 0.50, 0.70):
    v = float(iso.predict([c])[0]); l = a0 + b0*c
    print(f"  {c:>7.2f} {v:>10.1f} {l:>8.1f} {l-v:>+12.1f}")

# ---- FIGURE: show every observation, and resolve the high-cover end explicitly ----
# quantile bins collapse everything above ~0.39 into one point plotted at its MEAN cover,
# which hides the very data the curve is about. Bins here are quantile-spaced only up to
# 0.45 and FIXED above it, so each high-cover interval is drawn where it actually sits.
fig, ax = plt.subplots(2, 1, figsize=(9.2, 8.0), dpi=130, sharex=True,
                       gridspec_kw={"height_ratios": [3, 1]})
a0x = ax[0]
hb = a0x.hexbin(x, np.clip(y, -400, 400), gridsize=(70, 45), bins="log", cmap="Greys",
                mincnt=1, extent=(0, x.max(), -400, 400))
fig.colorbar(hb, ax=a0x, label="log10 cells", pad=0.01)
lowe = bs.quantile_edges(x[x <= 0.45], 10, first_edge=0.02)
edges = np.unique(np.concatenate([lowe[lowe <= 0.45], [0.45, 0.55, 0.65, 0.75, 0.85, 1.0]]))
ck = bs.binned_stats(x, y, edges, block=blk, min_n=5)
a0x.errorbar(ck.x, ck.y, yerr=ck.se, xerr=[ck.x-ck.lo, ck.hi-ck.x], fmt="o", ms=5, capsize=2,
             color="C3", lw=1.3, zorder=6, label="binned medians ± block SE (bars = bin span)")
xs = np.linspace(0, x.max(), 600)
a0x.plot(xs, iso.predict(xs), "-", lw=2.4, color="C2", zorder=7, label="isotonic (shape-free)")
a0x.fill_between(GRID, curve-se, curve+se, color="C2", alpha=.25, zorder=5)
a0x.plot(xs, a0 + b0*xs, "--", lw=1.4, color="C0", zorder=4,
         label=f"straight line {a0:+.0f}{b0:+.0f}·cover")
a0x.axhline(0, color="k", lw=.7)
a0x.set_ylim(-400, 400); a0x.set_xlim(0, x.max()*1.01)
a0x.set_ylabel("offset d (mm)   [gen1 − gen2; + = lower in 2021]")
a0x.set_title(f"offset vs cover, every cell shown — {TILE}\n"
              f"divides, |Laplacian|≤{A.curv_max:g}, incidence<{A.inc_max:g}°, {len(y):,} cells",
              fontsize=10)
a0x.legend(fontsize=8, loc="lower left"); a0x.grid(alpha=.3)
for lo_, hi_, n_ in zip(ck.lo, ck.hi, ck.n):
    if lo_ >= 0.45:
        a0x.annotate(f"n={n_}", (0.5*(lo_+hi_), -370), ha="center", fontsize=6.5, color="C3")
ax[1].bar(0.5*(ck.lo+ck.hi), ck.n, width=(ck.hi-ck.lo)*0.9, color="0.6")
ax[1].set_yscale("log"); ax[1].set_ylabel("cells per bin"); ax[1].set_xlabel("canopy cover (fraction)")
ax[1].grid(alpha=.3)
out = f"figures/refdatum/cover_offset_isotonic{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
print(f"\nbins actually drawn above cover 0.45:")
for lo_, hi_, xx_, yy_, se_, n_ in zip(ck.lo, ck.hi, ck.x, ck.y, ck.se, ck.n):
    if lo_ >= 0.45:
        print(f"   {lo_:.2f}-{hi_:.2f}: median {yy_:+8.1f} ± {se_:5.1f} mm   n={n_:,} cells")
