#!/usr/bin/env python3
"""CLI over :mod:`lidar_diff_icp.offset_model`: predict the median gen1-vs-gen2 offset
jointly from surface SLOPE and CANOPY COVER on one tile, with the covariance between
those predictors kept in view.

Basis: the 2026-08-24 slope x cover disentangling (originally run inline) -- per-cell
medians of the per-beam offset, near-planar restriction, additive vs interaction fits,
partial correlations, and matched-band tables in both directions. The reusable maths now
lives in the package; this script is the tile-level driver, report, and figure.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/offset_model_slope_cover.py --tile data/derived/elbaext \
        --curv-max 0.015 --ridge
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

from lidar_diff_icp.offset_model import (fit_offset_model, matched_band_effects,
                                         median_surface, partial_correlations,
                                         predictor_covariance, surface_centres)

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--curv-max", type=float, default=0.015, help="|Laplacian| cut (near-planar cells)")
ap.add_argument("--ridge", action="store_true", help="restrict to ridge_mask cells (divides: no overland flow)")
ap.add_argument("--min-n", type=int, default=3, help="minimum returns per cell for a cell median")
ap.add_argument("--min-cells", type=int, default=30, help="cells needed before a grid box is reported")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
# the curvature cut MUST be in the name: different cuts are different analyses and
# gave materially different coefficients, so one filename for all of them loses work.
TAG = (("" if TILE == "elba_fulldensity" else f"_{TILE}")
       + f"_curv{A.curv_max:g}" + ("_ridge" if A.ridge else ""))

SLOPE_EDGES = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 45], float)
COVER_EDGES = np.array([0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 1.01])

# ---------------------------------------------------------------- data: per-cell medians
df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", "d_mm", "slope", "canopy_cover", "curv_laplacian", "in_grid"])
sel = df.in_grid.to_numpy(bool) & np.isfinite(df.d_mm) & np.isfinite(df.slope) & np.isfinite(df.canopy_cover)
sel &= np.isfinite(df.curv_laplacian) & (df.curv_laplacian.abs().to_numpy() <= A.curv_max)
if A.ridge:
    rm = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
    sel &= rm[df.cell.to_numpy()]
d = df[sel]
g = d.groupby("cell")
cell_n = g.size(); keep = cell_n >= A.min_n
med = g.d_mm.median()[keep].to_numpy(float)
slp = g.slope.first()[keep].to_numpy(float)          # slope/cover are per-cell constants
cov = g.canopy_cover.first()[keep].to_numpy(float)
n = med.size
lab = f"{TILE}{' RIDGELINES' if A.ridge else ''}, |Laplacian|<={A.curv_max}, cells with >={A.min_n} returns"
print("=" * 92)
print(f"MEDIAN OFFSET as f(slope, canopy cover)   [{lab}]")
print(f"{len(d):,} returns -> {n:,} cell medians")
print("=" * 92)

# ---------------------------------------------------------------- 1. the covariance itself
r_sc, vif = predictor_covariance(slp, cov)
print(f"\n1. COVARIANCE OF THE PREDICTORS")
print(f"   corr(slope, cover) = {r_sc:+.3f}    VIF = {vif:.2f}   (1 = orthogonal)")
print(f"   slope {slp.min():.1f}-{slp.max():.1f} deg (median {np.median(slp):.1f});"
      f"  cover {cov.min():.2f}-{cov.max():.2f} (median {np.median(cov):.2f})")

# ---------------------------------------------------------------- 2. per-cell fits
sz = (slp - slp.mean()) / slp.std(); cz = (cov - cov.mean()) / cov.std()
za = fit_offset_model(sz, cz, med, interaction=False)
zi = fit_offset_model(sz, cz, med, interaction=True)
pa = fit_offset_model(slp, cov, med, interaction=False)
pi = fit_offset_model(slp, cov, med, interaction=True)
print(f"\n2. PER-CELL MODEL (response = per-cell median offset, mm)")
print(f"   z-scored, additive    : d = {za.coeffs[0]:+7.1f} {za.coeffs[1]:+7.1f}*slope_z "
      f"{za.coeffs[2]:+7.1f}*cover_z                       R2 = {za.r2:.4f}")
print(f"   z-scored, +interaction: d = {zi.coeffs[0]:+7.1f} {zi.coeffs[1]:+7.1f}*slope_z "
      f"{zi.coeffs[2]:+7.1f}*cover_z {zi.coeffs[3]:+7.1f}*(slope x cover)_z   R2 = {zi.r2:.4f}")
print(f"   physical units        : {pi}")
print(f"   NB R2 at per-cell scale is ~1e-3 by construction: per-cell scatter is ~100x the")
print(f"   systematic term. The median surface (5b) is the scale on which prediction lives.")

# ---------------------------------------------------------------- 3. partial correlations
print(f"\n3. PARTIAL CORRELATION (each predictor NET of the other)")
for nm, (rp, rs) in partial_correlations(slp, cov, med).items():
    print(f"   d ~ {nm:14s} Pearson {rp:+.3f}   Spearman {rs:+.3f}")

# ---------------------------------------------------------------- 4. matched bands
print(f"\n4a. COVER effect with SLOPE held in a narrow band  (d d/d cover, mm per full cover unit)")
for lo, hi, nb, grad, vmin, vmax in matched_band_effects(slp, cov, med, SLOPE_EDGES):
    if not np.isfinite(grad): print(f"    {lo:4.0f}-{hi:<4.0f} deg: n={nb:>6,}  (sparse)"); continue
    print(f"    {lo:4.0f}-{hi:<4.0f} deg: {grad:+7.0f} mm/unit   (n={nb:>6,}, cover {vmin:.2f}-{vmax:.2f})")
print(f"\n4b. SLOPE effect with COVER held in a narrow band  (d d/d slope, mm/deg)")
for lo, hi, nb, grad, vmin, vmax in matched_band_effects(cov, slp, med, COVER_EDGES):
    if not np.isfinite(grad): print(f"    cover {lo:.2f}-{hi:.2f}: n={nb:>6,}  (sparse)"); continue
    print(f"    cover {lo:.2f}-{hi:.2f}: {grad:+7.2f} mm/deg   (n={nb:>6,}, slope {vmin:.0f}-{vmax:.0f} deg)")

# ---------------------------------------------------------------- 5. the 2-D surface
grid, cnt = median_surface(slp, cov, med, SLOPE_EDGES, COVER_EDGES, min_cells=A.min_cells)
NS, NC = grid.shape
print(f"\n5. EMPIRICAL MEDIAN OFFSET (mm) on the (slope x cover) grid -- the prediction table.")
print(f"   '.' = fewer than {A.min_cells} cells: the terrain does not supply that combination,")
print(f"   so it is UNSUPPORTED, not extrapolated.  n cells in parentheses.")
print("   slope\\cover " + "".join(f"{COVER_EDGES[j]:.2f}-{COVER_EDGES[j+1]:.2f}".rjust(14) for j in range(NC)))
for i in range(NS):
    row = f"   {SLOPE_EDGES[i]:4.0f}-{SLOPE_EDGES[i+1]:<4.0f}   "
    for j in range(NC):
        row += (f"{grid[i,j]:+7.0f}({cnt[i,j]:>5,})" if np.isfinite(grid[i, j])
                else f"{'.':>7}({cnt[i,j]:>5,})").rjust(14)
    print(row)
ok = np.isfinite(grid)
print(f"   supported boxes: {ok.sum()}/{grid.size} ({100*ok.sum()/grid.size:.0f}%) "
      f"-- the empty ones ARE the slope-cover covariance.")

SC, CC = surface_centres(SLOPE_EDGES, COVER_EDGES)
print(f"   per-cell interaction model vs this table: "
      f"RMS {np.sqrt(np.mean((grid[ok]-pi.predict(SC[ok], CC[ok]))**2)):.1f} mm, "
      f"max |resid| {np.max(np.abs(grid[ok]-pi.predict(SC[ok], CC[ok]))):.1f} mm")

# ---------------------------------------------------------------- 5b. the median-surface fit
w = cnt[ok].astype(float)
ma = fit_offset_model(SC[ok], CC[ok], grid[ok], weights=w, interaction=False)
mi = fit_offset_model(SC[ok], CC[ok], grid[ok], weights=w, interaction=True)
print(f"\n5b. MEDIAN-SURFACE FIT (n-weighted over the {ok.sum()} supported boxes) -- THE PREDICTOR")
print(f"   additive    : {ma}")
print(f"   +interaction: {mi}")
for s_at in (5, 25):
    print(f"   d(offset)/d(cover) at slope {s_at:2d} deg: {mi.d_dcover(s_at):+.0f} mm/unit")
for c_at in (0.02, 0.60):
    print(f"   d(offset)/d(slope) at cover {c_at:.2f} : {mi.d_dslope(c_at):+.2f} mm/deg")

# ---------------------------------------------------------------- 6. figures
# PRIMARY figure follows the offset_vs_angle house style: per-beam density with a binned
# median on the left, medians split by canopy cover on the right -- here at the model's
# full cover resolution, with the fitted surface overlaid so model and data are read
# against each other in one place.
sl_ret = d.slope.to_numpy(float); d_ret = d.d_mm.to_numpy(float); cc_ret = d.canopy_cover.to_numpy(float)
XLIM = (0, 45)
fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.6), dpi=130)

hb = ax[0].hexbin(sl_ret, d_ret, gridsize=(60, 60), bins="log", cmap="viridis",
                  extent=(XLIM[0], XLIM[1], -300, 300), mincnt=1)
fig.colorbar(hb, ax=ax[0], label="log10 count")
be = np.arange(XLIM[0], XLIM[1] + 1.0, 1.0)
bc, bm = [], []
for lo, hi in zip(be[:-1], be[1:]):
    m = (sl_ret >= lo) & (sl_ret < hi)
    if m.sum() >= 50:
        bc.append(0.5*(lo+hi)); bm.append(np.median(d_ret[m]))
ax[0].plot(bc, bm, "o-", color="crimson", ms=4, lw=1.4, label="binned median")
ax[0].axhline(0, color="k", lw=.6); ax[0].set_xlim(*XLIM); ax[0].set_ylim(-300, 300)
ax[0].set_xlabel("surface slope (deg)")
ax[0].set_ylabel("offset d (mm) = gen1 − gen2   (+ = ground lower in 2021)")
ax[0].set_title("PRIMARY: per-beam offset vs slope (all returns)"); ax[0].legend()

cen_s = 0.5*(SLOPE_EDGES[:-1]+SLOPE_EDGES[1:])
allc, allm = [], []
for lo, hi in zip(SLOPE_EDGES[:-1], SLOPE_EDGES[1:]):
    m = (sl_ret >= lo) & (sl_ret < hi)
    if m.sum() >= 50: allc.append(0.5*(lo+hi)); allm.append(np.median(d_ret[m]))
ax[1].plot(allc, allm, "o-", color="0.35", lw=2.0, ms=5, label="all returns", zorder=5)
for j in range(NC):
    c_c = 0.5*(COVER_EDGES[j]+COVER_EDGES[j+1]); good = np.isfinite(grid[:, j])
    if good.sum() < 2: continue
    line, = ax[1].plot(cen_s[good], grid[good, j], "o-", ms=4, lw=1.3,
                       label=f"cover {COVER_EDGES[j]:.2f}-{COVER_EDGES[j+1]:.2f} (n={cnt[good, j].sum():,})")
    ax[1].plot(cen_s[good], mi.predict(cen_s[good], c_c), "--", lw=1.0, color=line.get_color(), alpha=.75)
ax[1].axhline(0, color="k", lw=.6); ax[1].set_xlim(*XLIM)
ax[1].set_xlabel("surface slope (deg)"); ax[1].set_ylabel("median offset d (mm)   [gen1 − gen2; + = lower in 2021]")
ax[1].set_title("median offset vs slope, by canopy cover\n(solid = observed, dashed = joint model)")
ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
fig.suptitle(f"gen1 offset vs slope and canopy cover — {lab}", y=1.0)
out = f"figures/refdatum/offset_model_slope_cover{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")

# SECONDARY: the (slope x cover) support grid -- where the terrain does and does not sample
fig2, bx = plt.subplots(figsize=(7.6, 5.6), dpi=130)
vm = np.nanmax(np.abs(grid))
im = bx.imshow(grid, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vm, vmax=vm,
               extent=[0, NC, 0, NS])
for i in range(NS):
    for j in range(NC):
        bx.text(j+0.5, i+0.5, f"{grid[i,j]:+.0f}" if np.isfinite(grid[i,j]) else "·",
                ha="center", va="center", fontsize=7.5,
                color="k" if np.isfinite(grid[i,j]) else "0.55")
bx.set_xticks(np.arange(NC)+0.5)
bx.set_xticklabels([f"{COVER_EDGES[j]:.2f}-{COVER_EDGES[j+1]:.2f}" for j in range(NC)], fontsize=7, rotation=30)
bx.set_yticks(np.arange(NS)+0.5)
bx.set_yticklabels([f"{SLOPE_EDGES[i]:.0f}-{SLOPE_EDGES[i+1]:.0f}" for i in range(NS)], fontsize=7)
bx.set_xlabel("canopy cover (PyForestScan fraction)"); bx.set_ylabel("surface slope (deg)")
bx.set_title(f"median offset (mm) on the slope x cover grid — {lab}\n'·' = unsupported by the terrain", fontsize=9)
fig2.colorbar(im, ax=bx, label="median offset d (mm)   [gen1 − gen2; + = lower in 2021]")
out2 = f"figures/refdatum/offset_model_grid{TAG}.png"
fig2.savefig(out2, bbox_inches="tight"); plt.close(fig2)
print(f"wrote {out2}")
