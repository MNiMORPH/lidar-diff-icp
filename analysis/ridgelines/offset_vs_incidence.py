#!/usr/bin/env python3
"""PRIMARY PRODUCT: gen1 per-beam elevation offset vs gen2 (d_mm) as a function of the
beam INCIDENCE angle to the local surface -- the physically-correct axis that unifies
slope + aspect + scan angle.

Reads the canonical per-beam table (beam_offset_table.parquet). Reports, in this order:
  1. PRIMARY  -- d_mm vs incidence, all in-grid returns: robust binned median + NMAD + n,
                 plus a hexbin figure. This is the headline relationship.
  2. CLEAN    -- within-cell residual (d_resid_mm) vs incidence: per-cell effects (real
                 change, datum, canopy height) differenced out; leverage lives in
                 flight-line-overlap cells, so we also report the overlap-weighted slope.
  3. STRATIFY -- the primary split by canopy cover, to see the leaf-state offset that sits
                 on top of the geometry term.
  4. CHECK    -- per-cell scatter vs within-cell incidence spread: a falsifiable prediction
                 (scatter should grow with incidence spread if incidence drives the offset).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/offset_vs_incidence.py [tile_dir] [curv_max]

Optional curv_max restricts to near-planar cells (|curv_laplacian| <= curv_max) to
suppress hillslope-diffusion / convex-concave real change; the figure form is identical,
only the point subset changes (output filename and title carry the threshold).

No fitting is imposed beyond robust bin statistics and simple correlations.
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr

TILE = sys.argv[1] if len(sys.argv) > 1 else "data/derived/elba_fulldensity"
CURV_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else None
INC_EDGES = np.arange(0, 48, 3)          # incidence bins: 3 deg, 0-45 (held fixed for this axis)
MIN_N = 300                              # drop bins below this (matches prior slope analysis)

df = pd.read_parquet(f"{TILE}/beam_offset_table.parquet")
df = df[df.in_grid.values].copy()
lab, suffix = "all curvatures", ""
if CURV_MAX is not None:
    keep = (df.curv_laplacian.abs() <= CURV_MAX).to_numpy()
    frac = 100 * keep.mean(); df = df[keep].copy()
    lab = f"|Laplacian| <= {CURV_MAX:g}  ({len(df):,} returns, {frac:.0f}% kept)"
    suffix = f"_curv{CURV_MAX:g}"
inc = df.incidence.to_numpy(float); d = df.d_mm.to_numpy(float)
print(f"{len(df):,} returns  [{lab}]\n")

def nmad(x): return 1.4826 * np.median(np.abs(x - np.median(x)))

def binned(xv, yv, edges=INC_EDGES, min_n=MIN_N):
    """robust median + NMAD + n per bin; returns (centers, med, nmad, n)."""
    c, m, s, k = [], [], [], []
    for i in range(len(edges) - 1):
        b = (xv >= edges[i]) & (xv < edges[i + 1]) & np.isfinite(yv)
        if b.sum() < min_n: continue
        c.append((edges[i] + edges[i + 1]) / 2); m.append(np.median(yv[b]))
        s.append(nmad(yv[b])); k.append(int(b.sum()))
    return np.array(c), np.array(m), np.array(s), k

def table(title, xv, yv):
    c, m, s, k = binned(xv, yv)
    print(title)
    print(f"  {'incid(deg)':>10s} {'median':>9s} {'NMAD':>8s} {'n':>12s}")
    for a, mm, ss, nn in zip(c, m, s, k):
        print(f"  {a:10.0f} {mm:+9.1f} {ss:8.1f} {nn:12,d}")
    fin = np.isfinite(yv)
    pr = np.corrcoef(xv[fin], yv[fin])[0, 1]; sp = spearmanr(xv[fin], yv[fin]).statistic
    print(f"  corr(offset, incidence): Pearson {pr:+.3f}  Spearman {sp:+.3f}\n")
    return c, m, s

# ---------------------------------------------------------------- 1. PRIMARY
c, m, s = table("1. PRIMARY  d_mm vs incidence (all in-grid returns):", inc, d)

# ---------------------------------------------------------------- 2. CLEAN (within-cell)
dr = df.d_resid_mm.to_numpy(float)
table("2. CLEAN  within-cell residual d_resid_mm vs incidence (per-cell effects removed):", inc, dr)
# overlap-powered slope: regress residual on within-cell incidence deviation, weighting by
# how much incidence actually varies in the cell (cell_inc_std). Cells with ~0 spread carry ~0 info.
cim = df.groupby("cell")["incidence"].transform("mean").to_numpy(float)
inc_dev = inc - cim                      # within-cell incidence deviation
w = df.cell_inc_std.to_numpy(float); ok = np.isfinite(w) & np.isfinite(inc_dev) & np.isfinite(dr) & (w > 0)
W = w[ok]
slope_wc = np.sum(W * inc_dev[ok] * dr[ok]) / np.sum(W * inc_dev[ok] ** 2)
print(f"   within-cell slope d(offset)/d(incidence) = {slope_wc:+.2f} mm/deg "
      f"(overlap-weighted, {ok.sum():,} returns; leverage from high cell_inc_std cells)\n")

# ---------------------------------------------------------------- 3. STRATIFY by canopy cover
cc = df.canopy_cover.to_numpy(float)
open_m = np.isfinite(cc) & (cc < 0.10)   # provisional (forest/open threshold not yet calibrated)
for_m  = np.isfinite(cc) & (cc > 0.50)
print(f"3. STRATIFY by canopy cover (PROVISIONAL: open cc<0.10 n={open_m.sum():,}, "
      f"forest cc>0.50 n={for_m.sum():,}; threshold calibration is an open item):")
co, mo, so = table("   -- OPEN (cc<0.10):", inc[open_m], d[open_m])
cf, mf, sf = table("   -- FOREST (cc>0.50):", inc[for_m], d[for_m])

# ---------------------------------------------------------------- 4. CHECK scatter vs incidence spread
per_cell = df.groupby("cell").agg(inc_std=("incidence", "std"), std_d=("d_mm", "std")).dropna()
print("4. CHECK  per-cell scatter vs within-cell incidence spread "
      "(prediction: scatter grows with spread if incidence drives offset):")
sp_edges = [0, 0.5, 1, 2, 4, 8, 20]
print(f"  {'inc_std(deg)':>14s} {'median cell_std_d(mm)':>22s} {'n cells':>10s}")
xs = per_cell.inc_std.to_numpy(); ys = per_cell.std_d.to_numpy()
for i in range(len(sp_edges) - 1):
    b = (xs >= sp_edges[i]) & (xs < sp_edges[i + 1])
    if b.sum() < 50: continue
    print(f"  {sp_edges[i]:5.1f}-{sp_edges[i+1]:<5.1f}    {np.median(ys[b]):18.1f} {b.sum():12,d}")
print(f"  corr(cell scatter, incidence spread): Spearman {spearmanr(xs, ys).statistic:+.3f}\n")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
view = np.isfinite(d) & (np.abs(d) <= 300)
hb = ax[0].hexbin(inc[view], d[view], gridsize=60, bins="log", cmap="viridis", mincnt=1)
ax[0].plot(c, m, "w-", lw=2.5); ax[0].plot(c, m, "C3o-", lw=1.5, label="binned median")
ax[0].axhline(0, color="k", lw=.6); ax[0].set_xlim(0, 45); ax[0].set_ylim(-300, 300)
ax[0].set_xlabel("incidence angle to surface (deg)"); ax[0].set_ylabel("offset d (mm), gen1 vs gen2")
ax[0].set_title("PRIMARY: per-beam offset vs incidence (all returns)")
ax[0].legend(loc="upper right"); fig.colorbar(hb, ax=ax[0], label="log10 count")
ax[1].plot(c, m, "C0o-", label="all returns");
if len(co): ax[1].plot(co, mo, "C2s-", label="open (cc<0.10)")
if len(cf): ax[1].plot(cf, mf, "C1^-", label="forest (cc>0.50)")
ax[1].axhline(0, color="k", lw=.6); ax[1].set_xlim(0, 45)
ax[1].set_xlabel("incidence angle to surface (deg)"); ax[1].set_ylabel("median offset d (mm)")
ax[1].set_title("median offset vs incidence, by canopy cover"); ax[1].legend(); ax[1].grid(alpha=.3)
fig.suptitle(f"gen1 per-beam offset vs beam incidence angle (elba) — {lab}", y=1.0)
fig.savefig(f"figures/refdatum/offset_vs_incidence{suffix}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"wrote figures/refdatum/offset_vs_incidence{suffix}.png")
