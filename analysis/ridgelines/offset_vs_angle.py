#!/usr/bin/env python3
"""gen1 per-beam elevation offset vs gen2 (d_mm) as a function of a beam/terrain ANGLE axis:
  --x incidence  (default) : beam incidence to the local surface -- unifies slope+aspect+scan
  --x scan_angle           : |scan angle| (UNSIGNED; direction is already in the incidence recon)
  --x slope                : local surface slope

Reads the canonical per-beam table (beam_offset_table.parquet). Reports, in this order:
  1. PRIMARY  -- d_mm vs x, all in-grid returns: robust binned median + NMAD + n, plus a hexbin.
  2. CLEAN    -- within-cell residual (d_resid_mm) vs x: per-cell effects (real change, datum,
                 canopy height) differenced out; leverage lives in flight-line-overlap cells,
                 so we also report the overlap-weighted within-cell slope.
  3. STRATIFY -- the primary split by canopy cover, to see the leaf-state offset atop geometry.
  4. CHECK    -- per-cell scatter vs within-cell x-spread: a falsifiable prediction (scatter
                 grows with spread if x drives the offset).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/offset_vs_angle.py [--tile DIR] [--x incidence|scan_angle|slope] [--curv-max T]

Optional --curv-max restricts to near-planar cells (|curv_laplacian| <= T) to suppress
hillslope-diffusion / convex-concave real change; the figure form is identical, only the
point subset changes (output filename and title carry the threshold). No fitting is imposed
beyond robust bin statistics and simple correlations.
"""
import argparse, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr

XCFG = {  # values: how to derive the x axis from the table; hdr: table header; fname: file/label token
    "incidence":  dict(values=lambda df: df.incidence,          edges=np.arange(0, 48, 3),
                       label="incidence + δ (deg; δ = unknown nadir offset)", short="incidence", hdr="incid(deg)", xlim=(0, 45)),
    "scan_angle": dict(values=lambda df: df.scan_angle.abs(),   edges=np.arange(0, 20, 2),
                       label="|scan angle| (deg)",               short="|scan angle|", hdr="|scan|(deg)", xlim=(0, 18)),
    "scan_angle_signed": dict(values=lambda df: df.scan_angle,  edges=np.arange(-17, 18, 2),
                       label="scan angle (deg, signed)",         short="scan angle", hdr="scan(deg)", xlim=(-17, 17)),
    "slope":      dict(values=lambda df: df.slope,              edges=np.arange(0, 48, 3),
                       label="surface slope (deg)",              short="slope", hdr="slope(deg)", xlim=(0, 45)),
}
MIN_N = 300                              # drop bins below this (matches prior slope analysis)

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--x", default="incidence", choices=list(XCFG))
ap.add_argument("--curv-max", type=float, default=None)
ap.add_argument("--bin", type=float, default=None, help="override bin width (deg)")
ap.add_argument("--ridge", action="store_true", help="restrict to ridge_mask cells (divides: no overland flow)")
ap.add_argument("--cover-bands", action="store_true",
                help="stratify the right panel into canopy-cover PERCENT bands (one line per "
                     "band) instead of the two open/forest classes, so the progression with "
                     "forest fraction is visible rather than just its end members")
ap.add_argument("--offset", default="raw", choices=("raw", "corr"),
                help="raw = d_mm as measured (pre-registration); corr = d_mm_corr, with the "
                     "geoid, lateral tie, per-swath alignment and along-track drift applied")
ap.add_argument("--boresight", type=float, default=None,
                help="subtract this roll (mm/deg) * scan_angle from d_mm before analysis, and "
                     "recompute the within-cell residual -- to show the boresight asymmetry removed")
A = ap.parse_args()
cfg = XCFG[A.x]; XNAME = cfg["short"]; XLIM = cfg["xlim"]
EDGES = np.arange(cfg["edges"][0], cfg["edges"][-1] + A.bin, A.bin) if A.bin else cfg["edges"]

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet")
df = df[df.in_grid.values].copy()
lab, suffix = "all curvatures", ""
if A.curv_max is not None:
    keep = (df.curv_laplacian.abs() <= A.curv_max).to_numpy()
    frac = 100 * keep.mean(); df = df[keep].copy()
    lab = f"|Laplacian| <= {A.curv_max:g}  ({len(df):,} returns, {frac:.0f}% kept)"
    suffix = f"_curv{A.curv_max:g}"
if A.offset == "corr":                            # registration-corrected offset, then re-derive residual
    # d_resid_mm and the per-cell scatter in the table were built from the RAW d_mm, so they
    # must be re-derived here or the "within-cell" panels would mix two different offsets.
    df["d_mm"] = df.d_mm_corr
    df["d_resid_mm"] = df.d_mm - df.groupby("cell")["d_mm"].transform("mean")
    lab += "; REGISTRATION-CORRECTED"; suffix += "_reg"
if A.boresight is not None:                       # remove boresight roll, then re-derive residual
    df["d_mm"] = df.d_mm - A.boresight * df.scan_angle
    df["d_resid_mm"] = df.d_mm - df.groupby("cell")["d_mm"].transform("mean")
    lab += f"; boresight-corrected {A.boresight:g} mm/deg"; suffix += f"_bore{A.boresight:g}"
if A.ridge:                                       # divides only: zero contributing area -> no overland flow
    rm = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
    df = df[rm[df.cell.to_numpy()]].copy()
    lab += "; RIDGELINES"; suffix += "_ridge"
df["_x"] = cfg["values"](df).to_numpy(float)
xv = df["_x"].to_numpy(float); d = df.d_mm.to_numpy(float)
print(f"{len(df):,} returns  [x = {A.x}; {lab}]\n")

def nmad(x): return 1.4826 * np.median(np.abs(x - np.median(x)))

def binned(xa, ya, edges=EDGES, min_n=MIN_N):
    """robust median + NMAD + n per bin; returns (centers, med, nmad, n)."""
    c, m, s, k = [], [], [], []
    for i in range(len(edges) - 1):
        b = (xa >= edges[i]) & (xa < edges[i + 1]) & np.isfinite(ya)
        if b.sum() < min_n: continue
        c.append((edges[i] + edges[i + 1]) / 2); m.append(np.median(ya[b]))
        s.append(nmad(ya[b])); k.append(int(b.sum()))
    return np.array(c), np.array(m), np.array(s), k

def table(title, xa, ya):
    c, m, s, k = binned(xa, ya)
    print(title)
    print(f"  {cfg['hdr']:>10s} {'median':>9s} {'NMAD':>8s} {'n':>12s}")
    for a, mm, ss, nn in zip(c, m, s, k):
        print(f"  {a:10.0f} {mm:+9.1f} {ss:8.1f} {nn:12,d}")
    fin = np.isfinite(ya)
    pr = np.corrcoef(xa[fin], ya[fin])[0, 1]; sp = spearmanr(xa[fin], ya[fin]).statistic
    print(f"  corr(offset, {XNAME}): Pearson {pr:+.3f}  Spearman {sp:+.3f}\n")
    return c, m, s

# ---------------------------------------------------------------- 1. PRIMARY
c, m, s = table(f"1. PRIMARY  d_mm vs {XNAME} (all in-grid returns):", xv, d)

# ---------------------------------------------------------------- 2. CLEAN (within-cell)
dr = df.d_resid_mm.to_numpy(float)
table(f"2. CLEAN  within-cell residual d_resid_mm vs {XNAME} (per-cell effects removed):", xv, dr)
# overlap-powered slope: regress residual on within-cell x-deviation, weighting by how much
# x actually varies in the cell (within-cell x spread). Cells with ~0 spread carry ~0 info.
cxm = df.groupby("cell")["_x"].transform("mean").to_numpy(float)
x_dev = xv - cxm
w = df.groupby("cell")["_x"].transform("std").to_numpy(float)
ok = np.isfinite(w) & np.isfinite(x_dev) & np.isfinite(dr) & (w > 0)
W = w[ok]
slope_wc = np.sum(W * x_dev[ok] * dr[ok]) / np.sum(W * x_dev[ok] ** 2)
print(f"   within-cell slope d(offset)/d({XNAME}) = {slope_wc:+.2f} mm/deg "
      f"(overlap-weighted, {ok.sum():,} returns; leverage from high within-cell {XNAME} spread)\n")

# ---------------------------------------------------------------- 3. STRATIFY by canopy cover
cc = df.canopy_cover.to_numpy(float)
open_m = np.isfinite(cc) & (cc < 0.10)   # provisional (forest/open threshold not yet calibrated)
for_m  = np.isfinite(cc) & (cc > 0.50)
COVER_EDGES = np.array([0, .05, .10, .20, .35, .50, 1.01])   # same bands as the slope x cover model
bands = []
if A.cover_bands:
    print(f"3. STRATIFY by canopy-cover band (percent forest fraction):")
    for lo, hi in zip(COVER_EDGES[:-1], COVER_EDGES[1:]):
        bm = np.isfinite(cc) & (cc >= lo) & (cc < hi)
        if bm.sum() < 500:
            print(f"   -- cover {100*lo:.0f}-{100*hi:.0f}%: n={bm.sum():,} (sparse, skipped)"); continue
        cb, mb, sb = table(f"   -- cover {100*lo:.0f}-{100*min(hi,1.0):.0f}% (n={bm.sum():,}):",
                           xv[bm], d[bm])
        bands.append((lo, min(hi, 1.0), cb, mb, int(bm.sum())))
else:
    print(f"3. STRATIFY by canopy cover (PROVISIONAL: open cc<0.10 n={open_m.sum():,}, "
          f"forest cc>0.50 n={for_m.sum():,}; threshold calibration is an open item):")
    co, mo, so = table("   -- OPEN (cc<0.10):", xv[open_m], d[open_m])
    cf, mf, sf = table("   -- FOREST (cc>0.50):", xv[for_m], d[for_m])

# ---------------------------------------------------------------- 4. CHECK scatter vs x-spread
per_cell = df.groupby("cell").agg(x_std=("_x", "std"), std_d=("d_mm", "std")).dropna()
print(f"4. CHECK  per-cell scatter vs within-cell {XNAME} spread "
      f"(prediction: scatter grows with spread if {XNAME} drives offset):")
sp_edges = [0, 0.5, 1, 2, 4, 8, 20]
print(f"  {XNAME+' std(deg)':>16s} {'median cell_std_d(mm)':>22s} {'n cells':>10s}")
xs = per_cell.x_std.to_numpy(); ys = per_cell.std_d.to_numpy()
for i in range(len(sp_edges) - 1):
    b = (xs >= sp_edges[i]) & (xs < sp_edges[i + 1])
    if b.sum() < 50: continue
    print(f"  {sp_edges[i]:6.1f}-{sp_edges[i+1]:<6.1f}  {np.median(ys[b]):18.1f} {b.sum():12,d}")
print(f"  corr(cell scatter, {XNAME} spread): Spearman {spearmanr(xs, ys).statistic:+.3f}\n")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
view = np.isfinite(d) & (np.abs(d) <= 300)
hb = ax[0].hexbin(xv[view], d[view], gridsize=60, bins="log", cmap="viridis", mincnt=1)
ax[0].plot(c, m, "w-", lw=2.5); ax[0].plot(c, m, "C3o-", lw=1.5, label="binned median")
ax[0].axhline(0, color="k", lw=.6); ax[0].set_xlim(*XLIM); ax[0].set_ylim(-300, 300)
ax[0].set_xlabel(cfg["label"]); ax[0].set_ylabel("offset d (mm) = gen1 − gen2   (+ = ground lower in 2021)")
ax[0].set_title(f"PRIMARY: per-beam offset vs {XNAME} (all returns)")
ax[0].legend(loc="upper right"); fig.colorbar(hb, ax=ax[0], label="log10 count")
if A.cover_bands:
    ax[1].plot(c, m, "k-", lw=2.2, alpha=.55, label="all returns", zorder=1)
    cmap = plt.get_cmap("viridis")
    for i, (lo, hi, cb, mb, nb) in enumerate(bands):
        if not len(cb): continue
        col = cmap(0.08 + 0.84 * i / max(len(bands) - 1, 1))
        ax[1].plot(cb, mb, "o-", ms=4, lw=1.6, color=col,
                   label=f"{100*lo:.0f}\u2013{100*hi:.0f}% cover (n={nb:,})")
    ttl = f"median offset vs {XNAME}, by forest fraction"
else:
    ax[1].plot(c, m, "C0o-", label="all returns")
    if len(co): ax[1].plot(co, mo, "C2s-", label="open (cc<0.10)")
    if len(cf): ax[1].plot(cf, mf, "C1^-", label="forest (cc>0.50)")
    ttl = f"median offset vs {XNAME}, by canopy cover"
ax[1].axhline(0, color="k", lw=.6); ax[1].set_xlim(*XLIM)
ax[1].set_xlabel(cfg["label"]); ax[1].set_ylabel("median offset d (mm)   [gen1 − gen2; + = lower in 2021]")
ax[1].set_title(ttl); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
import os as _os
_tile = _os.path.basename(A.tile.rstrip("/"))                       # tag figures by tile
_tt = "" if _tile == "elba_fulldensity" else f"_{_tile}"
if A.bin is not None: suffix += f"_bin{A.bin:g}"   # bin width changes the curve: keep it in the name
if A.cover_bands: suffix += "_cbands"
fig.suptitle(f"gen1 per-beam offset vs {XNAME} ({_tile}) — {lab}", y=1.0)
fig.savefig(f"figures/refdatum/offset_vs_{A.x}{suffix}{_tt}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"wrote figures/refdatum/offset_vs_{A.x}{suffix}{_tt}.png")
