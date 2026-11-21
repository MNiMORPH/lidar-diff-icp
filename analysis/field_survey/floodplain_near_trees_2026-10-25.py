#!/usr/bin/env python
"""The stacked 2021 return record near and around the trees (Andy, 2026-11-21: the zone with more nettle and goldenrod).
Cells: floodplain-mask cells of the field extent that are not tree (intensity_2021_1m_tree.tif), binned by distance to
the nearest tree cell: crown edge 0-3 m, near 3-10 m, open > 10 m (bands MINE); and the same split north/south of
N 4885480. Returns of both 2021 clouds in those cells, height = z minus the official USGS 1 m DEM of the cell.
Outputs  data/field/2026-10-25/intensity/near_trees_log.txt; figures/rtk_2026-10-25_floodplain_near_trees.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/near_trees_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = 380, 480; NSPLIT = 4885480.0
with rasterio.open(f"{OUT}/intensity_2021_1m_tree.tif") as r: TREE = r.read(1) == 1
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_field.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan
z10 = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); xs, ys, fpm = z10["xs"], z10["ys"], z10["fpm"].astype(bool)
cy = N1 - (np.arange(ny) + 0.5) * RES; cx = E0 + (np.arange(nx) + 0.5) * RES
FP = fpm[np.clip(np.round((ys[0] - cy) / 10).astype(int), 0, len(ys) - 1)[:, None], np.clip(np.round((cx - xs[0]) / 10).astype(int), 0, len(xs) - 1)[None, :]]
dist = distance_transform_edt(~TREE) * RES; BAND = np.full((ny, nx), -1); BAND[(dist > 0) & (dist <= 3)] = 0; BAND[(dist > 3) & (dist <= 10)] = 1; BAND[dist > 10] = 2; BAND[~FP | TREE] = -1
NORTH = (cy >= NSPLIT)[:, None].repeat(nx, 1)
names = ["crown edge, 0-3 m from a tree cell", "near trees, 3-10 m", "open, > 10 m"]
say("cells per band (floodplain, non-tree): " + ", ".join(f"{names[k]}: {int((BAND == k).sum())} (north {int(((BAND == k) & NORTH).sum())}, south {int(((BAND == k) & ~NORTH).sum())})" for k in range(3)))
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); b = BAND[i, j]; k = b >= 0
            parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][k] - DEM[i[k], j[k]], inten=np.asarray(ch.intensity)[m][k].astype(float), band=b[k], north=NORTH[i[k], j[k]], nr=np.asarray(ch.number_of_returns)[m][k])))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
def stats(q):
    return dict(n=len(q), med=q.dz.median(), p05=q.dz.quantile(.05), p95=q.dz.quantile(.95), shelf=100 * (q.dz < -0.08).mean(), herb=100 * ((q.dz > 0.2) & (q.dz <= 1.0)).mean(), canopy=100 * (q.dz > 1.0).mean(), single=100 * (q.nr == 1).mean(),
                int_peak=q.inten[(q.dz > -0.03) & (q.dz < 0.03)].median(), int_shelf=q.inten[(q.dz < -0.08) & (q.dz > -0.25)].median(), int_herb=q.inten[(q.dz > 0.2) & (q.dz <= 1.0)].median())
say(f"{'set':42s} {'n':>8s} {'median':>7s} {'p05':>7s} {'p95':>7s} {'shelf<-0.08':>11s} {'0.2-1 m':>8s} {'>1 m':>6s} {'single':>7s} {'int peak':>8s} {'int shelf':>9s} {'int herb':>8s}")
for k in range(3):
    for lab, q in ((names[k], R[R.band == k]), (names[k] + " | north", R[(R.band == k) & R.north]), (names[k] + " | south", R[(R.band == k) & ~R.north])):
        s = stats(q); say(f"{lab:42s} {s['n']:8d} {s['med']:+7.3f} {s['p05']:+7.3f} {s['p95']:+7.3f} {s['shelf']:10.2f}% {s['herb']:7.2f}% {s['canopy']:5.2f}% {s['single']:6.1f}% {s['int_peak']:8.0f} {s['int_shelf']:9.0f} {s['int_herb']:8.0f}")
edges = np.arange(-0.50, 1.00 + 0.005, 0.01); mids = 0.5 * (edges[1:] + edges[:-1]); norm = Normalize(vmin=2000, vmax=5000); cmap = plt.get_cmap("viridis")
fig, axs = plt.subplots(2, 2, figsize=(16, 9.5))
for k in range(2):
    q = R[R.band == k]; h = np.histogram(q.dz, bins=edges)[0] / (len(q) * 0.01); b = np.digitize(q.dz, edges) - 1; med = np.array([q.inten[b == m].median() if (b == m).sum() >= 20 else np.nan for m in range(len(mids))])
    ax = axs[0, k]; ax.bar(mids, np.where(h > 0, h, np.nan), width=0.01, color=[cmap(norm(v)) if np.isfinite(v) else "0.85" for v in med], edgecolor="none"); ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.5, 1.0)
    ax.set_title(f"{names[k]} (n {len(q):,})", fontsize=10); ax.set_xlabel("2021 return minus official USGS 1 m DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.spines[["top", "right"]].set_visible(False)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs[0, :], shrink=0.8, label="median intensity in the bar")
COL = ["#d55e00", "#e69f00", "#0072b2"]
for col, (lab, sel) in zip(axs[1, :], (("all three bands", np.ones(len(R), bool)), ("south of N 4885480 only", ~R.north.to_numpy()))):
    for k in range(3):
        q = R[(R.band == k) & sel]; col.hist(q.dz, bins=edges, density=True, histtype="step", color=COL[k], lw=1.5, label=f"{names[k]} (n {len(q):,})")
    col.set_yscale("log"); col.axvline(0, color="0.4", lw=0.8, ls=":"); col.set_xlim(-0.5, 1.0); col.set_xlabel("2021 return minus official USGS 1 m DEM (m)"); col.set_ylabel("density (1/m), log"); col.set_title(lab, fontsize=10); col.legend(frameon=False, fontsize=8); col.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_near_trees.png", dpi=int(os.environ.get("FIGDPI", 160))); say("wrote figures/rtk_2026-10-25_floodplain_near_trees.png")
LOG.close()
