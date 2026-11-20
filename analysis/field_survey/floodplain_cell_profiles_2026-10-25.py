#!/usr/bin/env python
"""Does the vertical distribution of 2021 returns in each 1 m cell define vegetation cover? (Andy, 2026-11-20.)
Field-site extent as the intensity grid (E 579780-580160, N 4885280-4885760). All 2021 returns of both tiles (overlap
from the first). Ground proxy per 1 m cell = the lowest return within a 5 x 5 m window around it (MINE; the RTK points
say this proxy sits ~5 cm above the soil in dense sward and on the soil where the laser reached it). Per cell with at
least 8 returns (MINE): fractions of returns by height above the proxy in the bands < 0.05 (ground), 0.05-0.30 (low
sward), 0.30-1.0 (tall herb), 1-3 (shrub), > 3 m (tree); highest return; mean intensity and its CV. k-means (k = 6,
MINE; k = 4-8 inertia printed) on the standardized features; cluster table; cluster map; cross-tabulation against
the NAIP-2021 imagery types at 10 m; and at the 147 RTK points the proxy minus the surveyed ground.
Outputs  data/field/2026-10-25/intensity/cell_profiles_1m.csv (one row per cell), cell_profiles_cluster.tif, cell_profiles_log.txt
         figures/rtk_2026-10-25_floodplain_cell_profiles.png  cluster map + the mean band profile of each cluster
"""
import os, json
import numpy as np, pandas as pd, laspy, rasterio
from rasterio.transform import from_origin
from scipy.ndimage import minimum_filter
from sklearn.cluster import KMeans
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/cell_profiles_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = int((E1 - E0) / RES), int((N1 - N0) / RES); K = int(os.environ.get("K", 6)); NMIN = 8
BANDS = [(-np.inf, 0.05, "ground"), (0.05, 0.30, "low sward"), (0.30, 1.0, "tall herb"), (1.0, 3.0, "shrub"), (3.0, np.inf, "tree")]
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if m.any(): parts.append(np.column_stack([x[m], y[m], np.asarray(ch.z)[m], np.asarray(ch.intensity)[m].astype(float)]))
R = np.vstack(parts); x, y, z, inten = R.T; del parts
i = np.clip(((N1 - y) / RES).astype(int), 0, ny - 1); j = np.clip(((x - E0) / RES).astype(int), 0, nx - 1); cell = i * nx + j
zmin = np.full(ny * nx, np.inf); np.minimum.at(zmin, cell, z); zmin[~np.isfinite(zmin)] = np.nan
g = minimum_filter(np.where(np.isfinite(zmin), zmin, np.inf).reshape(ny, nx), size=5, mode="nearest").ravel(); g[~np.isfinite(g)] = np.nan
h = z - g[cell]; n = np.bincount(cell, minlength=ny * nx)
feat = {}
for lo, hi, name in BANDS: feat[name] = np.bincount(cell[(h >= lo) & (h < hi)], minlength=ny * nx) / np.maximum(n, 1)
hmax = np.full(ny * nx, -np.inf); np.maximum.at(hmax, cell, h); hmax[~np.isfinite(hmax)] = np.nan
with np.errstate(invalid="ignore", divide="ignore"):
    mi = np.bincount(cell, weights=inten, minlength=ny * nx) / n; cv = np.sqrt(np.maximum(np.bincount(cell, weights=inten ** 2, minlength=ny * nx) / n - mi ** 2, 0)) / mi
C = pd.DataFrame(dict(cell=np.arange(ny * nx), n=n, ground_proxy_m=g, hmax=hmax, intensity=mi, cv=cv, **feat)); C["E"] = E0 + (C.cell % nx + 0.5) * RES; C["N"] = N1 - (C.cell // nx + 0.5) * RES
ok = (C.n >= NMIN) & np.isfinite(C.intensity) & np.isfinite(C.cv); say(f"returns {len(R):,}; cells {ny * nx}, with >= {NMIN} returns {int(ok.sum())}")
# RTK check of the proxy
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")].copy()
fi = np.clip(((N1 - F.northing_tied.values) / RES).astype(int), 0, ny - 1); fj = np.clip(((F.easting_tied.values - E0) / RES).astype(int), 0, nx - 1); F["proxy_minus_rtk"] = g[fi * nx + fj] - F.navd88_tied_m.values
say(f"ground proxy (5 m lowest return) minus surveyed ground at the {len(F)} floodplain FIX points: median {F.proxy_minus_rtk.median():+.3f} m, p10 {F.proxy_minus_rtk.quantile(.1):+.3f}, p90 {F.proxy_minus_rtk.quantile(.9):+.3f}; by rover " + ", ".join(f"{k.replace('floodplain_', '')} {v:+.3f}" for k, v in F.groupby("group").proxy_minus_rtk.median().items()))
# clustering
X = C.loc[ok, ["ground", "low sward", "tall herb", "shrub", "tree", "intensity", "cv"]].copy(); X["log_hmax"] = np.log10(C.loc[ok, "hmax"].clip(lower=0.01) + 0.05)
Z = (X - X.mean()) / X.std()
say("k-means inertia by k (standardized features ground/low/tall/shrub/tree fractions, intensity, CV, log top height): " + ", ".join(f"k {k}: {KMeans(n_clusters=k, n_init=4, random_state=0).fit(Z).inertia_:.0f}" for k in range(3, 9)))
km = KMeans(n_clusters=K, n_init=10, random_state=0).fit(Z); C.loc[ok, "cluster"] = km.labels_
# order clusters by mean height band for readability
order = C[ok].groupby("cluster").apply(lambda q: (q["low sward"] * 0.15 + q["tall herb"] * 0.6 + q["shrub"] * 2 + q["tree"] * 8).mean()).sort_values().index.tolist(); remap = {c: k for k, c in enumerate(order)}; C.loc[ok, "cluster"] = C.loc[ok, "cluster"].map(remap)
T = C[ok].groupby("cluster").agg(cells=("cell", "size"), ground=("ground", "mean"), low_sward=("low sward", "mean"), tall_herb=("tall herb", "mean"), shrub=("shrub", "mean"), tree=("tree", "mean"), hmax_med=("hmax", "median"), intensity=("intensity", "median"), cv=("cv", "median"), n_med=("n", "median"))
pd.set_option("display.width", 250); say(f"\ncluster table (k = {K}; means of the band fractions, medians otherwise):"); say(T.round(3).to_string())
# imagery types at 10 m
z10 = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); names = list(z10["names"]); xs, ys, fpm = z10["xs"], z10["ys"], z10["fpm"].astype(bool); st = z10["stack"]
Ty = pd.read_csv("data/field/2026-10-25/imagery/floodplain_types_points.csv")
ri = np.clip(np.round((ys[0] - C.N.values) / 10).astype(int), 0, len(ys) - 1); rj = np.clip(np.round((C.E.values - xs[0]) / 10).astype(int), 0, len(xs) - 1)
C["floodplain"] = fpm[ri, rj]; C["naip21_ndvi"] = st[names.index("naip_2021-08_NDVI")][ri, rj]; C["ndvi_apr26"] = st[names.index("s2_2026-04-10_NDVI")][ri, rj]; C["ndvi_oct26"] = st[names.index("s2_2026-10-17_NDVI")][ri, rj]; C["nelike21"] = st[names.index("nelike_p_s2_20211018")][ri, rj]
say("\nper cluster, medians of the optical layers at 10 m (floodplain cells): "); say(C[ok & C.floodplain].groupby("cluster")[["naip21_ndvi", "ndvi_apr26", "ndvi_oct26", "nelike21"]].median().round(3).to_string())
say("\ncluster share north (N >= 4885480) vs south of the split, floodplain cells (%):"); say(pd.crosstab(C.loc[ok & C.floodplain, "cluster"], np.where(C.loc[ok & C.floodplain, "N"] >= 4885480, "north", "south"), normalize="columns").mul(100).round(1).to_string())
# clusters at the RTK points
F["cluster"] = C.cluster.values[fi * nx + fj]; say("\nclusters at the RTK points by rover:"); say(pd.crosstab(F.group.str.replace("floodplain_", ""), F.cluster).to_string())
C.to_csv(f"{OUT}/cell_profiles_1m.csv", index=False, float_format="%.4f")
tr = from_origin(E0, N1, RES, RES)
with rasterio.open(f"{OUT}/cell_profiles_cluster.tif", "w", driver="GTiff", height=ny, width=nx, count=1, dtype="uint8", nodata=255, crs="EPSG:6344", transform=tr, compress="deflate") as dst:
    dst.write(np.where(ok, C.cluster.fillna(255), 255).reshape(ny, nx).astype("uint8"), 1)
say(f"wrote {OUT}/cell_profiles_1m.csv and cell_profiles_cluster.tif")
cols = plt.get_cmap("tab10"); fig = plt.figure(figsize=(18, 8)); ax = fig.add_subplot(1, 2, 1)
img = np.where(ok, C.cluster, np.nan).reshape(ny, nx); ax.imshow(img, extent=(E0, E1, N0, N1), cmap=ListedColormap([cols(k) for k in range(K)]), vmin=-0.5, vmax=K - 0.5, interpolation="nearest"); ax.plot(F.easting_tied, F.northing_tied, ".", color="k", ms=2); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
from matplotlib.patches import Patch; ax.legend(handles=[Patch(color=cols(k), label=f"cluster {k}: top {T.hmax_med[k]:.2f} m, intensity {T.intensity[k]:.0f}, tree {100 * T.tree[k]:.0f} %") for k in range(K)], fontsize=7, loc="lower left")
ax2 = fig.add_subplot(1, 2, 2); w = 0.8 / K; labels = [b[2] for b in BANDS]
for k in range(K): ax2.bar(np.arange(len(BANDS)) + (k - K / 2 + 0.5) * w, T.loc[k, ["ground", "low_sward", "tall_herb", "shrub", "tree"]].values, w, color=cols(k), label=f"cluster {k} (n {int(T.cells[k])})")
ax2.set_xticks(range(len(BANDS))); ax2.set_xticklabels(labels); ax2.set_ylabel("mean fraction of the cell's returns"); ax2.legend(frameon=False, fontsize=8); ax2.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_cell_profiles.png", dpi=int(os.environ.get("FIGDPI", 170))); say("wrote figures/rtk_2026-10-25_floodplain_cell_profiles.png")
LOG.close()
