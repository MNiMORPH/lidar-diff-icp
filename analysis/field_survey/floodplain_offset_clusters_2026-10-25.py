#!/usr/bin/env python
"""Does the lidar-minus-RTK offset cluster by location on the floodplain? (Andy, 2026-11-20: "see if there is
location-based clustering on offset to help define the binning", for densities re-binned by place, not by rover.)

Per floodplain FIX point (all three rovers, 147 points), from the per-return tables (floodplain_returns_2026-10-25.py,
benchmark-tied frame) and the imagery covariates:
  med21_sub_m    median height above the RTK ground of the 2021 returns within 1 m that lie below +0.50 m (the sward and
                 ground population; the 0.50 m cut, MINE, removes the tree canopy at the lot edge and the CORS line)
  ground_frac    share of the 2021 returns within 1 m below +0.05 m (MINE: the ground shelf seen in the pooled densities)
  canopy_frac    share of the 2021 returns within 1 m above +0.50 m
  p10_21_m       10th percentile of all 2021 returns within 1 m
  med08_3m_m     median height of the 2008 returns within 3 m (star frame; 1 m holds ~2 returns)
  product offsets gen2_minus_rtk_m, gen1_minus_rtk_m; NDVI metrics from imagery/floodplain_points_covariates.csv
Then (1) the empirical semivariogram of med21_sub_m and ground_frac over point pairs (lag bins 0-10, 10-20, ... m) to
see whether near points share their offset; (2) Ward agglomerative clustering with a spatial connectivity graph
(each point linked to its 8 nearest neighbors, MINE) on standardized (med21_sub_m, ground_frac), for k = 2..7, with
the fraction of variance explained per k; (3) maps and the per-cluster table. The cluster count to adopt is a
PROPOSAL; the script prints the table for every k so that the choice can be read, and writes labels for k = 2..7.
Outputs  data/field/2026-10-25/returns/floodplain_offset_points.csv      the per-point table with cluster_k2..k7
         data/field/2026-10-25/returns/floodplain_offset_clusters_log.txt
         figures/rtk_2026-10-25_floodplain_offset_map.png                 points colored by med21_sub_m and ground_frac, on NAIP 2021
         figures/rtk_2026-10-25_floodplain_offset_variogram.png
         figures/rtk_2026-10-25_floodplain_offset_clusters_k<k>.png       cluster map for the proposed k (and K env var)
"""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_offset_clusters_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
DPI = int(os.environ.get("FIGDPI", 170)); LOT = (579979.0, 4885340.0); SUB = 0.50; GROUND = 0.05; KNN = 8
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); C = pd.read_csv("data/field/2026-10-25/imagery/floodplain_points_covariates.csv")
rows = []
for tag, group in (("rsplus_b", "floodplain_RSplus_B"), ("rsplus_a", "floodplain_RSplus_A"), ("cors", "floodplain_CORS_Wickert")):
    t21 = pd.read_csv(f"{D}/{tag}_returns_2021.csv.gz"); t08 = pd.read_csv(f"{D}/{tag}_returns_2008.csv")
    pts = P[(P.group == group) & (P.status == "FIX")]
    for _, r in pts.iterrows():
        a = t21[(t21.point == r.point) & t21.radius_1m].dz_m; s = a[a < SUB]; b = t08[(t08.point == r.point) & t08.radius_3m].dz_star_m
        rows.append(dict(group=group, rover=group.replace("floodplain_", "").replace("_Wickert", ""), point=int(r.point), E=r.easting_tied, N=r.northing_tied,
                         d_lot=float(np.hypot(r.easting_tied - LOT[0], r.northing_tied - LOT[1])), n21=len(a), med21_sub_m=s.median(), p10_21_m=a.quantile(.1),
                         ground_frac=(a < GROUND).mean(), canopy_frac=(a > SUB).mean(), n08_3m=len(b), med08_3m_m=b.median() if len(b) else np.nan,
                         gen2_minus_rtk_m=r.gen2_minus_rtk_m, gen1_minus_rtk_m=r.gen1_minus_rtk_m, lowveg_gen2=r.lowveg_gen2))
X = pd.DataFrame(rows).merge(C[["group", "point", "metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "s2_2021-05-01"]], on=["group", "point"], how="left")
say(f"points: {len(X)} floodplain FIX points; sub-canopy cut {SUB} m, ground cut {GROUND} m, kNN {KNN} (all MINE)")
say("per-rover medians:", X.groupby("rover")[["med21_sub_m", "ground_frac", "canopy_frac", "med08_3m_m", "metric_ndvi_apr", "metric_ndvi_oct"]].median().round(3).to_string())

# ---- semivariogram
def variogram(v, lags=np.arange(0, 220, 10)):
    E, N = X.E.to_numpy(), X.N.to_numpy(); out = []
    i, j = np.triu_indices(len(X), 1); d = np.hypot(E[i] - E[j], N[i] - N[j]); g = 0.5 * (v[i] - v[j]) ** 2; ok = np.isfinite(g)
    for a, b in zip(lags[:-1], lags[1:]):
        m = ok & (d >= a) & (d < b); out.append((0.5 * (a + b), g[m].mean() if m.sum() > 5 else np.nan, int(m.sum())))
    return pd.DataFrame(out, columns=["lag_m", "semivariance", "pairs"])
V1 = variogram(X.med21_sub_m.to_numpy()); V2 = variogram(X.ground_frac.to_numpy()); V3 = variogram(X.med08_3m_m.to_numpy())
say("semivariogram of med21_sub_m (m^2) and ground_frac, by lag (sill = half the variance over all pairs: "
    f"{0.5 * np.nanvar(X.med21_sub_m) * 2:.5f} m^2 i.e. variance {np.nanvar(X.med21_sub_m):.5f}, ground_frac variance {np.nanvar(X.ground_frac):.5f}):")
say(pd.concat([V1, V2.semivariance.rename("semivar_ground_frac"), V3.semivariance.rename("semivar_med08_3m")], axis=1).round(5).to_string())
say(f"2008 median within 3 m: variance over all points {np.nanvar(X.med08_3m_m):.5f} m^2")
fig, axs = plt.subplots(1, 3, figsize=(14, 3.8))
axs[2].plot(V3.lag_m, V3.semivariance, "o-", color="#009e73"); axs[2].axhline(np.nanvar(X.med08_3m_m), color="0.5", lw=0.8, ls="--"); axs[2].set_xlabel("separation (m)"); axs[2].set_ylabel("semivariance of 2008 median within 3 m (m$^2$)")
axs[0].plot(V1.lag_m, V1.semivariance, "o-", color="#0072b2"); axs[0].axhline(np.nanvar(X.med21_sub_m), color="0.5", lw=0.8, ls="--"); axs[0].set_xlabel("separation (m)"); axs[0].set_ylabel("semivariance of sub-canopy median (m$^2$)")
axs[1].plot(V2.lag_m, V2.semivariance, "o-", color="#d55e00"); axs[1].axhline(np.nanvar(X.ground_frac), color="0.5", lw=0.8, ls="--"); axs[1].set_xlabel("separation (m)"); axs[1].set_ylabel("semivariance of ground fraction")
for ax in axs: ax.spines[["top", "right"]].set_visible(False); ax.set_ylim(bottom=0)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_offset_variogram.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_offset_variogram.png (dashed = variance over all points)")

# ---- spatially connected Ward clustering
F = X[["med21_sub_m", "ground_frac"]].to_numpy(); F = (F - F.mean(0)) / F.std(0)
conn = kneighbors_graph(X[["E", "N"]].to_numpy(), KNN, include_self=False)
tot = ((F - F.mean(0)) ** 2).sum()
say(f"Ward clustering on standardized (med21_sub_m, ground_frac) with a {KNN}-nearest-neighbor spatial connectivity graph; variance explained by k:")
for k in range(2, 8):
    lab = AgglomerativeClustering(n_clusters=k, connectivity=conn, linkage="ward").fit_predict(F); X[f"cluster_k{k}"] = lab
    within = sum(((F[lab == c] - F[lab == c].mean(0)) ** 2).sum() for c in range(k)); say(f"   k {k}: explained {100 * (1 - within / tot):5.1f} %  sizes {np.bincount(lab).tolist()}")
for name, cols in (("med08", ["med08_3m_m"]), ("all3", ["med21_sub_m", "ground_frac", "med08_3m_m"])):
    G = X[cols].to_numpy(float); ok = np.isfinite(G).all(1); G = (G[ok] - G[ok].mean(0)) / G[ok].std(0); conn2 = kneighbors_graph(X.loc[ok, ["E", "N"]].to_numpy(), KNN, include_self=False); tot2 = ((G - G.mean(0)) ** 2).sum()
    say(f"Ward clustering on {cols} ({int(ok.sum())} points with values), same spatial connectivity; variance explained by k:")
    for k in range(2, 7):
        lab = AgglomerativeClustering(n_clusters=k, connectivity=conn2, linkage="ward").fit_predict(G); X.loc[ok, f"cluster_{name}_k{k}"] = lab
        within = sum(((G[lab == c] - G[lab == c].mean(0)) ** 2).sum() for c in range(k)); say(f"   k {k}: explained {100 * (1 - within / tot2):5.1f} %  sizes {np.bincount(lab).tolist()}")
X.to_csv(f"{D}/floodplain_offset_points.csv", index=False, float_format="%.4f")
def table(k, col=None):
    g = X.groupby(col or f"cluster_k{k}")
    T = g.agg(n=("point", "size"), rovers=("rover", lambda v: dict(v.value_counts())), N_range=("N", lambda v: f"{v.min():.0f}-{v.max():.0f}"), E_range=("E", lambda v: f"{v.min():.0f}-{v.max():.0f}"),
              d_lot=("d_lot", "median"), med21_sub=("med21_sub_m", "median"), p10_21=("p10_21_m", "median"), ground_frac=("ground_frac", "median"), canopy_frac=("canopy_frac", "median"),
              med08_3m=("med08_3m_m", "median"), ndvi_apr=("metric_ndvi_apr", "median"), ndvi_sum=("metric_ndvi_sum", "median"), ndvi_oct=("metric_ndvi_oct", "median"), naip21=("naip_2021", "median"))
    return T.round(3)
for k in (2, 3, 4, 5):
    say(f"\ncluster table (2021 features), k = {k}:"); say(table(k).to_string())
for k in (2, 3, 4):
    say(f"\ncluster table (2008 median within 3 m), k = {k}:"); say(table(k, f"cluster_med08_k{k}").to_string())
for k in (3, 4):
    say(f"\ncluster table (all three features), k = {k}:"); say(table(k, f"cluster_all3_k{k}").to_string())

# ---- maps
naip = None
try:
    import planetary_computer as pc, rasterio
    from pystac_client import Client
    from pyproj import Transformer
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds
    E0, E1, N0, N1 = 579860, 580120, 4885320, 4885710
    tr = Transformer.from_crs("EPSG:6344", "EPSG:4326", always_xy=True); lon0, lat0 = tr.transform(E0, N0); lon1, lat1 = tr.transform(E1, N1)
    cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
    it = [i for i in cat.search(collections=["naip"], bbox=[lon0, lat0, lon1, lat1]).items() if i.datetime.year == 2021][0]
    with rasterio.open(it.assets["image"].href) as rr:
        w = from_bounds(*transform_bounds("EPSG:6344", rr.crs, E0, N0, E1, N1), transform=rr.transform); naip = np.moveaxis(rr.read([1, 2, 3], window=w), 0, -1)
    say(f"NAIP {it.id} fetched for the maps")
except Exception as e:
    E0, E1, N0, N1 = 579860, 580120, 4885320, 4885710; say("NAIP not fetched:", repr(e)[:100])
def base(ax):
    if naip is not None: ax.imshow(naip, extent=(E0, E1, N0, N1), alpha=0.85)
    ax.set_xlim(E0, E1); ax.set_ylim(N0, N1); ax.set_aspect("equal"); ax.plot([LOT[0]], [LOT[1]], "w+", ms=12, mew=2)
    ax.tick_params(labelsize=7)
fig, axs = plt.subplots(1, 3, figsize=(16, 7.5))
for ax, col, cm, lab, vr in ((axs[0], "med21_sub_m", "viridis", "2021 sub-canopy median above RTK ground (m)", (0.08, 0.22)), (axs[1], "ground_frac", "magma", "share of 2021 returns within 5 cm of the ground", (0, 0.15)), (axs[2], "med08_3m_m", "viridis", "2008 median within 3 m above RTK ground (m)", (0.1, 0.5))):
    base(ax); sc = ax.scatter(X.E, X.N, c=X[col], cmap=cm, vmin=vr[0], vmax=vr[1], s=22, edgecolors="k", linewidths=0.3); plt.colorbar(sc, ax=ax, shrink=0.6, label=lab)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_offset_map.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_offset_map.png")
for col, k, png in (("cluster_k3", 3, "figures/rtk_2026-10-25_floodplain_offset_clusters_k3.png"), ("cluster_k4", 4, "figures/rtk_2026-10-25_floodplain_offset_clusters_k4.png"), ("cluster_med08_k3", 3, "figures/rtk_2026-10-25_floodplain_offset_clusters_med08_k3.png"), ("cluster_med08_k4", 4, "figures/rtk_2026-10-25_floodplain_offset_clusters_med08_k4.png"), ("cluster_all3_k4", 4, "figures/rtk_2026-10-25_floodplain_offset_clusters_all3_k4.png")):
    fig, ax = plt.subplots(figsize=(6.5, 8)); base(ax); cols = plt.get_cmap("tab10")
    for c in range(k):
        q = X[X[col] == c]; ax.scatter(q.E, q.N, color=cols(c), s=26, edgecolors="k", linewidths=0.3, label=f"cluster {c} (n {len(q)}, 2021 {q.med21_sub_m.median():.3f} m, ground {100 * q.ground_frac.median():.0f} %, 2008 {q.med08_3m_m.median():.3f} m)")
    ax.legend(fontsize=7, loc="lower left"); fig.tight_layout(); fig.savefig(png, dpi=DPI); plt.close(fig); say("wrote", png)
LOG.close()
