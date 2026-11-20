#!/usr/bin/env python
"""Is there a "full shape factor" of the lidar return-height distribution? A test on the 147 floodplain FIX points
(Andy, 2026-11-20). For each point the 2021 returns within 3 m (~490 per point; 1 m gives ~54, too few for shape)
below +0.50 m (the sward-and-ground population; the 0.50 m cut, 3 m radius and all bins are MINE):
  summary shape      mode height (1 cm bins), quantiles q05/q10/q25/q50/q75/q90/q95 above the RTK ground, IQR,
                     skewness, excess kurtosis, Shannon entropy of the 1 cm histogram, ground share (< +0.05 m),
                     canopy share (> +0.50 m, of all returns), median intensity, single-return share
  two-parameter      Weibull (shape k, scale lambda) fitted by maximum likelihood to the heights above -0.05 m
                     (shifted so that the fit domain starts at the lowest ground return; the shift is MINE)
  full shape         the 1 cm density on -0.05..+0.45 m (50 bins) per point; principal components across points
                     (functional PCA): the scores are the "full shape factor", the loadings say what shape varies.
Then the shape metrics by location bin (floodplain_offset_points.csv, cluster_med08_k4 -> south zone / open field /
swale stretch) and by the imagery types already assigned to the points (floodplain_types_points.csv: cls_naip_2021_k5,
cls_s2_4dates_k6), and the Spearman correlation of each metric with the NDVI metrics. Measurement only; no model.
Outputs  data/field/2026-10-25/returns/floodplain_return_shape_points.csv, _log.txt
         figures/rtk_2026-10-25_floodplain_return_shape.png   mean density with +/- PC1, PC2; PC1 vs PC2 scores by bin;
                                                               Weibull k vs lambda by bin
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_return_shape_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
SUB, GROUND, LO, HI = 0.50, 0.05, -0.05, 0.45; EDGES = np.arange(LO, HI + 0.005, 0.01); MID = 0.5 * (EDGES[1:] + EDGES[:-1])
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); OP = pd.read_csv(f"{D}/floodplain_offset_points.csv")
BIN = {0: "open field", 1: "south zone", 2: "south zone", 3: "swale stretch"}; OP["loc_bin"] = OP.cluster_med08_k4.map(BIN)
T = pd.read_csv("data/field/2026-10-25/imagery/floodplain_types_points.csv")[["group", "point", "cls_naip_2021_k5", "cls_s2_4dates_k6"]]
C = pd.read_csv("data/field/2026-10-25/imagery/floodplain_points_covariates.csv")[["group", "point", "metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "gen1_spread_mm"]]
rows, dens = [], []
for tag, group in (("rsplus_b", "floodplain_RSplus_B"), ("rsplus_a", "floodplain_RSplus_A"), ("cors", "floodplain_CORS_Wickert")):
    t = pd.read_csv(f"{D}/{tag}_returns_2021.csv.gz")
    for _, r in P[(P.group == group) & (P.status == "FIX")].iterrows():
        a = t[(t.point == r.point) & t.radius_3m]; v = a.dz_m.to_numpy(); s = v[v < SUB]
        h, _ = np.histogram(s, bins=EDGES); d = h / max(h.sum(), 1) / 0.01; dens.append(d)
        p = d / max(d.sum(), 1e-12); ent = -float(np.sum(p[p > 0] * np.log(p[p > 0])))
        try:
            k, loc, lam = stats.weibull_min.fit(s[s > LO] - LO, floc=0)
        except Exception:
            k, lam = np.nan, np.nan
        q = np.quantile(s, [.05, .1, .25, .5, .75, .9, .95])
        rows.append(dict(group=group, point=int(r.point), E=r.easting_tied, N=r.northing_tied, n3m=len(v), n_sub=len(s), mode_m=MID[np.argmax(h)], q05=q[0], q10=q[1], q25=q[2], q50=q[3], q75=q[4], q90=q[5], q95=q[6],
                         iqr=q[4] - q[2], skew=stats.skew(s), kurt=stats.kurtosis(s), entropy=ent, ground_share=(v < GROUND).mean(), canopy_share=(v > SUB).mean(),
                         intensity_med=a.intensity.median(), single_share=(a.number_of_returns == 1).mean(), weibull_k=k, weibull_lambda=lam))
X = pd.DataFrame(rows).merge(OP[["group", "point", "loc_bin", "med08_3m_m"]], on=["group", "point"]).merge(T, on=["group", "point"], how="left").merge(C, on=["group", "point"], how="left")
Dm = np.array(dens); mean = Dm.mean(0); U, S, Vt = np.linalg.svd(Dm - mean, full_matrices=False); var = S ** 2 / (S ** 2).sum()
for i in range(4): X[f"pc{i + 1}"] = (Dm - mean) @ Vt[i]
X.to_csv(f"{D}/floodplain_return_shape_points.csv", index=False, float_format="%.4f")
say(f"points {len(X)}; 2021 returns within 3 m per point median {int(X.n3m.median())}, sub-canopy (< {SUB} m) median {int(X.n_sub.median())}; density bins 0.01 m on {LO}..{HI} m")
say("functional PCA of the per-point densities, variance explained: " + ", ".join(f"PC{i + 1} {100 * var[i]:.1f} %" for i in range(5)))
say("PC1 loading (density change per unit score) at heights 0.00/0.05/0.10/0.15/0.20/0.25 m: " + ", ".join(f"{Vt[0][np.argmin(np.abs(MID - h))]:+.3f}" for h in (0, .05, .1, .15, .2, .25)))
say("PC2 loading at the same heights: " + ", ".join(f"{Vt[1][np.argmin(np.abs(MID - h))]:+.3f}" for h in (0, .05, .1, .15, .2, .25)))
M = ["mode_m", "q10", "q50", "q90", "iqr", "skew", "kurt", "entropy", "ground_share", "canopy_share", "intensity_med", "single_share", "weibull_k", "weibull_lambda", "pc1", "pc2"]
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
say("\nmedians by location bin:"); say(X.groupby("loc_bin")[M].median().round(3).T.to_string())
say("\nKruskal-Wallis across the three location bins, p per metric:"); say(" ".join(f"{m} {stats.kruskal(*[X.loc[X.loc_bin == b, m].dropna() for b in X.loc_bin.unique()]).pvalue:.3g};" for m in M))
for col in ("cls_naip_2021_k5", "cls_s2_4dates_k6"):
    g = X.dropna(subset=[col]); say(f"\nmedians by imagery type {col} (n per type {g[col].value_counts().sort_index().to_dict()}):"); say(g.groupby(col)[M].median().round(3).T.to_string())
    say("Kruskal-Wallis p: " + " ".join(f"{m} {stats.kruskal(*[g.loc[g[col] == c, m].dropna() for c in g[col].unique()]).pvalue:.3g};" for m in M))
say("\nSpearman rho of each shape metric with the imagery metrics (147 points; |rho| > 0.16 is p < 0.05 for n 147):")
rho = pd.DataFrame({c: [X[[m, c]].dropna().corr(method="spearman").iloc[0, 1] for m in M] for c in ("metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "med08_3m_m", "gen1_spread_mm")}, index=M)
say(rho.round(2).to_string())
within = X[X.loc_bin == "open field"]; say(f"\nwithin the open field only (n {len(within)}), Spearman with April NDVI: " + " ".join(f"{m} {within[[m, 'metric_ndvi_apr']].dropna().corr(method='spearman').iloc[0, 1]:+.2f};" for m in ("q50", "iqr", "ground_share", "intensity_med", "weibull_k", "pc1", "pc2")))
COL = {"south zone": "#d55e00", "open field": "#0072b2", "swale stretch": "#009e73"}
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
axs[0].plot(MID, mean, "k", lw=1.8, label="mean density")
for i, ls in ((0, "--"), (1, ":")):
    sd = X[f"pc{i + 1}"].std(); axs[0].plot(MID, mean + sd * Vt[i], color="#0072b2", ls=ls, lw=1.2, label=f"mean + 1 sd of PC{i + 1} ({100 * var[i]:.0f} %)"); axs[0].plot(MID, mean - sd * Vt[i], color="#d55e00", ls=ls, lw=1.2, label=f"mean - 1 sd of PC{i + 1}")
axs[0].axvline(0, color="0.4", lw=0.8, ls=":"); axs[0].set_xlabel("2021 return minus RTK ground (m)"); axs[0].set_ylabel("density (1/m)"); axs[0].legend(frameon=False, fontsize=7)
for b, c in COL.items():
    q = X[X.loc_bin == b]; axs[1].scatter(q.pc1, q.pc2, color=c, s=18, label=b, edgecolors="k", linewidths=0.3); axs[2].scatter(q.weibull_k, q.weibull_lambda, color=c, s=18, label=b, edgecolors="k", linewidths=0.3)
axs[1].set_xlabel(f"PC1 score ({100 * var[0]:.0f} %)"); axs[1].set_ylabel(f"PC2 score ({100 * var[1]:.0f} %)"); axs[1].legend(frameon=False, fontsize=8)
axs[2].set_xlabel("Weibull shape k"); axs[2].set_ylabel("Weibull scale (m above -0.05 m)")
for ax in axs: ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_return_shape.png", dpi=int(os.environ.get("FIGDPI", 170))); say("wrote figures/rtk_2026-10-25_floodplain_return_shape.png")
LOG.close()
