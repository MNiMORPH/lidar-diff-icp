#!/usr/bin/env python
"""The hillslope RTK points: 2021 returns within 1 m (and 3 m) of each FIX point of the two RS3 traverses, as the
slope-normal distance from the surveyed ground (Andy, 2026-11-21: the bluff counterpart of the floodplain's 1 m stacks).
Reference per point: a plane through the tied RTK point with the official 1 m DEM's local gradient (dem3dep/
usgs_1m_dem_bluff_wide.tif, central differences); height = (z - plane(x, y)) x cos(slope). Pooled over the points, drawn
as the bluff-wide figure was: symmetric-log height axis ground to canopy top with intensity-colored bars, a
near-ground zoom with the DEM-referenced bluff curve, and the two teams separately. Log: shares by height band,
intensity by band, per-point ground share, DEM minus RTK.
Outputs  data/field/2026-10-25/intensity/hillslope_points_log.txt; figures/rtk_2026-10-25_hillslope_points_returns.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import map_coordinates
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/hillslope_points_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_bluff_wide.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan; tr = r.transform; E0, N1 = tr.c, tr.f
gy, gx = np.gradient(DEM, 1.0); gy = -gy
def at(a, x, y): return map_coordinates(a, [((N1 - y) - 0.5), ((x - E0) - 0.5)], order=1, mode="nearest")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy()
assert (H.easting_tied.between(E0 + 2, E0 + DEM.shape[1] - 2) & H.northing_tied.between(N1 - DEM.shape[0] + 2, N1 - 2)).all(), "a hillslope point lies outside the DEM window"
H["gx"] = at(gx, H.easting_tied.values, H.northing_tied.values); H["gy"] = at(gy, H.easting_tied.values, H.northing_tied.values); H["slope"] = np.degrees(np.arctan(np.hypot(H.gx, H.gy))); H["dem"] = at(DEM, H.easting_tied.values, H.northing_tied.values)
say(f"hillslope FIX points: {len(H)} ({H.group.value_counts().to_dict()}); DEM slope at the points median {H.slope.median():.1f} deg; official DEM minus RTK (vertical) median {(H.dem - H.navd88_tied_m).median():+.3f}, slope-normal median {((H.dem - H.navd88_tied_m) * np.cos(np.radians(H.slope))).median():+.3f}")
bbox = (H.easting_tied.min() - 4, H.easting_tied.max() + 4, H.northing_tied.min() - 4, H.northing_tied.max() + 4); parts = []
with laspy.open("data/after/valley/4342-29-64_3dep_full.laz") as r:
    for ch in r.chunk_iterator(5_000_000):
        x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])
        if m.any(): parts.append(pd.DataFrame(dict(x=x[m], y=y[m], z=np.asarray(ch.z)[m], inten=np.asarray(ch.intensity)[m].astype(float), nr=np.asarray(ch.number_of_returns)[m], rn=np.asarray(ch.return_number)[m])))
G = pd.concat(parts, ignore_index=True); rows = []
for _, p in H.iterrows():
    d = np.hypot(G.x - p.easting_tied, G.y - p.northing_tied); w = G[d <= 3].copy(); w["dist"] = d[d <= 3]
    plane = p.navd88_tied_m + p.gx * (w.x - p.easting_tied) + p.gy * (w.y - p.northing_tied); w["h"] = (w.z - plane) * np.cos(np.radians(p.slope)); w["group"] = p.group; w["point"] = p.point; w["slope"] = p.slope; rows.append(w)
T = pd.concat(rows, ignore_index=True); T1 = T[T.dist <= 1]
def bands(q, lab):
    say(f"{lab}: n {len(q):,}; within +-0.05 m {100 * (q.h.abs() < 0.05).mean():.1f} %, below -0.05 {100 * (q.h < -0.05).mean():.1f} %, 0.05-0.5 m {100 * ((q.h >= 0.05) & (q.h < 0.5)).mean():.1f} %, 0.5-2 m {100 * ((q.h >= 0.5) & (q.h < 2)).mean():.1f} %, above 2 m {100 * (q.h >= 2).mean():.1f} %; "
        f"near-ground median {q.h[q.h.abs() < 0.5].median():+.3f}, NMAD {1.4826 * np.median(np.abs(q.h[q.h.abs() < 0.5] - q.h[q.h.abs() < 0.5].median())):.3f}; intensity: ground {q.inten[q.h.abs() < 0.05].median():.0f}, understory {q.inten[(q.h >= 0.05) & (q.h < 0.5)].median():.0f}, canopy {q.inten[q.h >= 2].median():.0f}; single {100 * (q.nr == 1).mean():.1f} %")
bands(T1, "within 1 m, both teams"); bands(T[T.dist <= 3], "within 3 m, both teams")
for g, q in T1.groupby("group"): bands(q, f"within 1 m, {g}")
pp = T1.groupby(["group", "point"]).apply(lambda q: pd.Series(dict(n=len(q), ground_share=(q.h.abs() < 0.05).mean(), canopy_share=(q.h >= 2).mean(), med_near=q.h[q.h.abs() < 0.5].median() if (q.h.abs() < 0.5).any() else np.nan)))
say("per point (1 m): returns median %d (%d-%d); ground share median %.2f, p10 %.2f, p90 %.2f; points with no return within 5 cm of the ground: %d of %d; near-ground median offset over points %+.3f m" % (pp.n.median(), pp.n.min(), pp.n.max(), pp.ground_share.median(), pp.ground_share.quantile(.1), pp.ground_share.quantile(.9), int((pp.ground_share == 0).sum()), len(pp), pp.med_near.median()))
F = pd.read_csv(f"{OUT}/hillslope_returns_log.txt", sep="\x00", header=None) if False else None
norm = Normalize(vmin=1000, vmax=4500); cmap = plt.get_cmap("viridis")
fig, axs = plt.subplots(1, 3, figsize=(21, 5.5), gridspec_kw=dict(width_ratios=[1.4, 1, 1]))
edges = np.concatenate([np.arange(-0.60, 0.30, 0.01), 0.30 * 1.15 ** np.arange(0, 40)]); edges = edges[edges <= T1.h.max() * 1.15]; mids = 0.5 * (edges[1:] + edges[:-1])
cnt = np.histogram(T1.h, bins=edges)[0]; dens = cnt / (len(T1) * np.diff(edges)); b = np.digitize(T1.h, edges) - 1; med = np.array([T1.inten[b == k].median() if (b == k).sum() >= 10 else np.nan for k in range(len(mids))])
ax = axs[0]; ax.bar(mids, np.where(dens > 0, dens, np.nan), width=np.diff(edges), color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in med], edgecolor="none", align="center"); ax.set_xscale("symlog", linthresh=0.3, linscale=1.0); ax.set_yscale("log"); ax.set_xlim(-0.6, edges[-1])
ax.set_xticks([-0.5, -0.2, 0, 0.2, 0.5, 1, 2, 5, 10, 20, 30]); ax.set_xticklabels(["-0.5", "-0.2", "0", "0.2", "0.5", "1", "2", "5", "10", "20", "30"]); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlabel("2021 return, slope-normal distance from the surveyed ground (m); symmetric-log axis"); ax.set_ylabel("density (1/m), log"); ax.set_title(f"within 1 m of the {len(H)} hillslope points (n {len(T1):,})", fontsize=10)
# near-ground zoom with the bluff-wide DEM-referenced curve recomputed from the bluff script's inputs is not stored; use the 3 m stack as the comparison
ax = axs[1]; e2 = np.arange(-0.5, 1.0 + 0.005, 0.01); m2 = 0.5 * (e2[1:] + e2[:-1]); c2 = np.histogram(T1.h, bins=e2)[0]; d2 = c2 / (len(T1) * 0.01); b2 = np.digitize(T1.h, e2) - 1; md2 = np.array([T1.inten[b2 == k].median() if (b2 == k).sum() >= 10 else np.nan for k in range(len(m2))])
ax.bar(m2, np.where(d2 > 0, d2, np.nan), width=0.01, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in md2], edgecolor="none"); ax.hist(T[T.dist <= 3].h, bins=e2, density=True, histtype="step", color="k", lw=1.0, label=f"within 3 m (n {int((T.dist <= 3).sum()):,})")
ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.5, 1.0); ax.set_xlabel("slope-normal distance from the surveyed ground (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=8); ax.set_title("near the ground", fontsize=10)
ax = axs[2]
for g, c in (("hillslope_RS3", "#0072b2"), ("hillslope_T4a", "#d55e00")):
    q = T1[T1.group == g]; ax.hist(q.h, bins=e2, density=True, histtype="step", color=c, lw=1.4, label=f"{g.replace('hillslope_', 'RS3 ')} (n {len(q):,}, {q.point.nunique()} points)")
ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.5, 1.0); ax.set_xlabel("slope-normal distance from the surveyed ground (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=8); ax.set_title("the two traverses, within 1 m", fontsize=10)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs[:2], shrink=0.8, label="median intensity in the bar")
for a in axs: a.spines[["top", "right"]].set_visible(False)
fig.savefig("figures/rtk_2026-10-25_hillslope_points_returns.png", dpi=int(os.environ.get("FIGDPI", 160)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_hillslope_points_returns.png")
LOG.close()
