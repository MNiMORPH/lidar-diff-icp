#!/usr/bin/env python3
"""Join the lidar site screen of every finished valley tile with the valley-wide imagery
layers, and rank within class (a POST-HOC VIEW under D4: every cut-off is printed).

Inputs: data/derived/valley/sites_valley_lidar.csv (analysis/valley_tile_sites.py --all,
re-run here first so new tiles are in); the Phase A continuous layers on the 10 m grid
(p_nelike_{2008,2021,2025,2026}.tif, spring_ndvi_2026.tif, oct_ndvi_2026.tif,
ndvi_drop_max_2026.tif); the public-land layers. Each lidar site is a cluster of 5 m cells
reported by its centroid; the imagery values are the median over the 3x3 window of 10 m
pixels around that centroid (MINE: the screener does not export the cluster masks).
Rank score within class = ha / (1 + km to the nearest WMA parking lot), as in Phase A (MINE).
Outputs: data/derived/valley/sites_valley_ranked.csv, figures/valley_lidar_sites.png.
"""
import json
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import shapely
from shapely.geometry import shape

V = "data/derived/valley"
subprocess.run(["env", "-u", "PROJ_DATA", "-u", "GDAL_DATA", "./lidar-icp/bin/python", "analysis/valley_tile_sites.py", "--all"], check=True, capture_output=True)
S = pd.read_csv(f"{V}/sites_valley_lidar.csv"); S["type"] = S.cls.str.split(" ").str[0]
print(f"lidar sites from the screener: {len(S)} on {S.tile.nunique()} finished tiles: {sorted(S.tile.unique())}")
def samp(path, band=1):
    with rasterio.open(path) as r:
        a = r.read(band); tr = r.transform; out = []
        for x, y in zip(S.easting, S.northing):
            c, rr = ~tr * (x, y); c, rr = int(c), int(rr); w = a[max(rr - 1, 0):rr + 2, max(c - 1, 0):c + 2]
            out.append(float(np.nanmedian(w)) if np.isfinite(w).any() else np.nan)
    return np.array(out)
for y in (2008, 2021, 2025, 2026): S[f"nelike_{y}"] = (100 * samp(f"{V}/p_nelike_{y}.tif")).round(0)
S["spring2026_ndvi"] = samp(f"{V}/spring_ndvi_2026.tif").round(2); S["oct2026_ndvi"] = samp(f"{V}/oct_ndvi_2026.tif").round(2); S["drop_max_2026"] = samp(f"{V}/ndvi_drop_max_2026.tif", 1).round(2); S["drop_doy_2026"] = samp(f"{V}/ndvi_drop_max_2026.tif", 2)
S["both_epochs_nelike"] = (S.nelike_2008 >= 50) & (S.nelike_2021 >= 50)
# WorldCover 2021 class at the centroid (mode of the 3x3 10 m window): the only present-day
# cover source off the imagery classes; it is what separates a 2008 standing crop from sward
with rasterio.open(f"{V}/worldcover.tif") as r:
    wc = r.read(1); trw = r.transform; names = {10: "tree", 20: "shrub", 30: "grass", 40: "crop", 50: "built", 60: "bare", 80: "water", 90: "wetland"}; cls_ = []
    for x, y in zip(S.easting, S.northing):
        c, rr = ~trw * (x, y); c, rr = int(c), int(rr); w = wc[max(rr - 1, 0):rr + 2, max(c - 1, 0):c + 2].ravel(); v, n = np.unique(w, return_counts=True); cls_.append(names.get(int(v[n.argmax()]), str(v[n.argmax()])))
S["worldcover"] = cls_
S["rank_score"] = S.ha / (1 + S.dist_parking_m / 1000.0); S = S.sort_values(["type", "rank_score"], ascending=[True, False])
S.to_csv(f"{V}/sites_valley_ranked.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
cols = ["tile", "site", "ha", "lat", "lon", "worldcover", "raw_dod_mm", "gen1_span_mm", "slope_deg", "penfail_pct", "lowveg", "fingerprint_score", "nelike_2008", "nelike_2021", "nelike_2026", "spring2026_ndvi", "drop_max_2026", "dist_parking_m", "dist_wma_road_m"]
print("\nimagery values = median of the 3x3 10 m window at the cluster centroid (MINE); nelike = NE-like probability %, per year; rank = ha / (1 + km to WMA parking) (MINE)")
for t, g in S.groupby("type", sort=True):
    print(f"\n{g.cls.iloc[0]}: {len(g)} sites; top 5:"); print(g[cols].head(5).to_string(index=False))
# map
feats = lambda n: json.load(open(f"{V}/public_land/{n}"))["features"]
units = [shape(f["geometry"]) for f in feats("wma_units.geojson")] + [shape(f["geometry"]) for f in feats("state_parks.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]
lots = [shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"]
fig, ax = plt.subplots(figsize=(11, 12)); col = {"F1": "#2a78d6", "F2": "#eb6834", "F3": "#1baf7a", "H2": "#e34948", "H3": "#eda100", "H4": "#4a3aa7"}
for g in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=0.7)
for r in rivers: ax.plot(*r.xy, color="#2a78d6", lw=0.7, alpha=0.6)
for p in lots: ax.plot(p.x, p.y, marker="P", color="yellow", ms=6, mec="black")
for t, g in S.groupby("type"): ax.scatter(g.easting, g.northing, s=18 + 4 * g.ha, c=col.get(t, "grey"), edgecolors="black", linewidths=0.4, label=f"{g.cls.iloc[0]} ({len(g)})", zorder=3)
tiles = json.load(open(f"{V}/public_tiles.json"))
for t in S.tile.unique():
    b = tiles[t]["bounds"]
    if b: ax.add_patch(plt.Rectangle((b[0], b[1]), b[2] - b[0], b[3] - b[1], fill=False, ec="grey", lw=0.8, ls=":"))
xs = S.easting; ys = S.northing; ax.set_xlim(xs.min() - 1500, xs.max() + 1500); ax.set_ylim(ys.min() - 1500, ys.max() + 1500); ax.set_aspect("equal"); ax.legend(loc="lower left", fontsize=8); ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"Lidar-screened RTK sites on public land, finished valley tiles ({S.tile.nunique()}; dotted = tile bounds); marker size ~ ha; yellow = WMA parking", fontsize=9)
fig.tight_layout(); fig.savefig("figures/valley_lidar_sites.png", dpi=110); print("\nfigures/valley_lidar_sites.png; table data/derived/valley/sites_valley_ranked.csv")
