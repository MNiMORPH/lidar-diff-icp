#!/usr/bin/env python3
"""One field day (Andy, 2026-10-20, D5): "We have one day in the field just now ... find us
the minimal set of places, ideally accessible and close, where we can survey RCG as well as
other land-cover types."

Inputs: the valley post-hoc clusters (data/derived/valley/valley_sites.csv: classes R1 NE-like
2025+2026, R2 one year, S late-green not NE-like, C bare, E forest edge; each with the
continuous values p NE-like per year, spring / October NDVI, largest 2026 drop); the Elba
lidar sites (data/derived/field_sites_rtk.csv: F1 sedge, F2 closed canopy, F3 bare, H2 steep
gen1-low, H3 open hillslope); WMA parking lots and roads; the 2021 project control marks on
public land (the datum tie); the class raster and MnGeo 2025 FSA imagery for the map.

Method: for each WMA parking lot, the sites within a WALKING RADIUS (straight line; the
radius is MINE, reported at 600 and 1000 m): the distinct cover types reached, whether a
"sward in both lidar epochs" cluster (NE-like >= 50% in 2008 and 2021; 50% mine) is among
them, and whether a control mark is within reach. Lots ranked by types reached, then by the
total walking distance to one site of each type. The minimal set = the best lot, plus the
nearest second lot only for a type the first cannot reach.

    ./lidar-icp/bin/python analysis/rcg/field_day_plan.py
Outputs data/derived/valley/field_day_sites.csv, figures/field_day.png
"""
import io
import json
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import rasterio
import shapely
from matplotlib.colors import ListedColormap
from PIL import Image
from shapely.geometry import shape

V = "data/derived/valley"; PL = f"{V}/public_land"
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
S = pd.read_csv(f"{V}/valley_sites.csv"); S["type"] = S.cls.str.split(" ").str[0]; S["src"] = "valley imagery"; S["site"] = [f"{t}-v{i}" for i, t in zip(S.index, S.type)]
S["both_epochs"] = (S.nelike_2008 >= 50) & (S.nelike_2021 >= 50)
E_ = pd.read_csv("data/derived/field_sites_rtk.csv"); E_["type"] = E_.cls.str.split(" ").str[0]; E_["src"] = "Elba lidar"; E_["both_epochs"] = False
A = pd.concat([S, E_], ignore_index=True, sort=False)
feats = lambda n: json.load(open(f"{PL}/{n}"))["features"]
lots = [(shape(f["geometry"]), f["properties"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"]
units = [shape(f["geometry"]) for f in feats("wma_units.geojson")] + [shape(f["geometry"]) for f in feats("state_parks.geojson")]
public = shapely.union_all(units); roads = [shape(f["geometry"]) for f in feats("wma_roads.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]
M = pd.read_csv("data/derived/control_marks_gen2_2021_control_LCP-NVA-VVA_marks.csv")
M = M[[public.contains(shapely.Point(x, y)) for x, y in zip(M.easting, M.northing)]].copy()
print(f"inputs: {len(S)} valley clusters, {len(E_)} Elba lidar sites, {len(lots)} WMA parking lots, {len(M)} control marks on public land: {M.point_id.tolist()}")
TYPES = ["R1", "R2", "S", "C", "E", "F1", "F2", "F3", "H2", "H3"]

def reach(lot, R):
    d = np.hypot(A.easting - lot.x, A.northing - lot.y); k = A[d <= R].assign(dist=d[d <= R])
    dm = np.hypot(M.easting - lot.x, M.northing - lot.y)
    return k, (M.point_id[dm <= R].tolist()), dm.min() if len(M) else np.nan

rows = []
for R in (600, 1000):
    for g, p in lots:
        k, marks, dmark = reach(g, R)
        types = sorted(set(k.type)); nearest = k.sort_values("dist").groupby("type").head(1)
        rows.append(dict(R=R, lot_id=p["objectid"], easting=round(g.x, 1), northing=round(g.y, 1), unit=p.get("habitat_mgmt_unit_name"), n_types=len(types), types=" ".join(types),
                         both_epochs=int(k.both_epochs.any()), marks_in_reach=" ".join(marks), nearest_mark_m=round(dmark, 0), walk_sum_m=round(nearest.dist.sum(), 0), n_sites=len(k)))
L = pd.DataFrame(rows).sort_values(["R", "n_types", "both_epochs", "walk_sum_m"], ascending=[True, False, False, True])
pd.set_option("display.width", 250)
for R in (600, 1000):
    print(f"\nPARKING LOTS ranked, walking radius {R} m (straight line):"); print(L[L.R == R].head(6).drop(columns="R").to_string(index=False))

# the minimal set at 1000 m: best lot, plus the nearest lot that adds a type the best cannot reach
best = L[L.R == 1000].iloc[0]; have = set(best.types.split())
k1, marks1, _ = reach(shapely.Point(best.easting, best.northing), 1000)
missing = [t for t in TYPES if t not in have and t in set(A.type)]
add = None
for _, r in L[L.R == 1000].iterrows():
    if r.lot_id == best.lot_id: continue
    new = set(r.types.split()) - have
    if new:
        add = r; break
print(f"\nMINIMAL SET: lot {best.lot_id} at E {best.easting} N {best.northing} ({best.unit}) reaches {best.n_types} types [{best.types}] within 1000 m, both-epoch sward {bool(best.both_epochs)}, marks {marks1}")
print(f"  types not reachable from it: {missing}" + (f"; the nearest lot adding some: lot {add.lot_id} at E {add.easting} N {add.northing}, adds {sorted(set(add.types.split()) - have)}, {np.hypot(add.easting - best.easting, add.northing - best.northing) / 1000:.1f} km away" if add is not None else ""))

chosen = [(best.lot_id, best.easting, best.northing)] + ([(add.lot_id, add.easting, add.northing)] if add is not None else [])
out = []
for lid, ex, ny_ in chosen:
    k, marks, _ = reach(shapely.Point(ex, ny_), 1000)
    for t, grp in k.sort_values("dist").groupby("type"):
        for _, r in grp.head(2).iterrows():
            out.append(dict(lot=lid, type=t, site=r.site, src=r.src, dist_m=round(r.dist, 0), ha=r.get("ha"), easting=r.easting, northing=r.northing, lat=round(r.lat, 6), lon=round(r.lon, 6),
                            nelike_2008=r.get("nelike_2008"), nelike_2021=r.get("nelike_2021"), nelike_2025=r.get("nelike_2025"), nelike_2026=r.get("nelike_2026"), both_epochs=r.both_epochs,
                            spring_ndvi=r.get("spring2026_ndvi"), oct_ndvi=r.get("oct2026_ndvi"), drop_max=r.get("drop_max"), drop_doy=r.get("drop_max_doy"), raw_dod_mm=r.get("raw_dod_mm"), gen1_span_mm=r.get("gen1_span_mm"), slope_deg=r.get("slope_deg")))
O = pd.DataFrame(out); O.to_csv(f"{V}/field_day_sites.csv", index=False)
print("\nFIELD-DAY SITES (nearest two of each type per lot, within 1000 m):"); print(O.to_string(index=False))
print("\nmarks on public land:"); print(M[["point_id", "easting", "northing", "point_type"]].to_string(index=False))

# map: 2025 FSA imagery (MnGeo WMS) under the class raster, around the chosen lot(s)
ex = [c[1] for c in chosen]; ny_ = [c[2] for c in chosen]; cx, cy = np.mean(ex), np.mean(ny_); half = max(1400, (max(ex) - min(ex)) / 2 + 1200, (max(ny_) - min(ny_)) / 2 + 1200)
X0, X1, Y0, Y1 = cx - half, cx + half, cy - half, cy + half; W = int(2 * half / 2)
url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH={W}&HEIGHT={W}&FORMAT=image/jpeg"
try: img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
except Exception as e: img = None; print("imagery not fetched:", e)
with rasterio.open(f"{V}/valley_classes.tif") as r:
    win = rasterio.windows.from_bounds(X0, Y0, X1, Y1, r.transform); cls = r.read(1, window=win)
fig, ax = plt.subplots(figsize=(13, 13))
if img is not None: ax.imshow(img, extent=[X0, X1, Y0, Y1])
cm = ListedColormap(["#8b0000", "#eb6834", "#4a3aa7", "#eda100", "#1baf7a"])
ax.imshow(np.where(cls > 0, cls, np.nan), extent=[X0, X1, Y0, Y1], cmap=cm, vmin=0.5, vmax=5.5, alpha=0.45, interpolation="nearest")
for g in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=1)
for rd in roads: ax.plot(*rd.xy, color="white", lw=1.2, ls="--")
for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=1.2)
for g, p in lots: ax.plot(g.x, g.y, marker="P", color="yellow", ms=12, mec="black")
for lid, x, y in chosen: ax.add_patch(plt.Circle((x, y), 1000, fill=False, color="yellow", lw=1.5, ls=":")); ax.annotate(f"lot {lid}", (x, y), xytext=(6, 6), textcoords="offset points", color="yellow", fontsize=10, weight="bold")
for _, r in O.iterrows():
    ax.plot(r.easting, r.northing, marker="o", ms=10, mfc="none", mec="white", mew=2); ax.annotate(f"{r.type}", (r.easting, r.northing), xytext=(5, 5), textcoords="offset points", color="white", fontsize=9, weight="bold")
for _, r in M.iterrows():
    if X0 < r.easting < X1 and Y0 < r.northing < Y1: ax.plot(r.easting, r.northing, marker="^", ms=12, color="cyan", mec="black"); ax.annotate(r.point_id, (r.easting, r.northing), xytext=(6, -12), textcoords="offset points", color="cyan", fontsize=9)
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xlabel("E (m)"); ax.set_ylabel("N (m)")
ax.set_title("One field day: 2025 FSA imagery; classes R1 dark red, R2 orange, S purple, C yellow, E green (post-hoc view); yellow P = WMA parking, dotted = 1 km; circles = sites; cyan = 2021 control marks", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_day.png", dpi=100); print("figure: figures/field_day.png; table: data/derived/valley/field_day_sites.csv")
