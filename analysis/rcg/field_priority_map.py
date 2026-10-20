#!/usr/bin/env python3
"""Priority map: the ten survey items numbered in priority order with one-to-three-word names
(Andy, 2026-10-21). Built on the route map's layers; FIELD_SITES_RTK.md carries the reasons.
Writes figures/field_priority.png.
"""
import io
import json
import re
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import shapely
from PIL import Image
from shapely.geometry import shape

V = "data/derived/valley"; from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
feats = lambda n: json.load(open(f"{V}/public_land/{n}"))["features"]
lots = {f["properties"]["objectid"]: shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"}
roads = [shape(f["geometry"]) for f in feats("wma_roads.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]; units = [shape(f["geometry"]) for f in feats("wma_units.geojson")]
def gpx_polys(path):
    root = ET.parse(path).getroot(); ns = {"g": "http://www.topografix.com/GPX/1/1"}; out = {}
    for t in root.findall("g:trk", ns):
        n = t.find("g:name", ns).text; pts = [from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in t.findall(".//g:trkpt", ns)]
        if len(pts) >= 4: out.setdefault(re.sub(r"_\d+$", "", n), []).append(shapely.Polygon(pts))
    return {k: shapely.union_all(v) for k, v in out.items()}
R = {k.replace("region_", ""): v for k, v in gpx_polys(f"{V}/field_day_waypoints.gpx").items() if k.startswith("region_")}
B = {k: v for k, v in gpx_polys(f"{V}/best_places.gpx").items() if k.startswith("1_")}; north = shapely.union_all(list(B.values()))
tx = np.array([from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in ET.parse(f"{V}/field_day_waypoints.gpx").getroot().iter("{http://www.topografix.com/GPX/1/1}trkpt")][:61])
def rp(k): g = next(v for kk, v in R.items() if kk.startswith(k)); return g.representative_point()
ITEMS = [(1, "Base + mark", (579637.0, 4885872.1)), (2, "Uncut sward", rp("R1-v140")), (3, "Bare control", rp("C-v52")), (4, "North sward", north.representative_point()),
         (5, "Transect", shapely.Point(tx[len(tx) // 2])), (6, "Cut sward", rp("R2-v237")), (7, "Meadow", rp("S-v324")), (8, "Zero penetration", rp("F2-1")), (8, "Zero penetration", rp("F2-3")),
         (9, "Slope patches", rp("gen2-high")), (9, "Slope control", rp("hillslope control")), (10, "Bluff stand", rp("bluff stand"))]
bb = (578700, 4884500, 580900, 4887900); X0, Y0, X1, Y1 = bb
url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH=1400&HEIGHT={int(1400 * (Y1 - Y0) / (X1 - X0))}&FORMAT=image/jpeg"
try: img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
except Exception as e: img = None; print("imagery not fetched:", e)
fig, ax = plt.subplots(figsize=(10, 10 * (Y1 - Y0) / (X1 - X0) + 0.6))
if img is not None: ax.imshow(img, extent=[X0, X1, Y0, Y1])
for g in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=0.8)
for rd in roads: ax.plot(*rd.xy, color="white", lw=1.5, ls="--")
for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=1.2)
for g in list(R.values()) + [north]:
    for pg in (g.geoms if hasattr(g, "geoms") else [g]): ax.fill(*pg.exterior.xy, color="white", alpha=0.2); ax.plot(*pg.exterior.xy, color="white", lw=1.2)
ax.plot(tx[:, 0], tx[:, 1], color="#eda100", lw=2.5)
for lot in (501, 2597, 204): L = lots[lot]; ax.plot(L.x, L.y, marker="P", ms=14, color="yellow", mec="black"); ax.annotate(f"lot {lot}", (L.x, L.y), xytext=(8, -14), textcoords="offset points", color="yellow", fontsize=11, weight="bold", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])
for n, name, p in ITEMS:
    x, y = (p.x, p.y) if hasattr(p, "x") else p
    ax.annotate(f"{n}", (x, y), ha="center", va="center", fontsize=13, weight="bold", color="black", bbox=dict(boxstyle="circle,pad=0.3", fc="white", ec="black", lw=1.5))
    ax.annotate(name, (x, y), xytext=(16, 0), textcoords="offset points", va="center", fontsize=11, weight="bold", color="white", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])
ax.plot([X0 + 100, X0 + 600], [Y0 + 100, Y0 + 100], color="white", lw=4); ax.annotate("500 m", (X0 + 350, Y0 + 130), color="white", ha="center", fontsize=11, weight="bold")
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([]); ax.set_title("Priority order", fontsize=12)
fig.tight_layout(); fig.savefig("figures/field_priority.png", dpi=110); print("figures/field_priority.png")
