#!/usr/bin/env python3
"""The field-day route: cover every terrain type with the least driving and parking (Andy,
2026-10-20: "Generate a plan to cover the terrain and minimize driving and parking. Build one
or more maps around these.").

Everything worth surveying lies along ONE road, the WMA road up the Main Branch north of
Elba, within 2.6 km: three parking lots, the base monument at the DNR headquarters junction
between them, and the survey regions. Stops, south to north:
  lot 204   R1-v140 the 11.5 ha RCG-type field (sward every year, uncut 2026), F2-1 inside
            it, R1-v149 beside the lot, the bluff lowveg stand and its gen2-high patches
  lot 501   base on DR9925 (310 m), mark 1079, C-v52 bare control, S-v324 meadow, R2-v237
            sward in 2008+2021 cut since, F2-3, the class transect
  lot 2597  best place 1: the 14 ha zero-penetration field, sward in every year, 283 m
            from the lot
Walking distances are exact (shapely) from each lot to the nearest edge of each region.
Driving legs are straight-line along the valley (the DNR road layer does not carry the
public road the lots sit on, so no routed length is claimed).
Outputs figures/field_route.png (overview on 2025 imagery), data/derived/valley/field_route.csv
(stop, region, walk m), field_route.gpx (the lots in order plus every region outline).
"""
import io
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import shapely
from PIL import Image
from shapely.geometry import shape

V = "data/derived/valley"; to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True); from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
feats = lambda n: json.load(open(f"{V}/public_land/{n}"))["features"]
lots = {f["properties"]["objectid"]: shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"}
roads = [shape(f["geometry"]) for f in feats("wma_roads.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]; units = [shape(f["geometry"]) for f in feats("wma_units.geojson")]
def gpx_polys(path, prefix=None):
    root = ET.parse(path).getroot(); ns = {"g": "http://www.topografix.com/GPX/1/1"}; out = {}
    for t in root.findall("g:trk", ns):
        n = t.find("g:name", ns).text
        if prefix is not None and not n.startswith(prefix): continue
        pts = [from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in t.findall(".//g:trkpt", ns)]
        if len(pts) >= 4: out.setdefault(re.sub(r"_\d+$", "", n) if prefix is None else "best1", []).append(shapely.Polygon(pts))
    return {k: shapely.union_all(v) for k, v in out.items()}
regions = {k.replace("region_", ""): v for k, v in gpx_polys(f"{V}/field_day_waypoints.gpx").items() if k.startswith("region_")}
regions["best place 1 (14 ha field)"] = gpx_polys(f"{V}/best_places.gpx", prefix="1_")["best1"]
transect = [t for t in ET.parse(f"{V}/field_day_waypoints.gpx").getroot().findall("{http://www.topografix.com/GPX/1/1}trk") if t.find("{http://www.topografix.com/GPX/1/1}name").text == "class_transect"][0]
tx = np.array([from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in transect.findall(".//{http://www.topografix.com/GPX/1/1}trkpt")])
# THE BASE STATION goes ON monument DR9925 (MnDOT "8508 Y 2023", 3/4-inch aluminium rod in a
# PVC access cover, 11.5 m SW of TH 74 at the DNR headquarters junction). Enter in Emlid Flow
# as a KNOWN base position, datum NAD 83(2011) epoch 2010.00, from the NGS datasheet (read
# 2026-10-20): 44 07 18.13863 N, 092 00 16.79318 W, ellipsoid height 188.669 m (network
# accuracy 0.61 cm horizontal, 1.49 cm height at 95%); NAVD 88 218.75 m is the GEOID18
# orthometric height of the same point, for the check. Antenna height is added on site.
BASE_LAT = 44 + 7 / 60 + 18.13863 / 3600; BASE_LON = -(92 + 0 / 60 + 16.79318 / 3600); BASE_H_ELL = 188.669
BASE_E, BASE_N = from_ll.transform(BASE_LON, BASE_LAT)
mons = {"DR9925 base": (BASE_E, BASE_N), "DG8385 check": (579729.0, 4885711.0), "mark 1079": (579659.8, 4885857.1)}
print(f"BASE on DR9925: {BASE_LAT:.8f} N {-BASE_LON:.8f} W, ellipsoid height {BASE_H_ELL} m, NAD 83(2011) epoch 2010.00; UTM 15N E {BASE_E:.2f} N {BASE_N:.2f}")
STOPS = [(204, ["R1-v140: RCG-type field, both epochs", "R1-v149", "F2-1: zero penetration", "bluff stand (green <= 25 deg)", "gen2-high patches", "hillslope control"]),
         (501, ["DR9925 base", "mark 1079", "DG8385 check", "C-v52: bare control", "S-v324: meadow", "R2-v237: sward 2008+2021, cut since", "F2-3: zero penetration", "transect"]),
         (2597, ["best place 1 (14 ha field)"])]
rows = []
for lot, items in STOPS:
    L = lots[lot]
    for it in items:
        if it in mons: d = np.hypot(mons[it][0] - L.x, mons[it][1] - L.y)
        elif it == "transect": d = shapely.LineString(tx).distance(L)
        else:
            key = next((k for k in regions if k.startswith(it.split(":")[0])), None); d = regions[key].distance(L) if key else np.nan
        rows.append(dict(stop=lot, item=it, walk_m=round(float(d))))
T = pd.DataFrame(rows); T.to_csv(f"{V}/field_route.csv", index=False); pd.set_option("display.width", 200); print(T.to_string(index=False))
legs = [(204, 501), (501, 2597)]; print("driving legs (straight line along the valley):", [(a, b, round(np.hypot(lots[a].x - lots[b].x, lots[a].y - lots[b].y) / 1000, 2)) for a, b in legs])
for lot in (204, 501, 2597): lo, la = to_ll.transform(lots[lot].x, lots[lot].y); print(f"lot {lot}: {la:.6f} N {-lo:.6f} W")
# map
bb = (578700, 4884500, 580900, 4887900); X0, Y0, X1, Y1 = bb
url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH=1400&HEIGHT={int(1400 * (Y1 - Y0) / (X1 - X0))}&FORMAT=image/jpeg"
try: img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
except Exception as e: img = None; print("imagery not fetched:", e)
fig, ax = plt.subplots(figsize=(9, 9 * (Y1 - Y0) / (X1 - X0) + 0.6))
if img is not None: ax.imshow(img, extent=[X0, X1, Y0, Y1])
for g in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=0.8)
for rd in roads: ax.plot(*rd.xy, color="white", lw=1.5, ls="--")
for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=1.2)
pe = [__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")]
col = {"R1": "white", "F2": "#eb6834", "R2": "#eb6834", "C": "#eda100", "S": "#4a3aa7", "E": "#1baf7a", "bluff stand (green": "#1baf7a", "steeper": "#e34948", "gen2-high": "#e87ba4", "hillslope control": "#eda100", "best": "white"}
for k, g in regions.items():
    c = next((v for kk, v in col.items() if k.startswith(kk)), "white")
    for pg in (g.geoms if hasattr(g, "geoms") else [g]): ax.fill(*pg.exterior.xy, color=c, alpha=0.25); ax.plot(*pg.exterior.xy, color=c, lw=1.6)
ax.plot(tx[:, 0], tx[:, 1], color="#eda100", lw=2.5)
for (lot, _), n in zip(STOPS, ("1", "2", "3")):
    L = lots[lot]; ax.plot(L.x, L.y, marker="P", ms=16, color="yellow", mec="black"); ax.annotate(f"stop {n}: lot {lot}", (L.x, L.y), xytext=(10, -4), textcoords="offset points", color="yellow", fontsize=12, weight="bold", path_effects=pe)
for k, (x, y) in mons.items():
    if k.startswith("mark"): ax.plot(x, y, marker="^", ms=11, mfc="none", mec="cyan", mew=2.2)
    else: ax.plot(x, y, marker="*", ms=16, color="magenta", mec="black")
ax.annotate("BASE STATION HERE: monument DR9925 (8508 Y)\nknown point, NAD 83(2011) epoch 2010.00\n44 07 18.13863 N  092 00 16.79318 W  h 188.669 m", mons["DR9925 base"], xytext=(-16, 14), textcoords="offset points", ha="right", color="magenta", fontsize=10, weight="bold", path_effects=pe)
# the driving line: lots in order along the road (straight segments, drawn thin)
P = [lots[l] for l, _ in STOPS]; ax.plot([p.x for p in P], [p.y for p in P], color="yellow", lw=1.2, ls=":")
ax.plot([X0 + 100, X0 + 600], [Y0 + 100, Y0 + 100], color="white", lw=4); ax.annotate("500 m", (X0 + 350, Y0 + 130), color="white", ha="center", fontsize=11, weight="bold")
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Field route: three lots on one road, 2.6 km end to end. Stop 1 lot 204 (RCG-type field, bluff), stop 2 lot 501 (base, mark, controls, cut sward, transect),\nstop 3 lot 2597 (best place 1, the 14 ha field). Fills = survey regions; orange line = transect; magenta = monuments; cyan = mark 1079. 2025 FSA imagery.", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_route.png", dpi=110); print("figures/field_route.png")
with open(f"{V}/field_route.gpx", "w") as g:
    g.write('<?xml version="1.0" encoding="UTF-8"?>\n<gpx version="1.1" creator="field_route_plan.py" xmlns="http://www.topografix.com/GPX/1/1">\n')
    for (lot, _), n in zip(STOPS, ("1", "2", "3")):
        lo, la = to_ll.transform(lots[lot].x, lots[lot].y); g.write(f'  <wpt lat="{la:.7f}" lon="{lo:.7f}"><name>stop{n}_lot{lot}</name></wpt>\n')
    g.write(f'  <wpt lat="{BASE_LAT:.8f}" lon="{BASE_LON:.8f}"><ele>{BASE_H_ELL}</ele><name>BASE_STATION_HERE_DR9925_8508Y</name><desc>known point NAD83(2011) epoch 2010.00; ellipsoid h 188.669 m; NAVD88 218.75 m (GEOID18)</desc></wpt>\n')
    for k, (x, y) in mons.items():
        if k.startswith("DR9925"): continue
        lo, la = to_ll.transform(x, y); g.write(f'  <wpt lat="{la:.7f}" lon="{lo:.7f}"><name>{escape(k)}</name></wpt>\n')
    g.write('  <rte><name>drive_order</name>\n' + "".join(f'    <rtept lat="{to_ll.transform(lots[l].x, lots[l].y)[1]:.7f}" lon="{to_ll.transform(lots[l].x, lots[l].y)[0]:.7f}"><name>lot {l}</name></rtept>\n' for l, _ in STOPS) + '  </rte>\n')
    for k, g_ in regions.items():
        for j, pg in enumerate(g_.geoms if hasattr(g_, "geoms") else [g_]):
            g.write(f'  <trk><name>{escape(k.split(":")[0])}_{j}</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in pg.exterior.coords) + '  </trkseg></trk>\n')
    g.write('  <trk><name>transect</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in tx) + '  </trkseg></trk>\n</gpx>\n')
print("data/derived/valley/field_route.gpx")
