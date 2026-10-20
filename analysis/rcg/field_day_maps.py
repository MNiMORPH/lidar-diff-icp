#!/usr/bin/env python3
"""Site maps for the one field day (Andy, 2026-10-20: "Make site maps of where we should park
and go"): one map per parking lot on 2025 FSA imagery (MnGeo WMS), with the WMA road, the lot,
every site to visit with its label, the walking routes, and a scale bar; plus a GPX of all
waypoints for a handheld GPS.

Lot 204 (44.116878 N 92.000471 W): the RCG-type field R1-v140 (and R1-v149 beside the lot),
the closed-canopy lidar point F2-1 inside it, the bluff lowveg stand with its least-steep
entry and its <= 25 deg part, the gen2-high hillslope patches, and the hillslope control.
Lot 501 (44.120576 N 92.008147 W): the bare control C-v52, the meadow S-v324, mark
1079_2021_MN, the cut sward R2-v237, the closed-canopy point F2-3, and the 300 m class
transect with its 25 m waypoints.

Inputs: data/derived/valley/field_day_sites.csv, hillslope_transect_lot501.csv,
hillslope_lowveg_patches.csv; the Elba grids (slope, lowveg) for the bluff stand; the DNR
layers. Outputs: figures/field_day_lot204.png, figures/field_day_lot501.png,
data/derived/valley/field_day_waypoints.gpx.
"""
import heapq
import io
import json
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import shapely
from PIL import Image
from scipy.ndimage import label, uniform_filter
from shapely.geometry import shape

V = "data/derived/valley"; PL = f"{V}/public_land"; S = "data/derived/elbaext_delong_star_nosurf"; X0, Y0, RES = 575600.0, 4882200.0, 5.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
feats = lambda n: json.load(open(f"{PL}/{n}"))["features"]
roads = [shape(f["geometry"]) for f in feats("wma_roads.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]
lots = {f["properties"]["objectid"]: shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"}
units = [shape(f["geometry"]) for f in feats("wma_units.geojson")]
F = pd.read_csv(f"{V}/field_day_sites.csv"); T = pd.read_csv(f"{V}/hillslope_transect_lot501.csv"); P = pd.read_csv(f"{V}/hillslope_lowveg_patches.csv")
sl = np.load(f"{S}/slope.npy"); lv = np.load(f"{S}/lowveg_grid.npy") if False else np.load("data/derived/elbaext_independent/lowveg_grid.npy"); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); z = np.load(f"{S}/z_after.npy")
ny, nx = sl.shape

# the bluff stand (as hillslope_lowveg_stands found it) and its least-steep entry from lot 204
lvs = uniform_filter(np.nan_to_num(lv), size=5); hill = np.isfinite(sl) & ~fp & (sl >= 10) & (sl <= 35) & np.isfinite(lv)
lab, _ = label(hill & (lvs >= 0.4), structure=np.ones((3, 3))); stand = lab == lab[int((4885218 - Y0) / RES), int((579768 - X0) / RES)]
lx, ly = int((lots[204].x - X0) / RES), int((lots[204].y - Y0) / RES)
cost = np.where(np.isfinite(sl), np.exp(sl / 8.0), 1e6); dist = np.full((ny, nx), np.inf); prev = {}; dist[ly, lx] = 0; pq = [(0, ly, lx)]; target = None
while pq:
    d, y, x = heapq.heappop(pq)
    if d > dist[y, x]: continue
    if stand[y, x]: target = (y, x); break
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == dx == 0: continue
            yy, xx = y + dy, x + dx
            if 0 <= yy < ny and 0 <= xx < nx:
                nd = d + cost[yy, xx] * np.hypot(dy, dx)
                if nd < dist[yy, xx]: dist[yy, xx] = nd; prev[(yy, xx)] = (y, x); heapq.heappush(pq, (nd, yy, xx))
path = [target]
while path[-1] != (ly, lx): path.append(prev[path[-1]])
path = path[::-1]; route = np.array([(X0 + (x + 0.5) * RES, Y0 + (y + 0.5) * RES) for y, x in path])

def wms(X0_, Y0_, X1_, Y1_, px=1400):
    url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0_},{Y0_},{X1_},{Y1_}&WIDTH={px}&HEIGHT={int(px * (Y1_ - Y0_) / (X1_ - X0_))}&FORMAT=image/jpeg"
    try: return np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
    except Exception as e: print("imagery not fetched:", e); return None

def base_map(ax, bb):
    img = wms(*bb)
    if img is not None: ax.imshow(img, extent=[bb[0], bb[2], bb[1], bb[3]])
    for g in units:
        for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=1)
    for rd in roads: ax.plot(*rd.xy, color="white", lw=2, ls="--")
    for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=1.5)
    for lid, g in lots.items():
        if bb[0] < g.x < bb[2] and bb[1] < g.y < bb[3]: ax.plot(g.x, g.y, marker="P", color="yellow", ms=16, mec="black"); ax.annotate(f"lot {lid}", (g.x, g.y), xytext=(8, 8), textcoords="offset points", color="yellow", fontsize=12, weight="bold")
    ax.plot([bb[0] + 60, bb[0] + 260], [bb[1] + 60, bb[1] + 60], color="white", lw=4); ax.annotate("200 m", (bb[0] + 160, bb[1] + 80), color="white", ha="center", fontsize=11, weight="bold")
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3]); ax.set_xticks([]); ax.set_yticks([])

def mark(ax, x, y, text, color="white", m="o"):
    ax.plot(x, y, marker=m, ms=12, mfc="none", mec=color, mew=2.5); ax.annotate(text, (x, y), xytext=(9, 6), textcoords="offset points", color=color, fontsize=10, weight="bold", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])

wps = []
# ---- map 1: lot 204
bb = (579400, 4884900, 580300, 4885900); fig, ax = plt.subplots(figsize=(11, 12)); base_map(ax, bb)
jy, jx = np.nonzero(stand); ax.scatter(X0 + (jx + 0.5) * RES, Y0 + (jy + 0.5) * RES, s=9, c=np.where(sl[stand] <= 25, "#1baf7a", "#e34948"), alpha=0.7, marker="s", lw=0)
ax.plot(route[:, 0], route[:, 1], color="#eda100", lw=3); mark(ax, route[-1, 0], route[-1, 1], "stand entry (least-steep route, max 21 deg)", "#eda100", "^")
for _, r in F[(F.lot == 204) & (F.type.isin(["R1", "F2"]))].iterrows(): mark(ax, r.easting, r.northing, f"{r.site} {r.type}", "white"); wps.append((r.site, r.easting, r.northing))
for _, r in P[(P.lot == "204") | (P.lot == 204)].iterrows():
    if r.cls.startswith("high lowveg (>= p75) AND"): mark(ax, r.E, r.N, f"gen2-high +{r.dod_med:.0f} mm, lowveg {r.lowveg:.2f}", "#e34948", "D"); wps.append((f"gen2high_{r.E:.0f}_{r.N:.0f}", r.E, r.N))
ctrl = P[(P.cls.str.startswith("low lowveg")) & ((P.lot == "204") | (P.lot == 204))].sort_values("lot_m").iloc[0]; mark(ax, ctrl.E, ctrl.N, f"hillslope control, lowveg {ctrl.lowveg:.2f}", "#1baf7a", "s"); wps.append(("hillslope_control", ctrl.E, ctrl.N))
wps += [("lot_204", lots[204].x, lots[204].y), ("stand_entry", route[-1, 0], route[-1, 1])]
ax.set_title("Stop 1, lot 204: RCG-type field R1-v140 (+ R1-v149 by the lot), closed-canopy point F2-1, the bluff lowveg stand (green <= 25 deg, red steeper),\nits least-steep entry (orange), gen2-high patches (red diamonds), hillslope control (green square). White dashed = WMA road; black = WMA boundary. 2025 FSA imagery.", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_day_lot204.png", dpi=110); print("figures/field_day_lot204.png")

# ---- map 2: lot 501
bb = (578900, 4885300, 579900, 4886300); fig, ax = plt.subplots(figsize=(11, 12)); base_map(ax, bb)
ax.plot(T.easting, T.northing, color="#eda100", lw=3); wp = T[np.isclose(T.dist_m % 25, 0)]
for _, r in wp.iterrows(): ax.plot(r.easting, r.northing, marker="o", ms=5, color="#eda100", mec="black"); wps.append((f"transect_{int(r.dist_m)}m", r.easting, r.northing))
ax.annotate("class transect, 300 m, bearing 300", (T.easting.iloc[0], T.northing.iloc[0]), xytext=(10, -14), textcoords="offset points", color="#eda100", fontsize=10, weight="bold")
for _, r in F[(F.lot == 501) & (F.type.isin(["C", "S", "R2", "F2"]))].sort_values("dist_m").groupby("type").head(1).iterrows(): mark(ax, r.easting, r.northing, f"{r.site} {r.type}", "white"); wps.append((r.site, r.easting, r.northing))
mark(ax, 579659.824, 4885857.071, "mark 1079_2021_MN", "cyan", "^"); wps += [("mark_1079_2021_MN", 579659.824, 4885857.071), ("lot_501", lots[501].x, lots[501].y)]
ax.set_title("Stop 2, lot 501: bare control C-v52, meadow S-v324, mark 1079, sward-both-epochs R2-v237 (cut since), closed-canopy F2-3,\nand the class transect (orange, dots every 25 m). White dashed = WMA road; black = WMA boundary. 2025 FSA imagery.", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_day_lot501.png", dpi=110); print("figures/field_day_lot501.png")

# ---- GPX
with open(f"{V}/field_day_waypoints.gpx", "w") as g:
    g.write('<?xml version="1.0" encoding="UTF-8"?>\n<gpx version="1.1" creator="lidar-diff-icp field_day_maps.py" xmlns="http://www.topografix.com/GPX/1/1">\n')
    for name, x, y in wps:
        lo, la = to_ll.transform(x, y); g.write(f'  <wpt lat="{la:.7f}" lon="{lo:.7f}"><name>{name}</name></wpt>\n')
    g.write('  <trk><name>class_transect</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in zip(T.easting, T.northing)) + '  </trkseg></trk>\n')
    g.write('  <trk><name>route_lot204_to_stand</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in route) + '  </trkseg></trk>\n</gpx>\n')
print(f"GPX: {V}/field_day_waypoints.gpx ({len(wps)} waypoints, 2 tracks)")
