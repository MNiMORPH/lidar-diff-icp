#!/usr/bin/env python3
"""Site maps for the one field day (Andy, 2026-10-20: "Make site maps of where we should park
and go"): one map per parking lot on 2025 FSA imagery (MnGeo WMS), with the WMA road, the lot,
every site to visit with its label, the walking routes, and a scale bar; plus a GPX of all
waypoints for a handheld GPS.

Lot 204 (44.116878 N 92.000471 W): the RCG-type field R1-v140 (and R1-v149 beside the lot),
the zero-penetration lidar point F2-1 inside it, the bluff lowveg stand with its least-steep
entry and its <= 25 deg part, the gen2-high hillslope patches, and the hillslope control.
Lot 501 (44.120576 N 92.008147 W): the bare control C-v52, the meadow S-v324, mark
1079_2021_MN, the cut sward R2-v237, the zero-penetration point F2-3, and the 300 m class
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

# monuments (NGS datasheets, read 2026-10-20): DR9925 = MnDOT "8508 Y" rod, NAD 83(2011) adjusted,
# the BASE STATION; DG8385 = USGS "9 DRL" leveled disk, the vertical CHECK shot
MONUMENTS = {"DR9925 (8508 Y): BASE": (-92.004665, 44.121705, (-14, 14), "right"), "DG8385 (9 DRL): vertical check": (-92.003536, 44.120247, (10, -16), "left")}
from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
def draw_monuments(ax, bb):
    for name, (lo, la, off, ha) in MONUMENTS.items():
        x, y = from_ll.transform(lo, la)
        if bb[0] < x < bb[2] and bb[1] < y < bb[3]:
            ax.plot(x, y, marker="*", ms=20, color="magenta", mec="black", mew=1.2); ax.annotate(name, (x, y), xytext=off, textcoords="offset points", ha=ha, color="magenta", fontsize=11, weight="bold", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])
            wps.append((name.split(" ")[0] + "_" + name.split("(")[1].split(")")[0].replace(" ", ""), x, y))

def mark(ax, x, y, text, color="white", m="o"):
    ax.plot(x, y, marker=m, ms=12, mfc="none", mec=color, mew=2.5); ax.annotate(text, (x, y), xytext=(9, 6), textcoords="offset points", color=color, fontsize=10, weight="bold", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])

wps = []
# ---- REGIONS (Andy, 2026-10-20: "Regions are the only ones worth visiting; hitting an exact
# point might be beyond the students' abilities"): each site is its full cell cluster, drawn
# as a filled outline, and its outline goes to the GPX as a track.
import rasterio
from scipy.ndimage import binary_dilation
def cells_to_poly(mask, x0, y0, res, dy=None):
    """cells -> polygon. Elba grids: row index increases NORTHWARD (y = y0 + (jy+0.5)*res).
    North-up GeoTIFFs: pass dy = transform.e (negative) and y0 = transform.f (the top edge)."""
    jy, jx = np.nonzero(mask); dy = res if dy is None else dy; pts = shapely.points(x0 + (jx + 0.5) * res, y0 + (jy + 0.5) * dy)
    return shapely.union_all(shapely.buffer(pts, res * 0.5 + 0.6, cap_style="square")).simplify(1.0)
with rasterio.open(f"{V}/valley_classes.tif") as r: cls10 = r.read(1); tr10 = r.transform
def valley_region(site):
    row = F[F.site == site].iloc[0]; v = {"R1": 1, "R2": 2, "S": 3, "C": 4, "E": 5}[row.type]
    lab_, _ = label(cls10 == v, structure=np.ones((3, 3))); c, rr = ~tr10 * (row.easting, row.northing); k = lab_[int(rr), int(c)]
    return cells_to_poly(lab_ == k, tr10.c, tr10.f, tr10.a, dy=tr10.e) if k else None
def elba_region(mask, x, y):
    lab_, _ = label(mask, structure=np.ones((3, 3))); k = lab_[int((y - Y0) / RES), int((x - X0) / RES)]
    return cells_to_poly(lab_ == k, X0, Y0, RES) if k else None
pf = np.load(f"{S}/penetration_failure.npy").astype(bool); flat = fp & (sl <= 2.0) & np.isfinite(sl)
units_u = shapely.union_all(units); jy_, jx_ = np.mgrid[0:ny, 0:nx]; public = shapely.contains_xy(units_u, (X0 + (jx_ + 0.5) * RES).ravel(), (Y0 + (jy_ + 0.5) * RES).ravel()).reshape(ny, nx)
REG = {}
for site in ("R1-v140", "R1-v149", "R2-v237", "C-v52", "S-v324"): REG[site] = valley_region(site)
REG["F2-1"] = elba_region(flat & public & pf, 579959.4, 4885620.6); REG["F2-3"] = elba_region(flat & public & pf, 579812.5, 4886050.8)
REG["bluff stand, <= 25 deg"] = cells_to_poly(stand & (sl <= 25), X0, Y0, RES); REG["bluff stand, steeper"] = cells_to_poly(stand & (sl > 25), X0, Y0, RES)
base_ = np.load(f"{S}/dod.npy") * 1000; base_ = base_ - np.nanmedian(base_[np.load(f"{S}/stable.npy").astype(bool) & np.isfinite(base_)]); lod_ = np.load(f"{S}/lod.npy") * 1000
hiq = np.nanpercentile(lv[hill], 75); loq = np.nanpercentile(lv[hill], 25)
def big_clusters(mask, min_cells):
    lab_, _ = label(mask, structure=np.ones((3, 3))); sz = np.bincount(lab_.ravel()); keep = np.isin(lab_, np.flatnonzero(sz >= min_cells)[1:] if sz[0] >= 0 else []); return keep & (lab_ > 0)
hillp = hill & public & np.isfinite(base_)
near204 = np.hypot(X0 + (jx_ + 0.5) * RES - lots[204].x, Y0 + (jy_ + 0.5) * RES - lots[204].y) < 400
REG["gen2-high patches"] = cells_to_poly(big_clusters(hillp & (lv >= hiq) & (base_ > lod_), 16) & near204, X0, Y0, RES)   # the >= 16-cell clusters of hillslope_lowveg_patches.csv
ctrl_mask = hillp & (lv < loq) & (np.abs(base_) < lod_ / 2); lab_c, _ = label(ctrl_mask, structure=np.ones((3, 3))); sz_c = np.bincount(lab_c.ravel())
kc = lab_c[int((4885442 - Y0) / RES), int((579728 - X0) / RES)]; print(f"  hillslope control cluster at the recorded cell: label {kc}, {sz_c[kc] if kc else 0} cells; largest control cluster within 400 m of lot 204: {max(sz_c[k] for k in np.unique(lab_c[near204 & ctrl_mask]) if k) if (near204 & ctrl_mask).any() else 0} cells")
kbest = max((k for k in np.unique(lab_c[near204 & ctrl_mask]) if k), key=lambda k: sz_c[k]); REG["hillslope control"] = cells_to_poly(lab_c == kbest, X0, Y0, RES)
for k, g in list(REG.items()):
    if g is not None and hasattr(g, "geoms"): REG[k] = shapely.union_all([pg for pg in g.geoms if pg.area >= 500.0]) if any(pg.area >= 500.0 for pg in g.geoms) else None   # fragments < 0.05 ha dropped (MINE)
for k, g in REG.items(): print(f"  region {k}: {g.area / 1e4:.2f} ha" if g is not None else f"  region {k}: NONE")
def draw_region(ax, g, color, text, alpha=0.3):
    if g is None: return
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.fill(*poly.exterior.xy, color=color, alpha=alpha, lw=0); ax.plot(*poly.exterior.xy, color=color, lw=2)
    c = g.representative_point(); ax.annotate(text, (c.x, c.y), color=color, fontsize=10, weight="bold", ha="center", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])
    for i, poly in enumerate(g.geoms if hasattr(g, "geoms") else [g]): TRACKS.append((f"{text.split(':')[0]}{'' if i == 0 else '_' + str(i)}", np.asarray(poly.exterior.coords)))
TRACKS = []
# ---- map 1: lot 204
bb = (579400, 4884900, 580300, 4885900); fig, ax = plt.subplots(figsize=(11, 12)); base_map(ax, bb); draw_monuments(ax, bb)
draw_region(ax, REG["R1-v140"], "white", "R1-v140: RCG-type field, both epochs"); draw_region(ax, REG["R1-v149"], "white", "R1-v149"); draw_region(ax, REG["F2-1"], "#eb6834", "F2-1: zero penetration")
draw_region(ax, REG["bluff stand, <= 25 deg"], "#1baf7a", "bluff stand (green <= 25 deg)"); draw_region(ax, REG["bluff stand, steeper"], "#e34948", "", 0.25)
draw_region(ax, REG["gen2-high patches"], "#e87ba4", "gen2-high patches"); draw_region(ax, REG["hillslope control"], "#eda100", "hillslope control")
ax.plot(route[:, 0], route[:, 1], color="#eda100", lw=3); ax.plot(route[-1, 0], route[-1, 1], marker="^", ms=12, mfc="none", mec="#eda100", mew=2.5)
wps += [("lot_204", lots[204].x, lots[204].y), ("stand_entry", route[-1, 0], route[-1, 1])]
ax.set_title("Stop 1, lot 204. Survey anywhere inside a region: RCG-type field R1-v140 and R1-v149 (white), zero-penetration F2-1 (orange),\nbluff lowveg stand (green <= 25 deg, red steeper; under tree canopy), gen2-high patches (pink), hillslope control (yellow).\nOrange line = least-steep route to the stand. Magenta stars = monuments (DR9925 base, DG8385 check). White dashed = WMA road. 2025 FSA imagery.", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_day_lot204.png", dpi=110); print("figures/field_day_lot204.png")

# ---- map 2: lot 501
bb = (578900, 4885300, 579900, 4886300); fig, ax = plt.subplots(figsize=(11, 12)); base_map(ax, bb); draw_monuments(ax, bb)
ax.plot(T.easting, T.northing, color="#eda100", lw=3); wp = T[np.isclose(T.dist_m % 25, 0)]
for _, r in wp.iterrows(): ax.plot(r.easting, r.northing, marker="o", ms=5, color="#eda100", mec="black"); wps.append((f"transect_{int(r.dist_m)}m", r.easting, r.northing))
ax.annotate("class transect, 300 m, bearing 300", (T.easting.iloc[0], T.northing.iloc[0]), xytext=(10, -14), textcoords="offset points", color="#eda100", fontsize=10, weight="bold")
draw_region(ax, REG["C-v52"], "#eda100", "C-v52: bare control"); draw_region(ax, REG["S-v324"], "#4a3aa7", "S-v324: meadow"); draw_region(ax, REG["R2-v237"], "#eb6834", "R2-v237: sward 2008+2021, cut since"); draw_region(ax, REG["F2-3"], "white", "F2-3: zero penetration")
ax.plot(579659.824, 4885857.071, marker="^", ms=12, mfc="none", mec="cyan", mew=2.5); ax.annotate("benchmark 1079_2021_MN", (579659.824, 4885857.071), xytext=(12, -22), textcoords="offset points", color="cyan", fontsize=10, weight="bold", path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")]); wps += [("mark_1079_2021_MN", 579659.824, 4885857.071), ("lot_501", lots[501].x, lots[501].y)]
ax.set_title("Stop 2, lot 501. Survey anywhere inside a region: bare control C-v52 (yellow), meadow S-v324 (purple),\nsward-both-epochs R2-v237 (orange), zero-penetration F2-3 (white). Exact points: benchmark 1079 (cyan), monuments DR9925 = 8508 Y BASE\nand DG8385 = 9 DRL check (magenta stars). Class transect: orange line, dots every 25 m. White dashed = WMA road. 2025 FSA imagery.", fontsize=9)
fig.tight_layout(); fig.savefig("figures/field_day_lot501.png", dpi=110); print("figures/field_day_lot501.png")

# ---- GPX
from xml.sax.saxutils import escape
with open(f"{V}/field_day_waypoints.gpx", "w") as g:
    g.write('<?xml version="1.0" encoding="UTF-8"?>\n<gpx version="1.1" creator="lidar-diff-icp field_day_maps.py" xmlns="http://www.topografix.com/GPX/1/1">\n')
    for name, x, y in wps:
        lo, la = to_ll.transform(x, y); g.write(f'  <wpt lat="{la:.7f}" lon="{lo:.7f}"><name>{escape(name)}</name></wpt>\n')
    g.write('  <trk><name>class_transect</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in zip(T.easting, T.northing)) + '  </trkseg></trk>\n')
    from xml.sax.saxutils import escape
    for name, xy in TRACKS: g.write(f'  <trk><name>region_{escape(name)}</name><trkseg>\n' + ''.join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in xy) + '  </trkseg></trk>\n')
    g.write('  <trk><name>route_lot204_to_stand</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in route) + '  </trkseg></trk>\n</gpx>\n')
print(f"GPX: {V}/field_day_waypoints.gpx ({len(wps)} waypoints, {len(TRACKS) + 2} tracks: region outlines, the transect, the bluff route)")
