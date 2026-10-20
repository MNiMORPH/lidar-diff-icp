#!/usr/bin/env python3
"""The best places to study on the field trip, valley-wide: zero-penetration floodplain sites
(F2: flat, public, zero gen1 penetration) that are GRASS or WETLAND today (WorldCover 2021),
with their full cluster outlines rebuilt from the tile layers, on an overview map and in a
GPX. Crop-covered F2 sites are listed but set aside: a zero-penetration cropped field in
November 2008 is standing crop, not sward (post-hoc view; WorldCover is the only cover source).
Reads data/derived/valley/sites_valley_ranked.csv and each tile's layers.
Writes figures/valley_best_places.png, data/derived/valley/best_places.gpx, best_places.csv.
"""
import json
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import shapely
from scipy.ndimage import label
from shapely.geometry import shape

V = "data/derived/valley"; to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
S = pd.read_csv(f"{V}/sites_valley_ranked.csv"); F2 = S[S.cls.str.startswith("F2")].copy()
F2["today"] = np.where(F2.worldcover.isin(["grass", "wetland"]), "grass/wetland", F2.worldcover)
F2 = F2.sort_values("ha", ascending=False)
print("F2 zero-penetration sites by today's cover:", F2.today.value_counts().to_dict())
keep = F2[F2.today == "grass/wetland"].head(10)
feats = lambda n: json.load(open(f"{V}/public_land/{n}"))["features"]
units = [shape(f["geometry"]) for f in feats("wma_units.geojson")] + [shape(f["geometry"]) for f in feats("state_parks.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]
lots = [shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"]
pub = shapely.union_all(units)
polys = {}
for _, r in keep.iterrows():
    D = f"data/derived/{r.tile}"; c = json.load(open(f"{D}/corrections.json")); X0, Y0 = c["bounds"][0], c["bounds"][1]; res = float(c.get("res") or c.get("res_m"))
    fp = np.load(f"{D}/floodplain_mask_clubb.npy").astype(bool); sl = np.load(f"{D}/slope.npy"); pf = np.load(f"{D}/penetration_failure.npy").astype(bool)
    ny, nx = sl.shape; jy, jx = np.mgrid[0:ny, 0:nx]; X = X0 + (jx + 0.5) * res; Y = Y0 + (jy + 0.5) * res
    m = fp & (sl <= 2.0) & np.isfinite(sl) & pf; lab, _ = label(m, structure=np.ones((3, 3))); k = lab[int((r.northing - Y0) / res), int((r.easting - X0) / res)]
    if k == 0:   # centroid off the cluster (concave shape): nearest cluster cell
        d = np.hypot(X - r.easting, Y - r.northing); d[lab == 0] = np.inf; k = lab[np.unravel_index(d.argmin(), d.shape)]
    mm = lab == k; pts = shapely.points(X[mm], Y[mm]); poly = shapely.union_all(shapely.buffer(pts, res * 0.5 + 0.6, cap_style="square")).simplify(1.0)
    polys[f"{r.tile}_{r.site}"] = (poly, r)
rows = []
for name, (poly, r) in polys.items():
    nl = min(np.hypot(g.x - r.easting, g.y - r.northing) for g in lots); rows.append(dict(name=name, ha=round(poly.area / 1e4, 1), lat=r.lat, lon=r.lon, dod_mm=r.raw_dod_mm, nelike_2008=r.nelike_2008, nelike_2021=r.nelike_2021, nelike_2026=r.nelike_2026, fingerprint=r.fingerprint_score, worldcover=r.worldcover, m_to_lot=round(nl), m_to_road=r.dist_wma_road_m))
B = pd.DataFrame(rows); B.to_csv(f"{V}/best_places.csv", index=False); pd.set_option("display.width", 220); print(B.to_string(index=False))
fig, ax = plt.subplots(figsize=(11, 12))
for g in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): ax.plot(*poly.exterior.xy, color="black", lw=0.7)
for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=0.8, alpha=0.7)
for p in lots: ax.plot(p.x, p.y, marker="P", color="yellow", ms=7, mec="black")
for i, (name, (poly, r)) in enumerate(polys.items()):
    for pg in (poly.geoms if hasattr(poly, "geoms") else [poly]): ax.fill(*pg.exterior.xy, color="#e34948", alpha=0.6); ax.plot(*pg.exterior.xy, color="#8b0000", lw=1)
    ax.annotate(f"{i + 1}", (r.easting, r.northing), xytext=(6, 6), textcoords="offset points", fontsize=11, weight="bold")
xs = keep.easting; ys = keep.northing; ax.set_xlim(xs.min() - 2000, xs.max() + 2000); ax.set_ylim(ys.min() - 2000, ys.max() + 2000); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Best places: zero-penetration floodplain sites (zero gen1 penetration) that are grass or wetland today, numbered by size; yellow = WMA parking", fontsize=9)
fig.tight_layout(); fig.savefig("figures/valley_best_places.png", dpi=110)
with open(f"{V}/best_places.gpx", "w") as g:
    g.write('<?xml version="1.0" encoding="UTF-8"?>\n<gpx version="1.1" creator="valley_best_places.py" xmlns="http://www.topografix.com/GPX/1/1">\n')
    for i, (name, (poly, r)) in enumerate(polys.items()):
        lo, la = to_ll.transform(r.easting, r.northing); g.write(f'  <wpt lat="{la:.7f}" lon="{lo:.7f}"><name>{i + 1}_{escape(name)}</name></wpt>\n')
        for j, pg in enumerate(poly.geoms if hasattr(poly, "geoms") else [poly]):
            g.write(f'  <trk><name>{i + 1}_{escape(name)}_{j}</name><trkseg>\n' + "".join(f'    <trkpt lat="{to_ll.transform(x, y)[1]:.7f}" lon="{to_ll.transform(x, y)[0]:.7f}"></trkpt>\n' for x, y in pg.exterior.coords) + '  </trkseg></trk>\n')
    g.write('</gpx>\n')
print("figures/valley_best_places.png, data/derived/valley/best_places.gpx")
