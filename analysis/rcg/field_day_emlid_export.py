#!/usr/bin/env python3
"""Field-day navigation files for Emlid Flow (free plan) and a phone map app.

Emlid Flow's free plan imports and stakes out POINTS only (CSV, DXF, Shapefile, KML; lines,
polygons and satellite background need the paid plan -- Emlid docs, 2026-10-20). So each
survey REGION becomes a ring of named points along its outline every ~20 m plus one CENTRE
point, which Flow draws on its map and can stake out; the phone app gets the polygons.

Reads data/derived/valley/field_day_waypoints.gpx (from field_day_maps.py). Writes, in
data/derived/valley/:
  field_day_points_emlid.csv   Name, Longitude, Latitude (WGS84) -- Flow's CSV import
  field_day_points_emlid.kml   the same points as KML (Flow imports KML points too)
  field_day_points_emlid.shp   the same as a point Shapefile (+ .shx .dbf .prj)
  field_day_regions_phone.kml  the region POLYGONS, centres and waypoints for a phone app
"""
import re
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
import pyproj
import shapefile
import shapely

V = "data/derived/valley"
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True); from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
root = ET.parse(f"{V}/field_day_waypoints.gpx").getroot(); ns = {"g": "http://www.topografix.com/GPX/1/1"}
rows = []; polys = {}
for w in root.findall("g:wpt", ns): rows.append((w.find("g:name", ns).text, float(w.get("lat")), float(w.get("lon")), "waypoint"))
for t in root.findall("g:trk", ns):
    name = t.find("g:name", ns).text; pts = [(float(p.get("lon")), float(p.get("lat"))) for p in t.findall(".//g:trkpt", ns)]
    if name.startswith("region_"):
        xy = np.array([from_ll.transform(lo, la) for lo, la in pts]); ring = shapely.LinearRing(xy); n = max(4, int(ring.length // 20))
        short = re.sub(r"[^A-Za-z0-9]+", "", name.replace("region_", ""))[:14]
        for i in range(n):
            p = ring.interpolate(i / n, normalized=True); lo, la = to_ll.transform(p.x, p.y); rows.append((f"{short}-{i+1:02d}", la, lo, "region outline"))
        c = shapely.Polygon(xy).representative_point(); lo, la = to_ll.transform(c.x, c.y); rows.append((f"{short}-CENTRE", la, lo, "region centre")); polys[name.replace("region_", "")] = pts
    elif name == "class_transect":
        for i, (lo, la) in enumerate(pts[::5]): rows.append((f"transect-{i*25:03d}m", la, lo, "transect"))
P = pd.DataFrame(rows, columns=["Name", "Latitude", "Longitude", "Note"]).drop_duplicates("Name")
P[["Name", "Longitude", "Latitude"]].to_csv(f"{V}/field_day_points_emlid.csv", index=False)
def kml(df, path, extra=""):
    with open(path, "w") as k:
        k.write('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Whitewater field day</name>\n')
        for _, r in df.iterrows(): k.write(f'<Placemark><name>{escape(r.Name)}</name><description>{escape(r.Note)}</description><Point><coordinates>{r.Longitude:.7f},{r.Latitude:.7f},0</coordinates></Point></Placemark>\n')
        k.write(extra + '</Document></kml>\n')
kml(P, f"{V}/field_day_points_emlid.kml")
poly_kml = "".join(f'<Placemark><name>{escape(n)}</name><Style><LineStyle><color>ff00ffff</color><width>3</width></LineStyle><PolyStyle><color>4000ffff</color></PolyStyle></Style><Polygon><outerBoundaryIs><LinearRing><coordinates>' + " ".join(f"{lo:.7f},{la:.7f},0" for lo, la in pts) + '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n' for n, pts in polys.items())
kml(P[P.Note != "region outline"], f"{V}/field_day_regions_phone.kml", poly_kml)
with shapefile.Writer(f"{V}/field_day_points_emlid", shapeType=shapefile.POINT) as w:
    w.field("Name", "C", 40); w.field("Note", "C", 20)
    for _, r in P.iterrows(): w.point(r.Longitude, r.Latitude); w.record(r.Name, r.Note)
open(f"{V}/field_day_points_emlid.prj", "w").write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]]')
print(f"points: {len(P)} ({P.Note.value_counts().to_dict()}); regions as polygons: {len(polys)}")
print("written: field_day_points_emlid.csv / .kml / .shp(.shx .dbf .prj); field_day_regions_phone.kml")
