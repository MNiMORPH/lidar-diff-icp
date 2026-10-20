#!/usr/bin/env python3
"""Zoomed-in looks at the field-survey sites with every valley tile and the imagery in hand
(Andy, 2026-10-20: "create a zoomed-in look at our field-survey sites" and "produce a lidar
difference map of the area(s), highlighting likely spurious 'deposition' and 'erosion'").

For each window (stop 1 lot 204, stop 2 lot 501, the valley-wide best place 1 by the 2021
photo site, and the whole reach) three panels:
  A  2025 FSA imagery (MnGeo WMS) with the survey REGION outlines, lots, monuments, transect;
  B  the 2021 - 2008 difference, mosaicked from the Phase B tiles (5 m, star frame, no
     correction surface, stable-ground median removed per tile), +/- 300 mm, grey = no data;
  C  the same difference with the LIKELY SPURIOUS clusters hatched and outlined (a 3x3
     opening drops isolated cells; the json carries both the raw and the drawn counts):
       erosion, dashed //   = zero 2008 ground return (penetration_failure) on the flat
                           floodplain: the 2008 "ground" is a canopy top, so the lowering
                           is not ground (the F2 mechanism);
       erosion, dotted ..  = flat-floodplain cells that the vegetation LoD silences: the 2008
                           return span (p90-p10) makes the detection insignificant
                           (change.npy detected, change_veg.npy not);
       deposition, solid \\ = hillslope (10-35 deg) cells above +LoD whose gen2 near-ground
                           vegetation share (lowveg) is in the tile's top quartile or that
                           WorldCover 2021 calls tree: gen2 (May 2021, in leaf) reads the
                           vegetation top as ground (the gen2-vegetation-offset mechanism).
     Everything else that exceeds the LoD is left as it is: a candidate for real change.
All three rules are post-hoc views of layers already in the tile directories; the
quartile is the tile's own; nothing is fitted here. Outlines only, no text on the maps.

Outputs figures/zoom_<name>_{imagery,dod,spurious}.png and a combined
figures/zoom_<name>.png per window; data/derived/valley/spurious_<name>.json with the
cell counts behind each outline.
"""
import io
import json
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import rasterio
import shapely
from matplotlib.colors import TwoSlopeNorm
from PIL import Image
from scipy.ndimage import binary_dilation
from shapely.geometry import shape

V = "data/derived/valley"; RES = 5.0
from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
TILES = {t: v["bounds"] for t, v in json.load(open(f"{V}/public_tiles.json")).items()}
WINDOWS = {"stop1_lot204": (579250, 4884900, 580700, 4885950), "stop2_lot501": (578900, 4885250, 580100, 4886450),
           "bestplace1_photo": (579000, 4888000, 580900, 4889900), "reach": (578000, 4883500, 581500, 4890500)}

def tile_layer(t, name):
    p = f"data/derived/{t}/{name}"
    return np.load(p) if __import__("os").path.exists(p) else None

def mosaic(bb, name, dtype=float, fill=np.nan):
    X0, Y0, X1, Y1 = bb; nx, ny = int((X1 - X0) / RES), int((Y1 - Y0) / RES); out = np.full((ny, nx), fill, dtype)
    for t, b in TILES.items():
        c = f"data/derived/{t}/corrections.json"
        if not __import__("os").path.exists(c) or tile_layer(t, "dod.npy") is None: continue
        tb = json.load(open(c))["bounds"]
        if tb[2] <= X0 or tb[0] >= X1 or tb[3] <= Y0 or tb[1] >= Y1: continue
        a = tile_layer(t, name)
        if a is None: continue
        if name == "dod.npy":
            st = tile_layer(t, "stable.npy").astype(bool); a = a * 1000.0; a = a - np.nanmedian(a[st & np.isfinite(a)])
        elif name in ("lod.npy", "lod_veg.npy"): a = a * 1000.0
        ix0 = int(max(tb[0], X0) - X0) // int(RES); iy0 = int(max(tb[1], Y0) - Y0) // int(RES); jx0 = int(max(tb[0], X0) - tb[0]) // int(RES); jy0 = int(max(tb[1], Y0) - tb[1]) // int(RES)
        w = int(min(tb[2], X1) - max(tb[0], X0)) // int(RES); h = int(min(tb[3], Y1) - max(tb[1], Y0)) // int(RES)
        sub = a[jy0:jy0 + h, jx0:jx0 + w]; tgt = out[iy0:iy0 + h, ix0:ix0 + w]
        m = np.isfinite(sub) if sub.dtype.kind == "f" else np.ones(sub.shape, bool)
        if tgt.dtype.kind == "f": fillm = ~np.isfinite(tgt) & m
        else: fillm = m
        tgt[fillm] = sub[fillm]
    return out

def lowveg_q75(bb):
    """the top-quartile edge of lowveg on hillslope cells, per tile, mosaicked as a threshold grid"""
    X0, Y0, X1, Y1 = bb; nx, ny = int((X1 - X0) / RES), int((Y1 - Y0) / RES); out = np.full((ny, nx), np.nan)
    for t, b in TILES.items():
        c = f"data/derived/{t}/corrections.json"
        if not __import__("os").path.exists(c) or tile_layer(t, "lowveg_grid.npy") is None: continue
        tb = json.load(open(c))["bounds"]
        if tb[2] <= X0 or tb[0] >= X1 or tb[3] <= Y0 or tb[1] >= Y1: continue
        lv = tile_layer(t, "lowveg_grid.npy"); sl = tile_layer(t, "slope.npy"); fp = tile_layer(t, "floodplain_mask_clubb.npy").astype(bool)
        q = np.nanpercentile(lv[(sl >= 10) & (sl <= 35) & ~fp & np.isfinite(lv)], 75)
        ix0 = int(max(tb[0], X0) - X0) // int(RES); iy0 = int(max(tb[1], Y0) - Y0) // int(RES); w = int(min(tb[2], X1) - max(tb[0], X0)) // int(RES); h = int(min(tb[3], Y1) - max(tb[1], Y0)) // int(RES)
        tgt = out[iy0:iy0 + h, ix0:ix0 + w]; tgt[~np.isfinite(tgt)] = q
    return out

# overlays: survey regions, transect, lots, monuments (field_day_waypoints.gpx), best places (best_places.gpx)
def gpx_tracks(path):
    root = ET.parse(path).getroot(); ns = {"g": "http://www.topografix.com/GPX/1/1"}
    tr = [(t.find("g:name", ns).text, np.array([from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in t.findall(".//g:trkpt", ns)])) for t in root.findall("g:trk", ns)]
    wp = {w.find("g:name", ns).text: from_ll.transform(float(w.get("lon")), float(w.get("lat"))) for w in root.findall("g:wpt", ns)}
    return tr, wp
tracks, wpts = gpx_tracks(f"{V}/field_day_waypoints.gpx"); btracks, bwpts = gpx_tracks(f"{V}/best_places.gpx")
feats = lambda n: json.load(open(f"{V}/public_land/{n}"))["features"]
roads = [shape(f["geometry"]) for f in feats("wma_roads.geojson")]; rivers = [shape(f["geometry"]) for f in feats("rivers_whitewater.geojson")]
lots = [shape(f["geometry"]) for f in feats("wma_parking.geojson") if f["properties"].get("facility_subgroup") == "Parking Lot"]
with rasterio.open(f"{V}/worldcover.tif") as r: WC = r.read(1); WCT = r.transform

def wms(bb, px=1400):
    X0, Y0, X1, Y1 = bb; url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH={px}&HEIGHT={int(px * (Y1 - Y0) / (X1 - X0))}&FORMAT=image/jpeg"
    try: return np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
    except Exception as e: print("imagery not fetched:", e); return None

def overlays(ax, bb, lw=1.3):
    for name, xy in tracks:
        if name.startswith("region_"): ax.plot(xy[:, 0], xy[:, 1], color="black", lw=lw)
        elif name == "class_transect": ax.plot(xy[:, 0], xy[:, 1], color="black", lw=lw, ls="--")
    for name, xy in btracks: ax.plot(xy[:, 0], xy[:, 1], color="black", lw=lw)
    for rd in roads: ax.plot(*rd.xy, color="white", lw=1.0, ls=":")
    for rv in rivers: ax.plot(*rv.xy, color="#2a78d6", lw=1.0)
    for g in lots:
        if bb[0] < g.x < bb[2] and bb[1] < g.y < bb[3]: ax.plot(g.x, g.y, marker="P", ms=11, color="yellow", mec="black")
    for k, (x, y) in wpts.items():
        if (k.startswith("DR9925") or k.startswith("DG8385")) and bb[0] < x < bb[2] and bb[1] < y < bb[3]: ax.plot(x, y, marker="*", ms=14, color="magenta", mec="black")
    ax.plot([bb[0] + 0.04 * (bb[2] - bb[0]), bb[0] + 0.04 * (bb[2] - bb[0]) + 200], [bb[1] + 0.04 * (bb[3] - bb[1])] * 2, color="black", lw=4)
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3]); ax.set_xticks([]); ax.set_yticks([])

from scipy.ndimage import binary_opening
def outline(ax, mask, bb, color, ls, hatch):
    """clusters only (a 3x3 opening drops isolated cells, which drew as a speckle of dots), as a
    hatched wash with its outline"""
    X0, Y0, X1, Y1 = bb; m = binary_opening(mask, structure=np.ones((3, 3)))
    if m.any():
        xs = np.linspace(X0 + RES / 2, X1 - RES / 2, m.shape[1]); ys = np.linspace(Y0 + RES / 2, Y1 - RES / 2, m.shape[0])
        ax.contourf(xs, ys, m.astype(float), levels=[0.5, 1.5], colors="none", hatches=[hatch]); ax.contour(xs, ys, m.astype(float), levels=[0.5], colors=color, linewidths=0.9, linestyles=ls)
    return int(m.sum())

cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#bdbdbd"); norm = TwoSlopeNorm(vmin=-300, vcenter=0, vmax=300)
for name, bb in WINDOWS.items():
    X0, Y0, X1, Y1 = bb; ext = [X0, X1, Y0, Y1]
    dod = mosaic(bb, "dod.npy"); lod = mosaic(bb, "lod.npy"); lodv = mosaic(bb, "lod_veg.npy"); sl = mosaic(bb, "slope.npy")
    pf = mosaic(bb, "penetration_failure.npy", bool, False); fp = mosaic(bb, "floodplain_mask_clubb.npy", bool, False); ch = mosaic(bb, "change.npy", bool, False); chv = mosaic(bb, "change_veg.npy", bool, False); lv = mosaic(bb, "lowveg_grid.npy"); q75 = lowveg_q75(bb)
    ny, nx = dod.shape; jy, jx = np.mgrid[0:ny, 0:nx]; X = X0 + (jx + 0.5) * RES; Y = Y0 + (jy + 0.5) * RES
    c, rr = ~WCT * (X.ravel(), Y.ravel()); tree = (WC[np.clip(rr.astype(int), 0, WC.shape[0] - 1), np.clip(c.astype(int), 0, WC.shape[1] - 1)] == 10).reshape(ny, nx)
    flat = fp & (sl <= 2.0); hill = (sl >= 10) & (sl <= 35) & ~fp
    sp_ero_pf = flat & pf & (dod < 0); sp_ero_veg = flat & ~pf & ch & ~chv & (dod < 0); sp_dep = hill & (dod > lod) & ((lv >= q75) | tree)
    counts = {k: int(v.sum()) for k, v in (("cells", np.isfinite(dod)), ("beyond_lod", np.abs(dod) > lod), ("spurious_erosion_zero_penetration", sp_ero_pf), ("spurious_erosion_vegetation_lod", sp_ero_veg), ("spurious_deposition_hillslope_vegetation", sp_dep))}
    counts["remaining_beyond_lod_after_outlines"] = int(np.sum((np.abs(dod) > lod) & ~sp_ero_pf & ~sp_ero_veg & ~sp_dep)); json.dump(counts, open(f"{V}/spurious_{name}.json", "w"), indent=1); print(name, counts)
    img = wms(bb); fig, axs = plt.subplots(1, 3, figsize=(21, 7.4 * (Y1 - Y0) / (X1 - X0) + 0.8))
    if img is not None: axs[0].imshow(img, extent=ext)
    overlays(axs[0], bb)
    im = axs[1].imshow(np.ma.masked_invalid(dod), extent=ext, origin="lower", cmap=cmap, norm=norm, interpolation="nearest"); overlays(axs[1], bb)
    axs[2].imshow(np.ma.masked_invalid(dod), extent=ext, origin="lower", cmap=cmap, norm=norm, interpolation="nearest", alpha=0.55); overlays(axs[2], bb, lw=0.8)
    plt.rcParams["hatch.linewidth"] = 0.6
    counts["drawn_clusters_cells"] = dict(zero_penetration=outline(axs[2], sp_ero_pf, bb, "black", "dashed", "////"), vegetation_lod=outline(axs[2], sp_ero_veg, bb, "black", "dotted", "...."), deposition=outline(axs[2], sp_dep, bb, "black", "solid", "\\\\\\\\"))
    json.dump(counts, open(f"{V}/spurious_{name}.json", "w"), indent=1)
    cb = fig.colorbar(im, ax=axs[1:], shrink=0.7, pad=0.01); cb.set_label("2021 - 2008 (mm)")
    fig.savefig(f"figures/zoom_{name}.png", dpi=110, bbox_inches="tight"); plt.close(fig); print(f"figures/zoom_{name}.png")
