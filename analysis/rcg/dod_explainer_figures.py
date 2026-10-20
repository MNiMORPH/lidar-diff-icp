#!/usr/bin/env python3
"""Minimal DoD figures for the field day (Andy, 2026-10-20): "showing the lidar differences for
explaining why they are nonphysical. 2021 - 2008 ... Minimal: locations or area outlines only."

DoD = gen2 (3DEP, 1 May 2021) - gen1 (MN DNR, 25 Nov 2008), the star-only product
(data/derived/elbaext_delong_star_nosurf/dod.npy, m -> mm), stable-ground median removed as
the gauge; diverging blue (2021 lower) / red (2021 higher), +/- 300 mm, grey = no data.
Overlays: survey-region outlines and the class transect from field_day_waypoints.gpx, the two
parking lots, the two monuments; no text on the maps.
  figures/dod_explainer_tile.png     the whole Elba tile with a hillshade under the DoD
  figures/dod_explainer_lot204.png   zoom, stop 1
  figures/dod_explainer_lot501.png   zoom, stop 2
"""
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyproj
from matplotlib.colors import LightSource, TwoSlopeNorm

S = "data/derived/elbaext_delong_star_nosurf"; X0, Y0, RES = 575600.0, 4882200.0, 5.0
dod = np.load(f"{S}/dod.npy") * 1000.0; st = np.load(f"{S}/stable.npy").astype(bool); z = np.load(f"{S}/z_after.npy")
dod = dod - np.nanmedian(dod[st & np.isfinite(dod)]); ny, nx = dod.shape; ext = [X0, X0 + nx * RES, Y0, Y0 + ny * RES]
from_ll = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
root = ET.parse("data/derived/valley/field_day_waypoints.gpx").getroot(); ns = {"g": "http://www.topografix.com/GPX/1/1"}
tracks = [(t.find("g:name", ns).text, np.array([from_ll.transform(float(p.get("lon")), float(p.get("lat"))) for p in t.findall(".//g:trkpt", ns)])) for t in root.findall("g:trk", ns)]
wpts = {w.find("g:name", ns).text: from_ll.transform(float(w.get("lon")), float(w.get("lat"))) for w in root.findall("g:wpt", ns)}
lots = [wpts[k] for k in ("lot_204", "lot_501")]; mons = [v for k, v in wpts.items() if k.startswith("DR9925") or k.startswith("DG8385")]
hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(z, nan=np.nanmean(z)), vert_exag=1.5, dx=RES, dy=RES)
cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#bdbdbd"); norm = TwoSlopeNorm(vmin=-300, vcenter=0, vmax=300)

def panel(bb, path, figsize, bar=200, lw=1.2):
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(np.ma.masked_invalid(dod), extent=ext, origin="lower", cmap=cmap, norm=norm, alpha=0.75, interpolation="nearest")
    for name, xy in tracks:
        if name.startswith("region_"): ax.plot(xy[:, 0], xy[:, 1], color="black", lw=lw)
        elif name == "class_transect": ax.plot(xy[:, 0], xy[:, 1], color="black", lw=lw, ls="--")
    for x, y in lots: ax.plot(x, y, marker="P", ms=11, color="yellow", mec="black")
    for x, y in mons: ax.plot(x, y, marker="*", ms=14, color="magenta", mec="black")
    ax.plot([bb[0] + 0.04 * (bb[2] - bb[0]), bb[0] + 0.04 * (bb[2] - bb[0]) + bar], [bb[1] + 0.04 * (bb[3] - bb[1])] * 2, color="black", lw=4)
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3]); ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.01); cb.set_label("2021 - 2008 (mm)")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig); print(path)

panel((X0, Y0, X0 + nx * RES, Y0 + ny * RES), "figures/dod_explainer_tile.png", (11, 10), bar=500, lw=0.8)
panel((579400, 4884900, 580300, 4885900), "figures/dod_explainer_lot204.png", (9, 9.5))
panel((578900, 4885300, 579900, 4886300), "figures/dod_explainer_lot501.png", (9, 9.5))
