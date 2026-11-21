#!/usr/bin/env python
"""Map of the three points with returns below the surveyed ground (RS+ B 3, RS+ B 100, RS+ A 13): for each, a 50 m
window on NAIP 2021 (Planetary Computer, skipped if offline), the 1 m mean-intensity grid and the 1 m structure-class
grid, with the point in red, its line neighbors numbered, and an overview of all three on the field. (Andy, 2026-11-21.)
Outputs  figures/rtk_2026-10-25_floodplain_three_points_map.png
"""
import os
import numpy as np, pandas as pd, rasterio
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")]
cases = [("floodplain_RSplus_B", 3, "RS+ B 3: low intensity overall"), ("floodplain_RSplus_B", 100, "RS+ B 100: high intensity only below the surveyed ground"), ("floodplain_RSplus_A", 13, "RS+ A 13: high intensity overall")]
with rasterio.open("data/field/2026-10-25/intensity/intensity_2021_1m_mean.tif") as r: INT = r.read(1).astype(float); INT[INT == r.nodata] = np.nan; ext = (r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top)
with rasterio.open("data/field/2026-10-25/intensity/cell_profiles_cluster.tif") as r: CL = r.read(1).astype(float); CL[CL == 255] = np.nan
naip = None
try:
    import planetary_computer as pc
    from pystac_client import Client
    from pyproj import Transformer
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds
    E0, E1, N0, N1 = ext[0], ext[1], ext[2], ext[3]; tr = Transformer.from_crs("EPSG:6344", "EPSG:4326", always_xy=True); lon0, lat0 = tr.transform(E0, N0); lon1, lat1 = tr.transform(E1, N1)
    cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace); it = [i for i in cat.search(collections=["naip"], bbox=[lon0, lat0, lon1, lat1]).items() if i.datetime.year == 2021][0]
    with rasterio.open(it.assets["image"].href) as rr:
        w = from_bounds(*transform_bounds("EPSG:6344", rr.crs, E0, N0, E1, N1), transform=rr.transform); naip = np.moveaxis(rr.read([1, 2, 3], window=w), 0, -1)
except Exception as e: print("NAIP not fetched:", repr(e)[:80])
tab = ListedColormap([plt.get_cmap("tab10")(k) for k in range(6)]); H = 25
fig, axs = plt.subplots(3, 4, figsize=(20, 14.5))
for i, (g, p, title) in enumerate(cases):
    r = F[(F.group == g) & (F.point == p)].iloc[0]; e, n = r.easting_tied, r.northing_tied; win = (e - H, e + H, n - H, n + H)
    line = F[(F.group == g) & F.point.between(p - 4, p + 4)]
    panels = [(naip, "NAIP 2021", {}), (INT, "2021 mean intensity, 1 m", dict(cmap="gray", vmin=1500, vmax=4800)), (CL, "1 m structure class", dict(cmap=tab, vmin=-0.5, vmax=5.5, interpolation="nearest"))]
    for j, (img, lab, kw) in enumerate(panels):
        ax = axs[i, j]
        if img is not None: ax.imshow(img, extent=ext, **kw)
        ax.plot(F.easting_tied, F.northing_tied, ".", color="w", ms=5, mec="k", mew=0.4)
        for _, q in line.iterrows(): ax.annotate(str(int(q.point)), (q.easting_tied, q.northing_tied), xytext=(3, 3), textcoords="offset points", fontsize=8, color="w" if j < 2 else "k", path_effects=None)
        ax.plot([e], [n], "o", ms=11, mfc="none", mec="red", mew=2); ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
        if j == 0: ax.set_title(title, fontsize=10, loc="left")
        else: ax.set_title(lab, fontsize=9)
    ax = axs[i, 3]; ax.imshow(INT, extent=ext, cmap="gray", vmin=1500, vmax=4800); ax.plot(F.easting_tied, F.northing_tied, ".", color="#e69f00", ms=2)
    for gg, pp, _ in cases:
        rr_ = F[(F.group == gg) & (F.point == pp)].iloc[0]; ax.plot([rr_.easting_tied], [rr_.northing_tied], "o", ms=9, mfc="none", mec="red", mew=1.5)
    ax.plot([e], [n], "o", ms=14, mfc="none", mec="cyan", mew=2); ax.set_xlim(579860, 580100); ax.set_ylim(4885330, 4885700); ax.set_aspect("equal"); ax.set_title("overview (this point in cyan)", fontsize=9); ax.tick_params(labelsize=7)
from matplotlib.patches import Patch
names = ["smooth sward", "low sward", "tall herb", "sward under crowns", "shrub", "tree"]
axs[2, 2].legend(handles=[Patch(color=plt.get_cmap("tab10")(k), label=names[k]) for k in range(6)], fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_three_points_map.png", dpi=int(os.environ.get("FIGDPI", 150))); print("wrote figures/rtk_2026-10-25_floodplain_three_points_map.png")
