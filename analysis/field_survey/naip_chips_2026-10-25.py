"""NAIP over the two hillslope RTK lines (Planetary Computer COGs read through rasterio.merge so a box that crosses
a tile edge is mosaicked, not shifted): RGB and color-infrared chips with the tied points drawn, for leaf-on
2021-08-15 (the lidar year), 2025-09-28/29 (latest) and the two autumn acquisitions 2015-11-03 and 2017-10-16;
and, per point, the NAIP NDVI averaged within 3 m, summarized by team and slope position (upper > 275 m NAVD88,
mid 245-275, toe < 245). In November and mid-October green under a bare canopy is conifer or shrub layer."""
import pandas as pd, numpy as np, planetary_computer as pc, rasterio
from rasterio.merge import merge
from pystac_client import Client
from pyproj import Transformer
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "data/field/2026-10-25/imagery"
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy(); H["team"] = H.group.map({"hillslope_T4a": "T4a", "hillslope_RS3": "T3a"})
H["pos"] = np.where(H.navd88_tied_m > 275, "upper", np.where(H.navd88_tied_m < 245, "toe", "mid"))
E0, N0 = H.easting_tied.mean(), H.northing_tied.mean(); half = 220
lon0, lat0 = H.lon.mean(), H.lat.mean(); dlat = half / 111320; dlon = half / (111320 * np.cos(np.radians(lat0)))
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
items = list(cat.search(collections=["naip"], bbox=[lon0 - dlon, lat0 - dlat, lon0 + dlon, lat0 + dlat]).items())
rows = []
for want in ("2015-11-03", "2017-10-16", "2021-08-15", "2025-09-29", "2025-09-28"):
    its = [i for i in items if i.properties["datetime"][:10] == want]
    if want == "2025-09-29": its += [i for i in items if i.properties["datetime"][:10] == "2025-09-28"]
    if want == "2025-09-28" or not its: continue
    srcs = [rasterio.open(i.assets["image"].href) for i in its]
    crs = srcs[0].crs; T = Transformer.from_crs("EPSG:6344", crs, always_xy=True)
    cx, cy = T.transform(E0, N0); b = (cx - half, cy - half, cx + half, cy + half)
    img, tr = merge(srcs, bounds=b); [s.close() for s in srcs]
    ext = [tr.c, tr.c + tr.a * img.shape[2], tr.f + tr.e * img.shape[1], tr.f]
    px, py = T.transform(H.easting_tied.values, H.northing_tied.values)
    print(want, "crs", crs, "pixel", round(tr.a, 2), "m, chip", img.shape, "from", len(its), "tile(s)")
    # NDVI within 3 m of each point
    r, nir = img[0].astype(float), img[3].astype(float); ndvi = (nir - r) / np.where(nir + r == 0, np.nan, nir + r)
    yy, xx = np.indices(ndvi.shape); gx = tr.c + (xx + 0.5) * tr.a; gy = tr.f + (yy + 0.5) * tr.e
    for k, (idx, h) in enumerate(H.iterrows()):
        m = (gx - px[k]) ** 2 + (gy - py[k]) ** 2 <= 3.0 ** 2
        rows.append(dict(date=want, team=h.team, pos=h.pos, point=h.point, ndvi3m=float(np.nanmean(ndvi[m])), n=int(m.sum())))
    for name, arr in (("rgb", np.moveaxis(img[:3], 0, -1)), ("cir", np.moveaxis(img[[3, 0, 1]], 0, -1))):
        fig, ax = plt.subplots(figsize=(9, 9)); ax.imshow(arr, extent=ext)
        for tm, mk, col in (("T4a", "o", "#ff7f0e"), ("T3a", "s", "#00c3ff")):
            s = H[H.team == tm]; ax.scatter(px[(H.team == tm).values], py[(H.team == tm).values], s=40, marker=mk, facecolors="none", edgecolors=col, linewidths=1.4, label=tm)
        ax.legend(loc="lower right"); ax.set_xticks([]); ax.set_yticks([])
        fn = f"{OUT}/naip_{want}_{name}.png"; fig.savefig(fn, dpi=150, bbox_inches="tight"); plt.close(fig); print("  wrote", fn)
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/naip_ndvi_3m_per_point.csv", index=False)
print("\nNAIP NDVI within 3 m of each point, median by team and position [n points]:")
from scipy.stats import mannwhitneyu
for d, s in R.groupby("date"):
    line = f"  {d}: " + " | ".join(f"{tm} {ps} {g.ndvi3m.median():.3f} [{len(g)}]" for (tm, ps), g in s.groupby(["team", "pos"]))
    a = s[(s.team == "T4a") & (s.pos == "upper")].ndvi3m; bb = s[(s.team == "T3a") & (s.pos == "mid")].ndvi3m
    print(line + f" || T4a upper minus T3a mid {a.median() - bb.median():+.3f}, Mann-Whitney p {mannwhitneyu(a, bb).pvalue:.3f}")
