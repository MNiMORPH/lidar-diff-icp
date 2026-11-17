"""Addendum to surface_from_imagery_2026-10-25.py: (1) NAIP chips for the fall leaf-off-era acquisitions (2015-11-03,
2017-10-16) and the lidar year (2021-08-15), points drawn; (2) a per-point test of the late-October NDVI difference
between T4a's upper slope and T3a's mid-slope (Mann-Whitney, two-sided), per scene."""
import pandas as pd, numpy as np, planetary_computer as pc, odc.stac
from pystac_client import Client
from scipy.stats import mannwhitneyu
OUT = "data/field/2026-10-25/imagery"
T = pd.read_csv(f"{OUT}/sentinel2_ndvi_per_point.csv"); T = T[T.clear]
print("T4a upper (n 14) vs T3a mid (n 17), per-point NDVI: median, p16-p84, Mann-Whitney p")
for (w, d), s in T.groupby(["window", "date"]):
    a = s[(s.team == "T4a") & (s.pos == "upper")].ndvi; b = s[(s.team == "T3a") & (s.pos == "mid")].ndvi
    if len(a) < 10 or len(b) < 10: continue
    p = mannwhitneyu(a, b).pvalue
    print(f"  {w:20s} {d}  T4a upper {a.median():.3f} ({a.quantile(.16):.3f}-{a.quantile(.84):.3f})  T3a mid {b.median():.3f} ({b.quantile(.16):.3f}-{b.quantile(.84):.3f})  diff {a.median()-b.median():+.3f}  p {p:.3f}")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy(); H["team"] = H.group.map({"hillslope_T4a": "T4a", "hillslope_RS3": "T3a"})
lon0, lat0 = H.lon.mean(), H.lat.mean(); buf = 300; dlat = buf / 111320; dlon = buf / (111320 * np.cos(np.radians(lat0)))
bbox = [lon0 - dlon, lat0 - dlat, lon0 + dlon, lat0 + dlat]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
items = list(cat.search(collections=["naip"], bbox=bbox).items())
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
for want in ("2015-11-03", "2017-10-16", "2021-08-15"):
    its = [i for i in items if i.properties["datetime"][:10] == want]
    if not its: print("no NAIP", want); continue
    ds = odc.stac.load(its[:1], bands=["image"], bbox=bbox, chunks={}); img = ds.image.squeeze().values; xs = ds.x.values; ys = ds.y.values
    for name, arr in (("rgb", np.moveaxis(img[:3], 0, -1)), ("cir", np.moveaxis(img[[3, 0, 1]], 0, -1))):
        fig, ax = plt.subplots(figsize=(9, 9)); ax.imshow(arr, extent=[xs[0], xs[-1], ys[-1], ys[0]])
        for tm, mk, col in (("T4a", "o", "#ff7f0e"), ("T3a", "s", "#1f77b4")):
            s = H[H.team == tm]; ax.scatter(s.easting_tied, s.northing_tied, s=34, marker=mk, facecolors="none", edgecolors=col, linewidths=1.3, label=tm)
        ax.legend(loc="lower right"); ax.set_xticks([]); ax.set_yticks([])
        fn = f"{OUT}/naip_{want}_{name}.png"; fig.savefig(fn, dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", fn)
