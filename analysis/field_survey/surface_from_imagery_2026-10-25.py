"""Characterize the surfaces under the two hillslope RTK lines from satellite / aerial imagery (Andy, 2026-11-17).
Sentinel-2 L2A (10 m, Planetary Computer, same access path as mn-3dep-audit): NDVI at every hillslope FIX point
for every clear scene in three windows: around the survey (2026-10-10 .. 2026-11-10, canopy mostly down, so
green = conifers + shrub layer), spring leaf-off (2026-04-01 .. 2026-05-05) and midsummer (2026-07-01 .. 2026-08-15).
NAIP (0.6 m aerial, Planetary Computer): the most recent chip over the lines saved as PNG with the points drawn.
Positions: upper > 275 m NAVD88, mid 245-275, toe < 245 (the classes used in the README).
Output: data/field/2026-10-25/imagery/ (CSV of per-point NDVI, PNG chips) and this log on stdout."""
import pandas as pd, numpy as np, datetime as dt, sys, os
import planetary_computer as pc
from pystac_client import Client
import odc.stac
OUT = "data/field/2026-10-25/imagery"
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy()
H["pos"] = np.where(H.navd88_tied_m > 275, "upper", np.where(H.navd88_tied_m < 245, "toe", "mid"))
H["team"] = H.group.map({"hillslope_T4a": "T4a", "hillslope_RS3": "T3a"})
lon0, lat0 = H.lon.mean(), H.lat.mean(); buf = 400
dlat = buf / 111320; dlon = buf / (111320 * np.cos(np.radians(lat0)))
bbox = [lon0 - dlon, lat0 - dlat, lon0 + dlon, lat0 + dlat]
print("INPUT rtk_points_all.csv hillslope FIX points:", len(H), "| bbox", [round(b, 5) for b in bbox])
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
rows = []
for label, t0, t1 in (("survey_2026-10-25", "2026-10-10", "2026-11-10"), ("leafoff_spring_2026", "2026-04-01", "2026-05-05"), ("summer_2026", "2026-07-01", "2026-08-15")):
    items = list(cat.search(collections=["sentinel-2-l2a"], bbox=bbox, datetime=f"{t0}/{t1}", query={"eo:cloud_cover": {"lt": 40}}).items())
    print(f"\n{label}: {len(items)} Sentinel-2 L2A scenes with cloud < 40 %")
    for it in sorted(items, key=lambda i: i.properties["datetime"]):
        ds = odc.stac.load([it], bands=["red", "nir", "SCL"], bbox=bbox, resolution=10, chunks={})
        red = ds.red.astype("float32").squeeze().values; nir = ds.nir.astype("float32").squeeze().values; scl = ds.SCL.squeeze().values
        ndvi = (nir - red) / (nir + red)
        clear = ~np.isin(scl, [0, 1, 3, 8, 9, 10, 11])
        # sample at the tied positions: the loaded grid is in the item's UTM CRS (WGS84 UTM 15N); NAD83(2011) differs by ~1 m, below the 10 m pixel
        xs = ds.x.values; ys = ds.y.values
        j = np.clip(np.searchsorted(xs, H.easting_tied.values) - 1, 0, len(xs) - 1); i = np.clip(len(ys) - np.searchsorted(ys[::-1], H.northing_tied.values), 0, len(ys) - 1)
        v = ndvi[i, j]; c = clear[i, j]
        nclear = int(c.sum())
        print(f"  {it.properties['datetime'][:10]} cloud {it.properties.get('eo:cloud_cover', np.nan):5.1f} %  clear at {nclear}/{len(H)} points", end="")
        if nclear >= len(H) // 2:
            for (tm, ps), s in H.assign(ndvi=v, clear=c).groupby(["team", "pos"]):
                s = s[s.clear]
                if len(s): print(f" | {tm} {ps} {np.median(s.ndvi):.3f} (n {len(s)})", end="")
        print()
        for k, (idx, r) in enumerate(H.iterrows()):
            rows.append(dict(window=label, scene=it.id, date=it.properties["datetime"][:10], cloud=it.properties.get("eo:cloud_cover"), team=r.team, pos=r.pos, point=r.point, ndvi=float(v[k]), clear=bool(c[k]), navd88=r.navd88_tied_m, slope=r.slope_deg))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/sentinel2_ndvi_per_point.csv", index=False)
print("\nPER-SCENE SUMMARY, median NDVI by team and position, clear pixels only, scenes clear at >= half the points:")
S = T[T.clear].groupby(["window", "date", "team", "pos"]).ndvi.agg(["median", "count"]).reset_index()
pd.set_option("display.width", 250); print(S.pivot_table(index=["window", "date"], columns=["team", "pos"], values="median").round(3).to_string())
# NAIP chip
items = list(cat.search(collections=["naip"], bbox=bbox).items())
print(f"\nNAIP items over the lines: {len(items)}:", sorted({i.properties['datetime'][:10] for i in items}))
if items:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for it in sorted(items, key=lambda i: i.properties["datetime"])[-2:]:
        ds = odc.stac.load([it], bands=["image"], bbox=bbox, chunks={})
        img = ds.image.squeeze().values  # bands x y x
        xs = ds.x.values; ys = ds.y.values
        rgb = np.moveaxis(img[:3], 0, -1); cir = np.moveaxis(img[[3, 0, 1]], 0, -1)
        for name, arr in (("rgb", rgb), ("cir", cir)):
            fig, ax = plt.subplots(figsize=(9, 9)); ax.imshow(arr, extent=[xs[0], xs[-1], ys[-1], ys[0]])
            for tm, mk, col in (("T4a", "o", "#ff7f0e"), ("T3a", "s", "#1f77b4")):
                s = H[H.team == tm]; ax.scatter(s.easting_tied, s.northing_tied, s=28, marker=mk, facecolors="none", edgecolors=col, linewidths=1.2, label=tm)
            ax.legend(loc="lower right", frameon=True); ax.set_xticks([]); ax.set_yticks([])
            fn = f"{OUT}/naip_{it.properties['datetime'][:10]}_{name}.png"; fig.savefig(fn, dpi=150, bbox_inches="tight"); plt.close(fig); print("  wrote", fn)
