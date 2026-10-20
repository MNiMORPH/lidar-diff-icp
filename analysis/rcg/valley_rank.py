#!/usr/bin/env python3
"""Phase A, stage 2 (no network): classify and rank candidate RTK sites across the Whitewater
valley from imagery alone, on public land, near access.

Inputs (from valley_fetch.py, data/derived/valley/): ten-band Sentinel-2 mosaics for
2026-10-17, 2025-10-27, 2021-10-18; NDVI mosaics for 2026-04-10/04-15/05-15 and the 2026
cut-history dates; Landsat-5 2008-10-29; WorldCover; WMA units, state parks, parking,
roads, river centrelines.

METHOD, each step as it was done on the Elba tile (RCG_PHENOLOGY.md sections N, P, U-W):
  valley floor   the 800 m centreline buffer, the plan's proxy (D1 is Andy's); 'open' =
                 not WorldCover tree / built / water inside it.
  presence       per October scene, two-class maximum likelihood on all bands, trained on
                 the NE patch's pixels (the Elba tile's lidar-found field, its mask carried
                 onto this grid) against the rest of the open valley floor; p > 0.5 =
                 NE-like THAT YEAR. Bands are normalised per scene by the open valley floor's
                 median so the classifier sees relative reflectance. Confusion check by
                 WorldCover printed.
  density        spring 2026 vigour = mean NDVI of 04-10, 04-15, 05-15 (section W).
  cut history    per pixel over the 2026 Jun-Oct NDVI series: a CUT = a drop of >= 0.25
                 between consecutive clear dates (section V: a cut halves NDVI by the next
                 pass); the last cut date is kept.
  classes        R1 NE-like in 2026 AND 2025; R2 NE-like in exactly one; S late-green
                 herbaceous not NE-like (Oct-2026 NDVI >= 0.30, WorldCover grassland or
                 herbaceous wetland); C bare/cropped control (Oct-2026 NDVI < 0.25 and
                 summer-2026 max NDVI < 0.45, i.e. not a crop in leaf... crops ARE in leaf
                 in summer, so C = Oct NDVI < 0.25 only, flagged with its summer max);
                 E forest edge (open pixel within 30 m of WorldCover tree). All cut-offs
                 are MINE (D4) and printed.
  sites          connected clusters >= 50 pixels (0.5 ha) on PUBLIC land; per cluster:
                 centroid (UTM, lat/lon), ha, unit, % of the cluster inside the TRAINING
                 patch (a site there is resubstitution, not a find), class memberships
                 2008/2021/2025/2026 (2008 NaN until the Landsat scene is on disk),
                 spring vigour, last 2026 cut date, distance to nearest WMA parking and
                 road; ranked within class by hectares / (1 + distance to parking in km).
Outputs: data/derived/valley/valley_sites.csv, valley_classes.tif (uint8 class map),
figures/valley_sites.png, and the per-year p rasters.

    ./lidar-icp/bin/python analysis/rcg/valley_rank.py
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import rioxarray
import shapely
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation, label, distance_transform_edt
from shapely.geometry import shape

V = "data/derived/valley"; PL = f"{V}/public_land"; FIG = "figures"; S = "data/derived/elbaext_delong_star_nosurf"
BUF_M = 800.0; MIN_PX = 50; CUT_DROP = 0.25
S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]; LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
ref = rioxarray.open_rasterio(f"{V}/s2_20261017_bands.tif"); xs, ys = ref.x.values, ref.y.values; ny, nx = len(ys), len(xs)
E, N = np.meshgrid(xs, ys); print(f"grid {ny} x {nx} at 10 m")
geo = lambda name: [shape(f["geometry"]) for f in json.load(open(f"{PL}/{name}"))["features"]]
props = lambda name: [f["properties"] for f in json.load(open(f"{PL}/{name}"))["features"]]
wma, wmap = geo("wma_units.geojson"), props("wma_units.geojson"); parks, parkp = geo("state_parks.geojson"), props("state_parks.geojson")
parking = [g for g, p in zip(geo("wma_parking.geojson"), props("wma_parking.geojson")) if p.get("facility_subgroup") == "Parking Lot"]; roads = geo("wma_roads.geojson"); rivers = geo("rivers_whitewater.geojson")
units = [(g, p.get("name") or p.get("unit_name")) for g, p in zip(wma, wmap)] + [(g, p.get("name") or p.get("area_name") or p.get("unit_name") or "state park") for g, p in zip(parks, parkp)]
public = shapely.union_all([g for g, _ in units]); valley = shapely.union_all(rivers).buffer(BUF_M)
def rasterize(geom): return shapely.contains_xy(geom, E.ravel(), N.ravel()).reshape(ny, nx)
pub = rasterize(public); val = rasterize(valley)
wc = rioxarray.open_rasterio(f"{V}/worldcover.tif").squeeze(drop=True).reindex(x=xs, y=ys, method="nearest").values
open_ = val & ~np.isin(wc, [10, 50, 80])
print(f"valley-floor proxy (800 m buffer): {val.sum()*0.01:,.0f} ha; public inside it: {(pub & val).sum()*0.01:,.0f} ha; open inside it: {open_.sum()*0.01:,.0f} ha")
# the NE patch (Elba) onto this grid
dodT = np.load(f"{S}/dod.npy"); newT = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fpT = np.load(f"{S}/floodplain_mask.npy").astype(bool); lodT = np.load(f"{S}/lod.npy") * 1e3
below = fpT & np.isfinite(newT) & (newT < -lodT); labT, nT = label(below, structure=np.ones((3, 3))); sz = np.bincount(labT.ravel()); NE = labT == (np.argsort(sz[1:])[::-1][0] + 1)
TX0, TY0, tres = 575600.0, 4882200.0, 5.0; jy, jx = np.nonzero(NE); nepts = shapely.points(TX0 + (jx + 0.5) * tres, TY0 + (jy + 0.5) * tres)
ne_poly = shapely.union_all(shapely.buffer(nepts, 3.6, cap_style="square")); ne = rasterize(ne_poly); print(f"NE patch on the valley grid: {ne.sum()} px")

def logpdf(Z, C):
    mu = C.mean(0); cov = np.cov(C.T) + np.eye(C.shape[1]) * 1e-6; inv = np.linalg.inv(cov); dz = Z - mu
    return -0.5 * np.einsum("ij,jk,ik->i", dz, inv, dz) - 0.5 * np.linalg.slogdet(cov)[1]
def classify(tif, bn, tag):
    da = rioxarray.open_rasterio(tif)
    if not (np.allclose(da.x.values, xs) and np.allclose(da.y.values, ys)): da = da.reindex(x=xs, y=ys, method="nearest")
    X = np.stack([da.values[i] for i in range(len(bn))], -1); fin = np.all(np.isfinite(X), -1)
    refm = np.nanmedian(X[open_ & fin], 0); Xr = X / refm
    c1 = Xr[ne & fin]; c0 = Xr[open_ & ~ne & fin]
    p = np.full((ny, nx), np.nan); Zf = Xr[fin]; p[fin] = 1 / (1 + np.exp(np.clip(logpdf(Zf, c0) - logpdf(Zf, c1), -50, 50)))
    like = p > 0.5
    print(f"  {tag}: trained on {len(c1):,} NE px vs {len(c0):,} open-valley px; NE resubstitution {100*np.nanmean(like[ne & fin]):.0f}%; NE-like on open valley floor {like[open_ & fin].sum()*0.01:,.0f} ha ({100*np.nanmean(like[open_ & fin]):.1f}%); on public open valley {like[open_ & pub & fin].sum()*0.01:,.0f} ha")
    print("    confusion by WorldCover (whole box):", ", ".join(f"{nm} {100*np.nanmean(like[(wc == k) & fin]):.1f}%" for k, nm in ((10, "tree"), (30, "grass"), (40, "crop"), (90, "wetland")) if ((wc == k) & fin).sum() > 100))
    da.isel(band=0).copy(data=p.astype("float32")).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{V}/p_nelike_{tag}.tif", compress="deflate")
    nd = (X[..., bn.index("B08") if "B08" in bn else bn.index("nir08")] - X[..., bn.index("B04") if "B04" in bn else bn.index("red")]) / (X[..., bn.index("B08") if "B08" in bn else bn.index("nir08")] + X[..., bn.index("B04") if "B04" in bn else bn.index("red")])
    return like, nd
print("PRESENCE classifiers:")
like26, ndvi26 = classify(f"{V}/s2_20261017_bands.tif", S2B, "2026")
like25, _ = classify(f"{V}/s2_20251027_bands.tif", S2B, "2025")
like21, _ = classify(f"{V}/s2_20211018_bands.tif", S2B, "2021")
like08 = np.full((ny, nx), np.nan)   # NaN until the Landsat 2008 scene is on disk: absent is not "not NE-like"
if os.path.exists(f"{V}/ls_20081029_bands.tif"):
    like08, _ = classify(f"{V}/ls_20081029_bands.tif", LSB, "2008"); like08 = like08.astype(float)
else: print("  2008: Landsat scene not on disk; nelike_2008 is NaN")
# density and cut history
nd = lambda f: rioxarray.open_rasterio(f).squeeze(drop=True).reindex(x=xs, y=ys, method="nearest").values
spring = np.nanmean(np.stack([nd(f"{V}/s2_{d}_ndvi.tif") for d in ("20260410", "20260415", "20260515") if os.path.exists(f"{V}/s2_{d}_ndvi.tif")]), 0)
cutfiles = sorted(f for f in glob.glob(f"{V}/s2_2026*_ndvi.tif") if os.path.basename(f)[3:11] >= "20260615")
dates = [os.path.basename(f)[3:11] for f in cutfiles]; series = np.stack([nd(f) for f in cutfiles]) if cutfiles else np.zeros((0, ny, nx))
summer_max = np.nanmax(series, 0) if len(series) else np.full((ny, nx), np.nan); last_cut = np.full((ny, nx), "", dtype=object); last_cut_doy = np.full((ny, nx), np.nan)
prev = None; prev_d = None
for k, d in enumerate(dates):
    cur = series[k]
    if prev is not None:
        drop = (prev - cur) >= CUT_DROP; last_cut[drop] = d; last_cut_doy[drop] = pd.Timestamp(d).dayofyear
    both = np.isfinite(cur)
    prev = np.where(both, cur, prev if prev is not None else cur); prev_d = d
print(f"cut-history dates used: {dates}; pixels with a 2026 cut on open valley floor: {100*np.mean(np.isfinite(last_cut_doy[open_])):.1f}%")
# classes
near_tree = binary_dilation(wc == 10, structure=np.ones((3, 3)), iterations=3) & (wc != 10)
cls = np.zeros((ny, nx), np.uint8)   # 0 none, 1 R1, 2 R2, 3 S, 4 C, 5 E
R1 = open_ & like26 & like25; R2 = open_ & (like26 ^ like25); Sc = open_ & ~like26 & ~like25 & (ndvi26 >= 0.30) & np.isin(wc, [30, 90]); Cc = open_ & (ndvi26 < 0.25) & ~like26; Ec = open_ & near_tree & ~R1 & ~R2
for v, m in ((5, Ec), (4, Cc), (3, Sc), (2, R2), (1, R1)): cls[m] = v
ref.isel(band=0).copy(data=cls.astype("float32")).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{V}/valley_classes.tif")
names = {1: "R1 NE-like 2025 and 2026", 2: "R2 NE-like one of 2025/2026", 3: "S late-green herbaceous, not NE-like", 4: "C bare / low October NDVI", 5: "E forest edge"}
print("\nCLASS AREAS (ha): open valley floor / public open valley floor")
for v, nm in names.items(): print(f"  {nm:<40} {(cls == v).sum()*0.01:>8,.0f} / {((cls == v) & pub).sum()*0.01:>8,.0f}")
# sites on public land
dpark = distance_transform_edt(~rasterize(shapely.union_all([g.buffer(5) for g in parking]))) * 10.0; droad = distance_transform_edt(~rasterize(shapely.union_all([g.buffer(5) for g in roads]))) * 10.0
rows = []
for v, nm in names.items():
    m = (cls == v) & pub; lab_, n_ = label(m, structure=np.ones((3, 3))); szs = np.bincount(lab_.ravel())
    for k in np.flatnonzero(szs >= MIN_PX):
        if k == 0: continue
        mm = lab_ == k; e0, n0 = E[mm].mean(), N[mm].mean(); lon, lat = to_ll.transform(e0, n0); pt = shapely.Point(e0, n0)
        unit = next((u for g, u in units if g.contains(pt)), "?"); cuts = last_cut[mm]; cuts = cuts[cuts != ""]
        rows.append(dict(cls=nm, unit=unit, ha=round(mm.sum() * 0.01, 2), easting=round(e0, 1), northing=round(n0, 1), lat=round(lat, 6), lon=round(lon, 6),
                         training_patch_pct=round(100 * ne[mm].mean(), 0), nelike_2008=round(100 * np.nanmean(like08[mm]), 0), nelike_2021=round(100 * np.nanmean(like21[mm]), 0), nelike_2025=round(100 * np.nanmean(like25[mm]), 0), nelike_2026=round(100 * np.nanmean(like26[mm]), 0),
                         oct2026_ndvi=round(float(np.nanmedian(ndvi26[mm])), 2), spring2026_ndvi=round(float(np.nanmedian(spring[mm])), 2), summer2026_max_ndvi=round(float(np.nanmedian(summer_max[mm])), 2),
                         cut_2026=(pd.Series(cuts).mode().iloc[0] if len(cuts) else ""), cut_share=round(100 * len(cuts) / mm.sum(), 0), dist_parking_m=round(float(dpark[mm].min()), 0), dist_road_m=round(float(droad[mm].min()), 0)))
F = pd.DataFrame(rows); F["rank_score"] = F.ha / (1 + F.dist_parking_m / 1000.0); F = F.sort_values(["cls", "rank_score"], ascending=[True, False]); F.to_csv(f"{V}/valley_sites.csv", index=False)
print(f"\nCANDIDATE SITES on public land (clusters >= {MIN_PX} px = 0.5 ha): {len(F)}; per class: {F.cls.value_counts().to_dict()}")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
for v, nm in names.items():
    s = F[F.cls == nm].head(6)
    if len(s): print(f"\n  top {nm}:"); print(s.drop(columns=["cls", "rank_score"]).to_string(index=False))
# map
fig, ax = plt.subplots(figsize=(14, 13)); cm = ListedColormap(["#f2f2f0", "#8b0000", "#eb6834", "#4a3aa7", "#eda100", "#1baf7a"])
ax.imshow(np.where(val, cls, np.nan), extent=[xs[0], xs[-1], ys[-1], ys[0]], origin="upper", cmap=cm, vmin=-0.5, vmax=5.5, interpolation="nearest")
for g, u in units:
    for poly in (g.geoms if hasattr(g, "geoms") else [g]): px, py = poly.exterior.xy; ax.plot(px, py, color="black", lw=0.8)
for r in rivers: px, py = r.xy; ax.plot(px, py, color="#2a78d6", lw=0.8)
for p in parking: ax.plot(p.x, p.y, marker="P", color="yellow", ms=7, mec="black")
for _, r in F.groupby("cls").head(4).iterrows(): ax.plot(r.easting, r.northing, marker="o", ms=9, mfc="none", mec="black", mew=1.5)
ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[-1], ys[0]); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Whitewater valley, imagery-only classes on the 800 m valley-floor proxy: dark red R1, orange R2, purple S, yellow C, green E; black = public land; yellow = WMA parking; circles = top-4 sites per class", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FIG}/valley_sites.png", dpi=110); print(f"\nfigure: {FIG}/valley_sites.png; table: {V}/valley_sites.csv; classes: {V}/valley_classes.tif")
