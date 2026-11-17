"""Land cover of the floodplain survey field, its history, and the point-to-point structure of lidar minus RTK
(Andy, 2026-11-17). INPUT: data/field/2026-10-25/rtk_points_all.csv, floodplain FIX points (tied, pole-corrected).
COVARIATES PER POINT: lowveg_gen2 (2021 lidar, from the table); gen1_spread_mm (2008 lidar per-cell return spread,
nearest 5 m cell of data/derived/<tile>/gen1_spread_mm.npy, the project's 2008 vegetation proxy); NAIP NDVI within
3 m for every NAIP year on Planetary Computer (2010-2025; uncalibrated digital numbers, relative within a date);
Sentinel-2 L2A NDVI (10 m) per clear scene around the 2021 lidar flight (2021-04-15..05-20; 3DEP flown 2021-05-01),
summer 2021, April 2026 and around the survey (2026-10-10..11-10). FIELD-LEVEL: Landsat C2 L2 NDVI (30 m) median over
the field for 2008 (the 2008 lidar was flown in November 2008) and 2021, by month. RESIDUALS: gen2_minus_rtk_m and
gen1_minus_rtk_m (lidar minus RTK); structure is compared after removing each rover's median (level is the rover's
tie; structure is the ground's). OUTPUT: data/field/2026-10-25/imagery/floodplain_* (CSV, PNG) and this log.
Written incrementally. No threshold is chosen here."""
import pandas as pd, numpy as np, os, json, glob, datetime as dt
import planetary_computer as pc, rasterio, odc.stac
from rasterio.merge import merge
from pystac_client import Client
from pyproj import Transformer
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
OUT = "data/field/2026-10-25/imagery"
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")].copy().reset_index(drop=True)
F["rover"] = F.group.map({"floodplain_CORS_Wickert": "CORS", "floodplain_RSplus_A": "RS+ A", "floodplain_RSplus_B": "RS+ B"})
F["r21_mm"] = 1000 * F.gen2_minus_rtk_m; F["r08_mm"] = 1000 * F.gen1_minus_rtk_m
for c in ("r21_mm", "r08_mm"): F[c + "_rel"] = F[c] - F.groupby("rover")[c].transform("median")
print("INPUT floodplain FIX points:", len(F), dict(F.rover.value_counts()))
print("COLUMNS: r21_mm / r08_mm = 2021 / 2008 lidar minus tied RTK; *_rel = minus the rover's median (structure only)")
# gen1 spread at the points
tiles = {}
for c in glob.glob("data/derived/4*/corrections.json"):
    t = os.path.basename(os.path.dirname(c)); j = json.load(open(c)); tiles[t] = (j["bounds"], float(j.get("res") or j.get("res_m")))
def nearest(name, x, y):
    for t, (b, res) in tiles.items():
        if b[0] <= x < b[2] and b[1] <= y < b[3]:
            p = f"data/derived/{t}/{name}"
            if not os.path.exists(p): return np.nan
            a = np.load(p, mmap_mode="r"); return float(a[int((y - b[1]) / res), int((x - b[0]) / res)])
    return np.nan
F["gen1_spread_mm"] = [nearest("gen1_spread_mm.npy", x, y) for x, y in zip(F.easting_tied, F.northing_tied)]
print("gen1_spread_mm at the points: median %.0f, p10 %.0f, p90 %.0f, NaN %d" % (F.gen1_spread_mm.median(), F.gen1_spread_mm.quantile(.1), F.gen1_spread_mm.quantile(.9), F.gen1_spread_mm.isna().sum()))
E0, E1, N0, N1 = F.easting_tied.min() - 60, F.easting_tied.max() + 60, F.northing_tied.min() - 60, F.northing_tied.max() + 60
to_ll = Transformer.from_crs("EPSG:6344", "EPSG:4326", always_xy=True); lo0, la0 = to_ll.transform(E0, N0); lo1, la1 = to_ll.transform(E1, N1); bbox = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
# ---- NAIP by year: per-point NDVI within 3 m, and an RGB chip per year
items = list(cat.search(collections=["naip"], bbox=bbox).items()); years = sorted({i.properties["datetime"][:4] for i in items})
print("\nNAIP years over the field:", years)
chips = {}
for yr in years:
    its = [i for i in items if i.properties["datetime"][:4] == yr]; dates = sorted({i.properties["datetime"][:10] for i in its})
    srcs = [rasterio.open(i.assets["image"].href) for i in its]; crs = srcs[0].crs; T = Transformer.from_crs("EPSG:6344", crs, always_xy=True)
    x0, y0 = T.transform(E0, N0); x1, y1 = T.transform(E1, N1)
    img, tr = merge(srcs, bounds=(x0, y0, x1, y1)); [s.close() for s in srcs]
    r, nir = img[0].astype(float), img[3].astype(float); ndvi = (nir - r) / np.where(nir + r == 0, np.nan, nir + r)
    yy, xx = np.indices(ndvi.shape); gx = tr.c + (xx + 0.5) * tr.a; gy = tr.f + (yy + 0.5) * tr.e
    px, py = T.transform(F.easting_tied.values, F.northing_tied.values)
    F[f"naip_{yr}"] = [float(np.nanmean(ndvi[(gx - a) ** 2 + (gy - b) ** 2 <= 9.0])) for a, b in zip(px, py)]
    chips[yr] = (np.moveaxis(img[:3], 0, -1), np.moveaxis(img[[3, 0, 1]], 0, -1), [tr.c, tr.c + tr.a * img.shape[2], tr.f + tr.e * img.shape[1], tr.f], px, py, dates)
    print(f"  {yr} {dates}: pixel {tr.a:.1f} m; NDVI at points median {F[f'naip_{yr}'].median():.3f} (p10 {F[f'naip_{yr}'].quantile(.1):.3f}, p90 {F[f'naip_{yr}'].quantile(.9):.3f}); field-chip median {np.nanmedian(ndvi):.3f}")
    F.to_csv(f"{OUT}/floodplain_points_covariates.csv", index=False)
# ---- Sentinel-2 per point
def s2(label, t0, t1):
    its = list(cat.search(collections=["sentinel-2-l2a"], bbox=bbox, datetime=f"{t0}/{t1}", query={"eo:cloud_cover": {"lt": 40}}).items())
    print(f"\nSentinel-2 {label} ({t0}..{t1}): {len(its)} scenes with cloud < 40 %")
    for it in sorted(its, key=lambda i: i.properties["datetime"]):
        ds = odc.stac.load([it], bands=["red", "nir", "SCL"], bbox=bbox, resolution=10, chunks={})
        red = ds.red.astype("float32").squeeze().values; nir = ds.nir.astype("float32").squeeze().values; scl = ds.SCL.squeeze().values
        ndvi = (nir - red) / (nir + red); clear = ~np.isin(scl, [0, 1, 3, 8, 9, 10, 11])
        xs = ds.x.values; ys = ds.y.values
        j = np.clip(np.searchsorted(xs, F.easting_tied.values) - 1, 0, len(xs) - 1); i = np.clip(len(ys) - np.searchsorted(ys[::-1], F.northing_tied.values), 0, len(ys) - 1)
        v = np.where(clear[i, j], ndvi[i, j], np.nan); d = it.properties["datetime"][:10]
        if np.isfinite(v).sum() >= len(F) // 2:
            F[f"s2_{d}"] = v; print(f"  {d} cloud {it.properties.get('eo:cloud_cover', np.nan):5.1f} %: NDVI at points median {np.nanmedian(v):.3f} (p10 {np.nanpercentile(v, 10):.3f}, p90 {np.nanpercentile(v, 90):.3f}), clear {int(np.isfinite(v).sum())}/{len(F)}")
        else: print(f"  {d} cloud {it.properties.get('eo:cloud_cover', np.nan):5.1f} %: clear at only {int(np.isfinite(v).sum())} points, skipped")
    F.to_csv(f"{OUT}/floodplain_points_covariates.csv", index=False)
s2("around the 2021 lidar flight (2021-05-01)", "2021-04-15", "2021-05-20"); s2("summer 2021", "2021-07-01", "2021-08-31"); s2("April 2026", "2026-04-01", "2026-05-05"); s2("around the survey", "2026-10-10", "2026-11-10")
# ---- Landsat field-level, 2008 and 2021
print("\nLandsat C2 L2 (30 m), median NDVI over clear pixels of the field box, by scene:")
L = []
for yr in ("2008", "2021"):
    its = list(cat.search(collections=["landsat-c2-l2"], bbox=bbox, datetime=f"{yr}-04-01/{yr}-12-31", query={"eo:cloud_cover": {"lt": 60}}).items())
    for it in sorted(its, key=lambda i: i.properties["datetime"]):
        try:
            ds = odc.stac.load([it], bands=["red", "nir08", "qa_pixel"], bbox=bbox, resolution=30, chunks={})
        except Exception as e: print("  ", it.id, "load failed:", str(e)[:60]); continue
        red = ds.red.astype("float32").squeeze().values; nir = ds.nir08.astype("float32").squeeze().values; qa = ds.qa_pixel.squeeze().values.astype(int)
        clear = (qa & (1 << 6)) > 0   # bit 6: clear (C2 QA_PIXEL)
        red = red * 2.75e-5 - 0.2; nir = nir * 2.75e-5 - 0.2; ndvi = (nir - red) / (nir + red)
        v = ndvi[clear & np.isfinite(ndvi)]
        if v.size >= 10:
            L.append(dict(date=it.properties["datetime"][:10], platform=it.properties.get("platform"), ndvi_median=float(np.median(v)), n=int(v.size), total=int(ndvi.size)))
            print(f"  {it.properties['datetime'][:10]} {it.properties.get('platform'):12s} NDVI median {np.median(v):.3f} over {v.size}/{ndvi.size} clear pixels")
pd.DataFrame(L).to_csv(f"{OUT}/floodplain_landsat_field_ndvi.csv", index=False)
# ---- structure vs covariates
print("\nSPEARMAN of lidar minus RTK with each covariate: per rover, and pooled after removing each rover's median (rho, p, n):")
covs = ["lowveg_gen2", "gen1_spread_mm"] + [c for c in F.columns if c.startswith("naip_") or c.startswith("s2_")]
for resid in ("r21_mm", "r08_mm"):
    print(f" {resid}:")
    for c in covs:
        parts = []
        for rv, s in F.groupby("rover"):
            s = s[[resid, c]].dropna()
            if len(s) >= 10: r = spearmanr(s[resid], s[c]); parts.append(f"{rv} {r.correlation:+.2f} (p {r.pvalue:.2f}, n {len(s)})")
        s = F[[resid + "_rel", c]].dropna(); r = spearmanr(s[resid + "_rel"], s[c])
        print(f"   {c:16s} pooled-rel {r.correlation:+.2f} (p {r.pvalue:.3f}, n {len(s)}) | " + "; ".join(parts))
F.to_csv(f"{OUT}/floodplain_points_covariates.csv", index=False)
# ---- figures
yrs = [y for y in years]; n = len(yrs); cols = 4; rows_ = int(np.ceil(n / cols))
fig, axs = plt.subplots(rows_, cols, figsize=(4 * cols, 4.6 * rows_)); axs = np.atleast_1d(axs).ravel()
for k, yr in enumerate(yrs):
    rgb, cir, ext, px, py, dates = chips[yr]; ax = axs[k]; ax.imshow(rgb, extent=ext); ax.scatter(px, py, s=3, c="k"); ax.set_xticks([]); ax.set_yticks([]); ax.text(0.02, 0.98, ", ".join(dates), transform=ax.transAxes, va="top", fontsize=8, color="w", bbox=dict(fc="k", alpha=0.5, lw=0))
for ax in axs[n:]: ax.axis("off")
fig.savefig(f"{OUT}/floodplain_naip_years.png", dpi=130, bbox_inches="tight"); plt.close(fig); print("wrote", f"{OUT}/floodplain_naip_years.png")
rgb, cir, ext, px, py, dates = chips["2021"]
fig, axs = plt.subplots(1, 2, figsize=(13, 7.5))
for ax, col, vmax, lab in ((axs[0], "r21_mm", 250, "2021 lidar minus RTK (mm)"), (axs[1], "r08_mm", 600, "2008 lidar minus RTK (mm)")):
    ax.imshow(cir, extent=ext); ok = F[col].notna().values
    sc = ax.scatter(px[ok], py[ok], c=F[col].values[ok], cmap=cmc.bam_r, vmin=-vmax, vmax=vmax, s=22, edgecolors="k", linewidths=0.3)
    ax.scatter(px[~ok], py[~ok], s=14, facecolors="none", edgecolors="k", linewidths=0.4); ax.set_xticks([]); ax.set_yticks([]); plt.colorbar(sc, ax=ax, shrink=0.7, label=lab)
fig.savefig(f"{OUT}/floodplain_residual_map_naip2021.png", dpi=140, bbox_inches="tight"); plt.close(fig); print("wrote", f"{OUT}/floodplain_residual_map_naip2021.png")
