#!/usr/bin/env python3
"""Phase 0 of PLAN_RCG_REMOTE_SENSING.md: fetch the Sentinel-2 and Landsat scenes over the
Elba tile as NDVI GeoTIFFs in EPSG:26915, with their cloud/snow masks, and index them.

Scenes are named EXPLICITLY (the STAC query of 2026-10-19 listed them with cloud cover) so
the inputs of everything downstream are fixed and stated; summer scenes, which that query
did not cover, are searched here with a stated cloud cap and the chosen ones are written
to the index. Source: Microsoft Planetary Computer STAC (collections sentinel-2-l2a,
landsat-c2-l2). Landsat L2 surface reflectance is scaled 0.0000275 x DN - 0.2.

Masks. Sentinel-2: SCL classes 0,1,3,8,9,10,11 (nodata, saturated, shadow, cloud x3, snow)
are dropped -- the audit repo's rule. Landsat: QA_PIXEL bits 0 fill, 1 dilated cloud,
3 cloud, 4 shadow, 5 snow are dropped; the snow bit is REPORTED separately because the
2008-10-29 scene's value hinges on it. For every scene the clear fraction over the tile
and over the floodplain mask is printed; a scene is not silently used if it is mostly
masked.

    ./lidar-icp/bin/python analysis/rcg/fetch_scenes.py
"""
import json
import os

import numpy as np
import odc.stac
import planetary_computer as pc
import rioxarray  # noqa: F401  (registers .rio)
from pystac_client import Client

S = "data/derived/elbaext_delong_star_nosurf"; OUT = f"{S}/rcg"; os.makedirs(OUT, exist_ok=True)
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0          # the tile (corrections.json bounds)
BUF = 300.0
import pyproj
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lon0, lat0 = to_ll.transform(X0 - BUF, Y0 - BUF); lon1, lat1 = to_ll.transform(X1 + BUF, Y1 + BUF)
BBOX = [lon0, lat0, lon1, lat1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)

S2_DATES = ["2024-04-05", "2024-04-15", "2024-05-05", "2024-05-10", "2024-09-27", "2024-10-02", "2024-10-17",
            "2024-10-22", "2024-10-27", "2024-11-01", "2023-10-23"]          # tile T15TWJ; ids carry a processing suffix, so select by date
LS_IDS = ["LE07_L2SP_026029_20080919_02_T1", "LE07_L2SP_026029_20081021_02_T1", "LT05_L2SP_026029_20081029_02_T1",
          "LT05_L2SP_026029_20080404_02_T1", "LT05_L2SP_026029_20080506_02_T1", "LE07_L2SP_026029_20080514_02_T1",
          "LT05_L2SP_026029_20080522_02_T1"]
SUMMER = [("sentinel-2-l2a", "2024-06-25/2024-08-25", 10.0, 3), ("landsat-c2-l2", "2008-06-20/2008-08-31", 15.0, 3)]
print(f"AOI: tile {X0:.0f}-{X1:.0f} E, {Y0:.0f}-{Y1:.0f} N (+{BUF:.0f} m); bbox lon/lat {np.round(BBOX, 4).tolist()}")
fpm = np.load(f"{S}/floodplain_mask.npy").astype(bool)


def by_id(coll, ids):
    out = []
    for i in ids:
        it = list(cat.search(collections=[coll], ids=[i]).items())
        if not it: print(f"  NOT FOUND: {i}"); continue
        out.append(it[0])
    return out


def s2_by_date(dates):
    out = []
    for d in dates:
        its = [i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{d}T00:00:00Z/{d}T23:59:59Z").items()
               if "T15TWJ" in i.id]
        if not its: print(f"  NOT FOUND: Sentinel-2 {d} T15TWJ"); continue
        its.sort(key=lambda i: i.id); out.append(its[-1])          # latest processing baseline
    return out


def summer(coll, dt, maxc, n):
    its = list(cat.search(collections=[coll], bbox=BBOX, datetime=dt, query={"eo:cloud_cover": {"lt": maxc}}).items())
    its = [i for i in its if coll != "landsat-c2-l2" or i.id.split("_")[2] == "026029"]
    its.sort(key=lambda i: i.properties.get("eo:cloud_cover", 100))
    print(f"  summer {coll} {dt} cloud<{maxc:g}: {len(its)} candidates; taking {[i.id for i in its[:n]]}")
    return its[:n]


def floodplain_fraction(da_mask):
    """share of floodplain-mask area that is CLEAR in this scene (sampled at 5 m cell centres)."""
    yy, xx = np.nonzero(fpm); E = X0 + (xx + 0.5) * 5.0; N = Y0 + (yy + 0.5) * 5.0
    import xarray as xr
    v = da_mask.sel(x=xr.DataArray(E, dims="pts"), y=xr.DataArray(N, dims="pts"), method="nearest").values
    return float(np.nanmean(v))


index = []
for coll, items in (("sentinel-2-l2a", s2_by_date(S2_DATES) + summer(*SUMMER[0])),
                    ("landsat-c2-l2", by_id("landsat-c2-l2", LS_IDS) + summer(*SUMMER[1]))):
    for it in items:
        p = it.properties; date = p["datetime"][:10]; plat = p.get("platform", "")
        if coll == "sentinel-2-l2a":
            ds = odc.stac.load([it], bands=["B04", "B08", "SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
            red = ds["B04"].astype("float32") / 10000.0; nir = ds["B08"].astype("float32") / 10000.0
            scl = ds["SCL"].squeeze(); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11]); snow = scl == 11; res = 10
        else:
            ds = odc.stac.load([it], bands=["red", "nir08", "qa_pixel"], bbox=BBOX, crs="EPSG:26915", resolution=30, chunks={})
            red = ds["red"].astype("float32") * 0.0000275 - 0.2; nir = ds["nir08"].astype("float32") * 0.0000275 - 0.2
            qa = ds["qa_pixel"].squeeze().astype("uint16"); bad = 0
            for b in (0, 1, 3, 4, 5): bad = bad | ((qa >> b) & 1)
            clear = (bad == 0) & (ds["red"].squeeze() > 0); snow = ((qa >> 5) & 1) == 1; res = 30
        ndvi = ((nir - red) / (nir + red)).squeeze().where(clear).astype("float32")
        ndvi = ndvi.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan)
        tag = f"{'s2' if coll.startswith('sent') else 'ls'}_{date.replace('-', '')}"
        ndvi.rio.to_raster(f"{OUT}/{tag}_ndvi.tif", compress="deflate")
        clear.astype("uint8").rio.write_crs("EPSG:26915").rio.to_raster(f"{OUT}/{tag}_clear.tif", compress="deflate")
        cf = float(clear.mean().values); sf = float(snow.mean().values); fpf = floodplain_fraction(clear.astype("float32"))
        rec = dict(id=it.id, collection=coll, date=date, platform=plat, cloud_catalogue_pct=p.get("eo:cloud_cover"),
                   res_m=res, clear_frac_aoi=round(cf, 3), clear_frac_floodplain=round(fpf, 3), snow_frac_aoi=round(sf, 4), file=f"{tag}_ndvi.tif")
        index.append(rec)
        print(f"  {date} {plat:<12} cloud(cat) {p.get('eo:cloud_cover', float('nan')):5.1f}%  clear AOI {cf:.3f}  clear floodplain {fpf:.3f}  snow {sf:.4f}  -> {tag}_ndvi.tif", flush=True)
json.dump(index, open(f"{OUT}/scenes.json", "w"), indent=1)
print(f"\nwrote {len(index)} scenes to {OUT}/ and {OUT}/scenes.json")
