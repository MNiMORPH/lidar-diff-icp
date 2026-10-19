#!/usr/bin/env python3
"""Fetch FULL reflective band sets for the October scenes, the way the Wisconsin DNR method
needs them (Bernthal & Willis 2004; Hatch & Bernthal 2008): all six Landsat reflective bands
on one mid-October date, no NDVI; band 5 (SWIR 1.55-1.75 um) carried the separation.

Landsat-5 2008-10-29 (gap-free, 0% cloud): blue, green, red, nir08, swir16, swir22 (+qa).
Landsat-7 2008-10-21 and 2008-09-19 likewise. Sentinel-2 2024-10-17 and 2023-10-23: B02,
B03, B04, B05, B06, B07, B08, B8A, B11, B12 (+SCL) at 10 m (20 m bands resampled by odc).
Written as multi-band GeoTIFFs in EPSG:26915 with surface-reflectance scaling applied and
masked pixels NaN.

    ./lidar-icp/bin/python analysis/rcg/fetch_bands.py
"""
import numpy as np
import odc.stac
import planetary_computer as pc
import pyproj
import rioxarray  # noqa: F401
import xarray as xr
from pystac_client import Client

S = "data/derived/elbaext_delong_star_nosurf"; OUT = f"{S}/rcg"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; BUF = 300.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - BUF, Y0 - BUF); lo1, la1 = to_ll.transform(X1 + BUF, Y1 + BUF); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
LS = ["LT05_L2SP_026029_20081029_02_T1", "LE07_L2SP_026029_20081021_02_T1", "LE07_L2SP_026029_20080919_02_T1", "LT05_L2SP_026029_20080810_02_T1"]
LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]
for i in LS:
    it = list(cat.search(collections=["landsat-c2-l2"], ids=[i]).items())[0]
    ds = odc.stac.load([it], bands=LSB + ["qa_pixel"], bbox=BBOX, crs="EPSG:26915", resolution=30, chunks={})
    qa = ds["qa_pixel"].squeeze(drop=True).astype("uint16"); bad = 0
    for b in (0, 1, 3, 4, 5): bad = bad | ((qa >> b) & 1)
    clear = (bad == 0) & (ds["red"].squeeze(drop=True) > 0)
    arr = xr.concat([(ds[b].squeeze(drop=True).astype("float32") * 0.0000275 - 0.2).where(clear) for b in LSB], dim="band").assign_coords(band=LSB)
    arr.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{OUT}/ls_{it.properties['datetime'][:10].replace('-', '')}_bands.tif", compress="deflate")
    print(f"  {i}: bands {LSB}, clear {float(clear.mean()):.3f}")
S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
for d in ("2024-10-17", "2023-10-23", "2024-10-02", "2024-08-03"):
    its = sorted([i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{d}T00:00:00Z/{d}T23:59:59Z").items() if "T15TWJ" in i.id], key=lambda i: i.id)
    it = its[-1]
    ds = odc.stac.load([it], bands=S2B + ["SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
    scl = ds["SCL"].squeeze(drop=True); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11])
    arr = xr.concat([(ds[b].squeeze(drop=True).astype("float32") / 10000.0).where(clear) for b in S2B], dim="band").assign_coords(band=S2B)
    arr.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{OUT}/s2_{d.replace('-', '')}_bands.tif", compress="deflate")
    print(f"  {it.id}: bands {S2B}, clear {float(clear.mean()):.3f}")
print("done")
