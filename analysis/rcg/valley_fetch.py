#!/usr/bin/env python3
"""Phase A, stage 1: FETCH everything the valley-wide imagery ranking needs, resumably.

Every output is written to a temporary name and renamed on completion, and every output
that already exists is skipped, so a dropped connection costs only the scene in flight.
Re-run until it prints "all present".

AOI: the Whitewater valley box E 552000-592000, N 4862000-4900000 (EPSG:26915), from
PLAN_VALLEY_FIELD_SITES.md section 1. Rasters are 10 m on that box; Landsat 30 m.

Fetched into data/derived/valley/:
  public_land/   WMA units, state-park management-unit POLYGONS (layer 2), WMA parking
                 (points) and roads (lines), DNR river centrelines with kittle_name
                 containing WHITEWATER -- Geospatial Commons feature services.
  s2_<date>_bands.tif   ten-band Sentinel-2 L2A mosaics of MGRS 15TWJ + 15TWK for the
                 PRESENCE dates (2026-10-17 current; 2025-10-27; 2021-10-18 the gen2 year),
                 SCL-masked. 
  s2_<date>_ndvi.tif    NDVI-only mosaics for the DENSITY dates (2026-04-10, 04-15, 05-15)
                 and the CUT-HISTORY dates (every scene 2026-06-15..10-20 with catalogue
                 cloud < 30% on both tiles).
  ls_20081029_bands.tif Landsat-5 2008-10-29 six bands over the box, QA-masked, with its
                 footprint coverage of the box printed (the plan marks this UNVERIFIED).
  scenes.json    the exact item ids used, per date and tile, with clear fractions.

    ./lidar-icp/bin/python analysis/rcg/valley_fetch.py
"""
import json
import os
import subprocess

import numpy as np
import odc.stac
import planetary_computer as pc
import pyproj
import rioxarray  # noqa: F401
import xarray as xr
from pystac_client import Client

V = "data/derived/valley"; PL = f"{V}/public_land"; os.makedirs(PL, exist_ok=True)
X0, Y0, X1, Y1 = 552000.0, 4862000.0, 592000.0, 4900000.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0, Y0); lo1, la1 = to_ll.transform(X1, Y1); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]; LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]
idx_path = f"{V}/scenes.json"; index = json.load(open(idx_path)) if os.path.exists(idx_path) else {}
def save_index(): json.dump(index, open(idx_path + ".tmp", "w"), indent=1); os.replace(idx_path + ".tmp", idx_path)
def atomic_raster(da, out, **kw):
    tmp = out + ".tmp.tif"; da.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(tmp, compress="deflate", **kw); os.replace(tmp, out)

# ---- 1. public land and hydrography (small; atomic via curl -o tmp + rename)
BB = f"{X0:.0f},{Y0:.0f},{X1:.0f},{Y1:.0f}"
Q = f"geometry={BB}&geometryType=esriGeometryEnvelope&inSR=26915&spatialRel=esriSpatialRelIntersects&outFields=*&outSR=26915&f=geojson"
BASE = "https://enterprise.gisdata.mn.gov/aghost/rest/services/us_mn_state_dnr"
for name, url in (("wma_units.geojson", f"{BASE}/bdry_dnr_wildlife_mgmt_areas_pub/FeatureServer/0/query?where=1%3D1&{Q}"),
                  ("state_parks.geojson", f"{BASE}/bdry_dnr_lrs_prk/FeatureServer/2/query?where=1%3D1&{Q}"),
                  ("wma_parking.geojson", f"{BASE}/struc_dnr_wld_mgmt_areas_pub/FeatureServer/0/query?where=1%3D1&{Q}"),
                  ("wma_roads.geojson", f"{BASE}/struc_dnr_wld_mgmt_areas_pub/FeatureServer/1/query?where=1%3D1&{Q}"),
                  ("rivers_whitewater.geojson", f"{BASE}/water_dnr_hydrography/FeatureServer/0/query?where=kittle_name%20LIKE%20%27%25WHITEWATER%25%27&{Q}")):
    out = f"{PL}/{name}"
    if os.path.exists(out): continue
    subprocess.run(["curl", "-s", "-m", "300", url, "-o", out + ".tmp"], check=True)
    try:
        d = json.load(open(out + ".tmp")); n = len(d.get("features", []))
        if "error" in d or n == 0: raise ValueError(d.get("error", "no features"))
    except Exception as e:
        os.remove(out + ".tmp"); print(f"  {name}: FAILED ({e})"); continue
    os.replace(out + ".tmp", out); print(f"  {name}: {n} features")

# ---- 2. Sentinel-2 helpers: both MGRS tiles for one date, mosaicked on the box grid
def s2_items(date):
    its = [i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{date}T00:00:00Z/{date}T23:59:59Z").items() if ("T15TWJ" in i.id or "T15TWK" in i.id)]
    best = {}
    for i in its:
        t = "T15TWJ" if "T15TWJ" in i.id else "T15TWK"
        if t not in best or i.id > best[t].id: best[t] = i
    return best
def s2_load(items, bands):
    ds = odc.stac.load(list(items.values()), bands=bands + ["SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, groupby="solar_day", chunks={})
    scl = ds["SCL"].isel(time=0); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11]) & (scl != 0)
    return ds, clear
def fetch_bands(date):
    out = f"{V}/s2_{date.replace('-', '')}_bands.tif"
    if os.path.exists(out): return
    items = s2_items(date)
    if not items: print(f"  {date}: no Sentinel-2 items"); return
    ds, clear = s2_load(items, S2B)
    arr = xr.concat([(ds[b].isel(time=0).astype("float32") / 10000.0).where(clear) for b in S2B], dim="band").assign_coords(band=S2B)
    atomic_raster(arr, out); cf = float(clear.mean())
    index[f"s2_{date}"] = {t: i.id for t, i in items.items()} | {"clear_frac": cf, "bands": S2B}; save_index(); print(f"  {date} bands: tiles {sorted(items)} clear {cf:.3f}", flush=True)
def fetch_ndvi(date):
    out = f"{V}/s2_{date.replace('-', '')}_ndvi.tif"
    if os.path.exists(out): return
    items = s2_items(date)
    if not items: print(f"  {date}: no Sentinel-2 items"); return
    ds, clear = s2_load(items, ["B04", "B08"])
    r = ds["B04"].isel(time=0).astype("float32") / 1e4; n_ = ds["B08"].isel(time=0).astype("float32") / 1e4
    atomic_raster(((n_ - r) / (n_ + r)).where(clear), out); cf = float(clear.mean())
    index[f"s2_{date}_ndvi"] = {t: i.id for t, i in items.items()} | {"clear_frac": cf}; save_index(); print(f"  {date} ndvi: tiles {sorted(items)} clear {cf:.3f}", flush=True)

print("PRESENCE scenes (ten bands):")
for d in ("2026-10-17", "2025-10-27", "2021-10-18"): fetch_bands(d)
print("DENSITY scenes (NDVI):")
for d in ("2026-04-10", "2026-04-15", "2026-05-15"): fetch_ndvi(d)
print("CUT-HISTORY scenes 2026-06-15..10-20, catalogue cloud < 30% on both tiles (NDVI):")
cands = {}
for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime="2026-06-15/2026-10-20", query={"eo:cloud_cover": {"lt": 30}}).items():
    if "T15TWJ" in i.id or "T15TWK" in i.id: cands.setdefault(i.properties["datetime"][:10], set()).add("T15TWJ" if "T15TWJ" in i.id else "T15TWK")
for d in sorted(d for d, ts in cands.items() if len(ts) == 2): fetch_ndvi(d)

# ---- 3. Landsat-5 2008-10-29 over the box, with coverage check
out = f"{V}/ls_20081029_bands.tif"
if not os.path.exists(out):
    it = list(cat.search(collections=["landsat-c2-l2"], ids=["LT05_L2SP_026029_20081029_02_T1"]).items())[0]
    ds = odc.stac.load([it], bands=LSB + ["qa_pixel"], bbox=BBOX, crs="EPSG:26915", resolution=30, chunks={})
    qa = ds["qa_pixel"].squeeze(drop=True).astype("uint16"); bad = 0
    for b in (0, 1, 3, 4, 5): bad = bad | ((qa >> b) & 1)
    clear = (bad == 0) & (ds["red"].squeeze(drop=True) > 0)
    arr = xr.concat([(ds[b].squeeze(drop=True).astype("float32") * 0.0000275 - 0.2).where(clear) for b in LSB], dim="band").assign_coords(band=LSB)
    atomic_raster(arr, out); cov = float((ds["red"].squeeze(drop=True) > 0).mean()); cf = float(clear.mean())
    index["ls_2008-10-29"] = {"id": it.id, "box_covered_frac": cov, "clear_frac": cf}; save_index(); print(f"  Landsat 2008-10-29: box covered {cov:.3f}, clear {cf:.3f}")
missing = [f for f in ([f"{V}/s2_{d.replace('-', '')}_bands.tif" for d in ("2026-10-17", "2025-10-27", "2021-10-18")] + [f"{V}/s2_{d.replace('-', '')}_ndvi.tif" for d in ("2026-04-10", "2026-04-15", "2026-05-15")] + [out]) if not os.path.exists(f)]
print("all present" if not missing else f"MISSING (re-run when online): {missing}")
