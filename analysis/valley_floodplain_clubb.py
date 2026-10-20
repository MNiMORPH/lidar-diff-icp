#!/usr/bin/env python3
"""D1 (Andy, 2026-10-20): the valley floor from the OFFICIAL gen2 DEM by Clubb et al. 2017.

Mosaics the USGS 3DEP 1 m DEM tiles (project MN_SE_Driftless_2021_B21, data/after/usgs_1m_dem/)
over a box, resamples to --res (average), writes ENVI bil, and runs lsdtt-valley-metrics
(LSDTopoTools2 0.9, conda env `lsdtt`) with the floodplain switch. Every lsdtt parameter is
written to the driver file from the dict below so the run is reproducible from the repo;
the defaults are the tool's own (its README.csv) and Clubb 2017's (threshold_SO 3). Any value
Andy changes is a decision; any value I change is a proposal and must be said.

    ./lidar-icp/bin/python analysis/valley_floodplain_clubb.py --name elba_test \\
        --bounds 570000 4880000 580000 4890000 --res 5          # one tile, the tool-chain check
    ./lidar-icp/bin/python analysis/valley_floodplain_clubb.py --name valley5m --res 5   # the box

Outputs under data/derived/valley/clubb/<name>/: <name>.bil (+.hdr), <name>_valley.bil (the floodplain
mask as 0.9 names it; 1 = floodplain), <name>_filtered/_slope/_channel_relief.bil, <name>_qq_relief.txt,
<name>_qq_slope.txt, driver <name>.driver, log <name>.log, and floodplain_<name>.tif (the FP
mask as GeoTIFF for the rest of the pipeline).
"""
import argparse
import glob
import os
import subprocess
import sys
import time

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

LSDTT = os.path.expanduser("~/anaconda3/envs/lsdtt/bin/lsdtt-valley-metrics")
DEM_DIR = "data/after/usgs_1m_dem"; OUT = "data/derived/valley/clubb"
BOX = (552000.0, 4862000.0, 592000.0, 4900000.0)

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True)
ap.add_argument("--bounds", type=float, nargs=4, default=BOX, metavar=("X0", "Y0", "X1", "Y1"))
ap.add_argument("--res", type=float, default=5.0, help="grid spacing in m for the lsdtt run (PROPOSAL: 5)")
ap.add_argument("--threshold-SO", type=int, default=3, help="Clubb 2017 used 3 at every site")
ap.add_argument("--skip-lsdtt", action="store_true")
ap.add_argument("--corridor", type=float, default=None, help="if set, NoData outside this buffer (m) of the DNR Whitewater centrelines, so the quantile-quantile thresholds describe the VALLEY (hillslope vs floor), not the plateau. PROPOSAL (mine); Clubb 2017 assumes a dissected landscape")
A = ap.parse_args()
X0, Y0, X1, Y1 = A.bounds; D = f"{OUT}/{A.name}"; os.makedirs(D, exist_ok=True)

# ---- 1. mosaic the 1 m tiles that intersect the bounds, resample to res by averaging
tiles = []
for f in sorted(glob.glob(f"{DEM_DIR}/USGS_1M_15_*.tif")):
    with rasterio.open(f) as r:
        b = r.bounds
        if b.left < X1 and b.right > X0 and b.bottom < Y1 and b.top > Y0: tiles.append(f)
if not tiles: sys.exit("no DEM tiles on disk intersect the bounds")
print(f"DEM tiles intersecting the bounds: {len(tiles)}"); [print("  ", os.path.basename(t)) for t in tiles]
bil = f"{D}/{A.name}.bil"
if not os.path.exists(bil):
    # resample EACH 1 m tile into the destination grid (average), never a 1 m mosaic of the
    # box: the box at 1 m is 1.52e9 cells = 6 GB, the shared laptop's swap is full
    ny, nx = int(round((Y1 - Y0) / A.res)), int(round((X1 - X0) / A.res))
    dst = np.full((ny, nx), np.nan, "float32"); dtr = rasterio.transform.from_origin(X0, Y1, A.res, A.res); crs = None
    for t in tiles:
        with rasterio.open(t) as r:
            crs = r.crs; a = r.read(1).astype("float32"); a[a == r.nodata] = np.nan
            part = np.full((ny, nx), np.nan, "float32")
            reproject(a, part, src_transform=r.transform, src_crs=r.crs, dst_transform=dtr, dst_crs=r.crs, src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.average)
            dst = np.where(np.isfinite(part), part, dst); del a, part
        print(f"  resampled {os.path.basename(t)}; filled {np.isfinite(dst).mean():.4f} of the grid so far")
    if A.corridor:
        import json, shapely
        from shapely.geometry import shape
        riv = shapely.union_all([shape(f["geometry"]) for f in json.load(open("data/derived/valley/public_land/rivers_whitewater.geojson"))["features"]])
        jy, jx = np.mgrid[0:ny, 0:nx]; cx, cy = rasterio.transform.xy(dtr, jy.ravel(), jx.ravel())
        inside = shapely.contains_xy(riv.buffer(A.corridor), np.asarray(cx), np.asarray(cy)).reshape(ny, nx)
        dst[~inside] = np.nan; print(f"corridor: kept {inside.mean():.3f} of the grid within {A.corridor:.0f} m of the Whitewater centrelines")
    nod_out = -9999.0; out = np.where(np.isfinite(dst), dst, nod_out).astype("float32")
    with rasterio.open(bil, "w", driver="ENVI", height=ny, width=nx, count=1, dtype="float32", crs=crs, transform=dtr, nodata=nod_out) as w: w.write(out, 1)
    print(f"wrote {bil}: {ny} x {nx} at {A.res} m, nodata {np.mean(out == nod_out):.4f}, z {np.nanmin(dst):.1f}-{np.nanmax(dst):.1f} m")
    del dst, out

# ---- 2. the lsdtt driver: every parameter stated
P = {
    "read path": os.path.abspath(D), "write path": os.path.abspath(D), "read fname": A.name, "write fname": A.name,
    "dem read extension": "bil", "dem write extension": "bil",
    "min_slope_for_fill": "0.0001", "carve_before_fill": "true",   # roads cross the channels
    "surface_fitting_radius": "30",                 # tool default (m)
    "threshold_contributing_pixels": str(int(round(1e6 / A.res ** 2))),   # MINE: 1 km2 of drainage area, stated in pixels at this res
    "threshold_SO": str(A.threshold_SO),
    "floodplain_raster": "NULL", "get_terraces": "false",
    "slope_lower_percentile": "25", "slope_upper_percentile": "75", "relief_lower_percentile": "25", "relief_upper_percentile": "75",
    "QQ_threshold": "0.005", "use_absolute_thresholds": "false", "minimum_patch_size": "1000", "fill_floodplain": "true",
    "print_channels_to_csv": "true", "print_junctions_to_csv": "false",
    "write_hillshade": "false",
}
drv = f"{D}/{A.name}.driver"
with open(drv, "w") as f: f.write("# lsdtt-valley-metrics driver written by analysis/valley_floodplain_clubb.py\n" + "".join(f"{k}: {v}\n" for k, v in P.items()))
print("driver parameters:"); [print(f"  {k}: {v}") for k, v in P.items()]
if A.skip_lsdtt: sys.exit(0)

# ---- 3. run, timed, logged
t0 = time.time()
with open(f"{D}/{A.name}.log", "w") as lg:
    rc = subprocess.run([LSDTT, os.path.abspath(D) + "/", f"{A.name}.driver"], stdout=lg, stderr=subprocess.STDOUT).returncode
print(f"lsdtt-valley-metrics exit {rc} in {time.time() - t0:.0f} s; log {D}/{A.name}.log")
fp = f"{D}/{A.name}_valley.bil"   # LSDTopoTools2 0.9 names the floodplain mask _valley, not _FP (verified on the smoke test)
if os.path.exists(fp):
    with rasterio.open(fp) as r:
        m = r.read(1); prof = r.profile
    val = m[m != prof.get("nodata", -9999)] if prof.get("nodata") is not None else m
    print(f"floodplain raster {fp}: values {np.unique(m)[:6]}, floodplain share {np.mean(m == 1):.4f} = {np.sum(m == 1) * A.res ** 2 / 1e4:,.0f} ha")
    prof.update(driver="GTiff", compress="deflate")
    with rasterio.open(f"{D}/floodplain_{A.name}.tif", "w", **prof) as w: w.write(m, 1)
    for q in ("qq_relief", "qq_slope"):
        fq = f"{D}/{A.name}_{q}.txt"
        if os.path.exists(fq): print(f"  {q}: " + open(fq).read().strip().replace(chr(10), " | ")[:200])
else:
    print("NO _FP raster written; read the log"); print(open(f"{D}/{A.name}.log").read()[-3000:])
