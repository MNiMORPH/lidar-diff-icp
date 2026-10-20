#!/usr/bin/env python3
"""Sensitivity of the Clubb floodplain to the algorithm's parameters (Andy, 2026-10-20: "the
Clubb value is probably more realistic. Test against Clubb-algorithm parameters.").

Baseline = the adopted run corr1500_10m (10 m, 1500 m Whitewater corridor, tool defaults:
QQ_threshold 0.005, slope/relief percentiles 25-75, threshold_SO 3, surface_fitting_radius
30 m, minimum_patch_size 1000 px, threshold_contributing_pixels = 1 km2 at 10 m). Each
variant changes ONE parameter and reuses the baseline DEM. Reported per run: floodplain ha,
ha within 800 m of the Whitewater centrelines, Jaccard with the baseline, Jaccard with the
Elba pipeline mask, the QQ thresholds the tool chose, and the per-tile valley top the
Phase B rule derives (p95 of the DEM over the floodplain within 500 m of the river) for
three tiles. Writes data/derived/valley/clubb/sensitivity.csv.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import rasterio
import shapely
from shapely.geometry import shape

OUT = "data/derived/valley/clubb"; BASE = "corr1500_10m"
VARIANTS = {"qq0.002": ["QQ_threshold=0.002"], "qq0.01": ["QQ_threshold=0.01"], "pct10-90": ["slope_lower_percentile=10", "slope_upper_percentile=90", "relief_lower_percentile=10", "relief_upper_percentile=90"],
            "SO2": ["threshold_SO=2"], "SO4": ["threshold_SO=4"], "fit15": ["surface_fitting_radius=15"], "fit60": ["surface_fitting_radius=60"],
            "patch250": ["minimum_patch_size=250"], "patch4000": ["minimum_patch_size=4000"], "contrib0.25km2": ["threshold_contributing_pixels=2500"], "contrib4km2": ["threshold_contributing_pixels=40000"]}
riv = shapely.union_all([shape(f["geometry"]) for f in json.load(open("data/derived/valley/public_land/rivers_whitewater.geojson"))["features"]])
with rasterio.open(f"{OUT}/{BASE}/{BASE}.bil") as r: dem = r.read(1); dem = np.where(dem == r.nodata, np.nan, dem); tr = r.transform
ny, nx = dem.shape; jy, jx = np.mgrid[0:ny, 0:nx]; X, Y = tr * (jx + 0.5, jy + 0.5)
near800 = shapely.contains_xy(riv.buffer(800), X.ravel(), Y.ravel()).reshape(ny, nx); near500 = shapely.contains_xy(riv.buffer(500), X.ravel(), Y.ravel()).reshape(ny, nx)
m5 = np.load("data/derived/elbaext_delong_star_nosurf/floodplain_mask.npy").astype(bool); EX0, EY0 = 575600.0, 4882200.0
ey, ex = np.mgrid[0:m5.shape[0], 0:m5.shape[1]]; c, rr = ~tr * (EX0 + (ex.ravel() + 0.5) * 5, EY0 + (ey.ravel() + 0.5) * 5); ii = (np.clip(rr.astype(int), 0, ny - 1), np.clip(c.astype(int), 0, nx - 1))
tiles = {t: v["bounds"] for t, v in json.load(open("data/derived/valley/public_tiles.json")).items() if t in ("4342-29-64", "4342-28-64", "4358-26-03")}
def load_fp(name):
    with rasterio.open(f"{OUT}/{name}/{name}_valley.bil") as r: return r.read(1) == 1
def thresholds(name):
    t = {}
    for line in open(f"{OUT}/{name}/{name}.log"):
        if line.startswith("relief_threshold:"): parts = line.split(); t = {"relief_thr": float(parts[1]), "slope_thr": float(parts[3])}
    return t
base = load_fp(BASE); rows = []
def metrics(name, fp):
    row = dict(run=name, ha=fp.sum() * 0.01, ha_within_800m=np.sum(fp & near800) * 0.01, jaccard_base=np.sum(fp & base) / max(np.sum(fp | base), 1), elba_jaccard=np.sum(m5 & fp[ii].reshape(m5.shape)) / max(np.sum(m5 | fp[ii].reshape(m5.shape)), 1), **thresholds(name))
    for t, B in tiles.items():
        inside = (X >= B[0]) & (X < B[2]) & (Y >= B[1]) & (Y < B[3]); k = inside & near500 & fp & np.isfinite(dem)
        row[f"vtop_{t}"] = float(np.nanpercentile(dem[k], 95)) if k.sum() >= 50 else np.nan
    return row
rows.append(metrics(BASE, base))
for name, sets in VARIANTS.items():
    run = f"sens_{name}"
    if not os.path.exists(f"{OUT}/{run}/{run}_valley.bil"):
        cmd = ["nice", "-n", "10", "./lidar-icp/bin/python", "analysis/valley_floodplain_clubb.py", "--name", run, "--res", "10", "--corridor", "1500", "--dem-from", BASE] + sum([["--set", s_] for s_ in sets], [])
        print("RUN", run, sets, flush=True); r = subprocess.run(cmd, capture_output=True, text=True)
        if not os.path.exists(f"{OUT}/{run}/{run}_valley.bil"): print(f"  {run} FAILED (exit {r.returncode}); tail:", r.stdout[-300:], r.stderr[-300:]); rows.append(dict(run=run, failed=True)); continue
    rows.append(metrics(run, load_fp(run))); print(pd.DataFrame(rows[-1:]).round(3).to_string(index=False), flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/sensitivity.csv", index=False); pd.set_option("display.width", 250)
print("\nSENSITIVITY (one parameter changed per row from the baseline corr1500_10m):"); print(D.round(3).to_string(index=False))
