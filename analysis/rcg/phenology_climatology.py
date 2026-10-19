#!/usr/bin/env python3
"""Day-of-year phenology climatology of the floodplain classes from every clear Sentinel-2
scene 2016-2025, with the two lidar flight days placed on it.

Andy, 2026-10-19: "how do we handle the difference in time of year that we have the two
lidar surveys? Can we generate any day-of-year-based expectations?"

For each Sentinel-2 L2A scene over the tile (catalogue cloud < 40%) the per-class median
NDVI is computed from pixels the scene's own SCL mask calls clear (a class enters a date
only if >= 70% of its pixels are clear). Classes: the NE patch; the 2008 imagery class
(open floodplain, NE excluded); F1 sedge cells; F3 cropped cells; open floodplain in no
class; forest hillslope. Then a climatology: median and 25-75% band of NDVI in 15-day
day-of-year bins across all years, per class, with the year-to-year spread shown as the
band. Flight days: gen2 2021-05-01 (DOY 121; from the point gps_time, mn-3dep-audit);
gen1 2008-11-25 for Elba's own lines (DOY 330; analysis/LIDAR_DOCUMENTATION_MINE.md), the block 18-24 Nov 2008 per the NOAA metadata (DOY 323-329).

    ./lidar-icp/bin/python analysis/rcg/phenology_climatology.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import odc.stac
import pandas as pd
import planetary_computer as pc
import pyproj
import rioxarray
from pystac_client import Client
from scipy.ndimage import label

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res = 5.0; BUF = 100.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - BUF, Y0 - BUF); lo1, la1 = to_ll.transform(X1 + BUF, Y1 + BUF); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); forest = np.load(f"{S}/forest_pfs.npy").astype(bool); D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3; ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif")
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"]); open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
F1 = flat & (sp >= 400) & (dod < -lod) & (D > 0) & ~pf; F3 = flat & (sp < 140) & (np.abs(dod) < 30)
classes5 = {"NE patch": NE & flat, "2008 class (open, no NE)": flat & open5 & ~NE & (p08 > 0.5), "F1 sedge cells": F1 & open5, "F3 cropped cells": F3, "open floodplain, no class": flat & open5 & ~NE & ~(p08 > 0.5) & ~F1 & ~F3, "forest hillslope": forest & (sl > 10)}
def frac10(mask5, xs, ys):
    out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - 5 - Y0) / res); iy1 = iy0 + 2
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - 5 - X0) / res); ix1 = ix0 + 2
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = mask5[iy0:iy1, ix0:ix1].mean()
    return out
cache = f"{R}/phenology_climatology_scenes.csv"
if os.path.exists(cache):
    T = pd.read_csv(cache)
else:
    its = [i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime="2016-01-01/2025-12-31", query={"eo:cloud_cover": {"lt": 40}}).items() if "T15TWJ" in i.id]
    seen = set(); its = [i for i in sorted(its, key=lambda i: i.id) if not (i.properties["datetime"][:10] in seen or seen.add(i.properties["datetime"][:10]))]
    its.sort(key=lambda i: i.properties["datetime"]); print(f"{len(its)} candidate scenes 2016-2025 with catalogue cloud < 40%", flush=True)
    G = None; rows = []
    for k, it in enumerate(its):
        d = it.properties["datetime"][:10]
        try:
            ds = odc.stac.load([it], bands=["B04", "B08", "SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
        except Exception as e:
            print(f"  {d}: load failed {e}"); continue
        if G is None:
            xs, ys = ds.x.values, ds.y.values; G = {kk: frac10(v, xs, ys) > 0.5 for kk, v in classes5.items()}
        clear = (~ds["SCL"].squeeze(drop=True).isin([0, 1, 3, 8, 9, 10, 11])).values; r = (ds["B04"].squeeze(drop=True).astype("float32") / 1e4).values; n_ = (ds["B08"].squeeze(drop=True).astype("float32") / 1e4).values
        nd = np.where(clear & (n_ + r > 0), (n_ - r) / (n_ + r + 1e-9), np.nan)
        row = dict(date=d, year=int(d[:4]), doy=pd.Timestamp(d).dayofyear, cloud=it.properties["eo:cloud_cover"])
        for kk, m in G.items():
            cf = np.nanmean(clear[m]); row[kk] = np.nanmedian(nd[m]) if cf >= 0.7 else np.nan; row[kk + "_clear"] = cf
        rows.append(row)
        if k % 25 == 0: print(f"  {k}/{len(its)} {d}", flush=True)
    T = pd.DataFrame(rows); T.to_csv(cache, index=False)
print(f"scenes used: {len(T)}; per class, dates with >= 70% clear: " + ", ".join(f"{kk} {T[kk].notna().sum()}" for kk in classes5))
# climatology: 15-day bins
bins = np.arange(1, 367, 15); T["bin"] = np.digitize(T.doy, bins)
print("\nDOY CLIMATOLOGY (NDVI median across years in 15-day bins; 25-75% across scenes in parentheses)")
hdr = f"  {'DOY bin':<10}" + "".join(f"{kk[:16]:>22}" for kk in classes5); print(hdr)
clim = {}
for b in range(1, len(bins) + 1):
    s = T[T.bin == b]
    if len(s) == 0: continue
    line = f"  {bins[b-1]:>3}-{min(bins[b-1]+14, 366):<3}   "; clim[b] = {}
    for kk in classes5:
        v = s[kk].dropna(); clim[b][kk] = (np.median(v) if len(v) else np.nan, np.percentile(v, 25) if len(v) else np.nan, np.percentile(v, 75) if len(v) else np.nan, len(v))
        line += f"{np.median(v) if len(v) else np.nan:>8.2f} ({np.percentile(v,25) if len(v) else np.nan:4.2f}-{np.percentile(v,75) if len(v) else np.nan:4.2f}) n{len(v):<3}"
    print(line)
# figure
fig, ax = plt.subplots(figsize=(14, 6.5)); cols = {"NE patch": "#e34948", "2008 class (open, no NE)": "#eb6834", "F1 sedge cells": "#4a3aa7", "F3 cropped cells": "#eda100", "open floodplain, no class": "#9e9e9e", "forest hillslope": "#008300"}
for kk, col in cols.items():
    bx = [bins[b - 1] + 7 for b in clim]; med = [clim[b][kk][0] for b in clim]; lo = [clim[b][kk][1] for b in clim]; hi = [clim[b][kk][2] for b in clim]
    ax.fill_between(bx, lo, hi, color=col, alpha=0.12, lw=0); ax.plot(bx, med, color=col, lw=2, marker="o", ms=3, label=kk)
ax.axvline(121, color="black", lw=1.2, ls="--"); ax.text(123, 0.02, "gen2 flown 2021-05-01 (DOY 121)", fontsize=9, rotation=90, va="bottom")
ax.axvspan(323, 330, color="black", alpha=0.10); ax.axvline(330, color="black", lw=1.2, ls="--"); ax.text(332, 0.02, "gen1 flown 2008-11-25 (DOY 330; block 18-24 Nov)", fontsize=9, rotation=90, va="bottom")
ax.set_xlabel("day of year"); ax.set_ylabel("NDVI, class median per scene; line = median across 2016-2025, band = 25-75%"); ax.set_xlim(1, 366); ax.set_ylim(0, 0.95)
ax.set_title("Sentinel-2 phenology climatology of the floodplain classes, 2016-2025, with the two lidar flight days", fontsize=11); ax.legend(fontsize=8, frameon=False, loc="upper left")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_phenology_climatology.png", dpi=120); print(f"\nfigure: {FIG}/rcg_phenology_climatology.png; table: {cache}")
