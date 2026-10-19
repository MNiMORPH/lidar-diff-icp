#!/usr/bin/env python3
"""Sentinel-2 timeline, 2016-2025: when did the NE field (and the cells flagged only in
2021) leave the bright senesced sward state?

Andy, 2026-10-19: the land was likely converted at one time; use the earliest Sentinel-2
imagery to check whether the sward was still in the indicated field then.

One late-October scene per year (clearest in Sep 25 - Nov 8, widened for years with none),
all ten reflective bands. Each scene is normalised to its own flat open floodplain (band /
floodplain median), which removes sensor/atmosphere/date level so years are comparable in
the way the 2008-vs-2025 comparison was. STATE per year, per group (NE patch; cells in the
class only in 2021; cells in the class only in 2008): the relative spectrum, and the share
of the group's pixels assigned to the BRIGHT state by ONE classifier -- two-class maximum
likelihood on relative bands trained on the 2021-10-18 scene (NE patch = bright state vs
rest of open floodplain), applied unchanged to every year. Also: the class's area on the
floodplain mask per year.

    ./lidar-icp/bin/python analysis/rcg/s2_timeline.py
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
import xarray as xr
from pystac_client import Client
from scipy.ndimage import label

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res = 5.0; BUF = 300.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - BUF, Y0 - BUF); lo1, la1 = to_ll.transform(X1 + BUF, Y1 + BUF); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
YEARS = {2016: "2016-09-29", 2017: "2017-10-09", 2018: None, 2019: "2019-10-29", 2020: "2020-10-13", 2021: "2021-10-18", 2022: "2022-10-18", 2023: "2023-10-23", 2024: "2024-10-17", 2025: "2025-10-27"}
if YEARS[2018] is None:
    its = [i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime="2018-09-15/2018-11-20", query={"eo:cloud_cover": {"lt": 40}}).items() if "T15TWJ" in i.id]
    its.sort(key=lambda i: (i.properties["eo:cloud_cover"], abs((pd.Timestamp(i.properties["datetime"][:10]) - pd.Timestamp("2018-10-20")).days)))
    YEARS[2018] = its[0].properties["datetime"][:10] if its else None; print(f"2018 widened search: {[(i.properties['datetime'][:10], round(i.properties['eo:cloud_cover'])) for i in its[:5]]} -> {YEARS[2018]}")
def fetch(d):
    out = f"{R}/s2_{d.replace('-', '')}_bands.tif"
    if os.path.exists(out): return out
    its = sorted([i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{d}T00:00:00Z/{d}T23:59:59Z").items() if "T15TWJ" in i.id], key=lambda i: i.id); it = its[-1]
    ds = odc.stac.load([it], bands=S2B + ["SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
    scl = ds["SCL"].squeeze(drop=True); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11])
    xr.concat([(ds[b].squeeze(drop=True).astype("float32") / 10000.0).where(clear) for b in S2B], dim="band").assign_coords(band=S2B).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out, compress="deflate")
    print(f"  fetched {it.id} clear {float(clear.mean()):.3f}"); return out
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"]); open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); p21 = sample(f"{R}/nelike_p_ls_20211017.tif")
dom = flat & open5 & ~NE; groups = {"NE patch": NE & flat, "2021-only cells": dom & ~(p08 > 0.5) & (p21 > 0.5), "2008-only cells": dom & (p08 > 0.5) & ~(p21 > 0.5), "both-years cells": dom & (p08 > 0.5) & (p21 > 0.5), "open floodplain, neither": dom & ~(p08 > 0.5) & ~(p21 > 0.5)}
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
# load all years as RELATIVE bands on the common 10 m grid
rel = {}; grids = None
for yr, d in YEARS.items():
    if d is None: continue
    da = rioxarray.open_rasterio(fetch(d)); xs, ys = da.x.values, da.y.values
    if grids is None:
        grids = (xs, ys); G = {k: frac10(v, xs, ys) > 0.5 for k, v in groups.items()}; fl = frac10(flat & open5, xs, ys) > 0.5; fpx = frac10(fp, xs, ys) > 0.5
    assert np.allclose(xs, grids[0]) and np.allclose(ys, grids[1])
    X = np.stack([da.values[i] for i in range(len(S2B))], -1); fin = np.all(np.isfinite(X), -1); ref = np.nanmedian(X[fl & fin], 0); rel[yr] = (X / ref, fin)
# one classifier, trained on 2021 relative bands: NE patch (bright state) vs rest of open floodplain
Xt, fint = rel[2021]; c1 = Xt[G["NE patch"] & fint]; c0 = Xt[fl & ~G["NE patch"] & fint]
def logpdf(Z, Cc):
    mu = Cc.mean(0); cov = np.cov(Cc.T) + np.eye(Cc.shape[1]) * 1e-6; inv = np.linalg.inv(cov); dz = Z - mu
    return -0.5 * np.einsum("ij,jk,ik->i", dz, inv, dz) - 0.5 * np.linalg.slogdet(cov)[1]
# a SECOND state: the dense GREEN canopy the field shows in 2016-2017 (relative NIR ~2.4, red ~0.4-0.8); classifier trained on 2016
Xg, fing = rel[2016]; g1 = Xg[G["NE patch"] & fing]; g0 = Xg[fl & ~G["NE patch"] & fing]
print("\nTWO-STATE TIMELINE: share of each group's pixels in the BRIGHT (2021-trained) and in the GREEN-DENSE (2016-trained) state; and each group's relative B04 / B08")
print(f"  {'year':<6}{'date':<12}" + "".join(f"{k[:14]:>31}" for k in G) )
print(f"  {'':<18}" + "".join(f"{'bright%':>8}{'green%':>7}{'B04':>7}{'B08':>7}  " for k in G))
for yr, (Xr, fin) in rel.items():
    Zf = Xr[fin]; pb = np.full(fin.shape, np.nan); pb[fin] = 1 / (1 + np.exp(np.clip(logpdf(Zf, c0) - logpdf(Zf, c1), -50, 50)))
    pg = np.full(fin.shape, np.nan); pg[fin] = 1 / (1 + np.exp(np.clip(logpdf(Zf, g0) - logpdf(Zf, g1), -50, 50)))
    line = f"  {yr:<6}{YEARS[yr]:<12}"
    for k, m in G.items():
        mm = m & fin; line += f"{100*np.nanmean(pb[mm] > 0.5):>8.0f}{100*np.nanmean(pg[mm] > 0.5):>7.0f}{np.nanmedian(Xr[..., 2][mm]):>7.2f}{np.nanmedian(Xr[..., 6][mm]):>7.2f}  "
    print(line)
rows = []
print(f"\nSTATE PER YEAR (Sentinel-2, late October; relative to each scene's open flat floodplain). Bright-state classifier trained on 2021.")
print(f"  {'year':<6}{'date':<12}{'NE: B04':>8}{'B08':>6}{'B11':>6}{'NDVI':>7}   {'share in BRIGHT state: NE':>26}{'2021-only':>11}{'2008-only':>11}{'both':>7}{'neither':>9}   {'class ha on floodplain':>22}")
for yr, (Xr, fin) in rel.items():
    Zf = Xr[fin]; p = np.full(fin.shape, np.nan); p[fin] = 1 / (1 + np.exp(np.clip(logpdf(Zf, c0) - logpdf(Zf, c1), -50, 50))); like = p > 0.5
    ne = G["NE patch"] & fin; b04, b08, b11 = (np.nanmedian(Xr[..., i][ne]) for i in (2, 6, 8)); nd = np.nanmedian(((Xr[..., 6] - Xr[..., 2]) / (Xr[..., 6] + Xr[..., 2]))[ne]) * 0 + np.nanmedian(((rioxarray.open_rasterio(fetch(YEARS[yr])).values[6] - rioxarray.open_rasterio(fetch(YEARS[yr])).values[2]) / (rioxarray.open_rasterio(fetch(YEARS[yr])).values[6] + rioxarray.open_rasterio(fetch(YEARS[yr])).values[2]))[ne])
    shares = {k: 100 * np.nanmean(like[m & fin]) for k, m in G.items()}; ha = like[fpx & fin].sum() * 0.01
    print(f"  {yr:<6}{YEARS[yr]:<12}{b04:>8.2f}{b08:>6.2f}{b11:>6.2f}{nd:>7.2f}   {shares['NE patch']:>25.0f}%{shares['2021-only cells']:>10.0f}%{shares['2008-only cells']:>10.0f}%{shares['both-years cells']:>6.0f}%{shares['open floodplain, neither']:>8.0f}%   {ha:>22.1f}")
    rows.append(dict(year=yr, date=YEARS[yr], ne_rel_B04=b04, ne_rel_B08=b08, ne_rel_B11=b11, ne_ndvi=nd, **{f"bright_{k}": v for k, v in shares.items()}, class_ha=ha))
pd.DataFrame(rows).to_csv(f"{R}/s2_timeline.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(17, 5.5)); wl = [490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190]; cm = plt.get_cmap("viridis")
good = [yr for yr in rel if yr != 2018]     # 2018-11-18 is 26% clear and reads as snow/cloud over the field; shown in the table, not the spectrum
for k, yr in enumerate(good):
    Xr, fin = rel[yr]; axes[0].plot(wl, np.nanmedian(Xr[G["NE patch"] & fin], 0), color=cm(k / (len(good) - 1)), marker="o", ms=3, lw=1.6, label=f"{yr} ({YEARS[yr][5:]})")
axes[0].set_ylim(0, 3)
axes[0].axhline(1, color="#52514e", lw=0.8); axes[0].set_xlabel("wavelength (nm)"); axes[0].set_ylabel("NE patch reflectance / floodplain median"); axes[0].set_title("NE patch relative spectrum, late October, by year", fontsize=10); axes[0].legend(fontsize=7, frameon=False, ncol=2)
df = pd.DataFrame(rows)
for k, col in (("NE patch", "#e34948"), ("2021-only cells", "#00bcd4"), ("2008-only cells", "#eda100"), ("both-years cells", "#4a3aa7"), ("open floodplain, neither", "#9e9e9e")):
    axes[1].plot(df.year, df[f"bright_{k}"], marker="o", lw=1.8, color=col, label=k)
axes[1].set_xlabel("year"); axes[1].set_ylabel("share of group's pixels in the BRIGHT senesced-sward state (%)"); axes[1].set_title("Bright-state membership by year (classifier trained on 2021)", fontsize=10); axes[1].legend(fontsize=8, frameon=False)
for ax in axes:
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_s2_timeline.png", dpi=110); print(f"\nfigure: {FIG}/rcg_s2_timeline.png; table: {R}/s2_timeline.csv")
