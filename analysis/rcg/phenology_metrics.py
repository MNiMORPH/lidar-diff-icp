#!/usr/bin/env python3
"""Phase 1 of PLAN_RCG_REMOTE_SENSING.md: present-day phenology metrics from Sentinel-2
2024 (replicated on 2023-10-23 for October), and the NULL CHECK on populations whose
state is known before the populations in question are read.

METRICS, per 10 m pixel, from the clear-masked NDVI rasters of fetch_scenes.py (MINE,
proposals; each is printed with its definition):
    ndvi_apr    mean of 2024-04-05 and 04-15          early green-up
    ndvi_may    mean of 05-05 and 05-10
    ndvi_sum    max of 08-03 and 08-08                 summer maximum
    ndvi_oct    mean of 10-02 and 10-17                mid-October greenness (the Wisconsin window)
    ndvi_lateoct  10-27 where clear, else 10-22        late October
    retained    ndvi_oct / ndvi_sum                    greenness retained into October
    greenup     ndvi_apr / ndvi_sum                    greenness reached by early April
    senesc_doy  first autumn date (09-27..11-01) with NDVI < 0.5 * ndvi_sum; 999 if never

POPULATIONS. Known state (the null check): ESA WorldCover 2021 tree cover (10), cropland
(40), grassland (30), herbaceous wetland (90), each restricted to the tile; forest
hillslopes (product forest mask, slope > 10 deg); F3 bare floodplain controls. In
question: F1 sedge cells (flat floodplain, pooled span >= 400 mm, read erosion, corrected),
F2 closed-canopy cells (flat, penetration failure), and the NE patch alone. Product layers
(5 m) are read at the Sentinel-2 pixel centres; population statistics are over Sentinel-2
PIXELS, not 5 m cells, so the counts are honest.

    ./lidar-icp/bin/python analysis/rcg/phenology_metrics.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import odc.stac
import pandas as pd
import planetary_computer as pc
import pyarrow.parquet as pq
import rioxarray
import xarray as xr
from pystac_client import Client
from scipy.ndimage import label, uniform_filter

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
ld = lambda tag: rioxarray.open_rasterio(f"{R}/{tag}_ndvi.tif").squeeze(drop=True)
dates = {"0405": "2024-04-05", "0415": "2024-04-15", "0505": "2024-05-05", "0510": "2024-05-10", "0803": "2024-08-03", "0808": "2024-08-08",
         "0927": "2024-09-27", "1002": "2024-10-02", "1017": "2024-10-17", "1022": "2024-10-22", "1027": "2024-10-27", "1101": "2024-11-01"}
N = {k: ld(f"s2_{v.replace('-', '')}") for k, v in dates.items()}; N["2023oct"] = ld("s2_20231023")
grid = N["1017"]; xs, ys = grid.x.values, grid.y.values
apr = (N["0405"] + N["0415"]) / 2; may = (N["0505"] + N["0510"]) / 2; summ = xr.concat([N["0803"], N["0808"]], "t").max("t")
octm = (N["1002"] + N["1017"]) / 2; lateoct = N["1027"].where(N["1027"].notnull(), N["1022"])
retained = octm / summ; greenup = apr / summ
autumn = [("0927", 271), ("1002", 276), ("1017", 291), ("1022", 296), ("1027", 301), ("1101", 306)]
sen = xr.full_like(summ, 999.0)
for k, doy in reversed(autumn):
    sen = sen.where(~(N[k] < 0.5 * summ), doy)
M = dict(ndvi_apr=apr, ndvi_may=may, ndvi_sum=summ, ndvi_oct=octm, ndvi_lateoct=lateoct, retained=retained, greenup=greenup, senesc_doy=sen, ndvi_oct2023=N["2023oct"])
for k, v in M.items():
    v.astype("float32").rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{R}/metric_{k}.tif", compress="deflate")
print("METRICS written:", ", ".join(f"metric_{k}.tif" for k in M))
print("DEFINITIONS: apr=mean(04-05,04-15); may=mean(05-05,05-10); sum=max(08-03,08-08); oct=mean(10-02,10-17); lateoct=10-27 else 10-22;")
print("  retained=oct/sum; greenup=apr/sum; senesc_doy=first of 09-27..11-01 with NDVI<0.5*sum (999 never); oct2023=2023-10-23")

# ---- ESA WorldCover at the same grid
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
import pyproj
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - 300, Y0 - 300); lo1, la1 = to_ll.transform(X1 + 300, Y1 + 300)
wc_items = list(cat.search(collections=["esa-worldcover"], bbox=[lo0, la0, lo1, la1]).items())
wc_items.sort(key=lambda i: i.properties.get("start_datetime", ""))
wc = odc.stac.load([wc_items[-1]], bands=["map"], bbox=[lo0, la0, lo1, la1], crs="EPSG:26915", resolution=10, chunks={})["map"].squeeze(drop=True)
wc = wc.reindex_like(grid, method="nearest"); wc.astype("uint8").rio.write_crs("EPSG:26915").rio.to_raster(f"{R}/worldcover.tif")
print(f"WorldCover item: {wc_items[-1].id}; classes on the tile: {dict(zip(*np.unique(wc.values, return_counts=True)))}")

# ---- product layers (5 m) at Sentinel-2 pixel centres
def at_px(arr5):
    ix = np.clip(((xs - X0) / res5).astype(int), 0, arr5.shape[1] - 1); iy = np.clip(((ys - Y0) / res5).astype(int), 0, arr5.shape[0] - 1)
    out = arr5[np.ix_(iy, ix)].astype(float); out[(ys[:, None] < Y0) | (ys[:, None] >= Y1) | (xs[None, :] < X0) | (xs[None, :] >= X1)] = np.nan
    return out
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); forest = np.load(f"{S}/forest_pfs.npy").astype(bool)
ny, nx = dod.shape
LO, HI, DZ = -5000.0, 5000.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H = H.reshape(ny, nx, NB)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9); C = np.cumsum(Hp, 2); tot = C[..., -1]
q = lambda p: ctr[np.argmax(C >= p * tot[..., None], axis=2)]; span_pool = np.where(tot > 0, q(0.9) - q(0.1), np.nan); del H, Hp, C
flat = fp & (sl <= 2) & np.isfinite(dod)
F1 = flat & (span_pool >= 400) & (dod < -lod) & (D > 0) & ~pf; F2 = flat & pf; F3 = flat & (sp < 140) & (np.abs(dod) < 30)
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
P = {k: at_px(v.astype(np.uint8)) == 1 for k, v in dict(F1=F1, F2=F2, F3=F3, NE=NE, forest_hs=forest & (sl > 10), flat_fp=flat).items()}
wcv = wc.values
pops = [("KNOWN: WorldCover tree cover (10)", wcv == 10), ("KNOWN: WorldCover cropland (40)", wcv == 40), ("KNOWN: WorldCover grassland (30)", wcv == 30),
        ("KNOWN: WorldCover herb. wetland (90)", wcv == 90), ("KNOWN: forest hillslope (product)", P["forest_hs"]), ("KNOWN: F3 bare floodplain control", P["F3"]),
        ("all flat floodplain", P["flat_fp"]), ("F1 sedge cells (corrected)", P["F1"]), ("F2 closed canopy (penfail)", P["F2"]), ("NE patch", P["NE"])]
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
keys = ["ndvi_apr", "ndvi_may", "ndvi_sum", "ndvi_oct", "ndvi_lateoct", "retained", "greenup", "senesc_doy", "ndvi_oct2023"]
rows = []
print("\nPOPULATION MEDIANS (NMAD) over Sentinel-2 pixels -- read the KNOWN rows first")
print(f"  {'population':<38}{'px':>7}" + "".join(f"{k:>13}" for k in keys))
for name, m in pops:
    m = m & np.isfinite(M["ndvi_sum"].values)
    vals = {k: M[k].values[m] for k in keys}; r = dict(population=name, pixels=int(m.sum()))
    line = f"  {name:<38}{m.sum():>7,}"
    for k in keys:
        v = vals[k]; md = np.nanmedian(v); nm = nmad(v); r[k] = round(float(md), 3); r[k + "_nmad"] = round(float(nm), 3)
        line += f"{md:>8.2f}({nm:4.2f})" if k != "senesc_doy" else f"{md:>8.0f}({nm:4.0f})"
    print(line); rows.append(r)
pd.DataFrame(rows).to_csv(f"{R}/phenology_populations_2024.csv", index=False)

# per-patch: F2 clusters >= 5 pixels
labF2, nF2 = label(P["F2"], structure=np.ones((3, 3))); sz = np.bincount(labF2.ravel())
print(f"\nF2 (closed canopy) patches with >= 5 Sentinel-2 pixels: {(sz[1:] >= 5).sum()} -- per patch median retained / greenup / senesc_doy")
for k in np.argsort(sz[1:])[::-1][:10] + 1:
    if sz[k] < 5: break
    m = labF2 == k; print(f"   {sz[k]:>4} px  E {xs[np.nonzero(m)[1]].mean():.0f} N {ys[np.nonzero(m)[0]].mean():.0f}  retained {np.nanmedian(M['retained'].values[m]):.2f}  greenup {np.nanmedian(M['greenup'].values[m]):.2f}  senesc {np.nanmedian(M['senesc_doy'].values[m]):.0f}  oct {np.nanmedian(M['ndvi_oct'].values[m]):.2f}")

# figure: class phenology curves + retained-greenness map over the floodplain
fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
order = ["0405", "0415", "0505", "0510", "0803", "0808", "0927", "1002", "1017", "1022", "1027", "1101"]
doys = [95, 105, 125, 130, 215, 220, 271, 276, 291, 296, 301, 306]
cols = ["#008300", "#eda100", "#1baf7a", "#4a3aa7", "#52514e", "#2a78d6", "#9e9e9e", "#2a78d6", "#eb6834", "#e34948"]
for (name, m), col in zip(pops, cols):
    ls = "--" if name.startswith("KNOWN") else "-"
    axes[0].plot(doys, [np.nanmedian(N[k].values[m]) for k in order], color=col, ls=ls, lw=2, marker="o", ms=3, label=f"{name} (n={m.sum():,})")
axes[0].set_xlabel("day of year 2024"); axes[0].set_ylabel("median NDVI"); axes[0].legend(fontsize=7, frameon=False); axes[0].set_title("Phenology by population, Sentinel-2 2024", fontsize=10)
for s_ in ("top", "right"): axes[0].spines[s_].set_visible(False)
im = axes[1].imshow(np.where(at_px(fp.astype(np.uint8)) == 1, M["retained"].values, np.nan), extent=[xs[0], xs[-1], ys[-1], ys[0]], cmap="YlGn", vmin=0.3, vmax=1.1, interpolation="nearest")
plt.colorbar(im, ax=axes[1], fraction=0.03, pad=0.02)
yy, xx = np.mgrid[0:len(ys), 0:len(xs)]
axes[1].contour(xs[xx], ys[yy], P["F2"].astype(float), levels=[0.5], colors="black", linewidths=0.6)
axes[1].contour(xs[xx], ys[yy], P["F1"].astype(float), levels=[0.5], colors="#2a78d6", linewidths=0.6)
axes[1].set_title("Retained greenness Oct/summer 2024 on the floodplain mask; black = closed-canopy cells, blue = sedge-corrected cells", fontsize=9)
axes[1].set_xticks([]); axes[1].set_yticks([]); axes[1].set_xlim(X0, X1); axes[1].set_ylim(Y0, Y1)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_phenology_2024.png", dpi=110); print(f"\nfigure: {FIG}/rcg_phenology_2024.png; table: {R}/phenology_populations_2024.csv")
