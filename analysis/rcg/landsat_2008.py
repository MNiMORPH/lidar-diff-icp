#!/usr/bin/env python3
"""Phase 2 of PLAN_RCG_REMOTE_SENSING.md: the 2008 state from Landsat 5/7, at 30 m, the
same metrics and the same KNOWN-first null check as Phase 1, then the question: were the
surviving erosion patches late-green in the autumn the gen1 lidar flew?

METRICS per 30 m pixel (MINE, stated):
    ndvi_spr08   mean of clear values among 2008-04-04, 05-06, 05-14, 05-22   spring green-up
    ndvi_sum08   max of 08-02, 08-10, 08-18                                   summer maximum
    ndvi_sep08   09-19 (L7)
    ndvi_oct08   10-29 (L5, gap-free, 0% cloud) -- THE scene; 10-21 (L7) as a check
    retained08   ndvi_oct08 / ndvi_sum08
    greenup08    mean(04-04, 05-06) / ndvi_sum08  (the two earliest dates)
Cross-sensor NDVI levels are NOT comparable to Sentinel-2; only within-scene ranks among
populations are read. Populations as Phase 1, read at the 30 m pixel centres; per-patch
rows need >= 5 Landsat pixels (0.45 ha) and report the count.

    ./lidar-icp/bin/python analysis/rcg/landsat_2008.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
import xarray as xr
from scipy.ndimage import label, uniform_filter

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
ld = lambda tag: rioxarray.open_rasterio(f"{R}/{tag}_ndvi.tif").squeeze(drop=True)
L = {k: ld(f"ls_2008{k}") for k in ("0404", "0506", "0514", "0522", "0802", "0810", "0818", "0919", "1021", "1029")}
grid = L["1029"]; xs, ys = grid.x.values, grid.y.values
spr = xr.concat([L[k] for k in ("0404", "0506", "0514", "0522")], "t").mean("t", skipna=True)
summ = xr.concat([L[k] for k in ("0802", "0810", "0818")], "t").max("t", skipna=True)
early = xr.concat([L["0404"], L["0506"]], "t").mean("t", skipna=True)
M = dict(ndvi_spr08=spr, ndvi_sum08=summ, ndvi_sep08=L["0919"], ndvi_oct08=L["1029"], ndvi_oct21_08=L["1021"], retained08=L["1029"] / summ, greenup08=early / summ)
for k, v in M.items():
    v.astype("float32").rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{R}/metric2008_{k}.tif", compress="deflate")
print("METRICS 2008 written:", ", ".join(f"metric2008_{k}.tif" for k in M))
print("DEFINITIONS: spr=mean(04-04,05-06,05-14,05-22); sum=max(08-02,08-10,08-18); sep=09-19; oct=10-29 (L5); oct21=10-21 (L7); retained=oct/sum; greenup=mean(04-04,05-06)/sum")

wc = rioxarray.open_rasterio(f"{R}/worldcover.tif").squeeze(drop=True)
def at_px(arr, xs_, ys_, x0, y0, r):
    ix = np.clip(((xs_ - x0) / r).astype(int), 0, arr.shape[1] - 1); iy = np.clip(((ys_ - y0) / r).astype(int), 0, arr.shape[0] - 1)
    out = arr[np.ix_(iy, ix)].astype(float); out[(ys_[:, None] < y0) | (ys_[:, None] >= y0 + arr.shape[0] * r) | (xs_[None, :] < x0) | (xs_[None, :] >= x0 + arr.shape[1] * r)] = np.nan
    return out
wcv = at_px(wc.values, xs, ys, float(wc.x[0]) - 5, float(wc.y[-1]) - 5, 10.0)   # WorldCover 10 m sampled at 30 m centres (y ascending index handled below)
# WorldCover raster y is descending; build it as an ascending-origin array for at_px
wca = wc.values[::-1]; wcv = at_px(wca, xs, ys, float(wc.x[0]) - 5, float(wc.y[-1]) - 5, 10.0)
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
# 30 m pixel = majority of the 36 5-m cells it covers, for the lidar-defined classes (so a pixel is "F2" only if mostly F2)
def frac_in(mask5):
    out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - 15 - Y0) / res5); iy1 = iy0 + 6
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - 15 - X0) / res5); ix1 = ix0 + 6
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = mask5[iy0:iy1, ix0:ix1].mean()
    return out
P = {k: frac_in(v) for k, v in dict(F1=F1, F2=F2, F3=F3, NE=NE, forest_hs=forest & (sl > 10), flat_fp=flat).items()}
pops = [("KNOWN: WorldCover tree cover (10)", wcv == 10), ("KNOWN: WorldCover cropland (40)", wcv == 40), ("KNOWN: WorldCover grassland (30)", wcv == 30),
        ("KNOWN: forest hillslope (>50% of px)", P["forest_hs"] > 0.5), ("KNOWN: F3 bare control (>50%)", P["F3"] > 0.5),
        ("all flat floodplain (>50%)", P["flat_fp"] > 0.5), ("F1 sedge cells (>50%)", P["F1"] > 0.5), ("F2 closed canopy (>50%)", P["F2"] > 0.5), ("NE patch (>50%)", P["NE"] > 0.5)]
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
keys = list(M)
rows = []
print("\n2008 POPULATION MEDIANS (NMAD) over 30 m Landsat pixels -- KNOWN rows first; levels not comparable to Sentinel-2, ranks are")
print(f"  {'population':<40}{'px':>6}" + "".join(f"{k:>14}" for k in keys))
for name, m in pops:
    m = m & np.isfinite(M["ndvi_sum08"].values) & np.isfinite(M["ndvi_oct08"].values)
    r = dict(population=name, pixels=int(m.sum())); line = f"  {name:<40}{m.sum():>6,}"
    for k in keys:
        v = M[k].values[m]; md = np.nanmedian(v); nm = nmad(v); r[k] = round(float(md), 3); r[k + "_nmad"] = round(float(nm), 3); line += f"{md:>9.2f}({nm:4.2f})"
    print(line); rows.append(r)
pd.DataFrame(rows).to_csv(f"{R}/phenology_populations_2008.csv", index=False)
labF2, nF2 = label(P["F2"] > 0.5, structure=np.ones((3, 3))); sz = np.bincount(labF2.ravel())
print(f"\nF2 patches with >= 5 Landsat pixels (0.45 ha): {(sz[1:] >= 5).sum()} -- 2008 per patch: retained / greenup / oct / sep / summer")
for k in np.argsort(sz[1:])[::-1][:10] + 1:
    if sz[k] < 5: break
    m = labF2 == k; jy, jx = np.nonzero(m)
    print(f"   {sz[k]:>3} px  E {xs[jx].mean():.0f} N {ys[jy].mean():.0f}  retained {np.nanmedian(M['retained08'].values[m]):.2f}  greenup {np.nanmedian(M['greenup08'].values[m]):.2f}  oct {np.nanmedian(M['ndvi_oct08'].values[m]):.2f}  sep {np.nanmedian(M['ndvi_sep08'].values[m]):.2f}  sum {np.nanmedian(M['ndvi_sum08'].values[m]):.2f}")
# rank test: where do F2 / NE sit in the floodplain's distribution of Oct-2008 retained greenness?
fpm = P["flat_fp"] > 0.5; r_fp = M["retained08"].values[fpm & np.isfinite(M["retained08"].values)]
for name, m in (("F2 closed canopy", P["F2"] > 0.5), ("NE patch", P["NE"] > 0.5), ("F1 sedge", P["F1"] > 0.5), ("F3 bare", P["F3"] > 0.5)):
    v = M["retained08"].values[m & np.isfinite(M["retained08"].values)]
    print(f"  percentile of {name:<18} median retained08 within the flat floodplain: {100*np.mean(r_fp < np.nanmedian(v)):.0f}th")
fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
order = ["0404", "0506", "0514", "0522", "0802", "0810", "0818", "0919", "1021", "1029"]; doys = [95, 127, 135, 143, 215, 223, 231, 263, 295, 303]
cols = ["#008300", "#eda100", "#1baf7a", "#52514e", "#2a78d6", "#9e9e9e", "#2a78d6", "#eb6834", "#e34948"]
for (name, m), col in zip(pops, cols):
    axes[0].plot(doys, [np.nanmedian(L[k].values[m]) for k in order], color=col, ls="--" if name.startswith("KNOWN") else "-", lw=2, marker="o", ms=3, label=f"{name} (n={m.sum():,})")
axes[0].set_xlabel("day of year 2008"); axes[0].set_ylabel("median NDVI (Landsat 5/7)"); axes[0].legend(fontsize=7, frameon=False); axes[0].set_title("Phenology by population, 2008; gen1 lidar flown November", fontsize=10)
for s_ in ("top", "right"): axes[0].spines[s_].set_visible(False)
im = axes[1].imshow(np.where(P["flat_fp"] > 0.2, M["retained08"].values, np.nan), extent=[xs[0] - 15, xs[-1] + 15, ys[-1] - 15, ys[0] + 15], cmap="YlGn", vmin=0.2, vmax=1.0, interpolation="nearest")
plt.colorbar(im, ax=axes[1], fraction=0.03, pad=0.02); yy, xx = np.mgrid[0:len(ys), 0:len(xs)]
axes[1].contour(xs[xx], ys[yy], (P["F2"] > 0.5).astype(float), levels=[0.5], colors="black", linewidths=0.8)
axes[1].contour(xs[xx], ys[yy], (P["F1"] > 0.5).astype(float), levels=[0.5], colors="#2a78d6", linewidths=0.8)
axes[1].set_title("Retained greenness 29 Oct 2008 / summer 2008 on the floodplain; black = closed-canopy, blue = sedge-corrected", fontsize=9)
axes[1].set_xticks([]); axes[1].set_yticks([]); axes[1].set_xlim(X0, X1); axes[1].set_ylim(Y0, Y1)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_phenology_2008.png", dpi=110); print(f"\nfigure: {FIG}/rcg_phenology_2008.png; table: {R}/phenology_populations_2008.csv")
