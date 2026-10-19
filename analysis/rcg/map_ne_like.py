#!/usr/bin/env python3
"""Map the NE-patch spectral class in 2008 and 2025 and difference them.

Andy, 2026-10-19: "Can you identify the regions with reed canary grass at each time? ... Can
you then identify the difference, where it is now (and was not before), and where it was
before (but is not now)?"

WHAT IS MAPPED, honestly: pixels whose October reflectance is like the NE patch's IN THAT
YEAR -- in 2008 a bright, senesced, dense sward; in 2025 a green, high-NIR, late sward.
Both are states reed canary grass takes, and the Wisconsin reports list the confusions
(grazed pasture, cut hay, cranberry). Without field truth this is a spectral-class map,
not a species map, and it is labelled "NE-like".

METHOD (the Wisconsin shape, minus their field truth): two-class maximum likelihood on
all bands per scene. Class 1 = the NE patch's pixels (the training reference, by the
lidar's definition of the patch); class 0 = the rest of the flat floodplain. Gaussian
densities with full covariance, equal priors (MINE, stated); posterior p(NE-like); a
pixel is NE-like at p > 0.5. Domain for the map: the floodplain mask (the valley floor).
CONFUSION CHECK: the same classifier is run over every WorldCover class on the whole tile
and the share classified NE-like is reported -- if cropland or grassland come out mostly
NE-like, the class is not specific and the map says so.

CHANGE: Landsat-5 2008-10-29 against Landsat-9 2025-10-04 on the same 30 m grid:
  gained = NE-like 2025 and not 2008;  lost = NE-like 2008 and not 2025;  kept = both.
Sentinel-2 2025-10-02 gives the current map at 10 m. Public land (WMA) overlaid.

    ./lidar-icp/bin/python analysis/rcg/map_ne_like.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyproj
import rioxarray
import shapely
from matplotlib.colors import ListedColormap
from matplotlib.image import imread
from scipy.ndimage import label, uniform_filter
from shapely.geometry import shape

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
pf = np.load(f"{S}/penetration_failure.npy").astype(bool); ny, nx = dod.shape; flat = fp & (sl <= 2) & np.isfinite(dod)
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
wc = rioxarray.open_rasterio(f"{R}/worldcover.tif").squeeze(drop=True); wca = wc.values[::-1]; wcx0, wcy0 = float(wc.x[0]) - 5, float(wc.y[-1]) - 5
wma = shapely.union_all([shape(f["geometry"]) for f in json.load(open("data/derived/public_land/whitewater_wma_units.geojson"))["features"]])
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
WCN = {10: "tree", 20: "shrub", 30: "grassland", 40: "cropland", 50: "built", 60: "bare", 80: "water", 90: "herb. wetland"}


def frac_in(mask5, xs, ys, r):
    k = int(round(r / res5)); out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - r / 2 - Y0) / res5); iy1 = iy0 + k
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - r / 2 - X0) / res5); ix1 = ix0 + k
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = mask5[iy0:iy1, ix0:ix1].mean()
    return out


def at_px(arr, xs_, ys_, x0, y0, r):
    ix = np.clip(((xs_ - x0) / r).astype(int), 0, arr.shape[1] - 1); iy = np.clip(((ys_ - y0) / r).astype(int), 0, arr.shape[0] - 1)
    out = arr[np.ix_(iy, ix)].astype(float); out[(ys_[:, None] < y0) | (ys_[:, None] >= y0 + arr.shape[0] * r) | (xs_[None, :] < x0) | (xs_[None, :] >= x0 + arr.shape[1] * r)] = np.nan
    return out


def ml_map(tag, bn, r):
    da = rioxarray.open_rasterio(f"{R}/{tag}_bands.tif"); xs, ys = da.x.values, da.y.values; X = np.stack([da.values[i] for i in range(len(bn))], -1)
    fin = np.all(np.isfinite(X), -1); fpx = frac_in(fp, xs, ys, r) > 0.5; flx = frac_in(flat, xs, ys, r) > 0.5; nex = frac_in(NE, xs, ys, r) > 0.5
    wcv = at_px(wca, xs, ys, wcx0, wcy0, 10.0)
    c1 = X[nex & fin]; c0 = X[flx & ~nex & fin]
    def logpdf(Z, C):
        mu = C.mean(0); cov = np.cov(C.T) + np.eye(C.shape[1]) * 1e-6; inv = np.linalg.inv(cov); dz = Z - mu
        return -0.5 * np.einsum("ij,jk,ik->i", dz, inv, dz) - 0.5 * np.linalg.slogdet(cov)[1]
    Zf = X[fin]; l1 = logpdf(Zf, c1); l0 = logpdf(Zf, c0); p = np.full(fin.shape, np.nan); p[fin] = 1 / (1 + np.exp(np.clip(l0 - l1, -50, 50)))
    like = (p > 0.5)
    print(f"\n{tag}: trained on {len(c1):,} NE-patch pixels vs {len(c0):,} other flat-floodplain pixels, {len(bn)} bands")
    print(f"  resubstitution: NE patch classified NE-like {100*np.nanmean(like[nex & fin]):.0f}%; other flat floodplain {100*np.nanmean(like[flx & ~nex & fin]):.0f}%; whole floodplain mask {100*np.nanmean(like[fpx & fin]):.0f}%")
    print("  CONFUSION CHECK, share classified NE-like, by WorldCover class over the whole tile:")
    for k_, nm in WCN.items():
        m = (wcv == k_) & fin
        if m.sum() > 50: print(f"    {nm:<14} {m.sum():>8,} px   NE-like {100*np.nanmean(like[m]):5.1f}%")
    px_ha = r * r / 1e4
    out = like & fpx; print(f"  NE-like within the floodplain mask: {out.sum():,} px = {out.sum()*px_ha:.1f} ha of {fpx.sum()*px_ha:.1f} ha")
    da2 = da.isel(band=0).copy(data=p.astype("float32")); da2.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{R}/nelike_p_{tag}.tif", compress="deflate")
    return xs, ys, p, like, fpx, nex, wcv


LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]; S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
xs, ys, p08, like08, fpx, nex, wcv = ml_map("ls_20081029", LSB, 30.0)
xs2, ys2, p25, like25, fpx2, nex2, wcv2 = ml_map("ls_20251004", LSB, 30.0)
xs3, ys3, p25s, like25s, fpx3, nex3, wcv3 = ml_map("s2_20251002", S2B, 10.0)
assert np.allclose(xs, xs2) and np.allclose(ys, ys2)
dom = fpx & np.isfinite(p08) & np.isfinite(p25)
kept = dom & like08 & like25; gained = dom & ~like08 & like25; lost = dom & like08 & ~like25; neither = dom & ~like08 & ~like25
ha = 0.09
print(f"\nCHANGE on the floodplain mask, 30 m, Landsat 2008-10-29 -> 2025-10-04 (NE-like class each year):")
print(f"  kept   {kept.sum():>5,} px = {kept.sum()*ha:6.1f} ha\n  gained {gained.sum():>5,} px = {gained.sum()*ha:6.1f} ha  (NE-like now, not in 2008)\n  lost   {lost.sum():>5,} px = {lost.sum()*ha:6.1f} ha  (NE-like in 2008, not now)\n  neither{neither.sum():>5,} px = {neither.sum()*ha:6.1f} ha")
rows = []
for nm, m in (("gained", gained), ("lost", lost), ("kept", kept)):
    l_, n_ = label(m, structure=np.ones((3, 3))); sz = np.bincount(l_.ravel()); order = np.argsort(sz[1:])[::-1][:6] + 1
    print(f"  largest {nm} patches (>= 5 px = 0.45 ha):")
    for k in order:
        if sz[k] < 5: break
        mm = l_ == k; jy, jx = np.nonzero(mm); e, n0 = xs[jx].mean(), ys[jy].mean(); lon, lat = to_ll.transform(e, n0); pub = bool(wma.contains(shapely.Point(e, n0)))
        wcm = pd.Series(wcv[mm]).mode().iloc[0] if np.isfinite(wcv[mm]).any() else np.nan
        print(f"    {sz[k]:>4} px {sz[k]*ha:5.1f} ha  E {e:.0f} N {n0:.0f}  lat {lat:.5f} lon {lon:.5f}  WMA {pub}  WorldCover {WCN.get(int(wcm), '?') if np.isfinite(wcm) else '?'}  p08 {np.nanmedian(p08[mm]):.2f} p25 {np.nanmedian(p25[mm]):.2f}")
        rows.append(dict(change=nm, px=int(sz[k]), ha=round(sz[k] * ha, 2), easting=round(e, 1), northing=round(n0, 1), lat=round(lat, 6), lon=round(lon, 6), wma=pub, worldcover=WCN.get(int(wcm), "?") if np.isfinite(wcm) else "?"))
pd.DataFrame(rows).to_csv(f"{R}/nelike_change_patches.csv", index=False)
chg = np.full(dom.shape, np.nan); chg[neither] = 0; chg[kept] = 1; chg[gained] = 2; chg[lost] = 3
da = rioxarray.open_rasterio(f"{R}/ls_20081029_bands.tif").isel(band=0).copy(data=chg.astype("float32")); da.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{R}/nelike_change_2008_2025_30m.tif")

img = imread("data/derived/ne_patch_imagery/elbaext_fsa2023_2m.jpg")
fig, axes = plt.subplots(1, 3, figsize=(24, 8.5))
ext30 = [xs[0] - 15, xs[-1] + 15, ys[-1] - 15, ys[0] + 15]; ext10 = [xs3[0] - 5, xs3[-1] + 5, ys3[-1] - 5, ys3[0] + 5]
for ax, (ttl, arr, ext, cm, vmin, vmax) in zip(axes, (("p(NE-like), Landsat-5 29 Oct 2008", np.where(fpx, p08, np.nan), ext30, "magma", 0, 1), ("p(NE-like), Sentinel-2 2 Oct 2025, 10 m", np.where(fpx3, p25s, np.nan), ext10, "magma", 0, 1),
                                                     ("change 2008 -> 2025 (Landsat, 30 m)", chg, ext30, ListedColormap(["#d9d9d9", "#4a3aa7", "#1baf7a", "#e34948"]), -0.5, 3.5))):
    ax.imshow(img, extent=[X0, X1, Y0, Y1], origin="upper", alpha=0.55)
    im = ax.imshow(arr, extent=ext, origin="upper", cmap=cm, vmin=vmin, vmax=vmax, interpolation="nearest", alpha=0.9)
    for geom in (wma.geoms if hasattr(wma, "geoms") else [wma]):
        xg, yg = geom.exterior.xy; ax.plot(xg, yg, color="white", lw=1.0)
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_title(ttl, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    if ttl.startswith("change"): cb.set_ticks([0, 1, 2, 3]); cb.set_ticklabels(["neither", "kept", "gained (now, not 2008)", "lost (2008, not now)"])
fig.suptitle("'NE-like' = spectrally like the NE patch in that year's October scene (two-class maximum likelihood, all bands); floodplain mask; white = WMA", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_nelike_change.png", dpi=100); print(f"\nfigure: {FIG}/rcg_nelike_change.png; layers: {R}/nelike_p_*.tif, nelike_change_2008_2025_30m.tif; table: {R}/nelike_change_patches.csv")
