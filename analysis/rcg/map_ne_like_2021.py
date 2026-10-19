#!/usr/bin/env python3
"""NE-like map for the gen2 year: Sentinel-2 2021-10-18 (0.1% cloud) and Landsat-8
2021-10-17 (0%), the October before/after gen2's May 2021 flight. Same classifier as
map_ne_like.py (two-class maximum likelihood, all bands, NE patch vs rest of flat
floodplain, equal priors, p > 0.5). Writes rcg/nelike_p_{s2_20211018,ls_20211017}.tif and
the band files.

    ./lidar-icp/bin/python analysis/rcg/map_ne_like_2021.py
"""
import os

import numpy as np
import odc.stac
import planetary_computer as pc
import pyproj
import rioxarray
import xarray as xr
from pystac_client import Client
from scipy.ndimage import label

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0; BUF = 300.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - BUF, Y0 - BUF); lo1, la1 = to_ll.transform(X1 + BUF, Y1 + BUF); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]; S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
out = f"{R}/ls_20211017_bands.tif"
if not os.path.exists(out):
    it = list(cat.search(collections=["landsat-c2-l2"], ids=["LC08_L2SP_026029_20211017_02_T1"]).items())[0]
    ds = odc.stac.load([it], bands=LSB + ["qa_pixel"], bbox=BBOX, crs="EPSG:26915", resolution=30, chunks={})
    qa = ds["qa_pixel"].squeeze(drop=True).astype("uint16"); bad = 0
    for b in (0, 1, 3, 4, 5): bad = bad | ((qa >> b) & 1)
    clear = (bad == 0) & (ds["red"].squeeze(drop=True) > 0)
    xr.concat([(ds[b].squeeze(drop=True).astype("float32") * 0.0000275 - 0.2).where(clear) for b in LSB], dim="band").assign_coords(band=LSB).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out, compress="deflate")
    print(f"fetched Landsat-8 2021-10-17 clear {float(clear.mean()):.3f}")
out = f"{R}/s2_20211018_bands.tif"
if not os.path.exists(out):
    its = sorted([i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime="2021-10-18T00:00:00Z/2021-10-18T23:59:59Z").items() if "T15TWJ" in i.id], key=lambda i: i.id); it = its[-1]
    ds = odc.stac.load([it], bands=S2B + ["SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
    scl = ds["SCL"].squeeze(drop=True); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11])
    xr.concat([(ds[b].squeeze(drop=True).astype("float32") / 10000.0).where(clear) for b in S2B], dim="band").assign_coords(band=S2B).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out, compress="deflate")
    print(f"fetched {it.id} clear {float(clear.mean()):.3f}")

dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod)
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
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
for tag, bn, r in (("ls_20211017", LSB, 30.0), ("s2_20211018", S2B, 10.0)):
    da = rioxarray.open_rasterio(f"{R}/{tag}_bands.tif"); xs, ys = da.x.values, da.y.values; X = np.stack([da.values[i] for i in range(len(bn))], -1)
    fin = np.all(np.isfinite(X), -1); flx = frac_in(flat, xs, ys, r) > 0.5; nex = frac_in(NE, xs, ys, r) > 0.5; fpx = frac_in(fp, xs, ys, r) > 0.5
    c1 = X[nex & fin]; c0 = X[flx & ~nex & fin]
    def logpdf(Z, C):
        mu = C.mean(0); cov = np.cov(C.T) + np.eye(C.shape[1]) * 1e-6; inv = np.linalg.inv(cov); dz = Z - mu
        return -0.5 * np.einsum("ij,jk,ik->i", dz, inv, dz) - 0.5 * np.linalg.slogdet(cov)[1]
    Zf = X[fin]; p = np.full(fin.shape, np.nan); p[fin] = 1 / (1 + np.exp(np.clip(logpdf(Zf, c0) - logpdf(Zf, c1), -50, 50)))
    like = p > 0.5
    print(f"{tag}: NE patch resubstitution {100*np.nanmean(like[nex & fin]):.0f}%; other flat floodplain {100*np.nanmean(like[flx & ~nex & fin]):.0f}%; NE-like on floodplain mask {like[fpx & fin].sum() * r * r / 1e4:.1f} ha")
    # the patch's 2021 state, relative to the floodplain, for the record
    ref = np.nanmedian(X[flx & fin], 0); print("  NE patch / floodplain, " + ", ".join(f"{b} {v:.2f}" for b, v in zip(bn, np.nanmedian(X[nex & fin], 0) / ref)))
    da.isel(band=0).copy(data=p.astype("float32")).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(f"{R}/nelike_p_{tag}.tif", compress="deflate")
