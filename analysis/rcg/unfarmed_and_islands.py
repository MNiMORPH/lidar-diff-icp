#!/usr/bin/env python3
"""Two tests Andy asked for (2026-10-20).

1. UNFARMED, BARE-IN-BOTH-EPOCHS SURFACES. Not a land-cover class but a time series: a
   cell is 'persistently bare' if its NDVI never exceeds 0.30 in any clear summer scene
   2016-2025 (one July/August Sentinel-2 scene per year) AND was below 0.30 in Landsat
   August 2008, and it is not water (WorldCover 80) and has no gen2 return above 2 m. Crops
   fail this every year; gravel, roads, lots, bare bars pass. The DoD on those cells is the
   two-epoch null on surfaces nobody farmed. (Thresholds MINE and stated.)

2. ISLANDS. Small 2008-class patches on open flat floodplain, fully ringed by non-class
   open floodplain. For each patch, gen1's own surface (z2 - DoD) is interpolated across
   the patch from a 2-cell ring (plane fit, least squares) and compared INSIDE the patch
   with the gen2 surface (the reference), beside the raw gen1 reading and the section-R
   percentile estimator (p12 of 3x3-pooled returns). Scored by patch size.

    ./lidar-icp/bin/python analysis/rcg/unfarmed_and_islands.py
"""
import os

import numpy as np
import odc.stac
import pandas as pd
import planetary_computer as pc
import pyarrow.parquet as pq
import pyproj
import rioxarray
import xarray as xr
from pystac_client import Client
from scipy.ndimage import binary_dilation, label, uniform_filter

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res = 5.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True); lo0, la0 = to_ll.transform(X0 - 100, Y0 - 100); lo1, la1 = to_ll.transform(X1 + 100, Y1 + 100); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
z2 = np.load(f"{S}/z_after.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"]); open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
wc = sample(f"{R}/worldcover.tif")

# ---------------- 1. persistently bare
print("1. PERSISTENTLY BARE: one clear July/August Sentinel-2 scene per year 2016-2025 (max NDVI per cell) + Landsat Aug 2008")
mx = np.full((ny, nx), -np.inf); used = []
for yr in range(2016, 2026):
    out = f"{R}/s2_summer_{yr}_ndvi.tif"
    if not os.path.exists(out):
        its = [i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{yr}-07-01/{yr}-08-31", query={"eo:cloud_cover": {"lt": 10}}).items() if "T15TWJ" in i.id]
        if not its: print(f"   {yr}: no clear summer scene"); continue
        it = sorted(its, key=lambda i: i.properties["eo:cloud_cover"])[0]
        ds = odc.stac.load([it], bands=["B04", "B08", "SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={}); clear = ~ds["SCL"].squeeze(drop=True).isin([0, 1, 3, 8, 9, 10, 11])
        r_ = ds["B04"].squeeze(drop=True).astype("float32") / 1e4; n_ = ds["B08"].squeeze(drop=True).astype("float32") / 1e4; ((n_ - r_) / (n_ + r_)).where(clear).rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out)
        used.append(it.properties["datetime"][:10])
    v = sample(out); mx = np.fmax(mx, np.where(np.isfinite(v), v, -np.inf))
print("   scenes:", used if used else "(cached)")
ls08 = np.fmax.reduce([sample(f"{R}/ls_2008{d}_ndvi.tif") for d in ("0802", "0810", "0818")])
bare = (mx < 0.30) & (ls08 < 0.30) & (wc != 80) & np.isfinite(dod) & (sl <= 5)
bare_open = bare & open5
print(f"   cells with max summer NDVI < 0.30 every year 2016-2025 AND Aug-2008 NDVI < 0.30, not water, slope <= 5: {bare.sum():,} (no gen2 return above 2 m: {bare_open.sum():,}); on the floodplain mask {int((bare & fp).sum()):,}")
for lab_, m in (("all persistently bare", bare), ("persistently bare, open", bare_open), ("persistently bare, on floodplain mask", bare & fp), ("same, flat (<= 2 deg)", bare & flat)):
    if m.sum() >= 10: print(f"   {lab_:<40} n {m.sum():>6,}  DoD median {np.median(dod[m]):>+6.1f}  NMAD {nmad(dod[m]):>5.1f}  p25/p75 {np.percentile(dod[m],25):>+5.0f}/{np.percentile(dod[m],75):>+5.0f}  WorldCover: {dict(zip(*np.unique(wc[m][np.isfinite(wc[m])].astype(int), return_counts=True)))}")
lb, nb = label(bare & flat, structure=np.ones((3, 3))); sb = np.bincount(lb.ravel())
print("   largest flat persistently-bare patches:")
for k in np.argsort(sb[1:])[::-1][:6] + 1:
    if sb[k] < 10: break
    m = lb == k; jy, jx = np.nonzero(m); print(f"     {sb[k]:>4} cells  E {E[m].mean():.0f} N {N[m].mean():.0f}  DoD med {np.median(dod[m]):+6.1f}  NMAD {nmad(dod[m]):5.1f}  WorldCover {int(pd.Series(wc[m]).mode().iloc[0]) if np.isfinite(wc[m]).any() else '?'}")

# ---------------- 2. islands
print("\n2. ISLANDS: 2008-class patches on open flat floodplain, ringed by non-class open floodplain; gen1 surface interpolated across each from its 2-cell ring")
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); dom = flat & open5 & ~NE; cls = dom & (p08 > 0.5); nonc = dom & ~(p08 > 0.5)
LO, HI, DZ = -1500.0, 1500.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool); cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); Hp = np.rint(uniform_filter(H.reshape(ny, nx, NB).astype(np.float32), size=(3, 3, 1), mode="nearest") * 9).reshape(-1, NB)
C = np.cumsum(Hp, 1); tot = C[:, -1]; p12 = np.where(tot > 0, ctr[np.argmax(C >= 0.12 * tot[:, None], axis=1)], np.nan).reshape(ny, nx)     # gen1 p12 of pooled returns, mm vs gen2
z1 = z2 - dod / 1e3                                                                                                                     # gen1 surface (m)
lc, ncl = label(cls, structure=np.ones((3, 3))); sz = np.bincount(lc.ravel())
rows = []
for k in range(1, ncl + 1):
    m = lc == k
    if sz[k] < 4 or sz[k] > 400: continue
    RING_IT = int(os.environ.get("RING_IT", "2")); RING_PURITY = float(os.environ.get("RING_PURITY", "0.8"))
    ring = binary_dilation(m, structure=np.ones((3, 3)), iterations=RING_IT) & ~m
    if (ring & ~nonc).sum() > (1 - RING_PURITY) * ring.sum() or (ring & nonc).sum() < 8: continue   # ring must be (>= RING_PURITY) non-class open floodplain
    rm = ring & nonc; A = np.column_stack([np.ones(rm.sum()), E[rm] - E[m].mean(), N[rm] - N[m].mean()]); coef, *_ = np.linalg.lstsq(A, z1[rm], rcond=None)
    pred = coef[0] + coef[1] * (E[m] - E[m].mean()) + coef[2] * (N[m] - N[m].mean())
    rows.append(dict(size=int(sz[k]), ring_n=int(rm.sum()), ring_dod=np.median(dod[rm]), interp_resid=np.median((pred - z2[m]) * 1e3), raw_resid=np.median(-dod[m]), p12_resid=np.nanmedian(p12[m]), penfail=pf[m].mean(), plane_rms=float(np.sqrt(np.mean((A @ coef - z1[rm]) ** 2)) * 1e3)))
F = pd.DataFrame(rows); print(f"   islands accepted: {len(F)} (size 4-400 cells, ring of {RING_IT} cells, >= {100*RING_PURITY:.0f}% non-class open floodplain)")
print(f"   {'size bin (cells)':<18}{'n':>4}{'ring DoD':>10}{'RAW gen1-gen2':>15}{'INTERP-gen2':>13}{'p12-gen2':>10}{'plane rms':>11}{'penfail':>9}   (medians over islands; NMAD in parentheses; mm)")
for lo, hi in ((4, 10), (10, 25), (25, 60), (60, 150), (150, 401)):
    s = F[(F["size"] >= lo) & (F["size"] < hi)]
    if len(s) == 0: continue
    print(f"   {f'{lo}-{hi-1}':<18}{len(s):>4}{np.median(s.ring_dod):>+10.0f}{np.median(s.raw_resid):>+10.0f} ({nmad(s.raw_resid):3.0f}){np.median(s.interp_resid):>+8.0f} ({nmad(s.interp_resid):3.0f}){np.median(s.p12_resid):>+6.0f} ({nmad(s.p12_resid):3.0f}){np.median(s.plane_rms):>11.0f}{100*np.median(s.penfail):>8.0f}%")
F.to_csv(f"{R}/islands_interpolation.csv", index=False); print(f"   table: {R}/islands_interpolation.csv")
