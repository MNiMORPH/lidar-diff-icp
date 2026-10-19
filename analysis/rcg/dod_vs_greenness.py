#!/usr/bin/env python3
"""Phase 4 of PLAN_RCG_REMOTE_SENSING.md: the DoD against an INDEPENDENT vegetation
measure for the gen1 epoch -- October-2008 Landsat greenness -- on the 30 m grid.

Like with like: the 5 m DoD, corrected DoD, gen1 span and penetration-failure flag are
AGGREGATED to each 30 m Landsat pixel (mean over the 36 cells; a pixel enters only if
> 50% of its cells are flat floodplain). Then, per pixel: DoD vs ndvi_oct08, vs
retained08, vs greenup08, and -- for contrast -- vs the 2024 Sentinel-2 October metric
aggregated the same way. Reported: Pearson r, slope (mm per unit), and DoD medians by
quartile of the greenness. x is not derived from the lidar, so the circularity of
section I of FLOODPLAIN_RETURN_PDF.md does not apply.

    ./lidar-icp/bin/python analysis/rcg/dod_vs_greenness.py
"""
import numpy as np
import pandas as pd
import rioxarray
from scipy.ndimage import label

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
g = rioxarray.open_rasterio(f"{R}/metric2008_ndvi_oct08.tif").squeeze(drop=True); xs, ys = g.x.values, g.y.values
M = {k: rioxarray.open_rasterio(f"{R}/metric2008_{k}.tif").squeeze(drop=True).values for k in ("ndvi_oct08", "retained08", "greenup08", "ndvi_sep08", "ndvi_sum08")}
S2 = {k: rioxarray.open_rasterio(f"{R}/metric_{k}.tif").squeeze(drop=True) for k in ("ndvi_oct", "retained", "greenup")}
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{S}/slope.npy"); sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); lod = np.load(f"{S}/lod.npy") * 1e3
ny, nx = dod.shape; flat = fp & (sl <= 2) & np.isfinite(dod)
def agg(arr5, fn=np.nanmean):
    out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - 15 - Y0) / res5); iy1 = iy0 + 6
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - 15 - X0) / res5); ix1 = ix0 + 6
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = fn(arr5[iy0:iy1, ix0:ix1])
    return out
flatfrac = agg(flat.astype(float)); A = dict(dod=agg(np.where(flat, dod, np.nan)), corr=agg(np.where(flat, new, np.nan)), span=agg(np.where(flat, sp, np.nan)),
                                           penfail=agg(np.where(flat, pf.astype(float), np.nan)), lod=agg(np.where(flat, lod, np.nan)))
# Sentinel-2 10 m -> 30 m by 3x3 block mean, on the Landsat grid
def s2_on_ls(da):
    v = da.values; sx, sy = da.x.values, da.y.values; out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        jy = np.flatnonzero(np.abs(sy - y) <= 15)
        for i, x in enumerate(xs):
            ix = np.flatnonzero(np.abs(sx - x) <= 15); out[j, i] = np.nanmean(v[np.ix_(jy, ix)]) if jy.size and ix.size else np.nan
    return out
for k in S2: M["s2_" + k] = s2_on_ls(S2[k])
sel = (flatfrac > 0.5) & np.isfinite(A["dod"]) & np.isfinite(M["ndvi_oct08"]) & np.isfinite(M["retained08"])
print(f"30 m pixels > 50% flat floodplain with 2008 metrics: {sel.sum():,}")
print("\nINDEPENDENT-x REGRESSIONS, per 30 m pixel (y = aggregated DoD, mm)")
print(f"  {'x':<18}{'y':<12}{'r':>7}{'slope':>14}{'DoD by x-quartile Q1..Q4 (mm)':>40}")
for xk in ("ndvi_oct08", "retained08", "greenup08", "ndvi_sep08", "s2_ndvi_oct", "s2_retained", "s2_greenup"):
    x = M[xk][sel]
    for yk in ("dod", "corr", "span", "penfail"):
        y = A[yk][sel]; okm = np.isfinite(x) & np.isfinite(y); x_, y_ = x[okm], y[okm]
        r = np.corrcoef(x_, y_)[0, 1]; s_ = np.polyfit(x_, y_, 1)[0]
        qe = np.percentile(x_, [25, 50, 75]); qs = [np.median(y_[(x_ >= lo) & (x_ < hi)]) for lo, hi in zip([-9] + list(qe), list(qe) + [9])]
        unit = "mm/unit" if yk != "penfail" else "frac/unit"
        if yk in ("dod", "corr") or xk in ("ndvi_oct08", "retained08"):
            print(f"  {xk:<18}{yk:<12}{r:>7.3f}{s_:>10.0f} {unit:<6}" + "".join(f"{v:>10.1f}" for v in qs))
print("\n  y: dod = raw star-only DoD; corr = vegetation-corrected DoD (section I); span = gen1 p90-p10 (mm); penfail = share of cells with no gen1 return at gen2")
print("  2008 x are the gen1 epoch's own autumn; the s2_ rows are 2024 and are the wrong epoch, shown for contrast")
pd.DataFrame({k: v[sel] for k, v in {**{"x_" + k: M[k] for k in M}, **{"y_" + k: A[k] for k in A}}.items()}).to_csv(f"{R}/dod_vs_greenness_30m.csv", index=False)
print(f"table: {R}/dod_vs_greenness_30m.csv")
