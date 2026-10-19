#!/usr/bin/env python3
"""Phase 4 redone with the FULL October-2008 band set: does any Landsat-5 2008-10-29 band,
SWIR index or multi-band combination predict the floodplain DoD at 30 m? (NDVI did not:
r -0.013. The Wisconsin method found the information in band 5 / all bands, not NDVI.)

y = 5 m DoD aggregated to each 30 m pixel (> 50% flat floodplain); x = reflectance bands,
brightness (sum of the six), NDMI, NBR, NDVI. Single-x: Pearson r, slope, DoD by
x-quartile. Multi-x: ordinary least squares on the six bands with 5-fold CV by 500 m block
(the same guard as the transfer function), reporting CV R^2 -- an INDEPENDENT-x model of
the gen1 floodplain bias, if it holds. Also the same for the vegetation-corrected DoD, for
what remains after the lidar's own correction.

    ./lidar-icp/bin/python analysis/rcg/dod_vs_bands.py
"""
import numpy as np
import pandas as pd
import rioxarray

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
da = rioxarray.open_rasterio(f"{R}/ls_20081029_bands.tif"); xs, ys = da.x.values, da.y.values
BN = ["blue", "green", "red", "nir08", "swir16", "swir22"]; B = {b: da.values[i] for i, b in enumerate(BN)}
B["bright"] = sum(B[b] for b in BN); B["NDVI"] = (B["nir08"] - B["red"]) / (B["nir08"] + B["red"]); B["NDMI"] = (B["nir08"] - B["swir16"]) / (B["nir08"] + B["swir16"]); B["NBR"] = (B["nir08"] - B["swir22"]) / (B["nir08"] + B["swir22"])
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{S}/slope.npy"); sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool)
ny, nx = dod.shape; flat = fp & (sl <= 2) & np.isfinite(dod)
def agg(arr5):
    out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - 15 - Y0) / res5); iy1 = iy0 + 6
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - 15 - X0) / res5); ix1 = ix0 + 6
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = np.nanmean(arr5[iy0:iy1, ix0:ix1])
    return out
ff = agg(flat.astype(float)); Y = dict(dod=agg(np.where(flat, dod, np.nan)), corr=agg(np.where(flat, new, np.nan)), span=agg(np.where(flat, sp, np.nan)), penfail=agg(np.where(flat, pf.astype(float), np.nan)))
sel = (ff > 0.5) & np.all([np.isfinite(B[b]) for b in BN], axis=0) & np.isfinite(Y["dod"])
print(f"30 m pixels > 50% flat floodplain, all bands finite: {sel.sum():,}")
jj, ii = np.nonzero(sel); block = (jj // 17) * 1000 + (ii // 17)       # ~500 m blocks on the 30 m grid
ub, binv = np.unique(block, return_inverse=True); fold = np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]
print(f"\nSINGLE-x, y = raw DoD (mm) and vegetation-corrected DoD; Landsat-5 2008-10-29 surface reflectance")
print(f"  {'x':<8}{'r(dod)':>8}{'slope':>12}{'DoD by x-quartile Q1..Q4':>36}{'r(corr)':>9}{'r(span)':>9}{'r(penfail)':>11}")
for k in BN + ["bright", "NDVI", "NDMI", "NBR"]:
    x = B[k][sel]; y = Y["dod"][sel]; r = np.corrcoef(x, y)[0, 1]; s_ = np.polyfit(x, y, 1)[0]
    qe = np.percentile(x, [25, 50, 75]); qs = [np.median(y[(x >= lo) & (x < hi)]) for lo, hi in zip([-9] + list(qe), list(qe) + [9])]
    print(f"  {k:<8}{r:>8.3f}{s_:>9.0f} mm" + "".join(f"{v:>9.1f}" for v in qs) + f"{np.corrcoef(x, Y['corr'][sel])[0,1]:>9.3f}{np.corrcoef(x, Y['span'][sel])[0,1]:>9.3f}{np.corrcoef(x, Y['penfail'][sel])[0,1]:>11.3f}")
print("\nMULTI-x OLS on the six bands, 5-fold CV by ~500 m block")
Xm = np.column_stack([B[b][sel] for b in BN]); Xm = np.column_stack([np.ones(Xm.shape[0]), Xm])
for yk in ("dod", "corr", "span"):
    y = Y[yk][sel]; pred = np.empty_like(y)
    for f in range(5):
        tr = fold != f; beta, *_ = np.linalg.lstsq(Xm[tr], y[tr], rcond=None); pred[~tr] = Xm[~tr] @ beta
    res = y - pred; beta, *_ = np.linalg.lstsq(Xm, y, rcond=None)
    print(f"  y = {yk:<6} CV R2 {1 - res.var() / y.var():.3f}   NMAD resid {1.4826*np.median(np.abs(res - np.median(res))):.1f}   in-sample R2 {1 - (y - Xm @ beta).var() / y.var():.3f}   coefficients (mm per unit reflectance): " + ", ".join(f"{b} {c:+.0f}" for b, c in zip(BN, beta[1:])))
pd.DataFrame({**{"x_" + k: B[k][sel] for k in BN + ["bright", "NDVI", "NDMI", "NBR"]}, **{"y_" + k: v[sel] for k, v in Y.items()}}).to_csv(f"{R}/dod_vs_bands_30m.csv", index=False)
print(f"table: {R}/dod_vs_bands_30m.csv")
