#!/usr/bin/env python3
"""Oak forest ground error structure (Battle Creek characterization site).

The top rung of the vegetation ladder (surface_roughness.py): a deciduous
leaf-off oak block. Unlike that script — which used gen1 last-return and so read
canopy, not ground, in tall vegetation — here gen1 is CSF-classified GROUND, for
a fair ground-vs-ground comparison, plus the ground-return DENSITY (penetration)
that drives forest error and the spatial structure of the ground roughness.

KEY RESULT: with proper CSF ground, gen1 and gen2 oak-forest ground roughness
differ by only ~1.6× (the sensor ratio), NOT the 8.7× seen with gen1 last-return
on prairie. The vegetation "divergence" is chiefly a ground-DEFINITION artifact
(last-return grabs canopy); CSF closes it. Both epochs penetrate leaf-off oak
comparably (~150 ground returns / 25 m^2 cell). Contrast Cook (conifer, leaf-on):
that forest breaks down because of poor penetration, a different regime.

The expensive step (CSF on the 24 M-point gen1 cloud + gen2 class-2 count) is
cached to ``data/derived/oak_forest.npz``; delete it to recompute. Roughness
rasters come from surface_roughness.py's cache.

Run:  env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/roughness_characterization/oak_forest.py
"""
import os
import numpy as np
import laspy
import pandas as pd
from scipy.stats import spearmanr

FRAME = (498135.0, 4975136.0, 499365.0, 4976876.0)
RES = 5.0
FOREST_BOX = (498150, 498850, 4975200, 4976050)   # oak block, NW-central
GEN1 = "data/before_battlecreek/gen1_4tile.laz"
GEN2 = "data/after_battlecreek/gen2_4tile.laz"
CACHE = "data/derived/oak_forest.npz"
CHAR = "data/derived/bc_roughness_char.npz"        # r2 (gen2 roughness), ND, from surface_roughness.py


def rough_count(path, class2, x0, y0, nx, ny, res=RES):
    """NMAD plane-residual roughness and point count per cell (ny x nx)."""
    f = laspy.read(path)
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    m = (np.asarray(f.classification) == 2) if class2 else \
        (np.asarray(f.return_number) == np.asarray(f.number_of_returns))
    x, y, z = x[m], y[m], z[m]
    inb = (x >= x0) & (x < x0 + nx * res) & (y >= y0) & (y < y0 + ny * res)
    x, y, z = x[inb], y[inb], z[inb]
    ix = np.floor((x - x0) / res).astype(np.int64); iy = np.floor((y - y0) / res).astype(np.int64)
    f2 = iy * nx + ix; N = nx * ny
    u = x - (x0 + (ix + 0.5) * res); v = y - (y0 + (iy + 0.5) * res)
    S = lambda w: np.bincount(f2, w, minlength=N)
    n = S(np.ones_like(z)); Su = S(u); Sv = S(v); Sz = S(z)
    Suu = S(u * u); Suv = S(u * v); Svv = S(v * v); Suz = S(u * z); Svz = S(v * z); Szz = S(z * z)
    M = np.empty((N, 3, 3))
    M[:, 0, 0] = Suu; M[:, 0, 1] = Suv; M[:, 0, 2] = Su
    M[:, 1, 0] = Suv; M[:, 1, 1] = Svv; M[:, 1, 2] = Sv
    M[:, 2, 0] = Su;  M[:, 2, 1] = Sv;  M[:, 2, 2] = n
    r = np.stack([Suz, Svz, Sz], axis=1)
    valid = (n >= 6) & (np.abs(np.linalg.det(M)) > 1e-6)
    beta = np.zeros((N, 3)); beta[valid] = np.linalg.solve(M[valid], r[valid])
    pred = beta[f2, 0] * u + beta[f2, 1] * v + beta[f2, 2]; resid = z - pred
    med = pd.Series(resid).groupby(f2).transform("median")
    nm = pd.Series(np.abs(resid - med.values)).groupby(f2).median()
    rough = np.full(N, np.nan); rough[nm.index.values] = 1.4826 * nm.values; rough[~valid] = np.nan
    return rough.reshape(ny, nx), n.reshape(ny, nx).astype(float)


def corr_length(field, mask, maxlag=8):
    """1/e spatial autocorrelation length (cells) of a masked anomaly field,
    averaged over row and column lags among in-mask pairs."""
    a = np.where(mask, field - np.nanmedian(field[mask]), np.nan)
    var = np.nanvar(a[mask])
    prof = [1.0]
    for lag in range(1, maxlag):
        vals = []
        for A, Bx in [(a[:, :-lag], a[:, lag:]), (a[:-lag, :], a[lag:, :])]:
            ok = np.isfinite(A) & np.isfinite(Bx)
            if ok.sum() > 100:
                vals.append(np.mean(A[ok] * Bx[ok]))
        prof.append(np.mean(vals) / var if vals else np.nan)
    return np.array(prof)


def main():
    x0, y0, x1, y1 = FRAME
    nx = int(round((x1 - x0) / RES)); ny = int(round((y1 - y0) / RES))
    d = np.load(CHAR); r2 = d["r2"]; ND = d["ND"]
    if os.path.exists(CACHE):
        c = np.load(CACHE); r1g = c["r1g"]; n1 = c["n1"]; n2 = c["n2"]
    else:
        print("CSF gen1 (whole frame) + gen2 count ... (~10 min)", flush=True)
        from lidar_diff_icp.ground import classify_ground_csf
        import shutil
        csf = classify_ground_csf(GEN1); r1g, n1 = rough_count(csf, True, x0, y0, nx, ny)
        shutil.rmtree(os.path.dirname(csf), ignore_errors=True)
        _, n2 = rough_count(GEN2, True, x0, y0, nx, ny)
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        np.savez(CACHE, r1g=r1g, n1=n1, n2=n2)

    ix = np.arange(nx); iy = np.arange(ny); EE, NN = np.meshgrid(x0 + (ix + 0.5) * RES, y0 + (iy + 0.5) * RES)
    b = FOREST_BOX
    forest = ((EE >= b[0]) & (EE <= b[1]) & (NN >= b[2]) & (NN <= b[3]) & (ND > 0.5)
              & np.isfinite(r1g) & np.isfinite(r2) & np.isfinite(n1) & np.isfinite(n2))
    print(f"oak forest cells = {forest.sum()}")
    print(f"GROUND PENETRATION (returns/cell): gen1(CSF) median={np.median(n1[forest]):.0f} "
          f"p10={np.percentile(n1[forest],10):.0f} | gen2 median={np.median(n2[forest]):.0f}")
    print(f"GROUND ROUGHNESS: gen1(CSF)={np.median(r1g[forest]):.3f}  gen2(class2)={np.median(r2[forest]):.3f} m "
          f"(ratio {np.median(r1g[forest])/np.median(r2[forest]):.2f} ~ sensor ratio)")
    print(f"roughness vs density (Spearman): gen1 r1g~n1={spearmanr(r1g[forest],n1[forest])[0]:+.2f} | "
          f"gen2 r2~n2={spearmanr(r2[forest],n2[forest])[0]:+.2f}")
    # spatial structure needs a CONTIGUOUS region: use the forest box (finite r2),
    # not the NDVI-scattered mask, so adjacent-cell pairs exist.
    boxm = ((EE >= b[0]) & (EE <= b[1]) & (NN >= b[2]) & (NN <= b[3])
            & (ND > 0.4) & np.isfinite(r2))
    prof = corr_length(r2, boxm)
    print(f"gen2 ground-roughness autocorr by lag (5 m cells), forest box n={boxm.sum()}: {np.round(prof,2)}")


if __name__ == "__main__":
    main()
