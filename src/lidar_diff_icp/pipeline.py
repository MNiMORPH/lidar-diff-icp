"""End-to-end differencing of 2008-era Minnesota lidar against USGS 3DEP.

One function, ``difference_dem``, runs the whole validated workflow and returns a
DEM of Difference plus its corrections and uncertainty. It encodes the lessons
from the Elba pilot; see the module notes below and the project README.

Lessons baked in
----------------
* **Bare earth = last return, ``return_number == number_of_returns`` INCLUDING
  single returns.** Dropping singles empties flat/open ground (they dominate it).
* **Ground = a LOW PERCENTILE (10th) per cell, never mean/median.** On rough,
  vegetated, or sloping cells the true ground sits at the bottom of the return
  distribution; any central-tendency estimate rides above it and, because the
  two epochs differ greatly in density, that offset becomes COHERENT false change
  (~16-32% of convex hillslopes falsely depositional). A low percentile tracks
  the ground and removes it (~4%). This is the single most important choice.
* **Correct in the acquisition-honest order, per point, BEFORE gridding:**
  per-swath internal alignment (translation) -> spatially varying quadratic tie
  -> DeLong 400 m correction surface on flats (TPI floodplain buffer) -> per-swath
  along-track GNSS-drift spline ``f(gps_time)``. The drift is the deterministic,
  reusable core: the same early-lidar failure statewide, only coefficients differ.
* **TPI, not flow accumulation, buffers the floodplain** out of the stable set
  (flow routing is unreliable on flats).
* Convention: DoD is always ``after - before`` (positive = deposition); plot red =
  erosion, blue = deposition; standard NW (315/45) hillshade.
"""
from __future__ import annotations

import numpy as np
import laspy
import pandas as pd
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt as edt

from . import io, coreg


def read_last_return(path, bounds=None):
    """Read a LAZ and return last-return points (rn == nr, singles included).

    Returns dict of arrays: x, y, z, point_source_id, gps_time (gps_time is zeros
    if the file lacks it). ``bounds`` = (minx, miny, maxx, maxy) clips if given.
    """
    f = laspy.read(str(path))
    rn = np.asarray(f.return_number); nr = np.asarray(f.number_of_returns)
    m = rn == nr
    x = np.asarray(f.x)[m]; y = np.asarray(f.y)[m]; z = np.asarray(f.z)[m]
    ps = np.asarray(f.point_source_id)[m]
    try:
        gt = np.asarray(f.gps_time)[m]
    except Exception:
        gt = np.zeros_like(z)
    if bounds is not None:
        X0, Y0, X1, Y1 = bounds
        k = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)
        x, y, z, ps, gt = x[k], y[k], z[k], ps[k], gt[k]
    return dict(x=x, y=y, z=z, point_source_id=ps, gps_time=gt)


def rasterize(x, y, value, bounds, res=5.0, agg="median"):
    """Grid a per-point attribute (e.g. change or error) to a raster by per-cell
    ``agg`` ("median" or "mean"). Returns an ny x nx array (NaN where empty).
    Use to turn a point-based change product (m3c2 + lod dims) into GeoTIFFs."""
    X0, Y0, X1, Y1 = bounds
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))
    ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & np.isfinite(value)
    gb = pd.Series(value[ok]).groupby(iy[ok] * nx + ix[ok])
    s = gb.mean() if agg == "mean" else gb.median()
    out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
    return out.reshape(ny, nx)


def heteroscedastic_lod(dod, slope_deg, abs_curv, stable, *, z=1.96):
    """Per-cell level of detection from a calibrated error model (xdem / Hugonnet
    et al., 2022). Models the stable-ground DoD dispersion (NMAD) as a function of
    slope and curvature, then predicts ``z * sigma`` everywhere -- so the error is
    calibrated to the *actual* stable-ground scatter, honestly rising with slope
    rather than inflated by intra-cell relief. Returns None if xdem is
    unavailable (import needs PROJ_DATA unset, as pip rasterio bundles PROJ)."""
    try:
        import xdem.spatialstats as ss
    except Exception:
        return None
    m = stable & np.isfinite(dod) & np.isfinite(slope_deg) & np.isfinite(abs_curv)
    if m.sum() < 500:
        return None
    try:  # the model fit can fail on degenerate inputs (e.g. constant curvature)
        _, errfun = ss._estimate_model_heteroscedasticity(
            dod[m], [slope_deg[m], abs_curv[m]], list_var_names=["slope", "curv"])
        sig = errfun((slope_deg.ravel(), abs_curv.ravel())).reshape(dod.shape)
    except Exception:
        return None
    return z * sig


def difference_dem(before_laz, after_last_laz, bounds, *, res=5.0, ground_q=0.10,
                   correction_surface=True, along_track_drift=True,
                   before_crs=io.MN_2008_CRS):
    """Corrected bare-earth DEM of Difference (after - before).

    ``before_laz``  : 2008-era MN lidar tile (retains point_source_id + gps_time).
    ``after_last_laz``: 3DEP last-return cloud over the same bbox, same CRS.
    ``bounds``      : (minx, miny, maxx, maxy) in the working CRS (EPSG:26915).
    ``ground_q``    : ground percentile (0.10 default; lower = less slope bias,
                      slightly more noise).

    Returns dict: dod, lod (ny x nx arrays), z_after (for hillshade), corrections
    (JSON-serialisable), stable_sigma (empirical 1-sigma on stable ground, m), and
    grid meta (bounds, res, nx, ny).
    """
    X0, Y0, X1, Y1 = bounds
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    def cellstat(x, y, z, how, q=ground_q):
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        gb = pd.Series(z[ok]).groupby(iy[ok] * nx + ix[ok])
        if how == "ground":
            s = gb.quantile(q)
        elif how == "spread":
            s = 1.4826 * (gb.quantile(0.75) - gb.quantile(0.25)) / 1.349
        else:  # count
            s = gb.size()
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
        return out.reshape(ny, nx)

    # --- after (reference) ground + within-cell spread/count ---
    A = read_last_return(after_last_laz, bounds)
    Z21 = cellstat(A["x"], A["y"], A["z"], "ground")
    s21 = cellstat(A["x"], A["y"], A["z"], "spread")
    n21 = cellstat(A["x"], A["y"], A["z"], "count")

    # terrain masks from the reference ground
    Zf = Z21.copy(); nanm = np.isnan(Zf)
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    tpi = Z21 - uniform_filter(Zf, size=int(2 * 300 / res), mode="nearest")
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zf, 2.0), res)[0])
    Zsm = gaussian_filter(Zf, 50 / res / 2)
    lap = (np.gradient(np.gradient(Zsm, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(Zsm, res, axis=1), res, axis=1))
    convex = (sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0)
    stable = ((sdeg < 3) & (tpi > -2)) | convex
    floodplain = np.isfinite(Z21) & (tpi < -2)

    # --- before: align -> tie -> correction surface -> along-track drift ---
    f = laspy.read(str(before_laz))
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); gt8 = np.asarray(f.gps_time)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    be = rn8 == nr8
    pc = io.PointCloud(x8, y8, z8, ps8, np.asarray(f.classification),
                       np.zeros_like(z8), np.zeros_like(ps8), before_crs)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz

    tie = coreg.tie_polynomial(Z21, cellstat(xc[be], yc[be], zc[be], "ground"),
                               res, X0, Y0, order=2)
    xc += coreg.eval_poly_field(tie["a"], xc, yc, tie["norm"], 2)
    yc += coreg.eval_poly_field(tie["b"], xc, yc, tie["norm"], 2)
    zc += coreg.eval_poly_field(tie["c"], xc, yc, tie["norm"], 2)

    if correction_surface:
        C = coreg.correction_surface(Z21, cellstat(xc[be], yc[be], zc[be], "ground"),
                                     res, X0, Y0, radius=400.0, exclude=floodplain)["C"]
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        Cpt = C[iyp, ixp]; zc[np.isfinite(Cpt)] += Cpt[np.isfinite(Cpt)]

    curves = {}
    if along_track_drift:
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        resid = Z21 - cellstat(xc[be], yc[be], zc[be], "ground")
        chg = resid[iyp, ixp]
        stab_pt = be & stable[iyp, ixp] & np.isfinite(chg) & (np.abs(chg) < 0.15)
        drift, curves = coreg.fit_along_track_drift(gt8, chg, stab_pt, ps8)
        zc += drift

    # --- final gridded ground DoD + per-cell LoD ---
    Z08c = cellstat(xc[be], yc[be], zc[be], "ground")
    s08 = cellstat(xc[be], yc[be], zc[be], "spread")
    n08 = cellstat(xc[be], yc[be], zc[be], "count")
    dod = Z21 - Z08c
    r = dod[stable & np.isfinite(dod)]
    sigma = float(1.4826 * np.median(np.abs(r - np.median(r))))
    # LoD: calibrated heteroscedastic model (xdem/Hugonnet 2022) if available,
    # else a within-cell spread proxy (relief-inflated on slopes -- fallback only).
    abs_curv = np.abs(np.gradient(np.gradient(gaussian_filter(Zf, 1.0), res, axis=0), res, axis=0)
                      + np.gradient(np.gradient(gaussian_filter(Zf, 1.0), res, axis=1), res, axis=1))
    lod = heteroscedastic_lod(dod, sdeg, abs_curv, stable)
    lod_method = "xdem heteroscedastic (slope,curv), calibrated on stable ground"
    if lod is None:
        lod = 1.96 * np.sqrt(np.nan_to_num(s08**2 / np.maximum(n08, 1))
                             + np.nan_to_num(s21**2 / np.maximum(n21, 1)))
        lod_method = "within-cell spread proxy (fallback; relief-inflated on slopes)"

    corrections = {
        "epochs": "after - before (positive = deposition)",
        "crs": "EPSG:26915", "res_m": res, "ground_percentile": ground_q,
        "bounds": [float(b) for b in bounds], "stable_1sigma_m": round(sigma, 4),
        "lod_method": lod_method,
        "per_swath_internal_alignment_dxdydz_m":
            {str(k): [round(float(v), 4) for v in val] for k, val in corr.items()},
        "cross_epoch_tie_order2_coef": {
            "dx": [round(float(v), 6) for v in tie["a"]],
            "dy": [round(float(v), 6) for v in tie["b"]],
            "dz": [round(float(v), 6) for v in tie["c"]],
            "norm_xm_xhr_ym_yhr": [round(float(v), 3) for v in tie["norm"]]},
        "along_track_drift_gpsTime_to_m":
            {str(p): {"gps_time": [round(t, 3) for t in c[0]],
                      "drift_m": [round(d, 4) for d in c[1]]} for p, c in curves.items()},
    }
    return dict(dod=dod, lod=lod, z_after=Z21, corrections=corrections,
                stable_sigma=sigma, bounds=tuple(bounds), res=res, nx=nx, ny=ny)
