"""End-to-end differencing of First-Generation (2008-2012) Minnesota lidar against Second-Generation USGS 3DEP.

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
* **Correct in the ACQUISITION frame, per point, BEFORE gridding:** per-swath
  internal alignment (translation) -> spatially varying quadratic tie -> per-swath
  along-track GNSS-drift spline ``f(gps_time)``. The residual warp and real
  localized change share the same ~100-400 m scale, so no data-driven interpolator
  on the elevation residual can separate them; only the acquisition geometry can.
  The drift uses it (per-swath, time-ordered), is deterministic and reusable, and
  cannot absorb a localized deposit. A DeLong 400 m correction surface is available
  (``correction_surface=True``) for legacy data lacking ``gps_time``, but it is a
  data-driven IDW that absorbs localized flat change up to its dz threshold, adds
  only ~4 mm here, and is OFF by default.
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
from .ground import classify_ground_csf


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


def heteroscedastic_lod(dod, slope_deg, abs_curv, stable, *, density=None, z=1.96):
    """Per-cell level of detection from a calibrated error model (xdem / Hugonnet
    et al., 2022). Models the stable-ground DoD dispersion (NMAD) as a function of
    slope, curvature, and -- if given -- **ground-return DENSITY**, then predicts
    ``z * sigma`` everywhere, calibrated to the *actual* stable-ground scatter.

    ``density`` is a raw physical covariate: the number of lidar pulses that
    reached the ground per cell (the limiting epoch). It is NOT a land-cover class
    -- nothing is classified. Where the sparser survey penetrates poorly (dense
    canopy), the ground estimate rests on few returns and is noisier, so sigma
    honestly rises there. That drops the forest 'speckle' (which is
    penetration-noise, not ground change) as a modeled *error*, not a mask, and it
    flows straight into detection via ``perror = lod/1.96``. Returns None if xdem
    is unavailable (import needs PROJ_DATA unset, as pip rasterio bundles PROJ)."""
    try:
        import xdem.spatialstats as ss
    except Exception:
        return None
    covs = [slope_deg, abs_curv]; names = ["slope", "curv"]
    if density is not None:
        covs = covs + [density]; names = names + ["density"]
    m = stable & np.isfinite(dod)
    for c in covs:
        m = m & np.isfinite(c)
    if m.sum() < 500:
        return None
    try:  # the model fit can fail on degenerate inputs (e.g. a constant covariate)
        _, errfun = ss._estimate_model_heteroscedasticity(
            dod[m], [c[m] for c in covs], list_var_names=names)
        sig = errfun(tuple(c.ravel() for c in covs)).reshape(dod.shape)
    except Exception:
        return None
    return z * sig


def _stream_ground(path, bounds, res, nx, ny, q, *, plane=None, chunk=8_000_000,
                   coarse_bins=120, bw=0.02, down=3.0, up=2.0):
    """Per-cell low-q ground, spread, and count by STREAMING the cloud in chunks,
    so peak RAM is O(cells), not O(points) -- for statewide runs where the dense
    3DEP cloud will not fit in memory.

    Reads ``path`` in chunks (never holds the whole cloud). ``plane`` = flat
    per-cell (Z_reg, dz_deast, dz_dnorth) turns this into the slope-normal residual
    ground (``ground="slope_normal"``): low-q of ``z - regional plane``, plus the
    plane back. Blunder-robust via a coarse-histogram anchor + a downward-widened
    fine window, then read the q-th percentile off the per-cell CDF with in-bin
    interpolation. Matches an exact ``groupby.quantile`` to ~mm on well-sampled
    cells; SPARSE cells (few points) can differ (histogram vs the exact's linear
    interpolation across large gaps) and should be dropped by a min-count mask.
    Returns (ground, spread, count) as ny x nx arrays.
    """
    X0, Y0, X1, Y1 = bounds
    N = nx * ny

    def chunks():
        with laspy.open(str(path)) as fh:
            for pts in fh.chunk_iterator(chunk):
                rn = np.asarray(pts.return_number); nr = np.asarray(pts.number_of_returns)
                last = rn == nr
                x = np.asarray(pts.x)[last]; y = np.asarray(pts.y)[last]; z = np.asarray(pts.z)[last]
                ix = ((x - X0) / res).astype(np.int64); iy = ((y - Y0) / res).astype(np.int64)
                ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
                f = iy[ok] * nx + ix[ok]; v = z[ok]
                if plane is not None:
                    Zc, dzde, dzdn = plane
                    dxe = x[ok] - (X0 + (ix[ok] + 0.5) * res)
                    dyn = y[ok] - (Y0 + (iy[ok] + 0.5) * res)
                    v = v - (Zc[f] + dxe * dzde[f] + dyn * dzdn[f])
                yield f, v

    lo = np.full(N, np.inf); hi = np.full(N, -np.inf)              # pass 1: min/max
    for f, v in chunks():
        gmn = pd.Series(v).groupby(f).min(); gmx = pd.Series(v).groupby(f).max()
        np.minimum.at(lo, gmn.index.values, gmn.values)
        np.maximum.at(hi, gmx.index.values, gmx.values)
    span = np.where(np.isfinite(lo) & (hi > lo), hi - lo, 1.0)
    CB = coarse_bins; chist = np.zeros(N * CB, np.int64)           # pass 2: coarse anchor
    for f, v in chunks():
        b = np.clip(((v - lo[f]) / span[f] * CB).astype(np.int64), 0, CB - 1)
        chist += np.bincount(f * CB + b, minlength=N * CB)
    chist = chist.reshape(N, CB); ccdf = np.cumsum(chist, axis=1)
    ntot = ccdf[:, -1].astype(float)

    with np.errstate(invalid="ignore"):
        def coarse_pct(p):
            b = np.argmax(ccdf >= (p * ntot)[:, None], axis=1)
            return np.where(np.isfinite(lo), lo + (b + 0.5) * span / CB, np.nan)
        anchor = coarse_pct(q)
        spread = 1.4826 * (coarse_pct(0.75) - coarse_pct(0.25)) / 1.349
    flo = anchor - down; SPAN = down + up; FB = int(round(SPAN / bw))  # pass 3: fine window
    below = np.zeros(N, np.int64); fhist = np.zeros(N * FB, np.int64)
    for f, v in chunks():
        d = v - flo[f]; inw = (d >= 0) & (d < SPAN)
        below += np.bincount(f[d < 0], minlength=N)
        ff = f[inw]; b = np.clip((d[inw] / bw).astype(np.int64), 0, FB - 1)
        fhist += np.bincount(ff * FB + b, minlength=N * FB)
    fhist = fhist.reshape(N, FB); fcdf = below[:, None] + np.cumsum(fhist, axis=1)
    tgt = q * ntot
    bf = np.argmax(fcdf >= tgt[:, None], axis=1)
    cprev = np.where(bf > 0, fcdf[np.arange(N), np.clip(bf - 1, 0, FB - 1)], below).astype(float)
    hb = fhist[np.arange(N), bf].astype(float)
    frac = np.divide(tgt - cprev, hb, out=np.zeros(N), where=hb > 0)
    g = np.where(ntot > 0, flo + (bf + frac) * bw, np.nan)
    if plane is not None:
        g = np.where(np.isfinite(g), plane[0] + g, np.nan)
    cnt = np.where(ntot > 0, ntot, np.nan)
    return g.reshape(ny, nx), spread.reshape(ny, nx), cnt.reshape(ny, nx)


def difference_dem(before_laz, after_last_laz, bounds, *, res=5.0, ground_q=0.10,
                   correction_surface=False, along_track_drift=True,
                   ground="slope_normal", sn_smooth_cells=1.2, stream=False,
                   ground_source="csf", csf_pdal=None, robust_stable=True,
                   before_crs=io.MN_GEN1_CRS):
    """Corrected bare-earth DEM of Difference (after - before).

    ``before_laz``  : first-generation (gen1) MN lidar tile (retains point_source_id + gps_time).
    ``after_last_laz``: gen2 3DEP last-return cloud over the same bbox, same CRS.
    ``bounds``      : (minx, miny, maxx, maxy) in the working CRS (EPSG:26915).
    ``ground_q``    : ground percentile (0.10 default; lower = less slope bias,
                      slightly more noise).
    ``ground``      : ground GRIDDING estimator. "slope_normal" (default) = low
                      percentile of the residual to a common smoothed regional
                      surface (both epochs), which removes the downhill bias a
                      horizontal low-pick has on a slope (it necessarily selects the
                      downhill-lowest points); the shared surface cancels in the
                      difference. "low_q" = low percentile of raw z per horizontal
                      cell (the older heuristic). ``sn_smooth_cells`` sets the
                      regional-slope smoothing (in cells).
    ``ground_source``: how the before-epoch bare-earth is obtained. "csf" (default)
                      runs PDAL CSF (tuned for sparse steep/wooded terrain) for a
                      cleaner, more general ground -- SLOW (min/tile) and needs PDAL.
                      "last_return" uses the raw last-return heuristic (fast, no
                      dependency; near-identical DoD, so choose it to skip CSF).
                      ``csf_pdal`` optionally points to the PDAL binary.
    ``robust_stable``: if True (default), the stable-ground mask used to REPORT
                      uncertainty (stable_sigma, the LoD calibration) is refined by
                      an iterative 3-NMAD sigma-clip of the DoD, removing real
                      change that the geometric mask admits. This matters where the
                      geometric heuristic fails -- e.g. a valley wider than the
                      600 m TPI window, whose flat floodplain interior reads as
                      "stable" and pulls its own change into the calibration (~37%
                      of stable cells on the MN River valley pilot; <5% at Elba).
                      The CORRECTIONS (tie, drift) are already robust to this (the
                      tie fits only sloped cells with its own NMAD rejection; the
                      drift gates on |change|), so this only cleans the reporting
                      layer and never alters the DoD surface.

    Returns dict: dod, lod (ny x nx arrays), z_after (for hillshade), stable (the
    reporting stable mask), corrections (JSON-serialisable), stable_sigma (empirical
    1-sigma on stable ground, m), and grid meta (bounds, res, nx, ny).
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
    # stream=True grids the (dense) after cloud in chunks (O(cells) RAM) for
    # statewide runs; else load it in memory. A stays None in streaming mode.
    A = None
    if stream:
        Z21, s21, n21 = _stream_ground(after_last_laz, bounds, res, nx, ny, ground_q)
    else:
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

    # ground estimator: "low_q" (horizontal low percentile) or "slope_normal"
    # (low percentile of the residual to a common smoothed regional surface, which
    # removes the downhill bias of a horizontal low-pick on a slope). The shared
    # surface Zreg is the smoothed reference ground, so it cancels in after - before.
    if ground == "slope_normal":
        Zreg = gaussian_filter(Zf, sn_smooth_cells)
        dzde = np.gradient(Zreg, res, axis=1).ravel()   # d/deast (columns)
        dzdn = np.gradient(Zreg, res, axis=0).ravel()   # d/dnorth (rows; iy grows north)
        Zreg_f = Zreg.ravel()

    def groundg(x, y, z):
        if ground != "slope_normal":
            return cellstat(x, y, z, "ground")
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        f = iy[ok] * nx + ix[ok]
        dxe = x[ok] - (X0 + (ix[ok] + 0.5) * res)
        dyn = y[ok] - (Y0 + (iy[ok] + 0.5) * res)
        resid = z[ok] - (Zreg_f[f] + dxe * dzde[f] + dyn * dzdn[f])
        s = pd.Series(resid).groupby(f).quantile(ground_q)
        out = np.full(nx * ny, np.nan)
        out[s.index.values] = Zreg_f[s.index.values] + s.values
        return out.reshape(ny, nx)

    # after-epoch reference ground in the chosen estimator (== Z21 for low_q)
    if ground == "slope_normal":
        Zref = (_stream_ground(after_last_laz, bounds, res, nx, ny, ground_q,
                               plane=(Zreg_f, dzde, dzdn))[0] if stream
                else groundg(A["x"], A["y"], A["z"]))
    else:
        Zref = Z21

    # --- before: (CSF ground classification) -> align -> tie -> drift ---
    # ground_source="csf" (default) runs PDAL CSF on the before cloud first for a
    # cleaner, more general bare-earth (removes structures/understory); "last_return"
    # skips it and uses the raw last-return heuristic. CSF is slow (min/tile).
    _csf_tmp = None
    if ground_source == "csf":
        _csf_tmp = classify_ground_csf(before_laz, pdal=csf_pdal)
        before_laz = _csf_tmp
    f = laspy.read(str(before_laz))
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); gt8 = np.asarray(f.gps_time)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    be = rn8 == nr8
    if _csf_tmp is not None:
        import shutil, os
        shutil.rmtree(os.path.dirname(_csf_tmp), ignore_errors=True)
    pc = io.PointCloud(x8, y8, z8, ps8, np.asarray(f.classification),
                       np.zeros_like(z8), np.zeros_like(ps8), before_crs)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz

    tie = coreg.tie_polynomial(Zref, groundg(xc[be], yc[be], zc[be]),
                               res, X0, Y0, order=2)
    xc += coreg.eval_poly_field(tie["a"], xc, yc, tie["norm"], 2)
    yc += coreg.eval_poly_field(tie["b"], xc, yc, tie["norm"], 2)
    zc += coreg.eval_poly_field(tie["c"], xc, yc, tie["norm"], 2)

    if correction_surface:
        C = coreg.correction_surface(Zref, groundg(xc[be], yc[be], zc[be]),
                                     res, X0, Y0, radius=400.0, exclude=floodplain)["C"]
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        Cpt = C[iyp, ixp]; zc[np.isfinite(Cpt)] += Cpt[np.isfinite(Cpt)]

    curves = {}
    if along_track_drift:
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        resid = Zref - groundg(xc[be], yc[be], zc[be])
        chg = resid[iyp, ixp]
        stab_pt = be & stable[iyp, ixp] & np.isfinite(chg) & (np.abs(chg) < 0.15)
        drift, curves = coreg.fit_along_track_drift(gt8, chg, stab_pt, ps8)
        zc += drift

    # --- final gridded ground DoD + per-cell LoD ---
    Z08c = groundg(xc[be], yc[be], zc[be])
    s08 = cellstat(xc[be], yc[be], zc[be], "spread")
    n08 = cellstat(xc[be], yc[be], zc[be], "count")
    dod = Zref - Z08c
    # Reporting stable mask. The geometric `stable` admits real change where its
    # heuristics fail (a floodplain wider than the TPI window reads as flat-stable),
    # which inflates sigma and the LoD calibration. Refine it by an iterative
    # 3-NMAD sigma-clip of the DoD so the reported error is the true stable-ground
    # error, not the change bleeding into it. Corrections are untouched (already
    # robust), so the DoD surface is identical either way.
    stable_rep = stable & np.isfinite(dod)
    stable_geom_n = int(stable_rep.sum())
    if robust_stable:
        for _ in range(8):
            v = dod[stable_rep]
            med = np.median(v); nm = 1.4826 * np.median(np.abs(v - med))
            keep = stable_rep & (np.abs(dod - med) < 3.0 * max(nm, 1e-3))
            if keep.sum() == stable_rep.sum():
                break
            stable_rep = keep
    stable_clip_frac = (1.0 - stable_rep.sum() / stable_geom_n) if stable_geom_n else 0.0
    r = dod[stable_rep]
    sigma = float(1.4826 * np.median(np.abs(r - np.median(r))))
    # LoD: calibrated heteroscedastic model (xdem/Hugonnet 2022) if available,
    # else a within-cell spread proxy (relief-inflated on slopes -- fallback only).
    abs_curv = np.abs(np.gradient(np.gradient(gaussian_filter(Zf, 1.0), res, axis=0), res, axis=0)
                      + np.gradient(np.gradient(gaussian_filter(Zf, 1.0), res, axis=1), res, axis=1))
    # ground-return DENSITY covariate: the LIMITING (sparser) epoch's ground-return
    # count per cell -- a raw physical measurement (no classification). Low density
    # (poor canopy penetration) => noisier ground estimate => honestly larger LoD,
    # which drops forest 'speckle' as modeled error rather than a mask.
    density = np.minimum(np.nan_to_num(n08, nan=0.0), np.nan_to_num(n21, nan=0.0))
    lod = heteroscedastic_lod(dod, sdeg, abs_curv, stable_rep, density=density)
    lod_method = "xdem heteroscedastic (slope,curv,ground-density), calibrated on stable ground"
    if lod is None:                                   # density model degenerate -> slope,curv only
        lod = heteroscedastic_lod(dod, sdeg, abs_curv, stable_rep)
        lod_method = "xdem heteroscedastic (slope,curv), calibrated on stable ground"
    if lod is None:
        lod = 1.96 * np.sqrt(np.nan_to_num(s08**2 / np.maximum(n08, 1))
                             + np.nan_to_num(s21**2 / np.maximum(n21, 1)))
        lod_method = "within-cell spread proxy (fallback; relief-inflated on slopes)"

    corrections = {
        "epochs": "after - before (positive = deposition)",
        "crs": "EPSG:26915", "res_m": res, "ground_percentile": ground_q,
        "ground_estimator": ground, "ground_source": ground_source,
        "bounds": [float(b) for b in bounds], "stable_1sigma_m": round(sigma, 4),
        "robust_stable": robust_stable,
        "stable_clip_fraction": round(float(stable_clip_frac), 4),
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
    return dict(dod=dod, lod=lod, z_after=Z21, stable=stable_rep,
                corrections=corrections, stable_sigma=sigma,
                bounds=tuple(bounds), res=res, nx=nx, ny=ny)
