"""End-to-end differencing of First-Generation (2008-2012) Minnesota lidar against Second-Generation USGS 3DEP.

One function, ``difference_dem``, runs the whole validated workflow and returns a
DEM of Difference plus its corrections and uncertainty. It encodes the lessons
from the Elba pilot; see the module notes below and the project README.

Lessons baked in
----------------
* **Bare earth = last return, ``return_number == number_of_returns`` INCLUDING
  single returns.** Dropping singles empties flat/open ground (they dominate it).
* **Ground = the MEDIAN per cell of the classified ground returns.** Once the CSF
  cloth has removed vegetation, the ground returns scatter symmetrically about the
  true surface, so the unbiased estimate is the central tendency (the median; robust
  to any residual high outlier). This is the single most important choice.
  (History: the earlier heuristic took a LOW PERCENTILE (10th) on RAW last-return
  points -- where true ground sits at the bottom of the return distribution and a low
  pick rejects canopy, the right call BEFORE classification. Kept on top of CSF it
  DOUBLE-COUNTS the cloth and biases the ground low by ~1.28*sigma of the cell
  roughness; because the two epochs differ in roughness/density that offset does not
  cancel -- it reappears as COHERENT slope-correlated false change, ~+8-20 mm,
  physically-impossible ridgetop "deposition". The median removes it.)
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

import warnings

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


def _laz_arrays(f, bounds):
    """(x, y, z, point_source_id, gps_time, classification) from a laspy object,
    plus an in-bounds mask (all-True if bounds is None)."""
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    ps = np.asarray(f.point_source_id)
    try:
        gt = np.asarray(f.gps_time)
    except Exception:
        gt = np.zeros_like(z)
    cls = np.asarray(f.classification)
    if bounds is not None:
        X0, Y0, X1, Y1 = bounds
        inb = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)
    else:
        inb = np.ones(z.shape, bool)
    return x, y, z, ps, gt, cls, inb


def read_after_ground(path, bounds=None, *, mode="class2", csf_pdal=None,
                      min_ground_frac=0.01):
    """gen2 (3DEP) bare-earth points, using the survey's OWN ground classification.

    ``mode="class2"`` (default) selects ASPRS ``Classification == 2`` -- a properly
    QC'd ground-return set produced by USGS from the full-density cloud. It is
    cleaner than the last-return heuristic, which keeps canopy/understory last hits
    (in our forest tiles 24-72% of last returns are NOT ground), so the low-percentile
    ground floats high and the within-cell spread is inflated by vegetation. Every
    class-2 point in 3DEP is itself a last/only return, so this is a strict, cleaner
    subset of last-return, not a different surface.

    Region-level fallback: if the tile carries no usable ground class over the
    requested region (an unclassified 3DEP delivery -- class-2 fraction below
    ``min_ground_frac``), fall back to CSF over the region -- the same cloth filter
    used for the gen1 cloud. gen1 data, which lacks a usable classification, should
    always take that fallback (handled by ``ground_source="csf"`` on the before path).

    ``mode="last_return"`` keeps the legacy ``rn == nr`` heuristic (singles included).

    Returns dict of arrays (x, y, z, point_source_id, gps_time) plus
    ``ground_mode`` = the method actually used ("class2", "csf_fallback", or
    "last_return").
    """
    f = laspy.read(str(path))
    x, y, z, ps, gt, cls, inb = _laz_arrays(f, bounds)
    if mode == "last_return":
        rn = np.asarray(f.return_number); nr = np.asarray(f.number_of_returns)
        m = (rn == nr) & inb
        used = "last_return"
    else:
        g = cls == 2
        nin = int(inb.sum())
        if int((g & inb).sum()) < min_ground_frac * max(nin, 1):
            warnings.warn(
                f"gen2 tile {path}: ASPRS ground (class 2) is "
                f"{100 * (g & inb).sum() / max(nin, 1):.2f}% of the region "
                f"(< {100 * min_ground_frac:.0f}%); treating as unclassified and "
                f"falling back to CSF ground over the region.")
            import os
            import shutil
            tmp = classify_ground_csf(path, pdal=csf_pdal)
            try:
                x, y, z, ps, gt, cls, inb = _laz_arrays(laspy.read(tmp), bounds)
            finally:
                shutil.rmtree(os.path.dirname(tmp), ignore_errors=True)
            m = (cls == 2) & inb          # CSF writes class-2 ground only
            used = "csf_fallback"
        else:
            m = g & inb
            used = "class2"
    return dict(x=x[m], y=y[m], z=z[m], point_source_id=ps[m], gps_time=gt[m],
                ground_mode=used)


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


def _stream_roughness(path, bounds, res, nx, ny, *, after_ground="class2",
                      chunk=8_000_000, min_n=6):
    """Detrended within-cell roughness (see :func:`cell_plane_roughness`) computed
    by STREAMING the gen2 cloud, accumulating the 9 per-cell plane-fit sufficient
    statistics across chunks -- O(cells) RAM, for the large-tile / statewide path.
    Ground selection matches :func:`_stream_ground` (``after_ground``)."""
    X0, Y0, X1, Y1 = bounds; N = nx * ny
    acc = {k: np.zeros(N) for k in
           ("n", "Su", "Sv", "Sz", "Suu", "Suv", "Svv", "Suz", "Svz", "Szz")}
    with laspy.open(str(path)) as fh:
        for pts in fh.chunk_iterator(chunk):
            if after_ground == "last_return":
                rn = np.asarray(pts.return_number); nr = np.asarray(pts.number_of_returns)
                sel = rn == nr
            else:
                sel = np.asarray(pts.classification) == 2
            x = np.asarray(pts.x)[sel]; y = np.asarray(pts.y)[sel]; z = np.asarray(pts.z)[sel]
            ix = np.floor((x - X0) / res).astype(np.int64)
            iy = np.floor((y - Y0) / res).astype(np.int64)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            ix, iy, x, y, z = ix[ok], iy[ok], x[ok], y[ok], z[ok]
            f = iy * nx + ix
            u = x - (X0 + (ix + 0.5) * res); v = y - (Y0 + (iy + 0.5) * res)
            def add(key, w): acc[key] += np.bincount(f, w, minlength=N)
            add("n", np.ones_like(z)); add("Su", u); add("Sv", v); add("Sz", z)
            add("Suu", u * u); add("Suv", u * v); add("Svv", v * v)
            add("Suz", u * z); add("Svz", v * z); add("Szz", z * z)
    n = acc["n"]
    M = np.empty((N, 3, 3))
    M[:, 0, 0] = acc["Suu"]; M[:, 0, 1] = acc["Suv"]; M[:, 0, 2] = acc["Su"]
    M[:, 1, 0] = acc["Suv"]; M[:, 1, 1] = acc["Svv"]; M[:, 1, 2] = acc["Sv"]
    M[:, 2, 0] = acc["Su"];  M[:, 2, 1] = acc["Sv"];  M[:, 2, 2] = n
    r = np.stack([acc["Suz"], acc["Svz"], acc["Sz"]], axis=1)
    valid = (n >= min_n) & (np.abs(np.linalg.det(M)) > 1e-6)
    rms = np.full(N, np.nan)
    if valid.any():
        beta = np.linalg.solve(M[valid], r[valid])
        rss = acc["Szz"][valid] - np.einsum("ij,ij->i", beta, r[valid])
        rms[valid] = np.sqrt(np.maximum(rss, 0.0) / n[valid])
    return rms.reshape(ny, nx)


def cell_plane_roughness(x, y, z, X0, Y0, res, nx, ny, *, min_n=6):
    """Detrended within-cell roughness: per-cell RMS of the residual to a plane
    fitted to the ground points in that cell.

    Raw within-cell spread of the ground returns is dominated by RELIEF -- a 5 m
    cell on a 30 deg slope spans ~2.9 m vertically from tilt alone -- so it cannot
    serve as an error covariate without removing the slope first (else slope
    masquerades as roughness and is collinear with the slope covariate already in
    the LoD model). We remove it exactly, at the analysis-cell scale, by fitting a
    plane per cell and taking the RMS of the residuals -- the standard geomorphic
    definition of surface roughness (Grohmann et al. 2011; Cavalli et al. 2008),
    one scale below the slope-normal ground fix. What survives is the
    slope-independent scatter: micro-topography plus the vegetation-penetration
    ambiguity that makes a cell's ground return set fuzzy -- the physical signal
    that a cell's ground elevation is poorly determined, hence the natural LoD
    covariate (Wheaton et al. 2010 use slope, point density AND roughness).

    Computed from streaming-friendly sufficient statistics (the 9 per-cell moments
    of a plane fit), so it needs no per-point storage. Fit ``z ~ a*u + b*v + c`` on
    the in-cell-centered (u, v); ``rms = sqrt(RSS / n)``. Returns an ny x nx array,
    NaN where a cell has < ``min_n`` points or is ill-conditioned (collinear).
    """
    ix = np.floor((x - X0) / res).astype(np.int64)
    iy = np.floor((y - Y0) / res).astype(np.int64)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    ix, iy = ix[ok], iy[ok]; x, y, z = x[ok], y[ok], z[ok]
    f = iy * nx + ix; N = nx * ny
    u = x - (X0 + (ix + 0.5) * res); v = y - (Y0 + (iy + 0.5) * res)
    S = lambda w: np.bincount(f, w, minlength=N)
    n = S(np.ones_like(z)); Su = S(u); Sv = S(v); Sz = S(z)
    Suu = S(u * u); Suv = S(u * v); Svv = S(v * v)
    Suz = S(u * z); Svz = S(v * z); Szz = S(z * z)
    M = np.empty((N, 3, 3))
    M[:, 0, 0] = Suu; M[:, 0, 1] = Suv; M[:, 0, 2] = Su
    M[:, 1, 0] = Suv; M[:, 1, 1] = Svv; M[:, 1, 2] = Sv
    M[:, 2, 0] = Su;  M[:, 2, 1] = Sv;  M[:, 2, 2] = n
    r = np.stack([Suz, Svz, Sz], axis=1)
    valid = (n >= min_n) & (np.abs(np.linalg.det(M)) > 1e-6)
    rms = np.full(N, np.nan)
    if valid.any():
        beta = np.linalg.solve(M[valid], r[valid])
        rss = Szz[valid] - np.einsum("ij,ij->i", beta, r[valid])
        rms[valid] = np.sqrt(np.maximum(rss, 0.0) / n[valid])
    return rms.reshape(ny, nx)


def heteroscedastic_lod(dod, slope_deg, abs_curv, stable, *, stderr=None,
                        density=None, roughness=None, z=1.96):
    """Per-cell level of detection from a calibrated error model (xdem / Hugonnet
    et al., 2022). Models the stable-ground DoD dispersion (NMAD) as a function of
    slope, curvature, and a within-cell ground-uncertainty covariate, then predicts
    ``z * sigma`` everywhere, calibrated to the *actual* stable-ground scatter.

    ``stderr`` is the BASIS covariate the pipeline uses: the cell's ground-estimate
    STANDARD ERROR, which combines the two distinct pieces of within-cell
    information -- surface variability and sample size -- in the form the statistics
    dictate, ``sqrt(sum_epoch roughness^2 / n)`` (per-epoch standard errors in
    quadrature; see :func:`difference_dem`). Detrended ROUGHNESS (:func:`cell_plane_roughness`,
    slope removed) is the numerator: the real internal variability of the surface
    being sampled. DENSITY (ground returns per cell) is the denominator: how well
    that variability is pinned down. They are genuinely different information
    (Aguilar et al. 2005 rank morphology and sampling density as separate,
    both-significant factors; Wheaton et al. 2010 keep density and roughness as
    distinct FIS inputs), so one scales the error while the other damps only its
    sampling part. Folding them into a single physical covariate uses both without
    fragmenting xdem's N-D bin space (a 4th independent covariate makes the Delaunay
    fit fail on small stable sets).

    ``density`` and ``roughness`` are OPTIONAL standalone covariates, retained for
    flexibility (e.g. reproducing the density-only model, or a roughness-only run);
    the pipeline passes ``stderr`` instead. Any given covariate is added to the fit.

    The result flows straight into detection via ``perror = lod/1.96``. Returns None
    if xdem is unavailable (import needs PROJ_DATA unset, as pip rasterio bundles PROJ)."""
    try:
        import xdem.spatialstats as ss
    except Exception:
        return None
    covs = [slope_deg, abs_curv]; names = ["slope", "curv"]
    for val, nm in ((stderr, "stderr"), (density, "density"), (roughness, "roughness")):
        if val is not None:
            covs = covs + [val]; names = names + [nm]
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
                   coarse_bins=120, bw=0.02, down=3.0, up=2.0, after_ground="class2"):
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
    ``after_ground`` selects the gen2 ground returns per chunk: "class2" (default)
    uses the survey's ASPRS ``Classification == 2`` (cleaner than last-return; see
    :func:`read_after_ground`); "last_return" uses the ``rn == nr`` heuristic. The
    CSF fallback for an unclassified tile lives in the in-memory path only -- a
    streamed statewide run assumes the 3DEP delivery is classified.

    Returns (ground, spread, count) as ny x nx arrays.
    """
    X0, Y0, X1, Y1 = bounds
    N = nx * ny

    def chunks():
        with laspy.open(str(path)) as fh:
            for pts in fh.chunk_iterator(chunk):
                if after_ground == "last_return":
                    rn = np.asarray(pts.return_number); nr = np.asarray(pts.number_of_returns)
                    sel = rn == nr
                else:
                    sel = np.asarray(pts.classification) == 2
                x = np.asarray(pts.x)[sel]; y = np.asarray(pts.y)[sel]; z = np.asarray(pts.z)[sel]
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


def difference_dem(before_laz, after_laz, bounds, *, res=5.0, ground_q=0.50,
                   correction_surface=False, along_track_drift=True,
                   ground="slope_normal", sn_smooth_cells=1.2, stream=False,
                   ground_source="csf", after_ground="class2", csf_pdal=None,
                   robust_stable=True, before_crs=io.MN_GEN1_CRS):
    """Corrected bare-earth DEM of Difference (after - before).

    ``before_laz``  : first-generation (gen1) MN lidar tile (retains point_source_id + gps_time).
    ``after_laz``   : gen2 3DEP cloud over the same bbox, same CRS, WITH its ASPRS
                      classification intact (pass the full delivery, not a
                      pre-filtered last-return file -- see ``after_ground``).
    ``bounds``      : (minx, miny, maxx, maxy) in the working CRS (EPSG:26915).
    ``ground_q``    : per-cell ground quantile for GRIDDING. **0.50 (median) by
                      default** -- the CSF/class-2 ground returns are already
                      vegetation-free, so they scatter symmetrically about the
                      surface and the median is the unbiased estimate. A low
                      percentile (e.g. the legacy 0.10) sits ~1.28*sigma below the
                      surface, double-counts the cloth, and biases the ground low by
                      an epoch-dependent, roughness-/slope-growing amount that does
                      not cancel in the difference.
    ``ground``      : ground GRIDDING estimator. "slope_normal" (default) = the
                      ``ground_q`` quantile of the residual to a common smoothed
                      regional surface (both epochs), which removes the downhill bias
                      a horizontal pick has on a slope (a low pick necessarily selects
                      the downhill-lowest points); the shared surface cancels in the
                      difference. "low_q" = the quantile of raw z per horizontal cell
                      (the older heuristic). ``sn_smooth_cells`` sets the
                      regional-slope smoothing (in cells).
    ``ground_source``: how the before-epoch bare-earth is obtained. "csf" (default)
                      runs PDAL CSF (tuned for sparse steep/wooded terrain) for a
                      cleaner, more general ground -- SLOW (min/tile) and needs PDAL.
                      "last_return" uses the raw last-return heuristic (fast, no
                      dependency; near-identical DoD, so choose it to skip CSF).
                      ``csf_pdal`` optionally points to the PDAL binary.
    ``after_ground``: how the gen2 (3DEP) bare-earth is obtained. "class2" (default)
                      uses the survey's OWN ASPRS ground classification, a QC'd
                      ground-return set that is cleaner than last-return (which keeps
                      canopy/understory last hits -- 24-72% of last returns are
                      non-ground in our forest tiles, inflating the low-percentile
                      ground and the within-cell spread). Region-level fallback to
                      CSF if the tile is unclassified. "last_return" keeps the legacy
                      ``rn == nr`` heuristic. See :func:`read_after_ground`.
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
        Z21, s21, n21 = _stream_ground(after_laz, bounds, res, nx, ny, ground_q,
                                       after_ground=after_ground)
        r21 = _stream_roughness(after_laz, bounds, res, nx, ny, after_ground=after_ground)
    else:
        A = read_after_ground(after_laz, bounds, mode=after_ground, csf_pdal=csf_pdal)
        print(f"  gen2 ground: {A['ground_mode']} ({A['x'].size} points)", flush=True)
        Z21 = cellstat(A["x"], A["y"], A["z"], "ground")
        s21 = cellstat(A["x"], A["y"], A["z"], "spread")
        n21 = cellstat(A["x"], A["y"], A["z"], "count")
        r21 = cell_plane_roughness(A["x"], A["y"], A["z"], X0, Y0, res, nx, ny)

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

    # ground estimator: "low_q" (horizontal ground_q quantile) or "slope_normal"
    # (ground_q quantile -- median by default -- of the residual to a common smoothed
    # regional surface, which removes the downhill bias of a horizontal pick on a
    # slope). The shared
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
        Zref = (_stream_ground(after_laz, bounds, res, nx, ny, ground_q,
                               plane=(Zreg_f, dzde, dzdn), after_ground=after_ground)[0]
                if stream else groundg(A["x"], A["y"], A["z"]))
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
    r08 = cell_plane_roughness(xc[be], yc[be], zc[be], X0, Y0, res, nx, ny)
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
    # Standard-error LoD covariate: the cell ground-estimate's DoD standard error,
    # combining the two DISTINCT pieces of within-cell information in the form the
    # statistics dictate -- detrended ROUGHNESS (surface variability, the numerator)
    # over ground-return DENSITY (sample size, the denominator): the per-epoch
    # standard errors sqrt(roughness^2 / n) added in quadrature. Roughness and
    # density are not redundant (Aguilar 2005; Wheaton 2010): one scales the error,
    # the other damps only its sampling part. Using the combined SE keeps both while
    # feeding xdem a single covariate (a 4th independent covariate fragments its N-D
    # bin space and the Delaunay fit fails on small stable sets). Roughness is
    # detrended (cell_plane_roughness), so this is the principled form of the old
    # relief-inflated within-cell-spread proxy.
    stderr = np.sqrt(np.nan_to_num(r08) ** 2 / np.maximum(n08, 1.0)
                     + np.nan_to_num(r21) ** 2 / np.maximum(n21, 1.0))
    stderr[~(np.isfinite(r08) | np.isfinite(r21))] = np.nan
    lod = heteroscedastic_lod(dod, sdeg, abs_curv, stable_rep, stderr=stderr)
    lod_method = ("xdem heteroscedastic (slope,curv,standard-error[roughness/sqrt(density)]), "
                  "calibrated on stable ground")
    if lod is None:                                   # stderr model degenerate -> slope,curv only
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
