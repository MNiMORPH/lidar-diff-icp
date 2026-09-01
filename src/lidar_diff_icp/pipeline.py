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
  internal alignment (translation) -> one constant Nuth & Kaeaeb lateral shift
  (order-0 tie) + a geoid-model vertical datum -> per-swath along-track GNSS-drift
  spline ``f(gps_time)``. (The earlier spatially-varying quadratic/parabola tie
  and the reference-plane fit were REMOVED -- geoid-only datum now; see git
  history.) The residual warp and real
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
                   correction_surface=False, along_track_drift=True, tie="reference",
                   ground="slope_normal", sn_smooth_cells=1.2, stream=False,
                   ground_source="csf", after_ground="class2", csf_pdal=None,
                   csf_cache=None, robust_stable=True, before_crs=io.MN_GEN1_CRS,
                   geoid_datum=None, correct_boresight=False,
                   boresight_roll_mm_per_deg=None, swath_tie="intercept",
                   absolute_datum=None):
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
    ``tie``         : must be "reference" (the only value). The cross-epoch datum is:
                      align to the lowest-numbered flightline, apply the Nuth-Kaeaeb
                      lateral (x,y) shift, then add the geoid-model z offset
                      N_gen1 - N_gen2 (auto-computed from the PROJ geoid grids by
                      :func:`lidar_diff_icp.references.geoid_difference` unless
                      ``geoid_datum`` is supplied). No plane is fitted to "stable"
                      surfaces; residual offsets are left for later analysis. The
                      reference_plane fit and the parabola tie were removed (git history).
    **The absolute level of the output depends on the ZERO LINE and is not a measured
    elevation.** ``align_swaths`` uses the lowest-numbered flight line as the ZERO LINE, so the
    mosaic inherits that line's own vertical error; re-gauging on another line shifts every
    elevation (measured span 44.60 mm at elbaext). Swath-to-swath DIFFERENCES are
    unaffected. To obtain an absolute elevation, apply a ground-control datum constant
    measured against this zero line -- ``ground_control/apply_datum.py``. The output records
    ``zero_line`` and leaves ``absolute_datum_mm`` None until one is supplied.

    ``absolute_datum`` : optional dict placing BOTH epochs on surveyed NAVD88, making the
                      output's elevation independent of the zero line. Keys:
                      ``gen1_mm``  constant to ADD to gen1 **as it sits in this DoD**, i.e.
                                   AFTER the geoid shift (= c1_own_frame − geoid_mm);
                      ``gen2_mm``  constant to ADD to gen2;
                      ``zero_line`` the ZERO LINE the gen1 constant was measured against -- the
                      flight line defined as zero when that tile's swath network
                      was solved (ground_control still names this key zero_line)
                                   — CHECKED against this run's zero line and raises on a
                                   mismatch, because a constant measured against another
                                   zero line belongs to a different product;
                      ``source``   where the constants came from.
                      Optional: ``gen1_sigma_mm``, ``gen2_sigma_mm``.
                      ``None`` (default) leaves the product exactly as before and records
                      ``absolute_datum_mm: None``. Constants are measured with
                      ``ground_control/`` and built by ``apply_datum.datum_for_pipeline``.

    ``swath_tie``   : how the VERTICAL offset of each flight-line pair is reduced from
                      their overlap. "intercept" (default) is the LAD (median-regression)
                      intercept at across-track position zero, i.e. at
                      ``tan(scan_ref) = tan(scan_src)``. "overlap_median" is the older
                      plain ``median(dh)``, which equals ``k + c*mean(dtan)`` whenever the
                      between-line difference has an across-track slope ``c`` and is
                      therefore EXTENT-DEPENDENT: two tiles covering different parts of
                      the same sidelap get different constants for the same pair of lines.
                      The LAD fit reduces exactly to the median when ``c = 0``. Recorded in
                      ``corrections.json`` as ``swath_tie``; pass "overlap_median" to
                      reproduce products built before 2026-08-26.
                      See :func:`lidar_diff_icp.coreg.across_track_tie`.
    ``geoid_datum`` : optional ``(const_m, b, c)`` geoid shift to ADD to gen1 (const +
                      E/N tilt in m/km about the bounds centroid). Auto-computed from the
                      geoid grids when None (default).
    ``correct_boresight`` : if True, remove a scanner boresight ROLL (elevation error
                      proportional to scan angle, a sensor constant) from gen1 per point,
                      ``z -= b*scan_angle``, BEFORE the empirical alignment -- an
                      instrumental term applied first. Default False (opt-in): it removes
                      the cross-track scan-angle asymmetry but its tile-wide DoD footprint
                      is small (self-cancels in overlap medians / absorbed by per-swath
                      alignment). See :mod:`lidar_diff_icp.boresight`.
    ``boresight_roll_mm_per_deg`` : the roll to apply (mm/deg). If None (and
                      ``correct_boresight``), it is self-calibrated from gen1 flight-line
                      overlap via :func:`coreg.estimate_boresight_roll`; pass a value to
                      reuse a lift-wide sensor constant. Recorded in ``corrections.json``.
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
    ``csf_cache``   : path to ARCHIVE and REUSE the raw CSF ground. If set and the file
                      exists, it is loaded (CSF skipped, ~20 s); if set and absent, CSF
                      runs once, is saved there, and KEPT. The move (lateral shift) and
                      tilt (datum) corrections and the drift are applied to this kept
                      cloud downstream, so iterating on the correction never re-runs the
                      slow deterministic classification.
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

    def _plane_ground(x, y, z, minpts=4):
        """Hill-normal ground: per-cell least-squares plane z = c + a*dE + b*dN fit to
        the cell's own points, read at the cell CENTRE (c). Unbiased under a tilt
        regardless of where the (sparse / occluded) points fall in the cell -- unlike a
        quantile of residuals, which lands at a tilt-correlated spot on a steep cell.
        Falls back to the cell median where too few points for a stable plane."""
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        f = iy[ok] * nx + ix[ok]
        u = x[ok] - (X0 + (ix[ok] + 0.5) * res)      # dEast from cell centre
        v = y[ok] - (Y0 + (iy[ok] + 0.5) * res)      # dNorth from cell centre
        zz = z[ok]; N = nx * ny
        n = np.bincount(f, minlength=N).astype(float)
        Su = np.bincount(f, u, N); Sv = np.bincount(f, v, N)
        Suu = np.bincount(f, u * u, N); Svv = np.bincount(f, v * v, N); Suv = np.bincount(f, u * v, N)
        Sz = np.bincount(f, zz, N); Suz = np.bincount(f, u * zz, N); Svz = np.bincount(f, v * zz, N)
        M = np.stack([np.stack([n, Su, Sv], 1), np.stack([Su, Suu, Suv], 1),
                      np.stack([Sv, Suv, Svv], 1)], 1)          # (N,3,3) normal equations
        rhs = np.stack([Sz, Suz, Svz], 1)                        # (N,3)
        valid = (n >= minpts) & (np.abs(np.linalg.det(M)) > 1e-9)
        c = np.full(N, np.nan)
        if valid.any():
            c[valid] = np.linalg.solve(M[valid], rhs[valid])[:, 0]
        med = cellstat(x, y, z, "ground").ravel()                # fallback where too sparse
        c[~valid] = med[~valid]
        return c.reshape(ny, nx)

    def _poly2_ground(x, y, z, minpts=18):
        """Windowed 2nd-order-polynomial ground: per cell fit z = a + b*u + c*v + d*u^2
        + e*v^2 + f*uv to the ground points in the 3x3 (15 m) window (u,v = offset from
        the cell centre, in cell units), and read the CONSTANT term a = surface value AT
        the cell centre. Curvature-UNBIASED, unlike the per-cell median (which carries
        the cell's curvature) or a plane (which has no curvature term). Robust over the
        window (gen2 ~60 pts / 3x3, gen1 ~115); falls back to the median where the window
        is too sparse or the normal matrix is singular. Windowed moments are accumulated
        via 9 shifts: a point contributes to each target cell with offset (u0-dj, v0-di)."""
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        ixf, iyf, xf, yf, zf2 = ix[ok], iy[ok], x[ok], y[ok], z[ok]
        u0 = (xf - (X0 + (ixf + 0.5) * res)) / res
        v0 = (yf - (Y0 + (iyf + 0.5) * res)) / res
        Nc = nx * ny; pairs = [(a, b) for a in range(6) for b in range(a, 6)]
        M = [np.zeros(Nc) for _ in range(21)]; Rr = [np.zeros(Nc) for _ in range(6)]
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ti = iyf + di; tj = ixf + dj
                mm = (ti >= 0) & (ti < ny) & (tj >= 0) & (tj < nx)
                t = ti[mm] * nx + tj[mm]; u = u0[mm] - dj; v = v0[mm] - di; zz = zf2[mm]
                phi = [np.ones_like(u), u, v, u * u, v * v, u * v]
                for k, (a, b) in enumerate(pairs): M[k] += np.bincount(t, phi[a] * phi[b], Nc)
                for k in range(6): Rr[k] += np.bincount(t, phi[k] * zz, Nc)
        a0 = np.full(Nc, np.nan); idx = np.where(M[0] >= minpts)[0]
        if len(idx):
            Mm = np.zeros((len(idx), 6, 6))
            for k, (a, b) in enumerate(pairs): Mm[:, a, b] = M[k][idx]; Mm[:, b, a] = M[k][idx]
            rhs = np.stack([Rr[k][idx] for k in range(6)], 1)
            good = np.abs(np.linalg.det(Mm)) > 1e-6
            if good.any(): a0[idx[good]] = np.linalg.solve(Mm[good], rhs[good])[:, 0]
        med = cellstat(x, y, z, "ground").ravel()                # fallback where sparse/singular
        a0[~np.isfinite(a0)] = med[~np.isfinite(a0)]
        return a0.reshape(ny, nx)

    def groundg(x, y, z):
        if ground == "poly2":
            return _poly2_ground(x, y, z)
        if ground == "plane":
            return _plane_ground(x, y, z)
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
    elif ground in ("plane", "poly2"):
        if stream:
            raise NotImplementedError(f"ground='{ground}' not supported with stream=True")
        Zref = groundg(A["x"], A["y"], A["z"])
    else:
        Zref = Z21

    # --- before: (CSF ground classification) -> align -> tie -> drift ---
    # ground_source="csf" (default) runs PDAL CSF on the before cloud first for a
    # cleaner, more general bare-earth (removes structures/understory); "last_return"
    # skips it and uses the raw last-return heuristic. CSF is slow (min/tile).
    # CSF ground is classified ONCE and kept: it is deterministic and slow (~min/tile),
    # and every downstream step (move/tilt/drift/grid) leaves it unchanged. With
    # csf_cache set, the raw classified cloud is archived there and reused on later runs
    # (move + tilt are applied to it below); without it, a temp is used and deleted.
    import os as _os, shutil as _sh
    _csf_tmp = None
    if ground_source == "csf":
        if csf_cache and _os.path.exists(csf_cache):
            before_laz = csf_cache                       # reuse the archived raw CSF
        else:
            _tmp = classify_ground_csf(before_laz, pdal=csf_pdal)
            if csf_cache:                                # archive the raw CSF, keep it
                _os.makedirs(_os.path.dirname(csf_cache) or ".", exist_ok=True)
                _sh.move(_tmp, csf_cache)
                _sh.rmtree(_os.path.dirname(_tmp), ignore_errors=True)
                before_laz = csf_cache
            else:
                _csf_tmp = _tmp; before_laz = _tmp
    f = laspy.read(str(before_laz))
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); gt8 = np.asarray(f.gps_time)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    cl8 = np.asarray(f.classification)
    try:                                             # LAS 1.4 scan angle is stored in 0.006 deg
        sa8 = np.asarray(f.scan_angle).astype(float) * 0.006
    except Exception:                                # older formats / missing -> no boresight term
        sa8 = np.zeros_like(z8)
    # before-epoch ground POINT selection. "csf"/"last_return" -> last returns of the
    # (CSF-classified or raw) cloud; "class2" -> the before survey's OWN ASPRS ground
    # class (a test path: gen1's 2008 vendor classification, which the CSF default was
    # chosen to replace -- see read_after_ground note).
    be = (cl8 == 2) if ground_source == "class2" else (rn8 == nr8)
    if _csf_tmp is not None:
        import shutil, os
        shutil.rmtree(os.path.dirname(_csf_tmp), ignore_errors=True)
    pc = io.PointCloud(x8, y8, z8, ps8, cl8,
                       np.zeros_like(z8), sa8, before_crs)
    # INSTRUMENTAL boresight roll (opt-in), removed per point BEFORE the empirical swath
    # alignment. Calibrated from gen1 self-overlap (gen2-free), so it is decoupled from the
    # gen1-vs-gen2 lateral tie and needs no iteration.
    boresight_used = None
    if correct_boresight:
        boresight_used = (boresight_roll_mm_per_deg if boresight_roll_mm_per_deg is not None
                          else coreg.estimate_boresight_roll(pc, res).b)
        z8 = z8 - boresight_used * sa8 / 1000.0      # mm/deg * deg -> mm -> m
        pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8), sa8, before_crs)
    # ZERO LINE. align_swaths solves a FREE NETWORK and only then subtracts the reference
    # swath's value, so this choice does not touch any swath-to-swath DIFFERENCE -- but it
    # DOES set the absolute level the whole mosaic inherits, because that level becomes the
    # reference line's own error. Measured on elbaext: the six per-swath dz span 44.60 mm
    # (-22.60 .. +22.00), so re-gauging on a different line moves every elevation by up to
    # that much. The product's absolute level therefore DEPENDS ON THE ZERO LINE until a
    # ground-control datum constant is applied; see `zero_line` and
    # `absolute_datum_mm` in corrections.json, and ground_control/apply_datum.py.
    # The ZERO LINE: the flight line defined as zero when this tile's swath network
    # is solved. Arbitrary and per-tile; it sets only the level the tile inherits,
    # and an absolute datum cancels it exactly. Recorded so two products can be
    # related: see scripts/backfill_zero_line.py.
    zero_line = int(ps8.min())
    corr, _, _ = coreg.align_swaths(pc, ref=zero_line, tie=swath_tie)
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz

    # cross-epoch vertical datum, in the principled order: get x,y right, then z.
    # (1) HORIZONTAL: one constant Nuth & Kaeaeb lateral shift from the full topography
    #     (order-0 tie) so the two DEMs are registered in x,y before z is touched.
    #     Drainage divides / ridgelines do not move to first order, so registering the
    #     whole DEM recovers the lateral shift; the aspect-DIPOLE it fits is EROSION-
    #     robust (diffuse erosion has no aspect dependence). (2) VERTICAL: the geoid-model
    #     datum shift N_gen1 - N_gen2 (e.g. GEOID03 - GEOID18) -- a REQUIRED geodetic
    #     offset, auto-computed from the PROJ geoid grids if not supplied (no hard-coded
    #     constant, no arbitrary plane fitted to "stable" surfaces). Residual offsets are
    #     left for later analysis, not baked into the datum.
    from . import references
    if tie != "reference":
        raise ValueError(
            f"tie={tie!r} is not supported. The only cross-epoch datum is the geoid "
            "difference applied after the lateral shift; the reference_plane fit and the "
            "parabola tie were removed (see git history if ever needed).")
    # ORDER 0 = a single CONSTANT (dx, dy) lateral shift (Nuth & Kaeaeb order-0 tie), NOT the
    # removed order-2 parabola. tie_polynomial is a general polynomial fit; only order 0 is used
    # here, and only its horizontal shift (hs["a"][0], hs["b"][0]) is applied to xc, yc.
    hs = coreg.tie_polynomial(Zref, groundg(xc[be], yc[be], zc[be]),
                              res, X0, Y0, order=0)
    xc += coreg.eval_poly_field(hs["a"], xc, yc, hs["norm"], 0)
    yc += coreg.eval_poly_field(hs["b"], xc, yc, hs["norm"], 0)
    if geoid_datum is None:                          # auto-compute from the geoid grids
        geoid_datum = references.geoid_difference(bounds, 26915)
    gc, gb, gcc = geoid_datum        # (const_m, b East, c North) m,m/km of (N_gen1 - N_gen2), ADD to gen1
    cxg = 0.5*(bounds[0]+bounds[2]); cyg = 0.5*(bounds[1]+bounds[3])
    zc += gc + gb*(xc-cxg)/1000.0 + gcc*(yc-cyg)/1000.0
    print(f"  geoid-difference datum: const {1000*gc:+.1f} mm, tilt "
          f"({1000*gb:+.3f},{1000*gcc:+.3f}) mm/km; lateral shift "
          f"({100*hs['a'][0]:+.1f},{100*hs['b'][0]:+.1f}) cm", flush=True)
    tie_info = {"method": "geoid_difference", "const_m": gc, "tilt_b_m_per_km": gb,
                "tilt_c_m_per_km": gcc, "centroid": [cxg, cyg],
                "horizontal_shift_m": [round(float(hs["a"][0]),4), round(float(hs["b"][0]),4)]}

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

    # ABSOLUTE DATUM. Placing each epoch on surveyed NAVD88 removes the zero-line dependence
    # exactly: re-gauging by d shifts z by +d and the measured constant by -d, so the
    # corrected surface is unchanged. Applied to BOTH epochs, so the DoD moves by the
    # DIFFERENCE of the two constants and true change on stable ground goes to zero.
    datum_applied = None
    if absolute_datum is not None:
        need = {"gen1_mm", "gen2_mm", "zero_line", "source"}
        missing = need - set(absolute_datum)
        if missing:
            raise ValueError(f"absolute_datum is missing {sorted(missing)}")
        if int(absolute_datum["zero_line"]) != zero_line:
            raise ValueError(
                f"absolute_datum was measured against zero line "
                f"{absolute_datum['zero_line']} but this run's zero line is "
                f"{zero_line}. A constant belongs to the product it was measured "
                f"against; re-express it with ground_control apply_datum.on_zero_line().")
        g1 = float(absolute_datum["gen1_mm"]); g2 = float(absolute_datum["gen2_mm"])
        Z21 = Z21 + g2 / 1000.0                      # gen2 onto NAVD88
        dod = dod + (g2 - g1) / 1000.0               # DoD moves by the DIFFERENCE
        datum_applied = {
            "gen1_mm": round(g1, 3), "gen2_mm": round(g2, 3),
            "dod_shift_mm": round(g2 - g1, 3),
            "gen1_sigma_mm": absolute_datum.get("gen1_sigma_mm"),
            "gen2_sigma_mm": absolute_datum.get("gen2_sigma_mm"),
            "zero_line": zero_line, "source": absolute_datum["source"]}

    corrections = {
        "epochs": "after - before (positive = deposition)",
        "crs": "EPSG:26915", "res_m": res, "ground_percentile": ground_q,
        "ground_estimator": ground, "ground_source": ground_source,
        "bounds": [float(b) for b in bounds], "stable_1sigma_m": round(sigma, 4),
        "robust_stable": robust_stable,
        "stable_clip_fraction": round(float(stable_clip_frac), 4),
        "lod_method": lod_method,
        "boresight_roll_mm_per_deg": (round(float(boresight_used), 3)
                                      if boresight_used is not None else None),
        "swath_tie": swath_tie,
        "zero_line": zero_line,
        "absolute_level_depends_on_zero_line": True,
        "absolute_datum_mm": datum_applied,
        "absolute_datum_note": (
            "The absolute level of this product is the ZERO LINE's own error, not a "
            "measured elevation: using another line as zero shifts every elevation (44.60 "
            "mm across the six lines at elbaext). Apply a ground-control datum constant "
            "measured against THIS zero line to make the result independent of it -- "
            "corrected = z + c, and if the zero line moves by d then z moves by +d and c "
            "by -d, so the corrected surface is unchanged. See ground_control/apply_datum.py "
            "and ground_control/FRAME.md."),
        "per_swath_internal_alignment_dxdydz_m":
            {str(k): [round(float(v), 4) for v in val] for k, val in corr.items()},
        "cross_epoch_datum": tie_info,
        "along_track_drift_gpsTime_to_m":
            {str(p): {"gps_time": [round(t, 3) for t in c[0]],
                      "drift_m": [round(d, 4) for d in c[1]]} for p, c in curves.items()},
    }
    return dict(dod=dod, lod=lod, z_after=Z21, stable=stable_rep,
                corrections=corrections, stable_sigma=sigma,
                bounds=tuple(bounds), res=res, nx=nx, ny=ny)
