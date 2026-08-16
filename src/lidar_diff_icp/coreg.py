"""Nuth & Kaeaeb (2011) horizontal + vertical co-registration.

A horizontal misregistration between two elevation surfaces makes their
difference vary with terrain aspect, scaled by slope:

    dh / tan(slope) = a * cos(b - aspect) + c

where ``a`` is the horizontal-shift magnitude, ``b`` its direction, and ``c``
carries the vertical bias. Fitting the cosine and iterating (resample, refit)
converges to a rigid 3-D translation (dx, dy, dz) that brings ``z_src`` onto
``z_ref``.

Reference: Nuth, C. & Kaeaeb, A. (2011), *The Cryosphere* 5, 271-290.

Conventions: grids are row-major with row index increasing **north** (origin
"lower"), column index increasing **east**; ``res`` is the cell size in metres.
``dx`` is the eastward and ``dy`` the northward shift applied to ``z_src``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import map_coordinates


@dataclass
class Coreg:
    dx: float          # eastward shift of z_src -> z_ref (m)
    dy: float          # northward shift (m)
    dz: float          # vertical shift (m)
    dx_sigma: float    # 1-sigma, from the cosine-fit covariance
    dy_sigma: float
    dz_sigma: float
    n: int             # cells used in the final fit
    nmad_before: float # robust dh scatter before co-registration (m)
    nmad_after: float  # after (the achieved co-registration floor, m)
    n_iter: int
    converged: bool


def slope_aspect(z: np.ndarray, res: float):
    """Return (slope, aspect) in radians. Aspect is clockwise from north."""
    dz_dnorth = np.gradient(z, res, axis=0)  # +row = north
    dz_deast = np.gradient(z, res, axis=1)   # +col = east
    slope = np.arctan(np.hypot(dz_deast, dz_dnorth))
    # aspect: direction of steepest ascent, clockwise from north
    aspect = np.arctan2(dz_deast, dz_dnorth)
    return slope, aspect


def _shift_grid(z: np.ndarray, dx: float, dy: float, res: float) -> np.ndarray:
    """Bilinearly resample ``z`` so its content moves by (dx east, dy north).

    NaNs are propagated: an output cell is valid only if its full bilinear
    stencil was valid (interpolated weight ~ 1).
    """
    ny, nx = z.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    cx = xx - dx / res            # to move content +dx east, sample at x-dx
    cy = yy - dy / res            # +dy north
    valid = np.isfinite(z).astype(float)
    filled = np.where(np.isfinite(z), z, 0.0)
    num = map_coordinates(filled, [cy, cx], order=1, mode="constant", cval=0.0)
    wgt = map_coordinates(valid, [cy, cx], order=1, mode="constant", cval=0.0)
    out = np.where(wgt > 0.999, num / np.where(wgt == 0, 1.0, wgt), np.nan)
    return out


def _nmad(a: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def nuth_kaab(z_ref: np.ndarray, z_src: np.ndarray, res: float,
              slope_min_deg: float = 3.0, max_iter: int = 20,
              tol: float = 0.005) -> Coreg:
    """Co-register ``z_src`` to ``z_ref`` on a shared grid. See module docstring."""
    slope0 = slope_aspect(z_ref, res)[0]
    dh0 = z_ref - z_src
    nmad_before = _nmad(dh0[np.isfinite(dh0)])

    dx = dy = dz = 0.0
    tan_min = np.tan(np.radians(slope_min_deg))
    coef_cov = np.full((3, 3), np.nan)
    sigma = np.nan
    n_used = 0
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        z_shift = _shift_grid(z_src, dx, dy, res) + dz
        slope, aspect = slope_aspect(z_ref, res)
        dh = z_ref - z_shift
        tan_s = np.tan(slope)
        m = np.isfinite(dh) & (tan_s > tan_min)
        if m.sum() < 100:
            break
        # robust outlier rejection on dh
        med, nm = np.median(dh[m]), _nmad(dh[m])
        m &= np.abs(dh - med) < 3 * max(nm, 1e-3)
        y = (dh[m] - np.median(dh[m])) / tan_s[m]
        psi = aspect[m]
        X = np.c_[np.cos(psi), np.sin(psi), np.ones(psi.size)]
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        p1, p2, _ = coef
        a = np.hypot(p1, p2)
        b = np.arctan2(p2, p1)
        # covariance of the linear fit (for shift uncertainty)
        resid = y - X @ coef
        sigma = float(np.sqrt(np.sum(resid**2) / max(psi.size - 3, 1)))
        coef_cov = sigma**2 * np.linalg.inv(X.T @ X)
        n_used = int(psi.size)
        # The fit recovers the existing offset (east = a*sin b, north = a*cos b);
        # the correction is its negative.
        ddx, ddy = a * np.sin(b), a * np.cos(b)
        dx -= ddx
        dy -= ddy
        # update vertical shift as the median residual after horizontal move
        dz = float(np.nanmedian(z_ref - _shift_grid(z_src, dx, dy, res)))
        if np.hypot(ddx, ddy) < tol:
            converged = True
            break

    z_final = _shift_grid(z_src, dx, dy, res) + dz
    dh_final = z_ref - z_final
    nmad_after = _nmad(dh_final[np.isfinite(dh_final)])

    # Divergence guard (as in tie_polynomial): on gentle terrain the horizontal
    # shift (~1/tan slope) can run away and worsen the fit. Fall back to a rigid
    # vertical offset if the horizontal solution is implausible (> 10 m; real
    # airborne shifts are sub-metre) or fails to beat it.
    dz0 = float(np.nanmedian(dh0))
    nmad_rigid = _nmad((dh0 - dz0)[np.isfinite(dh0)])
    if (not np.isfinite(nmad_after)) or nmad_after > nmad_rigid or np.hypot(dx, dy) > 10.0:
        dx = dy = 0.0; dz = dz0; nmad_after = nmad_rigid; converged = True
        coef_cov = np.full((3, 3), np.nan)      # horizontal not estimated

    # Shift-component uncertainty from the cosine-fit covariance. The eastward
    # increment is a*sin(b) = p2 and the northward is a*cos(b) = p1, so their
    # variances are exactly coef_cov[1,1] (p2) and coef_cov[0,0] (p1).
    dx_sigma = float(np.sqrt(coef_cov[1, 1])) if np.isfinite(coef_cov[1, 1]) else np.nan
    dy_sigma = float(np.sqrt(coef_cov[0, 0])) if np.isfinite(coef_cov[0, 0]) else np.nan
    # dz uncertainty ~ standard error of the median residual
    dz_sigma = float(nmad_after / np.sqrt(max(n_used, 1)))

    return Coreg(dx, dy, dz, dx_sigma, dy_sigma, dz_sigma, n_used,
                 nmad_before, nmad_after, it, converged)


def tie_translation_tilt(z_ref, z_src, res, x_origin, y_origin,
                         slope_min_deg=3.0, max_iter=25, tol=0.003):
    """Best overall co-registration of ``z_src`` onto ``z_ref``: horizontal
    translation PLUS a planar vertical tilt, one robust joint fit over the DEM.

    First-order model, fit each iteration by Huber-IRLS with 3-sigma rejection
    of change outliers:

        dh = dx * dz/dE + dy * dz/dN + c0 + c1 * X + c2 * Y

    (dx,dy) is the horizontal shift; (c1,c2) the vertical tilt gradient; c0 the
    vertical offset. Returns a dict of the parameters plus tilt (mm/m), the total
    vertical ramp across the DEM (m), and the robust residual before/after.
    """
    ny, nx = z_ref.shape
    gy_i, gx_i = np.mgrid[0:ny, 0:nx]
    X = x_origin + (gx_i + 0.5) * res
    Y = y_origin + (gy_i + 0.5) * res
    fin0 = np.isfinite(z_ref)
    Xc = X - X[fin0].mean()
    Yc = Y - Y[fin0].mean()
    gx = np.gradient(z_ref, res, axis=1)      # dz/dE
    gy = np.gradient(z_ref, res, axis=0)      # dz/dN
    tanS = np.hypot(gx, gy)
    tan_min = np.tan(np.radians(slope_min_deg))
    nmad_before = _nmad((z_ref - z_src)[np.isfinite(z_ref - z_src)])

    dx = dy = c0 = c1 = c2 = 0.0
    converged = False
    for _ in range(max_iter):
        vert = c0 + c1 * Xc + c2 * Yc
        dh = z_ref - (_shift_grid(z_src, dx, dy, res) + vert)
        m = np.isfinite(dh) & (tanS > tan_min)
        med, nm = np.median(dh[m]), _nmad(dh[m])
        m &= np.abs(dh - med) < 3 * max(nm, 1e-3)
        A = np.c_[gx[m], gy[m], np.ones(m.sum()), Xc[m], Yc[m]]
        y = dh[m]
        w = np.ones(y.size); coef = np.zeros(5)
        for _ in range(5):                    # Huber IRLS
            Aw = A * w[:, None]
            coef, *_ = np.linalg.lstsq(Aw, y * w, rcond=None)
            r = y - A @ coef
            s = 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-9
            k = 1.345 * s; a = np.abs(r)
            w = np.where(a <= k, 1.0, k / a)
        dx -= coef[0]; dy -= coef[1]           # horizontal offset -> correction
        c0 += coef[2]; c1 += coef[3]; c2 += coef[4]
        if np.hypot(coef[0], coef[1]) < tol and abs(coef[2]) < tol:
            converged = True
            break
    dh_final = z_ref - (_shift_grid(z_src, dx, dy, res) + c0 + c1 * Xc + c2 * Yc)
    m = np.isfinite(dh_final)
    return dict(dx=float(dx), dy=float(dy), c0=float(c0), c1=float(c1), c2=float(c2),
                tilt_mm_per_m=float(np.hypot(c1, c2) * 1000),
                ramp_m=float(abs(c1) * (X[fin0].max() - X[fin0].min())
                             + abs(c2) * (Y[fin0].max() - Y[fin0].min())),
                nmad_before=float(nmad_before), nmad_after=float(_nmad(dh_final[m])),
                converged=converged, n=int(m.sum()))


def _poly_basis(xn, yn, order):
    b = [np.ones_like(xn), xn, yn]
    if order >= 2:
        b += [xn * xn, xn * yn, yn * yn]
    return b


def _warp_grid(z, dxf, dyf, res):
    """Resample z with a spatially varying shift (dxf, dyf per cell)."""
    ny, nx = z.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    cx = xx - dxf / res
    cy = yy - dyf / res
    valid = np.isfinite(z).astype(float)
    filled = np.where(np.isfinite(z), z, 0.0)
    num = map_coordinates(filled, [cy, cx], order=1, mode="constant", cval=0.0)
    wgt = map_coordinates(valid, [cy, cx], order=1, mode="constant", cval=0.0)
    return np.where(wgt > 0.999, num / np.where(wgt == 0, 1.0, wgt), np.nan)


def eval_poly_field(coef, x, y, norm, order):
    """Evaluate a fitted polynomial correction field at arbitrary (x, y)."""
    xm, xhr, ym, yhr = norm
    xn = (x - xm) / xhr; yn = (y - ym) / yhr
    basis = _poly_basis(xn, yn, order)
    return sum(coef[k] * basis[k] for k in range(len(basis)))


def tie_polynomial(z_ref, z_src, res, x_origin, y_origin, order=2,
                   slope_min_deg=3.0, max_iter=30, tol=0.003):
    """Spatially varying co-registration: dx, dy, dz each an order-`order`
    polynomial in (x, y), fit jointly by robust IRLS with change rejection.

    Removes a smooth (kilometre-scale) warp between the two surfaces while being
    too low-order to absorb finer geomorphic change. Returns the coefficient
    vectors (a=dx, b=dy, c=dz), the correction fields on the grid, the
    normalization for re-evaluating the fields at points, and residuals.
    """
    ny, nx = z_ref.shape
    gy_i, gx_i = np.mgrid[0:ny, 0:nx]
    X = x_origin + (gx_i + 0.5) * res
    Y = y_origin + (gy_i + 0.5) * res
    fin = np.isfinite(z_ref)
    xm = 0.5 * (X[fin].max() + X[fin].min()); xhr = 0.5 * (X[fin].max() - X[fin].min())
    ym = 0.5 * (Y[fin].max() + Y[fin].min()); yhr = 0.5 * (Y[fin].max() - Y[fin].min())
    Xn = (X - xm) / xhr; Yn = (Y - ym) / yhr
    basis = _poly_basis(Xn, Yn, order); K = len(basis)
    gx = np.gradient(z_ref, res, axis=1)
    gy = np.gradient(z_ref, res, axis=0)
    tanS = np.hypot(gx, gy); tan_min = np.tan(np.radians(slope_min_deg))
    nmad_before = _nmad((z_ref - z_src)[np.isfinite(z_ref - z_src)])

    a = np.zeros(K); b = np.zeros(K); c = np.zeros(K)
    fld = lambda coef: sum(coef[k] * basis[k] for k in range(K))
    converged = False
    for _ in range(max_iter):
        dh = z_ref - (_warp_grid(z_src, fld(a), fld(b), res) + fld(c))
        m = np.isfinite(dh) & (tanS > tan_min)
        if m.sum() < 3 * K:           # too few sloped cells to constrain the fit
            break
        med, nm = np.median(dh[m]), _nmad(dh[m])
        m &= np.abs(dh - med) < 3 * max(nm, 1e-3)
        cols = ([gx[m] * basis[k][m] for k in range(K)]
                + [gy[m] * basis[k][m] for k in range(K)]
                + [basis[k][m] for k in range(K)])
        A = np.column_stack(cols); y = dh[m]
        w = np.ones(y.size); coef = np.zeros(3 * K)
        for _ in range(5):
            coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
            r = y - A @ coef
            s = 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-9
            k = 1.345 * s; aa = np.abs(r); w = np.where(aa <= k, 1.0, k / aa)
        a -= coef[:K]; b -= coef[K:2 * K]; c += coef[2 * K:]
        if np.hypot(coef[0], coef[K]) < tol and abs(coef[2 * K]) < tol:
            converged = True
            break
    dxf, dyf, dzf = fld(a), fld(b), fld(c)
    dh_poly = z_ref - (_warp_grid(z_src, dxf, dyf, res) + dzf)
    mf = np.isfinite(dh_poly)
    nmad_poly = _nmad(dh_poly[mf]) if mf.any() else np.inf
    # Divergence guard. On gentle terrain few cells exceed slope_min, so the
    # Nuth & Kaeaeb horizontal shift (~1/tan(slope)) is unconstrained and can run
    # away (dx ~ 1000 m). Accept the polynomial tie ONLY if it beats a plain
    # rigid vertical offset; otherwise fall back to that offset -- so the tie is
    # never worse than doing nothing but removing a bias.
    diff = z_ref - z_src; md = np.isfinite(diff)
    dz0 = float(np.median(diff[md]))
    nmad_rigid = _nmad((diff - dz0)[md])
    # accept the polynomial tie only if it beats the rigid offset AND its
    # horizontal shift is physically plausible (< 10 m; real airborne shifts are
    # sub-metre) -- an over-fit tie can improve noise with an absurd shift.
    max_shift = max(np.ptp(dxf), np.ptp(dyf)) if mf.any() else np.inf
    if np.isfinite(nmad_poly) and nmad_poly <= nmad_rigid and max_shift < 10.0:
        fallback = False; nmad_after = nmad_poly
    else:
        fallback = True; converged = True
        a = np.zeros(K); b = np.zeros(K); c = np.zeros(K); c[0] = dz0
        dxf = fld(a); dyf = fld(b); dzf = fld(c); nmad_after = nmad_rigid
    return dict(a=a, b=b, c=c, order=order, norm=(xm, xhr, ym, yhr),
                dx_field=dxf, dy_field=dyf, dz_field=dzf,
                nmad_before=nmad_before, nmad_after=nmad_after,
                converged=converged, fallback=fallback, n=int(mf.sum()))


def correction_surface(z_ref, z_src, res, x_origin, y_origin, *,
                       slope_thresh_deg=3.0, dz_thresh=0.7, radius=400.0,
                       power=2.0, source_step=4, k=32, exclude=None):
    """Vertical error-correction surface after DeLong et al. (2022, ESS).

    From two gridded ground surfaces, form the raw vertical difference
    ``dz = z_ref - z_src``, mask cells likely to hold REAL change -- local slope
    > ``slope_thresh_deg`` or ``|dz| > dz_thresh`` -- and inverse-distance
    interpolate (power ``power``, search ``radius`` m) the remaining STABLE
    residual across every cell. The result is a smooth field of spurious
    vertical offset: add it to ``z_src`` (or subtract from ``dz``) to remove the
    error while leaving change that exceeds the masks. Unlike a low-order
    polynomial tie this is nonparametric, so it can follow a multi-lobe warp.

    Precondition (from the paper): valid only where stable ground is
    widespread. The paper also buffers streams / valley floors, since real
    low-slope deposition there can be absorbed. Pass such a mask as ``exclude``
    (a bool array, True = drop from the stable sources) -- e.g. a floodplain
    mask from a topographic position index, which flow accumulation cannot
    reliably give in a flat valley. Without it, ``dz_thresh`` is the only guard.

    ``source_step`` thins the stable cells used as IDW sources (every Nth cell
    each axis) for speed; the surface is smooth at ``radius`` so this is benign.

    Cells with NO stable source within ``radius`` are left ``NaN`` (uncorrected)
    rather than extrapolated: block cross-validation shows the interpolation has
    skill only across gaps up to ~the correlation length (here it removes
    23-41% of stable residual variance for gaps <=300 m, ~0 by 600 m, and
    INJECTS error beyond ~1 km). ``radius`` is therefore the honest reach of the
    correction, not just a search cutoff. The returned ``gap`` (distance to the
    nearest stable source) lets a caller see where the correction is trustworthy.
    Returns dict(C, gap, stable, dz, n_stable, frac_stable).
    """
    from scipy.spatial import cKDTree
    ny, nx = z_ref.shape
    dz = z_ref - z_src
    slope, _ = slope_aspect(z_ref, res)
    stable = (np.isfinite(dz) & (np.degrees(slope) < slope_thresh_deg)
              & (np.abs(dz) < dz_thresh))
    if exclude is not None:
        stable = stable & ~np.asarray(exclude, bool)   # e.g. floodplain buffer
    jj, ii = np.mgrid[0:ny, 0:nx]
    X = x_origin + (ii + 0.5) * res
    Y = y_origin + (jj + 0.5) * res
    sel = stable & (ii % source_step == 0) & (jj % source_step == 0)
    if sel.sum() < k:
        raise ValueError(f"too few stable source cells ({sel.sum()}) for k={k}")
    tree = cKDTree(np.c_[X[sel], Y[sel]])
    sdz = dz[sel]
    dist, idx = tree.query(np.c_[X.ravel(), Y.ravel()], k=k)
    w = 1.0 / np.maximum(dist, 1e-6) ** power
    w[dist > radius] = 0.0
    wsum = w.sum(1)
    C = np.where(wsum > 0, (w * sdz[idx]).sum(1) / np.where(wsum == 0, 1.0, wsum),
                 np.nan)                          # no stable within radius -> NaN
    gap = dist[:, 0]                              # distance to nearest stable source
    return dict(C=C.reshape(ny, nx), gap=gap.reshape(ny, nx), stable=stable,
                dz=dz, n_stable=int(stable.sum()), frac_stable=float(stable.mean()))


def fit_along_track_drift(gps_time, change_on_stable, is_stable, swath, *,
                          n_bins=120, s_frac=1.0, min_pts=2000):
    """Per-swath along-track (GNSS trajectory-drift) vertical correction.

    The dominant residual in early (2008) Minnesota lidar after per-swath
    translation and a smooth cross-epoch tie is a DETERMINISTIC along-track
    drift: a smooth vertical undulation as a function of ``gps_time`` within each
    flight line (``point_source_id``), ~constant across-track, no roll/scan-angle
    term. This models it directly, so it does not overfit terrain-correlated
    noise the way a spatial interpolator on the difference does.

    Per swath, the stable-ground ``change_on_stable`` (reference - 2008, ~0 where
    stable) is binned by ``gps_time``, robust-median per bin, smoothed, and
    evaluated at every point of the swath. Add the returned drift to the 2008
    elevations to correct them (new change = old change - drift).

    ``s_frac`` sets the smoothing spline stiffness (cross-validation on this data
    favours heavy smoothing, ~1; the drift is long-wavelength). Robust binned
    medians feed a cubic ``UnivariateSpline`` in ``gps_time`` -- a continuous
    trajectory-bias model, the standard strip-adjustment form.

    Returns (drift_per_point, curves) with curves[swath] = (gps_bin_centers,
    drift_curve). Universal form; only the per-tile coefficients differ.
    """
    from scipy.interpolate import UnivariateSpline
    gps_time = np.asarray(gps_time); change_on_stable = np.asarray(change_on_stable)
    drift = np.zeros(len(gps_time)); curves = {}
    for p in np.unique(swath):
        sw = swath == p
        s = sw & is_stable & np.isfinite(change_on_stable)
        if s.sum() < min_pts:
            continue
        t = gps_time[s]; c = change_on_stable[s]
        edges = np.linspace(t.min(), t.max(), n_bins + 1)
        cen = 0.5 * (edges[:-1] + edges[1:])
        bi = np.clip(np.digitize(t, edges) - 1, 0, n_bins - 1)
        prof = np.array([np.median(c[bi == k]) if (bi == k).sum() >= 5 else np.nan
                         for k in range(n_bins)])
        ok = np.isfinite(prof)
        if ok.sum() < 10:
            continue
        spl = UnivariateSpline(cen[ok], prof[ok], k=3,
                               s=s_frac * ok.sum() * np.var(prof[ok]))
        drift[sw] = spl(gps_time[sw])
        curves[int(p)] = (cen, spl(cen))
    return drift, curves


def align_swaths(pc, res: float = 2.0, exclude=(5, 6, 9), ref=None):
    """Free-network least-squares alignment of every swath into one frame.

    Runs Nuth & Kaeaeb on each overlapping swath pair, then solves for a
    per-swath 3-D shift (Dx, Dy, Dz) that makes all overlaps mutually
    consistent. The observation for edge (a, b) is ``c_b - c_a = s_ab`` where
    ``s_ab`` aligns b onto a.

    Gauge (choice of datum; does not change the *relative* solution):
    ``ref=None`` -> zero-mean (group offset spread evenly, absolute frame free);
    ``ref=<swath id>`` -> that swath is pinned to zero and becomes the local
    reference (all others measured relative to it). Either way the group's
    absolute offset from another epoch must be tied separately.

    Returns ``(corrections, edges, misclosure)`` where ``corrections`` maps
    swath id -> (Dx, Dy, Dz) m, ``edges`` lists the pairwise observations, and
    ``misclosure`` is the per-edge residual (~0 for a tree/chain; nonzero only
    where redundant overlaps -- loops -- exist to check consistency).
    """
    from itertools import combinations
    swaths = pc.swaths.tolist()
    idx = {s: k for k, s in enumerate(swaths)}
    edges = []
    for a, b in combinations(swaths, 2):
        try:
            c = coregister_swaths(pc, a, b, res, exclude)
        except ValueError:
            continue
        edges.append((a, b, c.dx, c.dy, c.dz, float(c.n)))
    if not edges:
        raise ValueError("no overlapping swath pairs")
    n, E = len(swaths), len(edges)
    A = np.zeros((E, n)); w = np.zeros(E); O = np.zeros((E, 3))
    for e, (a, b, dx, dy, dz, ww) in enumerate(edges):
        A[e, idx[a]] = -1.0; A[e, idx[b]] = 1.0
        w[e] = ww; O[e] = (dx, dy, dz)
    sw = np.sqrt(w)
    corr = np.zeros((n, 3)); mis = np.zeros((E, 3))
    for k in range(3):
        c, *_ = np.linalg.lstsq(A * sw[:, None], O[:, k] * sw, rcond=None)
        c -= c[idx[ref]] if ref is not None else c.mean()   # gauge choice
        corr[:, k] = c
        mis[:, k] = A @ c - O[:, k]
    corrections = {swaths[i]: tuple(float(v) for v in corr[i]) for i in range(n)}
    return corrections, edges, mis


def apply_alignment(pc, corrections):
    """Return copies of (x, y, z) with each swath shifted by its correction."""
    x = pc.x.copy(); y = pc.y.copy(); z = pc.z.copy()
    for s, (dx, dy, dz) in corrections.items():
        m = pc.point_source_id == s
        x[m] += dx; y[m] += dy; z[m] += dz
    return x, y, z


def coregister_swaths(pc, swath_ref: int, swath_src: int, res: float = 2.0,
                      exclude=(5, 6, 9)) -> Coreg:
    """Nuth & Kaeaeb co-registration of ``swath_src`` onto ``swath_ref`` over
    their overlap, using density-robust per-cell median-Z surfaces."""
    from .swathdiff import _median_grid
    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == swath_ref)
    mb = terr & (pc.point_source_id == swath_src)
    x, y, z = pc.x, pc.y, pc.z
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"swaths {swath_ref} and {swath_src} do not overlap")
    nx = int(np.ceil((x1 - x0) / res)); ny = int(np.ceil((y1 - y0) / res))
    z_ref = _median_grid(x[ma], y[ma], z[ma], res, x0, y0, nx, ny)
    z_src = _median_grid(x[mb], y[mb], z[mb], res, x0, y0, nx, ny)
    return nuth_kaab(z_ref, z_src, res)
