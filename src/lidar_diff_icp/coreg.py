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

    # Shift-component uncertainty from the cosine-fit covariance. The eastward
    # increment is a*sin(b) = p2 and the northward is a*cos(b) = p1, so their
    # variances are exactly coef_cov[1,1] (p2) and coef_cov[0,0] (p1).
    dx_sigma = float(np.sqrt(coef_cov[1, 1])) if np.isfinite(coef_cov[1, 1]) else np.nan
    dy_sigma = float(np.sqrt(coef_cov[0, 0])) if np.isfinite(coef_cov[0, 0]) else np.nan
    # dz uncertainty ~ standard error of the median residual
    dz_sigma = float(nmad_after / np.sqrt(max(n_used, 1)))

    return Coreg(dx, dy, dz, dx_sigma, dy_sigma, dz_sigma, n_used,
                 nmad_before, nmad_after, it, converged)


def align_swaths(pc, res: float = 2.0, exclude=(5, 6, 9)):
    """Free-network least-squares alignment of every swath into one frame.

    Runs Nuth & Kaeaeb on each overlapping swath pair, then solves for a
    per-swath 3-D shift (Dx, Dy, Dz) that makes all overlaps mutually
    consistent. The observation for edge (a, b) is ``c_b - c_a = s_ab`` where
    ``s_ab`` aligns b onto a; the system is solved per component with a
    **zero-mean gauge**, so the whole group's absolute offset (e.g. from the
    2021 3DEP) is left free and must be tied separately.

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
        c -= c.mean()                      # zero-mean gauge (group offset free)
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
