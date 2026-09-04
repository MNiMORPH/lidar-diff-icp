"""Per-cell ground estimators: how a cloud becomes one elevation per grid cell.

FOUR ESTIMATORS, ONE GRID. These were nested inside ``difference_dem``, closing over its
locals, so nothing could call them and nothing could test them in isolation -- and the
function was 643 lines partly because they lived in it. Moved out unchanged; the closure
over ``(X0, Y0, res, nx, ny)`` is now an explicit ``grid`` argument, and the ground
quantile, which used to be picked up from the enclosing ``_GQ_SCALAR``, is passed in.

  cellstat            per-cell quantile / robust spread / count
  plane_ground        least-squares plane per cell, read at the cell centre
  poly2_ground        2nd-order polynomial over a 3x3 window, read at the centre
  slope_normal_ground quantile of the residual to a common regional plane
  estimate_ground     the dispatcher over the four

WHICH QUANTILE. ``slope_normal_ground`` takes the plain scalar, never a calibrated per-cell
percentile: it grids BOTH epochs, and the ground-q calibration is gen2-only -- the gen1
curve came back flat at 0.219 with q = 0.50 biased +89.8 mm, which is a statement about the
2008 control's uncertain vertical reference rather than about the lidar
(analysis/GROUND_Q_FROM_CLASS2_SPREAD.md).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["cellstat", "plane_ground", "poly2_ground", "slope_normal_ground",
           "estimate_ground"]


def cellstat(x, y, z, how, grid, q=0.50):
    X0, Y0, res, nx, ny = grid
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

def plane_ground(x, y, z, grid, q=0.50, minpts=4):
    X0, Y0, res, nx, ny = grid
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
    med = cellstat(x, y, z, "ground", grid, q).ravel()                # fallback where too sparse
    c[~valid] = med[~valid]
    return c.reshape(ny, nx)

def poly2_ground(x, y, z, grid, q=0.50, minpts=18):
    X0, Y0, res, nx, ny = grid
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
    med = cellstat(x, y, z, "ground", grid, q).ravel()                # fallback where sparse/singular
    a0[~np.isfinite(a0)] = med[~np.isfinite(a0)]
    return a0.reshape(ny, nx)

def slope_normal_ground(x, y, z, grid, q, plane):
    """Quantile of the residual to a common smoothed regional plane, plus the plane
    back. ``plane`` is ``(Zreg_f, dzde, dzdn)``, flat per-cell arrays.

    Takes the PLAIN scalar quantile, never a calibrated per-cell percentile: this
    grids BOTH epochs and the ground-q calibration is gen2-only.
    """
    X0, Y0, res, nx, ny = grid
    Zreg_f, dzde, dzdn = plane
    ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    f = iy[ok] * nx + ix[ok]
    dxe = x[ok] - (X0 + (ix[ok] + 0.5) * res)
    dyn = y[ok] - (Y0 + (iy[ok] + 0.5) * res)
    resid = z[ok] - (Zreg_f[f] + dxe * dzde[f] + dyn * dzdn[f])
    # _GQ_SCALAR, never the calibrated per-cell q. groundg grids BOTH epochs -- gen1
    # via tie_polynomial below, gen2 on the non-stream path -- and the calibration is
    # GEN2-ONLY: it was fitted on gen2 control, and the gen1 curve came back flat at
    # 0.219 with q = 0.50 biased +89.8 mm, which is a statement about the 2008 control's
    # uncertain vertical reference rather than about the lidar
    # (analysis/GROUND_Q_FROM_CLASS2_SPREAD.md). gen1 keeps the plain median.
    s = pd.Series(resid).groupby(f).quantile(q)
    out = np.full(nx * ny, np.nan)
    out[s.index.values] = Zreg_f[s.index.values] + s.values
    return out.reshape(ny, nx)

def estimate_ground(x, y, z, grid, q, ground, plane=None):
    """Dispatch to the estimator named by ``ground``."""
    if ground == "poly2":
        return poly2_ground(x, y, z, grid, q)
    if ground == "plane":
        return plane_ground(x, y, z, grid, q)
    if ground != "slope_normal":
        return cellstat(x, y, z, "ground", grid, q)
    if plane is None:
        raise ValueError("ground='slope_normal' needs plane=(Zreg_f, dzde, dzdn)")
    return slope_normal_ground(x, y, z, grid, q, plane)
