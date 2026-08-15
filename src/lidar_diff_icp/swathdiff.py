"""Inter-swath (flight-line) difference surfaces — the self-calibration check.

Between two swaths flown minutes-to-hours apart there is no real land-surface
change, so their overlap difference is pure acquisition/navigation error. We
estimate a **density-robust** bare-earth surface per swath (per-cell median Z of
terrain returns, which — unlike per-cell *minimum* Z — is not biased by the two
swaths' differing point densities) and difference them on the shared grid.

Reported per pair:
- ``median_offset`` : robust vertical bias (m); the honest offset estimate.
- ``robust_std``    : normalized-MAD scatter (m); roughness + slope-aliased
                      horizontal error + residual vegetation.
- ``tilt``          : magnitude of a robustly fit planar gradient (mm/m).

The mean and an ordinary-least-squares tilt are intentionally avoided: the cell
difference distribution is right-skewed by canopy returns, which inflates both.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .io import PointCloud

# ASPRS classes excluded from the bare-earth proxy: high veg, building, water.
DEFAULT_EXCLUDE = (5, 6, 9)


@dataclass
class SwathDiff:
    swath_a: int
    swath_b: int
    n_cells: int
    median_offset: float
    robust_std: float
    tilt: float          # mm / m
    diff: np.ndarray     # 2-D difference grid (a - b), NaN where unpaired
    extent: tuple[float, float, float, float]  # (x0, x1, y0, y1)


def _median_grid(x, y, z, res, x0, y0, nx, ny) -> np.ndarray:
    """Per-cell median Z on a regular grid; NaN in empty cells."""
    ix = ((x - x0) / res).astype(int)
    iy = ((y - y0) / res).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    flat = iy[ok] * nx + ix[ok]
    med = pd.Series(z[ok]).groupby(flat).median()
    out = np.full(nx * ny, np.nan)
    out[med.index.values] = med.values
    return out.reshape(ny, nx)


def _robust_plane_fit(gx, gy, d, iters=5):
    """Huber-IRLS plane fit d ~ a*gx + b*gy + c. Returns (a, b, c)."""
    A = np.c_[gx, gy, np.ones_like(gx)]
    w = np.ones_like(d)
    coef = np.zeros(3)
    for _ in range(iters):
        Aw = A * w[:, None]
        coef, *_ = np.linalg.lstsq(Aw, d * w, rcond=None)
        r = d - A @ coef
        s = 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-9
        k = 1.345 * s
        a = np.abs(r)
        w = np.where(a <= k, 1.0, k / a)
    return coef


def swath_difference(pc: PointCloud, swath_a: int, swath_b: int,
                     res: float = 2.0,
                     exclude=DEFAULT_EXCLUDE) -> SwathDiff:
    """Density-robust bare-earth difference of two swaths over their overlap."""
    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == swath_a)
    mb = terr & (pc.point_source_id == swath_b)
    x, y, z = pc.x, pc.y, pc.z

    x0 = max(x[ma].min(), x[mb].min())
    x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min())
    y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"swaths {swath_a} and {swath_b} do not overlap")
    nx = int(np.ceil((x1 - x0) / res))
    ny = int(np.ceil((y1 - y0) / res))

    za = _median_grid(x[ma], y[ma], z[ma], res, x0, y0, nx, ny)
    zb = _median_grid(x[mb], y[mb], z[mb], res, x0, y0, nx, ny)
    d = za - zb
    m = np.isfinite(d)
    dd = d[m]

    med = float(np.median(dd))
    rstd = float(1.4826 * np.median(np.abs(dd - med)))
    gy, gx = np.mgrid[0:ny, 0:nx]
    a, b, _ = _robust_plane_fit(gx[m] * res, gy[m] * res, dd)
    tilt = float(np.hypot(a, b) * 1000.0)

    return SwathDiff(swath_a, swath_b, int(dd.size), med, rstd, tilt,
                     d, (x0, x1, y0, y1))
