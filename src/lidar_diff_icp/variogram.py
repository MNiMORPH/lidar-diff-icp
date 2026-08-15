"""Empirical (semi)variograms and correlated-error detection limits.

DEM/lidar errors are spatially correlated, so the uncertainty of a mean
elevation difference over a feature does **not** fall as 1/sqrt(N_pixels). It is
governed by the error correlation length. Here we estimate that length from an
empirical variogram of co-registration residuals and propagate it to a
detection limit following the correlated-error logic of Rolstad et al. (2009)
and Hugonnet et al. (2022): white (nugget) variance averages over all pixels,
correlated (sill) variance averages only over N_eff ~ area / (pi * range^2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit


def empirical_variogram(x, y, v, max_lag, n_lags=25, n_pairs=800_000,
                        direction=None, tol_deg=25.0, estimator="dowd", seed=0):
    """Random-pair empirical semivariogram.

    ``estimator``: ``"matheron"`` (classical mean of 0.5*dv^2, outlier-sensitive)
    or ``"dowd"`` (robust: gamma = 1.099 * median(|dv|)^2). Lidar difference
    fields are heavy-tailed (vegetation), so ``"dowd"`` is the default.

    If ``direction`` (degrees clockwise from north) is given, only pairs whose
    separation vector lies within ``tol_deg`` of that axis are used, giving a
    directional variogram (to expose flight-line anisotropy).
    Returns (lag_centers, gamma, counts).
    """
    rng = np.random.default_rng(seed)
    n = v.size
    i = rng.integers(0, n, n_pairs)
    j = rng.integers(0, n, n_pairs)
    dx = x[j] - x[i]
    dy = y[j] - y[i]
    h = np.hypot(dx, dy)
    keep = (h > 0) & (h <= max_lag)
    if direction is not None:
        ang = np.degrees(np.arctan2(dx, dy)) % 180.0   # 0 = N, 90 = E
        da = np.abs((ang - (direction % 180.0) + 90.0) % 180.0 - 90.0)
        keep &= da < tol_deg
    h = h[keep]
    adv = np.abs((v[j] - v[i])[keep])
    edges = np.linspace(0, max_lag, n_lags + 1)
    idx = np.clip(np.digitize(h, edges) - 1, 0, n_lags - 1)
    counts = np.bincount(idx, minlength=n_lags)
    centers = 0.5 * (edges[:-1] + edges[1:])
    gamma = np.full(n_lags, np.nan)
    for k in range(n_lags):
        sel = idx == k
        if counts[k] == 0:
            continue
        if estimator == "matheron":
            gamma[k] = 0.5 * np.mean(adv[sel] ** 2)
        elif estimator == "dowd":
            gamma[k] = 1.099 * np.median(adv[sel]) ** 2
        else:
            raise ValueError(f"unknown estimator {estimator!r}")
    return centers, gamma, counts


def _spherical(h, nugget, sill, rng):
    hr = np.clip(h / rng, 0, 1)
    return nugget + sill * (1.5 * hr - 0.5 * hr ** 3)


@dataclass
class VariogramModel:
    nugget: float
    sill: float          # partial sill (correlated variance)
    range_: float        # correlation length (m)

    @property
    def total_sill(self) -> float:
        return self.nugget + self.sill


def fit_spherical(centers, gamma, counts) -> VariogramModel:
    m = np.isfinite(gamma) & (counts > 0)
    h, g, w = centers[m], gamma[m], counts[m]
    s0 = np.nanmax(g)
    p0 = [0.2 * s0, 0.8 * s0, 0.3 * h.max()]
    bounds = ([0, 0, h.min()], [s0, 2 * s0, h.max()])
    popt, _ = curve_fit(_spherical, h, g, p0=p0, bounds=bounds,
                        sigma=1.0 / np.sqrt(w), maxfev=20000)
    return VariogramModel(*popt)


def detection_limit(model: VariogramModel, area_m2: float, cell_m2: float,
                    z: float = 1.96):
    """Uncertainty of a mean difference over ``area_m2`` and its detection limit.

    White (nugget) variance averages over all pixels; correlated (sill)
    variance averages over N_eff ~ area / (pi * range^2). Returns
    (sigma_mean, lod) in metres; ``lod = z * sigma_mean``.
    """
    n_pix = max(area_m2 / cell_m2, 1.0)
    n_eff = max(area_m2 / (np.pi * model.range_ ** 2), 1.0)
    var_mean = model.nugget / n_pix + model.sill / n_eff
    sigma_mean = float(np.sqrt(var_mean))
    return sigma_mean, z * sigma_mean
