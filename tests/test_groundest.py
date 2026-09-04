"""The four ground estimators, now that they can be called at all.

They were nested inside difference_dem, closing over its locals, so none of this could be
written before: nothing could reach them.
"""
import numpy as np
import pandas as pd
import pytest

from lidar_diff_icp import groundest

GRID = (0.0, 0.0, 5.0, 24, 20)          # X0, Y0, res, nx, ny


def _cloud(seed=0, n=60000, curv=0.0008):
    X0, Y0, res, nx, ny = GRID
    rng = np.random.default_rng(seed)
    x = rng.uniform(X0, X0 + nx * res, n)
    y = rng.uniform(Y0, Y0 + ny * res, n)
    z = 100 + 0.04 * x - 0.02 * y + curv * x * x + rng.normal(0, 0.05, n)
    return x, y, z


def test_cellstat_is_a_plain_groupby_quantile():
    X0, Y0, res, nx, ny = GRID
    x, y, z = _cloud()
    got = groundest.cellstat(x, y, z, "ground", GRID, 0.37)
    ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
    s = pd.Series(z).groupby(iy * nx + ix).quantile(0.37)
    exp = np.full(nx * ny, np.nan); exp[s.index.values] = s.values
    assert np.array_equal(got, exp.reshape(ny, nx), equal_nan=True)


def test_plane_ground_is_unbiased_under_a_tilt_where_a_quantile_is_not():
    """The reason this estimator exists: on a tilted cell a quantile of the points lands at
    a tilt-correlated spot, while the plane is read at the cell CENTRE."""
    X0, Y0, res, nx, ny = GRID
    x, y, z = _cloud(curv=0.0)                       # pure tilt, no curvature, no noise
    z = 100 + 0.04 * x - 0.02 * y
    p = groundest.plane_ground(x, y, z, GRID, 0.50)
    cx = X0 + (np.arange(nx) + 0.5) * res
    cy = Y0 + (np.arange(ny) + 0.5) * res
    truth = 100 + 0.04 * cx[None, :] - 0.02 * cy[:, None]
    assert np.nanmax(np.abs(p - truth)) < 1e-6


def test_poly2_ground_is_unbiased_under_CURVATURE_where_a_plane_is_not():
    x, y, z = _cloud()
    z = 100 + 0.0008 * x ** 2                        # curved, noiseless
    X0, Y0, res, nx, ny = GRID
    cx = X0 + (np.arange(nx) + 0.5) * res
    truth = np.repeat((100 + 0.0008 * cx ** 2)[None, :], ny, axis=0)
    p2 = groundest.poly2_ground(x, y, z, GRID, 0.50)
    pl = groundest.plane_ground(x, y, z, GRID, 0.50)
    assert np.nanmax(np.abs(p2 - truth)) < np.nanmax(np.abs(pl - truth))


def test_slope_normal_adds_the_plane_back():
    X0, Y0, res, nx, ny = GRID
    x, y, z = _cloud()
    plane = (np.full(nx * ny, 100.0), np.full(nx * ny, 0.04), np.full(nx * ny, -0.02))
    g = groundest.slope_normal_ground(x, y, z, GRID, 0.50, plane)
    assert np.isfinite(g).all() and np.nanmedian(g) > 90


def test_the_dispatcher_refuses_slope_normal_without_a_plane():
    x, y, z = _cloud()
    with pytest.raises(ValueError, match="needs plane"):
        groundest.estimate_ground(x, y, z, GRID, 0.50, "slope_normal")
