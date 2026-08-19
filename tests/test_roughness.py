"""Detrended within-cell roughness covariate.

The property that makes it a valid error covariate (and raw spread not one): it
recovers the true within-cell noise INDEPENDENT of slope, because the per-cell
plane fit removes relief. Raw IQR spread instead inflates with slope. Also pins
that the streaming accumulator matches the in-memory computation.
"""
import numpy as np
import pandas as pd
import laspy
import pytest

from lidar_diff_icp.pipeline import cell_plane_roughness, _stream_roughness


X0, Y0, res, nx, ny = 0.0, 0.0, 5.0, 6, 6
NOISE = 0.04


def _tilted(rng, slope, n=4000):
    x = rng.uniform(0, nx * res, n); y = rng.uniform(0, ny * res, n)
    z = 100.0 + slope * x - 0.2 * slope * y + rng.normal(0, NOISE, n)
    return x, y, z


@pytest.mark.parametrize("slope", [0.0, 0.1, 0.3, 0.6])   # up to ~31 deg
def test_roughness_is_slope_independent(slope):
    rng = np.random.default_rng(0)
    x, y, z = _tilted(rng, slope)
    rgh = cell_plane_roughness(x, y, z, X0, Y0, res, nx, ny)
    # detrended roughness ~ the true noise at every slope (within sampling)
    assert abs(np.nanmedian(rgh) - NOISE) < 0.01, np.nanmedian(rgh)
    # raw IQR spread inflates hard with slope -- the reason we must detrend
    ix = (x / res).astype(int); iy = (y / res).astype(int)
    gb = pd.Series(z).groupby(iy * nx + ix)
    raw = float((1.4826 * (gb.quantile(.75) - gb.quantile(.25)) / 1.349).median())
    if slope >= 0.3:
        assert raw > 5 * NOISE                        # >0.2 m from tilt alone
        assert np.nanmedian(rgh) < 0.5 * raw          # roughness immune to it


def test_roughness_rises_with_real_noise():
    rng = np.random.default_rng(1)
    lo = cell_plane_roughness(*_tilted_noise(rng, 0.02), X0, Y0, res, nx, ny)
    hi = cell_plane_roughness(*_tilted_noise(rng, 0.20), X0, Y0, res, nx, ny)
    assert np.nanmedian(hi) > 3 * np.nanmedian(lo)     # monotone in true scatter


def _tilted_noise(rng, sigma, n=4000):
    x = rng.uniform(0, nx * res, n); y = rng.uniform(0, ny * res, n)
    z = 100.0 + 0.3 * x + rng.normal(0, sigma, n)      # sloped, so detrending matters
    return x, y, z


def test_stream_roughness_matches_in_memory(tmp_path):
    rng = np.random.default_rng(2)
    x, y, z = _tilted(rng, 0.25, n=8000)
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [x.min(), y.min(), z.min()]; hdr.scales = [0.001, 0.001, 0.001]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.classification = np.full(len(x), 2, np.uint8)      # all ground (class 2)
    las.return_number = np.ones(len(x), np.uint8)
    las.number_of_returns = np.ones(len(x), np.uint8)
    p = tmp_path / "g.laz"; las.write(str(p))

    bounds = (X0, Y0, X0 + nx * res, Y0 + ny * res)
    # read back the quantized coords so the two paths see identical points
    r = laspy.read(str(p)); xr = np.asarray(r.x); yr = np.asarray(r.y); zr = np.asarray(r.z)
    mem = cell_plane_roughness(xr, yr, zr, X0, Y0, res, nx, ny)
    st = _stream_roughness(str(p), bounds, res, nx, ny, after_ground="class2")
    both = np.isfinite(mem) & np.isfinite(st)
    assert both.sum() > 20
    assert np.allclose(mem[both], st[both], atol=1e-6)
