"""End-to-end and utility tests for the differencing pipeline.

These use small synthetic last-return clouds (two overlapping swaths, gps_time,
a known bump in the "after" epoch) so the whole workflow -- read -> per-swath
align -> tie -> correction surface -> along-track drift -> gridded low-percentile
DoD -> LoD -- runs without any downloaded data.
"""
import numpy as np
import laspy
import pytest

from lidar_diff_icp.pipeline import difference_dem, rasterize, heteroscedastic_lod


X0, Y0, W = 1000.0, 2000.0, 200.0
BOUNDS = (X0, Y0, X0 + W, Y0 + W)
BUMP_XY = (X0 + 60.0, Y0 + 140.0)          # off the dome apex


def _ground(x, y):
    # a gentle dome: radial slopes (all aspects, < ~2.5 deg) so the ground is
    # "stable" everywhere yet WELL-CONDITIONED for the Nuth & Kaeaeb tie (a
    # near-flat or single-aspect surface leaves the horizontal shift unconstrained).
    cx, cy = X0 + W / 2, Y0 + W / 2
    return 105.0 - 1.5e-4 * ((x - cx) ** 2 + (y - cy) ** 2)


def _bump(x, y):  # a flat-topped 1 m "deposit" (radius 30 m) in the after epoch.
    # Flat-topped and above dz_thresh (0.7 m) everywhere within it, so the whole
    # patch is excluded from the correction surface's stable set and preserved --
    # a clearly detectable change the pipeline must keep. (A Gaussian bump has a
    # ring where 0.15 < |dz| < 0.7 that the CS would partly absorb.)
    r = np.hypot(x - BUMP_XY[0], y - BUMP_XY[1])
    return np.where(r < 30.0, 1.0, 0.0)


def _write_laz(path, x, y, z, psid, gps):
    hdr = laspy.LasHeader(point_format=1, version="1.2")   # format 1 carries gps_time
    hdr.offsets = [x.min(), y.min(), z.min()]; hdr.scales = [0.01, 0.01, 0.01]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.return_number = np.ones(len(x), np.uint8)          # single returns (rn==nr)
    las.number_of_returns = np.ones(len(x), np.uint8)
    las.point_source_id = psid.astype(np.uint16)
    las.gps_time = gps.astype(np.float64)
    las.classification = np.zeros(len(x), np.uint8)
    las.write(str(path))


def _make_tiles(tmp_path):
    rng = np.random.default_rng(0)
    # before: two swaths overlapping in x in [X0+80, X0+120], ~3 pts/m^2 (dense
    # enough for a stable per-cell low-10% ground)
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    zb = _ground(xb, yb) + rng.normal(0, 0.03, len(xb))
    _write_laz(tmp_path / "before.laz", xb, yb, zb, ps, yb)   # gps_time ~ along-track (y)
    # after: dense reference with the bump
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz(tmp_path / "after.laz", xa, ya, za, np.ones(na), ya)
    return str(tmp_path / "before.laz"), str(tmp_path / "after.laz")


def test_difference_dem_recovers_bump_and_zero_on_stable(tmp_path):
    before, after = _make_tiles(tmp_path)
    r = difference_dem(before, after, BOUNDS, res=5.0, ground_q=0.10)
    dod = r["dod"]; res = r["res"]
    ci = int((BUMP_XY[0] - X0) / res); ri = int((BUMP_XY[1] - Y0) / res)
    # the 1 m bump is recovered (above dz_thresh, so kept as real change)
    assert dod[ri, ci] > 0.7, f"bump not recovered: {dod[ri, ci]}"
    # away from the bump the difference is ~0 and tight
    far = np.ones_like(dod, bool)
    yy, xx = np.mgrid[0:dod.shape[0], 0:dod.shape[1]]
    d2 = (xx - ci) ** 2 + (yy - ri) ** 2
    far &= d2 > (15 ** 2)
    m = far & np.isfinite(dod)
    assert abs(np.median(dod[m])) < 0.05
    assert r["stable_sigma"] < 0.08


def test_rasterize_roundtrip():
    # two cells (10 m) each with two points; median per cell
    x = np.array([1., 2., 11., 12.]); y = np.array([1., 1., 1., 1.])
    v = np.array([10., 20., 30., 40.])
    g = rasterize(x, y, v, (0., 0., 20., 10.), res=10.0, agg="median")
    assert g.shape == (1, 2)
    assert abs(g[0, 0] - 15.0) < 1e-9 and abs(g[0, 1] - 35.0) < 1e-9


def test_heteroscedastic_lod_optional():
    """If xdem is importable (its import needs PROJ_DATA unset), the model must
    recover slope-scaled noise; otherwise the function returns None (fallback)."""
    try:
        import xdem  # noqa: F401
    except Exception:
        pytest.skip("xdem not importable in this environment (PROJ)")
    rng = np.random.default_rng(0); n = 200
    yy, xx = np.mgrid[0:n, 0:n]
    slope = (xx / n * 30.0).astype(float)          # 0..30 deg across x
    curv = (yy / n * 2.0).astype(float)            # 0..2 across y (non-degenerate)
    dod = rng.normal(0, 1, (n, n)) * (0.02 + 0.006 * slope)  # sigma grows with slope
    stable = np.ones((n, n), bool)
    lod = heteroscedastic_lod(dod, slope, curv, stable)
    assert lod is not None and lod.shape == (n, n)
    # LoD must increase from shallow to steep
    assert np.nanmedian(lod[:, n - 20:]) > 1.5 * np.nanmedian(lod[:, :20])
