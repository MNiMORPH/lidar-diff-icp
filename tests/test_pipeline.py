"""End-to-end and utility tests for the differencing pipeline.

These use small synthetic last-return clouds (two overlapping swaths, gps_time,
a known bump in the "after" epoch) so the whole workflow -- read -> per-swath
align -> tie -> correction surface -> along-track drift -> gridded low-percentile
DoD -> LoD -- runs without any downloaded data.
"""
import numpy as np
import pandas as pd
import laspy
import pytest

from lidar_diff_icp.pipeline import (difference_dem, rasterize, heteroscedastic_lod,
                                     read_last_return, _stream_ground)
from lidar_diff_icp.ground import classify_ground_csf, find_pdal


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


def _write_laz14(path, x, y, z, psid, gps, scan_deg):
    """LAS 1.4 / point format 6 writer that carries a per-point scan angle (0.006 deg units,
    what the pipeline reads for the boresight term)."""
    hdr = laspy.LasHeader(point_format=6, version="1.4")
    hdr.offsets = [x.min(), y.min(), z.min()]; hdr.scales = [0.01, 0.01, 0.01]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.return_number = np.ones(len(x), np.uint8); las.number_of_returns = np.ones(len(x), np.uint8)
    las.point_source_id = psid.astype(np.uint16); las.gps_time = gps.astype(np.float64)
    las.classification = np.zeros(len(x), np.uint8)
    las.scan_angle = np.round(scan_deg / 0.006).astype(np.int16)
    las.write(str(path))


def test_boresight_correction_recovers_injected_roll(tmp_path):
    """End-to-end: inject a known scanner roll into two overlapping gen1 swaths, then check
    difference_dem(correct_boresight=True) recovers it, records None when off, and keeps the
    real bump. Bites: without the wiring the recorded roll is None and the roll tilt survives."""
    rng = np.random.default_rng(1)
    ROLL = 3.0                                             # mm/deg injected
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    sc = np.concatenate([(x1 - (X0 + 60)) * 0.3, (x2 - (X0 + 140)) * 0.3])   # opposite, x-varying
    zb = _ground(xb, yb) + ROLL * sc / 1000.0 + rng.normal(0, 0.02, len(xb))  # inject roll tilt
    _write_laz14(tmp_path / "before.laz", xb, yb, zb, ps, yb, sc)
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz14(tmp_path / "after.laz", xa, ya, za, np.ones(na), ya, np.zeros(na))
    before = str(tmp_path / "before.laz"); after = str(tmp_path / "after.laz")
    # A synthetic flat-ish tile has no valley; state the cut rather than let anything
    # compute one. The caller ALWAYS says which (Andy, 2026-09-04).
    kw = dict(res=5.0, ground_q=0.10, ground="low_q", ground_source="last_return",
              after_ground="last_return", geoid_datum=(0.0, 0.0, 0.0),
              valley_top_m=-1e9)
    r_off = difference_dem(before, after, BOUNDS, correct_boresight=False, **kw)
    r_on = difference_dem(before, after, BOUNDS, correct_boresight=True, **kw)
    assert r_off["corrections"]["boresight_roll_mm_per_deg"] is None
    b = r_on["corrections"]["boresight_roll_mm_per_deg"]
    assert b is not None and 2.0 < b < 4.0, f"injected 3.0 mm/deg, recovered {b}"
    # applying it flattens the roll-induced cross-swath disagreement (bites on the apply step)
    assert r_on["stable_sigma"] < 0.7 * r_off["stable_sigma"], \
        f"correction did not flatten the roll: {r_off['stable_sigma']:.4f} -> {r_on['stable_sigma']:.4f}"
    ci = int((BUMP_XY[0] - X0) / 5.0); ri = int((BUMP_XY[1] - Y0) / 5.0)
    assert r_on["dod"][ri, ci] > 0.7, "boresight correction ate the real bump"


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
    # test the deterministic core without the PDAL/CSF dependency (the synthetic
    # clouds are last-return with no ASPRS classification, so opt into the
    # last-return heuristic for both epochs rather than 3DEP's class 2)
    # synthetic clouds sit at fake coordinates with no geoid-grid coverage, so pass an
    # explicit zero geoid datum rather than let the datum step compute (nan) from PROJ.
    r = difference_dem(before, after, BOUNDS, res=5.0, ground_q=0.10,
                       ground="low_q", ground_source="last_return",
                       after_ground="last_return", geoid_datum=(0.0, 0.0, 0.0),
                       valley_top_m=-1e9)
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


def test_read_last_return_keeps_singles(tmp_path):
    """Bare earth = last return (rn == nr) INCLUDING single returns. Dropping
    singles (rn==nr & nr>1) empties flat open ground -- the bug this guards."""
    rn = np.array([1, 2, 1], np.uint8)         # single, last-of-2, first-of-2
    nr = np.array([1, 2, 2], np.uint8)
    x = np.array([1., 2., 3.]); y = np.array([1., 1., 1.]); z = np.array([10., 11., 12.])
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [0, 0, 0]; hdr.scales = [.01, .01, .01]
    las = laspy.LasData(hdr); las.x, las.y, las.z = x, y, z
    las.return_number = rn; las.number_of_returns = nr
    las.point_source_id = np.ones(3, np.uint16); las.gps_time = np.zeros(3)
    p = tmp_path / "multi.laz"; las.write(str(p))
    r = read_last_return(p)
    assert len(r["z"]) == 2                     # single + last-of-2 kept
    assert 10.0 in r["z"] and 11.0 in r["z"] and 12.0 not in r["z"]


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


def test_stream_ground_matches_exact(tmp_path):
    """The streaming (O(cells) RAM) low-percentile ground must match the exact
    per-cell groupby.quantile to ~cm on well-sampled cells -- it never holds the
    whole cloud, so it enables statewide runs. Sparse cells are excluded (cnt>50)
    since the histogram cannot reproduce the exact's linear interpolation there."""
    rng = np.random.default_rng(0)
    X0, Y0, res, nx, ny = 0.0, 0.0, 5.0, 10, 10
    bounds = (X0, Y0, X0 + nx * res, Y0 + ny * res)
    n = 200 * nx * ny
    x = rng.uniform(X0, X0 + nx * res, n); y = rng.uniform(Y0, Y0 + ny * res, n)
    z = 100.0 + 0.05 * x - 0.03 * y + rng.exponential(0.2, n)      # ground + one-sided veg
    z[rng.integers(0, n, 20)] -= 15.0                              # low blunders (must not corrupt it)
    _write_laz(tmp_path / "c.laz", x, y, z, np.ones(n), np.zeros(n))
    g, spread, cnt = _stream_ground(str(tmp_path / "c.laz"), bounds, res, nx, ny, 0.10,
                                    after_ground="last_return")   # synthetic cloud has no class 2
    ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
    ex = pd.Series(z).groupby(iy * nx + ix).quantile(0.10)
    Ge = np.full(nx * ny, np.nan); Ge[ex.index.values] = ex.values; Ge = Ge.reshape(ny, nx)
    m = np.isfinite(g) & np.isfinite(Ge) & (cnt > 50)
    assert m.sum() > 50
    assert np.median(np.abs(g[m] - Ge[m])) < 0.02                 # cm agreement, blunder-robust


def test_classify_ground_csf_optional(tmp_path):
    """CSF ground classification via PDAL: removes a high (building/canopy) cluster
    and keeps ground, preserving point attributes. Skipped if PDAL isn't installed."""
    try:
        find_pdal()
    except Exception:
        pytest.skip("PDAL (filters.csf) not available")
    rng = np.random.default_rng(0); n = 40000
    x = rng.uniform(0, 100, n); y = rng.uniform(0, 100, n)
    z = 100.0 + 0.02 * x + rng.normal(0, 0.03, n)                 # gentle ground
    hi = rng.integers(0, n, 800); z[hi] += 8.0                    # a cluster to remove
    _write_laz(tmp_path / "c.laz", x, y, z, np.ones(n), np.zeros(n))
    out = classify_ground_csf(str(tmp_path / "c.laz"), resolution=2.0, iterations=100)
    g = laspy.read(out)
    assert 0 < len(g.x) < n                                       # filtered, not pass-through
    assert float(np.max(g.z)) < 105.0                             # the +8 m cluster is gone
    assert "gps_time" in g.point_format.dimension_names           # attributes preserved
