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
    # correction_surface is PINNED OFF here, though it is the default since 2026-09-11.
    # This test is about the BORESIGHT mechanism, and a scanner roll is a spatially
    # varying field, so a DeLong IDW absorbs most of it before boresight sees it: with
    # the surface on, stable_sigma goes 0.0098 -> 0.0087 (ratio 0.88) and the 0.7 ratio
    # below stops biting, not because the correction stopped working but because there
    # was little left for it to remove. Isolating the mechanism means pinning the other
    # one off. The absorption itself is asserted separately below, so it stays visible.
    kw = dict(res=5.0, ground_q=0.10, ground="low_q", ground_source="last_return",
              after_ground="last_return", geoid_datum=(0.0, 0.0, 0.0),
              valley_top_m=-1e9, correction_surface=False, along_track_drift=True)
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

    # The DeLong correction surface ABSORBS this injected roll -- recorded, not hidden.
    # It is why the ratio assertion above needs the surface pinned off, and it is a real
    # cost of the 2026-09-11 default: with the surface on, a scanner roll is removed from
    # the product but is no longer VISIBLE as a roll, so it cannot be diagnosed or
    # physically modelled. Bites if someone turns the surface on inside this test.
    cs_kw = dict(kw); cs_kw["correction_surface"] = True; cs_kw["along_track_drift"] = False
    r_cs = difference_dem(before, after, BOUNDS, correct_boresight=False, **cs_kw)
    assert r_cs["stable_sigma"] < 0.5 * r_off["stable_sigma"], (
        f"the correction surface should absorb most of the injected roll: "
        f"{r_off['stable_sigma']:.4f} -> {r_cs['stable_sigma']:.4f}")
    # ... and it must NOT eat the real 1 m bump; the |dz| > 0.7 m mask excludes the bump
    # from the IDW's stable sources. Measured 0.994 of 1.0 m, BETTER than the drift's
    # 0.978, which is what refuted the concern that an IDW would absorb real change.
    assert r_cs["dod"][ri, ci] > 0.9, (
        f"the correction surface ate real localized change: {r_cs['dod'][ri, ci]:.3f} of 1.0 m")


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


def test_the_correction_surface_absorbs_a_wrong_geoid_and_hides_it(tmp_path):
    """A wrong gen1 geoid is INVISIBLE in the product once correction_surface is on.

    Andy, 2026-09-11: "DeLong's purely relative would become absolute by necessity because
    gen2 is considered absolute. So it's the same thing." That is right, and this is the
    mechanism -- and its cost. The surface is fit on the stable residual, so ANY smooth
    offset between the epochs is interpolated away: gen1 is pulled onto gen2's frame
    whatever datum was applied first.

    Measured here with the +54.87 mm battlecreek error (GEOID03 applied to a GEOID09
    survey). Surface OFF: the DoD moves by -54.870 mm, mean and max -- that is the bug, and
    it was found only because nothing absorbed it. Surface ON: 0.000 mm, mean and max.

    So with the surface on, the geoid term is no longer load-bearing for the product's
    NUMBERS. It stays as a physically meaningful, transferable decomposition -- but it can
    no longer be validated from the product, because the product is now insensitive to it.
    Same pattern as the scanner roll the surface absorbs: better product, blind diagnosis.
    That is the argument for keeping the geoid and the control marks as INDEPENDENT checks
    OUTSIDE the pipeline, which is exactly where they are being moved (task #68).

    Bites if someone makes the surface respect the datum, or drops the geoid term believing
    it was never needed -- it IS needed whenever the surface is off.
    """
    rng = np.random.default_rng(0)
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    zb = _ground(xb, yb) + rng.normal(0, 0.02, len(xb))
    _write_laz14(tmp_path / "before.laz", xb, yb, zb, ps, yb, np.zeros(len(xb)))
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz14(tmp_path / "after.laz", xa, ya, za, np.ones(na), ya, np.zeros(na))
    before = str(tmp_path / "before.laz"); after = str(tmp_path / "after.laz")

    kw = dict(res=5.0, ground_q=0.10, ground="low_q", ground_source="last_return",
              after_ground="last_return", valley_top_m=-1e9, along_track_drift=False)
    GEOID_ERR = 0.05487                       # the measured battlecreek error, in metres

    def dod(cs, g):
        # apply_geoid FORCED ON in both arms. Since 2026-09-11 the delong route defaults it
        # off, and this test is about whether the surface ABSORBS a geoid that WAS applied
        # -- not about whether the route applies one. Without the force, both arms would
        # skip the geoid and the test would pass by measuring nothing.
        return difference_dem(before, after, BOUNDS, correction_surface=cs,
                              apply_geoid=True, geoid_datum=(g, 0.0, 0.0), **kw)["dod"]

    off = dod(False, GEOID_ERR) - dod(False, 0.0)
    on = dod(True, GEOID_ERR) - dod(True, 0.0)
    m_off, m_on = np.isfinite(off), np.isfinite(on)

    # OFF: the error passes straight through, to the millimetre
    assert abs(1000 * np.nanmean(off[m_off]) + 54.87) < 0.5, (
        f"without the surface a geoid error must pass through: "
        f"{1000 * np.nanmean(off[m_off]):+.3f} mm")
    # ON: absorbed entirely -- this is what makes the product insensitive, and blind
    assert 1000 * np.nanmax(np.abs(on[m_on])) < 1.0, (
        f"the correction surface should absorb a wrong geoid completely: "
        f"max {1000 * np.nanmax(np.abs(on[m_on])):.3f} mm survived")


def test_both_routes_run_and_stay_distinct(tmp_path):
    """The fork is real: `independent` and `delong` both run, and differ.

    Andy, 2026-09-11: "Ensure that we keep our old pipeline intact so we can run it.
    And: set up the possibility of working through the DeLong method."

    A route is only a coherent set of defaults; every switch still works alone, so any
    combination stays reachable. This pins three things:

      * `independent` reproduces the pre-2026-09-11 pipeline -- drift on, geoid applied,
        no correction surface -- and still produces a DoD.
      * `delong` runs with no gen1_geoid at all. Under `independent` that same call is a
        REFUSAL, because a defaulted geoid silently added +54.87 mm at Battle Creek.
      * the two give DIFFERENT answers. If they ever agree to the millimetre, the route
        switch has stopped being wired to anything.

    Bites if a route's defaults are collapsed into one another, if the geoid guard is
    dropped from the independent route, or if `route=` stops reaching apply_datum.
    """
    rng = np.random.default_rng(3)
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    zb = _ground(xb, yb) + rng.normal(0, 0.02, len(xb))
    _write_laz14(tmp_path / "before.laz", xb, yb, zb, ps, yb, np.zeros(len(xb)))
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz14(tmp_path / "after.laz", xa, ya, za, np.ones(na), ya, np.zeros(na))
    before = str(tmp_path / "before.laz"); after = str(tmp_path / "after.laz")
    kw = dict(res=5.0, ground_q=0.10, ground="low_q", ground_source="last_return",
              after_ground="last_return", valley_top_m=-1e9)
    ci = int((BUMP_XY[0] - X0) / 5.0); ri = int((BUMP_XY[1] - Y0) / 5.0)

    # the old pipeline, unchanged, still runs -- it just has to be asked for
    r_ind = difference_dem(before, after, BOUNDS, route="independent",
                           geoid_datum=(0.0, 0.0, 0.0), **kw)
    assert r_ind["corrections"]["route"] == "independent"
    assert r_ind["corrections"]["along_track_drift"] is True
    assert r_ind["corrections"]["correction_surface"] is False
    assert r_ind["corrections"]["geoid_applied"] is True
    assert r_ind["dod"][ri, ci] > 0.7

    # the DeLong route needs no geoid at all
    r_del = difference_dem(before, after, BOUNDS, route="delong", **kw)
    assert r_del["corrections"]["route"] == "delong"
    assert r_del["corrections"]["correction_surface"] is True
    assert r_del["corrections"]["along_track_drift"] is False
    assert r_del["corrections"]["geoid_applied"] is False
    assert r_del["corrections"]["cross_epoch_datum"]["const_m"] is None
    assert r_del["dod"][ri, ci] > 0.7

    # the geoid GUARD survives on the route that uses it
    with pytest.raises(ValueError, match="gen1_geoid is required"):
        difference_dem(before, after, BOUNDS, route="independent", **kw)

    # and the two routes are actually different products
    d = r_del["dod"] - r_ind["dod"]
    m = np.isfinite(d)
    assert np.nanmax(np.abs(d[m])) > 0.001, (
        "the two routes gave the same DoD to the millimetre -- the route switch is "
        "no longer wired to anything")

    # an explicit switch still overrides its route, so every combination stays reachable
    r_mix = difference_dem(before, after, BOUNDS, route="delong", along_track_drift=True,
                           **kw)
    assert r_mix["corrections"]["along_track_drift"] is True
    assert r_mix["corrections"]["correction_surface"] is True


def test_the_chain_refuses_orders_the_measurements_rule_out():
    """`run_chain` enforces the two orderings that are method, not taste.

    Both are refusals rather than warnings because the wrong order does not raise on its
    own -- it produces a quietly WORSE product, which is the hardest kind of error to
    notice later.
    """
    from lidar_diff_icp import chain as C
    ctx = C.Ctx(x=None, y=None, z=None, ground=None, Zref=None, ground_of=None,
                grid=None, bounds=None)

    # drift before the surface converts a spatial field into per-line offsets, and the
    # surface cannot undo them: 19.91 mm mean per-line step that way vs 5.49 the other
    with pytest.raises(ValueError, match="BEFORE CorrectionSurface"):
        C.run_chain([C.LateralShift(), C.AlongTrackDrift(), C.CorrectionSurface()], ctx)

    # z before x,y lets terrain slope leak into the elevation difference
    with pytest.raises(ValueError, match="not 'lateral_shift'"):
        C.run_chain([C.CorrectionSurface(), C.LateralShift()], ctx)

    # the two shipped orders are accepted (they fail later, on the None arrays, not on order)
    for good in ([C.LateralShift(), C.CorrectionSurface()],
                 [C.LateralShift(), C.GeoidConversion(gen1_geoid="x"), C.AlongTrackDrift()]):
        with pytest.raises(Exception) as e:
            C.run_chain(good, ctx)
        assert "BEFORE CorrectionSurface" not in str(e.value)
        assert "not 'lateral_shift'" not in str(e.value)


def test_the_differenced_gen2_surface_is_shipped_and_reconstructs_the_dod(tmp_path):
    """`z_after_differenced` IS the surface the DoD was taken against; `z_after` is not.

    `z_after` is the q = 0.50 grid. It is the hillshade backdrop AND the reference surface
    groundq.reference_surface hangs the near-ground columns off, so it must stay the median
    grid -- the cover correction measures percentiles against it. That means it is NOT
    what gets differenced once any gen2-side step runs, and shipping it alone let a
    corrected product carry a gen2 raster that was not what it differenced.

    Caught 2026-09-11 by comparing two cover runs' z_after, finding them byte-identical,
    and nearly concluding the correction had not applied -- when the DoD had in fact moved
    on 696,131 cells.

    Bites if z_after_differenced stops being returned, or is quietly aliased to z_after.
    """
    rng = np.random.default_rng(11)
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    zb = _ground(xb, yb) + rng.normal(0, 0.02, len(xb))
    _write_laz14(tmp_path / "before.laz", xb, yb, zb, ps, yb, np.zeros(len(xb)))
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz14(tmp_path / "after.laz", xa, ya, za, np.ones(na), ya, np.zeros(na))

    r = difference_dem(str(tmp_path / "before.laz"), str(tmp_path / "after.laz"), BOUNDS,
                       route="delong", res=5.0, ground_q=0.10, ground="low_q",
                       ground_source="last_return", after_ground="last_return",
                       valley_top_m=-1e9)
    assert "z_after_differenced" in r, "the differenced gen2 surface must be shipped"
    zd = r["z_after_differenced"]
    assert zd.shape == r["dod"].shape

    # With no gen2 step the two must coincide numerically, so a reader of an uncorrected
    # product is not misled by there being two names. (They may even be the SAME object:
    # with ground="low_q" the pipeline does Zref = Z21 outright. Identity is not the
    # property under test -- agreement is.)
    d = zd - r["z_after"]
    m = np.isfinite(d)
    assert not m.any() or np.nanmax(np.abs(d[m])) < 1e-9, (
        f"with no gen2 correction the two grids must coincide; "
        f"max|diff| {np.nanmax(np.abs(d[m])):.3e} m")

    # And a gen2 step must REBIND Zref, which is what makes the two diverge and is the
    # whole reason the second grid has to be shipped. Tested on the chain directly rather
    # than through a pipeline run, which would need a fitted curve and a canopy raster.
    from lidar_diff_icp import chain as C

    class _Lower(C.Gen2Correction):
        name = "test_lower"

        def apply(self, ctx):
            ctx.Zref = ctx.Zref - 0.05          # rebinds, as CoverPercentile does

    g2 = C.Gen2Ctx(Zref=r["z_after"].copy(), Z21=r["z_after"], after_laz="", grid=None)
    before = g2.Zref
    C.run_gen2_chain([_Lower()], g2)
    assert g2.Zref is not before, "a gen2 step must rebind Zref, not mutate it in place"
    dd = g2.Zref - r["z_after"]
    mm = np.isfinite(dd)
    assert abs(np.nanmedian(dd[mm]) + 0.05) < 1e-9
    assert g2.record["gen2_chain"] == ["test_lower"]
