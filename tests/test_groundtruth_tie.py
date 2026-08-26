"""The tie estimator, against surfaces whose answer is known analytically.

The regression test that matters is ``test_plane_fit_is_radius_dependent_on_a_local_high``
paired with ``test_order_two_removes_the_radius_dependence``: the first reproduces the
pathology measured on real checkpoint 2210 (plane fit runs -200 -> -589 -> -1169 mm as
the radius goes 5 -> 10 -> 20 m), the second shows the order-2 surface removes it. The
first is the "without the fix" case, so the pair proves the fix bites.
"""
import numpy as np
import pytest

from lidar_diff_icp.groundtruth import checkpoints as C
from lidar_diff_icp.groundtruth import tie as T


def _cp(elev, e=0.0, n=0.0, ptype="NVA", geoid="GEOID18"):
    return C.Checkpoint(point_id="SYN", point_type=ptype, easting=e, northing=n,
                        elevation=elev, elevation_units="m", horizontal_crs="EPSG:26915",
                        vertical_datum="NAVD88", geoid_model=geoid)


def _returns(x, y, z, line=1):
    return T.GroundReturns(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float),
                           np.zeros(np.size(x)), np.full(np.size(x), line, int),
                           source="synthetic", origin="test")


def _disc(n=4000, r=30.0, seed=0):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, n)
    rr = r * np.sqrt(rng.uniform(0, 1, n))
    return rr * np.cos(t), rr * np.sin(t)


# ------------------------------------------------------------------ analytic surfaces

def test_tilted_plane_gives_the_surveyed_elevation_exactly():
    """A checkpoint ON a tilted plane: the tie is zero at every radius, any order."""
    x, y = _disc()
    z = 100.0 + 0.10 * x - 0.05 * y                 # 10% east, 5% south
    g = _returns(x, y, z)
    for order in (1, 2):
        t = T.estimate_tie(_cp(100.0), g, line=1, geoid_shift_m=0.0, surface_order=order)
        assert t.tie_mm == pytest.approx(0.0, abs=1e-6)
        for r in t.curve:
            assert r.ok and r.z_lidar_m == pytest.approx(100.0, abs=1e-6)
        assert t.radius_spread_mm == pytest.approx(0.0, abs=1e-6)
        assert t.curve[0].slope_deg == pytest.approx(np.degrees(np.arctan(np.hypot(0.10, 0.05))))


def test_a_known_offset_is_recovered():
    x, y = _disc()
    z = 100.0 + 0.10 * x - 0.05 * y - 0.123          # lidar reads 123 mm low
    t = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0)
    assert t.tie_mm == pytest.approx(123.0, abs=1e-3)


def test_the_geoid_term_enters_with_the_sign_that_lifts_gen1():
    x, y = _disc()
    z = 100.0 + 0.0 * x
    t0 = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0)
    t1 = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.067)
    assert t0.tie_mm == pytest.approx(0.0, abs=1e-6)
    assert t1.tie_mm == pytest.approx(-67.0, abs=1e-3)   # gen1 lifted 67 mm -> tie drops


def test_the_swath_shift_is_applied_before_estimating():
    x, y = _disc()
    z = 100.0 + 0.10 * x
    t = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0,
                       swath_shift_m=(0.0, 0.0, -0.020))
    assert t.tie_mm == pytest.approx(20.0, abs=1e-3)


# ------------------------------------------------------- the radius pathology, and the fix

def _local_high(curv=-0.004, r=30.0, spacing=0.25):
    """A convex dome, z = 100 + curv*(x^2+y^2) with curv < 0: a road crown / shoulder,
    the geometry ELBAEXT2_SCOPE measured at checkpoint 2210 (the surveyed height sits at
    the p95 of returns within 5 m, with ground falling away inside 10 m).

    Sampled on a REGULAR grid, not at random, so the analytic plane bias below is checked
    against the estimator and not against Monte-Carlo noise."""
    g = np.arange(-r, r + spacing, spacing)
    X, Y = np.meshgrid(g, g)
    m = (X * X + Y * Y) <= r * r
    x = X[m]; y = Y[m]
    return x, y, 100.0 + curv * (x * x + y * y)


def test_plane_fit_is_radius_dependent_on_a_local_high():
    """WITHOUT the fix (surface_order=1): the estimate walks with the radius,
    reproducing the real 2210 curve. This is the failing case the fix must remove."""
    x, y, z = _local_high()
    t = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0,
                       surface_order=1)
    by_r = {r.radius_m: r for r in t.curve}
    # Analytic: over a disc of radius R the LS plane's constant term is the MEAN of
    # z = 100 + c*r^2, and mean(r^2) = R^2/2, so it reads 100 + c*R^2/2 -- quadratic in R.
    # (The median residual adds nothing: r^2 is uniform on [0, R^2], so its median is
    # also R^2/2.) Checked where the window holds enough points to determine it.
    for R in (10.0, 15.0, 20.0, 25.0):
        assert by_r[R].z_lidar_m == pytest.approx(100.0 - 0.004 * R * R / 2.0, abs=2e-3)
    zs = [by_r[R].z_lidar_m for R in T.radius_ladder(5.0)]
    assert all(b < a for a, b in zip(zs, zs[1:]))          # sinks monotonically with R
    assert t.radius_spread_mm > 100.0            # a decimetre-scale walk over 2.5-10 m
    ok, why = t.verdict(50.0, tolerance_source="synthetic test")
    assert ok is False and "radius spread" in why


def test_order_two_removes_the_radius_dependence():
    """WITH the fix: a 2nd-order surface carries the curvature term, so the read at the
    mark is the same at every radius."""
    x, y, z = _local_high()
    t = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0,
                       surface_order=2)
    for r in t.curve:
        assert r.ok and r.z_lidar_m == pytest.approx(100.0, abs=1e-6)
    assert t.radius_spread_mm == pytest.approx(0.0, abs=1e-6)
    assert t.tie_mm == pytest.approx(0.0, abs=1e-6)
    assert t.verdict(50.0, tolerance_source="synthetic test")[0] is True


def test_the_radius_curve_is_always_returned():
    x, y, z = _local_high()
    t = T.estimate_tie(_cp(100.0), _returns(x, y, z), line=1, geoid_shift_m=0.0)
    assert [r.radius_m for r in t.curve] == list(T.radius_ladder(5.0))
    assert len(t.table_rows()) == len(t.curve)
    for c in ("R_m", "tie_mm", "fit_rms_mm", "n_lines", "scan_deg"):
        assert c in t.table_columns()


def test_radius_ladder_is_tied_to_the_grid_resolution():
    assert T.radius_ladder(5.0) == (2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0)
    assert T.radius_ladder(2.0) == (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0)


# -------------------------------------------------------------------------- refusals

def test_too_few_returns_is_reported_not_guessed():
    x = np.array([0.0, 1.0, -1.0]); y = np.array([0.0, 1.0, -1.0]); z = np.zeros(3)
    zh, info = T.ground_elevation_at(x, y, z, 0.0, 0.0, 5.0, surface_order=2)
    assert np.isnan(zh) and "6 coefficients" in info["note"]
    zh1, info1 = T.ground_elevation_at(x, y, z, 0.0, 0.0, 5.0, surface_order=1)
    assert np.isfinite(zh1) and info1["n"] == 3


def test_empty_window_yields_an_unusable_verdict():
    x, y = _disc(n=200, r=2.0)
    t = T.estimate_tie(_cp(100.0), _returns(x + 500.0, y, np.full(x.size, 100.0)),
                       line=1, geoid_shift_m=0.0)
    ok, why = t.verdict(1000.0, tolerance_source="synthetic test")
    assert ok is False and "no radius produced an estimate" in why


def test_an_unknown_geoid_cannot_reach_a_tie():
    x, y = _disc()
    with pytest.raises(C.UnknownDatumError):
        T.estimate_tie(_cp(100.0, geoid="") , _returns(x, y, np.full(x.size, 100.0)),
                       line=1, geoid_shift_m=0.0)


def test_a_vva_checkpoint_is_kept_and_flagged_not_dropped():
    x, y = _disc()
    t = T.estimate_tie(_cp(100.0, ptype="VVA"), _returns(x, y, np.full(x.size, 100.0)),
                       line=1, geoid_shift_m=0.0)
    assert t.point_type == "VVA"
    assert any("VVA" in n for n in t.notes)
    assert np.isfinite(t.tie_mm)


def test_mixed_flight_lines_are_flagged():
    x, y = _disc()
    z = np.full(x.size, 100.0)
    g = T.GroundReturns(x, y, z, np.zeros(x.size),
                        np.where(x > 0, 1, 2), source="synthetic", origin="test")
    mixed = T.estimate_tie(_cp(100.0), g, line=None, geoid_shift_m=0.0)
    assert any("mixes them" in n for n in mixed.notes)
    single = T.estimate_tie(_cp(100.0), g, line=1, geoid_shift_m=0.0)
    assert single.notes == []


def test_params_carry_their_origin():
    x, y = _disc()
    t = T.estimate_tie(_cp(100.0), _returns(x, y, np.full(x.size, 100.0)), line=1,
                       geoid_shift_m=0.0)
    names = {p.name: p for p in t.params}
    assert names["ground_quantile"].value == 0.50
    assert "corrections.json" in names["ground_quantile"].why
    assert names["surface_order"].value == 2
    assert all(p.src in ("andy", "repo", "MINE") for p in t.params)
    assert all(p.why for p in t.params)
