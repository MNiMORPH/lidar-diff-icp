"""The cross-epoch datum: lateral shift, then geoid.

This is the only stage that reads across the epochs. It was inline in difference_dem, so
its two load-bearing properties were untestable: that the geoid term is COMPUTED from the
PROJ grids rather than carried as a constant, and that an unsupported `tie` refuses instead
of silently doing something else. The removed reference_plane and parabola ties are exactly
the sort of thing that must not come back by accident.
"""
import numpy as np
import pytest

from lidar_diff_icp import pipeline, groundest


# a corner of the real elba tile, so references.geoid_difference reads real geoid grids
BOUNDS = (577492.8, 4882737.6, 578492.8, 4883237.6)
RES = 5.0
NX = int((BOUNDS[2] - BOUNDS[0]) / RES)
NY = int((BOUNDS[3] - BOUNDS[1]) / RES)
GRID = (BOUNDS[0], BOUNDS[1], RES, NX, NY)


def _scene(n=120_000, seed=11, gen2_above_m=0.030):
    """A gen1 cloud and a gen2 reference grid sitting a known distance above it."""
    X0, Y0 = BOUNDS[0], BOUNDS[1]
    rng = np.random.default_rng(seed)
    x = rng.uniform(BOUNDS[0], BOUNDS[2], n)
    y = rng.uniform(BOUNDS[1], BOUNDS[3], n)
    z = (300.0 + 6.0 * np.sin((x - X0) / 70.0) + 6.0 * np.cos((y - Y0) / 60.0)
         + 0.01 * (x - X0) + rng.normal(0.0, 0.02, n))
    ground_of = lambda a, b, c: groundest.cellstat(a, b, c, "ground", GRID, 0.50)
    Zref = ground_of(x, y, z) + gen2_above_m
    return dict(x=x, y=y, z=z, ground=np.ones(n, bool), Zref=Zref, ground_of=ground_of,
                gps_time=np.linspace(0.0, 900.0, n),
                source_id=(((x - X0) // 250).astype(np.int32) + 10),
                stable=np.isfinite(Zref), floodplain=np.zeros((NY, NX), bool))


def _run(s, **kw):
    return pipeline.apply_datum(s["x"].copy(), s["y"].copy(), s["z"].copy(), s["ground"],
                                s["Zref"], s["ground_of"], GRID, BOUNDS,
                                gps_time=s["gps_time"], source_id=s["source_id"],
                                stable=s["stable"], floodplain=s["floodplain"],
                                verbose=False, **kw)


def test_an_unsupported_tie_refuses(): 
    """reference_plane and the parabola were REMOVED, not deprecated. A tie that is not
    the geoid difference must raise, not fall through to it -- a silent fallback would
    reintroduce a fitted datum under a name that says it is geodetic."""
    s = _scene(n=2000)
    for bad in ("reference_plane", "parabola", "median", ""):
        with pytest.raises(ValueError, match="not supported"):
            _run(s, tie=bad)


def test_the_geoid_term_is_computed_not_assumed():
    """No hard-coded constant: with geoid_datum=None the term comes from the PROJ grids.
    At elba that is the GEOID03 -> GEOID18 difference, tens of mm, and it must be non-zero
    and finite or the whole cross-epoch level is wrong."""
    from lidar_diff_icp import references
    t = _run(_scene(n=20_000))["tie_info"]
    assert t["method"] == "geoid_difference"
    # compared against the source itself rather than a window I would have had to invent
    gc, gb, gcc = references.geoid_difference(BOUNDS, 26915)
    assert (t["const_m"], t["tilt_b_m_per_km"], t["tilt_c_m_per_km"]) == (gc, gb, gcc)
    assert np.isfinite([gc, gb, gcc]).all() and gc != 0.0


def test_a_stated_geoid_datum_is_applied_exactly():
    """When the caller states the datum, that is what is applied -- no re-derivation."""
    s = _scene(n=20_000)
    z0 = s["z"].copy()
    r = _run(s, geoid_datum=(0.0674, 0.0, 0.0))
    assert r["tie_info"]["const_m"] == 0.0674
    assert np.allclose(r["z"] - z0, 0.0674, atol=1e-12)      # flat: no tilt term


def test_the_geoid_tilt_is_referenced_to_the_bounds_centroid():
    """The tilt is m/km about the centre of the tile, so a point AT the centroid gets the
    constant alone. Referencing it to a corner instead would put a whole-tile bias in."""
    s = _scene(n=20_000)
    z0 = s["z"].copy()
    r = _run(s, geoid_datum=(0.05, 1.0, -2.0))               # 1 and -2 m per km
    cx, cy = r["tie_info"]["centroid"]
    assert cx == pytest.approx(0.5 * (BOUNDS[0] + BOUNDS[2]))
    assert cy == pytest.approx(0.5 * (BOUNDS[1] + BOUNDS[3]))
    dz = r["z"] - z0
    expect = 0.05 + 1.0 * (r["x"] - cx) / 1000.0 + (-2.0) * (r["y"] - cy) / 1000.0
    assert np.allclose(dz, expect, atol=1e-12)


def test_it_mutates_in_place_and_returns_the_same_arrays():
    """Documented behaviour, pinned so it cannot change silently: at statewide scale these
    are tens of millions of points and three copies is real memory. A caller that needs the
    registered cloud afterwards must copy BEFORE calling."""
    s = _scene(n=5000)
    x, y, z = s["x"], s["y"], s["z"]
    r = pipeline.apply_datum(x, y, z, s["ground"], s["Zref"], s["ground_of"], GRID, BOUNDS,
                             geoid_datum=(0.05, 0.0, 0.0), gps_time=s["gps_time"],
                             source_id=s["source_id"], stable=s["stable"],
                             floodplain=s["floodplain"], verbose=False)
    assert r["z"] is z and r["x"] is x and r["y"] is y


def test_drift_curves_are_empty_unless_the_drift_is_fitted():
    """An empty dict and a fitted-but-zero curve are different states, and corrections.json
    records which."""
    s = _scene(n=40_000)
    assert _run(s)["drift_curves"] == {}
    assert _run(s, along_track_drift=True)["drift_curves"] != {}
