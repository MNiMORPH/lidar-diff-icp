"""q2(cover): the matching relation between gen2's percentile and gen1's median.

The recipe is documented and load-bearing (analysis/ridgelines/Q2_COVER_RELATION.md), and
the failure that motivated these tests is that a DEFAULT changed underneath it -- refcells
began excluding the valley -- and the published fit silently stopped reproducing. So the
constants are pinned here, and the behaviours that would go wrong quietly are covered.
"""
import numpy as np
import pytest

from lidar_diff_icp import q2cover as Q


def test_the_documented_recipe_constants_are_what_the_doc_says():
    """If any of these move, the relation moves. They are not tuning knobs."""
    assert Q.BIN_WIDTH == 0.05
    assert Q.MIN_GEN1_RETURNS == 5
    assert Q.MIN_GEN2_RETURNS == 10
    assert Q.Q1 == 0.50
    assert Q.Q2_AT_ZERO == 0.50


def test_q2_at_is_the_relation_and_stays_a_percentile():
    assert Q.q2_at(0.0, -0.1871) == pytest.approx(0.50)
    assert Q.q2_at(0.5, -0.20) == pytest.approx(0.40)
    # linear extrapolation would leave (0,1); a percentile cannot
    assert 0.0 < Q.q2_at(10.0, -0.20) < 1.0
    assert 0.0 < Q.q2_at(-10.0, -0.20) < 1.0


def test_hist_quantile_recovers_a_known_quantile():
    """A flat histogram over -1..1 m: the median must land at 0, p25 at -0.5 m."""
    nb = 100; dz = 0.02; zlo = -1.0
    C = np.cumsum(np.ones((1, nb)), 1)
    ntot = np.array([float(nb)])
    assert Q.hist_quantile(C, ntot, 0.50, zlo, dz)[0] == pytest.approx(0.0, abs=25)
    assert Q.hist_quantile(C, ntot, 0.25, zlo, dz)[0] == pytest.approx(-500.0, abs=25)


def test_ragged_quantile_matches_numpy_per_cell():
    rng = np.random.default_rng(0)
    cell = np.repeat(np.arange(4), 25)
    val = rng.normal(size=100)
    vs, off, n = Q.ragged_sorted(cell, val, 4)
    got = Q.ragged_quantile(vs, off, n, 0.5, np.arange(4))
    want = [np.quantile(val[cell == c], 0.5) for c in range(4)]
    np.testing.assert_allclose(got, want)


def _synthetic(cover, true_q2, ncell=400, nb=150, dz=0.02, zlo=-1.0, seed=0):
    """Cells whose gen2 column is uniform, with gen1's median placed at a KNOWN percentile."""
    C = np.cumsum(np.ones((ncell, nb)), 1)
    ntot = np.full(ncell, float(nb))
    g1 = Q.hist_quantile(C, ntot, true_q2, zlo, dz)
    return np.full(ncell, cover), g1, C, ntot, zlo, dz


def test_solve_recovers_a_planted_percentile():
    cov, g1, C, ntot, zlo, dz = _synthetic(0.30, 0.42)
    bins, un = Q.solve_q2_per_bin(cov, g1, C, ntot, zlo, dz)
    assert not un
    assert len(bins) == 1
    assert bins[0]["q2"] == pytest.approx(0.42, abs=0.01)


def test_an_unmatchable_bin_is_reported_not_raised():
    """gen1's median outside gen2's whole column. This HAPPENS -- elba, valley excluded,
    bin 0.70-0.75, n=15 -- and it must not crash the fit or vanish."""
    cov, g1, C, ntot, zlo, dz = _synthetic(0.30, 0.50)
    g1 = g1 - 10_000.0                      # push gen1 far below gen2's column
    bins, un = Q.solve_q2_per_bin(cov, g1, C, ntot, zlo, dz)
    assert bins == []
    assert len(un) == 1
    assert un[0]["side"] == "below"
    assert un[0]["n"] == 400


def test_fit_slope_passes_through_the_imposed_anchor():
    bins = [dict(cover=c, q2=0.50 - 0.19 * c, n=1000) for c in (0.1, 0.2, 0.3, 0.4)]
    slope, se = Q.fit_slope(bins)
    assert slope == pytest.approx(-0.19, abs=1e-9)
    assert Q.q2_at(0.0, slope) == pytest.approx(0.50)


def test_free_intercept_recovers_a_line_that_does_not_pass_through_the_anchor():
    """The check that caught the bare-ground disagreement between tiles."""
    bins = [dict(cover=c, q2=0.57 - 0.35 * c, n=1000) for c in (0.05, 0.1, 0.2, 0.3, 0.4)]
    a, b, sa, sb = Q.free_intercept(bins)
    assert a == pytest.approx(0.57, abs=1e-6)
    assert b == pytest.approx(-0.35, abs=1e-6)


def test_fit_tile_names_the_missing_producer_rather_than_failing_obscurely():
    with pytest.raises(FileNotFoundError, match="nearground_class_split"):
        Q.fit_tile("tests")
