"""An unmeasurable variance term is not a zero one.

SE^2 = SE08^2 + SE21^2. cell_plane_roughness returns NaN where a cell holds fewer than
min_n points, so a missing term marks the THINNEST cells -- exactly the ones that deserve
the widest detection limit. Substituting zero there gave them the most optimistic error
available.
"""
import numpy as np

from lidar_diff_icp.pipeline import cell_plane_roughness


def _cell(n, spread=0.0, seed=0):
    """n points inside one 5 m cell, with a given vertical spread."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.5, 4.5, n)
    y = rng.uniform(0.5, 4.5, n)
    z = rng.normal(0.0, spread, n) if spread else np.zeros(n)
    return x, y, z


def test_roughness_is_nan_below_min_n_and_finite_above():
    """The condition that creates the missing term in the first place."""
    for n, want_finite in ((3, False), (5, False), (6, True), (30, True)):
        x, y, z = _cell(n, spread=0.05, seed=n)
        r = cell_plane_roughness(x, y, z, 0.0, 0.0, 5.0, 1, 1)
        assert np.isfinite(r[0, 0]) == want_finite, f"n={n}"


def test_an_unmeasurable_epoch_yields_no_standard_error():
    """The propagation itself: NaN in either term must survive to the total."""
    r08 = np.array([[0.05, 0.05, np.nan, np.nan]])
    r21 = np.array([[0.03, np.nan, 0.03, np.nan]])
    n08 = np.array([[10.0, 10.0, 2.0, 2.0]])
    n21 = np.array([[10.0, 3.0, 10.0, 3.0]])

    stderr = np.sqrt(r08**2 / np.maximum(n08, 1.0) + r21**2 / np.maximum(n21, 1.0))

    assert np.isfinite(stderr[0, 0]), "both measured -> a standard error"
    assert not np.isfinite(stderr[0, 1]), "gen2 unmeasurable -> no total"
    assert not np.isfinite(stderr[0, 2]), "gen1 unmeasurable -> no total"
    assert not np.isfinite(stderr[0, 3]), "neither measured -> no total"


def test_zero_filling_a_missing_term_understates_the_error():
    """Why it matters, in the direction it matters."""
    r08, n08 = 0.05, 9.0
    r21_unknown, n21 = np.nan, 4.0

    honest = np.sqrt(r08**2 / n08 + r21_unknown**2 / n21)
    zero_filled = np.sqrt(r08**2 / n08 + np.nan_to_num(r21_unknown)**2 / n21)

    assert not np.isfinite(honest)
    assert np.isfinite(zero_filled)
    # and the zero-filled value is exactly the one-epoch error -- the smallest it could be
    assert zero_filled == np.sqrt(r08**2 / n08)
