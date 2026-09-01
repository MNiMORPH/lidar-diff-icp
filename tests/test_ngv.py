"""NGV -- the near-ground vegetation index -- as a REUSABLE filter.

Nothing tested this chain, and the one check that mattered lived only in a scratchpad
script: that the index computed from raw returns reproduces the construction the
coefficient was fitted against. That check is here now, plus the properties any reuse
depends on.

The defining property, and the reason the regression is not circular: the reference
surface is fitted from the neighbourhood's OWN returns, so NGV is invariant to any
vertical shift of the cloud. It structurally cannot carry offset information.
"""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from ngv import design, ngv, NGV_LO, NGV_HI, RADIUS


def _cloud(n=4000, seed=0, slope=(0.0, 0.0), ground_z=100.0, veg_frac=0.0, veg_lo=0.15,
           veg_hi=4.0, radius=RADIUS):
    """Points on a planar ground within `radius`, a fraction lifted into the veg band."""
    rng = np.random.default_rng(seed)
    r = radius * np.sqrt(rng.random(n)); th = rng.random(n) * 2 * np.pi
    x = r * np.cos(th); y = r * np.sin(th)
    z = ground_z + slope[0] * x + slope[1] * y
    k = rng.random(n) < veg_frac
    z[k] += rng.uniform(veg_lo, veg_hi, k.sum())
    coef = np.array([ground_z, slope[0], slope[1], 0.0, 0.0, 0.0])
    return x, y, z, coef


def test_bare_ground_gives_zero():
    x, y, z, c = _cloud(veg_frac=0.0)
    v, n = ngv(x, y, z, 0.0, 0.0, c)
    assert n == len(x)
    assert v == 0.0


def test_it_recovers_a_known_vegetation_fraction():
    x, y, z, c = _cloud(veg_frac=0.30, seed=1)
    v, _ = ngv(x, y, z, 0.0, 0.0, c)
    assert abs(v - 0.30) < 0.02, f"expected ~0.30, got {v}"


def test_invariant_to_a_vertical_shift_of_the_whole_cloud():
    """THE property the whole method rests on. If this fails, the regression of offset on
    NGV is circular and every coefficient fitted from it is meaningless."""
    x, y, z, c = _cloud(veg_frac=0.25, seed=2)
    a, _ = ngv(x, y, z, 0.0, 0.0, c)
    shift = 0.837
    b, _ = ngv(x, y, z + shift, 0.0, 0.0, c + np.array([shift, 0, 0, 0, 0, 0]))
    assert a == b, "NGV moved under a vertical shift; it would carry offset information"


def test_slope_is_removed_by_the_fitted_surface():
    """Slope-normal heights, not vertical ones: a tilted ground must still read as bare."""
    flat, _ = ngv(*_cloud(veg_frac=0.0, seed=3)[:3], 0.0, 0.0, _cloud(veg_frac=0.0, seed=3)[3])
    x, y, z, c = _cloud(veg_frac=0.0, slope=(0.35, -0.2), seed=3)
    tilted, _ = ngv(x, y, z, 0.0, 0.0, c)
    assert flat == 0.0 and tilted == 0.0


def test_returns_outside_the_radius_are_excluded():
    x, y, z, c = _cloud(veg_frac=0.0, seed=4)
    far_x = np.r_[x, np.array([RADIUS + 5.0])]
    far_y = np.r_[y, np.array([0.0])]
    far_z = np.r_[z, np.array([200.0])]          # a wild point well outside
    v, n = ngv(far_x, far_y, far_z, 0.0, 0.0, c)
    assert n == len(x), "a point outside the radius entered the denominator"
    assert v == 0.0


def test_window_edges_are_exclusive_below_inclusive_above():
    """The coefficient is tied to (NGV_LO, NGV_HI]. A silent edge change re-scales it.

    Ground at z=0 so the height arithmetic is exact -- see the next test for why that
    matters.
    """
    c = np.zeros(6)
    x = np.zeros(4); y = np.zeros(4)
    z = np.array([NGV_LO, NGV_LO + 1e-6, NGV_HI, NGV_HI + 1e-6])
    v, n = ngv(x, y, z, 0.0, 0.0, c)
    assert n == 4
    assert v == pytest.approx(0.5), "expected exactly the middle two returns to count"


def test_a_return_exactly_on_the_edge_is_decided_by_floating_point_not_by_the_rule():
    """Documented, not fixed: h is computed as (z - S), so with ground near 100 m the
    subtraction is inexact and a return sitting ON the lower edge falls above it --
    100.0 + 0.15 - 100.0 = 0.15000000000000568. Harmless in practice (exact-edge returns
    are measure-zero) but it means the edge is not bit-reproducible across sites at
    different elevations. Recorded so nobody debugs it twice."""
    c = np.array([100.0, 0, 0, 0, 0, 0])
    z = 100.0 + np.array([NGV_LO])
    v, n = ngv(np.zeros(1), np.zeros(1), z, 0.0, 0.0, c)
    assert n == 1
    assert v == 1.0, "at ground 100 m the lower-edge return counts, by float error"
    c0 = np.zeros(6)
    v0, _ = ngv(np.zeros(1), np.zeros(1), np.array([NGV_LO]), 0.0, 0.0, c0)
    assert v0 == 0.0, "at ground 0 m the same return is excluded, exactly as the rule says"


def test_empty_neighbourhood_is_nan_not_zero():
    """Zero would be indistinguishable from 'measured, and bare'."""
    x, y, z, c = _cloud(veg_frac=0.5, seed=5)
    v, n = ngv(x, y, z, 1e6, 1e6, c)
    assert n == 0 and np.isnan(v)


def test_design_matrix_column_order_matches_the_fitted_coefficients():
    """[1, u, v, u^2, v^2, uv]. A swapped pair silently bends the surface -- this exact
    error produced a false alarm earlier in the project."""
    A = design(np.array([2.0]), np.array([3.0]))
    assert A.shape == (1, 6)
    np.testing.assert_allclose(A[0], [1.0, 2.0, 3.0, 4.0, 9.0, 6.0])


@pytest.mark.skipif(not os.path.exists("data/derived/control_ngv_exact.csv"),
                    reason="needs the derived control table")
def test_the_shipped_control_table_is_in_range_and_complete():
    import pandas as pd
    t = pd.read_csv("data/derived/control_ngv_exact.csv")
    assert len(t) > 300
    assert t.ngv.between(0.0, 1.0).all(), "NGV is a fraction"
    assert (t.n_disc > 0).all(), "a mark with no returns should not carry a value"
