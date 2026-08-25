"""Tests for the joint slope x cover offset model (lidar_diff_icp.offset_model).

The model exists to predict an epoch-difference offset from two CORRELATED predictors, so
the tests are built around that confounding: synthetic data with a known slope-cover
correlation and known coefficients, checked for recovery, for parameterization invariance,
and for the covariance controls behaving as advertised.
"""
import numpy as np
import pytest

from lidar_diff_icp import offset_model as om

SLOPE_EDGES = np.array([0, 5, 10, 15, 20, 25, 30, 45], float)
COVER_EDGES = np.array([0.0, 0.1, 0.25, 0.5, 1.01])

TRUE = dict(a=-80.0, b=1.0, c=-10.0, e=-2.3)      # intercept, slope, cover, slope*cover


def synth(n=20000, rho=0.6, noise=5.0, seed=0):
    """Correlated (slope, cover) with an offset built from TRUE, plus noise."""
    rng = np.random.default_rng(seed)
    z1 = rng.normal(size=n); z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.normal(size=n)
    slope = np.clip(10 + 7 * z1, 0, 45)
    cover = np.clip(0.25 + 0.15 * z2, 0, 1)
    d = (TRUE["a"] + TRUE["b"] * slope + TRUE["c"] * cover + TRUE["e"] * slope * cover
         + rng.normal(scale=noise, size=n))
    return slope, cover, d


def test_predictor_covariance_reports_the_confounding():
    slope, cover, _ = synth()
    r, vif = om.predictor_covariance(slope, cover)
    assert 0.4 < r < 0.75, f"synthetic corr {r} not in the intended range"
    assert vif == pytest.approx(1 / (1 - r ** 2), rel=1e-12)
    r0, vif0 = om.predictor_covariance([0, 1, 2, 3], [1, 0, 0, 1])   # orthogonal by construction
    assert abs(r0) < 1e-12 and vif0 == pytest.approx(1.0)


def test_fit_recovers_known_coefficients():
    slope, cover, d = synth()
    m = om.fit_offset_model(slope, cover, d, interaction=True)
    assert m.coeffs[0] == pytest.approx(TRUE["a"], abs=1.0)
    assert m.coeffs[1] == pytest.approx(TRUE["b"], abs=0.05)
    assert m.coeffs[2] == pytest.approx(TRUE["c"], abs=2.0)
    assert m.coeffs[3] == pytest.approx(TRUE["e"], abs=0.1)
    # R2 here is governed by the noise budget (systematic spread ~ the 5 mm noise), not by
    # fit quality, so the meaningful check is that the residual equals the noise put in.
    assert m.rms == pytest.approx(5.0, rel=0.1)


def test_z_scored_and_physical_fits_predict_identically():
    """Regression test: the physical-unit coefficients must come from a fit in physical
    units, never from hand-converting z-scored ones. The interaction term makes that
    conversion easy to get wrong (b3*sz*cz also contributes to the slope and cover terms),
    and a version of this analysis carried an error of tens of mm from exactly that."""
    slope, cover, d = synth()
    sz = (slope - slope.mean()) / slope.std(); cz = (cover - cover.mean()) / cover.std()
    zi = om.fit_offset_model(sz, cz, d, interaction=True)
    pi = om.fit_offset_model(slope, cover, d, interaction=True)
    assert np.max(np.abs(zi.predict(sz, cz) - pi.predict(slope, cover))) < 1e-6

    # the WRONG conversion (dropping the cross terms) must be visibly different
    b = zi.coeffs; s_sd, c_sd, s_mu, c_mu = slope.std(), cover.std(), slope.mean(), cover.mean()
    k_s, k_c = b[1] / s_sd, b[2] / c_sd; k_sc = b[3] / (s_sd * c_sd)
    k_0 = b[0] - k_s * s_mu - k_c * c_mu + k_sc * s_mu * c_mu
    bad = k_0 + k_s * slope + k_c * cover + k_sc * slope * cover
    assert np.max(np.abs(bad - pi.predict(slope, cover))) > 1.0, \
        "the naive conversion should be materially wrong -- if not, this test proves nothing"


def test_interaction_beats_additive_when_the_truth_interacts():
    slope, cover, d = synth()
    add = om.fit_offset_model(slope, cover, d, interaction=False)
    inter = om.fit_offset_model(slope, cover, d, interaction=True)
    assert inter.r2 > add.r2 + 0.01
    assert inter.rms < add.rms
    assert len(add.coeffs) == 3 and len(inter.coeffs) == 4


def test_sensitivities_follow_the_interaction():
    slope, cover, d = synth()
    m = om.fit_offset_model(slope, cover, d, interaction=True)
    assert m.d_dcover(25) < m.d_dcover(5)                    # cover bites harder on steep ground
    assert m.d_dslope(0.6) < m.d_dslope(0.02)                # slope reverses under canopy
    assert m.d_dcover(0) == pytest.approx(m.coeffs[2], abs=1e-9)
    flat = om.fit_offset_model(slope, cover, d, interaction=False)
    assert flat.d_dcover(25) == pytest.approx(flat.d_dcover(5))   # additive: no slope dependence


def test_median_surface_marks_unsupported_boxes():
    slope, cover, d = synth()
    grid, cnt = om.median_surface(slope, cover, d, SLOPE_EDGES, COVER_EDGES, min_cells=30)
    assert grid.shape == (len(SLOPE_EDGES) - 1, len(COVER_EDGES) - 1) == cnt.shape
    assert np.all(np.isnan(grid[cnt < 30])), "under-populated boxes must be NaN, not a number"
    assert np.all(np.isfinite(grid[cnt >= 30]))
    assert cnt.sum() == slope.size
    # steep-and-bare is rare under a positive slope-cover correlation: some box must be empty
    assert np.isnan(grid).any(), "synthetic covariance should leave at least one box unsupported"


def test_median_surface_values_are_medians_not_means():
    slope = np.array([1.0] * 40); cover = np.array([0.02] * 40)
    d = np.concatenate([np.zeros(39), [1e6]])                # one wild outlier
    grid, cnt = om.median_surface(slope, cover, d, SLOPE_EDGES, COVER_EDGES, min_cells=10)
    assert cnt[0, 0] == 40 and grid[0, 0] == pytest.approx(0.0)


def test_matched_bands_break_the_confounding():
    """With slope and cover correlated, the MARGINAL cover gradient is contaminated by
    slope; holding slope inside a band must recover the true local cover effect."""
    slope, cover, d = synth()
    # marginal gradient = one band spanning everything, via the same code path
    marg = om.matched_band_effects(np.zeros_like(slope), cover, d, np.array([-1.0, 1.0]))[0][3]
    rows = om.matched_band_effects(slope, cover, d, SLOPE_EDGES, min_n=200)
    for lo, hi, n, grad, vmin, vmax in rows:
        if not np.isfinite(grad):
            continue
        centre = 0.5 * (lo + min(hi, 45.0))
        assert grad == pytest.approx(TRUE["c"] + TRUE["e"] * centre, abs=12.0)
    low = next(r for r in rows if r[0] == 0.0)          # 0-5 deg: local truth ~ TRUE["c"]
    assert abs(marg - TRUE["c"]) > 5.0, "the marginal gradient should be visibly confounded"
    assert abs(low[3] - TRUE["c"]) < abs(marg - TRUE["c"]), \
        "holding slope fixed must land CLOSER to the local truth than the marginal does"


def test_matched_band_reports_sparse_bands_as_nan():
    slope, cover, d = synth(n=500)
    rows = om.matched_band_effects(slope, cover, d, SLOPE_EDGES, min_n=100000)
    assert all(not np.isfinite(g) for _, _, _, g, _, _ in rows)
    assert sum(n for _, _, n, _, _, _ in rows) == 500


def test_partial_correlations_isolate_each_predictor():
    rng = np.random.default_rng(3)
    n = 8000
    slope = rng.uniform(0, 30, n)
    cover = np.clip(0.02 * slope + rng.normal(scale=0.05, size=n), 0, 1)   # strongly confounded
    d = 2.0 * slope + rng.normal(scale=1.0, size=n)                        # cover has NO effect
    p = om.partial_correlations(slope, cover, d)
    assert p["slope|cover"][0] > 0.5, "slope's own effect must survive the control"
    assert abs(p["cover|slope"][0]) < 0.2, "cover has no true effect; partial must be ~0"
    assert set(p) == {"cover|slope", "slope|cover"}
    for v in p.values():
        assert len(v) == 2 and all(np.isfinite(v))


def test_weighted_fit_follows_the_weights():
    slope = np.array([0.0, 10.0, 20.0, 30.0]); cover = np.zeros(4)
    d = np.array([0.0, 10.0, 20.0, 300.0])                  # last point is an outlier
    heavy = om.fit_offset_model(slope, cover, d, weights=[1, 1, 1, 1000], interaction=False)
    light = om.fit_offset_model(slope, cover, d, weights=[1000, 1000, 1000, 1], interaction=False)
    assert heavy.coeffs[1] > light.coeffs[1]
    assert light.coeffs[1] == pytest.approx(1.0, abs=0.05)
    assert heavy.weighted and not om.fit_offset_model(slope, cover, d).weighted


def test_cell_reduce_medians_and_min_n():
    cell = np.array([5, 5, 5, 7, 7, 9])
    val = np.array([1.0, 2.0, 300.0, 4.0, 6.0, 99.0])
    cells, med, cnt = om.cell_reduce(cell, val, min_n=2)
    assert list(cells) == [5, 7]                            # cell 9 has one return, dropped
    assert med[0] == pytest.approx(2.0) and med[1] == pytest.approx(5.0)
    assert list(cnt) == [3, 2]
