"""The control residual as a spatial field -- checked against answers known in advance.

The load-bearing checks here are the ones that would let a wrong number out quietly:

* ``test_loo_shortcut_equals_a_full_refit_*`` -- the leave-one-out cross-validation in
  ``analysis/control_residual_field.py`` uses a one-inverse identity instead of refitting
  963 times. If that identity were wrong every cross-validation number in
  ``analysis/CONTROL_RESIDUAL_FIELD.md`` would be wrong and nothing would look odd.
* ``test_kriging_is_exact_at_a_datum_when_there_is_no_nugget`` and
  ``test_kriging_reproduces_a_constant_field`` -- the two identities any kriging
  implementation must satisfy.
* ``test_the_two_prediction_sds_differ_by_exactly_the_nugget`` -- the module reports two
  standard deviations with different meanings, and the whole point of reporting both is
  that they are not interchangeable.
* ``test_every_tunable_is_required`` -- the project rule is that no threshold, bin width
  or block size has a default. This test fails the moment one acquires one.
"""
import numpy as np
import pytest

from lidar_diff_icp.groundtruth import residual_field as RF
from lidar_diff_icp.variogram import VariogramModel


@pytest.fixture(scope="module")
def control():
    return RF.load_residuals()


def _synthetic(n=40, seed=3):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 10_000, n)
    y = rng.uniform(0, 10_000, n)
    v = 50.0 + 0.004 * x - 0.002 * y + rng.normal(0, 10, n)
    return x, y, v


# ------------------------------------------------------------------ the bundled file

def test_the_bundled_control_has_963_marks_in_1004_rows(control):
    assert control.n_rows_in == 1004
    assert len(control) == 963
    assert control.n_dup_rows == 41
    assert control.n_dup_groups == 39


def test_dnr_error_is_control_minus_surface_and_not_the_other_way(control):
    n, n_cms, mx_cms, n_smc, mx_smc = RF.check_sign_convention(tol_m=1e-9)
    assert (n, n_cms, n_smc) == (1004, 1004, 0)
    assert mx_cms < 1e-9
    assert mx_smc > 1.0


def test_stratify_refuses_a_cover_class_that_does_not_exist(control):
    with pytest.raises(ValueError, match="unknown cover class"):
        RF.stratify(control, ("L1O", "L9Z"))
    assert RF.stratify(control, ("L1O",)).sum() == int((control.cover == "L1O").sum())


# ------------------------------------------------------------------------- the model

def test_the_spherical_model_is_zero_at_zero_and_flat_beyond_the_range():
    m = VariogramModel(nugget=100.0, sill=400.0, range_=5000.0)
    assert RF._gamma(0.0, m) == 0.0
    assert RF._gamma(1e-9, m) == pytest.approx(100.0, abs=1e-3)
    assert RF._gamma(5000.0, m) == pytest.approx(500.0)
    assert RF._gamma(50_000.0, m) == pytest.approx(500.0)
    h = np.linspace(1.0, 20_000.0, 50)
    assert np.all(np.diff(RF._gamma(h, m)) >= -1e-9)


def test_kriging_is_exact_at_a_datum_when_there_is_no_nugget():
    x, y, v = _synthetic()
    m = VariogramModel(nugget=0.0, sill=300.0, range_=6000.0)
    r = RF.krige(x, y, v, m, x[7], y[7])
    assert r.value_mm == pytest.approx(v[7], abs=1e-6)
    assert r.sd_new_mark_mm == pytest.approx(0.0, abs=1e-4)


def test_kriging_reproduces_a_constant_field():
    x, y, _ = _synthetic()
    v = np.full(x.size, -37.5)
    m = VariogramModel(nugget=90.0, sill=300.0, range_=4000.0)
    r = RF.krige(x, y, v, m, 4321.0, 8765.0)
    assert r.value_mm == pytest.approx(-37.5, abs=1e-8)


def test_the_two_prediction_sds_differ_by_exactly_the_nugget():
    x, y, v = _synthetic()
    m = VariogramModel(nugget=120.0, sill=300.0, range_=4000.0)
    r = RF.krige(x, y, v, m, 4321.0, 8765.0)
    assert r.sd_new_mark_mm > r.sd_field_mm
    assert r.sd_new_mark_mm ** 2 - r.sd_field_mm ** 2 == pytest.approx(120.0, rel=1e-9)


def test_krige_many_matches_krige_one_at_a_time():
    x, y, v = _synthetic()
    m = VariogramModel(nugget=50.0, sill=300.0, range_=4000.0)
    tx = np.array([1000.0, 5000.0, 9000.0])
    ty = np.array([2000.0, 5000.0, 1000.0])
    vals, vars_ = RF.krige_many(x, y, v, m, tx, ty)
    for k in range(tx.size):
        one = RF.krige(x, y, v, m, tx[k], ty[k])
        assert vals[k] == pytest.approx(one.value_mm, rel=1e-10)
        assert vars_[k] == pytest.approx(one.sd_new_mark_mm ** 2, rel=1e-10)


# ---------------------------------------------------------------- the LOO shortcut

def test_loo_shortcut_equals_a_full_refit_ordinary_kriging():
    x, y, v = _synthetic()
    m = VariogramModel(nugget=60.0, sill=300.0, range_=4000.0)
    de, dv = RF.verify_loo_shortcut(x, y, v, m, range(x.size))
    assert de < 1e-8
    assert dv < 1e-6


def test_loo_shortcut_equals_a_full_refit_with_a_cover_drift():
    x, y, v = _synthetic(n=48, seed=11)
    cover = np.array(["L1O", "L2T", "L5U"] * 16)
    v = v + np.where(cover == "L2T", -90.0, 0.0)
    X, labels, at = RF.cover_design(cover, ("L1O", "L2T", "L5U"))
    m = VariogramModel(nugget=60.0, sill=300.0, range_=4000.0)
    de, dv = RF.verify_loo_shortcut(x, y, v, m, range(x.size), X=X)
    assert de < 1e-7
    assert dv < 1e-5


def test_loo_of_a_constant_field_has_zero_error():
    x, y, _ = _synthetic()
    v = np.full(x.size, 12.0)
    m = VariogramModel(nugget=60.0, sill=300.0, range_=4000.0)
    err, var = RF.loo_errors(x, y, v, m)
    assert np.abs(err).max() < 1e-8
    assert np.all(var > 0)


# --------------------------------------------------------------------- the drift

def test_cover_design_absorbs_the_first_class_into_the_constant():
    cover = np.array(["L1O", "L5U", "L2T", "L1O"])
    X, labels, at = RF.cover_design(cover, ("L1O", "L2T", "L5U"))
    assert labels[0] == "const(=L1O)"
    assert np.allclose(X[:, 0], 1.0)
    assert np.allclose(at("L1O"), [1, 0, 0])
    assert np.allclose(at("L2T"), [1, 1, 0])
    assert np.allclose(at("L5U"), [1, 0, 1])
    assert np.linalg.matrix_rank(X) == 3


def test_a_pure_cover_offset_is_absorbed_by_the_drift_not_the_field():
    """Two classes, one offset by a known amount, no spatial structure at all.

    Kriging with the cover drift must return the reference class's level at a new
    location, not the pooled mean of the two.
    """
    x, y, _ = _synthetic(n=60, seed=5)
    cover = np.where(np.arange(60) % 2 == 0, "L1O", "L2T")
    v = np.where(cover == "L1O", 20.0, -80.0)
    X, labels, at = RF.cover_design(cover, ("L1O", "L2T"))
    m = VariogramModel(nugget=100.0, sill=200.0, range_=3000.0)
    r = RF.krige(x, y, v, m, 5000.0, 5000.0, X=X, x0_drift=at("L1O"), drift_labels=labels)
    assert r.value_mm == pytest.approx(20.0, abs=1e-6)
    r2 = RF.krige(x, y, v, m, 5000.0, 5000.0, X=X, x0_drift=at("L2T"), drift_labels=labels)
    assert r2.value_mm == pytest.approx(-80.0, abs=1e-6)


# ------------------------------------------------------------------ the null model

def test_the_constant_null_uses_only_the_training_fold():
    v = np.array([1.0, 3.0, 5.0, 7.0])
    fold = np.array([0, 0, 1, 1])
    em, ed = RF.constant_null_errors(v, fold)
    assert em[0] == pytest.approx(6.0 - 1.0)      # mean(5,7) - 1
    assert em[2] == pytest.approx(2.0 - 5.0)      # mean(1,3) - 5
    assert ed[0] == pytest.approx(6.0 - 1.0)


def test_block_cv_scores_every_mark_and_names_its_variogram_target():
    x, y, v = _synthetic(n=80, seed=9)
    err, bid, nb = RF.block_cv(x, y, v, block_m=4000.0, max_lag_m=8000.0, n_lags=8,
                               n_pairs=20_000, estimator="dowd", seed=0,
                               refit_variogram=True, variogram_on="raw")
    assert err.size == v.size
    assert np.isfinite(err).all()
    assert 1 < nb <= v.size
    with pytest.raises(ValueError, match="variogram_on"):
        RF.block_cv(x, y, v, block_m=4000.0, max_lag_m=8000.0, n_lags=8, n_pairs=20_000,
                    estimator="dowd", seed=0, refit_variogram=True, variogram_on="detrend")


def test_blocked_cv_is_harder_than_leave_one_out_on_a_smooth_field():
    """Holding out a whole block removes the near neighbours LOO leaves in place."""
    x, y, v = _synthetic(n=120, seed=13)
    m, *_ = RF.fit_field(x, y, v, max_lag_m=8000.0, n_lags=10, n_pairs=50_000,
                         estimator="dowd", seed=0)
    loo, _ = RF.loo_errors(x, y, v, m)
    blk, _, _ = RF.block_cv(x, y, v, block_m=5000.0, max_lag_m=8000.0, n_lags=10,
                            n_pairs=50_000, estimator="dowd", seed=0,
                            refit_variogram=False, variogram_on="raw")
    assert np.sqrt(np.mean(blk ** 2)) > np.sqrt(np.mean(loo ** 2))


# ------------------------------------------------------- no tunable has a default

@pytest.mark.parametrize("call", [
    lambda x, y, v: RF.fit_field(x, y, v),
    lambda x, y, v: RF.fit_field(x, y, v, max_lag_m=1000.0),
    lambda x, y, v: RF.fit_field(x, y, v, max_lag_m=1000.0, n_lags=5, n_pairs=1000),
    lambda x, y, v: RF.block_cv(x, y, v, block_m=1000.0),
    lambda x, y, v: RF.check_sign_convention(),
])
def test_every_tunable_is_required(call):
    x, y, v = _synthetic(n=12)
    with pytest.raises(TypeError):
        call(x, y, v)
