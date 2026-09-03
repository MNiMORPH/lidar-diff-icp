"""The calibrated ground percentile: the curve, and the refusals that keep it honest.

Every test here pins a behaviour that, if it silently reverted, would put an unmeasured
bias into elevations rather than raise. That is the failure mode this route exists to avoid,
so the refusals are tested as carefully as the arithmetic.
"""
import os

import numpy as np
import pytest

from lidar_diff_icp import groundq


def _write_curve(tmp_path, epoch="gen2_2021_control", log_sd=None, q=None):
    p = tmp_path / f"ground_q_vs_class2sd_{epoch}.npz"
    log_sd = np.log([20.0, 60.0, 120.0, 400.0]) if log_sd is None else log_sd
    q = np.array([0.506, 0.445, 0.395, 0.101]) if q is None else q
    np.savez(p, log_sd_mm=log_sd, q=q, n_marks=519, set=epoch,
             fitted_on="test", response="test", covariate="test",
             shape="isotonic", cv="test", known_limits="test")
    return str(p)


def test_curve_refuses_an_epoch_it_was_not_calibrated_on(tmp_path):
    """A gen1 curve must never be applied to gen2. 2008 is leaf-off November, 2021 is
    green-up May, and the deliveries used different classifiers -- a silent mismatch would
    bias every elevation on the tile."""
    p = _write_curve(tmp_path, epoch="gen1_2008_control")
    with pytest.raises(ValueError, match="valid only for its own epoch"):
        groundq.load_curve(p, expect_epoch="gen2_2021_control")
    # and the matching case loads
    assert groundq.load_curve(p, expect_epoch="gen1_2008_control")["epoch"] == "gen1_2008_control"


def test_curve_checks_the_epoch_INSIDE_the_file_not_the_filename(tmp_path):
    """Renaming a file must not smuggle the wrong calibration through."""
    p = _write_curve(tmp_path, epoch="gen1_2008_control")
    renamed = str(tmp_path / "ground_q_vs_class2sd_gen2_2021_control.npz")
    os.rename(p, renamed)
    with pytest.raises(ValueError, match="calibrated on 'gen1_2008_control'"):
        groundq.load_curve(renamed, expect_epoch="gen2_2021_control")


def test_missing_curve_refuses_rather_than_defaulting(tmp_path):
    with pytest.raises(FileNotFoundError, match="no default curve and no fallback"):
        groundq.load_curve(str(tmp_path / "absent.npz"))


def test_q_interpolates_in_log_spread(tmp_path):
    c = groundq.load_curve(_write_curve(tmp_path))
    # a spread exactly at a knot returns that knot's q
    assert groundq.q_from_spread([60.0], c)[0] == pytest.approx(0.445, abs=1e-9)
    # between knots it interpolates in LOG spread, not linear spread
    mid_log = float(np.exp(0.5 * (np.log(60.0) + np.log(120.0))))
    assert groundq.q_from_spread([mid_log], c)[0] == pytest.approx(0.5 * (0.445 + 0.395),
                                                                   abs=1e-9)


def test_q_is_held_at_the_ends_and_never_leaves_zero_one(tmp_path):
    """Beyond the calibrated range the curve is HELD, not extrapolated. A linear
    continuation ran q out of [0, 1] on 8-11% of cells -- the failure that dogged every
    cover-relation version of this correction."""
    c = groundq.load_curve(_write_curve(tmp_path))
    q = groundq.q_from_spread([1.0, 1e6], c)
    assert q[0] == pytest.approx(0.506)      # below the first knot: held
    assert q[1] == pytest.approx(0.101)      # above the last: held
    assert np.all((q >= 0) & (q <= 1))


def test_q_declines_to_estimate_rather_than_defaulting(tmp_path):
    """NaN means 'this method will not estimate here'. It is deliberately NOT 0.50: a cell
    that was never corrected must not end up looking corrected."""
    c = groundq.load_curve(_write_curve(tmp_path))
    q = groundq.q_from_spread([np.nan, 0.0, -5.0, 60.0], c)
    assert np.isnan(q[:3]).all()
    assert np.isfinite(q[3])
    # and the count gate
    q2 = groundq.q_from_spread([60.0, 60.0], c, min_count=20, count=[5, 50])
    assert np.isnan(q2[0]) and np.isfinite(q2[1])


def test_curve_is_monotone_non_increasing_as_calibrated(tmp_path):
    """More contamination cannot mean a HIGHER ground rank. The isotonic fit enforces it;
    this pins that the loader does not reorder or resample it away."""
    c = groundq.load_curve(_write_curve(tmp_path))
    assert np.all(np.diff(c["q"]) <= 0)
    sd = np.array([10.0, 30.0, 60.0, 100.0, 200.0, 500.0])
    assert np.all(np.diff(groundq.q_from_spread(sd, c)) <= 1e-12)


def test_difference_dem_refuses_calibrated_on_the_in_memory_path():
    """The in-memory path grids with pandas groupby.quantile, which takes ONE quantile for
    all cells. Refusing beats falling back to 0.50 silently."""
    from lidar_diff_icp.pipeline import difference_dem
    with pytest.raises(ValueError, match="needs stream=True"):
        difference_dem("nonexistent_before.laz", "nonexistent_after.laz",
                       (0.0, 0.0, 100.0, 100.0), stream=False)


def test_difference_dem_rejects_an_unknown_ground_q_string():
    from lidar_diff_icp.pipeline import difference_dem
    with pytest.raises(ValueError, match="must be a float or 'calibrated'"):
        difference_dem("nonexistent_before.laz", "nonexistent_after.laz",
                       (0.0, 0.0, 100.0, 100.0), ground_q="median", stream=True)
