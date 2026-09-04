"""The calibrated ground percentile: the curve, and the refusals that keep it honest.

Every test here pins a behaviour that, if it silently reverted, would put an unmeasured
bias into elevations rather than raise. That is the failure mode this route exists to avoid,
so the refusals are tested as carefully as the arithmetic.
"""
import os

import numpy as np
import pytest

from lidar_diff_icp import groundq


def _write_curve(tmp_path, epoch="gen2_2021_control", log_sd=None, q=None,
                 point_types="NVA"):
    p = tmp_path / f"ground_q_vs_class2sd_{epoch}_{point_types}.npz"
    log_sd = np.log([20.0, 60.0, 120.0, 400.0]) if log_sd is None else log_sd
    q = np.array([0.506, 0.445, 0.395, 0.101]) if q is None else q
    kw = {} if point_types is None else {"point_types": point_types}
    np.savez(p, log_sd_mm=log_sd, q=q, n_marks=519, set=epoch,
             fitted_on="test", response="test", covariate="test",
             shape="isotonic", cv="test", known_limits="test", **kw)
    return str(p)


def test_curve_without_point_types_is_refused(tmp_path):
    """NVA, VVA and LCP are three different populations -- the class-2 median sits -3.5 mm
    from truth at NVA, +103.3 mm at VVA (sited under vegetation BY DESIGN) and -23.1 mm at
    LCP. A pooled curve's falling limb is the VVA marks, so a curve that cannot say what
    went into it cannot be read as a vegetation correction for ordinary ground."""
    p = _write_curve(tmp_path, point_types=None)
    with pytest.raises(ValueError, match="records no point_types"):
        groundq.load_curve(p)


def test_difference_dem_refuses_calibrated_without_a_named_curve():
    """There is no default curve, on purpose: on open ground the calibrated curve measured
    WORSE than ground_q = 0.50 (RMS 52.5 vs 49.1 mm, held out), so it must be asked for by
    name rather than arrived at by default."""
    from lidar_diff_icp.pipeline import difference_dem
    with pytest.raises(ValueError, match="requires gen2_curve"):
        difference_dem("nonexistent_before.laz", "nonexistent_after.laz",
                       (0.0, 0.0, 100.0, 100.0), ground_q="calibrated", stream=True)


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
    renamed = str(tmp_path / "ground_q_vs_class2sd_gen2_2021_control_NVA.npz")
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


def test_difference_dem_refuses_calibrated_on_the_in_memory_path(tmp_path):
    """The in-memory path grids with pandas groupby.quantile, which takes ONE quantile for
    all cells. Refusing beats falling back to 0.50 silently."""
    from lidar_diff_icp.pipeline import difference_dem
    with pytest.raises(ValueError, match="needs stream=True"):
        difference_dem("nonexistent_before.laz", "nonexistent_after.laz",
                       (0.0, 0.0, 100.0, 100.0), ground_q="calibrated",
                       gen2_curve=_write_curve(tmp_path), stream=False)


def test_difference_dem_rejects_an_unknown_ground_q_string():
    from lidar_diff_icp.pipeline import difference_dem
    with pytest.raises(ValueError, match="must be a float or 'calibrated'"):
        difference_dem("nonexistent_before.laz", "nonexistent_after.laz",
                       (0.0, 0.0, 100.0, 100.0), ground_q="median", stream=True)


def test_stream_ground_refuses_an_unresolved_ground_q_string():
    """difference_dem resolves 'calibrated' into a scalar + curve before calling. A string
    reaching _stream_ground means a call site was missed -- which happened: one of two was
    patched and the other died 300 lines away with a numpy ufunc error naming no cause."""
    from lidar_diff_icp.pipeline import _stream_ground
    with pytest.raises(TypeError, match="this call site was missed"):
        _stream_ground("nonexistent.laz", (0.0, 0.0, 10.0, 10.0), 5.0, 2, 2, "calibrated")


def _one_cell_histogram(hg, zlo=-1.0, zhi=2.0, dz=0.02):
    nz = int(round((zhi - zlo) / dz))
    idx = np.floor((np.asarray(hg) - zlo) / dz).astype(int)
    keep = (idx >= 0) & (idx < nz)
    H = np.zeros((1, nz), np.int32)
    np.add.at(H, (np.zeros(int(keep.sum()), int), idx[keep]), 1)
    return H, zlo, dz


def _ground_and_mat(seed=0):
    """A cell that looks like a real one: a tight ground return plus a vegetation mat."""
    rng = np.random.default_rng(seed)
    return np.r_[rng.normal(0.0, 0.05, 400), rng.normal(0.35, 0.20, 200)]


def test_tile_spread_is_the_same_statistic_the_curve_was_calibrated_on():
    """THE GUARD THIS MODULE EXISTS FOR. The curve is indexed by the plain standard
    deviation of the ground-class column, measured at the marks by mark_statistics. On a
    tile the same number has to come out of the histogram, or the curve is being read at
    the wrong place on its own x-axis."""
    hg = _ground_and_mat()
    H, zlo, dz = _one_cell_histogram(hg)
    sd_hist, count = groundq.spread_from_histogram(H, zlo, dz, min_count=20)
    sd_marks = groundq.mark_statistics(hg, 0.0)["sd"] * 1000.0
    assert count[0] == hg.size
    assert abs(float(sd_hist[0]) - sd_marks) < 3.0, (float(sd_hist[0]), sd_marks)


def test_a_robust_spread_is_a_different_number_and_would_misread_the_curve():
    """Why the test above is not trivially satisfied by any spread. The pipeline briefly
    indexed the curve with 1.4826*(p75-p25)/1.349 -- a robust spread of a residual to a
    different plane -- and on Whitewater that pushed the ground down by ~1 m in the worst
    tenth of cells. On a ground-plus-mat column the two disagree by more than the whole
    calibrated range is wide."""
    hg = _ground_and_mat()
    plain = groundq.mark_statistics(hg, 0.0)["sd"] * 1000.0
    robust = 1.4826 * (np.percentile(hg, 75) - np.percentile(hg, 25)) / 1.349 * 1000.0
    assert abs(robust - plain) > 20.0, (robust, plain)


def test_ground_at_q_reads_the_column_like_a_quantile():
    hg = _ground_and_mat()
    H, zlo, dz = _one_cell_histogram(hg)
    for q in (0.10, 0.25, 0.50, 0.75):
        got = float(groundq.ground_at_q(H, zlo, dz, np.array([q]))[0]) / 1000.0
        assert abs(got - float(np.quantile(hg, q))) < dz, (q, got)


def test_ground_at_q_propagates_a_declined_cell_as_nan():
    """A cell the curve declined to estimate must not come back looking corrected."""
    H, zlo, dz = _one_cell_histogram(_ground_and_mat())
    assert np.isnan(groundq.ground_at_q(H, zlo, dz, np.array([np.nan]))[0])


def test_spatial_folds_split_whole_blocks_not_marks():
    """Neighbouring marks share a flight line and a phenology; splitting them at random
    would let the curve be scored on a mark it effectively already saw."""
    e = np.array([100., 200., 10_100., 10_200., 20_100., 20_200.])
    n = np.zeros(6)
    f, blocks = groundq.spatial_folds(e, n, n_folds=3)
    assert len(blocks) == 3
    assert f[0] == f[1] and f[2] == f[3] and f[4] == f[5]
    assert len(set(f.tolist())) == 3
