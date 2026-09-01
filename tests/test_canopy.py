"""Ground penetration: the per-cell ground fraction, and the distinction that matters.

The contract this pins is the one a stored product violated: a cell with NO returns is
NaN, and a cell with returns but none of them ground is 0.0. Those two are not the same
thing -- 0.0 says "measured, nothing reached the ground", which every downstream forest
cut (pen < 0.25) reads as maximally closed canopy. Filling the first with the second
turns absent data into a measurement.
"""
import numpy as np
import pytest

laspy = pytest.importorskip("laspy")

from lidar_diff_icp.canopy import ground_penetration, leafon_slope_flag, inflate_lod

BOUNDS = (0.0, 0.0, 20.0, 20.0)          # 4 x 4 cells at res 5
RES, NX, NY = 5.0, 4, 4


def _write_las(path, pts):
    """pts: iterable of (x, y, classification)."""
    hdr = laspy.LasHeader(version="1.2", point_format=1)
    hdr.scales = [0.001, 0.001, 0.001]
    hdr.offsets = [0.0, 0.0, 0.0]
    las = laspy.LasData(hdr)
    x, y, c = (np.asarray(a) for a in zip(*pts))
    las.x = x.astype(float)
    las.y = y.astype(float)
    las.z = np.zeros(len(x), float)
    las.classification = c.astype(np.uint8)
    las.write(str(path))
    return path


def _cell_centre(col, row):
    return BOUNDS[0] + (col + 0.5) * RES, BOUNDS[1] + (row + 0.5) * RES


def test_ground_fraction_is_ground_over_non_noise(tmp_path):
    x0, y0 = _cell_centre(0, 0)
    pts = [(x0, y0, 2)] * 3 + [(x0, y0, 5)]          # 3 ground of 4 -> 0.75
    p = ground_penetration(_write_las(tmp_path / "a.laz", pts), BOUNDS, RES, NX, NY)
    assert p.shape == (NY, NX)
    assert p[0, 0] == pytest.approx(0.75)


def test_a_cell_with_no_returns_is_nan_not_zero(tmp_path):
    """The defect this exists to catch: absent must not read as measured-and-empty."""
    x0, y0 = _cell_centre(0, 0)
    x1, y1 = _cell_centre(2, 1)
    pts = [(x0, y0, 2), (x0, y0, 5),                 # cell (0,0) has returns
           (x1, y1, 5), (x1, y1, 5)]                 # cell (2,1) has returns, NO ground
    p = ground_penetration(_write_las(tmp_path / "b.laz", pts), BOUNDS, RES, NX, NY)

    assert p[1, 2] == 0.0, "returns present but none ground -> a measured 0.0"
    empty = np.isnan(p)
    assert empty.sum() == NX * NY - 2, "every cell without returns must be NaN"
    assert not np.isnan(p[0, 0]) and not np.isnan(p[1, 2])

    # and the two states must be distinguishable by the cut every consumer uses
    forest = p < 0.25                                 # NaN compares False, which is correct
    assert forest[1, 2], "a genuinely bare-of-ground cell is forest by this cut"
    assert not forest[0, 1], "an EMPTY cell must not be swept into forest"


def test_noise_is_excluded_from_both_numerator_and_denominator(tmp_path):
    x0, y0 = _cell_centre(1, 1)
    pts = [(x0, y0, 2), (x0, y0, 5), (x0, y0, 7), (x0, y0, 7)]
    p = ground_penetration(_write_las(tmp_path / "c.laz", pts), BOUNDS, RES, NX, NY)
    assert p[1, 1] == pytest.approx(0.5), "1 ground of 2 non-noise, not 1 of 4"


def test_a_cell_of_pure_noise_has_no_measurement(tmp_path):
    x0, y0 = _cell_centre(3, 3)
    xk, yk = _cell_centre(0, 0)
    p = ground_penetration(_write_las(tmp_path / "d.laz",
                                      [(x0, y0, 7), (x0, y0, 7), (xk, yk, 2)]),
                           BOUNDS, RES, NX, NY)
    assert np.isnan(p[3, 3]), "noise-only cell is unmeasured, not zero-penetration"


def test_points_outside_the_bounds_are_dropped(tmp_path):
    xin, yin = _cell_centre(0, 0)
    pts = [(xin, yin, 2), (xin, yin, 5),
           (-1.0, 5.0, 2), (25.0, 5.0, 2), (5.0, -3.0, 2), (5.0, 40.0, 2),
           (20.0, 5.0, 2)]                            # the far edge is exclusive
    p = ground_penetration(_write_las(tmp_path / "e.laz", pts), BOUNDS, RES, NX, NY)
    assert p[0, 0] == pytest.approx(0.5)
    assert np.isnan(p).sum() == NX * NY - 1


def test_leafon_flag_never_fires_on_an_unmeasured_cell():
    pen = np.array([[np.nan, 0.10], [0.10, 0.90]])
    slope = np.full((2, 2), 30.0)
    flag = leafon_slope_flag(pen, slope)
    assert not flag[0, 0], "NaN penetration is not evidence of closed canopy"
    assert flag[0, 1] and flag[1, 0]
    assert not flag[1, 1], "open ground is not flagged however steep"


def test_leafon_flag_needs_both_low_penetration_and_slope():
    pen = np.full((2, 2), 0.10)
    slope = np.array([[0.0, 30.0], [11.9, 12.1]])
    flag = leafon_slope_flag(pen, slope)
    assert list(flag.ravel()) == [False, True, False, True]


def test_inflate_lod_widens_only_flagged_cells_and_does_not_mutate():
    lod = np.array([[10.0, 20.0], [30.0, 40.0]])
    before = lod.copy()
    flag = np.array([[True, False], [False, True]])
    out = inflate_lod(lod, flag, factor=3.0)
    assert np.array_equal(lod, before), "must return a copy"
    assert np.array_equal(out, np.array([[30.0, 20.0], [30.0, 120.0]]))
