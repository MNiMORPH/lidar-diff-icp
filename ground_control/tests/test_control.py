"""Tests for the epoch-agnostic control adapter.

Each test pins a fact that was re-derived from the bundled CSVs, not quoted from a
report.  Two of them exist specifically to keep a known trap from reopening.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ground_control"))
sys.path.insert(0, str(ROOT / "src"))

import control  # noqa: E402


# ------------------------------------------------------------- sign convention

def test_gen1_residual_is_surveyed_minus_lidar_on_every_row():
    (c,) = control.verify_sign_convention("gen1", tol_m=1e-9)
    assert c.n_rows_checked == 1004
    assert c.n_exact_surveyed_minus_lidar == 1004
    assert c.n_exact_lidar_minus_surveyed == 0
    # the reverse ordering is not a near-miss; it is wrong by a metre
    assert c.worst_miss_lidar_minus_surveyed_m > 1.0
    assert c.is_surveyed_minus_lidar


@pytest.mark.parametrize("surface,n", [("ql1_dem", 238), ("ql1_laz", 238),
                                       ("ql0_dem", 157), ("ql0_laz", 157)])
def test_gen2_every_surface_is_surveyed_minus_lidar(surface, n):
    (c,) = control.verify_sign_convention("gen2", tol_m=1e-9, surface=surface)
    assert c.n_rows_checked == n
    assert c.n_exact_surveyed_minus_lidar == n
    assert c.n_exact_lidar_minus_surveyed == 0
    assert c.is_surveyed_minus_lidar


# ------------------------------------------------------------- de-duplication

def test_gen1_dedup_1004_rows_to_963_marks():
    L = control.load_control("gen1")
    assert L.n_rows_in == 1004
    assert L.n_marks_out == 963
    assert L.n_dup_rows == 41
    assert L.n_dup_groups == 39


# --------------------------------------------------- the L1O / L10 trap (regression)

def test_open_marks_come_from_the_column_not_the_point_id_prefix():
    """21 gen1 marks have a point_id starting 'L10' with a DIGIT ZERO.

    Reading cover from the point_id prefix silently drops them: 209 open marks instead
    of the CSV's own 230.  This test fails if load_control is ever changed to parse the
    prefix.  Proven to bite -- see REPORT.md.
    """
    L = control.load_control("gen1")
    r = L.residuals
    from_column = int((r.cover == "L1O").sum())
    from_prefix = int(sum(1 for p in r.point_id
                          if control.cover_from_point_id_prefix(p) == "L1O"))
    digit_zero = int(sum(1 for p in r.point_id
                         if control.cover_from_point_id_prefix(p) == "L10"))
    assert from_column == 230, "the CSV's point_type column is the correct source"
    assert from_prefix == 209, "prefix parsing loses marks"
    assert digit_zero == 21
    assert from_column - from_prefix == digit_zero


# ------------------------------------------------- gen2 has four surfaces, not one

def test_gen2_requires_an_explicit_surface():
    with pytest.raises(ValueError, match="requires surface"):
        control.load_control("gen2")


def test_gen1_rejects_a_surface_argument():
    with pytest.raises(ValueError, match="ONE surface"):
        control.load_control("gen1", surface="ql1_dem")


def test_the_four_gen2_surfaces_are_not_interchangeable():
    """Different marks, different counts, different answers -- so `surface` is a choice."""
    got = {s: control.load_control("gen2", surface=s) for s in control.GEN2_SURFACES}
    assert got["ql1_dem"].n_marks_out == 238
    assert got["ql0_dem"].n_marks_out == 157
    means = {s: float(np.mean(L.residuals.resid_mm)) for s, L in got.items()}
    assert max(means.values()) - min(means.values()) > 5.0, means


# ------------------------------------------------------------------ the LCP fact

def test_lcps_carry_no_residual_and_so_never_enter_a_field():
    """All 143 LCPs are absent from every surface -- a fact of the table, not a filter."""
    for s in control.GEN2_SURFACES:
        L = control.load_control("gen2", surface=s)
        assert "calibration" not in L.roles_present, s
        assert set(L.roles_present) == {"check"}, s
        assert "LCP" not in L.covers_present, s


def test_control_residuals_is_the_shape_residual_field_already_consumes():
    """The adapter must hand `residual_field`'s own dataclass to its own estimators."""
    from lidar_diff_icp.groundtruth.residual_field import ControlResiduals, stratify
    L = control.load_control("gen1")
    assert isinstance(L.residuals, ControlResiduals)
    m = stratify(L.residuals, ("L1O",))
    assert int(np.asarray(m).sum()) == 230
