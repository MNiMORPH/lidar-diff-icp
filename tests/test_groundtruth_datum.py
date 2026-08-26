"""Combining ties into a datum constant, and keeping common-mode error at full size."""
import json

import numpy as np
import pytest

from lidar_diff_icp.groundtruth.datum import BudgetTerm, DatumConstant, combine_ties


def _terms(**kw):
    d = dict(name="t", value_mm=10.0, kind="common",
             source="measured in analysis/MISSION_TIME_DRIFT.md section 4")
    d.update(kw)
    return BudgetTerm(**d)


# --------------------------------------------------------------------------- the term

def test_budget_term_refuses_an_unknown_kind():
    with pytest.raises(ValueError, match="kind must be one of"):
        _terms(kind="probably fine")


def test_budget_term_refuses_a_term_with_no_source():
    """An error term with no provenance is a number someone typed."""
    with pytest.raises(ValueError, match="needs a source"):
        _terms(source="")
    with pytest.raises(ValueError, match="needs a source"):
        _terms(source="guess")


def test_budget_term_refuses_a_negative_or_nan_magnitude():
    with pytest.raises(ValueError, match="finite and >= 0"):
        _terms(value_mm=-1.0)
    with pytest.raises(ValueError, match="finite and >= 0"):
        _terms(value_mm=float("nan"))


# --------------------------------------------------------------------- the combination

def test_inverse_variance_weighting_pulls_towards_the_tighter_mark():
    """Elba's own case: +21.3 +/- 12.4 and +28.9 +/- 27.0."""
    d = combine_ties([("2210", 21.3, 12.4), ("2036", 28.9, 27.0)], [])
    assert 21.3 < d.value_mm < 28.9
    assert abs(d.value_mm - 21.3) < abs(d.value_mm - 28.9)      # nearer the tighter mark
    w = np.array([1 / 12.4 ** 2, 1 / 27.0 ** 2])
    assert d.value_mm == pytest.approx(float(np.dot(w, [21.3, 28.9]) / w.sum()))


def test_spread_is_reported_and_is_not_the_uncertainty():
    """The two chains agree to 7.5 mm; that is smaller than either sigma, so it must not
    become the error bar."""
    d = combine_ties([("2210", 21.3, 12.4), ("2036", 28.9, 27.0)], [])
    assert d.spread_mm == pytest.approx(7.6, abs=0.1)
    assert d.sigma_total_mm > d.spread_mm


def test_a_mark_without_a_sigma_is_kept_and_noted_not_dropped():
    d = combine_ties([("a", 10.0, 5.0), ("b", 20.0, float("nan"))], [])
    assert [p for p, _, _ in d.ties] == ["a", "b"]
    assert any("Kept, not dropped" in n for n in d.notes)


def test_combine_refuses_a_plain_number_as_a_budget_term():
    with pytest.raises(TypeError, match="must be BudgetTerm"):
        combine_ties([("a", 1.0, 1.0)], [12.4])


# ------------------------------------------------------- THE regression test (see below)

def test_common_mode_error_does_not_average_down_with_more_marks():
    """A term shared by every tie is not reduced by adding ties.

    This is the bug the module exists to prevent: the two Elba chains share one
    extrapolated lateral shift and one alignment estimator, so quoting sigma/sqrt(n)
    over all terms would shrink an error that is common to both. Reverting
    ``DatumConstant.common_mm`` to divide by the tie count (or to sqrt(k)) makes this
    fail -- which is how it bites.
    """
    common = _terms(name="lateral extrapolation", value_mm=20.0, kind="common",
                    source="docs/groundtruth.md section 7, measured per-tie effect")
    two = combine_ties([("a", 10.0, 12.0), ("b", 12.0, 12.0)], [common])
    ten = combine_ties([(f"m{i}", 10.0 + i, 12.0) for i in range(10)], [common])
    assert two.common_mm == pytest.approx(20.0)
    assert ten.common_mm == pytest.approx(20.0)                 # unchanged by n
    assert ten.sigma_total_mm >= 20.0                           # never below the floor
    # the RANDOM part does average down, so the contrast is real, not a no-op assertion
    assert ten.random_mm < two.random_mm


def test_random_terms_do_average_down():
    r = _terms(name="link error", value_mm=9.0, kind="random",
               source="chain.solve_link per-link sigma, docs/groundtruth.md section 5")
    two = combine_ties([("a", 10.0, 12.0), ("b", 12.0, 12.0)], [r])
    eight = combine_ties([(f"m{i}", 10.0, 12.0) for i in range(8)], [r])
    assert eight.random_mm < two.random_mm


def test_unmodelled_terms_are_reported_but_never_folded_into_the_total():
    """Along-track drift has a scale but no distribution; folding it into a quadrature
    sigma would dress a knowledge gap as a measurement."""
    u = _terms(name="along-track drift", value_mm=50.0, kind="unmodelled",
               source="pipeline.fit_along_track_drift measured 11-29 mm/km on elbaext")
    d = combine_ties([("a", 10.0, 12.0)], [u])
    assert d.unmodelled_mm == pytest.approx(50.0)
    assert d.sigma_total_mm < 50.0                              # not absorbed
    assert d.sigma_total_mm == pytest.approx(12.0)


# -------------------------------------------------------------------------- the outputs

def test_every_budget_row_carries_its_source_and_kind():
    d = combine_ties([("a", 10.0, 12.0)],
                     [_terms(name="x", kind="common", value_mm=5.0,
                             source="measured somewhere real")])
    cols = DatumConstant.table_columns()
    assert set(cols) == {"term", "kind", "mm", "applies_to", "source"}
    rows = d.table_rows()
    assert all(len(r) == len(cols) for r in rows)
    assert any(r[0] == "x" and r[4] == "measured somewhere real" for r in rows)
    assert any(r[0].startswith("== sigma_total") for r in rows)


def test_sidecar_records_the_sign_convention_and_round_trips(tmp_path):
    d = combine_ties([("2210", 21.3, 12.4), ("2036", 28.9, 27.0)],
                     [_terms(name="x", value_mm=12.4, kind="common",
                             source="MISSION_TIME_DRIFT.md section 4 dz repeatability")])
    p = tmp_path / "sidecar.json"
    d.to_json(p)
    got = json.loads(p.read_text())
    assert "ADD to gen1" in got["sign_convention"]
    assert got["datum_constant_mm"] == pytest.approx(d.value_mm)
    assert got["sigma_common_mm"] == pytest.approx(12.4)
    assert [t["point_id"] for t in got["ties"]] == ["2210", "2036"]
    assert got["uncertainty_budget"][0]["source"].startswith("MISSION_TIME_DRIFT")


def test_combine_needs_at_least_one_tie():
    with pytest.raises(ValueError, match="at least one tie"):
        combine_ties([], [])
