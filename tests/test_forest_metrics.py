"""The ground-cover filter's mask logic: it is reused by ten scripts and had no tests.

These cover the parts that can silently go wrong without anyone noticing: the two masks
are not a partition, an inverted threshold pair would make a cell both forest and open,
and NaN cover must belong to neither.
"""
import sys, os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from forest_metrics_pfs import classify, report


def test_masks_are_not_a_partition_and_the_accounting_closes():
    cover = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.49, 0.5, 0.9])
    forest, openg, a = classify(cover, 0.5, 0.1)
    assert forest.sum() == 2                      # 0.5 and 0.9
    assert openg.sum() == 3                       # 0.0, 0.05, 0.1 (inclusive)
    assert a["n_between"] == 3                    # 0.2, 0.3, 0.49 belong to NEITHER
    assert a["n_forest"] + a["n_open"] + a["n_between"] == a["n_cover"]
    assert not (forest & openg).any(), "a cell must not be both forest and open"


def test_nan_cover_is_in_neither_mask_and_is_counted_separately():
    cover = np.array([np.nan, 0.0, np.nan, 0.9])
    forest, openg, a = classify(cover, 0.5, 0.1)
    assert not forest[0] and not openg[0]
    assert a["n_nocover"] == 2
    assert a["n_cover"] == 2
    assert a["n_cells"] == 4


def test_inverted_thresholds_raise_rather_than_overlap():
    """open_cover above forest_cover would make mid-cover cells BOTH, which downstream
    code has no way to detect. It must fail loudly instead."""
    with pytest.raises(ValueError, match="must be below"):
        classify(np.array([0.0, 0.3, 0.9]), 0.1, 0.5)


def test_equal_thresholds_also_raise():
    with pytest.raises(ValueError):
        classify(np.array([0.0, 0.5, 0.9]), 0.5, 0.5)


def test_boundaries_are_inclusive_on_both_sides():
    forest, openg, _ = classify(np.array([0.1, 0.5]), 0.5, 0.1)
    assert openg[0], "cover == open_cover must count as open"
    assert forest[1], "cover == forest_cover must count as forest"


def test_report_names_the_unclassified_fraction_and_the_two_caveats():
    _, _, a = classify(np.array([0.0, 0.3, 0.9]), 0.5, 0.1)
    txt = report(a, 0.5, 0.1, np.array([0.0, 0.3, 0.9]))
    assert "NEITHER" in txt, "the unclassified fraction must be stated, not dropped"
    assert "not calibrated" in txt, "the thresholds must be declared as uncalibrated"
    assert "blind to undergrowth" in txt, "the canopy/undergrowth limit must travel with it"


def test_a_tile_of_all_mid_cover_classifies_nothing_and_says_so():
    """The Elba case in miniature: most cells between the thresholds. The filter must
    report that it classified nothing rather than return two empty masks quietly."""
    cover = np.full(100, 0.3)
    forest, openg, a = classify(cover, 0.5, 0.1)
    assert forest.sum() == 0 and openg.sum() == 0
    assert a["n_between"] == 100
    assert "100" in report(a, 0.5, 0.1, cover)
