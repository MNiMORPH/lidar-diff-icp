"""The reference-project resolver must never downgrade the epoch silently.

Measured on the mnrv tile 2026-09-04: MN_SEDriftless_5_2021 misses 0.0064% of the bbox --
a ~50 x 20 m sliver at one corner -- and the resolver returned
USGS_LPC_MN_Phase1_LeSueurCO_2010_LAS_2016 instead, with no warning. A 2010 acquisition
offered as the gen2 half of a 2008-to-2021 difference is a different experiment, and
nothing downstream would have said so.
"""
import pytest
from shapely.geometry import box

from lidar_diff_icp import threedep


def test_bbox_covered_is_true_only_when_nothing_is_left_out():
    """The predicate itself was never the bug -- shapely's `contains` is already True for an
    identical box, so switching to `covers` changed nothing about the mnrv case. Recorded
    here so the next reader does not credit the wrong fix. What matters is that a REAL
    shortfall, however small, still reads as not-covered; the mnrv sliver was 0.0064% of the
    bbox."""
    bb = (0.0, 0.0, 1.0, 1.0)
    assert threedep.bbox_covered(box(0, 0, 1, 1), bb)          # exact
    assert threedep.bbox_covered(box(-1, -1, 2, 2), bb)        # roomy
    assert not threedep.bbox_covered(box(0, 0, 1, 0.999), bb)  # a sliver missing


def test_resolver_refuses_to_fall_back_to_an_older_project(monkeypatch):
    newest = dict(name="MN_New_2021", url="u1", latest=2021, is_mosaic=False,
                  geom=box(0, 0, 1, 0.999))          # misses a sliver at the top
    older = dict(name="MN_Old_2010_LAS_2016", url="u2", latest=2016, is_mosaic=False,
                 geom=box(-1, -1, 2, 2))             # covers everything
    monkeypatch.setattr(threedep, "find_projects", lambda *a, **k: [newest, older])
    with pytest.raises(LookupError, match="Refusing to downgrade the epoch"):
        threedep.resolve_reference(0.5, 0.5, (0.0, 0.0, 1.0, 1.0))


def test_resolver_returns_the_newest_when_it_does_cover(monkeypatch):
    newest = dict(name="MN_New_2021", url="u1", latest=2021, is_mosaic=False,
                  geom=box(-1, -1, 2, 2))
    older = dict(name="MN_Old_2010_LAS_2016", url="u2", latest=2016, is_mosaic=False,
                 geom=box(-1, -1, 2, 2))
    monkeypatch.setattr(threedep, "find_projects", lambda *a, **k: [newest, older])
    got = threedep.resolve_reference(0.5, 0.5, (0.0, 0.0, 1.0, 1.0))
    assert got["name"] == "MN_New_2021" and got["covers"] is True
