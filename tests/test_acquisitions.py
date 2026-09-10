"""gen1 is four surveys on two geoid models, and the frame is never defaulted.

These exist because the failure they cover was real and shipped. Until 2026-09-07
``references.geoid_difference`` defaulted ``before_geoid`` to GEOID03 and
``pipeline.difference_dem`` never overrode it, so three of the six sites were differenced
in the wrong vertical frame: +54.87 mm at battlecreek, +27.75 at carlton, +26.39 at cook.
Nothing inside a tile reveals that -- a wrong geoid does not look like an error, it looks
like erosion.
"""
import inspect

import pytest

from lidar_diff_icp import acquisitions as A
from lidar_diff_icp import references
from lidar_diff_icp.pipeline import apply_datum
from lidar_diff_icp.sites import SITES


def test_every_site_names_a_survey_we_have_read_the_metadata_for():
    """A site with no survey cannot be given a frame, so it must not silently get one."""
    for name, s in SITES.items():
        assert s.gen1_project, f"{name} names no gen1 survey"
        acq = A.for_project(s.gen1_project)          # raises if unrecorded
        assert acq.geoid_grid                        # raises if the model has no grid
        assert acq.datum_quote and acq.metadata_page, (
            f"{name}: the geoid is a DATASET-level assertion and must carry the page and "
            f"the sentence that asserts it, or it cannot be checked")


def test_the_three_geoid09_sites_are_recorded_as_such():
    """The specific fact the bug turned on. If these ever flip back to GEOID03 silently,
    the DoD moves by tens of mm at three sites with nothing else changing."""
    assert A.for_project(SITES["battlecreek"].gen1_project).geoid_model == "GEOID09"
    assert A.for_project(SITES["carlton"].gen1_project).geoid_model == "GEOID09"
    assert A.for_project(SITES["cook"].gen1_project).geoid_model == "GEOID09"
    assert A.for_project(SITES["elba"].gen1_project).geoid_model == "GEOID03"
    assert A.for_project(SITES["whitewater"].gen1_project).geoid_model == "GEOID03"
    assert A.for_project(SITES["mnrv"].gen1_project).geoid_model == "GEOID03"


def test_geoid_difference_has_no_before_default():
    """THE REGRESSION. A default here is what made the bug invisible: every caller that
    forgot the argument got GEOID03 and no error. Restoring the default makes this fail."""
    p = inspect.signature(references.geoid_difference).parameters["before_geoid"]
    assert p.default is inspect.Parameter.empty, (
        "before_geoid has a default again. It must not: a caller that omits it would "
        "silently get one survey's frame applied to another's data.")


def test_apply_datum_refuses_rather_than_defaulting_the_frame():
    """THE REGRESSION, at the layer that shipped the wrong numbers. apply_datum is the
    only stage that reads across the epochs, so it is where the frame is chosen. The
    refusal must come BEFORE any work -- passing junk arrays proves nothing is touched."""
    with pytest.raises(ValueError, match="gen1_geoid is required"):
        apply_datum(None, None, None, None, None, None, None,
                    (0.0, 0.0, 100.0, 100.0))


def test_an_unknown_survey_raises_and_says_what_to_do():
    with pytest.raises(KeyError, match="no acquisition recorded"):
        A.for_project("lidar_not_a_survey_2099")
    with pytest.raises(KeyError, match="no acquisition recorded for county"):
        A.for_county("nowhere")


def test_county_lookup_agrees_with_the_project_lookup():
    """The two indexes are built from one table; if they ever disagree, one of them is
    reading a stale copy."""
    for county, acq in A._BY_COUNTY.items():
        assert A.for_county(county) is A.for_project(acq.project_id)


def test_every_site_resolves_its_own_survey_control():
    """A Site's marks come from ITS survey, never from a default.

    Both loaders used to mean the SE-Minnesota 2008 table: gen1_datum.load_control()
    resolves DEFAULT_CONTROL at 18 bare call sites, and residual_field.GEN1_CSV hard-codes
    the same file. So four of the six sites' marks -- bundled 2026-09-09 -- were reachable
    but never reached. Same shape as the geoid default that cost +54.87 mm at Ramsey.
    """
    from lidar_diff_icp.groundtruth import gen1_datum as G
    seen = {}
    for name, s in SITES.items():
        acq = A.for_project(s.gen1_project)
        assert acq.control_set, f"{name}: no control set recorded for {s.gen1_project}"
        cs = G.control_for_survey(s.gen1_project)
        assert len(cs) > 0
        assert all(m.checkpoint.project_id == s.gen1_project for m in cs.marks), (
            f"{name}: control_for_survey returned another survey's marks")
        seen[name] = len(cs)
    # the four non-SE-MN sites must NOT be getting the SE-MN table
    semn = seen["elba"]
    for name in ("mnrv", "cook", "carlton", "battlecreek"):
        assert seen[name] != semn, f"{name} is still being handed the SE-MN control"


def test_a_point_id_is_unique_only_within_a_survey():
    """THE REGRESSION. Two bundled files carry TWO surveys each, and a point_id repeats
    across them: L1O-1108 is a lesueur monument AND a ramsey one, at different positions.

    load_control's uniqueness guard is right to refuse that, so the survey filter must be
    applied BEFORE it. Filtering afterwards can never run -- the load raises first. This
    fails if `project` is ever moved to a post-load filter.
    """
    from lidar_diff_icp.groundtruth import gen1_datum as G
    with pytest.raises(ValueError, match="two different positions"):
        G.load_control("mn_dnr_control_swmn2010_metro2011")          # whole file
    for proj, n in (("lidar_swmn2010", 100), ("lidar_metro2011", 108)):
        assert len(G.load_control("mn_dnr_control_swmn2010_metro2011",
                                  project=proj)) == n


def test_an_unrecorded_control_set_refuses_rather_than_substituting():
    from lidar_diff_icp.groundtruth import gen1_datum as G
    with pytest.raises(KeyError, match="no acquisition recorded"):
        G.control_for_survey("lidar_not_a_survey_2099")


def test_control_reach_counts_and_names_without_judging():
    """The gen1 mirror of `completeness`. It exists because mnrv built completely, passed
    every other step, and could not carry a datum -- zero marks of any cover inside its
    single tile -- with nothing in the graph saying so.

    REPORTED, NEVER JUDGED: no threshold anywhere, and marks excluded by the cover rule
    are COUNTED in a note rather than silently dropped.
    """
    from lidar_diff_icp import control_reach as CR
    r = CR.reach_for_site(SITES["elba"])
    assert r.project_id == "lidar_semn2008"
    # independently cross-checked: same_line.marks_in_tiles finds 29 open marks there
    assert r.marks_reachable == 29
    assert r.n_marks_considered <= r.n_marks_in_survey
    assert any("excluded by cover" in n for n in r.notes)
    # mnrv is the case that motivated it. Assert the INVARIANT, not the count: this test
    # first asserted marks_reachable == 0, which the very next act -- fetching mnrv's 19
    # tiles -- made false. A test pinned to how much data happens to be on disk fails for
    # the wrong reason and teaches nothing.
    m = CR.reach_for_site(SITES["mnrv"])
    assert m.project_id == "lidar_swmn2010"
    assert m.marks_reachable + len(m.missing) >= m.n_marks_in_radius or not m.missing, (
        "every considered mark is either reachable or accounted for by a missing tile")
    assert (m.marks_reachable == m.n_marks_in_radius) == (not m.missing), (
        "nothing missing iff every mark is reachable")


def test_control_reach_uses_the_site_s_own_survey():
    """A mark from the wrong survey is not a weaker tie, it is a different datum."""
    from lidar_diff_icp import control_reach as CR
    for name in ("cook", "carlton", "battlecreek"):
        r = CR.reach_for_site(SITES[name])
        assert r.project_id == SITES[name].gen1_project
        assert A.for_project(r.project_id).geoid_model == "GEOID09"


def test_the_swath_cache_key_tracks_the_alignment_CODE():
    """THE REGRESSION. The key was tile|res|tie|exclude -- nothing about the METHOD.

    On 2026-09-10 that let a bridge run mix constants from two alignments: entries cached
    before commit 9e78f4a (swath network weighted by cell count) served alongside entries
    recomputed after it (1/variance). Four of 29 marks moved, one by 45.1 mm, and the
    bridge mean went -4.04 -> -5.52 mm with nothing announcing a method change.

    Editing coreg must change the key, or a stale constant is served silently again.
    """
    import inspect
    from lidar_diff_icp import coreg
    from lidar_diff_icp.groundtruth import reconstruct as R

    k = R._alignment_code_key()
    assert k and len(k) == 12

    real = inspect.getsource
    try:                                   # a changed coreg must give a changed key
        inspect.getsource = lambda o: real(o) + "\n# an edit to the alignment code\n"
        assert R._alignment_code_key() != k, (
            "the cache key does not track coreg; a superseded alignment would be served")
    finally:
        inspect.getsource = real
    assert R._alignment_code_key() == k    # and restoring it restores the key


def test_the_swath_cache_is_not_inside_package_data():
    """It is written at RUNTIME. It lived beside the bundled control CSVs until
    2026-09-10, because CACHE was module-relative and the module was promoted."""
    from lidar_diff_icp.groundtruth import reconstruct as R
    from lidar_diff_icp.groundtruth import residual_field as RF
    assert RF._DATA not in R.CACHE.parents, "runtime cache inside package data"
    assert R.CACHE.parent.name == "data" and R.CACHE.parent.parent.name != "groundtruth"
