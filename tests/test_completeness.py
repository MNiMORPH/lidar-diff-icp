"""The completeness gate: unknown is not a pass, and no threshold is invented."""
import json

import pytest

from lidar_diff_icp import completeness as C


def test_an_unmeasured_site_raises_rather_than_passing(tmp_path):
    """THE POINT OF THE MODULE. A site whose completeness has never been measured must not
    be indistinguishable from one measured and found whole. Whitewater's truncated gen2
    averaged 11.39 returns/m2 over a tile that was 15.45 west and 5.52 east of a seam --
    nothing inside the tile showed it."""
    with pytest.raises(C.CompletenessUnknown, match="UNKNOWN IS NOT A PASS"):
        C.check(tmp_path)


def test_it_can_be_asked_to_proceed_while_saying_so(tmp_path):
    c = C.check(tmp_path, require=False)
    assert c["gen2"]["known"] is False
    assert "UNKNOWN" in C.summary_line(tmp_path, "gen2", c)


def test_the_ratio_is_stored_with_both_its_inputs(tmp_path):
    """A bare ratio cannot be rechecked, and which of the two numbers moved is the whole
    diagnosis -- a fetch that got shorter and a source that got bigger are different
    problems with the same ratio."""
    C.write(tmp_path, epoch="gen2", cloud="x.laz", points_in_file=148_050_625,
            points_available=140_000_000.0, measured_by="test")
    r = json.load(open(C.record_path(tmp_path)))["epochs"]["gen2"]
    assert r["points_in_file"] == 148_050_625
    assert r["points_available_in_bbox"] == 140_000_000.0
    assert r["ratio"] == pytest.approx(148_050_625 / 140_000_000.0)
    assert r["measured_by"] == "test" and r["cloud"] == "x.laz"


def test_a_ratio_against_zero_is_refused_not_recorded(tmp_path):
    """Dividing by an unmeasured availability would produce a number that looks like a
    measurement and is not one."""
    with pytest.raises(ValueError, match="not a measurement"):
        C.write(tmp_path, epoch="gen2", cloud="x.laz", points_in_file=10,
                points_available=0.0, measured_by="test")


def test_no_threshold_is_applied_anywhere(tmp_path):
    """DELIBERATE. A low ratio is REPORTED, never judged: what counts as complete enough for
    a particular question is a decision this module must not make silently. mnrv sits at
    0.67 and that is a fact to act on, not a failure to swallow."""
    C.write(tmp_path, epoch="gen2", cloud="x.laz", points_in_file=67,
            points_available=100.0, measured_by="test")
    c = C.check(tmp_path, epochs=("gen2",))
    assert c["gen2"]["known"] is True
    assert c["gen2"]["ratio"] == pytest.approx(0.67)
    # the source carries no comparison operator against a ratio
    import inspect
    src = inspect.getsource(C)
    assert "ratio >" not in src and "ratio <" not in src


def test_a_second_epoch_does_not_overwrite_the_first(tmp_path):
    C.write(tmp_path, epoch="gen2", cloud="a.laz", points_in_file=1, points_available=2.0,
            measured_by="t")
    C.write(tmp_path, epoch="gen1", cloud="b.laz", points_in_file=3, points_available=4.0,
            measured_by="t")
    rec = C.read(tmp_path)
    assert set(rec["epochs"]) == {"gen1", "gen2"}
    assert rec["epochs"]["gen2"]["cloud"] == "a.laz"
