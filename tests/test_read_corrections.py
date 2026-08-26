"""``read_corrections`` must prefer the geoid sidecar when a tile carries both.

A tile rebuilt onto the geoid-only datum writes ``corrections_geoid.json`` and leaves the
older ``reference_plane`` ``corrections.json`` in place beside it (elbaext carries both on
disk today). Anything that opens ``corrections.json`` by name therefore reads the OBSOLETE
datum on exactly the tiles that have been brought up to date, and does so silently.
"""
import json

from lidar_diff_icp import registration as reg

GEOID = {"cross_epoch_datum": {"method": "geoid_difference", "const_m": 0.0667},
         "bounds": [0.0, 0.0, 100.0, 100.0], "res_m": 5.0}
PLANE = {"cross_epoch_datum": {"method": "reference_plane", "const_m": -0.0849},
         "bounds": [0.0, 0.0, 100.0, 100.0], "res_m": 5.0}


def _tile(tmp_path, **files):
    for name, body in files.items():
        (tmp_path / name).write_text(json.dumps(body))
    return str(tmp_path)


def test_geoid_sidecar_wins_when_both_are_present(tmp_path):
    d = _tile(tmp_path, **{"corrections.json": PLANE, "corrections_geoid.json": GEOID})
    assert reg.read_corrections(d)["cross_epoch_datum"]["method"] == "geoid_difference"
    # and the datum accessor must agree with the whole-file reader, or callers mixing the
    # two would silently combine terms from different datums
    assert reg.read_cross_epoch_datum(d)["method"] == "geoid_difference"


def test_falls_back_to_the_plain_name_when_that_is_all_there_is(tmp_path):
    d = _tile(tmp_path, **{"corrections.json": PLANE})
    assert reg.read_corrections(d)["cross_epoch_datum"]["method"] == "reference_plane"


def test_missing_sidecar_raises_rather_than_returning_empty(tmp_path):
    import pytest
    with pytest.raises(FileNotFoundError):
        reg.read_corrections(str(tmp_path))
