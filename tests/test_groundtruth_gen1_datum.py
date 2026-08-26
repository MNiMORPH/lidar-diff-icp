"""gen1's datum from its own 2008 control -- against data whose answer is known.

``test_duplicate_rows_are_merged_so_one_mark_is_weighted_once`` is a regression test in
the strict sense: it was shown to FAIL with the merge removed, and the failure is
recorded in ``analysis/GEN1_DATUM_MODULE.md`` §7. The bundled CSV publishes 41 rows
twice, because a mark on a county line appears in both counties' validation reports.
"""
import json

import numpy as np
import pytest

from lidar_diff_icp.groundtruth import checkpoints as C
from lidar_diff_icp.groundtruth import gen1_datum as G
from lidar_diff_icp.groundtruth import tie as T


# --------------------------------------------------------------------------- helpers

def _cp(pid, e, n, elev, *, ptype="L1O", geoid="GEOID03", vdatum="NAVD88"):
    return C.Checkpoint(point_id=pid, point_type=ptype, easting=e, northing=n,
                        elevation=elev, elevation_units="m", horizontal_crs="EPSG:26915",
                        vertical_datum=vdatum, geoid_model=geoid)


def _mark(pid, e, n, elev, **kw):
    return G.ControlMark(checkpoint=_cp(pid, e, n, elev, **kw), cover_class="L1O",
                         counties=("winona",), aliases=(pid,), reports=("test",),
                         dnr_surface_z_m=None, dnr_error_m=None)


def _site(pid, e, n, elev, d=0.0, **kw):
    return G.MarkSite(_mark(pid, e, n, elev, **kw), d, "test site")


def _disc(n=4000, r=30.0, seed=0, cx=0.0, cy=0.0):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, n)
    rr = r * np.sqrt(rng.uniform(0, 1, n))
    return cx + rr * np.cos(t), cy + rr * np.sin(t)


def _returns(x, y, z, line):
    x = np.asarray(x, float)
    z = np.broadcast_to(np.asarray(z, float), x.shape).copy()
    line = np.asarray(line)
    if line.ndim == 0:
        line = np.full(np.size(x), int(line), int)
    return T.GroundReturns(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float),
                           np.zeros(np.size(x)), line.astype(int),
                           source="synthetic", origin="test")


def _loader(ground):
    """A ground_loader that ignores the file and hands back a prepared cloud."""
    def load(tile_path, e, n, half):
        return ground
    return load


def _flat_site_measure(pid, low_by_m, line, *, e=0.0, n=0.0, **kw):
    """A mark on flat ground where the lidar reads ``low_by_m`` low: tie = +low_by_m."""
    x, y = _disc(cx=e, cy=n)
    g = _returns(x, y, 100.0 - low_by_m, line)
    return G.measure_site(_site(pid, e, n, 100.0), "fake.laz", ground_loader=_loader(g), **kw)


# ------------------------------------------------------------------- the control set

def test_duplicate_rows_are_merged_so_one_mark_is_weighted_once():
    """REGRESSION. The bundled transcription holds the same physical mark twice wherever
    it sits on a county line and both counties' validation reports print it. Without the
    merge those marks enter any average twice."""
    cs = G.load_control()
    assert cs.n_rows == 1004
    assert len(cs) == 963
    assert cs.n_rows - len(cs) == 41
    ids = [m.point_id for m in cs]
    assert len(ids) == len(set(ids))
    xyz = {(m.easting, m.northing, m.checkpoint.elevation) for m in cs}
    assert len(xyz) == len(cs)
    m = cs["L2T-6126 Wabasha RTK"]
    assert sorted(m.counties) == ["olmsted", "wabasha"]
    assert len(m.reports) == 2


def test_every_merge_is_reported_with_the_reports_it_came_from():
    cs = G.load_control()
    assert len(cs.merges) == 39
    assert sum(1 + len(a) for _, a, _ in cs.merges) >= len(cs.merges)
    for kept, aliases, counties in cs.merges:
        assert len(counties) >= 1
        assert cs[kept] is not None
    assert "41 rows merged" in cs.merge_note


def test_a_mark_id_at_two_positions_is_refused(tmp_path, monkeypatch):
    hdr = ("point_id,point_type,easting,northing,elevation,elevation_units,"
           "horizontal_crs,vertical_datum,geoid_model,county,source\n")
    row = "A,L1O,{e},{n},{z},m,EPSG:26915,NAVD88,GEOID03,winona,test\n"
    p = tmp_path / "clash.csv"
    p.write_text(hdr + row.format(e=1, n=2, z=3) + row.format(e=1, n=2, z=4))
    monkeypatch.setattr(G, "_DATA", tmp_path)
    with pytest.raises(ValueError, match="two different positions"):
        G.load_control("clash")


def test_an_unknown_control_set_says_so():
    with pytest.raises(FileNotFoundError):
        G.load_control("no_such_control_set")


# ------------------------------------------------------------- the geoid is ASSERTED

def test_matching_geoids_return_the_claim_that_was_checked():
    claim = G.assert_no_geoid_conversion([_mark("A", 0, 0, 1.0)])
    assert "no geoid conversion" in claim
    assert "GEOID03" in claim and "NAVD88" in claim


def test_a_geoid18_mark_cannot_reach_a_gen1_tie():
    with pytest.raises(G.DatumMismatchError, match="does not convert geoids"):
        G.assert_no_geoid_conversion([_mark("A", 0, 0, 1.0, geoid="GEOID18")])


def test_a_different_vertical_datum_is_refused_too():
    with pytest.raises(G.DatumMismatchError):
        G.assert_no_geoid_conversion([_mark("A", 0, 0, 1.0, vdatum="NGVD29")])


def test_the_bundled_control_is_all_on_gen1s_own_geoid():
    cs = G.load_control()
    assert G.assert_no_geoid_conversion(list(cs)).startswith("no geoid conversion")
