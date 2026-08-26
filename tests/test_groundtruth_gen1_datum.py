"""gen1's datum from its own 2008 control -- against data whose answer is known.

Three of these are regression tests in the strict sense: each was shown to FAIL with the
code under test removed, and the failure recorded in ``analysis/GEN1_DATUM_MODULE.md`` §7.

* ``test_line_assignment_comes_from_the_returns_not_the_centreline`` -- build a mark that
  sits closest to line A's fitted centreline while every ground return under it belongs
  to line B. Assigning by centreline gives A; the module must give B.
* ``test_the_clustered_se_does_not_shrink_when_marks_are_added_to_one_line`` -- adding
  marks to a line that is already represented must not shrink the datum's SE, because
  they share that swath's constant. The per-mark SE does shrink, and is returned beside
  it so the difference is visible.
* ``test_duplicate_rows_are_merged_so_one_mark_is_weighted_once`` -- the bundled CSV
  publishes 41 rows twice (marks on county lines appear in both counties' reports).
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


def test_the_assertion_runs_inside_a_measurement():
    x, y = _disc()
    g = _returns(x, y, 100.0, 7)
    with pytest.raises(G.DatumMismatchError):
        G.measure_site(_site("A", 0, 0, 100.0, geoid="GEOID18"), "fake.laz",
                       ground_loader=_loader(g))


def test_the_bundled_control_is_all_on_gen1s_own_geoid():
    cs = G.load_control()
    assert G.assert_no_geoid_conversion(list(cs)).startswith("no geoid conversion")


# ------------------------------------------------------------------------- discovery

def test_the_search_radius_has_no_default():
    cs = G.load_control()
    with pytest.raises(TypeError):
        G.discover_near_point(cs, 579705.72, 4883677.71)


def test_marks_come_back_nearest_first_and_inside_the_radius():
    cs = G.load_control()
    e, n = 579705.72, 4883677.71                    # the Elba reference point
    got = G.discover_near_point(cs, e, n, 10_000.0)
    d = [s.distance_m for s in got]
    assert d == sorted(d)
    assert max(d) <= 10_000.0
    assert all(np.hypot(s.mark.easting - e, s.mark.northing - n) == pytest.approx(s.distance_m)
               for s in got)
    assert len(G.discover_near_point(cs, e, n, 5_000.0)) < len(got)


def test_a_line_search_measures_perpendicular_distance_to_the_track():
    cs = G.ControlSet([_mark("ON", 100.0, 500.0, 1.0), _mark("OFF", 5000.0, 500.0, 1.0)],
                      origin="test", n_rows=2)
    lines = {137: ((0.0, 0.0), (0.0, 10_000.0))}    # a due-north track on easting 0
    got = G.discover_near_lines(cs, lines, 730.0)
    assert [s.point_id for s in got] == ["ON"]
    assert got[0].distance_m == pytest.approx(100.0)
    assert got[0].nearest_feature == "137"


def test_a_line_search_uses_the_nearest_of_several_tracks():
    cs = G.ControlSet([_mark("M", 900.0, 500.0, 1.0)], origin="test", n_rows=1)
    lines = {136: ((0.0, 0.0), (0.0, 10_000.0)), 137: ((1000.0, 0.0), (1000.0, 10_000.0))}
    got = G.discover_near_lines(cs, lines, 730.0)
    assert got[0].nearest_feature == "137"
    assert got[0].distance_m == pytest.approx(100.0)


# ------------------------------------------------------------------ tile resolution

def test_tiles_are_reported_as_on_disk_or_to_fetch_and_never_fetched(tmp_path, monkeypatch):
    calls = []

    def fake_find_tile(e, n, **kw):
        calls.append((e, n))
        return "4342-29-64" if e < 100 else "9999-99-99"

    import lidar_diff_icp.tiles as tiles
    monkeypatch.setattr(tiles, "find_tile", fake_find_tile)
    monkeypatch.setattr(tiles, "download_tile",
                        lambda *a, **k: pytest.fail("resolve_tiles must never download"))
    (tmp_path / "4342-29-64.laz").write_bytes(b"x")

    sites = [_site("HERE", 10.0, 0.0, 1.0), _site("AWAY", 1000.0, 0.0, 1.0)]
    res = G.resolve_tiles(sites, [tmp_path])
    assert len(calls) == 2
    assert [n.tile for n in res.on_disk] == ["4342-29-64"]
    assert [n.tile for n in res.to_fetch] == ["9999-99-99"]
    assert res.path_for("HERE").endswith("4342-29-64.laz")
    assert res.path_for("AWAY") is None
    assert res.per_mark == {"HERE": "4342-29-64", "AWAY": "9999-99-99"}
    assert set(res.table_columns()) >= {"tile", "on_disk", "n_marks"}


def test_a_missing_tile_is_skipped_with_its_reason_or_raises(tmp_path, monkeypatch):
    import lidar_diff_icp.tiles as tiles
    monkeypatch.setattr(tiles, "find_tile", lambda e, n, **kw: "9999-99-99")
    res = G.resolve_tiles([_site("AWAY", 1.0, 0.0, 1.0)], [tmp_path])
    meas, skipped = G.measure_sites([_site("AWAY", 1.0, 0.0, 1.0)], res)
    assert meas == []
    assert skipped == [("AWAY", "9999-99-99", "tile not on disk")]
    with pytest.raises(FileNotFoundError, match="never downloads"):
        G.measure_sites([_site("AWAY", 1.0, 0.0, 1.0)], res, on_missing="raise")


# -------------------------------------------------- flight lines FROM THE RETURNS

def test_line_assignment_comes_from_the_returns_not_the_centreline():
    """REGRESSION, proven to fail without the returns-based assignment.

    The mark sits 100 m from line 136's fitted track and 900 m from line 137's, so a
    centreline assignment calls it 136. Every ground return under it carries
    ``point_source_id`` 137 -- the acquisition's own record of which sortie lit the
    ground. The module must say 137.
    """
    e, n = 100.0, 500.0
    tracks = {136: ((0.0, 0.0), (0.0, 10_000.0)), 137: ((1000.0, 0.0), (1000.0, 10_000.0))}
    nearest_track = min(tracks, key=lambda k: G._point_segment_distance(
        e, n, *tracks[k][0], *tracks[k][1]))
    assert nearest_track == 136                       # what the centreline rule would say

    x, y = _disc(cx=e, cy=n)
    got = G.assign_line_from_returns(_returns(x, y, 100.0, 137), e, n, 7.5)
    assert got.dominant == 137
    assert got.n_lines == 1
    assert not got.mixed
    assert got.dominant_fraction == pytest.approx(1.0)


def test_a_mark_lit_by_two_lines_is_flagged_and_counted_not_dropped():
    x, y = _disc(n=2000)
    lines = np.where(x < 0, 136, 137)
    a = G.assign_line_from_returns(_returns(x, y, 100.0, lines), 0.0, 0.0, 7.5)
    assert a.mixed and a.n_lines == 2
    assert set(a.counts) == {136, 137}
    assert a.n == sum(a.counts.values())
    assert 0.0 < a.dominant_fraction < 1.0

    m = G.measure_site(_site("MIX", 0.0, 0.0, 100.0), "fake.laz",
                       ground_loader=_loader(_returns(x, y, 100.0, lines)))
    assert m.line.mixed
    assert any("flight lines at the mark" in s for s in m.notes)
    assert set(m.per_line_tie_mm) == {136, 137}


def test_an_empty_window_gives_no_line_rather_than_a_guess():
    a = G.assign_line_from_returns(_returns([1000.0], [1000.0], [1.0], 5), 0.0, 0.0, 7.5)
    assert a.dominant is None and a.counts == {} and a.n_lines == 0


def test_the_assignment_radius_is_the_report_radius_of_the_estimator():
    x, y = _disc(cx=0.0, cy=0.0, r=30.0)
    lines = np.where(np.hypot(x, y) < 7.5, 137, 999)   # 999 only outside the report radius
    m = G.measure_site(_site("R", 0.0, 0.0, 100.0), "fake.laz",
                       ground_loader=_loader(_returns(x, y, 100.0, lines)), res=5.0)
    assert m.line.radius_m == pytest.approx(7.5)
    assert m.line.dominant == 137 and not m.line.mixed


# ------------------------------------------------------------- measurement + screen

def test_a_known_offset_is_recovered_with_no_geoid_and_no_lateral_term():
    m = _flat_site_measure("A", 0.123, 137)
    assert m.tie_mm == pytest.approx(123.0, abs=1e-3)
    assert m.swath_shift_m == (0.0, 0.0, 0.0)
    p = {q.name: q for q in m.params}
    assert p["geoid_shift_m"].value == 0.0
    assert p["lateral_shift_m"].value is None
    assert "no place in gen1 vs its own control" in p["lateral_shift_m"].why


def test_the_screen_reports_the_statistics_and_applies_no_cut():
    x, y = _disc()
    z = 100.0 + 0.10 * x - 0.05 * y                 # a 10%/5% plane, no relief beyond it
    m = G.measure_site(_site("S", 0.0, 0.0, 100.0), "fake.laz",
                       ground_loader=_loader(_returns(x, y, z, 137)))
    s = m.screen
    assert s.n > 0
    assert s.slope_deg == pytest.approx(np.degrees(np.arctan(np.hypot(0.10, 0.05))), abs=1e-6)
    assert s.fit_rms_mm == pytest.approx(0.0, abs=1e-6)
    assert s.radius_spread_mm == pytest.approx(0.0, abs=1e-6)
    assert s.relief_mm > 0
    assert set(G.SitingScreen.table_columns()) == {
        "n", "slope_deg", "relief_mm", "fit_rms_mm", "radius_spread_mm"}


def test_the_crop_half_width_is_derived_and_cannot_truncate_the_ladder():
    m = _flat_site_measure("A", 0.010, 137)
    p = {q.name: q for q in m.params}
    assert p["crop_half_width_m"].value == pytest.approx(25.0)      # 5 * res
    assert "DERIVED" in p["crop_half_width_m"].why
    x, y = _disc()
    with pytest.raises(ValueError, match="smaller than the largest fitting radius"):
        G.measure_site(_site("A", 0.0, 0.0, 100.0), "fake.laz",
                       ground_loader=_loader(_returns(x, y, 100.0, 137)),
                       crop_half_width_m=10.0)


def test_a_bigger_crop_changes_nothing_for_a_vendor_class_read():
    """The derived half-width is the smallest square holding every fitting window, so a
    larger one cannot move the answer -- which is why it is derived and not chosen."""
    x, y = _disc(r=200.0, n=20000)
    g = _returns(x, y, 100.0 + 0.03 * x, 137)
    a = G.measure_site(_site("A", 0.0, 0.0, 100.0), "f.laz", ground_loader=_loader(g))
    b = G.measure_site(_site("A", 0.0, 0.0, 100.0), "f.laz", ground_loader=_loader(g),
                       crop_half_width_m=300.0)
    assert a.tie_mm == pytest.approx(b.tie_mm, abs=1e-9)


def test_the_swath_constant_moves_the_tie_by_its_dz_and_is_recorded():
    plain = _flat_site_measure("A", 0.0, 137)
    shifted = _flat_site_measure("A", 0.0, 137,
                                 swath_constants={137: (0.0, 0.0, -0.020)},
                                 swath_constants_source="test corrections.json")
    assert plain.tie_mm == pytest.approx(0.0, abs=1e-6)
    assert shifted.tie_mm == pytest.approx(20.0, abs=1e-3)
    assert shifted.swath_constant_source == "test corrections.json"


def test_a_line_with_no_constant_is_flagged_not_dropped():
    m = _flat_site_measure("A", 0.0, 999, swath_constants={137: (0.0, 0.0, -0.02)},
                           swath_constants_source="test")
    assert m.swath_constant_source == ""
    assert any("no swath constant for line 999" in s for s in m.notes)


def test_the_gen2_lateral_shift_is_opt_in_and_recorded_when_taken():
    x, y = _disc()
    g = _returns(x, y, 100.0 + 0.10 * x, 137)         # 10% east slope: a shift moves z
    off = G.measure_site(_site("A", 0.0, 0.0, 100.0), "f.laz", ground_loader=_loader(g))
    on = G.measure_site(_site("A", 0.0, 0.0, 100.0), "f.laz", ground_loader=_loader(g),
                        lateral_shift_m=(-0.75, -0.19))
    assert off.tie_mm != pytest.approx(on.tie_mm)
    # the returns move 0.75 m west, so the surface read at the mark climbs 0.10*0.75 m
    # and the tie (surveyed - lidar) falls by that much
    assert on.tie_mm - off.tie_mm == pytest.approx(-0.10 * 0.75 * 1000.0, rel=1e-6)
    p = {q.name: q for q in on.params}
    assert p["lateral_shift_m"].value == (-0.75, -0.19)
    assert p["lateral_shift_m"].src == "andy"
    assert any("cross-epoch term" in s for s in on.notes)


def test_swath_constants_are_read_from_a_corrections_json(tmp_path):
    p = tmp_path / "corrections.json"
    p.write_text(json.dumps({"per_swath_internal_alignment_dxdydz_m":
                             {"135": [0, 0, 0], "136": [0.32, -0.08, -0.018]},
                             "swath_tie": "intercept", "ground_source": "csf", "res_m": 5.0}))
    const, src = G.swath_constants_from_corrections(p)
    assert const == {135: (0.0, 0.0, 0.0), 136: (0.32, -0.08, -0.018)}
    assert "swath_tie='intercept'" in src and str(p) in src
    (tmp_path / "empty.json").write_text("{}")
    with pytest.raises(KeyError):
        G.swath_constants_from_corrections(tmp_path / "empty.json")


# ------------------------------------------------------------------- the combination

def _measurements(spec, **kw):
    """spec: [(id, tie_mm, line)] -> measurements whose ties are exactly those."""
    return [_flat_site_measure(pid, t / 1000.0, line, e=100.0 * k, **kw)
            for k, (pid, t, line) in enumerate(spec)]


def test_per_line_averages_within_line_then_over_lines():
    m = _measurements([("a", 0.0, 136), ("b", 100.0, 136), ("c", 200.0, 137)])
    est = G.combine_datum(m, mode="per_line")
    assert est.n_marks == 3 and est.n_lines == 2
    assert [g.line for g in est.groups] == [136, 137]
    assert est.groups[0].mean_mm == pytest.approx(50.0, abs=1e-3)
    assert est.value_mm == pytest.approx(125.0, abs=1e-3)            # (50 + 200) / 2
    assert est.mean_over_marks_mm == pytest.approx(100.0, abs=1e-3)  # the naive answer
    assert est.value_mm != pytest.approx(est.mean_over_marks_mm)


def test_the_se_says_what_it_is_the_se_of():
    m = _measurements([("a", 0.0, 136), ("b", 100.0, 137), ("c", 200.0, 138)])
    est = G.combine_datum(m, mode="per_line")
    assert est.se_of.startswith("SE of the mean over flight lines of the within-line mean tie")
    assert "unit of replication" in est.se_of or "unit of replication" in \
        {p.name: p.value for p in est.params}.get("unit_of_replication", "flight line")
    line_means = np.array([0.0, 100.0, 200.0])
    assert est.se_mm == pytest.approx(line_means.std(ddof=1) / np.sqrt(3), abs=1e-3)
    assert "flight line, not the\nmark" in est.se_of.replace("  ", " ") or \
        "flight line, not the" in est.se_of


def test_the_clustered_se_does_not_shrink_when_marks_are_added_to_one_line():
    """REGRESSION, proven to fail if the marks are pooled as if independent.

    Four marks on two lines; then eight marks on the same two lines, each duplicated.
    The extra marks carry NO new information about the datum -- they share their swath's
    constant -- so the SE over lines must not move. The per-mark SE does shrink, and the
    design effect records by how much.
    """
    few = _measurements([("a", 0.0, 136), ("b", 40.0, 136),
                         ("c", 200.0, 137), ("d", 240.0, 137)])
    many = _measurements([("a", 0.0, 136), ("b", 40.0, 136), ("e", 0.0, 136), ("f", 40.0, 136),
                          ("c", 200.0, 137), ("d", 240.0, 137), ("g", 200.0, 137),
                          ("h", 240.0, 137)])
    e1 = G.combine_datum(few, mode="per_line")
    e2 = G.combine_datum(many, mode="per_line")
    assert e1.value_mm == pytest.approx(e2.value_mm, abs=1e-6)
    assert e1.se_mm == pytest.approx(e2.se_mm, abs=1e-6)
    assert e2.se_over_marks_mm < e1.se_over_marks_mm
    assert e2.n_marks == 8 and e2.n_lines == 2


def test_the_anova_reports_that_the_scatter_is_organised_by_line():
    tight = _measurements([("a", 0.0, 136), ("b", 2.0, 136),
                           ("c", 300.0, 137), ("d", 302.0, 137)])
    est = G.combine_datum(tight, mode="per_line")
    assert est.anova_F > 100.0 and est.anova_p < 1e-3
    assert est.anova_df == (1, 2)
    assert est.icc > 0.9
    mixed = _measurements([("a", 0.0, 136), ("b", 300.0, 136),
                           ("c", 2.0, 137), ("d", 302.0, 137)])
    assert G.combine_datum(mixed, mode="per_line").anova_F < 1.0


def test_one_line_only_cannot_separate_the_datum_from_that_swaths_constant():
    est = G.combine_datum(_measurements([("a", 0.0, 136), ("b", 100.0, 136)]),
                          mode="per_line")
    assert est.n_lines == 1
    assert not np.isfinite(est.se_mm)
    assert any("not separable" in s for s in est.notes)


def test_common_datum_needs_the_constants_and_returns_the_line_residuals():
    const = {136: (0.0, 0.0, 0.0), 137: (0.0, 0.0, -0.100)}
    m = _measurements([("a", 0.0, 136), ("c", 0.0, 137)],
                      swath_constants=const, swath_constants_source="test corrections.json")
    est = G.combine_datum(m, mode="common_datum", swath_constants_source="test corrections.json")
    assert est.mode == "common_datum"
    assert est.groups[0].mean_mm == pytest.approx(0.0, abs=1e-3)
    assert est.groups[1].mean_mm == pytest.approx(100.0, abs=1e-3)
    assert est.value_mm == pytest.approx(50.0, abs=1e-3)
    assert est.line_residual_mm[136] == pytest.approx(-50.0, abs=1e-3)
    assert est.line_residual_mm[137] == pytest.approx(+50.0, abs=1e-3)
    assert est.swath_constant_source == "test corrections.json"


def test_the_mode_is_checked_against_how_the_marks_were_measured():
    plain = _measurements([("a", 0.0, 136), ("c", 0.0, 137)])
    shifted = _measurements([("a", 0.0, 136), ("c", 0.0, 137)],
                            swath_constants={136: (0, 0, 0), 137: (0, 0, -0.1)},
                            swath_constants_source="test")
    with pytest.raises(ValueError, match="mode='per_line'"):
        G.combine_datum(shifted, mode="per_line")
    with pytest.raises(ValueError, match="no measurement survived"):
        G.combine_datum(plain, mode="common_datum")     # none is in the common frame


def test_a_mark_the_common_frame_cannot_hold_is_excluded_with_its_reason():
    const = {136: (0.0, 0.0, 0.0)}
    m = _measurements([("a", 0.0, 136), ("z", 0.0, 999)],
                      swath_constants=const, swath_constants_source="test")
    est = G.combine_datum(m, mode="common_datum", swath_constants_source="test")
    assert est.n_marks == 1
    assert [p for p, _ in est.excluded] == ["z"]
    assert "no swath constant for line 999" in est.excluded[0][1]


def test_an_unknown_mode_is_refused():
    with pytest.raises(ValueError, match="mode must be one of"):
        G.combine_datum(_measurements([("a", 0.0, 136)]), mode="whatever")


def test_the_estimate_serialises_with_its_convention_and_its_se_of():
    est = G.combine_datum(_measurements([("a", 0.0, 136), ("b", 100.0, 137)]),
                          mode="per_line")
    d = est.to_dict()
    assert d["sign_convention"].startswith("constant to ADD to gen1")
    assert d["se_of"] == est.se_of
    assert d["mode"] == "per_line"
    assert [g["line"] for g in d["line_groups"]] == [136, 137]
    json.dumps(d)                                   # must be serialisable
    assert "SE is:" in est.summary() and "ANOVA" in est.summary()


def test_every_returned_number_carries_a_definition():
    cols = G.MarkMeasurement.table_columns()
    for key in ("tie_mm", "sigma_mm", "line", "n_lines", "slope_deg", "relief_mm",
                "fit_rms_mm", "radius_spread_mm"):
        assert key in cols and len(cols[key].split()) >= 3
    est = G.combine_datum(_measurements([("a", 0.0, 136), ("b", 100.0, 137)]),
                          mode="per_line")
    assert set(est.table_columns()) >= {"line", "n", "mean_mm", "resid_mm"}
    assert len(est.table_rows()) == 2
