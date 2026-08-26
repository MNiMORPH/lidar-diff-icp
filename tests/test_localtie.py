"""Ties measured at a stated place: the window, the ladder, the chain, and what breaks.

The synthetic tile writer is deliberately the one from ``test_groundtruth_chain.py``
(long-wavelength terrain with real slope and varied aspect, per-line Gaussian noise at
the measured 17 mm gen1 within-cell precision), extended with a **scan angle**, because
``tie="intercept"`` reads it and ``chain.py``'s pair cloud does not carry it.
"""
import numpy as np
import pytest

from lidar_diff_icp import localtie as LT


# --------------------------------------------------------------- synthetic tile writing

def _terrain(x0, x1, y0, y1, spacing=1.0):
    """Smooth surface with slope above Nuth & Kaeaeb's 3 deg floor and varied aspect."""
    gx = np.arange(x0, x1, spacing)
    gy = np.arange(y0, y1, spacing)
    X, Y = np.meshgrid(gx, gy)
    Z = 8.0 * np.sin(X / 60.0) + 8.0 * np.cos(Y / 70.0) + 0.10 * X - 0.06 * Y
    return X.ravel(), Y.ravel(), Z.ravel()


def _flat(x0, x1, y0, y1, spacing=1.0):
    gx = np.arange(x0, x1, spacing)
    gy = np.arange(y0, y1, spacing)
    X, Y = np.meshgrid(gx, gy)
    return X.ravel(), Y.ravel(), np.full(X.size, 300.0)


def _write_tile(path, lines, *, spacing=1.0, noise=0.017, surface=_terrain,
                nadir=None, angle_per_m=0.2, tan_error_m=0.0, dz_field=None):
    """``lines`` maps psid -> (x0, x1, y0, y1, dz_m).

    ``nadir`` maps psid -> ``(easting of that line's nadir track, flight direction)``,
    from which the integer ``scan_angle_rank`` is written as
    ``direction * angle_per_m * (x - nadir)`` degrees. The direction matters and is not
    decoration: adjacent lines are flown there-and-back, so at the middle of their
    sidelap both see the ground at the SAME aircraft-relative angle and
    ``dtan = tan(scan_ref) - tan(scan_src)`` passes through zero -- which is the position
    ``coreg.across_track_tie`` reports its intercept at. Give both lines the same
    direction and dtan never reaches zero, and the intercept becomes an extrapolation.
    ``tan_error_m`` adds ``tan_error_m * tan(scan_angle)`` to EVERY line's z -- a common
    across-track error, so the between-line difference is ``dz + tan_error_m * dtan`` and
    the intercept tie and the overlap-median tie have different, known, answers.
    ``dz_field`` (callable of x, y) adds a spatially varying offset to the last line, so
    the tie has a different true value in different windows.
    """
    import laspy
    xs, ys, zs, ps, ang = [], [], [], [], []
    last = list(lines)[-1]
    for sid, (x0, x1, y0, y1, dz) in lines.items():
        x, y, z = surface(x0, x1, y0, y1, spacing)
        a = (np.zeros(x.size) if nadir is None
             else nadir[sid][1] * angle_per_m * (x - nadir[sid][0]))
        z = z + dz + np.random.default_rng(sid).normal(0.0, noise, z.size)
        z = z + tan_error_m * np.tan(np.radians(np.round(a)))
        if dz_field is not None and sid == last:
            z = z + dz_field(x, y)
        xs.append(x); ys.append(y); zs.append(z)
        ps.append(np.full(x.size, sid, np.uint16))
        ang.append(np.round(a))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = np.array([x.min(), y.min(), z.min()])
    hdr.scales = np.array([0.001, 0.001, 0.001])
    las = laspy.LasData(hdr)
    las.x = x; las.y = y; las.z = z
    las.point_source_id = np.concatenate(ps)
    las.classification = np.full(x.size, 2, np.uint8)
    las.scan_angle_rank = np.concatenate(ang).astype(np.int8)
    las.write(str(path))
    return str(path)


@pytest.fixture
def pair_tile(tmp_path):
    """Lines 10 and 11 over x 0-100 and 60-160; line 11 sits 40 mm HIGH."""
    return _write_tile(tmp_path / "pair.las",
                       {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                        11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                       nadir={10: (50.0, +1), 11: (110.0, -1)})


EX = (5, 6, 9)          # the repo's terrain proxy, named at every call site on purpose


# ------------------------------------------------------------------------ the window

def test_the_window_is_the_size_asked_for_and_the_shape_asked_for(pair_tile):
    sq = LT.window_cloud([pair_tile], easting=80.0, northing=100.0,
                         half_width_m=20.0, shape="square")
    assert np.abs(sq.pc.x - 80.0).max() <= 20.0
    assert np.abs(sq.pc.y - 100.0).max() <= 20.0
    dk = LT.window_cloud([pair_tile], easting=80.0, northing=100.0,
                         half_width_m=20.0, shape="disk")
    assert np.hypot(dk.pc.x - 80.0, dk.pc.y - 100.0).max() <= 20.0
    assert len(dk.pc) < len(sq.pc)                      # a disk is smaller than its box
    with pytest.raises(ValueError):
        LT.window_cloud([pair_tile], easting=80.0, northing=100.0,
                        half_width_m=20.0, shape="hexagon")


def test_the_window_carries_the_scan_angle_in_degrees(pair_tile):
    """chain.py zeroes this field, which silently turns tie='intercept' into the median
    tie. The whole point of this module's own reader is that it does not."""
    w = LT.window_cloud([pair_tile], easting=80.0, northing=100.0,
                        half_width_m=20.0, shape="square")
    a10 = w.pc.scan_angle[w.pc.point_source_id == 10]
    a11 = w.pc.scan_angle[w.pc.point_source_id == 11]
    # both lines see this ground on the same side of their own track, because they are
    # flown in opposite directions -- so both angles are positive and dtan crosses zero
    assert a10.min() > 0 and a11.min() > 0
    assert np.ptp(a10) > 1.0 and np.ptp(a11) > 1.0      # degrees, not raw counts
    # and at the sidelap centre they are equal, which is where the intercept is read
    mid = LT.window_cloud([pair_tile], easting=80.0, northing=100.0,
                          half_width_m=2.0, shape="square")
    assert np.median(mid.pc.scan_angle[mid.pc.point_source_id == 10]) == \
        pytest.approx(np.median(mid.pc.scan_angle[mid.pc.point_source_id == 11]))


def test_there_is_no_default_window_size_or_tie_mode():
    """Every cut is a caller-supplied argument. Omitting one must fail loudly."""
    with pytest.raises(TypeError):
        LT.window_cloud(["x.las"], easting=0.0, northing=0.0)     # no half_width/shape
    with pytest.raises(TypeError):
        LT.pair_tie_at(["x.las"], 10, 11, easting=0.0, northing=0.0,
                       half_width_m=50.0, shape="square")          # no res/tie/exclude


# --------------------------------------------------------------------- the pair tie

def test_a_known_offset_is_recovered_in_both_tie_modes(pair_tile):
    for tie in ("overlap_median", "intercept"):
        t = LT.pair_tie_at([pair_tile], 10, 11, easting=80.0, northing=100.0,
                           half_width_m=95.0, shape="square", res_m=2.0, tie=tie,
                           exclude=EX)
        assert t.tie_mode == tie
        # dz is added to SRC to reach REF: line 11 sits 40 mm high, so it is pushed down
        assert t.dz_m == pytest.approx(-0.040, abs=0.006), tie
        assert not t.degenerate


def test_the_two_tie_modes_split_a_common_across_track_error(tmp_path):
    """With a common error ``a*tan(scan)`` on both lines the between-line difference is
    ``dz + a*dtan``: the intercept tie reads ``dz`` at dtan = 0, the overlap median reads
    it at the sampled mean dtan. They must therefore DIFFER, in the known direction."""
    p = _write_tile(tmp_path / "tan.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (130.0, -1)}, tan_error_m=0.5)
    kw = dict(easting=80.0, northing=100.0, half_width_m=95.0, shape="square",
              res_m=2.0, exclude=EX)
    med = LT.pair_tie_at([p], 10, 11, tie="overlap_median", **kw)
    itc = LT.pair_tie_at([p], 10, 11, tie="intercept", **kw)
    # line 11's nadir is at x = 130 while the pair only overlaps over x 60-100, so the
    # sampled dtan straddles zero ASYMMETRICALLY: the median tie is read at a mean dtan
    # that is not zero and the intercept tie is read at dtan = 0.
    assert abs(med.dz_m - itc.dz_m) > 0.010
    assert itc.dz_m == pytest.approx(-0.040, abs=0.012)


def test_the_across_track_diagnostic_reproduces_coreg_own_intercept(tmp_path):
    """``k_check_m`` is computed by this module's own mirror of ``coregister_swaths``'s
    intercept branch. If the mirror drifts from the original, the two stop agreeing --
    so pin them together here rather than assuming the mirroring holds."""
    p = _write_tile(tmp_path / "chk.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (130.0, -1)}, tan_error_m=0.5)
    kw = dict(easting=80.0, northing=100.0, half_width_m=95.0, shape="square",
              res_m=2.0, exclude=EX)
    itc = LT.pair_tie_at([p], 10, 11, tie="intercept", **kw)
    assert itc.k_check_m == pytest.approx(itc.dz_m, abs=1e-12)
    med = LT.pair_tie_at([p], 10, 11, tie="overlap_median", **kw)
    assert med.k_check_m != pytest.approx(med.dz_m)        # a different estimator
    assert itc.c_mm_per_tan == pytest.approx(500.0, rel=0.25)   # the 0.5 m/tan we wrote
    assert not itc.extrapolated                            # dtan straddles zero here


def test_a_window_that_does_not_reach_dtan_zero_is_flagged_as_extrapolated(tmp_path):
    """The intercept tie is read at dtan = 0. A window cut near one edge of the sidelap
    may not sample dtan = 0 at all, and then the 'extent-invariant' tie is an
    extrapolation. It is reported, not corrected and not dropped."""
    p = _write_tile(tmp_path / "edge.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)})
    edge = LT.pair_tie_at([p], 10, 11, easting=98.0, northing=100.0, half_width_m=12.0,
                          shape="square", res_m=2.0, tie="intercept", exclude=EX)
    assert edge.dtan_min > 0.0 and edge.extrapolated is True
    mid = LT.pair_tie_at([p], 10, 11, easting=80.0, northing=100.0, half_width_m=12.0,
                         shape="square", res_m=2.0, tie="intercept", exclude=EX)
    assert mid.dtan_min < 0.0 < mid.dtan_max and mid.extrapolated is False


def test_the_tie_is_local__two_windows_read_their_own_offset(tmp_path):
    """The module's whole purpose. Line 11 carries an offset that varies with northing;
    a tie measured in the south window must read the SOUTH value and one in the north
    window the NORTH value. A tie fitted over the whole extent reads neither."""
    p = _write_tile(tmp_path / "vary.las",
                    {10: (0.0, 100.0, 0.0, 600.0, 0.0),
                     11: (60.0, 160.0, 0.0, 600.0, 0.0)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)},
                    dz_field=lambda x, y: 0.100 * (y / 600.0))     # 0 -> 100 mm south to north
    kw = dict(half_width_m=80.0, shape="square", res_m=2.0, tie="overlap_median",
              exclude=EX)
    south = LT.pair_tie_at([p], 10, 11, easting=80.0, northing=90.0, **kw)
    north = LT.pair_tie_at([p], 10, 11, easting=80.0, northing=510.0, **kw)
    assert south.dz_m == pytest.approx(-0.015, abs=0.008)
    assert north.dz_m == pytest.approx(-0.085, abs=0.008)
    assert south.dz_m - north.dz_m == pytest.approx(0.070, abs=0.010)


def test_a_flat_window_is_flagged_degenerate_and_the_median_read_survives(tmp_path):
    """REGRESSION. ``coreg.nuth_kaab`` abandons its fit when fewer than 100 cells clear
    its 3 deg slope floor, and in that branch it returns ``dz = 0.0`` EXACTLY with
    ``n = 0`` -- not NaN, not the overlap median. On a whole tile that never happens; on
    a small window over flat ground it does, and a silent 0.0 mm tie would be the worst
    possible failure for a module whose job is small windows. The flag and the
    independent overlap-median column are what make it visible."""
    p = _write_tile(tmp_path / "flat.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)}, surface=_flat)
    t = LT.pair_tie_at([p], 10, 11, easting=80.0, northing=100.0, half_width_m=95.0,
                       shape="square", res_m=2.0, tie="overlap_median", exclude=EX)
    assert t.n_nk_cells == 0
    assert t.dz_m == 0.0                       # coreg's abandoned-fit value, not the truth
    assert t.degenerate is True                # ... and the caller is told
    assert t.dz_overlap_median_m == pytest.approx(-0.040, abs=0.005)  # the truth is beside it
    assert t.row()[-1] == "YES"


def test_lines_that_do_not_overlap_in_the_window_raise(pair_tile):
    with pytest.raises(ValueError):
        LT.pair_tie_at([pair_tile], 10, 11, easting=20.0, northing=100.0,
                       half_width_m=15.0, shape="square", res_m=2.0,
                       tie="overlap_median", exclude=EX)


# ------------------------------------------------------------------- the window ladder

def test_the_ladder_reports_every_requested_window_and_its_spread(pair_tile):
    lad = LT.window_ladder([pair_tile], 10, 11, easting=80.0, northing=100.0,
                           half_widths_m=[30.0, 60.0, 95.0], shape="square", res_m=2.0,
                           tie="overlap_median", exclude=EX)
    assert [t.half_width_m for t in lad.ties] == [30.0, 60.0, 95.0]
    assert [t.n_overlap_cells for t in lad.ties] == sorted(
        [t.n_overlap_cells for t in lad.ties])       # bigger window, more overlap cells
    assert np.isfinite(lad.spread_mm) and lad.spread_mm >= 0.0


def test_a_window_that_cannot_be_tied_is_reported_not_silently_dropped(pair_tile):
    """Nothing in this module drops a measurement: 'too small to tie here' IS the answer
    to the question the ladder asks."""
    lad = LT.window_ladder([pair_tile], 10, 11, easting=20.0, northing=100.0,
                           half_widths_m=[15.0, 95.0], shape="square", res_m=2.0,
                           tie="overlap_median", exclude=EX)
    assert len(lad.ties) == 2
    assert not np.isfinite(lad.ties[0].dz_m) and lad.ties[0].degenerate
    assert np.isfinite(lad.ties[1].dz_m)


def test_the_ladder_grows_when_the_offset_varies_across_the_window(tmp_path):
    """A tie that is a constant of the flight line does not care how big the window is.
    One that is local does. That difference is the ladder's whole content."""
    p = _write_tile(tmp_path / "vary2.las",
                    {10: (0.0, 100.0, 0.0, 600.0, 0.0),
                     11: (60.0, 160.0, 0.0, 600.0, 0.0)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)},
                    dz_field=lambda x, y: 0.200 * (y / 600.0))
    q = _write_tile(tmp_path / "const2.las",
                    {10: (0.0, 100.0, 0.0, 600.0, 0.0),
                     11: (60.0, 160.0, 0.0, 600.0, 0.100)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)})
    kw = dict(easting=80.0, northing=110.0, half_widths_m=[40.0, 110.0],
              shape="square", res_m=2.0, tie="overlap_median", exclude=EX)
    assert LT.window_ladder([p], 10, 11, **kw).spread_mm > \
           LT.window_ladder([q], 10, 11, **kw).spread_mm


# ----------------------------------------------------------------------- the network

def test_the_local_network_solves_every_line_into_the_gauge_line_frame(tmp_path):
    p = _write_tile(tmp_path / "three.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.000),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040),
                     12: (120.0, 220.0, 0.0, 200.0, 0.100)},
                    nadir={10: (50.0, +1), 11: (110.0, -1), 12: (170.0, +1)})
    net = LT.local_network([p], [10, 11, 12], easting=110.0, northing=100.0,
                           half_width_m=110.0, shape="square", res_m=2.0,
                           tie="overlap_median", exclude=EX, ref_line=10)
    c = net.constants_mm()
    assert c[10] == pytest.approx(0.0, abs=1e-9)            # the gauge
    # tolerances are 15 mm because the terrain here tilts 0.10 m/m eastward, so the
    # sub-decimetre horizontal shift Nuth & Kaeaeb solves trades directly against dz.
    assert c[11] == pytest.approx(-40.0, abs=15.0)          # 11 sits high -> pushed down
    assert c[12] == pytest.approx(-100.0, abs=20.0)


# ---------------------------------------------------------------- overlap point + chain

def test_the_nearest_overlap_point_is_in_the_overlap_and_is_the_nearest(pair_tile):
    op = LT.nearest_overlap_point([pair_tile], 10, 11, easting=0.0, northing=100.0,
                                  res_m=2.0, exclude=EX)
    assert 60.0 <= op.easting <= 100.0                     # inside the pair's overlap
    assert op.distance_m == pytest.approx(np.hypot(op.easting - 0.0,
                                                   op.northing - 100.0), rel=1e-9)
    assert op.easting == pytest.approx(61.0, abs=2.0)      # the closest edge of it
    assert op.n_overlap_cells > 0


def test_a_two_link_chain_accumulates_and_reports_how_far_it_reached(tmp_path):
    p = _write_tile(tmp_path / "chain.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.000),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040),
                     12: (120.0, 220.0, 0.0, 200.0, 0.100)},
                    nadir={10: (50.0, +1), 11: (110.0, -1), 12: (170.0, +1)})
    ch = LT.chain_local([p], easting=10.0, northing=100.0, source_line=10,
                        target_line=12, half_width_m=40.0, shape="square", res_m=2.0,
                        tie="overlap_median", exclude=EX)
    assert ch.nodes == [10, 11, 12]
    assert len(ch.links) == 2
    # line 10 sits 100 mm BELOW line 12, so +100 mm puts it in 12's frame
    assert ch.dz_total_mm == pytest.approx(100.0, abs=15.0)
    assert ch.dz_sigma_formal_m > 0
    # the far link cannot be solved at the mark; the module says how far away it was
    assert ch.max_solve_distance_m > 100.0
    # the walk runs from the TARGET inward, so the first link is the farthest from the mark
    assert [l.solve_distance_m for l in ch.links] == sorted(
        [l.solve_distance_m for l in ch.links], reverse=True)
    assert ch.degenerate_links == []


def test_the_chain_direction_flips_the_sign(tmp_path):
    p = _write_tile(tmp_path / "flip.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.000),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)})
    kw = dict(easting=80.0, northing=100.0, half_width_m=40.0, shape="square",
              res_m=2.0, tie="overlap_median", exclude=EX)
    a = LT.chain_local([p], source_line=10, target_line=11, **kw)
    b = LT.chain_local([p], source_line=11, target_line=10, **kw)
    assert a.dz_total_mm == pytest.approx(-b.dz_total_mm, abs=3.0)


def test_a_requested_path_that_does_not_run_source_to_target_raises(pair_tile):
    with pytest.raises(ValueError):
        LT.chain_local([pair_tile], easting=80.0, northing=100.0, source_line=10,
                       target_line=11, half_width_m=40.0, shape="square", res_m=2.0,
                       tie="overlap_median", exclude=EX, path=[11, 10])


def test_the_chain_ladder_gives_an_empirical_sigma_beside_the_formal_one(pair_tile):
    ch = LT.chain_local([pair_tile], easting=80.0, northing=100.0, source_line=10,
                        target_line=11, half_width_m=40.0, shape="square", res_m=2.0,
                        tie="overlap_median", exclude=EX,
                        ladder_half_widths_m=[25.0, 40.0, 60.0])
    assert ch.links[0].ladder is not None
    assert np.isfinite(ch.dz_sigma_window_m)
    assert ch.rows()[0][-1] != ""


def test_the_route_is_planned_by_groundtruth_chain_not_reimplemented(tmp_path, monkeypatch):
    from lidar_diff_icp.groundtruth import chain as gt_chain
    p = _write_tile(tmp_path / "route.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.000),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040),
                     12: (120.0, 220.0, 0.0, 200.0, 0.100)},
                    nadir={10: (50.0, +1), 11: (110.0, -1), 12: (170.0, +1)})
    seen = {}
    real = gt_chain.plan_path
    monkeypatch.setattr(gt_chain, "plan_path",
                        lambda *a, **k: seen.setdefault("called", real(*a, **k)))
    nodes = LT.plan_path_local([p], 10, 12, res_m=2.0, exclude=EX, cache=LT.TileCache())
    assert nodes == [10, 11, 12]
    assert "called" in seen


def test_no_route_is_a_reportable_error_not_a_silent_zero(tmp_path):
    p = _write_tile(tmp_path / "gap.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.0),
                     12: (300.0, 400.0, 0.0, 200.0, 0.1)},
                    nadir={10: (50.0, +1), 12: (350.0, -1)})
    with pytest.raises(ValueError, match="no overlap route"):
        LT.plan_path_local([p], 10, 12, res_m=2.0, exclude=EX, cache=LT.TileCache())


# ------------------------------------------------------- against an imported constant

def test_comparison_against_imported_constants_cancels_their_gauge(tmp_path):
    p = _write_tile(tmp_path / "cmp.las",
                    {10: (0.0, 100.0, 0.0, 200.0, 0.000),
                     11: (60.0, 160.0, 0.0, 200.0, 0.040)},
                    nadir={10: (50.0, +1), 11: (110.0, -1)})
    ch = LT.chain_local([p], easting=80.0, northing=100.0, source_line=10,
                        target_line=11, half_width_m=40.0, shape="square", res_m=2.0,
                        tie="overlap_median", exclude=EX)
    gauge_a = LT.compare_to_constants(ch, {10: 0.000, 11: -0.040})
    gauge_b = LT.compare_to_constants(ch, {10: 0.040, 11: 0.000})   # same set, regauged
    assert gauge_a.imported_mm == pytest.approx(gauge_b.imported_mm)
    assert gauge_a.difference_mm == pytest.approx(ch.dz_total_mm - 40.0)


# ------------------------------------------------------------------------ the cache

def test_the_disk_cache_round_trips_every_field(tmp_path, pair_tile):
    c1 = LT.TileCache(cache_dir=str(tmp_path / "cache"))
    a = c1.tile(pair_tile)
    c2 = LT.TileCache(cache_dir=str(tmp_path / "cache"))
    b = c2.tile(pair_tile)                      # from the npz this time
    for k in ("x", "y", "z"):
        assert np.abs(a[k] - b[k]).max() < 1e-3       # float32 offsets: ~1 mm over a tile
    assert np.array_equal(a["psid"], b["psid"])
    assert np.array_equal(a["cls"], b["cls"])
    assert np.array_equal(a["ang"], b["ang"])
    c1.release()
    assert c1._arrays == {}
