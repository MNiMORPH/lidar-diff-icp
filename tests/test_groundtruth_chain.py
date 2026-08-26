"""The swath chain: along-swath first, then the shortest chain, with known offsets."""
import numpy as np
import pytest

from lidar_diff_icp.groundtruth import chain as K


# --------------------------------------------------------------- synthetic tile writing

def _terrain(x0, x1, y0, y1, spacing=1.0):
    """A smooth surface with real slope AND varied aspect over the overlap.

    Nuth & Kaeaeb solves the horizontal shift from dh/tan(slope) against aspect, so it
    needs slope above its 3 deg floor and a spread of aspects to condition the cosine
    fit; and the wavelengths are long (~370-440 m) so that a 40 m overlap carries no
    periodic ambiguity for it to lock onto."""
    gx = np.arange(x0, x1, spacing)
    gy = np.arange(y0, y1, spacing)
    X, Y = np.meshgrid(gx, gy)
    Z = 8.0 * np.sin(X / 60.0) + 8.0 * np.cos(Y / 70.0) + 0.10 * X - 0.06 * Y
    return X.ravel(), Y.ravel(), Z.ravel()


def _write_tile(path, lines, spacing=1.0, noise=0.017):
    """``lines`` maps point_source_id -> (x0, x1, y0, y1, dz). Writes a LAS whose points
    are class 2 (kept by the (5, 6, 9) terrain proxy).

    ``noise`` is per-return Gaussian scatter, seeded per line so the two lines are
    independent; 17 mm is the measured gen1 within-cell precision, and without it the
    link residual is identically zero and its uncertainty is meaningless."""
    import laspy
    xs, ys, zs, ps = [], [], [], []
    for sid, (x0, x1, y0, y1, dz) in lines.items():
        x, y, z = _terrain(x0, x1, y0, y1, spacing)
        z = z + dz + np.random.default_rng(sid).normal(0.0, noise, z.size)
        xs.append(x); ys.append(y); zs.append(z)
        ps.append(np.full(x.size, sid, np.uint16))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = np.array([x.min(), y.min(), z.min()])
    hdr.scales = np.array([0.001, 0.001, 0.001])
    las = laspy.LasData(hdr)
    las.x = x; las.y = y; las.z = z
    las.point_source_id = np.concatenate(ps)
    las.classification = np.full(x.size, 2, np.uint8)
    las.write(str(path))
    return path


@pytest.fixture
def two_line_tile(tmp_path):
    """Lines 10 and 11 overlapping over x 60-100, line 11 sitting 40 mm HIGH."""
    p = _write_tile(tmp_path / "pair.las",
                    {10: (0.0, 100.0, 0.0, 120.0, 0.0),
                     11: (60.0, 160.0, 0.0, 120.0, 0.040)})
    return p


# ------------------------------------------------------------------- inventory + graph

def test_inventory_finds_the_lines_and_caches(tmp_path, two_line_tile):
    cache = tmp_path / "cache"
    inv = K.build_inventory([two_line_tile], cache_dir=str(cache),
                            inventory_json=str(tmp_path / "inv.json"))
    assert inv.lines == [10, 11]
    assert inv.tiles_for(10) == [str(two_line_tile)]
    assert inv.tiles_for_pair(10, 11) == [str(two_line_tile)]
    # second build reads the JSON rather than the LAS
    inv2 = K.build_inventory([two_line_tile], cache_dir=str(cache),
                             inventory_json=str(tmp_path / "inv.json"))
    assert inv2.lines == [10, 11]


def test_overlap_area_is_measured_not_assumed(two_line_tile):
    inv = K.build_inventory([two_line_tile])
    g = K.overlap_graph(inv, res=2.0)
    assert list(g) == [(10, 11)]
    # the shared strip is x 60-100 by y 0-120 = 4800 m^2
    assert g[(10, 11)]["area_km2"] == pytest.approx(4800 / 1e6, rel=0.05)


def test_overlap_cells_are_global_so_a_pair_is_not_counted_twice(tmp_path):
    """The same overlap present in two tiles must be measured once."""
    a = _write_tile(tmp_path / "a.las", {10: (0.0, 100.0, 0.0, 60.0, 0.0),
                                         11: (60.0, 160.0, 0.0, 60.0, 0.0)})
    b = _write_tile(tmp_path / "b.las", {10: (0.0, 100.0, 30.0, 90.0, 0.0),
                                         11: (60.0, 160.0, 30.0, 90.0, 0.0)})
    inv = K.build_inventory([a, b])
    g = K.overlap_graph(inv, res=2.0)
    # union of y 0-60 and y 30-90 is y 0-90, NOT 120
    assert g[(10, 11)]["area_km2"] == pytest.approx(40 * 90 / 1e6, rel=0.05)


def test_covering_lines_reports_returns_near_a_point(two_line_tile):
    inv = K.build_inventory([two_line_tile])
    only10 = K.covering_lines(inv, 20.0, 60.0, 10.0)
    assert list(only10) == [10]
    both = K.covering_lines(inv, 80.0, 60.0, 10.0)
    assert set(both) == {10, 11}
    assert K.covering_lines(inv, 5000.0, 60.0, 10.0) == {}


# ------------------------------------------------------------------------ link solving

def test_a_known_link_offset_is_recovered(two_line_tile):
    inv = K.build_inventory([two_line_tile])
    g = K.overlap_graph(inv, res=2.0)
    L = K.solve_link(inv, 10, 11, res=2.0, graph=g)
    # line 11 was written 40 mm high, so aligning it onto 10 must subtract 40 mm
    assert L.dz_m == pytest.approx(-0.040, abs=0.003)
    assert L.overlap_km2 == pytest.approx(4800 / 1e6, rel=0.05)
    assert L.nmad_after_m <= L.nmad_before_m
    assert L.n_cells > 0
    assert L.tiles == [str(two_line_tile)]


def test_a_link_reads_only_the_overlap(tmp_path, monkeypatch):
    """Only the pair's shared strip is gridded, not the tiles that contain it."""
    p = _write_tile(tmp_path / "wide.las", {10: (0.0, 400.0, 0.0, 120.0, 0.0),
                                            11: (360.0, 760.0, 0.0, 120.0, 0.0)})
    inv = K.build_inventory([p])
    pc, tiles, box = K._pair_cloud(inv, 10, 11, res=2.0, pad=4.0)
    assert box[0] >= 360.0 - 4.0 - 1e-6 and box[2] <= 400.0 + 4.0 + 1e-6
    assert len(pc) < 0.2 * inv.tiles[str(p)].n_terrain


# --------------------------------------------------------------------------- planning

def test_along_swath_is_tried_first_and_costs_zero_links(two_line_tile):
    """When a target line itself covers the mark there is no cross-swath transfer."""
    inv = K.build_inventory([two_line_tile])
    g = K.overlap_graph(inv)
    paths = K.plan_path(g, inv, source_lines=[11], target_lines=[10, 11])
    assert len(paths) == 1
    p = paths[0]
    assert p.along_swath and p.n_links == 0 and p.nodes == [11] and p.edges == []
    assert "along-swath" in p.note


def test_the_shortest_chain_is_chosen_by_LINK_COUNT():
    g = {(1, 2): {"area_km2": 1.4}, (2, 3): {"area_km2": 1.4}, (3, 4): {"area_km2": 1.4},
         (4, 5): {"area_km2": 1.4}, (1, 9): {"area_km2": 0.01}, (9, 5): {"area_km2": 0.01}}
    inv = K.SwathInventory(tiles={})
    paths = K.plan_path(g, inv, source_lines=[1], target_lines=[5])
    assert paths[0].nodes == [1, 9, 5]      # 2 links, not the 4-link route
    assert paths[0].n_links == 2
    assert not paths[0].along_swath
    # every route returned is a SHORTEST route: the longer one is not offered at all
    assert {p.n_links for p in paths} == {2}


def test_all_shortest_paths_are_returned_so_redundancy_is_visible():
    g = {(1, 2): {"area_km2": 1.0}, (2, 4): {"area_km2": 1.0},
         (1, 3): {"area_km2": 1.0}, (3, 4): {"area_km2": 1.0}}
    inv = K.SwathInventory(tiles={})
    paths = K.plan_path(g, inv, source_lines=[1], target_lines=[4])
    assert {tuple(p.nodes) for p in paths} == {(1, 2, 4), (1, 3, 4)}


def test_a_disconnected_graph_returns_no_path_rather_than_raising():
    g = {(1, 2): {"area_km2": 1.0}, (7, 8): {"area_km2": 1.0}}
    inv = K.SwathInventory(tiles={})
    assert K.plan_path(g, inv, source_lines=[1], target_lines=[8]) == []


def test_the_west_chain_shape_is_reproduced_from_a_line_geometry_graph():
    """The measured Elba adjacency: 128-129-...-138, plus the eastern run to 144."""
    g = {(a, a + 1): {"area_km2": 1.4} for a in range(128, 146)}
    inv = K.SwathInventory(tiles={})
    west = K.plan_path(g, inv, source_lines=[128], target_lines=[133, 134, 135, 136, 137, 138])
    assert west[0].nodes == [128, 129, 130, 131, 132, 133]
    assert west[0].n_links == 5
    assert {p.n_links for p in west} == {5}
    east = K.plan_path(g, inv, source_lines=[144], target_lines=[133, 134, 135, 136, 137, 138])
    assert east[0].nodes == [144, 143, 142, 141, 140, 139, 138]
    assert east[0].n_links == 6
    assert {p.n_links for p in east} == {6}


# ----------------------------------------------------------------- chain accumulation

def test_a_two_link_chain_accumulates_a_known_offset(tmp_path):
    """Lines 10 / 11 / 12, each 40 mm above the previous. Solving 12 into 10's frame must
    return -80 mm, and the per-link residuals must be visible."""
    p = _write_tile(tmp_path / "three.las",
                    {10: (0.0, 100.0, 0.0, 120.0, 0.000),
                     11: (60.0, 160.0, 0.0, 120.0, 0.040),
                     12: (120.0, 220.0, 0.0, 120.0, 0.080)})
    inv = K.build_inventory([p])
    g = K.overlap_graph(inv, res=2.0)
    assert set(g) == {(10, 11), (11, 12)}          # 10 and 12 do not overlap
    paths = K.plan_path(g, inv, source_lines=[12], target_lines=[10])
    assert paths[0].nodes == [12, 11, 10] and paths[0].n_links == 2
    sol = K.solve_chain(inv, paths[0], res=2.0, graph=g, reference="last")
    assert sol.reference_line == 10 and sol.far_line == 12
    assert sol.dz_total_m == pytest.approx(-0.080, abs=0.005)
    assert len(sol.links) == 2
    assert all(L.overlap_km2 > 0 for L in sol.links)
    assert all(L.nmad_after_m <= L.nmad_before_m for L in sol.links)
    assert sol.dz_sigma_m > 0
    rows = sol.table_rows()
    assert len(rows) == 2 and len(rows[0]) == len(K.ChainSolution.table_columns())


def test_reference_end_flips_the_sign(tmp_path):
    p = _write_tile(tmp_path / "two.las", {10: (0.0, 100.0, 0.0, 120.0, 0.0),
                                           11: (60.0, 160.0, 0.0, 120.0, 0.040)})
    inv = K.build_inventory([p])
    g = K.overlap_graph(inv, res=2.0)
    path = K.plan_path(g, inv, source_lines=[11], target_lines=[10])[0]
    last = K.solve_chain(inv, path, res=2.0, graph=g, reference="last")
    first = K.solve_chain(inv, path, res=2.0, graph=g, reference="first")
    assert last.dz_total_m == pytest.approx(-first.dz_total_m, abs=0.004)
    with pytest.raises(ValueError):
        K.solve_chain(inv, path, reference="middle")


def test_compare_paths_reports_disagreement_as_the_error_bar():
    class S:
        def __init__(self, dz, sig, ref, far, n):
            self.dz_total_m = dz; self.dz_sigma_m = sig
            self.reference_line = ref; self.far_line = far
            self.path = type("P", (), {"n_links": n})()
    out = K.compare_paths([S(0.010, 0.002, 133, 128, 5), S(0.031, 0.003, 133, 144, 6)])
    assert out["spread_mm"] == pytest.approx(21.0, abs=1e-6)
    assert out["n_paths"] == 2
    assert len(out["formal_sigma_mm"]) == 2
    assert np.isnan(K.compare_paths([S(0.010, 0.002, 133, 128, 5)])["spread_mm"])


def test_chain_params_name_the_repo_defaults(tmp_path):
    p = _write_tile(tmp_path / "two2.las", {10: (0.0, 100.0, 0.0, 120.0, 0.0),
                                            11: (60.0, 160.0, 0.0, 120.0, 0.040)})
    inv = K.build_inventory([p])
    g = K.overlap_graph(inv, res=2.0)
    path = K.plan_path(g, inv, source_lines=[11], target_lines=[10])[0]
    sol = K.solve_chain(inv, path, res=2.0, graph=g)
    names = {q.name: q for q in sol.params}
    assert names["link_res_m"].value == 2.0
    assert names["link_exclude_classes"].value == (5, 6, 9)
    assert "coregister_swaths" in names["link_exclude_classes"].why
    assert all(q.src == "repo" and q.why for q in sol.params)
