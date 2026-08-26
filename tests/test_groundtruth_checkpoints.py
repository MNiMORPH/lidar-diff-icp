"""Checkpoint ingestion: the datum must be stated, not assumed."""
import csv

import pytest

from lidar_diff_icp.groundtruth import checkpoints as C


def _row(**kw):
    r = dict(point_id="X1", point_type="NVA", easting="1.0", northing="2.0",
             elevation="100.0", elevation_units="m", horizontal_crs="EPSG:6344",
             vertical_datum="NAVD88", geoid_model="GEOID18")
    r.update(kw)
    return r


def _write(tmp_path, rows, name="cp.csv"):
    p = tmp_path / name
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def test_bundled_set_loads_with_full_datum():
    s = C.load_bundled()
    assert len(s) == 6
    for p in s:
        assert p.geoid_model == "GEOID18"
        assert p.vertical_datum == "NAVD88"
        assert p.elevation_units == "m"
        assert p.horizontal_crs == "EPSG:6344"
    # the four probed in ELBAEXT2_SCOPE, with their surveyed heights
    assert s["2210_2021_MN"].elevation_m == pytest.approx(349.288)
    assert s["3056_2021_MN"].point_type == "VVA"
    assert s["2024_2021_MN"].elevation_m == pytest.approx(344.735)
    assert s["2036_2021_MN"].elevation_m == pytest.approx(353.119)


def test_bundled_set_is_usable_and_names_its_origin():
    s = C.load_bundled()
    assert len(s.usable()) == 6
    assert s.origin.endswith(".csv")


@pytest.mark.parametrize("token", ["", "unknown", "N/A", "  ", "None"])
def test_unknown_geoid_model_refuses(tmp_path, token):
    p = _write(tmp_path, [_row(geoid_model=token)])
    s = C.read_checkpoint_csv(p)
    assert not s[0].datum_known
    with pytest.raises(C.UnknownDatumError, match="geoid"):
        s.usable()


def test_unknown_vertical_datum_refuses(tmp_path):
    p = _write(tmp_path, [_row(vertical_datum="")])
    with pytest.raises(C.UnknownDatumError):
        C.read_checkpoint_csv(p).usable()


def test_a_good_point_is_not_dropped_when_a_bad_one_is_present(tmp_path):
    """The bad point must RAISE, not vanish: a control point that cannot be used is a
    reported result, never a silent exclusion."""
    p = _write(tmp_path, [_row(point_id="good"), _row(point_id="bad", geoid_model="")])
    s = C.read_checkpoint_csv(p)
    assert s.ids == ["good", "bad"]          # both present after reading
    with pytest.raises(C.UnknownDatumError, match="bad"):
        s.usable()


def test_non_metre_units_refuse_rather_than_convert(tmp_path):
    p = _write(tmp_path, [_row(elevation_units="US Feet")])
    cp = C.read_checkpoint_csv(p)[0]
    with pytest.raises(ValueError, match="does not convert"):
        cp.elevation_m


def test_missing_column_is_an_error_not_a_default(tmp_path):
    rows = [_row()]
    del rows[0]["geoid_model"]
    p = _write(tmp_path, rows)
    with pytest.raises(ValueError, match="geoid_model"):
        C.read_checkpoint_csv(p)


def test_within_keeps_the_datum_metadata():
    s = C.load_bundled().within((570000, 4884000, 571000, 4885000))
    assert s.ids == ["2210_2021_MN", "3056_2021_MN"]
    assert all(p.geoid_model == "GEOID18" for p in s)


def test_shapefile_reader_requires_units_to_be_stated():
    """``source_ele_units`` has no default: the contractor .dbf mislabels metres as feet."""
    import inspect
    sig = inspect.signature(C.read_3dep_va_shapefile)
    p = sig.parameters["source_ele_units"]
    assert p.default is inspect.Parameter.empty
    assert p.kind is inspect.Parameter.KEYWORD_ONLY


def test_shapefile_reader_refuses_a_file_with_no_geoid_field(tmp_path):
    shapefile = pytest.importorskip("shapefile")
    p = tmp_path / "va"
    w = shapefile.Writer(str(p))
    w.field("unique_ind", "C")
    w.field("source_ele", "N", decimal=3)
    w.point(570000.0, 4884000.0)
    w.record("2210", 349.288)
    w.close()
    with pytest.raises(C.UnknownDatumError, match="geoid"):
        C.read_3dep_va_shapefile(p, source_ele_units="m")


def test_shapefile_reader_reads_the_geoid_field_when_present(tmp_path):
    shapefile = pytest.importorskip("shapefile")
    p = tmp_path / "va2"
    w = shapefile.Writer(str(p))
    w.field("unique_ind", "C")
    w.field("source_ele", "N", decimal=3)
    w.field("nva_vva", "C")
    w.field("geoid", "C")
    w.point(570492.1, 4884126.1)
    w.record("2210_2021_MN", 349.288, "NVA", "GEOID18")
    w.close()
    s = C.read_3dep_va_shapefile(p, source_ele_units="m")
    assert len(s) == 1
    cp = s.usable()[0]
    assert cp.point_id == "2210_2021_MN"
    assert cp.geoid_model == "GEOID18"
    assert cp.point_type == "NVA"
    assert cp.elevation_m == pytest.approx(349.288)
