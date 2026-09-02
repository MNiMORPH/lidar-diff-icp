"""read_tile must not invent a scan angle or a timestamp.

Zero is a MEASUREMENT for both: scan angle 0 is nadir, gps_time 0 is a time. The old
fallback wrote exactly those when the dimension was absent, and point format 6+ drops
`scan_angle_rank` in favour of a scaled `scan_angle` -- so every PF7 file, including the
CSF caches this project writes, read as "every return at nadir, at time zero".
"""
import numpy as np
import pytest

laspy = pytest.importorskip("laspy")

from lidar_diff_icp import io as L


def _write(path, pf, n=8, **dims):
    hdr = laspy.LasHeader(version="1.4" if pf >= 6 else "1.2", point_format=pf)
    hdr.scales = [0.001] * 3
    hdr.offsets = [0.0] * 3
    las = laspy.LasData(hdr)
    las.x = np.arange(n, dtype=float)
    las.y = np.zeros(n)
    las.z = np.zeros(n)
    las.point_source_id = np.ones(n, dtype=np.uint16)
    las.classification = np.full(n, 2, dtype=np.uint8)
    for k, v in dims.items():
        setattr(las, k, v)
    las.write(str(path))
    return path


def test_pf7_scan_angle_is_read_not_zeroed(tmp_path):
    """PF6+ stores scan_angle in 0.006-deg units and has no scan_angle_rank."""
    n = 8
    deg = np.linspace(-17.0, 17.0, n)
    p = _write(tmp_path / "pf7.las", 7, n=n,
               scan_angle=np.round(deg / 0.006).astype(np.int16),
               gps_time=np.arange(n, dtype=float) + 1000.0)
    pc = L.read_tile(p)
    sa = np.asarray(pc.scan_angle, float)
    assert np.isfinite(sa).all()
    assert not np.allclose(sa, 0.0), "the whole point: PF7 must not read as all-nadir"
    assert sa.min() == pytest.approx(-17.0, abs=0.01)
    assert sa.max() == pytest.approx(+17.0, abs=0.01)


def test_pf1_scan_angle_rank_is_already_degrees(tmp_path):
    n = 5
    rank = np.array([-17, -8, 0, 8, 17], dtype=np.int8)
    p = _write(tmp_path / "pf1.las", 1, n=n, scan_angle_rank=rank,
               gps_time=np.arange(n, dtype=float))
    pc = L.read_tile(p)
    assert np.allclose(np.asarray(pc.scan_angle, float), rank.astype(float))


def test_the_two_formats_agree_on_the_same_angles(tmp_path):
    """A PF1 and a PF7 file describing the same beam geometry must read the same."""
    n = 5
    deg = np.array([-17.0, -8.0, 0.0, 8.0, 17.0])
    a = L.read_tile(_write(tmp_path / "a.las", 1, n=n,
                           scan_angle_rank=deg.astype(np.int8)))
    b = L.read_tile(_write(tmp_path / "b.las", 7, n=n,
                           scan_angle=np.round(deg / 0.006).astype(np.int16)))
    assert np.allclose(np.asarray(a.scan_angle, float),
                       np.asarray(b.scan_angle, float), atol=0.01)


def test_absent_gps_time_is_nan_not_zero(tmp_path):
    """PF0 carries no gps_time; PF1 is the format that adds it."""
    p = _write(tmp_path / "pf0.las", 0, n=6)
    gt = np.asarray(L.read_tile(p).gps_time, float)
    assert np.isnan(gt).all(), "absent gps_time must be NaN, never 0 (= a timestamp)"


def test_every_standard_format_yields_a_real_scan_angle(tmp_path):
    """The NaN fallback for scan angle is a guard, not a live path.

    PF0-5 carry scan_angle_rank and PF6-10 carry scan_angle, so one of the two branches
    always fires for standard LAS. That matters for reading the original bug correctly: it
    was never "the dimension is missing", it was looking ONLY for scan_angle_rank, which
    PF6+ does not have -- so PF7 files silently took the zero fallback.
    """
    for pf in (0, 1, 6, 7):
        pc = L.read_tile(_write(tmp_path / f"f{pf}.las", pf, n=4))
        sa = np.asarray(pc.scan_angle, float)
        assert np.isfinite(sa).all(), f"PF{pf} should resolve a scan angle, not NaN"
