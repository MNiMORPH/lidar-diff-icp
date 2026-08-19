"""gen2 (3DEP) bare-earth selection: prefer the survey's ASPRS ground class,
fall back to CSF over an unclassified region.

These pin the behaviour that motivated the switch away from the last-return
heuristic: 3DEP's class 2 is a strict, cleaner subset of last-return (it drops
canopy/understory last hits), and an unclassified tile must transparently fall
back to CSF rather than silently returning a canopy-contaminated surface.
"""
import os

import numpy as np
import laspy
import pytest

from lidar_diff_icp import pipeline
from lidar_diff_icp.pipeline import read_after_ground


X0, Y0, W = 1000.0, 2000.0, 100.0
BOUNDS = (X0, Y0, X0 + W, Y0 + W)


def _write(path, x, y, z, cls, rn, nr):
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [x.min(), y.min(), z.min()]; hdr.scales = [0.01, 0.01, 0.01]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.classification = cls.astype(np.uint8)
    las.return_number = rn.astype(np.uint8)
    las.number_of_returns = nr.astype(np.uint8)
    las.point_source_id = np.ones(len(x), np.uint16)
    las.gps_time = y.astype(np.float64)
    las.write(str(path))


def _mixed_cloud(tmp_path):
    """500 class-2 ground last-returns + 300 class-5 veg last-returns + 200
    class-1 intermediate returns (rn<nr). last-return set = ground + veg (800);
    class-2 set = ground only (500)."""
    rng = np.random.default_rng(0)
    def xy(n): return rng.uniform(X0, X0 + W, n), rng.uniform(Y0, Y0 + W, n)
    gx, gy = xy(500); vx, vy = xy(300); ix, iy = xy(200)
    x = np.concatenate([gx, vx, ix]); y = np.concatenate([gy, vy, iy])
    z = np.concatenate([np.full(500, 100.0), np.full(300, 108.0), np.full(200, 105.0)])
    cls = np.concatenate([np.full(500, 2), np.full(300, 5), np.full(200, 1)])
    rn = np.concatenate([np.ones(500), np.ones(300), np.ones(200)])       # ground/veg are last (rn==nr)
    nr = np.concatenate([np.ones(500), np.ones(300), np.full(200, 2)])    # intermediates: rn=1<nr=2
    p = tmp_path / "gen2.laz"; _write(p, x, y, z, cls, rn, nr)
    return str(p)


def test_class2_is_clean_ground_subset(tmp_path):
    p = _mixed_cloud(tmp_path)
    g = read_after_ground(p, BOUNDS, mode="class2")
    lr = read_after_ground(p, BOUNDS, mode="last_return")
    assert g["ground_mode"] == "class2"
    assert g["x"].size == 500                       # only class-2 ground
    assert lr["x"].size == 800                       # ground + veg last returns
    assert np.allclose(g["z"], 100.0)                # ground only, no 108 m canopy
    assert g["x"].size < lr["x"].size                # class 2 is a strict subset


def test_unclassified_region_falls_back_to_csf(tmp_path, monkeypatch):
    """A tile with no class-2 ground triggers the region-level CSF fallback. We
    stub CSF (no PDAL in CI) with a ground-only file and assert it is used."""
    rng = np.random.default_rng(1)
    n = 400
    x = rng.uniform(X0, X0 + W, n); y = rng.uniform(Y0, Y0 + W, n)
    z = np.full(n, 100.0)
    _write(tmp_path / "unclassified.laz", x, y, z,
           np.ones(n), np.ones(n), np.ones(n))       # all class 1

    csf_dir = tmp_path / "csf"; csf_dir.mkdir()
    csf_out = csf_dir / "ground.las"
    _write(csf_out, x[:120], y[:120], np.full(120, 99.5),
           np.full(120, 2), np.ones(120), np.ones(120))   # CSF returns 120 class-2 pts

    def fake_csf(path, pdal=None):
        return str(csf_out)
    monkeypatch.setattr(pipeline, "classify_ground_csf", fake_csf)

    with pytest.warns(UserWarning, match="unclassified"):
        g = read_after_ground(str(tmp_path / "unclassified.laz"), BOUNDS, mode="class2")
    assert g["ground_mode"] == "csf_fallback"
    assert g["x"].size == 120
    assert np.allclose(g["z"], 99.5)                 # came from the CSF stub
    assert not os.path.exists(csf_out)               # temp CSF output cleaned up
