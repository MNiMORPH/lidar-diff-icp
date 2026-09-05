"""Gridding the reference (gen2) cloud: what `spread` and `rough` each mean, and what the
streaming route can and cannot hand back.

This stage was inline in difference_dem, so none of it could be reached without a real
tile. The property worth pinning is not that it runs but that its two roughness-like
outputs are DIFFERENT NUMBERS: `spread` is taken about the cell's own horizontal level and
therefore counts the slope itself, while `rough` is the residual to a fitted plane and does
not. The ground-q curve is indexed by `spread` and the LoD's standard error is built from
`rough`; swapping them would be silent and would bias every corrected elevation on a slope.
"""
import numpy as np
import laspy
import pytest

from lidar_diff_icp import pipeline


GRID = (0.0, 0.0, 5.0, 24, 20)          # X0, Y0, res, nx, ny
BOUNDS = (0.0, 0.0, 24 * 5.0, 20 * 5.0)


def _las(tmp_path, *, name="after.las", slope=0.0, noise=0.05, n=200_000, seed=7,
         classification=2):
    """A gen2-like ground cloud: a tilted plane plus Gaussian noise, all one class."""
    X0, Y0, res, nx, ny = GRID
    rng = np.random.default_rng(seed)
    x = rng.uniform(X0, X0 + nx * res, n)
    y = rng.uniform(Y0, Y0 + ny * res, n)
    z = 100.0 + slope * x + rng.normal(0.0, noise, n)
    hdr = laspy.LasHeader(version="1.4", point_format=6)
    hdr.offsets = [X0, Y0, 0.0]; hdr.scales = [0.001, 0.001, 0.001]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.classification = np.full(n, classification, np.uint8)
    las.return_number = np.ones(n, np.uint8)
    las.number_of_returns = np.ones(n, np.uint8)
    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / name
    las.write(str(p))
    return str(p)


def _ref(path, stream):
    return pipeline.grid_reference(path, BOUNDS, GRID, 0.50, stream=stream,
                                   after_ground="class2", verbose=False)


@pytest.mark.parametrize("stream", [True, False])
def test_spread_counts_the_slope_and_rough_does_not(tmp_path, stream):
    """THE GUARD THIS STAGE EXISTS FOR. On a 15% slope a 5 m cell spans 0.75 m, so the
    within-cell SD about a horizontal level is dominated by the slope, while the residual
    to a fitted plane still sees only the 50 mm of noise."""
    flat = _ref(_las(tmp_path, name="flat.las", slope=0.0, noise=0.05), stream)
    tilt = _ref(_las(tmp_path, name="tilt.las", slope=0.15, noise=0.05), stream)

    # spread grows with the slope ...
    assert np.nanmedian(tilt["spread"]) > 4 * np.nanmedian(flat["spread"])
    # ... rough does not: it stays at the noise level on both
    assert np.nanmedian(tilt["rough"]) == pytest.approx(np.nanmedian(flat["rough"]),
                                                        rel=0.35)
    assert np.nanmedian(tilt["rough"]) < 0.5 * np.nanmedian(tilt["spread"])


@pytest.mark.parametrize("stream", [True, False])
def test_count_is_the_number_of_ground_returns_in_the_cell(tmp_path, stream):
    p = _las(tmp_path, n=200_000)
    r = _ref(p, stream)
    inside = pipeline.read_after_ground(p, BOUNDS, mode="class2")["x"].size
    assert np.nansum(r["count"]) == inside       # every point in bounds is counted once
    assert inside == pytest.approx(200_000, abs=20)   # a few land exactly on the edge
    assert np.all(np.isfinite(r["count"]))       # every cell is populated at this density


@pytest.mark.parametrize("stream", [True, False])
def test_ground_recovers_a_known_level(tmp_path, stream):
    """q = 0.50 on symmetric noise about a flat plane returns the plane."""
    r = _ref(_las(tmp_path, slope=0.0, noise=0.05), stream)
    assert np.nanmedian(r["ground"]) == pytest.approx(100.0, abs=0.005)


def test_the_streaming_route_cannot_hand_back_the_cloud(tmp_path):
    """Why ground='plane'/'poly2' refuse under stream=True: those estimators need random
    access to the points, and the streaming route never holds them."""
    p = _las(tmp_path)
    assert _ref(p, stream=True)["cloud"] is None
    assert _ref(p, stream=True)["route"] == "stream"
    mem = _ref(p, stream=False)
    assert mem["route"] == "memory" and mem["cloud"] is not None
    assert mem["cloud"]["x"].size == pytest.approx(200_000, abs=20)


def test_the_two_routes_agree_on_the_ground_they_grid(tmp_path):
    """They are not bit-identical -- the streaming route bins the column -- but they must
    agree far inside the tolerance any conclusion rests on."""
    p = _las(tmp_path, slope=0.05)
    a, b = _ref(p, stream=True), _ref(p, stream=False)
    assert np.array_equal(a["count"], b["count"], equal_nan=True)
    d = np.abs(a["ground"] - b["ground"])
    assert np.nanmax(d) < 0.01, float(np.nanmax(d))       # under one bin
