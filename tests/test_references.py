"""Flat-hard stable-reference detection recovers a known vertical datum, keeps hard
flat surfaces, and rejects rough ones."""
import numpy as np

from lidar_diff_icp.references import (flat_hard_cells, datum_offset, datum_plane,
                                       eval_datum_correction)


def test_datum_plane_recovers_constant_and_tilt_rejects_outliers():
    rng = np.random.default_rng(5)
    n = 400; cx = cy = 1500.0
    x = rng.uniform(0, 3000, n); y = rng.uniform(0, 3000, n)
    a, b, c = -0.055, 0.010, -0.006                    # const (m), tilt (m/km)
    E = (x - cx) / 1000; N = (y - cy) / 1000
    off = a + b * E + c * N + rng.normal(0, 0.02, n)
    off[:20] += 0.30                                   # resurfacing outliers
    cells = dict(x=x, y=y, offset=off, roughness=np.full(n, 0.01))
    pl = datum_plane(cells)
    assert abs(pl["a"] - a) < 0.01                      # constant
    assert abs(pl["b"] - b) < 0.004 and abs(pl["c"] - c) < 0.004   # tilt
    assert pl["rejected"] >= 15                          # outliers dropped
    corr = eval_datum_correction(pl, x[20:], y[20:])     # apply to before -> stable ~0
    assert abs(np.median(off[20:] + corr)) < 0.005


def _flat_surface(x0, x1, y0, y1, z, sigma, seed, per_cell=12, res=2.0):
    """Dense points on a flat surface z (+ noise sigma) over a rectangle."""
    rng = np.random.default_rng(seed)
    nx = int((x1 - x0) / res); ny = int((y1 - y0) / res); n = nx * ny * per_cell
    x = rng.uniform(x0, x1, n); y = rng.uniform(y0, y1, n)
    zz = z + rng.normal(0, sigma, n)
    return x, y, zz


def test_recovers_known_datum_on_flat_hard_surface():
    bounds = (0, 0, 120, 120)
    dz_true = -0.055                                   # gen1 sits 55 mm below gen2
    # hard flat surface (low roughness), both epochs, offset by dz_true
    ax, ay, az = _flat_surface(10, 110, 10, 110, 10.00, 0.004, 1)
    bx, by, bz = _flat_surface(10, 110, 10, 110, 10.00 + dz_true, 0.004, 2)
    cells = flat_hard_cells(bx, by, bz, ax, ay, az, bounds)
    assert cells["x"].size > 500                       # plenty of flat-hard cells found
    d = datum_offset(cells)
    assert abs(d["raw"] - dz_true) < 0.005             # datum recovered to <5 mm
    assert abs(d["dz_before"] - (-dz_true)) < 0.005    # correction is +55 mm to gen1
    assert d["se"] < 0.003


def test_rejects_rough_surface():
    """A rough surface (high within-cell spread) is not a hard reference."""
    bounds = (0, 0, 80, 80)
    ax, ay, az = _flat_surface(5, 75, 5, 75, 10.0, 0.20, 3, per_cell=30)  # rough (20 cm) -> not hard
    bx, by, bz = _flat_surface(5, 75, 5, 75, 10.0, 0.20, 4, per_cell=30)
    cells = flat_hard_cells(bx, by, bz, ax, ay, az, bounds, max_rough_m=0.04)
    assert cells["x"].size < 3                            # essentially all rejected as non-hard
