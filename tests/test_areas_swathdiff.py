"""The two library modules that no test reached: areas.py and swathdiff.py.

Task 43. The question was "test them or retire them", and the answer is TEST: both are
load-bearing. coherence.py lazily imports areas.flow_corridor_counts inside a function --
which is why a plain grep for importers found nothing -- and coreg.coregister_swaths, the
heart of the swath network, uses swathdiff._median_grid. Nine files import swathdiff.

Being load-bearing AND untested is the worst of the two states, and it lasted this long
because a lazy import hides from the same scans that hid everything else this session.
"""
import numpy as np
import pytest

from lidar_diff_icp import areas, swathdiff
from lidar_diff_icp.io import PointCloud, MN_GEN1_CRS


# --- areas.flow_corridor_counts ----------------------------------------------------------

def _line_flow(ny, nx):
    """A single flow line down column 0: each cell drains to the one below it."""
    down = np.full(ny * nx, -1, int)
    up = np.full(ny * nx, -1, int)
    for r in range(ny - 1):
        down[r * nx] = (r + 1) * nx
        up[(r + 1) * nx] = r * nx
    return down, up


def test_the_corridor_follows_FLOW_not_a_disc():
    """THE POINT OF THIS MODULE. Wheaton coherence over a disc promotes a cell because its
    NEIGHBOURS agree; along flow it promotes because its CATCHMENT agrees. A gully is a
    line, so a disc footprint dilutes it -- that is why coherence erased real gullies under
    ~20-25 m until the flow corridor existed."""
    ny = nx = 9
    dod = np.zeros((ny, nx)); dod[:, 0] = 0.5          # deposition along the flow line only
    valid = np.ones((ny, nx), bool)
    down, up = _line_flow(ny, nx)
    nd, ne, n = areas.flow_corridor_counts(dod, valid, down, up, k=4, width=0)

    # a cell ON the line sees the whole corridor as same-sign deposition
    assert nd[4, 0] == n[4, 0] and n[4, 0] > 1
    # a cell OFF the line has a corridor of itself alone: no flow neighbours
    assert n[4, 5] == 1
    assert ne.sum() == 0, "nothing is erosion here"


def test_a_cell_not_built_gets_weight_zero_rather_than_a_nan():
    """`cells` restricts which footprints are built. An unbuilt cell keeps COUNT 0, and its
    footprint size stays at the nominal 2k+1 rather than 0 -- which looks wrong and is not.

    coherence.spatial_coherence_probability sets low = 0.80*n, up = 0.96*n and takes
    (ndepos - low) / (up - low). With count 0 that clips to weight 0, so the cell is not
    promoted, which is the documented contract. With n = 0 instead, low == up and the
    division is 0/0 -- every unbuilt cell would come back NaN. The initialisation is
    protective, and this test records that so nobody "fixes" it."""
    ny = nx = 7
    dod = np.full((ny, nx), 0.4)
    valid = np.ones((ny, nx), bool)
    down, up_ = _line_flow(ny, nx)
    only = np.zeros((ny, nx), bool); only[3, 0] = True
    nd, ne, n = areas.flow_corridor_counts(dod, valid, down, up_, k=3, cells=only)
    assert nd[3, 0] > 0, "the one requested footprint was built"
    assert nd.sum() == nd[3, 0], "and no other was"
    assert (n[nd == 0] == 2 * 3 + 1).all(), "unbuilt cells keep the nominal footprint size"
    # and that combination really does yield weight 0, not NaN
    low, hi = 0.80 * n, 0.96 * n
    w = np.clip((nd - low) / (hi - low), 0.0, 1.0)
    assert np.isfinite(w).all() and w[2, 3] == 0.0


def test_invalid_cells_are_not_counted():
    ny = nx = 7
    dod = np.full((ny, nx), 0.4)
    valid = np.ones((ny, nx), bool); valid[:, 0] = False      # the flow line is invalid
    down, up = _line_flow(ny, nx)
    nd, ne, n = areas.flow_corridor_counts(dod, valid, down, up, k=3)
    assert nd[:, 0].sum() == 0


# --- swathdiff.swath_difference ----------------------------------------------------------

def _two_swaths(dz=0.05, seed=0):
    """Two overlapping swaths over one surface, the second raised by a known dz."""
    rng = np.random.default_rng(seed)
    n = 20000
    xa = rng.uniform(0, 60, n); ya = rng.uniform(0, 60, n)
    xb = rng.uniform(40, 100, n); yb = rng.uniform(0, 60, n)
    f = lambda x, y: 100.0 + 0.05 * x + 0.02 * y
    za = f(xa, ya) + rng.normal(0, 0.01, n)
    zb = f(xb, yb) + dz + rng.normal(0, 0.01, n)
    x = np.r_[xa, xb]; y = np.r_[ya, yb]; z = np.r_[za, zb]
    ps = np.r_[np.full(n, 1), np.full(n, 2)].astype(np.int32)
    return PointCloud(x=x, y=y, z=z, point_source_id=ps,
                      classification=np.full(x.size, 2, np.uint8),
                      gps_time=np.zeros(x.size), scan_angle=np.zeros(x.size),
                      crs=MN_GEN1_CRS)


def test_a_known_offset_between_swaths_is_recovered():
    """The self-calibration check: between two swaths flown minutes apart there is no real
    land-surface change, so whatever the difference surface shows is the instrument."""
    d = swathdiff.swath_difference(_two_swaths(dz=0.05), 1, 2, res=2.0)
    assert d.median_offset == pytest.approx(-0.05, abs=0.005), d.median_offset
    assert d.n_cells > 0 and d.diff.shape[0] > 1
    assert d.robust_std < 0.02, "one surface plus 1 cm noise: the spread is the noise"


def test_non_overlapping_swaths_refuse_rather_than_return_nothing():
    pc = _two_swaths()
    far = PointCloud(x=pc.x + np.where(pc.point_source_id == 2, 1000.0, 0.0), y=pc.y,
                     z=pc.z, point_source_id=pc.point_source_id,
                     classification=pc.classification, gps_time=pc.gps_time,
                     scan_angle=pc.scan_angle, crs=pc.crs)
    with pytest.raises(ValueError, match="do not overlap"):
        swathdiff.swath_difference(far, 1, 2, res=2.0)


def test_the_median_grid_is_robust_to_a_wild_outlier():
    """_median_grid is what coreg.coregister_swaths builds every pairwise tie on, so its
    robustness is the swath network's robustness."""
    x = np.array([1.0, 1.1, 1.2, 1.3]); y = np.array([1.0, 1.1, 1.2, 1.3])
    z = np.array([10.0, 10.0, 10.0, 900.0])                  # one absurd return
    g = swathdiff._median_grid(x, y, z, 5.0, 0.0, 0.0, 1, 1)
    assert g[0, 0] == pytest.approx(10.0, abs=0.5), g[0, 0]
