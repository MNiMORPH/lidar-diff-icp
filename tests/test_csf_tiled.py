"""The tiled CSF's core assignment must PARTITION the cloud -- exactly once each.

Tiling CSF is a memory device, not a different computation: every point is classified
inside some tile's halo, and exactly one tile is allowed to keep it. Two failures are
possible and both are silent in the output LAS -- a seam GAP (points on a tile boundary
kept by nobody, so the ground cache has stripes of missing ground) and a seam DUPLICATE
(points kept twice, so those cells are over-weighted in every per-cell median downstream).

The half-open [cx0, cx1) rule handles interior seams; the ``i == 0`` / ``i == nx-1``
widenings handle the outer edges, where floating-point on ``x.min()``/``x.max()`` can put
a point just outside its own domain. Both are tested here.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis", "slope_bias"))
from csf_tiled import core_mask  # noqa: E402


def _partition_counts(x, y, nx, ny, overlap=150.0):
    """How many tile cores claim each point, under the real halo/core construction."""
    X0, X1, Y0, Y1 = x.min(), x.max(), y.min(), y.max()
    dx, dy = (X1 - X0) / nx, (Y1 - Y0) / ny
    claims = np.zeros(len(x), int)
    for j in range(ny):
        for i in range(nx):
            cx0, cx1 = X0 + i * dx, X0 + (i + 1) * dx
            cy0, cy1 = Y0 + j * dy, Y0 + (j + 1) * dy
            halo = ((x >= cx0 - overlap) & (x <= cx1 + overlap) &
                    (y >= cy0 - overlap) & (y <= cy1 + overlap))
            core = core_mask(x[halo], y[halo], i, j, nx, ny, cx0, cx1, cy0, cy1)
            idx = np.flatnonzero(halo)[core]
            claims[idx] += 1
    return claims


@pytest.mark.parametrize("nx,ny", [(1, 1), (2, 2), (3, 3), (3, 2)])
def test_cores_partition_a_random_cloud(nx, ny):
    rng = np.random.default_rng(0)
    x = 577_000 + rng.random(20_000) * 2_500
    y = 4_886_000 + rng.random(20_000) * 3_500
    assert np.all(_partition_counts(x, y, nx, ny) == 1)


@pytest.mark.parametrize("nx,ny", [(2, 2), (3, 3)])
def test_points_exactly_on_the_seams_are_claimed_once(nx, ny):
    """Seam coordinates are the failure case: an inclusive [cx0, cx1] rule double-counts
    them, and a doubly-exclusive one drops them."""
    X0, X1 = 577_000.0, 579_500.0
    Y0, Y1 = 4_886_000.0, 4_889_500.0
    dx, dy = (X1 - X0) / nx, (Y1 - Y0) / ny
    seam_x = np.array([X0 + i * dx for i in range(nx + 1)])
    seam_y = np.array([Y0 + j * dy for j in range(ny + 1)])
    gx, gy = np.meshgrid(seam_x, seam_y)
    # corners alone would not pin X0/X1 as the extent, so include them plus interior pts
    x = np.r_[gx.ravel(), X0, X1, (X0 + X1) / 2]
    y = np.r_[gy.ravel(), Y0, Y1, (Y0 + Y1) / 2]
    assert np.all(_partition_counts(x, y, nx, ny) == 1)


def test_an_inclusive_upper_edge_would_double_count():
    """Prove the test bites: relax the half-open rule and the seam points are claimed twice."""
    X0, X1, Y0, Y1 = 0.0, 100.0, 0.0, 100.0
    nx = ny = 2
    dx, dy = (X1 - X0) / nx, (Y1 - Y0) / ny
    x = np.array([0.0, 50.0, 100.0, 50.0, 25.0])
    y = np.array([0.0, 50.0, 100.0, 25.0, 50.0])
    claims = np.zeros(len(x), int)
    for j in range(ny):
        for i in range(nx):
            cx0, cx1 = X0 + i * dx, X0 + (i + 1) * dx
            cy0, cy1 = Y0 + j * dy, Y0 + (j + 1) * dy
            bad = (x >= cx0) & (x <= cx1) & (y >= cy0) & (y <= cy1)   # inclusive: WRONG
            claims += bad
    assert claims.max() > 1                      # the bug the real rule avoids
    assert np.all(_partition_counts(x, y, nx, ny) == 1)   # the real rule is clean
