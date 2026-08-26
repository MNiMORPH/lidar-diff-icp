"""The cross-line estimator must recover per-line coefficients a N-S pair cannot.

Two claims are load-bearing in `analysis/CROSS_LINE_FIT.md`, and both are tested against
data whose answer is known by construction rather than against the field data:

1. **A there-and-back N-S pair cannot separate the two lines.** Its two tangents sum to a
   near-constant, so ``D = k + p*(tA-tB) + q*(tA+tB)`` is near-singular in ``q`` and the
   individual ``c_A = p+q``, ``c_B = p-q`` are not estimable -- while the SUM ``p`` is.
   A cross pair, whose tangents are independent, recovers both.
2. **The slope-normal ground estimator removes the sub-cell tilt term**, so a between-line
   difference returns the injected per-line offsets and not the terrain gradient times
   however the two lines' returns happened to fall inside the cell.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from crossline_fit import fit_pair, cross_pair_rows, gen1_cellline_ground  # noqa: E402

C_A, C_B, K = 150.0, -40.0, 12.0          # mm per unit tangent, mm


def _pair_frame(ta, tb, rng, noise=8.0):
    D = K + C_A * ta - C_B * tb + rng.normal(0, noise, len(ta))
    n = len(ta)
    return pd.DataFrame(dict(D=D, ta=ta, tb=tb, dtan=ta - tb, stan=ta + tb,
                             blk=np.arange(n) // 25, cx=rng.random(n) * 2000,
                             cy=rng.random(n) * 2000))


def _ns_tangents(rng, n=6000, S=0.37):
    """There-and-back geometry: tan_A + tan_B is pinned to S = spacing/height."""
    ta = rng.uniform(0.05, 0.30, n)
    return ta, S - ta + rng.normal(0, 0.02, n)      # the measured sd(sum) is 0.017-0.024


def _cross_tangents(rng, n=6000):
    """A crossing line: its across-track coordinate is the other line's ALONG-track one."""
    return rng.uniform(-0.29, 0.29, n), rng.uniform(-0.31, 0.29, n)


def test_cross_pair_recovers_both_coefficients():
    rng = np.random.default_rng(1)
    ta, tb = _cross_tangents(rng)
    f = fit_pair(_pair_frame(ta, tb, rng))
    assert f["se_ratio"] < 1.5                       # the two are separable
    assert abs(f["c_A"] - C_A) < 4 * f["c_A_se"]
    assert abs(f["c_B"] - C_B) < 4 * f["c_B_se"]
    assert f["c_A_se"] < 15 and f["c_B_se"] < 15


def test_ns_pair_cannot_separate_them_but_gets_the_sum():
    rng = np.random.default_rng(2)
    ta, tb = _ns_tangents(rng)
    f = fit_pair(_pair_frame(ta, tb, rng))
    assert f["se_ratio"] > 4.0                       # the degeneracy, quantified
    # each line's OWN coefficient costs several times the standard error of their SUM,
    # which is the whole reason the N-S chain reports sums
    assert f["c_A_se"] / f["c_pair_se"] > 4.0
    assert f["c_B_se"] / f["c_pair_se"] > 4.0
    # ...while the pair SUM (c_A + c_B)/2 is measured tightly by the same data
    assert abs(f["c_pair"] - 0.5 * (C_A + C_B)) < 4 * f["c_pair_se"]


def test_the_two_designs_differ_only_in_the_tangent_geometry():
    """Same coefficients, same noise, same n -- only the flight geometry changes."""
    rng = np.random.default_rng(3)
    ns = fit_pair(_pair_frame(*_ns_tangents(rng), rng))
    cr = fit_pair(_pair_frame(*_cross_tangents(rng), rng))
    assert ns["c_A_se"] / cr["c_A_se"] > 5
    assert ns["se_ratio"] / cr["se_ratio"] > 5


def _synthetic_ground_las(path, *, offsets, slope=0.15, res=5.0, n_per_cell=14, seed=0):
    """A tilted plane sampled by several flight lines, each with a known vertical offset
    and its own scan-angle field, written as a class-2 LAS."""
    import laspy
    rng = np.random.default_rng(seed)
    xs, ys, zs, ps, sa = [], [], [], [], []
    ncell = 40
    for psid, (dz_mm, scan) in offsets.items():
        for i in range(ncell):
            for j in range(ncell):
                x = 577000 + (i + rng.random(n_per_cell)) * res
                y = 4886000 + (j + rng.random(n_per_cell)) * res
                z = 200.0 + slope * (x - 577000) + 0.05 * (y - 4886000) + dz_mm / 1000.0
                xs.append(x); ys.append(y); zs.append(z)
                ps.append(np.full(n_per_cell, psid))
                sa.append(np.full(n_per_cell, scan))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    hdr = laspy.LasHeader(point_format=7, version="1.4")
    hdr.scales = [0.001, 0.001, 0.001]
    hdr.offsets = [x.min(), y.min(), z.min()]
    las = laspy.LasData(hdr)
    las.x = x; las.y = y; las.z = z
    las.point_source_id = np.concatenate(ps).astype(np.uint16)
    las.scan_angle = (np.concatenate(sa) / 0.006).astype(np.int16)
    las.classification = np.full(len(x), 2, np.uint8)
    las.return_number = np.ones(len(x), np.uint8)
    las.number_of_returns = np.ones(len(x), np.uint8)
    las.write(path)


def test_ground_estimator_returns_the_injected_between_line_offset(tmp_path):
    """On a 15% slope, the between-line difference must be the injected offset difference.

    This is what the slope-normal residual buys: a per-cell median of RAW z would carry
    the horizontal sampling difference between the two lines times the terrain gradient.
    """
    p = str(tmp_path / "synth.las")
    _synthetic_ground_las(p, offsets={136: (0.0, +8.0), 137: (-60.0, -7.0)})
    clm, nx, X0, Y0, n = gen1_cellline_ground(p, 5.0, 1.2, 0.50)
    m = cross_pair_rows(clm, nx, 5.0, 50.0)
    assert len(m) > 1000
    assert abs(np.median(m.D) - 60.0) < 3.0          # 136 reads 60 mm above 137


def test_ground_estimator_refuses_a_zeroed_scan_angle(tmp_path):
    """The failure that has bitten this project twice: PDAL renames scan_angle_rank ->
    scan_angle on the PF1 -> PF7 promotion, and reading the wrong name gives silent zeros."""
    p = str(tmp_path / "flat.las")
    _synthetic_ground_las(p, offsets={136: (0.0, 0.0), 137: (-60.0, 0.0)})
    with pytest.raises(RuntimeError, match="identically zero"):
        gen1_cellline_ground(p, 5.0, 1.2, 0.50)
