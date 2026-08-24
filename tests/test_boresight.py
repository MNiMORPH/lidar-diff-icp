"""Boresight self-calibration: recover an injected roll from synthetic flight-line overlap.

The synthetic data carry a large per-cell shared component (terrain + datum, ~200 mm) that
MUST be differenced out, per-line registration offsets, and a known roll b0 proportional to
scan angle. Two opposing lines sample each cell at opposite scan angles that vary cross-track,
so the between-line scan difference has leverage. The estimator must return b0 (not 0, not the
shared component), apply must remove it, and with no roll it must return ~0.
"""
import numpy as np
import pytest

from lidar_diff_icp.boresight import estimate_boresight, apply_boresight


def _synth(b0=2.5, reg=(("10", -8.0), ("20", 8.0)), ncells=2000, per=5, noise=15.0, seed=0):
    reg = dict(reg)
    rng = np.random.default_rng(seed)
    xt = rng.uniform(-1, 1, ncells)                      # cross-track position of each cell
    shared = rng.normal(0, 200.0, ncells)               # terrain change + datum (differenced out)
    base = {"10": 15.0 * xt, "20": -13.0 * xt}          # scan angle vs cross-track, opposing lines
    cells, psid, scan, d = [], [], [], []
    for ci in range(ncells):
        for ln in reg:
            for _ in range(per):
                s = base[ln][ci] + rng.normal(0, 0.5)
                cells.append(ci); psid.append(ln); scan.append(s)
                d.append(shared[ci] + reg[ln] + b0 * s + rng.normal(0, noise))
    return (np.array(cells), np.array(psid, dtype=object),
            np.array(scan, float), np.array(d, float))


def test_recovers_injected_roll():
    cell, psid, scan, d = _synth(b0=2.5)
    sol = estimate_boresight(cell, psid, scan, d, min_cell_line=3, n_boot=50)
    assert abs(sol.b - 2.5) < 0.2, f"roll {sol.b} not ~2.5"
    # per-line registration recovered (mean-zero gauge): 20 - 10 ~ +16 mm
    assert abs((sol.registration["20"] - sol.registration["10"]) - 16.0) < 3.0
    assert sol.b_se < 0.5


def test_apply_removes_the_roll():
    cell, psid, scan, d = _synth(b0=2.5, noise=1.0)
    sol = estimate_boresight(cell, psid, scan, d, min_cell_line=3, n_boot=20)
    d_corr = d - apply_boresight(scan, sol.b)
    sol2 = estimate_boresight(cell, psid, scan, d_corr, min_cell_line=3, n_boot=20)
    assert abs(sol2.b) < 0.2, f"residual roll {sol2.b} after correction"


def test_no_false_roll_when_absent():
    """The bite: a scan-blind estimator would report the shared/registration structure as a
    roll. With b0=0 the estimate must be ~0; with b0=3 it must be clearly non-zero."""
    cell, psid, scan, d = _synth(b0=0.0, noise=1.0)
    assert abs(estimate_boresight(cell, psid, scan, d, n_boot=20).b) < 0.2
    cell, psid, scan, d = _synth(b0=3.0, noise=1.0)
    assert estimate_boresight(cell, psid, scan, d, n_boot=20).b > 1.5
