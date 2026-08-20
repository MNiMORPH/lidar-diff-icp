"""Mass-conserving DoD budget: the accumulation must conserve mass exactly, and the
routed continuity check must flag unphysical (unsupported) deposition while passing
physical erosion-then-deposition."""
import numpy as np
import pytest

from lidar_diff_icp.massbalance import weighted_accumulation, mass_balance


def _straight_down_props(ny, nx):
    """Hand-built D-infinity proportions: every cell sends all flow straight down
    (band 7 = 1.0); bottom row exits the grid. band0 = 0 (all valid)."""
    props = np.full((ny, nx, 9), -1.0)
    props[:, :, 0] = 0.0            # all resolved/valid
    props[:, :, 7] = 1.0            # DX[7]=0, DY[7]=+1 -> straight down
    props[:, :, [1, 2, 3, 4, 5, 6, 8]] = -1.0
    return props, np.ones((ny, nx), bool)


def test_weighted_accumulation_conserves_mass():
    props, valid = _straight_down_props(4, 3)
    w = np.zeros((4, 3)); w[0, 1] = 1.0                 # unit source, top-middle
    acc, exited = weighted_accumulation(w, props, valid)
    assert np.allclose(acc[:, 1], [1, 1, 1, 1])          # routes straight down its column
    assert np.allclose(acc[:, 0], 0) and np.allclose(acc[:, 2], 0)
    assert abs(exited - 1.0) < 1e-9                       # all mass leaves at the bottom edge

    # uniform weight: total exited == number of valid cells (nothing created/destroyed)
    acc2, ex2 = weighted_accumulation(valid.astype(float), props, valid)
    assert abs(ex2 - valid.sum()) < 1e-9


def test_mass_balance_flags_unsupported_deposition_passes_supported():
    """Column 1: erosion upstream then deposition downstream -> V_acc stays <= 0
    (supported, no flag). Column 2 gets deposition with NO upstream erosion ->
    V_acc climbs positive -> surplus flagged. (Columns are independent under
    straight-down flow, so they don't cross-contaminate.)"""
    ny, nx = 8, 4
    props, valid = _straight_down_props(ny, nx)
    res = 5.0; perror = np.full((ny, nx), 0.001)          # tiny error: don't mask the signal
    dod = np.zeros((ny, nx))
    dod[1:4, 1] = -0.20                                    # erosion upstream (col 1)
    dod[4:7, 1] = +0.10                                    # deposition downstream (col 1): supported
    dod[4:7, 2] = +0.10                                    # deposition (col 2) with no upstream erosion

    out = mass_balance(dod, perror, props, valid, res, z=1.96)
    V, surplus, contaminated = out["V_acc"], out["surplus"], out["contaminated"]

    # col 1 downstream budget never goes positive (erosion supplies the deposition)
    assert np.nanmax(V[:, 1]) <= 1e-9
    assert not surplus[:, 1].any()
    # col 2 accumulates unsupported deposition -> positive V_acc and a surplus flag,
    # among the cells the check actually evaluates (uncontaminated)
    evaluable = ~contaminated
    assert (surplus[:, 2] & evaluable[:, 2]).any() or contaminated[:, 2].all()
    if evaluable[4:7, 2].any():
        assert np.nanmax(V[4:7, 2]) > 0


def test_mass_balance_end_to_end_with_richdem():
    """Full path incl. D-infinity routing (RichDEM): a radial 'tent' DEM so flow
    originates at an interior peak (uncontaminated interior). A dug pit that only
    fills (deposition, no upstream erosion) must raise a surplus flag."""
    rd = pytest.importorskip("richdem")
    from lidar_diff_icp.massbalance import dinf_proportions
    n = 41; res = 5.0
    yy, xx = np.mgrid[0:n, 0:n]
    dem = 100 - np.hypot(yy - n // 2, xx - n // 2)          # cone, peak in the centre
    props, valid = dinf_proportions(dem, breach=True)
    perror = np.full((n, n), 0.002)
    dod = np.zeros((n, n))
    dod[n // 2 + 6, n // 2 + 6] = +0.5                       # pure deposition, no upstream erosion
    out = mass_balance(dod, perror, props, valid, res)
    # somewhere at/below that unsupported deposition the budget is positive & flagged
    assert np.nanmax(out["V_acc"]) > 0
    assert out["surplus"].any()
