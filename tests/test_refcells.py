"""The reference-cell population must apply each cut, and must NOT be offset-driven."""
import numpy as np
import pytest

from lidar_diff_icp.refcells import reference_cells


@pytest.fixture
def tile(tmp_path):
    """A 4x4 tile: all cells stable except one of each failure mode."""
    n = 16
    np.save(tmp_path / "slope.npy", np.full((4, 4), 5.0))
    np.save(tmp_path / "ridge_mask.npy", np.ones((4, 4), bool))
    np.save(tmp_path / "curv_laplacian.npy", np.zeros((4, 4)))
    np.save(tmp_path / "dod.npy", np.zeros((4, 4)))
    np.save(tmp_path / "floodplain_mask.npy", np.zeros((4, 4), bool))
    for g in ("gen1", "gen2"):
        np.savez(tmp_path / f"{g}_canopy_frac.npz", n_bldg=np.zeros((4, 4), int),
                 frac=np.full((4, 4), 0.5))
    return tmp_path


def _set(path, key, flat_i, value):
    a = np.load(path) if str(path).endswith(".npy") else None
    if a is not None:
        a.ravel()[flat_i] = value
        np.save(path, a)
    else:
        z = dict(np.load(path))
        z[key].ravel()[flat_i] = value
        np.savez(path, **z)


def test_all_stable_keeps_everything(tile):
    m, rep = reference_cells(tile)
    assert m.all() and rep["kept"] == 16


@pytest.mark.parametrize("fname,key,val,label", [
    ("slope.npy", None, 30.0, "slope >= 12 deg"),
    ("curv_laplacian.npy", None, 0.5, "|curv| > 0.015"),
    ("ridge_mask.npy", None, False, "not a divide cell"),
    ("gen1_canopy_frac.npz", "n_bldg", 3, "building returns"),
    ("dod.npy", None, 1.0, "|DoD| > 500 mm"),
])
def test_each_criterion_removes_its_cell(tile, fname, key, val, label):
    _set(tile / fname, key, 5, val)
    m, rep = reference_cells(tile)
    assert not m[5] and m.sum() == 15
    assert rep[label] == 1


def test_clearcut_is_one_sided_so_leaf_out_is_kept(tile):
    """Canopy GAINED between epochs is deciduous leaf-out, not change: keep it.
    Canopy LOST is a clear-cut: drop it."""
    _set(tile / "gen1_canopy_frac.npz", "frac", 5, 0.05)   # gen1 low, gen2 0.5: leaf-out
    _set(tile / "gen1_canopy_frac.npz", "frac", 6, 0.95)   # gen1 high, gen2 0.5: cut
    m, rep = reference_cells(tile)
    assert m[5], "a leaf-out gain must not be treated as change"
    assert not m[6]
    assert rep["clear-cut (frac drop > 0.3)"] == 1


def test_gross_change_cannot_remove_a_vegetation_scale_offset(tile):
    """The default gross-change guard sits far above the effect being calibrated."""
    _set(tile / "dod.npy", None, 5, -0.130)      # the largest dense-canopy offset seen
    m, _ = reference_cells(tile)
    assert m[5], "a 130 mm canopy-scale offset must survive the stability cut"


def test_cells_subset_and_report_order(tile):
    _set(tile / "slope.npy", None, 5, 30.0)
    m, rep = reference_cells(tile, cells=np.array([4, 5, 6]))
    assert m.tolist() == [True, False, True]
    assert rep["start"] == 3 and rep["kept"] == 2


def test_a_missing_floodplain_mask_refuses_rather_than_skipping_the_cut(tile):
    """The defect this prevents: the cut was skipped in SILENCE when the file was absent,
    so two tiles' reference populations differed by 39,038 cells with nothing to notice,
    and every comparison between them was invalid without saying so."""
    (tile / "floodplain_mask.npy").unlink()
    with pytest.raises(FileNotFoundError, match="not comparable"):
        reference_cells(tile, use_floodplain_mask=True)


def test_the_valley_cut_is_by_ELEVATION_not_the_tpi_mask(tile):
    """Andy, 2026-09-04: cut the floodplain on ELEVATION. The TPI mask is a
    topographic-position heuristic whose extent depends on the window width, and it keeps
    flat terrace ground sitting at valley level. The default must not silently apply it."""
    m, rep = reference_cells(tile)
    assert not any("floodplain mask" in k for k in rep), rep
    # The elevation cut fires only where the histogram is bimodal -- a tile with no valley
    # has nothing to remove, and this fixture is unimodal. On the real tiles it does fire:
    # elba cuts below 286.1 m and whitewater below 282.0 m.
    m2, rep2 = reference_cells(tile, use_floodplain_mask=True)
    assert any("floodplain mask" in k for k in rep2), rep2
    assert m2.sum() <= m.sum()


def test_working_without_the_mask_is_allowed_but_must_be_stated(tile):
    (tile / "floodplain_mask.npy").unlink()
    m, rep = reference_cells(tile, use_floodplain_mask=False)
    assert m.sum() > 0
    assert not any("floodplain" in k for k in rep)
