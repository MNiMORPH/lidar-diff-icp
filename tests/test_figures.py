"""The two standard figures, rebuildable from the saved products alone.

They used to live inside scripts/run_all_sites.py and take in-memory arrays, so redrawing
one meant re-running the whole site pipeline -- minutes of point-cloud work for a picture
whose every input was already on disk.
"""
import json
import os

import numpy as np
import pytest

from lidar_diff_icp import figures as F


def _tile(tmp_path, ny=20, nx=15, res=5.0, x0=1000.0, y0=2000.0, **arrays):
    d = tmp_path / "tile"
    d.mkdir()
    json.dump({"bounds": [x0, y0, x0 + nx * res, y0 + ny * res], "res_m": res},
              open(d / "corrections.json", "w"))
    rng = np.random.default_rng(0)
    base = {"z_after.npy": rng.normal(300.0, 5.0, (ny, nx)),
            "dod.npy": rng.normal(0.0, 0.05, (ny, nx)),
            "lod.npy": np.full((ny, nx), 0.08),
            "change.npy": rng.random((ny, nx)) > 0.7}
    base.update(arrays)
    for k, v in base.items():
        np.save(d / k, v)
    return str(d)


def test_grid_comes_from_the_tile_not_an_argument(tmp_path):
    d = _tile(tmp_path)
    assert F.grid_of(d) == (1000.0, 2000.0, 5.0, 15, 20)


def test_grid_metadata_is_required(tmp_path):
    d = tmp_path / "bare"
    d.mkdir()
    with pytest.raises(SystemExit, match="no grid metadata"):
        F.grid_of(str(d))


def test_both_figures_are_written_from_products_alone(tmp_path):
    d = _tile(tmp_path)
    fd = str(tmp_path / "figs")
    a = F.dod_lod_figure(d, fd)
    b = F.change_figure(d, fd)
    for p in (a, b):
        assert os.path.exists(p) and os.path.getsize(p) > 1000
    assert a.endswith("tile_dod_lod.png") and b.endswith("tile_change.png")


def test_a_missing_product_refuses_rather_than_drawing_something_else(tmp_path):
    d = _tile(tmp_path)
    os.remove(os.path.join(d, "lod.npy"))
    with pytest.raises(SystemExit, match="lod.npy is missing"):
        F.dod_lod_figure(d, str(tmp_path / "figs"))
    # the other figure does not need it, and still works
    assert os.path.exists(F.change_figure(d, str(tmp_path / "figs")))


def test_the_name_can_be_overridden_for_title_and_filename(tmp_path):
    d = _tile(tmp_path)
    out = F.change_figure(d, str(tmp_path / "figs"), name="carlton")
    assert out.endswith("carlton_change.png")
