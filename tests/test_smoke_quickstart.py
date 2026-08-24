"""Smoke test: the quick-start path (scripts/gridded_ground_dod.py -> difference_dem) is
compatible with the geoid-only pipeline.

The datum refactor removed allow_parabola / ref_polys / save_ref_cells / datum_tilt and
made the geoid the (auto-computed) datum. This catches a caller passing a now-removed
kwarg without running the heavy end-to-end pipeline.
"""
import inspect
from lidar_diff_icp.pipeline import difference_dem


def test_difference_dem_signature_is_geoid_only():
    p = set(inspect.signature(difference_dem).parameters)
    assert "geoid_datum" in p
    for removed in ("allow_parabola", "ref_polys", "save_ref_cells", "datum_tilt"):
        assert removed not in p, f"{removed} should have been removed"


def test_gridded_ground_dod_kwargs_are_valid():
    """Every kwarg scripts/gridded_ground_dod.py passes must be a current param."""
    used = {"res", "ground_q", "ground", "ground_source", "after_ground",
            "csf_pdal", "stream", "correction_surface", "along_track_drift"}
    p = set(inspect.signature(difference_dem).parameters)
    assert used <= p, f"stale kwargs: {used - p}"
