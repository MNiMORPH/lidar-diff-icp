"""Tests for the deterministic geoid-difference datum (references.geoid_difference).

The earlier reference_plane machinery (flat_hard_cells / datum_plane / datum_offset /
eval_datum_correction) was removed in the geoid-only datum refactor; this tests its
replacement — the GEOID03->GEOID18 model difference auto-computed from the PROJ grids.
"""
import pytest
from lidar_diff_icp import references

ELBA_BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6)   # EPSG:26915


def test_geoid_difference_reproduces_pilot():
    """N_gen1 - N_gen2 (GEOID03 - GEOID18) at Elba is independently ~+67 mm, with a
    small (sub-mm/km) planar tilt. Skips if the PROJ geoid grids can't be reached."""
    # GEOID03 is STATED, not defaulted: the default was removed 2026-09-07 after it was
    # silently applied to three sites whose gen1 is on GEOID09.
    from lidar_diff_icp import acquisitions
    grid = acquisitions.for_project("lidar_semn2008").geoid_grid
    try:
        a, b, c = references.geoid_difference(ELBA_BOUNDS, 26915, before_geoid=grid)
    except TypeError:
        raise                                    # a signature error is a BUG, never a skip
    except Exception as e:                       # noqa: BLE001 -- PROJ/network/grid availability
        pytest.skip(f"PROJ geoid grids unavailable (needs proj.db + GEOID03/18 grids): {e}")
    assert 0.060 < a < 0.075, f"const {a} m not ~+67 mm"
    assert abs(b) < 0.003 and abs(c) < 0.003, f"tilt ({b},{c}) m/km not small"
