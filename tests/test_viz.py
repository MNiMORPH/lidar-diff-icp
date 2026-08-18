"""GDAL shaded relief must be geographically NW-illuminated (not flipped)."""
import numpy as np
import pytest

from lidar_diff_icp.viz import hillshade, _find_gdaldem


def test_hillshade_is_geographically_NW():
    try:
        _find_gdaldem()
    except Exception:
        pytest.skip("gdaldem not available")
    n = 81
    y, x = np.mgrid[0:n, 0:n]
    z = 100 - 0.01 * ((x - 40) ** 2 + (y - 40) ** 2)     # central peak, origin='lower'
    hs = hillshade(z, 5.0, 500000.0, 4900000.0)
    nw = np.nanmean(hs[60:75, 5:20])                      # north (high row), west (low col)
    se = np.nanmean(hs[5:20, 60:75])
    assert nw > se + 0.2                                  # NW flank clearly lit (315 az)
    assert np.nanmin(hs) >= 0.0 and np.nanmax(hs) <= 1.0
