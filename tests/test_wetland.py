"""Return-statistic water/wetland flag: separate dropout (open water) and
attenuated-ground flat lows (wetland) from sloped dry land."""
import numpy as np

from lidar_diff_icp.wetland import wetland_flag, DRY, WETLAND, OPEN_WATER


def test_flag_separates_water_wetland_dry():
    # The flag assumes wetland is a MINORITY of the tile (as on real tiles): the
    # ground-fraction threshold is a low percentile, so a wetland-dominated tile
    # would need intensity or a raised percentile. Here water+wetland are minority.
    H = W = 60; res = 5.0
    Y, X = np.mgrid[0:H, 0:W]
    dem = np.where(Y < 40, 100.0 + 0.4 * X, 90.0)      # sloped upland; flat LOW basin below
    total = np.full((H, W), 250.0)
    ground = np.full((H, W), 200.0)                    # dry baseline: dense ground returns
    ground[45:55, 5:20] = 10.0                         # wetland strip: attenuated, flat low
    total[45:55, 40:55] = 3.0; ground[45:55, 40:55] = 0.0   # open-water strip: dropout

    flag = wetland_flag(ground, total, dem, res)
    assert (flag[:40, :] == DRY).mean() > 0.9          # sloped upland -> dry
    assert (flag[45:55, 5:20] == WETLAND).mean() > 0.6  # attenuated flat low -> wetland
    assert (flag[45:55, 40:55] == OPEN_WATER).mean() > 0.9  # dropout -> open water
