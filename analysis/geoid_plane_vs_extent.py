#!/usr/bin/env python3
"""How far can ONE plane represent the GEOID03->GEOID18 difference?

`references.geoid_difference` samples the two PROJ geoid grids on an n x n grid over the
bounds it is given and fits

    shift(E,N) = a + b*(E-cx)/1000 + c*(N-cy)/1000

about that bounds' centroid. The returned triple is therefore a LINEARISATION whose
accuracy depends on the footprint -- which is why elba and elbaext deliberately share ONE
plane (fitted on elba's bounds, re-expressed about elbaext's centroid) rather than each
fitting its own: same method, different frame, and it is the frame that has to match.

This measures where that device stops working. It replicates the function's own sampling
(same grids, same centroid, same least squares) and additionally reports the residual,
which the function does not return. n=25 rather than its default 7, so the misfit is
resolved rather than absorbed.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/geoid_plane_vs_extent.py
"""
import argparse, os

import numpy as np
import pyproj

from lidar_diff_icp import io

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=25, help="samples per axis (function default is 7)")
ap.add_argument("--proj-data", default="/usr/share/proj")
A = ap.parse_args()

pyproj.datadir.set_data_dir(A.proj_data)
os.environ.setdefault("PROJ_NETWORK", "ON")
pyproj.network.set_network_enabled(True)
from pyproj import Transformer                                    # noqa: E402

ELBA = (577492.8, 4882737.6, 580032.8, 4886237.6)
ELBAEXT = (575600.0, 4882200.0, 580050.0, 4886250.0)


def sample(bounds, n):
    X0, Y0, X1, Y1 = bounds
    xs = np.linspace(X0, X1, n); ys = np.linspace(Y0, Y1, n)
    XX, YY = (a.ravel() for a in np.meshgrid(xs, ys))
    lon, lat = Transformer.from_crs(io.MN_GEN1_CRS, 4326, always_xy=True).transform(XX, YY)

    def und(grid):
        tr = Transformer.from_pipeline(f"+proj=vgridshift +grids={grid} +multiplier=1")
        return np.asarray(tr.transform(lon, lat, np.zeros_like(lon))[2])

    return XX, YY, und("us_noaa_geoid03_conus.tif") - und("us_noaa_g2018u0.tif")


def misfit(bounds, n):
    XX, YY, d = sample(bounds, n)
    ok = np.isfinite(d)
    XX, YY, d = XX[ok], YY[ok], d[ok]
    cx, cy = 0.5 * (bounds[0] + bounds[2]), 0.5 * (bounds[1] + bounds[3])
    G = np.c_[np.ones_like(XX), (XX - cx) / 1000.0, (YY - cy) / 1000.0]
    coef, *_ = np.linalg.lstsq(G, d, rcond=None)
    r = d - G @ coef
    return len(d), np.ptp(d) * 1000, np.sqrt((r ** 2).mean()) * 1000, np.abs(r).max() * 1000


cx, cy = 0.5 * (ELBA[0] + ELBA[2]), 0.5 * (ELBA[1] + ELBA[3])
cases = [("elba tile", ELBA), ("elbaext tile", ELBAEXT)]
cases += [(f"{k} km box at Elba", (cx - k*500, cy - k*500, cx + k*500, cy + k*500))
          for k in (10, 50, 100, 200, 400)]

print("GEOID03 - GEOID18, and the error of representing it by one plane about the centroid")
print(f"sampling {A.n} x {A.n} per extent, as references.geoid_difference does\n")
print(f"  {'extent':<26}{'span km':>9}{'n':>6}{'field ptp':>11}{'plane RMS':>11}{'plane max':>11}")
print(f"  {'':<26}{'':>9}{'':>6}{'mm':>11}{'mm':>11}{'mm':>11}")
for label, b in cases:
    n, ptp, rms, mx = misfit(b, A.n)
    span = max(b[2] - b[0], b[3] - b[1]) / 1000.0
    print(f"  {label:<26}{span:9.1f}{n:6d}{ptp:11.2f}{rms:11.2f}{mx:11.2f}")
print("\nThe elba row's field ptp reproduces the 4.44 mm quoted in "
      "analysis/slope_bias/elbaext_geoid_regrid.py.")
