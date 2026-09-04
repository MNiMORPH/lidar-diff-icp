#!/usr/bin/env python3
"""Per-swath along-track drift correction for early Minnesota (2008) lidar.

The dominant residual after per-swath translation + a smooth cross-epoch tie is a
DETERMINISTIC along-track GNSS-trajectory drift: a smooth vertical undulation
along each flight line (a function of gps_time), roughly constant across-track,
with no roll/scan-angle component. This is the reusable core of a statewide 2008
correction -- the FORM is universal (per-swath f(gps_time)); only the coefficients
differ per tile. The 2008 data retains point_source_id + gps_time, which is all
the model needs; a modern reference (2021 3DEP) supplies stable-ground control.

Fit: on stable ground, the observed change (reference - 2008) is -drift. Per
swath, bin the change by gps_time, take a robust median, smooth it, and evaluate
that curve at every point of the swath. Correcting the 2008 elevations by this
smooth per-swath curve removes the drift while leaving everything shorter than the
along-track smoothing scale untouched (so it does NOT overfit terrain-correlated
noise -- unlike a spatial interpolator on the difference).

    python scripts/along_track_drift.py data/before/4342-29-64.laz \
        data/derived/change_core2008_robust.laz --bounds 577492.8 4882737.6 580035.0 4886238.3
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import pandas as pd
from scipy.ndimage import (gaussian_filter, uniform_filter,
                           distance_transform_edt as edt, gaussian_filter1d)

from lidar_diff_icp import coreg, terrain


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("before_laz"); ap.add_argument("change_laz",
                    help="difference product (reference - 2008) with an m3c2 dim")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--change-thresh", type=float, default=0.15)
    ap.add_argument("--valley-top", dest="valley_top", required=True,
                    help="valley top for the stable mask: an elevation in metres, "
                         "'registry', or 'histogram'. Never chosen for you.")
    ap.add_argument("--tile-dir", dest="tile_dir", default=None,
                    help="tile directory, needed by --valley-top registry/histogram")
    ap.add_argument("--out-drift-laz", default="data/derived/2008_along_track_drift.laz",
                    help="2008 last-return points with a per-point 'drift' dim")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = 5.0; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    def grid(x, y, v, q=None):
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & np.isfinite(v)
        gb = pd.Series(v[ok]).groupby(iy[ok] * nx + ix[ok])
        s = gb.quantile(q) if q is not None else gb.median()
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
        return out.reshape(ny, nx)

    # stable-ground mask + change field, gridded
    h = laspy.read(a.change_laz)
    D = grid(np.asarray(h.x), np.asarray(h.y), np.asarray(h.m3c2))
    Zt = grid(np.asarray(h.x), np.asarray(h.y), np.asarray(h.z))
    # ONE definition, from terrain.py. This was a hand copy of the pipeline's mask and
    # still carried the TPI valley cut after the pipeline moved to elevation.
    stable = terrain.terrain_masks(Zt, res, valley_top_m=a.valley_top,
                                   tile_dir=a.tile_dir)["stable"]

    # 2008 last returns + acquisition
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); gt = np.asarray(f.gps_time)
    rn = np.asarray(f.return_number); nr = np.asarray(f.number_of_returns); be = rn == nr
    x8, y8, z8, ps8, gt = x8[be], y8[be], z8[be], ps8[be], gt[be]
    ix = np.clip(((x8 - X0) / res).astype(int), 0, nx - 1)
    iy = np.clip(((y8 - Y0) / res).astype(int), 0, ny - 1)
    chg = D[iy, ix]
    is_stable = stable[iy, ix] & np.isfinite(chg) & (np.abs(chg) < a.change_thresh)

    drift, curves = coreg.fit_along_track_drift(gt, chg, is_stable, ps8)
    print(f"fitted along-track drift for {len(curves)} swaths; "
          f"drift NMAD {1.4826 * np.median(np.abs(drift - np.median(drift))):.3f} m; "
          f"per-swath amplitude (max-min): "
          + ", ".join(f"{p}:{(c[1].max() - c[1].min()):.3f}m" for p, c in curves.items()), flush=True)

    # write 2008 last returns with drift and corrected z (z + drift)
    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [x8.min(), y8.min(), z8.min()]; hdr.scales = [0.01, 0.01, 0.01]
    hdr.add_extra_dim(laspy.ExtraBytesParams(name="drift", type=np.float32))
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x8, y8, z8 + drift          # corrected elevation
    las.point_source_id = ps8.astype(np.uint16); las.gps_time = gt
    las.drift = drift.astype(np.float32)
    Path(a.out_drift_laz).parent.mkdir(parents=True, exist_ok=True)
    las.write(a.out_drift_laz)
    print(f"wrote {a.out_drift_laz} (z corrected by along-track drift; 'drift' dim retained)", flush=True)


if __name__ == "__main__":
    main()
