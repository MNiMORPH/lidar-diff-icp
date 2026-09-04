#!/usr/bin/env python3
"""Per-cell class-2 ground-surface spread (the ground Gaussian's SD), for a tile.

The covariate the correction is actually indexed by -- NOT a cover layer, no windows, no
external product. Written as its own producer because it was previously computed inside
whichever consumer needed it and thrown away, so a second consumer had to recompute it.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/class2_spread_grid.py \
        --tile data/derived/whitewater --gen2 data/after/3dep_4358_fulltile.laz
"""
import argparse, numpy as np
from lidar_diff_icp import groundq

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True)
ap.add_argument("--gen2", required=True)
ap.add_argument("--zlo", type=float, default=-1.0)
ap.add_argument("--zhi", type=float, default=2.0)
ap.add_argument("--dz", type=float, default=0.02)
ap.add_argument("--min-count", type=int, default=20)
ap.add_argument("--chunk", type=int, default=3_000_000)
A = ap.parse_args()

surf = groundq.reference_surface(A.tile)
H, n_in = groundq.column_histogram(A.gen2, surf, zlo=A.zlo, zhi=A.zhi, dz=A.dz,
                                   chunk=A.chunk)
sd, cnt = groundq.spread_from_histogram(H, A.zlo, A.dz, min_count=A.min_count)
np.save(f"{A.tile}/class2_sd_mm.npy", sd.reshape(surf["ny"], surf["nx"]))
np.save(f"{A.tile}/class2_n.npy", cnt.reshape(surf["ny"], surf["nx"]))
ok = np.isfinite(sd) & (sd > 0)
print(f"{n_in:,} class-2 returns in {A.zlo:+.2f}..{A.zhi:+.2f} m")
print(f"class-2 spread on {int(ok.sum()):,} cells (>= {A.min_count} returns): "
      f"median {np.nanmedian(sd):.1f} mm  p10 {np.nanpercentile(sd,10):.1f}  "
      f"p90 {np.nanpercentile(sd,90):.1f}")
print(f"wrote {A.tile}/class2_sd_mm.npy  and  {A.tile}/class2_n.npy")
