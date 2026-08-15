#!/usr/bin/env python3
"""Nuth & Kaeaeb co-registration of every overlapping swath pair in a tile.

Reports the horizontal (dx, dy) and vertical (dz) shift that aligns each source
swath onto a reference swath, with the formal fit sigma and the robust scatter
(NMAD) before and after. Note: the formal sigma assumes independent cells and
therefore *underestimates* the true uncertainty (DEM errors are spatially
correlated); treat it as a lower bound.

Example
-------
    python scripts/swath_coregister.py data/before/4342-29-64.laz
"""
import argparse
from itertools import combinations
from pathlib import Path

from lidar_diff_icp import io, coreg


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tile", help="path to a LAZ/LAS tile")
    p.add_argument("--res", type=float, default=2.0, help="grid size (m)")
    args = p.parse_args()

    pc = io.read_tile(args.tile)
    swaths = pc.swaths.tolist()
    print(f"{Path(args.tile).name}: swaths {swaths}")
    print(f"{'ref->src':>10} {'dx_m':>8} {'dy_m':>8} {'dz_m':>8} "
          f"{'nmad0':>7} {'nmad1':>7} {'conv':>5}")
    for a, b in combinations(swaths, 2):
        try:
            c = coreg.coregister_swaths(pc, a, b, res=args.res)
        except ValueError:
            continue  # non-overlapping
        print(f"{a:>4}->{b:<4} {c.dx:>8.3f} {c.dy:>8.3f} {c.dz:>8.3f} "
              f"{c.nmad_before:>7.3f} {c.nmad_after:>7.3f} {str(c.converged):>5}")


if __name__ == "__main__":
    main()
