#!/usr/bin/env python3
"""Inter-swath consistency for one tile: the 2008 self-calibration check.

Computes the density-robust bare-earth difference between every pair of
overlapping flight lines in a tile, prints the robust vertical offset, scatter,
and tilt, and (optionally) writes a difference-map figure per pair.

Example
-------
    python scripts/swath_consistency.py data/before/4342-29-64.laz --figdir figures
"""
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from lidar_diff_icp import io, swathdiff


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tile", help="path to a LAZ/LAS tile")
    p.add_argument("--res", type=float, default=2.0, help="grid size (m)")
    p.add_argument("--figdir", help="if set, write a diff map per overlapping pair")
    args = p.parse_args()

    pc = io.read_tile(args.tile)
    swaths = pc.swaths.tolist()
    print(f"{Path(args.tile).name}: {len(pc):,} points, swaths {swaths}")
    print(f"{'pair':>11} {'ncells':>8} {'offset_m':>9} {'rob_std_m':>10} {'tilt_mm/m':>10}")

    results = []
    for a, b in combinations(swaths, 2):
        try:
            d = swathdiff.swath_difference(pc, a, b, res=args.res)
        except ValueError:
            continue  # non-overlapping pair
        results.append(d)
        print(f"{a:>5}-{b:<5} {d.n_cells:>8,} {d.median_offset:>9.3f} "
              f"{d.robust_std:>10.3f} {d.tilt:>10.3f}")

    if args.figdir:
        _plot(results, args.figdir, Path(args.tile).stem)


def _plot(results, figdir, stem) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Path(figdir).mkdir(parents=True, exist_ok=True)
    for d in results:
        fig, ax = plt.subplots(figsize=(4, 6))
        im = ax.imshow(d.diff, origin="lower", extent=d.extent,
                       vmin=-0.3, vmax=0.3, cmap="RdBu_r", aspect="equal")
        ax.set_title(f"swath {d.swath_a} - {d.swath_b}\n"
                     f"offset {d.median_offset:+.3f} m, tilt {d.tilt:.2f} mm/m")
        ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
        fig.colorbar(im, ax=ax, label="Δz (m)", shrink=0.7)
        out = Path(figdir) / f"{stem}_swath_{d.swath_a}-{d.swath_b}.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
