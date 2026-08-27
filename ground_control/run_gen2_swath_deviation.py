"""Driver: gen2's per-swath vertical deviation over elbaext, and its sensitivity.

Run (conda env `lidar-icp`, which carries PDAL and laspy):

    ~/anaconda3/envs/lidar-icp/bin/python ground_control/run_gen2_swath_deviation.py \
        --copc data/after/elbaext_3dep_fulldensity.copc.laz \
        --pdal ~/anaconda3/envs/lidar-icp/bin/pdal \
        --res-m 2.0 --block-m 50.0 --half-width-m 200.0 \
        --ladder-half-widths-m 100 200 400 800 \
        --place-eastings 575900 576700 577500 578300 579100 579900 \
        --centre-easting 577825 --elba-easting 578762.8

There are no defaults for the window size, resolution, block size, ladder or places.  A
default would hide the answer to the question this script exists to ask, which is exactly
how much those choices move the number.  ``--nadir-*`` is likewise required: the flight
tracks are a measurement, and the script will not guess them.
"""

from __future__ import annotations

import argparse
import gc
import sys

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "ground_control")
import gen2_swath_deviation as G  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--copc", required=True)
    p.add_argument("--pdal", required=True)
    p.add_argument("--crs", required=True)
    p.add_argument("--res-m", type=float, required=True)
    p.add_argument("--block-m", type=float, required=True)
    p.add_argument("--half-width-m", type=float, required=True)
    p.add_argument("--ladder-half-widths-m", type=float, nargs="+", required=True)
    p.add_argument("--place-eastings", type=float, nargs="+", required=True)
    p.add_argument("--centre-easting", type=float, required=True)
    p.add_argument("--elba-easting", type=float, required=True)
    p.add_argument("--nadir-lines", type=int, nargs="+", required=True,
                   help="flight-line point_source_id values, in track order")
    p.add_argument("--nadir-northings", type=float, nargs="+", required=True,
                   help="the northing of each line's nadir track, same order")
    p.add_argument("--ladder-pair", type=int, nargs=2, required=True)
    p.add_argument("--keep-classes", type=int, nargs="+", required=True)
    return p.parse_args(argv)


def load(a, easting, northing, half_width_m):
    path = G.crop_window(a.copc, easting=easting, northing=northing,
                         half_width_m=half_width_m,
                         keep_classes=tuple(a.keep_classes), pdal_bin=a.pdal)
    return G.window_pointcloud(path, crs=a.crs)


def measure(a, pc, ref, src):
    grids = G.overlap_grids(pc, ref, src, res_m=a.res_m, exclude=())
    if grids is None:
        return None, None
    dh, dt, iy, ix, nx = grids
    bid = G.block_ids(iy, ix, nx, res_m=a.res_m, block_m=a.block_m)
    stats = G.pair_stats(pc, ref, src, res_m=a.res_m, block_m=a.block_m, exclude=())
    fit = G.across_track_ols(dh, dt, bid)
    return stats, fit


def main(argv=None):
    a = parse_args(argv)
    if len(a.nadir_lines) != len(a.nadir_northings):
        raise SystemExit("--nadir-lines and --nadir-northings must have equal length")
    nadir = dict(zip(a.nadir_lines, a.nadir_northings))

    print("PARAMETERS (source of each stated; none is a silent default)")
    print(f"  copc            = {a.copc}")
    print(f"  res_m           = {a.res_m}   [INHERITED: coreg.coregister_swaths default;"
          " res_m in analysis/LOCAL_TIE_CHAINING.md]")
    print(f"  block_m         = {a.block_m}  [INHERITED: 50 m cluster blocks, the repo"
          " default named in analysis/SWATH_ACROSS_TRACK_TEST.md]")
    print(f"  keep_classes    = {tuple(a.keep_classes)}     [FORCED: gen2 carries only"
          " classes 1 and 2; class 2 is gen2's vendor ground]")
    print(f"  exclude         = ()      [FORCED: nothing left to exclude after the"
          " class-2 read; coreg's (5,6,9) would remove NOTHING from gen2]")
    print(f"  half_width_m    = {a.half_width_m}  [CALLER]")
    print(f"  ladder          = {a.ladder_half_widths_m}  [CALLER]")
    print(f"  place eastings  = {a.place_eastings}  [CALLER]")
    print(f"  nadir tracks    = {nadir}  [CALLER; measured, not guessed]")
    print()

    ks = list(a.nadir_lines)
    seams = {(x, y): (nadir[x] + nadir[y]) / 2.0 for x, y in zip(ks, ks[1:])}

    hdr = ("%-11s %8s %6s %10s %9s %10s %9s | %10s %8s %10s %8s"
           % ("pair", "cells", "blks", "mean_mm", "SE_blk", "median_mm", "NMAD_mm",
              "OLS_k_mm", "SE_k", "OLS_c", "SE_c"))

    print("=== (a) EVERY ADJACENT PAIR, each at its own seam, easting %.0f, hw %.0f m ==="
          % (a.centre_easting, a.half_width_m))
    print(hdr)
    for (x, y), seam in seams.items():
        pc, _ = load(a, a.centre_easting, seam, a.half_width_m)
        s, f = measure(a, pc, x, y)
        if s is None:
            print("%-11s  (no overlap)" % ("%d-%d" % (x, y)))
        else:
            print("%-11s %8d %6d %+10.2f %9.2f %+10.2f %9.1f | %+10.2f %8.2f %+10.1f %8.1f"
                  % ("%d-%d" % (x, y), s.n_cells, s.n_blocks, s.mean_mm, s.se_block_mm,
                     s.median_mm, s.nmad_mm, f.k_mm, f.se_k_mm,
                     f.c_mm_per_tan, f.se_c_mm_per_tan))
        del pc; gc.collect()

    pair = tuple(a.ladder_pair)
    seam = seams[pair]
    print()
    print("=== (b) PLACE: pair %d-%d along its seam (N %.0f), hw %.0f m ==="
          % (pair[0], pair[1], seam, a.half_width_m))
    print(hdr.replace("pair", "east", 1))
    means = []
    for e in a.place_eastings:
        pc, _ = load(a, e, seam, a.half_width_m)
        s, f = measure(a, pc, *pair)
        means.append(s.mean_mm)
        print("%-11.0f %8d %6d %+10.2f %9.2f %+10.2f %9.1f | %+10.2f %8.2f %+10.1f %8.1f"
              % (e, s.n_cells, s.n_blocks, s.mean_mm, s.se_block_mm, s.median_mm,
                 s.nmad_mm, f.k_mm, f.se_k_mm, f.c_mm_per_tan, f.se_c_mm_per_tan))
        del pc; gc.collect()
    m = np.array(means)
    print("   --> over %d places: %+.2f to %+.2f mm, spread %.2f, sd %.2f"
          % (m.size, m.min(), m.max(), np.ptp(m), m.std(ddof=1)))

    print()
    print("=== (c) WINDOW SIZE: pair %d-%d at easting %.1f ===" % (*pair, a.elba_easting))
    print(hdr.replace("pair", "hw_m", 1))
    wm = []
    for hw in a.ladder_half_widths_m:
        pc, _ = load(a, a.elba_easting, seam, hw)
        s, f = measure(a, pc, *pair)
        wm.append(s.mean_mm)
        print("%-11.0f %8d %6d %+10.2f %9.2f %+10.2f %9.1f | %+10.2f %8.2f %+10.1f %8.1f"
              % (hw, s.n_cells, s.n_blocks, s.mean_mm, s.se_block_mm, s.median_mm,
                 s.nmad_mm, f.k_mm, f.se_k_mm, f.c_mm_per_tan, f.se_c_mm_per_tan))
        del pc; gc.collect()
    w = np.array(wm)
    print("   --> over %d window sizes: %+.2f to %+.2f mm, spread %.2f, sd %.2f"
          % (w.size, w.min(), w.max(), np.ptp(w), w.std(ddof=1)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
