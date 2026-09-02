#!/usr/bin/env python3
"""Per-cell ``lowveg`` for a tile: the near-ground vegetation COMPOSITION.

    lowveg = (returns with lo < h <= hi) / (returns in the whole near-ground window)

This is the metric of ``analysis/CONTROL_LOWVEG_OFFSET.md``, which fitted it against
SURVEYED CONTROL and found it outperforms canopy cover -- "It is a composition, not a
canopy density, which is why it outperforms canopy cover here." Nothing here is invented:
the band edges, the population and the bin-centre test are all read off that document and
``analysis/control_lowveg_offset.py``, so a coefficient fitted at the control marks and a
value computed here refer to the same quantity.

Three things make the tile cube and the control-mark boxes commensurate, and all three
were checked rather than assumed:

  window      ``nearground_cells_sn.npz`` stores -1.00..+2.00 m at 0.02 m; the control
              boxes store 150 bins over -1.00..+2.00 m. Same window, same bin width.
  population  ``H2`` is ALL gen2 returns -- nearground_cells.py: "no classification, no
              ground/vegetation decision, no height threshold beyond the window" -- which
              is the ``ng_all`` array ``control_lowveg_offset.lowveg()`` reads. The
              class-2 split lives in nearground_gen2_class_split.npz and is NOT used here.
  bin-centre  the band test is on bin CENTRES, so with 20 mm bins the (0.14, 0.16] bin is
              counted whole and the metric is really "fraction above 0.14 m". The document
              flags this explicitly and states that the -290 mm/unit coefficient is correct
              FOR THE BINNED METRIC and must not be applied to an exactly-computed one.
              Reproducing the binning is therefore required, not incidental.

COVERAGE. The cube holds only the cells nearground_cells.py selected, so lowveg is NaN
elsewhere. That is enough to FIT q2(lowveg), which is evaluated on those same cells. It is
NOT enough to APPLY a correction across the whole DoD -- that needs a grid-wide near-ground
pass over the cloud, which this script deliberately does not fake by interpolating.

    ./lidar-icp/bin/python scripts/make_lowveg.py --tile elba
    ./lidar-icp/bin/python scripts/make_lowveg.py --tile elba --check
"""
import argparse
import os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True, help="tile name under data/derived/, or a path")
ap.add_argument("--cube", default="nearground_cells_sn.npz")
# Band edges are CITED, not chosen here. CONTROL_LOWVEG_OFFSET.md, "lowveg, EXACT
# DEFINITION -- it must travel with any coefficient fitted from it": (0.15, 2.00] m.
ap.add_argument("--lo", type=float, default=0.15, help="band lower edge, m (doc: 0.15)")
ap.add_argument("--hi", type=float, default=2.00, help="band upper edge, m (doc: 2.00)")
ap.add_argument("--out", default="lowveg.npy")
ap.add_argument("--check", action="store_true",
                help="compare against the existing file and write nothing")
A = ap.parse_args()

D = A.tile if os.path.sep in A.tile else os.path.join("data", "derived", A.tile)
cube_path = os.path.join(D, A.cube)
if not os.path.exists(cube_path):
    raise SystemExit(f"{cube_path} is missing. It is the product of "
                     f"analysis/ridgelines/nearground_cells.py --tile {D} "
                     f"--out {A.cube}; lowveg is a slice of that cube and cannot be "
                     f"computed without it.")

z = np.load(cube_path)
cells, edges, H2 = z["cells"], z["edges"], z["H2"].astype(float)
mid = 0.5 * (edges[:-1] + edges[1:])
band = (mid > A.lo) & (mid <= A.hi)              # bin-CENTRE test, as at the control marks

tot = H2.sum(1)
num = H2[:, band].sum(1)
with np.errstate(invalid="ignore", divide="ignore"):
    lv_cells = np.where(tot > 0, num / tot, np.nan)

ny, nx = np.load(os.path.join(D, "z_after.npy")).shape
out = np.full(ny * nx, np.nan)
out[cells] = lv_cells
out = out.reshape(ny, nx)

fin = np.isfinite(out)
print(f"lowveg = (returns with {A.lo:g} < h <= {A.hi:g} m) / (returns in "
      f"{edges[0]:+.2f}..{edges[-1]:+.2f} m), bin-centre test on {edges[1]-edges[0]:.3f} m bins")
print(f"  population: H2 = ALL gen2 returns in the window (not class-2)")
print(f"  band covers {int(band.sum())} of {band.size} bins, "
      f"centres {mid[band][0]:+.3f}..{mid[band][-1]:+.3f} m")
print(f"  cells with a value: {int(fin.sum()):,} of {out.size:,} grid cells "
      f"({int((tot == 0).sum()):,} cube cells had no returns)")
if fin.any():
    q = np.nanpercentile(out, [10, 50, 90])
    print(f"  lowveg  p10 {q[0]:.4f}  median {q[1]:.4f}  p90 {q[2]:.4f}  max {np.nanmax(out):.4f}")

dst = os.path.join(D, A.out)
if A.check:
    if not os.path.exists(dst):
        raise SystemExit(f"--check: {dst} does not exist")
    old = np.load(dst)
    same = np.array_equal(np.nan_to_num(old, nan=-9e9), np.nan_to_num(out, nan=-9e9))
    print(f"  --check: {'IDENTICAL' if same else 'DIFFERS'} from {dst}; nothing written")
    raise SystemExit(0 if same else 1)

np.save(dst, out)
print(f"wrote {dst}")
