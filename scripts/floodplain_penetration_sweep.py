#!/usr/bin/env python3
"""Does gen1's ground keep dropping as you coarsen? The tussock-penetration test.

Andy, 2026-09-11: "Much of the floodplain 'erosion' is from gen1 lidar not penetrating
tussocky vegetation. Is it possible to correct for that by coarsening and looking for any
ground returns?"

THE MECHANISM UNDER TEST. In dense tussock a 5 m cell may contain NO pulse that reached
true ground. CSF's cloth then settles on the tussock tops and EVERY percentile of that
cell is biased high -- gen1 reads high, gen2 - gen1 reads negative, and apparent erosion
appears with no sediment moved. Widening the neighbourhood cannot lower a well-penetrated
cell, but it CAN find a rare penetrating shot where there was none. So this is
coarsen-to-find-penetration, not "take a lower quantile".

THE SIGNATURE, and it is the point of sweeping rather than measuring once:

    flat in block size          penetration is adequate; there is nothing lower to find
    falls, then plateaus        tussock-blocked; the PLATEAU SIZE is the penetration
                                length scale and the DROP is the bias

CONTROL: the same sweep on flat ground OUTSIDE the valley. If the drop is a tussock
effect it should be far weaker there; if both fall alike, we are watching the coarsening
artefact (a low percentile over a wider window drifts downhill) rather than penetration.

SCOPED BY SLOPE, AND THE CUT IS SWEPT. The tile's floodplain_mask is an ELEVATION cut --
everything below the 230 m valley top -- so at elbaext its median slope is 12.2 deg
against 5.9 deg elsewhere: it contains the valley WALLS, not just the floor. Coarsening on
12 deg ground reintroduces exactly the downhill low-percentile bias this is trying to
avoid. Only 7.3% of that mask is <= 1 deg and 18.9% is <= 3 deg. Which cut defines "the
flat floodplain" is a judgement, so every cut is reported rather than one chosen.

THE QUANTITY is `d_mm_corr` from the beam table: each gen1 CSF ground return's SLOPE-NORMAL
offset to the gen2 surface, in mm, after registration. Already detrended, so no plane fit
is needed and terrain slope does not enter the percentile. Its ABSOLUTE level carries
whatever the registration did; only its CHANGE WITH BLOCK SIZE is read here, which is
level-invariant.

    ./lidar-icp/bin/python scripts/floodplain_penetration_sweep.py --tile data/derived/elbaext_independent
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="data/derived/elbaext_independent")
ap.add_argument("--masks", default="data/derived/elbaext",
                help="where floodplain_mask.npy and slope.npy live")
ap.add_argument("--slope-cuts", type=float, nargs="+", default=[1.0, 2.0, 3.0, 5.0],
                help="SWEPT, not chosen: which cut defines 'flat floodplain' is a judgement")
ap.add_argument("--blocks", type=int, nargs="+", default=[1, 2, 4, 8, 16],
                help="block size in CELLS; at 5 m that is 5..80 m")
ap.add_argument("--min-returns", type=int, default=20,
                help="a block with fewer returns than this has no usable percentile and is "
                     "reported as dropped, never quietly averaged in")
ap.add_argument("--out", default="")
A = ap.parse_args()

fp = np.load(os.path.join(A.masks, "floodplain_mask.npy")).astype(bool)
slope = np.load(os.path.join(A.masks, "slope.npy"))
ny, nx = fp.shape
res = float(json.load(open(os.path.join(A.tile, "corrections.json")))["res_m"])

t = pq.read_table(os.path.join(A.tile, "beam_offset_table.parquet"),
                  columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]
d = t["d_mm_corr"].to_numpy()[g].astype(float)
iy, ix = np.divmod(cell, nx)
keep = (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx) & np.isfinite(d)
iy, ix, d = iy[keep], ix[keep], d[keep]
print(f"  {len(d):,} gen1 ground returns in grid, {ny}x{nx} at {res:g} m\n")


def sweep(sel_mask, label):
    """Per-block percentiles of d_mm, median over blocks, against block size."""
    sel = sel_mask[iy, ix]
    dv, yy, xx = d[sel], iy[sel], ix[sel]
    print(f"  {label}: {int(sel_mask.sum()):,} cells, {len(dv):,} returns")
    print(f"    {'block':>7}{'blocks':>9}{'dropped':>9}{'ret/blk':>9}"
          f"{'p05':>9}{'p25':>9}{'p50':>9}{'d(p05)':>9}")
    rows, first = [], None
    for b in A.blocks:
        key = (yy // b).astype(np.int64) * 100000 + (xx // b)
        s = pd.Series(dv).groupby(key)
        n = s.size()
        big = n >= A.min_returns
        if big.sum() < 5:
            print(f"    {b*res:>6.0f}m   (too few blocks with >= {A.min_returns} returns)")
            continue
        p05 = s.quantile(0.05)[big]; p25 = s.quantile(0.25)[big]; p50 = s.quantile(0.50)[big]
        m05, m50 = float(p05.median()), float(p50.median())
        if first is None:
            first = m05
        print(f"    {b*res:>6.0f}m{int(big.sum()):>9,}{int((~big).sum()):>9,}"
              f"{float(n[big].median()):>9.0f}{m05:>9.1f}{float(p25.median()):>9.1f}"
              f"{m50:>9.1f}{m05-first:>9.1f}")
        rows.append(dict(block_m=b*res, n_blocks=int(big.sum()), dropped=int((~big).sum()),
                         returns_per_block=float(n[big].median()), p05=m05,
                         p25=float(p25.median()), p50=m50, drop_from_first=m05-first))
    print()
    return rows


out = {}
for cut in A.slope_cuts:
    flat = slope <= cut
    out[f"floodplain_slope_le_{cut:g}"] = sweep(fp & flat, f"VALLEY  & slope <= {cut:g} deg")
    out[f"upland_slope_le_{cut:g}"] = sweep(~fp & flat, f"UPLAND  & slope <= {cut:g} deg  (CONTROL)")

print("  READ THE d(p05) COLUMN: how far the 5th percentile falls as blocks grow.")
print("  A real penetration deficit falls in the VALLEY and not in the UPLAND control.")
print("  Both falling together is the coarsening artefact, not tussocks.")
if A.out:
    json.dump(out, open(A.out, "w"), indent=1)
    print(f"\n  wrote {A.out}")
