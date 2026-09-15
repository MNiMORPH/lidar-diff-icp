#!/usr/bin/env python3
"""Rebuild every floodplain_mask.npy as a SIMPLE ELEVATION CUT, from a stated valley top.

Andy, 2026-09-15: "All the floodplain masks should be simple elevation cuts until we have a
reliable method to produce better ones."

WHY. `refcells.floodplain_by_elevation` already IS a simple cut -- ``isfinite(z) & (z < t)``,
nothing else. The bug is in ``t``. The producer
(``steps/convexity_dod_landcover.py``) defaults ``--valley-top`` to ``"histogram"``, and the
histogram method is WITHDRAWN in the docstring of the very function that consumes it: it
worked at elba only because that tile is genuinely bimodal, and at whitewater the histogram
simply decays so it "returned an arbitrary point in it, cutting 76% of the tile".

Measured at elbaext, the saved mask has median elevation 261.5 m and p90 328.9 m against a
registry valley top of 230.0 m: 214,460 of its 259,178 cells are ABOVE the valley top and
50,788 cells BELOW it are missing. 83% of the "floodplain" is upland. The withdrawal was
recorded in refcells but the producer's default was never changed to match.

THE VALLEY TOP IS NEVER INVENTED HERE. Only tiles with a value in
``refcells.VALLEY_TOP_M`` are rebuilt. Every other tile is REPORTED, with the elevation
summary needed to state one, and left untouched. A wrong valley top is exactly the failure
being fixed, so guessing one to finish the sweep would reproduce it.

Build VARIANTS of a registered tile (elbaext_independent, elbaext_delong) take that tile's
value by longest-prefix match, and the match is printed. They are deliberately NOT added to
the registry, which is keyed on real tiles.

The previous mask is copied to ``floodplain_mask_histogram.npy.bak`` before anything is
written, so a wrong call is recoverable.

    ./lidar-icp/bin/python scripts/rebuild_floodplain_masks.py            # report only
    ./lidar-icp/bin/python scripts/rebuild_floodplain_masks.py --write
"""
import argparse
import glob
import json
import os
import shutil

import numpy as np
from scipy.ndimage import distance_transform_edt

from lidar_diff_icp.refcells import VALLEY_TOP_M, floodplain_by_elevation

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--write", action="store_true", help="write the masks (default: report)")
ap.add_argument("--glob", default="data/derived/*/floodplain_mask.npy")
A = ap.parse_args()


def valley_top_for(tile):
    """Registry value, or the longest registered prefix (a build variant of a real tile)."""
    if tile in VALLEY_TOP_M:
        return float(VALLEY_TOP_M[tile]), tile
    cands = [k for k in VALLEY_TOP_M if tile.startswith(k + "_")]
    if cands:
        k = max(cands, key=len)
        return float(VALLEY_TOP_M[k]), k
    return None, None


done, blocked = [], []
print(f"{'tile':<24}{'valley top':>12}{'source':>22}{'old cells':>11}{'new cells':>11}{'old p50 m':>11}")
for p in sorted(glob.glob(A.glob)):
    D = os.path.dirname(p)
    tile = os.path.basename(D)
    zp = os.path.join(D, "z_after.npy")
    if not os.path.exists(zp):
        blocked.append((tile, "no z_after.npy")); continue
    z = np.load(zp)
    old = np.load(p).astype(bool)
    vt, src = valley_top_for(tile)
    zf = z.copy()
    nm = ~np.isfinite(zf)
    if nm.any():                       # same gap-fill the producer applies before cutting
        zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
    oldp50 = float(np.nanmedian(z[old & np.isfinite(z)])) if (old & np.isfinite(z)).any() else float("nan")
    if vt is None:
        blocked.append((tile, "no valley top in refcells.VALLEY_TOP_M"))
        print(f"{tile:<24}{'--':>12}{'BLOCKED':>22}{int(old.sum()):>11,}{'--':>11}{oldp50:>11.1f}")
        continue
    new, _ = floodplain_by_elevation(zf, vt)
    print(f"{tile:<24}{vt:>12.1f}{src:>22}{int(old.sum()):>11,}{int(new.sum()):>11,}{oldp50:>11.1f}")
    if A.write:
        bak = os.path.join(D, "floodplain_mask_histogram.npy.bak")
        if not os.path.exists(bak):
            shutil.copy(p, bak)
        np.save(p, new)
        with open(os.path.join(D, "floodplain_mask.json"), "w") as fh:
            json.dump({"method": "elevation cut: isfinite(z_after) & (z_after < valley_top_m)",
                       "valley_top_m": vt, "valley_top_source": f"refcells.VALLEY_TOP_M[{src!r}]",
                       "cells": int(new.sum()), "cells_before": int(old.sum()),
                       "rebuilt": "2026-09-15",
                       "why": "the previous mask used the WITHDRAWN histogram valley top; at "
                              "elbaext 83% of it was above the valley top"}, fh, indent=1)
            fh.write("\n")
    done.append(tile)

print(f"\n  {len(done)} rebuilt{'' if A.write else ' (report only; re-run with --write)'}, "
      f"{len(blocked)} BLOCKED")
if blocked:
    print("\n  BLOCKED -- these need a valley top STATED before they can be rebuilt:")
    for tile, why in blocked:
        D = f"data/derived/{tile}"
        zp = os.path.join(D, "z_after.npy")
        if os.path.exists(zp):
            z = np.load(zp); f = np.isfinite(z)
            print(f"    {tile:<22} {why}")
            print(f"      elevation p05 {np.percentile(z[f],5):.1f}  p25 {np.percentile(z[f],25):.1f}"
                  f"  p50 {np.percentile(z[f],50):.1f}  p75 {np.percentile(z[f],75):.1f}"
                  f"  p95 {np.percentile(z[f],95):.1f} m")
        else:
            print(f"    {tile:<22} {why}")
    print("\n  Add each to refcells.VALLEY_TOP_M, or pass one per tile. NOT invented here.")
