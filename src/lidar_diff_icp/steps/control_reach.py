#!/usr/bin/env python3
"""Report how many of a site's OWN survey's control marks sit in gen1 tiles on disk.

Writes ``control_reach.json`` into the tile directory. It is the gen1 mirror of
``data_completeness.json``: a site can build completely, pass all six pipeline steps, and
still be unable to carry a datum, and until this existed nothing in the graph said so.
mnrv is exactly that -- zero marks of any cover inside its single tile, discovered only by
a failed run on 2026-09-10.

REPORTED, NEVER JUDGED. No threshold is applied. The radius ladder is printed so the
choice of how far to reach is made against the counts, rather than a default being
inherited: how many marks a datum needs is a scientific decision.

It names the MISSING TILES rather than a region, because a mark's tie needs only a small
window (``5*res`` = 50 m; 300 m for the bridge's CSF) while MnGeo serves whole ~22 MB LAZ
with no COPC or .lax, so a tile is the download unit. At mnrv, 11 open marks within 25 km
sit in 11 tiles (~242 MB); the disc would be ~85 tiles (~1.9 GB). Fetch what the report
names, with scripts/fetch_tile.py.

    python -m lidar_diff_icp.steps.control_reach --site mnrv
    python -m lidar_diff_icp.steps.control_reach --site mnrv --radius-km 25 --write
"""
import argparse
import json
import os

from lidar_diff_icp import control_reach as CR
from lidar_diff_icp.sites import SITES, site as get_site

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--site", required=True, help="site name, or 'all'")
ap.add_argument("--radius-km", type=float, default=None,
                help="optional distance cut from the site centre. NO default: the ladder "
                     "below shows what each radius buys, and the choice is scientific.")
ap.add_argument("--covers", nargs="+", default=["L1O"],
                help="cover classes to count; default open ground only, the standing rule "
                     "for a datum. Pass ALL to count every cover.")
ap.add_argument("--write", action="store_true",
                help="write control_reach.json into the tile directory")
A = ap.parse_args()

names = list(SITES) if A.site == "all" else [A.site]
covers = None if [c.upper() for c in A.covers] == ["ALL"] else tuple(A.covers)

for name in names:
    s = get_site(name)
    r = CR.reach_for_site(s, covers=covers,
                          radius_m=None if A.radius_km is None else A.radius_km * 1000)
    print(CR.summary_line(r))
    for n in r.notes:
        print(f"    note: {n}")

    # the ladder: what each radius buys, so no radius is inherited as a default
    print(f"    {'radius':>9}{'marks':>8}{'on disk':>9}{'missing tiles':>15}{'~MB':>7}")
    for km in (5, 10, 15, 20, 25, 40, None):
        rr = CR.reach_for_site(s, covers=covers,
                               radius_m=None if km is None else km * 1000)
        lab = "all" if km is None else f"{km} km"
        print(f"    {lab:>9}{rr.n_marks_in_radius:>8}{rr.marks_reachable:>9}"
              f"{len(rr.missing):>15}{rr.approx_missing_mb:>7}")

    if r.missing:
        print(f"    missing tiles ({len(r.missing)}): {', '.join(r.missing[:12])}"
              f"{' ...' if len(r.missing) > 12 else ''}")
        print(f"    fetch with: python scripts/fetch_tile.py --tile <NAME> "
              f"--out {os.path.dirname(s.gen1)}")

    if A.write:
        out = os.path.join(s.tile_dir, "control_reach.json")
        os.makedirs(s.tile_dir, exist_ok=True)
        with open(out, "w") as fh:
            json.dump({k: (list(v) if isinstance(v, tuple) else v)
                       for k, v in r.__dict__.items()}, fh, indent=1)
            fh.write("\n")
        print(f"    wrote {out}")
