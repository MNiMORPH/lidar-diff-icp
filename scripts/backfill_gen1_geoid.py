#!/usr/bin/env python3
"""Record which GEOID each shipped corrections.json's gen1 was reduced in.

`difference_dem` writes `gen1_geoid_grid` since 2026-09-07. Files written before that do
not carry it, so a product cannot say what vertical frame it is on -- and that is exactly
the fact whose absence cost +54.87 mm at battlecreek, +27.75 at carlton and +26.39 at cook,
because `references.geoid_difference` defaulted `before_geoid` to GEOID03 for every survey.
Three sites (elba, whitewater, mnrv) were verified unaffected and deliberately NOT rebuilt,
so they are correct but silent. This labels them.

**The label is CERTIFIED, not asserted.** The geoid is not taken from the Site record and
written down. For each file the expected grid is resolved from `Site.gen1_project` through
`lidar_diff_icp.acquisitions`, `geoid_difference` is RE-RUN over that file's own bounds, and
the result is compared with the `cross_epoch_datum` the file already stores. Only a file
whose stored constant and tilts actually reproduce under that grid gets labelled. A file
that disagrees is REPORTED and left alone: it was built on a different frame, and it needs a
rebuild, not a sticker. That is the whole point -- a wrong geoid is invisible from inside a
tile, so a label that is merely assumed would make it invisible for good.

Idempotent: a file already carrying a matching value is verified and left alone; one
carrying a DIFFERENT value is reported as a conflict and never overwritten.

Tiles with no Site record (elbaext, elba_fulldensity, ...) are skipped and counted, because
nothing states which survey they came from. Guessing one is the error this script exists to
prevent.

    ./lidar-icp/bin/python scripts/backfill_gen1_geoid.py            # report only
    ./lidar-icp/bin/python scripts/backfill_gen1_geoid.py --write
"""
import argparse
import glob
import json
import os
import sys

from lidar_diff_icp import acquisitions, references          # noqa: E402
from lidar_diff_icp.sites import SITES                       # noqa: E402

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--glob", default="data/derived/*/corrections.json")
ap.add_argument("--write", action="store_true", help="write the key (default: report only)")
ap.add_argument("--tol-m", type=float, default=1e-6,
                help="agreement required between the file's stored cross-epoch constant "
                     "and the constant recomputed under the resolved geoid, in metres "
                     "(default 1e-6 = 0.001 mm). This is a REPRODUCTION gate, not a "
                     "physical tolerance: the same code over the same bounds with the same "
                     "grid should return the same number, so anything above numerical "
                     "noise means the file was built on a different frame. The achieved "
                     "difference is printed for every file, and none is excluded.")
A = ap.parse_args()

n_ok = n_written = n_skip = n_conflict = 0
print(f"{'tile':<20}{'survey':<21}{'geoid':<9}{'|stored-recomputed|':>21}  action")
for p in sorted(glob.glob(A.glob)):
    tile = os.path.basename(os.path.dirname(p))
    try:
        d = json.load(open(p))
    except Exception as e:
        print(f"{tile:<20}unreadable ({e})"); n_skip += 1; continue

    site = SITES.get(tile)
    if site is None or not site.gen1_project:
        print(f"{tile:<20}{'-':<21}{'-':<9}{'-':>21}  SKIP: no Site record naming a survey")
        n_skip += 1
        continue
    acq = acquisitions.for_project(site.gen1_project)      # raises on an unknown survey

    tie = d.get("cross_epoch_datum") or {}
    bounds = d.get("bounds")
    if not bounds or "const_m" not in tie:
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{'-':>21}"
              f"  SKIP: no bounds / cross_epoch_datum to certify against")
        n_skip += 1
        continue

    a, b, c = references.geoid_difference(bounds, 26915, before_geoid=acq.geoid_grid)
    diff = max(abs(a - tie["const_m"]),
               abs(b - tie.get("tilt_b_m_per_km", b)),
               abs(c - tie.get("tilt_c_m_per_km", c)))
    have = d.get("gen1_geoid_grid")

    if diff > A.tol_m:
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{diff:>21.3e}"
              f"  CONFLICT: does NOT reproduce -- built on another frame; REBUILD, not label")
        n_conflict += 1
    elif have == acq.geoid_grid:
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{diff:>21.3e}"
              f"  ok: already labelled, and it certifies")
        n_ok += 1
    elif have not in (None, acq.geoid_grid):
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{diff:>21.3e}"
              f"  CONFLICT: file says {have!r}; never overwritten")
        n_conflict += 1
    elif A.write:
        d["gen1_geoid_grid"] = acq.geoid_grid
        with open(p, "w") as fh:
            json.dump(d, fh, indent=1)
            fh.write("\n")
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{diff:>21.3e}  WROTE")
        n_written += 1
    else:
        print(f"{tile:<20}{acq.project_id:<21}{acq.geoid_model:<9}{diff:>21.3e}"
              f"  would write (re-run with --write)")
        n_written += 1

print(f"\n{n_ok} already correct, {n_written} {'written' if A.write else 'to write'}, "
      f"{n_skip} skipped, {n_conflict} CONFLICT")
if n_conflict:
    raise SystemExit(f"{n_conflict} file(s) do not reproduce under the survey's own geoid. "
                     f"They need rebuilding; labelling them would hide it.")
