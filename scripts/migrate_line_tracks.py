#!/usr/bin/env python3
"""Stamp a pre-2026-09-09 flight-line track file with the ACQUISITION it came from.

`groundtruth.lines.Pass` is keyed `project:psid.pass_index` since 2026-09-09, because
**psid is a per-project line number, not a global one**. gen1 is four acquisitions across
the six pilot sites, and cook uses psids 8-11 while battlecreek uses 1012-1014: keyed by
psid alone, a statewide TrackSet would eventually merge two different flight lines that
share a number, silently. `load_tracks` therefore REFUSES a file whose passes record no
project.

The project is NOT assumed. It is derived from the tile paths the file records, each
resolved to a Site and then to that Site's `gen1_project`, and the migration REFUSES unless
every tile agrees on one acquisition. A file spanning two surveys cannot be stamped with a
single project and must be re-derived per acquisition instead.

For `ground_control/data/gen1_line_tracks.json` the evidence is: 46 tiles, all under
`data/before/`, which is where elba's and whitewater's gen1 live, and both are
`lidar_semn2008`. Its psids (115-156, 1512, 1542, 10008-10012) are that survey's.

    ./lidar-icp/bin/python scripts/migrate_line_tracks.py <file>            # report
    ./lidar-icp/bin/python scripts/migrate_line_tracks.py <file> --write
"""
import argparse
import json
import os
import sys

from lidar_diff_icp.sites import SITES

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("path")
ap.add_argument("--write", action="store_true", help="write it (default: report only)")
A = ap.parse_args()

d = json.loads(open(A.path).read())
passes = d["passes"]
already = [q for q in passes if "project" in q]
if already and len(already) == len(passes):
    print(f"{A.path}: all {len(passes)} passes already record a project "
          f"({sorted({q['project'] for q in passes})}); nothing to do")
    raise SystemExit(0)
if already:
    raise SystemExit(f"{A.path}: {len(already)} of {len(passes)} passes record a project "
                     f"and the rest do not. Mixed state; not touching it.")

tiles = d.get("tiles") or []
if not tiles:
    raise SystemExit(f"{A.path}: records no tiles, so nothing identifies its acquisition. "
                     f"Re-derive it with derive_tracks(project=...) instead of guessing.")

# tile path -> the Site whose gen1 lives in that directory -> that Site's survey
by_dir = {}
for name, s in SITES.items():
    if s.gen1_project:
        by_dir.setdefault(os.path.dirname(s.gen1), set()).add(s.gen1_project)
found, unknown = set(), set()
for t in tiles:
    projs = by_dir.get(os.path.dirname(t))
    (found.update(projs) if projs else unknown.add(os.path.dirname(t)))

print(f"{A.path}: {len(passes)} passes, {len(tiles)} tiles")
print(f"  tile directories -> surveys: "
      f"{ {k: sorted(v) for k, v in by_dir.items() if k in {os.path.dirname(t) for t in tiles}} }")
if unknown:
    raise SystemExit(f"  REFUSED: {sorted(unknown)} matches no Site's gen1 directory, so "
                     f"the acquisition is not established for those tiles.")
if len(found) != 1:
    raise SystemExit(f"  REFUSED: tiles span {sorted(found)}. One file cannot carry one "
                     f"project; re-derive per acquisition.")
project = found.pop()
psids = sorted({q["psid"] for q in passes})
print(f"  every tile resolves to ONE survey: {project}")
print(f"  psids {psids[0]}-{psids[-1]} ({len(psids)} distinct) would be stamped with it")

if not A.write:
    print("  (report only; re-run with --write)")
    raise SystemExit(0)
for q in passes:
    q["project"] = project
d["migrated"] = {"project": project, "on": "2026-09-09",
                 "evidence": f"all {len(tiles)} tiles under "
                             f"{sorted({os.path.dirname(t) for t in tiles})}, which is the "
                             f"gen1 directory of Sites on {project}"}
with open(A.path, "w") as fh:
    json.dump(d, fh, indent=1)
    fh.write("\n")
print(f"  WROTE project={project!r} onto {len(passes)} passes, with the evidence recorded")
