"""Derive the gen1 flight-line tracks and COMMIT them into the repository.

    ./lidar-icp/bin/python ground_control/run_derive_tracks.py \
        --tiles 'data/before/*.laz' --exclude-substring merged \
        --out ground_control/data/gen1_line_tracks.json --chunk-size 2000000

The four method parameters are inherited from analysis/groundtruth/gen1_line_tracks.py,
not re-chosen, so these tracks are the object that module built.  Pass --stride etc. to
override and the banner will record it.

Nothing is downloaded.  Tiles are streamed, not read whole.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))

import lines as L  # noqa: E402
from trust.provenance import Run  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tiles", required=True, help="glob for gen1 tiles")
    p.add_argument("--exclude-substring", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--chunk-size", type=int, required=True)
    for k, v in L.INHERITED_PARAMS.items():
        p.add_argument(f"--{k.replace('_', '-')}", type=type(v), default=v)
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    tiles = sorted(glob.glob(a.tiles))
    if a.exclude_substring:
        tiles = [t for t in tiles if a.exclude_substring not in t]
    if not tiles:
        raise SystemExit(f"no tiles matched {a.tiles!r}")

    R = Run("Where does each gen1 flight-line PASS run, so that marks it could have "
            "seen can be found anywhere along it?")
    for t in tiles[:3]:
        R.input(t, role="gen1 2008 delivered tile; the near-nadir returns give the "
                        "flight track, and point_source_id names the pass")
    R.param("n_tiles", len(tiles), src="repo",
            why="every gen1 tile already on disk; nothing is downloaded")
    R.param("tiles_glob", a.tiles, src="andy")
    for k in L.INHERITED_PARAMS:
        got, base = getattr(a, k), L.INHERITED_PARAMS[k]
        R.param(k, got, src="repo" if got == base else "MINE",
                why=("inherited verbatim from analysis/groundtruth/gen1_line_tracks.py "
                     "so these tracks are the object that module built"
                     if got == base else
                     f"OVERRIDDEN on the command line; that module used {base}"))
    R.param("chunk_size", a.chunk_size, src="MINE",
            why="points per streamed read; affects peak memory on this shared laptop "
                "only, not the result")
    R.column("psid_pass", "point_source_id and pass index; a psid is NOT a flight line, "
                          "it is reused across missions and split at gps_time gaps")
    R.column("n_bins", "gps_time bins surviving min_bin_points, count")
    R.column("span_km", "straight-line distance between the pass's first and last "
                        "binned centroid, km")
    R.column("heading_deg", "bearing from first to last centroid, degrees clockwise "
                            "from north")
    R.column("gap_before_s", "gps_time gap to the previous pass of the SAME psid, "
                             "seconds; blank for the first pass")
    R.column("resid_med_m", "median perpendicular distance of the binned centroids to "
                            "the fitted straight track, m")
    R.column("resid_p95_m", "95th percentile of the same, m")
    R.notes.append("A gap is not by itself proof of a second pass: it is equally a "
                   "stretch of line under tiles we do not hold. The gap is printed so "
                   "the caller judges.")
    R.notes.append("These tracks TARGET marks. Which line actually illuminated a mark "
                   "is settled by gen1_datum.assign_line_from_returns from the returns, "
                   "never by distance to a track.")
    R.banner()

    def progress(i, n, path, kept):
        print(f"  [{i:>2}/{n}] {Path(path).name}  near-nadir kept {kept:,}", flush=True)

    P = {k: getattr(a, k) for k in L.INHERITED_PARAMS}
    ts = L.derive_tracks(tiles, chunk_size=a.chunk_size, progress=progress, **P)

    prev = {}
    rows = []
    for p in sorted(ts.passes, key=lambda q: (q.psid, q.pass_index)):
        gap = "" if p.psid not in prev else f"{p.t0 - prev[p.psid]:.0f}"
        prev[p.psid] = p.t1
        rows.append([p.key, p.n_bins, f"{p.span_km:.2f}", f"{p.heading_deg:.1f}",
                     gap, f"{p.resid_med_m:.1f}", f"{p.resid_p95_m:.1f}"])
    R.table(["psid_pass", "n_bins", "span_km", "heading_deg", "gap_before_s",
             "resid_med_m", "resid_p95_m"], rows)

    multi = sorted({p.psid for p in ts.passes
                    if len([q for q in ts.passes if q.psid == p.psid]) > 1})
    print()
    print(f"  psids present            : {len({p.psid for p in ts.passes})}")
    print(f"  passes                   : {len(ts.passes)}")
    print(f"  psids with >1 pass       : {len(multi)}  {multi}")
    print(f"    ^ each of these would be merged into one 'line' by any code keying on "
          f"point_source_id alone")

    out = L.save_tracks(ts, a.out)
    print(f"  wrote {out}  ({out.stat().st_size/1024:.0f} kB)")
    R.done(headline=f"{len(ts.passes)} passes over {len({p.psid for p in ts.passes})} "
                    f"psids from {ts.n_tiles_read} tiles; "
                    f"{len(multi)} psids carry more than one pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
