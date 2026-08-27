"""gen1's datum at a site with NO catchment radius: every mark in a tile, assigned by returns.

    ./lidar-icp/bin/python ground_control/run_datum_by_returns.py \
        --easting 578762.8 --northing 4884487.6 --psids 133 134 135 136 137 138 \
        --covers L1O L5U --tiles data/before --res 5.0 \
        --tracks ground_control/data/gen1_line_tracks.json --collinear-sigma 3

The catchment was only ever a compute bound on candidates -- assign_line_from_returns does
the assigning and can only reject -- so this bounds by "is the tile on disk" instead, which
is simpler, strictly more complete, and free of the cover-mix confound a radius introduces.

Tracks are still required, for ONE job: disambiguating a reused point_source_id. gps_time
cannot substitute (correlation with the collinearity sigma is -0.32, the wrong sign).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))

import lines as L  # noqa: E402
import same_line as S  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G  # noqa: E402
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--psids", type=int, nargs="+", required=True)
    p.add_argument("--covers", nargs="+", required=True)
    p.add_argument("--tiles", nargs="+", required=True)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--tracks", required=True)
    p.add_argument("--collinear-sigma", type=float, nargs="+", required=True)
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    ts = L.load_tracks(a.tracks)
    control = G.load_control()

    R = Run("gen1's datum at this site with NO catchment radius: every control mark "
            "inside a tile we hold, assigned to a line by its own ground returns.")
    R.input(control.origin, role="gen1's own 2008 control")
    R.input(a.tracks, role="flight-line tracks; used ONLY to disambiguate a reused "
                           "point_source_id, never to search or to assign")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="andy",
            why="the lines the site's tile is built from, per corrections.json")
    R.param("covers", tuple(a.covers), src="andy")
    R.param("catchment_radius_m", "NONE -- removed", src="andy",
            why="it was a compute bound, not a criterion; the returns assign, and a "
                "radius silently shifted the cover mix (at Elba all 4 marks it added "
                "were urban)")
    R.param("collinear_sigma", tuple(a.collinear_sigma), src="MINE",
            why="how many prediction-sd two passes of one psid may sit apart and still "
                "be one line; swept, no default")
    R.param("res_m", a.res, src="repo")
    R.column("sigma", "collinearity cut used for reused-psid disambiguation, unitless")
    R.column("candidates", "control marks of these covers inside a tile on disk, count")
    R.column("measured", "of those, marks that yielded a tie, count")
    R.column("kept", "marks the RETURNS place on one of --psids AND that survive "
                     "reused-psid disambiguation, count")
    R.column("n_lines", "distinct point_source_id among the kept marks, count")
    R.column("datum_mm", "constant to ADD to gen1 at the site, mm; mean over line means; "
                         "positive = gen1 reads LOW")
    R.column("se_mm", "SE of the mean over LINES, mm -- the line is the unit of "
                      "replication because marks under one line share its constant")
    R.column("max_km", "farthest kept mark from the site, km")
    R.banner()

    rows, keep = [], {}
    for sig in a.collinear_sigma:
        meas, kept, rejected, est = S.estimate_by_returns(
            ts, psids=a.psids, easting=a.easting, northing=a.northing,
            covers=a.covers, tile_dirs=a.tiles, res=a.res, collinear_sigma=sig,
            control=control)
        d = [np.hypot(m.site.mark.easting - a.easting,
                      m.site.mark.northing - a.northing) / 1000.0 for m in kept]
        keep[sig] = (meas, kept, rejected, est)
        rows.append([f"{sig:g}", len(S.marks_in_tiles(control, a.tiles, covers=a.covers)),
                     len(meas), len(kept), est.n_lines, f"{est.value_mm:+.2f}",
                     f"{est.se_mm:.2f}", f"{max(d):.1f}" if d else "--"])
    R.table(["sigma", "candidates", "measured", "kept", "n_lines", "datum_mm", "se_mm",
             "max_km"], rows)

    sig0 = a.collinear_sigma[0]
    meas, kept, rejected, est = keep[sig0]
    print(f"\n  at sigma={sig0:g}, why candidates were rejected:")
    import collections
    why = collections.Counter(r[2] for r in rejected)
    for k, v in why.most_common():
        print(f"    {v:3d}  {k}")
    print(f"\n  kept marks (line, km, tie):")
    for m in sorted(kept, key=lambda z: z.line_id):
        km = np.hypot(m.site.mark.easting - a.easting,
                      m.site.mark.northing - a.northing) / 1000.0
        print(f"    {m.point_id:<24s} {m.site.mark.cover_class}  line {m.line_id}  "
              f"{km:6.1f} km  {m.tie_mm:+8.1f} mm")
    if a.out:
        Path(a.out).write_text(json.dumps(dict(
            sigma=sig0, value_mm=est.value_mm, se_mm=est.se_mm, n_marks=est.n_marks,
            n_lines=est.n_lines,
            marks=[dict(point_id=m.point_id, cover=m.site.mark.cover_class,
                        line=int(m.line_id), tie_mm=m.tie_mm) for m in kept]), indent=1))
        print(f"\n  wrote {a.out}")
    R.done(headline=f"no-catchment datum {est.value_mm:+.2f} +/- {est.se_mm:.2f} mm "
                    f"({est.n_marks} marks, {est.n_lines} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
