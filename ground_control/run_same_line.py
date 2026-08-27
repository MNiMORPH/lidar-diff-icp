"""gen1's datum at a site from its OWN flight lines, and the psid-scope counterfactual.

    ./lidar-icp/bin/python ground_control/run_same_line.py \
        --easting 578762.8 --northing 4884487.6 \
        --tracks ground_control/data/gen1_line_tracks.json \
        --psids 133 134 135 136 137 138 --covers L1O L5U \
        --tiles data/before --res 5.0

Runs BOTH scopes.  They differ in one thing -- which tracks the mark search walks -- so
any difference between the two answers is the search and nothing else.
"""

from __future__ import annotations

import argparse
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


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--tracks", required=True)
    p.add_argument("--psids", type=int, nargs="+", required=True)
    p.add_argument("--covers", nargs="+", required=True,
                   help="MnDNR classes; no default -- pooling covers measures canopy")
    p.add_argument("--tiles", nargs="+", required=True)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--half-width-m", type=float, default=S.SEAM_HALF_SPACING_M)
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    ts = L.load_tracks(a.tracks)
    control = G.load_control()

    R = Run("What is gen1's datum at this site, from marks on the site's OWN flight "
            "lines, and does searching by point_source_id instead of by PASS change it?")
    R.input(a.tracks, role="gen1 flight-line tracks, one per PASS; TARGETS marks, never "
                           "assigns them")
    R.input(control.origin, role="gen1's own 2008 MnGeo/MnDNR control; dnr_error_m is "
                                 "Control Z - Surface Z")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="andy",
            why="the flight lines the site's tile is built from")
    R.param("covers", tuple(a.covers), src="andy",
            why="supplied on the command line; pooling covers measures canopy, not datum")
    R.param("half_width_m", a.half_width_m, src="repo",
            why="half the MEASURED gen1 line spacing (942/987/932/988/956 m, mean 961) "
                "= the vendor's class-12 seam at which bare-earth ground is cut, from "
                "analysis/GEN1_DATUM_MORE_MARKS.md section 1. Not a chosen radius")
    R.param("res_m", a.res, src="repo", why="corrections.json res_m; sets the radius ladder")
    R.param("scopes", ("pass", "psid"), src="MINE",
            why="both are run; they differ ONLY in which tracks the search walks, so the "
                "difference isolates the search")
    R.param("swath_constants", "none (mode=per_line)", src="repo",
            why="per_line treats each line's constant as unknown, which it is")
    R.column("scope", "which tracks the mark search walked: pass = only each psid's "
                      "nearest pass to the site; psid = every pass of every psid")
    R.column("tracks", "number of tracks searched")
    R.column("candidates", "marks within half_width_m of any searched track, count")
    R.column("measured", "candidates whose tile is on disk and which yielded a tie, count")
    R.column("on_site_psids", "of those, marks the RETURNS place on one of --psids, count")
    R.column("n_lines", "distinct point_source_id among them, count")
    R.column("datum_mm", "constant to ADD to gen1 at the site, mm; mean over line means; "
                         "positive = gen1 reads LOW")
    R.column("se_mm", "SE of the mean over LINES: sd of the line means / sqrt(n_lines). "
                      "The line, not the mark, is the unit of replication")
    R.column("max_km", "farthest measured mark from the site, km -- how far the search reached")
    R.banner()

    rows, keep = [], {}
    for scope in S.SCOPES:
        sc, sites, meas, skipped, est = S.estimate(
            ts, psids=a.psids, easting=a.easting, northing=a.northing, scope=scope,
            half_width_m=a.half_width_m, covers=a.covers, tile_dirs=a.tiles,
            res=a.res, control=control)
        on = S.marks_on_scope_psids(meas, a.psids)
        d = [np.hypot(m.site.mark.easting - a.easting,
                      m.site.mark.northing - a.northing) / 1000.0 for m in on]
        keep[scope] = (sc, meas, on, est)
        rows.append([scope, sc.n_tracks, len(sites), len(meas), len(on),
                     len({m.line_id for m in on}),
                     f"{est.value_mm:+.2f}", f"{est.se_mm:.2f}",
                     f"{max(d):.1f}" if d else "--"])
    R.table(["scope", "tracks", "candidates", "measured", "on_site_psids", "n_lines",
             "datum_mm", "se_mm", "max_km"], rows)

    print()
    for scope in S.SCOPES:
        sc, meas, on, est = keep[scope]
        print(f"  scope={scope}: {sc.note}")
        print(f"    tracks: {', '.join(sc.track_keys)}")
        if sc.dropped_track_keys:
            print(f"    dropped as other lines reusing the id: "
                  f"{', '.join(sc.dropped_track_keys)}")
    print()
    print("  per-mark ties, scope=pass (marks the RETURNS put on the site's psids):")
    print("    %-24s %5s %6s %9s %9s" % ("point_id", "cover", "line", "km", "tie_mm"))
    for m in sorted(keep["pass"][2], key=lambda z: z.line_id):
        km = np.hypot(m.site.mark.easting - a.easting,
                      m.site.mark.northing - a.northing) / 1000.0
        print("    %-24s %5s %6d %9.1f %+9.1f"
              % (m.point_id, m.site.mark.cover_class, m.line_id, km, m.tie_mm))

    extra = ({m.point_id for m in keep["psid"][2]} - {m.point_id for m in keep["pass"][2]})
    print(f"\n  marks the psid scope adds that the pass scope excludes: {len(extra)}")
    if extra:
        for m in keep["psid"][2]:
            if m.point_id in extra:
                km = np.hypot(m.site.mark.easting - a.easting,
                              m.site.mark.northing - a.northing) / 1000.0
                print(f"    {m.point_id:<24s} line {m.line_id}  {km:6.1f} km away  "
                      f"tie {m.tie_mm:+.1f} mm")

    e = keep["pass"][3]
    R.done(headline=f"gen1 at the site, scope=pass: {e.value_mm:+.2f} +/- {e.se_mm:.2f} mm "
                    f"over {e.n_lines} lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
