"""What is the catchment radius, and do the marks it admits really belong to the line?

The catchment is the perpendicular half-width of the SEARCH around a flight-line track.
It is a compute bound on candidates, not a physical criterion: every candidate is then
confirmed or rejected by ``gen1_datum.assign_line_from_returns``, which reads the actual
point_source_id of the ground returns at the mark.

This measures the two things that decide how wide it should be:

1. the LOCAL line spacing, from the track model itself -- half of it is the vendor's
   class-12 seam, where bare-earth ground is cut and reassigned to the neighbouring line;
2. for every mark a wider catchment admits, the DOMINANT FRACTION of its ground returns,
   i.e. how decisively the returns actually place it on that line.
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--tracks", required=True)
    p.add_argument("--psids", type=int, nargs="+", required=True)
    p.add_argument("--covers", nargs="+", required=True)
    p.add_argument("--tiles", nargs="+", required=True)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--half-widths-m", type=float, nargs="+", required=True)
    p.add_argument("--collinear-sigma", type=float, required=True)
    a = p.parse_args(argv)

    ts = L.load_tracks(a.tracks)
    ctl = G.load_control()

    R = Run("How wide should the catchment be, and do the marks a wider one admits "
            "really belong to the line the returns assign them to?")
    R.input(a.tracks, role="flight-line tracks; used for the SEARCH and to measure the "
                           "local line spacing")
    R.input(ctl.origin, role="gen1's own 2008 control")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="andy")
    R.param("covers", tuple(a.covers), src="andy")
    R.param("half_widths_m", tuple(a.half_widths_m), src="MINE",
            why="the two candidate catchments being compared; the narrow one is half the "
                "measured line spacing (the class-12 seam), the wide one is a bound that "
                "lets the RETURNS do the discriminating")
    R.param("collinear_sigma", a.collinear_sigma, src="MINE",
            why="passes merged into physical lines within this many prediction-sd")
    R.column("pair", "adjacent flight lines whose tracks are being separated")
    R.column("spacing_m", "perpendicular distance between the two tracks at the site's "
                          "northing, m")
    R.column("half_spacing_m", "half of it -- the class-12 seam, where the vendor cuts "
                               "bare-earth ground over to the neighbouring line, m")
    R.column("point_id", "control mark admitted only by the WIDER catchment")
    R.column("cover", "MnDNR land-cover class, unitless code")
    R.column("line", "point_source_id the RETURNS assign, not the nearest track")
    R.column("dist_to_track_m", "perpendicular distance to the nearest searched track, m")
    R.column("dominant_frac", "share of ground returns at the mark carrying that line, "
                              "0-1; how decisively the returns place it")
    R.column("n_lines_at_mark", "distinct point_source_id values in the mark's window, count")
    R.banner()

    # 1. local spacing between adjacent Elba lines, measured from the tracks
    near = {}
    for psid in a.psids:
        cands = ts.by_psid(psid)
        if cands:
            near[psid] = min(cands, key=lambda q: S._track_distance(
                a.easting, a.northing, q.vertices))
    ks = sorted(near)
    rows = []
    for x, y in zip(ks, ks[1:]):
        vx = np.asarray(near[x].vertices, float)
        vy = np.asarray(near[y].vertices, float)
        # easting of each track at the site's northing
        ex = np.interp(a.northing, vx[np.argsort(vx[:, 1]), 1],
                       vx[np.argsort(vx[:, 1]), 0])
        ey = np.interp(a.northing, vy[np.argsort(vy[:, 1]), 1],
                       vy[np.argsort(vy[:, 1]), 0])
        d = abs(ex - ey)
        rows.append([f"{x}-{y}", f"{d:.0f}", f"{d/2:.0f}"])
    R.table(["pair", "spacing_m", "half_spacing_m"], rows)
    sp = np.array([float(r[1]) for r in rows])
    print(f"\n  measured spacing at the site: {sp.min():.0f}-{sp.max():.0f} m, "
          f"mean {sp.mean():.0f} m")
    print(f"  => the class-12 seam runs {sp.min()/2:.0f}-{sp.max()/2:.0f} m, "
          f"mean {sp.mean()/2:.0f} m")

    # 2. what a wider catchment admits, and how firmly the returns place it
    got = {}
    for hw in a.half_widths_m:
        sc, sites, meas, skipped, est = S.estimate(
            ts, psids=a.psids, easting=a.easting, northing=a.northing, scope="track",
            half_width_m=hw, covers=a.covers, tile_dirs=a.tiles, res=a.res,
            control=ctl, collinear_sigma=a.collinear_sigma)
        got[hw] = {m.point_id: m for m in S.marks_on_scope_psids(meas, a.psids)}
        print(f"  half_width {hw:7.0f} m -> {len(got[hw])} returns-confirmed marks")
    lo, hi = min(a.half_widths_m), max(a.half_widths_m)
    add = sorted(set(got[hi]) - set(got[lo]))
    tracks = ts.as_search_tracks()
    keys = S.site_scope(ts, psids=a.psids, easting=a.easting, northing=a.northing,
                        scope="track", collinear_sigma=a.collinear_sigma).track_keys
    rows2 = []
    for pid in add:
        m = got[hi][pid]
        d = min(S._track_distance(m.site.mark.easting, m.site.mark.northing, tracks[k])
                for k in keys)
        rows2.append([pid, m.site.mark.cover_class, m.line_id, f"{d:.0f}",
                      f"{m.line.dominant_fraction:.2f}", len(m.line.counts)])
    R.table(["point_id", "cover", "line", "dist_to_track_m", "dominant_frac",
             "n_lines_at_mark"], rows2)
    R.done(headline=f"seam {sp.min()/2:.0f}-{sp.max()/2:.0f} m; the wider catchment adds "
                    f"{len(add)} marks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
