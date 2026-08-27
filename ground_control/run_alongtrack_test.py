"""Does a same-line mark 30 km away still say anything about the site?

The pass-scope estimate reaches 30.6 km along track.  The per-line model assumes a line
carries ONE constant along its whole length, but this repo has measured within-line
along-track drift (analysis/ABSOLUTE_BASIS_ELBA.md carries 15.5 mm/km on the elbaext
swaths, and kept its anchors within 0.33 km of the site for that reason).

The test: regress each mark's tie on its ALONG-TRACK distance from the site, measured as
arc length along that mark's own fitted flight-line track.

    ./lidar-icp/bin/python ground_control/run_alongtrack_test.py \
        --easting 578762.8 --northing 4884487.6 \
        --tracks ground_control/data/gen1_line_tracks.json \
        --psids 133 134 135 136 137 138 --covers L1O L5U --tiles data/before --res 5.0

Distance is confounded with line here -- the near marks and the far marks sit on
different lines -- so the pooled slope is reported AND the within-line slopes are
reported beside it.  With this many marks neither may resolve anything; that is a result,
not a failure, and the script says so rather than fitting harder.
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


def arclength_of(px, py, vertices):
    """Arc length along the polyline of the projection of (px, py) onto it, metres."""
    v = np.asarray(vertices, float)
    a, b = v[:-1], v[1:]
    ab = b - a
    seg = np.hypot(ab[:, 0], ab[:, 1])
    cum = np.r_[0.0, np.cumsum(seg)]
    denom = np.maximum((ab ** 2).sum(1), 1e-9)
    t = np.clip(((px - a[:, 0]) * ab[:, 0] + (py - a[:, 1]) * ab[:, 1]) / denom, 0.0, 1.0)
    c = a + t[:, None] * ab
    d = np.hypot(px - c[:, 0], py - c[:, 1])
    k = int(np.argmin(d))
    return float(cum[k] + t[k] * seg[k]), float(d[k])


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
    a = p.parse_args(argv)

    ts = L.load_tracks(a.tracks)
    control = G.load_control()

    R = Run("Does a mark on the site's own flight line, but tens of km along it, still "
            "carry the site's datum -- or has the line drifted in between?")
    R.input(a.tracks, role="gen1 flight-line tracks, one per pass; used here to measure "
                           "ARC LENGTH along the line, not to assign marks")
    R.input(control.origin, role="gen1's own 2008 control; dnr_error_m is Control Z - Surface Z")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="andy")
    R.param("covers", tuple(a.covers), src="andy")
    R.param("scope", "pass", src="repo",
            why="the defensible scope from ground_control/same_line.py; psid scope mixes "
                "different physical lines and would confound the distance test outright")
    R.param("reference_rate_mm_per_km", 15.5, src="repo",
            why="within-swath along-track drift measured on the elbaext swaths, "
                "analysis/ABSOLUTE_BASIS_ELBA.md uncertainty budget. Used only as a "
                "yardstick to compare the fitted slope against; nothing is fitted to it")
    R.column("point_id", "control mark id")
    R.column("line", "point_source_id from the ground returns at the mark")
    R.column("along_km", "arc length along that line's own fitted track between the "
                         "mark's projection and the site's projection, km")
    R.column("offtrack_m", "perpendicular distance of the mark from that track, m")
    R.column("tie_mm", "surveyed - lidar at the mark, mm; positive = gen1 reads LOW")
    R.column("group", "which regression a row enters: pooled, and its line's own")
    R.notes.append("Along-track distance is CONFOUNDED with line here: the near marks and "
                   "far marks sit on different lines. The pooled slope therefore mixes a "
                   "distance effect with a per-line offset and cannot separate them.")
    R.banner()

    sc, sites, meas, skipped, est = S.estimate(
        ts, psids=a.psids, easting=a.easting, northing=a.northing, scope="pass",
        half_width_m=S.SEAM_HALF_SPACING_M, covers=a.covers, tile_dirs=a.tiles,
        res=a.res, control=control)
    on = S.marks_on_scope_psids(meas, a.psids)
    tracks = {int(k.split(".")[0]): ts.as_search_tracks()[k] for k in sc.track_keys}

    rows, rec = [], []
    for m in sorted(on, key=lambda z: (z.line_id, z.point_id)):
        v = tracks[int(m.line_id)]
        s_mark, off = arclength_of(m.site.mark.easting, m.site.mark.northing, v)
        s_site, _ = arclength_of(a.easting, a.northing, v)
        along = abs(s_mark - s_site) / 1000.0
        rec.append((m.point_id, int(m.line_id), along, off, m.tie_mm))
        rows.append([m.point_id, m.line_id, f"{along:.2f}", f"{off:.0f}",
                     f"{m.tie_mm:+.1f}", f"pooled+line{m.line_id}"])
    R.table(["point_id", "line", "along_km", "offtrack_m", "tie_mm", "group"], rows)

    x = np.array([r[2] for r in rec]); y = np.array([r[4] for r in rec])
    print()
    print("  POOLED regression of tie on along-track distance (CONFOUNDED with line):")
    if x.size >= 3 and np.ptp(x) > 0:
        b, a0 = np.polyfit(x, y, 1)
        yh = a0 + b * x
        ss = float(((y - yh) ** 2).sum()); st = float(((y - y.mean()) ** 2).sum())
        r2 = 1 - ss / st if st > 0 else float("nan")
        se_b = float(np.sqrt(ss / max(x.size - 2, 1) / ((x - x.mean()) ** 2).sum()))
        print(f"    slope {b:+.2f} +/- {se_b:.2f} mm/km   intercept {a0:+.2f} mm   "
              f"R2 {r2:.3f}   n {x.size}")
        print(f"    |t| on the slope = {abs(b/se_b):.2f} on {x.size-2} d.f.")
        print(f"    reference within-swath rate, elbaext: 15.5 mm/km")
    print()
    print("  WITHIN-LINE slopes (the only version distance is not confounded in):")
    for ln in sorted({r[1] for r in rec}):
        g = [r for r in rec if r[1] == ln]
        xs = np.array([r[2] for r in g]); ys = np.array([r[4] for r in g])
        if xs.size >= 2 and np.ptp(xs) > 0:
            bb = np.polyfit(xs, ys, 1)[0] if xs.size >= 2 else float("nan")
            print(f"    line {ln}: n={xs.size}  along {xs.min():.1f}-{xs.max():.1f} km  "
                  f"ties {', '.join('%+.1f' % t for t in ys)}  slope {bb:+.2f} mm/km"
                  + ("   [2 points: a slope through 2 points has no residual d.f.]"
                     if xs.size == 2 else ""))
        else:
            print(f"    line {ln}: n={xs.size}  along {xs.min():.1f} km  "
                  f"tie {ys[0]:+.1f}  -- no slope from one mark")
    print()
    near = y[x <= np.median(x)]; far = y[x > np.median(x)]
    print(f"  Split at the median along-track distance ({np.median(x):.1f} km):")
    print(f"    near n={near.size} mean {near.mean():+.1f} mm | "
          f"far n={far.size} mean {far.mean():+.1f} mm | "
          f"difference {near.mean()-far.mean():+.1f} mm")
    R.done(headline=f"pooled slope on {x.size} marks; see within-line rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
