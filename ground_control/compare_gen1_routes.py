"""Every defensible route to gen1's datum at a site, side by side.

There is no single answer to pick from; there are routes with different assumptions, and
the choice between them is the caller's.  This computes all of them in one run so they
are produced the SAME WAY and are actually comparable, and states for each what it
assumes and how it can fail.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        ground_control/compare_gen1_routes.py --easting 578762.8 --northing 4884487.6 \
        --tracks ground_control/data/gen1_line_tracks.json --psids 133 134 135 136 137 138 \
        --covers L1O L5U --tiles data/before --res 5.0 \
        --field-product ground_control/products/elba__gen1__open__datum.json \
        --bridge-product ground_control/products/bridge_wide_L1O.json \
        --out ground_control/products/gen1_routes_elba.json
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

ASSUMES = {
    "residual_field": "the residual is a smooth spatial field with NO flight-line term, "
                      "so 963 marks anywhere inform this site",
    "same_line_pass": "only the site's own flight-line PASSES carry its datum; a line has "
                      "one constant along its length",
    "same_line_psid": "same, but point_source_id names a line -- which it does not",
}
FAILS_IF = {
    "residual_field": "the per-line structure is real (measured: F = 8.63, p < 0.001), in "
                      "which case this averages across lines that genuinely differ",
    "same_line_pass": "a line's constant drifts along track (tested: slope +0.74 +/- 2.37 "
                      "mm/km, and 15.5 mm/km excluded at 6.23 sigma), or 7 marks are too few",
    "same_line_psid": "always -- two passes sharing a psid sit 10.4-83.3 km apart and are "
                      "different physical lines",
}


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
    p.add_argument("--field-product", required=True)
    p.add_argument("--bridge-product", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    R = Run("What are gen1's defensible datum routes at this site, and how far apart are "
            "they?")
    R.input(a.field_product, role="kriged residual-field datum at this site, open ground")
    R.input(a.bridge_product, role="per-mark bridge; its population mean carries a "
                                   "delivered-surface constant onto OUR surface")
    R.input(a.tracks, role="flight-line tracks, one per PASS")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("covers", tuple(a.covers), src="andy")
    R.param("psids", tuple(a.psids), src="andy")
    R.param("routes", tuple(ASSUMES), src="MINE",
            why="every route reachable from committed code; none is preferred here")
    R.column("route", "which estimator produced the constant")
    R.column("n", "marks entering it, count")
    R.column("unit", "the unit of replication its uncertainty is over")
    R.column("delivered_mm", "constant to ADD to the DELIVERED gen1 surface, mm")
    R.column("uncert_mm", "its stated uncertainty, mm -- see 'unit' for what OF")
    R.column("bridged_mm", "the same constant carried onto OUR reconstructed surface, mm "
                           "(delivered + bridge)")
    R.column("reproducible", "whether committed code regenerates it")
    R.banner()

    bridge = np.array([m["bridge_mm"] for m in
                       json.loads(Path(a.bridge_product).read_text())["marks"]])
    b_mean = float(bridge.mean())
    b_se = float(bridge.std(ddof=1) / np.sqrt(bridge.size))

    fp = json.loads(Path(a.field_product).read_text())
    rows, recs = [], []

    def add(route, n, unit, val, unc, repro):
        rows.append([route, n, unit, f"{val:+.2f}", f"{unc:.2f}",
                     f"{val + b_mean:+.2f}", repro])
        recs.append(dict(route=route, n=n, unit=unit, delivered_mm=val, uncert_mm=unc,
                         bridged_mm=val + b_mean, assumes=ASSUMES[route],
                         fails_if=FAILS_IF[route], reproducible=repro))

    add("residual_field", fp["n_marks"], "prediction sd of the field",
        fp["constant_mm"], fp["sd_field_mm"], "yes")

    ts = L.load_tracks(a.tracks)
    control = G.load_control()
    for scope in ("pass", "psid"):
        sc, sites, meas, skipped, est = S.estimate(
            ts, psids=a.psids, easting=a.easting, northing=a.northing, scope=scope,
            half_width_m=S.SEAM_HALF_SPACING_M, covers=a.covers, tile_dirs=a.tiles,
            res=a.res, control=control)
        add(f"same_line_{scope}", est.n_marks, f"SE over {est.n_lines} lines",
            est.value_mm, est.se_mm, "yes")

    R.table(["route", "n", "unit", "delivered_mm", "uncert_mm", "bridged_mm",
             "reproducible"], rows)
    v = np.array([r["delivered_mm"] for r in recs])
    print()
    print(f"  spread across routes: {np.ptp(v):.2f} mm ({v.min():+.2f} .. {v.max():+.2f})")
    print(f"  bridge applied to all: {b_mean:+.2f} +/- {b_se:.2f} mm "
          f"({bridge.size} open marks)")
    print()
    for r in recs:
        print(f"  {r['route']}")
        print(f"     ASSUMES  {r['assumes']}")
        print(f"     FAILS IF {r['fails_if']}")
    Path(a.out).write_text(json.dumps(dict(routes=recs, bridge_mm=b_mean,
                                           bridge_se_mm=b_se), indent=1))
    print(f"\n  wrote {a.out}")
    R.done(headline=f"{len(recs)} routes spanning {np.ptp(v):.1f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
