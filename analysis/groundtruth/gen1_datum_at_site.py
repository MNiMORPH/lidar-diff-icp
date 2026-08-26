"""gen1's absolute vertical datum at any MN site, from the 2008 MnGeo control.

The driver for :mod:`lidar_diff_icp.groundtruth.gen1_datum`. It discovers the control
marks around a site, says which tiles are needed and which are already on disk, measures
every mark whose tile IS on disk, and combines them in whichever mode the caller asks
for. It never downloads a tile and never runs CSF.

Every table below is printed through ``trust/provenance.py``, so no column appears whose
meaning the script has not stated, and no parameter appears without where it came from.

    ./lidar-icp/bin/python analysis/groundtruth/gen1_datum_at_site.py \
        --easting 579705.72 --northing 4883677.71 --radius-km 20 \
        --tiles data/before --mode per_line

    ... --mode common_datum --corrections data/derived/elba_fulldensity/corrections.json
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from trust.provenance import Run                                        # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G                  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--easting", type=float, required=True)
    ap.add_argument("--northing", type=float, required=True)
    ap.add_argument("--radius-km", type=float, required=True,
                    help="how far to reach for control. No default: it is a statement "
                         "about how far the datum may be assumed constant.")
    ap.add_argument("--tiles", default="data/before",
                    help="directory searched for <tile>.laz; nothing is ever downloaded")
    ap.add_argument("--mode", choices=G.MODES, required=True)
    ap.add_argument("--corrections", default=None,
                    help="a tile's corrections.json, for --mode common_datum")
    ap.add_argument("--cover", default=None,
                    help="restrict to these MnDNR cover classes, comma separated "
                         "(L1O open / L2T weeds+crops / L3B brush / L4F forest / L5U urban). "
                         "There is NO default: pooling covers measures canopy, not datum, "
                         "and picking one is the caller's decision, recorded here.")
    ap.add_argument("--lateral-shift", default=None,
                    help="OPT-IN gen2-derived Nuth-Kaeaeb shift 'dx,dy' in metres. "
                         "Off by default: a cross-epoch term does not belong here.")
    ap.add_argument("--index-cache", default="data/mn_tile_centroids.csv")
    ap.add_argument("--out", default=None, help="write the estimate as JSON")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    R = Run("what is gen1's absolute vertical datum at this site, measured against the "
            "2008 MnGeo control the acquisition was itself validated on?")

    control = G.load_control()
    R.input(control.origin, role="the 2008 MnDNR control transcribed from the eight "
                                 "MnGeo county validation reports, NAVD88(GEOID03)")
    R.input(a.index_cache, role="MnGeo statewide tile centroid index, for naming the tile "
                                "each mark falls in")

    swath_const, swath_src = None, ""
    if a.mode == "common_datum":
        if not a.corrections:
            ap.error("--mode common_datum needs --corrections <tile>/corrections.json")
        R.input(a.corrections, role="per-swath alignment constants solved by "
                                    "coreg.align_swaths, gen1-internal")
        swath_const, swath_src = G.swath_constants_from_corrections(a.corrections)

    lat = None
    if a.lateral_shift:
        lat = tuple(float(v) for v in a.lateral_shift.split(","))

    R.param("site_easting_m", a.easting, src="andy", why="the site asked about")
    R.param("site_northing_m", a.northing, src="andy", why="the site asked about")
    R.param("search_radius_km", a.radius_km, src="andy",
            why="supplied on the command line; the module has no default radius")
    R.param("cover_class", a.cover or "ALL COVERS POOLED", src="andy",
            why="supplied on the command line; pooling covers measures canopy, not datum")
    R.param("mode", a.mode, src="andy", why="per_line = swath constants unknown; "
                                            "common_datum = align_swaths constants applied")
    R.param("swath_constants", swath_src or "none", src="repo",
            why=swath_src or "mode=per_line applies none")
    R.param("ground_source", "vendor class 2", src="repo",
            why="ASPRS bare earth as delivered; CSF is ~460 s/tile and the 2026-08-26 "
                "control run measured the ground-source effect at 6.5 mm median absolute "
                "over 16 marks")
    R.param("geoid_shift_m", 0.0, src="repo",
            why="both sides are NAVD88(GEOID03); assert_no_geoid_conversion checks it "
                "per mark and raises rather than converting")
    R.param("lateral_shift_m", lat, src="andy" if lat else "repo",
            why="gen2-derived Nuth-Kaeaeb shift, taken on explicit request" if lat else
                "NOT applied: a cross-epoch term has no place in gen1 vs its own control")
    R.param("res_m", 5.0, src="repo", why="corrections.json res_m; sets the radius ladder")
    R.param("threshold_on_the_screen", "none", src="repo",
            why="the screen statistics are printed for every mark and NOTHING is cut on "
                "them; the radius-spread screen was measured not to reduce site-to-site "
                "scatter (analysis/GEN1_DATUM_MODULE.md)")

    sites = G.discover_near_point(control, a.easting, a.northing, 1000.0 * a.radius_km)
    if a.cover:
        want = {c.strip() for c in a.cover.split(",")}
        sites = [s for s in sites if s.mark.cover_class in want]
    claim = G.assert_no_geoid_conversion([s.mark for s in sites])
    R.notes.append(control.merge_note)
    R.notes.append(claim)

    res = G.resolve_tiles(sites, [a.tiles], cache=a.index_cache)
    R.notes.append(f"{len(res.on_disk)} of {len(res.needs)} tiles are on disk; "
                   f"{len(res.to_fetch)} would have to be fetched and NONE was")

    for k, v in G.TileResolution.table_columns().items():
        R.column(k, v)
    mcols = ["point", "cover", "km", "tile", "line", "line_frac", "n_lines", "n",
             "slope_deg", "relief_mm", "fit_rms_mm", "radius_spread_mm", "tie_mm", "sigma_mm"]
    for k, v in G.MarkMeasurement.table_columns().items():
        R.column(k, v)
    for k, v in G.Gen1DatumEstimate.table_columns().items():
        if k not in R.columns:
            R.column(k, v)
    R.column("marks", "the control mark ids in that group")
    R.column("stat", "name of the statistic on this row")
    R.column("value", "its value, units in the name")
    R.banner()

    print(f"\n--- tiles ({len(res.needs)}) ---")
    R.table(list(G.TileResolution.table_columns()), res.table_rows())

    meas, skipped = G.measure_sites(sites, res, swath_constants=swath_const,
                                    swath_constants_source=swath_src, lateral_shift_m=lat)

    print(f"\n--- marks measured ({len(meas)}); skipped {len(skipped)} ---")
    R.table(mcols, [m.table_row(mcols) for m in meas])
    if skipped:
        print("  skipped:")
        for pid, tile, why in skipped:
            print(f"    {pid:32s} {tile or '--':12s} {why}")

    est = G.combine_datum(meas, mode=a.mode, swath_constants_source=swath_src)
    print(f"\n--- by flight line, from the RETURNS ({est.n_lines} lines) ---")
    R.table(list(G.Gen1DatumEstimate.table_columns()), est.table_rows())
    if est.excluded:
        print("  excluded from the combination:")
        for pid, why in est.excluded:
            print(f"    {pid:32s} {why}")

    print("\n--- the estimate ---")
    print(est.summary())

    ties = np.array([m.tie_mm for m in meas], float)
    spread = np.array([m.screen.radius_spread_mm for m in meas], float)
    print("\n--- does the radius-spread screen reduce the site-to-site scatter? ---")
    R.table(["stat", "value"], [
        ["marks measured", len(meas)],
        ["sd of the per-mark ties, mm", f"{ties.std(ddof=1):.1f}" if ties.size > 1 else "--"],
        ["sd of the line means, mm",
         f"{np.std([g.mean_mm for g in est.groups], ddof=1):.1f}" if est.n_lines > 1 else "--"],
    ])
    for cut in (15.0, 25.0, 50.0, 100.0, np.inf):
        k = spread <= cut
        if k.sum() > 1:
            R.table(["stat", "value"], [
                [f"marks with radius_spread <= {cut} mm", int(k.sum())],
                ["  sd of their ties, mm", f"{ties[k].std(ddof=1):.1f}"],
                ["  mean of their ties, mm", f"{ties[k].mean():+.1f}"],
            ])

    if a.out:
        with open(a.out, "w") as f:
            json.dump(dict(estimate=est.to_dict(),
                           marks=[dict(point_id=m.point_id, cover=m.site.mark.cover_class,
                                       km=m.site.distance_m / 1000.0, tile=m.tile,
                                       line=m.line_id, line_counts=m.line.counts,
                                       tie_mm=m.tie_mm, sigma_mm=m.sigma_mm,
                                       n=m.screen.n, slope_deg=m.screen.slope_deg,
                                       relief_mm=m.screen.relief_mm,
                                       fit_rms_mm=m.screen.fit_rms_mm,
                                       radius_spread_mm=m.screen.radius_spread_mm,
                                       notes=m.notes) for m in meas],
                           skipped=[list(s) for s in skipped],
                           tiles_to_fetch=[n.tile for n in res.to_fetch]), f, indent=1)
        print(f"\nwrote {a.out}")

    R.done(headline=f"gen1 datum {est.value_mm:+.1f} +/- {est.se_mm:.1f} mm "
                    f"({est.mode}, {est.n_marks} marks, {est.n_lines} lines)")


if __name__ == "__main__":
    main()
