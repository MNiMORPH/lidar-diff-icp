#!/usr/bin/env python3
"""Tie gen1's floating vertical datum at Elba to surveyed 3DEP checkpoints.

Runs the whole :mod:`lidar_diff_icp.groundtruth` path end to end on tiles already on
disk, downloading nothing:

1. load the surveyed checkpoints with their datum stated (bundled CSV);
2. inventory the gen1 tiles, measure every flight-line overlap;
3. for each checkpoint, test the ALONG-SWATH case first, then plan the shortest chain to
   the elbaext swath network and solve it link by link on the overlaps alone;
4. estimate the lidar ground at the mark over a ladder of radii, with the chain shift and
   the per-point geoid shift applied;
5. report each checkpoint's tie SEPARATELY with its own radius curve, then the three
   independent disagreements: two marks on one line (an estimator check), two links of
   one chain (a link check), and west versus east (a chain check).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/groundtruth/elba_absolute_tie.py [--ground csf|vendor] [--east/--no-east]

Every parameter printed is generated from the estimators' own Param records, not typed
here. Sign convention: ``tie`` is the constant to ADD to gen1 (already in swath 133's
frame and already geoid-shifted to GEOID18) to place it on surveyed NAVD88(GEOID18).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                    # noqa: E402
from lidar_diff_icp.groundtruth import chain as K                   # noqa: E402
from lidar_diff_icp.groundtruth import checkpoints as C             # noqa: E402
from lidar_diff_icp.groundtruth import tie as T                     # noqa: E402
from lidar_diff_icp.groundtruth.provenance import Param, declare    # noqa: E402

# gen1 row-29 tiles across the Elba latitude band -- the corridor the chain is solved in.
# Row 29 keeps every link at the same along-track position as Elba, so the per-swath
# along-track drift term stays common instead of accumulating (ELBAEXT2_SCOPE.md section 3).
WEST_TILES = ["data/before/4342-29-61.laz", "data/before/4342-29-62.laz",
              "data/before/4342-29-63.laz", "data/before/4342-29-64.laz"]
EAST_TILES = ["data/before/4358-29-01.laz", "data/before/4358-29-02.laz",
              "data/before/4358-29-03.laz"]
# Tile holding each checkpoint (may sit outside the corridor row; that is reported).
CP_TILE = {"2210_2021_MN": "data/before/4342-29-61.laz",
           "3056_2021_MN": "data/before/4342-29-61.laz",
           "2024_2021_MN": "data/before/4342-28-61.laz",
           "2036_2021_MN": "data/before/4358-29-03.laz"}
ELBAEXT = "data/derived/elbaext/corrections_geoid.json"
CACHE = "data/derived/groundtruth"

# gen1's own published vertical accuracy for the county containing Elba: per-county RMSEz
# = 0.161 m for Winona (InPort 68818, quoted in ABSOLUTE_ELEVATION_REFS.md section 1b).
# Used ONLY as the tolerance handed to TieEstimate.verdict -- a tie whose radius
# sensitivity exceeds the dataset's own accuracy cannot inform anything.
GEN1_WINONA_RMSE_MM = 161.0
TOL_SOURCE = ("gen1 per-county vertical RMSEz for Winona, InPort 68818, "
              "ABSOLUTE_ELEVATION_REFS.md section 1b")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ground", choices=("csf", "vendor"), default="csf",
                    help="ground source for the tie (repo default: csf)")
    ap.add_argument("--cover-radius", type=float, default=10.0,
                    help="radius within which a flight line counts as covering a mark, m")
    ap.add_argument("--csf-halfwidth", type=float, default=300.0,
                    help="half-width of the crop CSF is run on, m")
    ap.add_argument("--east", dest="east", action="store_true", default=True)
    ap.add_argument("--no-east", dest="east", action="store_false")
    ap.add_argument("--cross-check-ground", action="store_true", default=True,
                    help="also read the tie with the other ground source, to show the "
                         "choice does not move it")
    A = ap.parse_args()

    R = Run("what constant must be added to gen1 at Elba to place it on surveyed "
            "NAVD88(GEOID18), and do independent control points agree?")

    tiles = list(WEST_TILES) + (list(EAST_TILES) if A.east else [])
    for p in tiles:
        R.input(p, role="gen1 2008 MN DNR tile: the corridor the swath chain is solved in")
    for p in sorted(set(CP_TILE.values())):
        if p not in tiles:
            R.input(p, role="gen1 2008 MN DNR tile holding a checkpoint (outside the "
                            "corridor row; its along-track offset is reported)")
    R.input(ELBAEXT, role="elbaext datum + per-swath alignment; swath 133 is the "
                          "reference frame the tie is expressed in")

    cps = C.load_bundled()
    R.input(cps.origin, role="surveyed 3DEP QA checkpoints near Elba, NAVD88(GEOID18)")
    R.param("checkpoint_set", os.path.basename(cps.origin), src="repo",
            why="the six 3DEP vertical-accuracy marks within ~12 km of Elba, transcribed "
                "in ABSOLUTE_ELEVATION_REFS.md section 1a")
    R.param("verdict_tolerance_mm", GEN1_WINONA_RMSE_MM, src="repo", why=TOL_SOURCE)
    R.param("ground_source", A.ground, src="repo",
            why="data/derived/elbaext/corrections_geoid.json ground_source = csf")
    R.param("covering_line_radius_m", A.cover_radius, src="repo",
            why="ELBAEXT2_SCOPE.md section 2 counts covering returns within 10 m of each "
                "mark (173-201 class-2 returns there); it only decides which lines are "
                "candidate chain sources, not which returns enter the tie")
    R.param("csf_crop_halfwidth_m", A.csf_halfwidth, src="MINE",
            why="CSF is run on a crop, not a tile, to keep it seconds instead of minutes; "
                "300 m gives a 600 m box (~350k gen1 points) so the cloth is not dominated "
                "by the crop edge. It excludes nothing from the estimate: the widest "
                "fitting radius is 25 m.")

    elba = json.load(open(ELBAEXT))
    swath_corr = {int(k): v for k, v in
                  elba["per_swath_internal_alignment_dxdydz_m"].items()}
    ref_line = min(swath_corr, key=lambda k: abs(swath_corr[k][2]) + abs(swath_corr[k][0]))
    lat = elba["cross_epoch_datum"]["horizontal_shift_m"]
    target_lines = sorted(swath_corr)
    R.param("target_lines", target_lines, src="repo",
            why="the gen1 swaths of the elbaext product; swath "
                f"{ref_line} is its reference (0,0,0)")
    R.param("elbaext_lateral_shift_m", tuple(lat), src="repo",
            why="corrections_geoid.json cross_epoch_datum horizontal_shift_m -- the "
                "Nuth-Kaaab shift that puts gen1 in gen2's (and the checkpoints') "
                "horizontal frame; reported both applied and not")

    for k, v in T.TieEstimate.table_columns().items():
        R.column(k, v)
    for k, v in K.ChainSolution.table_columns().items():
        R.column(k, v)
    R.column("checkpoint", "surveyed 3DEP QA mark id")
    R.column("type", "NVA = open ground, VVA = under vegetation (3DEP convention)")
    R.column("line", "gen1 flight line (point_source_id) covering the mark")
    R.column("links", "cross-swath links between that line and the elbaext frame, count")
    R.column("n", "ground returns inside the report radius, count")
    R.column("tie_mm", "constant to ADD to gen1 in swath-133 frame to reach the surveyed "
                       "datum, mm")
    R.column("sigma_mm", "half the tie's spread across the pipeline-scale radii, mm")
    R.column("tie_med_mm", "the tie taken as the median over the pipeline-scale radii, mm")
    R.column("chain_mm", "chain-accumulated vertical shift into the elbaext frame, mm")
    R.column("chain_sig", "per-link sigmas of that shift in quadrature, mm")
    R.column("geoid_mm", "references.geoid_difference GEOID03->GEOID18 at the mark, mm")
    R.column("usable", "does the radius spread fall inside the verdict tolerance")
    R.column("pair", "which two results are being compared")
    R.column("what", "what a disagreement between them would mean")
    R.column("diff_mm", "their difference, mm")
    R.column("dN_km", "northing offset of the mark from the corridor band centre, km")
    R.column("swath", "gen1 flight line whose fitted along-track drift curve is summarised")
    R.column("p2p_mm", "peak-to-peak of that fitted drift curve, mm")
    R.column("span_km", "along-track distance the curve covers, km (gps_time x 76 m/s)")
    R.column("mm_per_km", "peak-to-peak divided by span -- the scale of the drift term "
                          "this module does NOT correct, mm/km")

    # ---------------------------------------------------------------- inventory + graph
    inv = K.build_inventory(tiles, cache_dir=CACHE,
                            inventory_json=os.path.join(CACHE, "inventory.json"))
    graph = K.overlap_graph(inv)
    declare(R, [Param("link_res_m", K.DEFAULT_RES, "repo",
                      "coreg.coregister_swaths / align_swaths default grid resolution"),
                Param("link_exclude_classes", K.DEFAULT_EXCLUDE, "repo",
                      "coreg.coregister_swaths terrain proxy ~isin(classification,(5,6,9)) "
                      "-- the VENDOR classification, so no CSF is needed for the chain")])
    R.banner()

    print("\n== flight-line overlap graph (measured, "
          f"{K.DEFAULT_RES:g} m cells) ==")
    R.column("area_km2", "measured overlap area of the pair on the solve grid, km^2")
    print("  " + "  ".join(f"{a}-{b}:{g['area_km2']:.2f}" for (a, b), g in sorted(graph.items())))
    print(f"  {len(graph)} adjacent pairs over {len(inv.lines)} lines "
          f"{inv.lines[0]}-{inv.lines[-1]}; no non-adjacent pair shares cells")

    # How big is the term this module cannot apply? The pipeline fits a per-swath
    # along-track drift spline against the gen2 grid, so it exists only where gen2 does --
    # i.e. over elbaext, not out at the checkpoints. Its measured size sets how far a
    # checkpoint may sit from the corridor band before the tie stops meaning anything.
    print("\n== the term this module does NOT apply: per-swath along-track drift ==")
    drows = []
    for sw, cur in sorted(elba["along_track_drift_gpsTime_to_m"].items()):
        t = np.asarray(cur["gps_time"]); z = np.asarray(cur["drift_m"])
        p2p = 1000.0 * (z.max() - z.min())
        span = (t.max() - t.min()) * 76.0 / 1000.0        # 76 m/s: elbaext grid / gps span
        drows.append([sw, f"{p2p:.1f}", f"{span:.2f}", f"{p2p/span:.0f}"])
    R.table(["swath", "p2p_mm", "span_km", "mm_per_km"], drows)
    drift_scale = float(np.median([float(r[3]) for r in drows]))
    print(f"  median {drift_scale:.0f} mm/km. pipeline.fit_along_track_drift regresses "
          "against the gen2 grid (pipeline.py:718), which does not exist at the "
          "checkpoints, so this term is UNMODELLED there.")

    band = [np.mean([t.lines[s]["bbox"][1] + t.lines[s]["bbox"][3]
                     for s in t.lines]) / 2.0 for t in inv.tiles.values()]
    band_centre = float(np.mean(band))

    # ------------------------------------------------------------------ per checkpoint
    results = []
    for cp in cps.usable():
        if cp.point_id not in CP_TILE:
            print(f"\n== {cp.point_id} ==\n  no tile on disk covers this mark "
                  f"(ELBAEXT2_SCOPE section 8: 2099 needs 4342-26-61, 3089 needs "
                  f"4358-26-02). NOT ATTEMPTED -- reported, not dropped.")
            results.append(dict(cp=cp, tie=None, sol=None, reason="no tile on disk"))
            continue
        print(f"\n{'='*78}\n== {cp.point_id}  ({cp.point_type})  surveyed "
              f"{cp.elevation_m:.3f} m {cp.vertical_datum}/{cp.geoid_model} ==")
        tile = CP_TILE[cp.point_id]

        # Step 1, always: which lines actually put returns on this mark?
        cov = K.covering_lines(K.build_inventory([tile], cache_dir=CACHE,
                                                 inventory_json=os.path.join(
                                                     CACHE, "cp_inventory.json")),
                               cp.easting, cp.northing, A.cover_radius)
        print(f"  covering lines within {A.cover_radius:g} m: " +
              ", ".join(f"{k} ({v} returns)" for k, v in cov.items()))
        paths = K.plan_path(graph, inv, source_lines=list(cov), target_lines=target_lines)
        if not paths:
            print("  NO PATH from a covering line to the elbaext frame in this corridor. "
                  "Reported, not dropped.")
            results.append(dict(cp=cp, tie=None, sol=None, reason="no path"))
            continue
        p = paths[0]
        if p.along_swath:
            print(f"  ALONG-SWATH: line {p.nodes[0]} is itself an elbaext swath -- "
                  "0 links, no cross-swath transfer.")
            sol = None
            dz_chain = dx_chain = dy_chain = 0.0
            sig_chain = 0.0
            far = p.nodes[0]
        else:
            print(f"  along-swath test FAILS (no covering line is an elbaext swath); "
                  f"shortest chain: {p}")
            if len(paths) > 1:
                print(f"  {len(paths)} equally short routes: " +
                      "; ".join(str(q) for q in paths))
            sol = K.solve_chain(inv, p, graph=graph, reference="last")
            print(f"  links solved on their overlaps only, in tiles: " +
                  ", ".join(sorted({os.path.basename(t) for L in sol.links for t in L.tiles})))
            R.table(list(K.ChainSolution.table_columns()), sol.table_rows())
            dz_chain, dx_chain, dy_chain = sol.dz_total_m, sol.dx_total_m, sol.dy_total_m
            sig_chain = sol.dz_sigma_m
            far = sol.far_line
            tgt = sol.reference_line
            cdx, cdy, cdz = swath_corr[tgt]
            if any(abs(v) > 0 for v in (cdx, cdy, cdz)):
                print(f"  + elbaext internal alignment of swath {tgt} into swath "
                      f"{ref_line}'s frame: ({cdx:+.4f}, {cdy:+.4f}, {cdz:+.4f}) m")
            dz_chain += cdz; dx_chain += cdx; dy_chain += cdy

        geoid = T.geoid_shift_for(cp)
        loader = (T.csf_ground_near if A.ground == "csf" else T.vendor_ground_near)
        kw = dict(cache_dir=os.path.join(CACHE, "csf")) if A.ground == "csf" else {}
        g = loader(tile, cp.easting, cp.northing, A.csf_halfwidth, **kw)
        print(f"  ground: {g.source}, {len(g):,} returns from {g.n_input:,} in the "
              f"{2*A.csf_halfwidth:.0f} m crop")

        est = T.estimate_tie(cp, g, line=far, geoid_shift_m=geoid,
                             swath_shift_m=(dx_chain + lat[0], dy_chain + lat[1], dz_chain))
        no_lat = T.estimate_tie(cp, g, line=far, geoid_shift_m=geoid,
                                swath_shift_m=(dx_chain, dy_chain, dz_chain))
        R.table(list(T.TieEstimate.table_columns()), est.table_rows())
        ok, why = est.verdict(GEN1_WINONA_RMSE_MM, tolerance_source=TOL_SOURCE)
        print(f"  TIE = {est.tie_mm:+.1f} +/- {est.sigma_mm:.1f} mm   "
              f"(n={est.n_report} at R={est.report_radius_m:g} m)")
        print(f"  median over the pipeline-scale radii: {est.tie_median_mm:+.1f} mm "
              f"({est.tie_median_mm - est.tie_mm:+.1f} mm from the R="
              f"{est.report_radius_m:g} m read)")
        print(f"  radius spread {est.radius_spread_mm:.1f} mm over the pipeline-scale "
              f"radii, {est.radius_spread_all_mm:.1f} mm over the whole ladder; "
              f"fit SE {est.fit_se_mm:.1f} mm (optimistic -- returns 1 m apart are not "
              f"independent)")
        print(f"  usable: {ok} -- {why}")
        print(f"  lateral term: with the elbaext Nuth-Kaaab shift {est.tie_mm:+.1f} mm, "
              f"without it {no_lat.tie_mm:+.1f} mm (difference "
              f"{est.tie_mm - no_lat.tie_mm:+.1f} mm)")
        for n in est.notes:
            print(f"  note: {n}")
        dN = (cp.northing - band_centre) / 1000.0
        if abs(dN) > 1.0:
            print(f"  CAVEAT: this mark sits {dN:+.2f} km from the corridor band the links "
                  "were solved in. At the measured drift scale above that is "
                  f"~{abs(dN) * drift_scale:.0f} mm of UNMODELLED along-track drift -- "
                  "the tie here carries that on top of everything in the table.")
        if A.cross_check_ground:
            other = ("vendor" if A.ground == "csf" else "csf")
            ldr = (T.vendor_ground_near if other == "vendor" else T.csf_ground_near)
            kw2 = {} if other == "vendor" else dict(cache_dir=os.path.join(CACHE, "csf"))
            g2 = ldr(tile, cp.easting, cp.northing, A.csf_halfwidth, **kw2)
            e2 = T.estimate_tie(cp, g2, line=far, geoid_shift_m=geoid,
                                swath_shift_m=(dx_chain + lat[0], dy_chain + lat[1], dz_chain))
            print(f"  ground-source check: {other} gives {e2.tie_mm:+.1f} mm "
                  f"({e2.tie_mm - est.tie_mm:+.1f} mm from {A.ground})")
        results.append(dict(cp=cp, tie=est, sol=sol, drift_mm=abs(dN) * drift_scale,
                            chain_mm=dz_chain * 1000.0,
                            chain_sig=sig_chain * 1000.0, geoid_mm=geoid * 1000.0,
                            line=far, links=p.n_links, ok=ok, dN=dN, reason=""))
        inv.release()

    # --------------------------------------------------------------------- the summary
    print(f"\n{'='*78}\n== every control point, separately ==")
    rows = []
    for r in results:
        cp = r["cp"]
        if r["tie"] is None:
            rows.append([cp.point_id, cp.point_type, "-", "-", "-", "-", "-", "-", "-",
                         "-", f"NOT ATTEMPTED: {r['reason']}"])
            continue
        rows.append([cp.point_id, cp.point_type, r["line"], r["links"], r["tie"].n_report,
                     f"{r['tie'].tie_mm:+.1f}", f"{r['tie'].sigma_mm:.1f}",
                     f"{r['tie'].tie_median_mm:+.1f}",
                     f"{r['chain_mm']:+.1f}", f"{r['geoid_mm']:+.1f}", str(r["ok"])])
    R.column("note", "why a checkpoint has no tie, where it has none")
    head = ["checkpoint", "type", "line", "links", "n", "tie_mm", "sigma_mm",
            "tie_med_mm", "chain_mm", "geoid_mm"]
    R.table(head + (["usable"] if all(r["tie"] is not None for r in results) else ["note"]),
            rows)

    print("\n== the three independent disagreements ==")
    got = {r["cp"].point_id: r for r in results if r["tie"] is not None}
    cmp_rows = []

    def add(a, b, what):
        if a in got and b in got:
            cmp_rows.append([f"{a} vs {b}", what,
                             f"{got[a]['tie'].tie_mm - got[b]['tie'].tie_mm:+.1f}"])

    add("2210_2021_MN", "3056_2021_MN",
        "two marks on ONE line: a disagreement is the ESTIMATOR, not the chain")
    add("2210_2021_MN", "2024_2021_MN",
        "two links of one chain: a disagreement is the LINK 128-129")
    add("2210_2021_MN", "2036_2021_MN",
        "west chain vs east chain: a disagreement is the CHAIN (the only real check)")
    if cmp_rows:
        R.table(["pair", "what", "diff_mm"], cmp_rows)

    # Two strata, both declared by the DATA and the METHOD before any tie was computed,
    # not chosen after seeing the answers:
    #   NVA vs VVA  -- the survey's own accuracy class (3DEP: 3.5 cm RMSE vs 27 cm @95%);
    #   on-band vs off-band -- whether the along-track drift term above is unmodelled.
    # Every checkpoint is reported in the full table above regardless of stratum.
    print("\n== strata (both declared before the ties were computed) ==")
    have = [r for r in results if r["tie"] is not None]
    for label, sel in (("all control points", have),
                       (f"on the corridor band (|dN| <= 1 km)",
                        [r for r in have if abs(r["dN"]) <= 1.0]),
                       ("on band AND NVA (open ground)",
                        [r for r in have if abs(r["dN"]) <= 1.0
                         and r["cp"].point_type.upper() == "NVA"])):
        if len(sel) < 2:
            continue
        v = [r["tie"].tie_mm for r in sel]
        print(f"  {label:34s} n={len(sel)}  spread {max(v) - min(v):6.1f} mm  "
              f"mean {np.mean(v):+7.1f} mm   [" +
              ", ".join(r["cp"].point_id.split("_")[0] for r in sel) + "]")
    for r in have:
        if abs(r["dN"]) > 1.0:
            print(f"  {r['cp'].point_id} sits {r['dN']:+.2f} km off band. At the measured "
                  f"{drift_scale:.0f} mm/km that is ~{r['drift_mm']:.0f} mm of unmodelled "
                  "drift -- PART of, not all of, its departure from the on-band ties.")

    ties = [r["tie"].tie_mm for r in results if r["tie"] is not None]
    if len(ties) > 1:
        print(f"\n  spread across {len(ties)} independent control points: "
              f"{max(ties) - min(ties):.1f} mm  (median {np.median(ties):+.1f} mm)")
        print("  The spread ACROSS MARKS is the error bar. It is not averaged away above,")
        print("  because a chain has no internal redundancy and its formal sigma "
              f"({np.hypot(*[r['chain_sig'] for r in results if r['tie'] is not None][:2]):.1f} "
              "mm scale) cannot see accumulated error.")
    R.done(headline=(f"gen1 tie at Elba from {len(ties)} control point(s): "
                     + ", ".join(f"{r['cp'].point_id}={r['tie'].tie_mm:+.0f}mm"
                                 for r in results if r["tie"] is not None)))


if __name__ == "__main__":
    main()
