"""gen1 against its OWN 2008 ground control, from the point cloud.

Question
--------
The Elba datum constant, ``+22.7 +/- 39.7 mm`` (``analysis/ABSOLUTE_BASIS_ELBA.md``), is
carried in from two 2021 3DEP marks by five- and six-link swath chains, across a
GEOID03 -> GEOID18 conversion.  gen1 has its own 2008 control -- same epoch, same geoid,
directly under gen1's flight lines -- and the MnDNR validation tables say the delivered
2008 *surface* sits ABOVE that control.  Opposite sign.  This script measures **our own**
residual at those marks, from the point cloud, with the estimator the anchor used, so the
two numbers can be compared without a proxy anywhere in between.

Why the vendor's ``Surface Z`` column is not used as the answer: it is the delivered 2008
DNR DEM, a different surface from our slope-normal reconstruction, and it may already
carry the (unpublished) bias adjustment of metadata process step 8.  It is reported
beside our number, never in place of it.

Sign convention -- one family throughout, and it is the SAME one on both sides:

    tie = surveyed - z_lidar     (groundtruth/tie.py; docs/groundtruth.md section 2)
    MnDNR Error = Control Z - Surface Z   (pinned by arithmetic on 1022/1022 rows,
                                           parse_mndnr_2008_control.py --check)

so **positive means the lidar reads BELOW the mark** in both.  The 2008 control is on
GEOID03 and so is the raw gen1 cloud, so ``geoid_shift_m = 0`` here is not an
approximation: it is the same frame on both sides.  And because the anchor's tie
subtracts the geoid shift it adds, ``tie21`` and ``tie08`` are the same quantity
(``H03_surveyed - H03_gen1``) -- the geoid cancels out of the comparison.

Usage
-----
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/groundtruth/gen1_own_control_tie.py [--ground csf|vendor]
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                        # noqa: E402
from lidar_diff_icp.groundtruth import checkpoints as C                 # noqa: E402
from lidar_diff_icp.groundtruth import tie as T                         # noqa: E402

REF = (579705.72, 4883677.71)            # the Elba reference point
ELBAEXT = "data/derived/elbaext/corrections_geoid.json"
CACHE = "data/derived/groundtruth"
CONTROL = "mn_dnr_2008_control_semn"
# The 2021 anchor this script exists to ARGUE AGAINST. It is read deliberately, not by
# oversight: the comparison IS the point, so it must keep reading the anchor's own numbers.
# That anchor was subsequently RETRACTED (z_before_absolute.RETRACTED.md, 2026-08-28) and
# the conclusion below is what replaced it -- gen1 measured against its OWN 2008 control,
# adopted in ground_control/products/ANSWER_gen1_elba.json (+62.74 +/- 23.38 mm delivered).
# Do NOT repoint this to the adopted product: it would compare that product with itself.
SIDECAR = "data/derived/elba_fulldensity/z_before_absolute.json"
ADOPTED = "ground_control/products/ANSWER_gen1_elba.json"
COVER_NAME = {"L1O": "open terrain", "L2T": "tall weeds/crops", "L3B": "brush/low trees",
              "L4F": "forested", "L5U": "urban", "other": "unclassed"}


def tile_bounds(paths):
    import laspy
    out = {}
    for p in paths:
        with laspy.open(p) as f:
            h = f.header
            out[p] = (h.mins[0], h.mins[1], h.maxs[0], h.maxs[1])
    return out


def siting(g, cp, radius):
    """Section 7.1 screen: local relief and where the surveyed value falls in it."""
    d = np.hypot(g.x - cp.easting, g.y - cp.northing)
    z = g.z[d <= radius]
    if z.size < 5:
        return np.nan, np.nan, int(z.size)
    p05, p95 = np.percentile(z, [5, 95])
    pct = 100.0 * float((z < cp.elevation_m).mean())
    return float(p95 - p05), pct, int(z.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ground", choices=("csf", "vendor"), default="csf")
    ap.add_argument("--csf-halfwidth", type=float, default=300.0)
    ap.add_argument("--siting-radius", type=float, default=5.0)
    ap.add_argument("--json", default=None)
    A = ap.parse_args()

    R = Run("does gen1's own 2008 ground control agree with the +22.7 mm datum anchor "
            "carried in from the 2021 3DEP marks, and with what sign?")

    tiles = sorted(p for p in glob.glob("data/before/*.laz") if "merged" not in p)
    for p in tiles:
        R.input(p, role="gen1 2008 MN DNR tile already on disk")
    cps = C.load_bundled(CONTROL)
    R.input(str(C._DATA / f"{CONTROL}.csv"),
            role="MnDNR 2008 control, 8 SE-MN counties, parsed from the MnGeo county "
                 "validation reports (NAVD88/GEOID03 -- gen1's own geoid)")
    R.input(ELBAEXT, role="elbaext datum + lateral shift (the Nuth-Kaaeb term)")
    marks21 = C.load_bundled()
    R.input(str(C._DATA / "mn_se_driftless_2021_ql1_near_elba.csv"),
            role="the 2021 3DEP marks the +22.7 mm anchor is built on, for a "
                 "same-estimator siting comparison")

    elba = json.load(open(ELBAEXT))
    lat = tuple(float(v) for v in elba["cross_epoch_datum"]["horizontal_shift_m"])
    res = float(elba["res_m"])

    R.param("ground_source", A.ground, src="repo",
            why=f"{ELBAEXT} ground_source; the same source the anchor's ties were read on")
    R.param("res_m", res, src="repo", why=f"{ELBAEXT} res_m")
    R.param("csf_crop_halfwidth_m", A.csf_halfwidth, src="MINE",
            why="CSF is run on a crop, not a tile, to keep it seconds instead of minutes. "
                "Copied unchanged from analysis/groundtruth/elba_absolute_tie.py so the "
                "two runs are the same method; its effect is measured below by re-reading "
                "every mark with the vendor class-2 ground, which uses no crop at all")
    R.param("siting_radius_m", A.siting_radius, src="repo",
            why="ADDITIONAL_GROUND_CONTROL.md section 7.1 states the screen at R = 5 m; "
                "kept at 5 m so this run is comparable to the 17 measurements there")
    R.param("geoid_shift_m", 0.0, src="repo",
            why="the 2008 control and the raw gen1 cloud are both NAVD88(GEOID03) "
                "(lidar_semn2008.html); no conversion enters this comparison")
    R.param("lateral_shift_m", lat, src="repo",
            why=f"{ELBAEXT} cross_epoch_datum.horizontal_shift_m: the elbaext Nuth-Kaaeb shift, reported applied AND "
                "withheld because the 2008 control is in gen1's own horizontal frame")
    R.param("swath_chain", "none", src="repo",
            why="every mark here lies under a gen1 flight line in a tile on disk, so no "
                "chain link is needed -- which is the whole point of this control set")

    R.column("tie_mm", "surveyed - z_lidar, mm; POSITIVE = gen1 reads BELOW the mark")
    R.column("dnr_mm", "the MnDNR table's own Error = Control Z - delivered 2008 DEM, mm; "
                       "same sign family, different surface")
    R.column("spread_mm", "p05-p95 of ground returns within the siting radius (section 7.1)")
    R.column("ctrl_pct", "percentile of the surveyed elevation within those returns; "
                         "near 50 is well sited, a tail is not")
    R.column("sigma_mm", "half the tie spread over the pipeline-scale radii (res/2..2*res)")
    R.column("point", "MnDNR 2008 control point name; the prefix is its land-cover class")
    R.column("cover", "land-cover class from the point name: L1O open, L2T tall "
                      "weeds/crops, L3B brush/low trees, L4F forested, L5U urban")
    R.column("what", "what the land-cover class means, from lidar_semn2008.html")
    R.column("km", "distance from the Elba reference point (579705.72, 4883677.71), km")
    R.column("tile", "the gen1 tile the mark falls in, already on disk")
    R.column("n", "ground returns inside the report radius (or the siting radius, in the "
                  "siting tables)")
    R.column("lines", "distinct flight lines contributing inside the report radius")
    R.column("tie_nolat_mm", "the same tie with the elbaext lateral shift WITHHELD")
    R.column("mean_tie_mm", "mean of tie_mm over the marks of this class")
    R.column("median_tie_mm", "median of tie_mm over the marks of this class")
    R.column("SE_mm", "standard error of the mean over the marks in the row")
    R.column("mean_spread_mm", "mean of spread_mm over the marks of this class")
    R.column("mean_dnr_mm", "mean of dnr_mm over the marks of this class")
    R.column("quantity", "which estimate of the gen1 datum constant this row holds")
    R.column("mm", "its value in mm, on the tie = surveyed - z_lidar convention")
    R.column("note", "how the row was obtained, or what qualifies it")
    R.column("line", "gen1 flight line (LAS point_source_id) the tie was read on alone")
    R.column("mixed_tie_mm", "the tie for this mark with all lines pooled, mm, for "
                             "comparison with the per-line reads")
    R.column("mark", "2021 3DEP vertical-accuracy mark id")
    R.column("type", "3DEP accuracy class: NVA open ground, VVA under vegetation")
    R.column("mean_mm", "mean of the MnDNR Error over the marks in the row, mm")
    R.column("median_mm", "median of the MnDNR Error over the marks in the row, mm")
    R.column("L1O_mean_mm", "mean MnDNR Error over the OPEN-cover marks of the row only")
    R.column("n_L1O", "how many of the row's marks are open cover")
    R.column("stat", "which statistic of (our point-cloud tie - the MnDNR table's own "
                     "Error) is on this row, with its units")
    R.column("value", "the value of the statistic named in the 'stat' column, in the "
                      "units that column states (mm, or dimensionless for a count or a "
                      "correlation)")
    R.column("stratum", "which marks the plane was fitted to")
    R.column("radius_km", "marks within this distance of the Elba reference point")
    R.column("intercept_mm", "the fitted plane's value AT the Elba reference point, mm")
    R.column("dE_mm_per_km", "east gradient of the fitted plane, mm/km")
    R.column("dN_mm_per_km", "north gradient of the fitted plane, mm/km")
    R.column("resid_rms_mm", "RMS of the marks about the fitted plane, mm")
    R.column("county", "MnGeo validation report the mark came from")
    for _r in T.radius_ladder(res):
        R.column(f"R={_r:g}", f"tie_mm read at radius {_r:g} m")
    R.banner()

    bounds = tile_bounds(tiles)
    loader = T.csf_ground_near if A.ground == "csf" else T.vendor_ground_near
    kw = dict(cache_dir=os.path.join(CACHE, "csf")) if A.ground == "csf" else {}

    # ---------------------------------------------------------------- the 2008 control
    todo, seen_xy = [], set()
    for cp in cps:
        key = (round(cp.easting, 3), round(cp.northing, 3))
        if key in seen_xy:
            continue
        seen_xy.add(key)
        for p, (x0, y0, x1, y1) in bounds.items():
            if x0 <= cp.easting <= x1 and y0 <= cp.northing <= y1:
                todo.append((cp, p))
                break
    todo.sort(key=lambda t: np.hypot(t[0].easting - REF[0], t[0].northing - REF[1]))
    print(f"\n{len(todo)} of {len(cps)} MnDNR 2008 control points fall inside the "
          f"{len(tiles)} gen1 tiles already on disk. No tile was downloaded.\n")

    rows, recs = [], []
    ladder = None
    for cp, tile in todo:
        g = loader(tile, cp.easting, cp.northing, A.csf_halfwidth, **kw)
        est = T.estimate_tie(cp, g, line=None, res=res, geoid_shift_m=0.0,
                             swath_shift_m=(lat[0], lat[1], 0.0))
        raw = T.estimate_tie(cp, g, line=None, res=res, geoid_shift_m=0.0)
        spread, pct, n5 = siting(g, cp, A.siting_radius)
        d = 1e-3 * float(np.hypot(cp.easting - REF[0], cp.northing - REF[1]))
        rep = next(r for r in est.curve if r.radius_m == est.report_radius_m)
        ladder = [r.radius_m for r in est.curve]
        recs.append(dict(cp=cp, est=est, raw=raw, spread=spread, pct=pct, n5=n5,
                         dist=d, tile=os.path.basename(tile), tile_path=tile,
                         n_lines=rep.n_lines))
        rows.append([cp.point_id, cp.point_type, f"{d:.2f}", os.path.basename(tile)[:10],
                     rep.n, rep.n_lines, f"{est.tie_mm:+.1f}", f"{est.sigma_mm:.1f}",
                     f"{raw.tie_mm:+.1f}", f"{1000 * cp_err(cp, cps):+.1f}",
                     f"{1000 * spread:.0f}", f"{pct:.0f}"])
        del g

    R.table(["point", "cover", "km", "tile", "n", "lines", "tie_mm", "sigma_mm",
             "tie_nolat_mm", "dnr_mm", "spread_mm", "ctrl_pct"], rows)

    print("\nFull radius ladder, tie in mm at each radius "
          f"(lateral shift applied; report radius {1.5 * res:g} m):")
    R.table(["point", "cover"] + [f"R={r:g}" for r in ladder],
            [[r["cp"].point_id, r["cp"].point_type] +
             [(f"{1000 * (r['cp'].elevation_m - e.z_lidar_m):+.0f}" if e.ok else "--")
              for e in r["est"].curve] for r in recs])

    # ------------------------------------------------------------------- by land cover
    print("\nBy land cover (our point-cloud tie, and the MnDNR DEM error beside it):")
    cov_rows = []
    for cvr in ("L1O", "L2T", "L3B", "L4F", "L5U"):
        sel = [r for r in recs if r["cp"].point_type == cvr]
        if not sel:
            continue
        t = np.array([r["est"].tie_mm for r in sel])
        sp = np.array([1000 * r["spread"] for r in sel])
        dn = np.array([1000 * cp_err(r["cp"], cps) for r in sel])
        cov_rows.append([cvr, COVER_NAME[cvr], len(sel), f"{t.mean():+.1f}",
                         f"{np.median(t):+.1f}",
                         f"{t.std(ddof=1) / np.sqrt(len(t)):.1f}" if len(t) > 1 else "--",
                         f"{sp.mean():.0f}", f"{dn.mean():+.1f}"])
    R.table(["cover", "what", "n", "mean_tie_mm", "median_tie_mm", "SE_mm",
             "mean_spread_mm", "mean_dnr_mm"], cov_rows)

    # --------------------------------------- the 2021 anchor marks, SAME estimator
    print("\nThe two anchor marks, read with this same estimator and NO chain, for the "
          "siting comparison (their tie is not meaningful without the chain):")
    a_rows = []
    for cp in marks21:
        hit = [p for p, (x0, y0, x1, y1) in bounds.items()
               if x0 <= cp.easting <= x1 and y0 <= cp.northing <= y1]
        if not hit:
            a_rows.append([cp.point_id, cp.point_type, "--", "--", "--", "--",
                           "no tile on disk"])
            continue
        g = loader(hit[0], cp.easting, cp.northing, A.csf_halfwidth, **kw)
        spread, pct, n5 = siting(g, cp, A.siting_radius)
        est = T.estimate_tie(cp, g, line=None, res=res, geoid_shift_m=0.0)
        a_rows.append([cp.point_id, cp.point_type, os.path.basename(hit[0])[:10], n5,
                       f"{1000 * spread:.0f}", f"{pct:.0f}", f"{est.sigma_mm:.1f}"])
        del g
    R.table(["mark", "type", "tile", "n", "spread_mm", "ctrl_pct", "sigma_mm"], a_rows)

    # ------------------------------------------------------------------ the whole set
    print("\nThe MnDNR table's own residual over all 8 SE-MN counties, by distance from "
          "Elba (Control Z - delivered 2008 DEM; POSITIVE = the DEM reads below control):")
    x = np.array([p.easting for p in cps]); y = np.array([p.northing for p in cps])
    err = np.array([cp_err(p, cps) for p in cps])
    cvr = np.array([p.point_type for p in cps])
    d = np.hypot(x - REF[0], y - REF[1]) / 1000.0
    band_rows = []
    for lo, hi in ((0, 5), (5, 10), (10, 20), (20, 40), (40, 200)):
        m = (d >= lo) & (d < hi)
        if not m.any():
            continue
        band_rows.append([f"{lo}-{hi}", int(m.sum()), f"{1000 * err[m].mean():+.1f}",
                          f"{1000 * np.median(err[m]):+.1f}",
                          f"{1000 * err[m].std(ddof=1) / np.sqrt(m.sum()):.1f}",
                          f"{1000 * (err[m & (cvr == 'L1O')].mean()):+.1f}"
                          if (m & (cvr == "L1O")).any() else "--",
                          int((m & (cvr == "L1O")).sum())])
    R.table(["km", "n", "mean_mm", "median_mm", "SE_mm", "L1O_mean_mm", "n_L1O"], band_rows)

    print("\nAnd by land cover over the whole 8-county set (the physical check on the "
          "sign: vegetation must make the SURFACE read HIGH, i.e. Control - Surface < 0):")
    R.table(["cover", "what", "n", "mean_mm", "median_mm", "SE_mm"],
            [[c, COVER_NAME[c], int((cvr == c).sum()),
              f"{1000 * err[cvr == c].mean():+.1f}",
              f"{1000 * np.median(err[cvr == c]):+.1f}",
              f"{1000 * err[cvr == c].std(ddof=1) / np.sqrt((cvr == c).sum()):.1f}"]
             for c in ("L1O", "L2T", "L3B", "L4F", "L5U", "other")
             if (cvr == c).any()])

    print("\nMarks with more than one gen1 flight line inside the report radius, read "
          "one line at a time -- mixing lines folds their relative offsets into the tie:")
    per_line = []
    for r in recs:
        if r["n_lines"] < 2:
            continue
        cp = r["cp"]
        g = loader(r["tile_path"], cp.easting, cp.northing, A.csf_halfwidth, **kw)
        d_ = np.hypot(g.x - cp.easting, g.y - cp.northing)
        for ln in sorted(set(int(v) for v in g.point_source_id[d_ <= 1.5 * res])):
            e_ = T.estimate_tie(cp, g, line=ln, res=res, geoid_shift_m=0.0,
                                swath_shift_m=(lat[0], lat[1], 0.0))
            rp = next(q for q in e_.curve if q.radius_m == e_.report_radius_m)
            per_line.append([cp.point_id, cp.point_type, ln, rp.n,
                             f"{e_.tie_mm:+.1f}", f"{e_.sigma_mm:.1f}",
                             f"{r['est'].tie_mm:+.1f}"])
        del g
    R.table(["point", "cover", "line", "n", "tie_mm", "sigma_mm", "mixed_tie_mm"], per_line)

    print("\nDoes a mark sitting high in its own local return distribution manufacture a "
          "positive tie? tie_mm against ctrl_pct over the 16 control marks:")
    pc = np.array([r["pct"] for r in recs], float)
    tm = np.array([r["est"].tie_mm for r in recs], float)
    k = np.isfinite(pc) & np.isfinite(tm)
    sl, ic = np.polyfit(pc[k], tm[k], 1)
    rr = np.corrcoef(pc[k], tm[k])[0, 1]
    res_ = tm[k] - (sl * pc[k] + ic)
    se_sl = (np.sqrt((res_ @ res_) / (k.sum() - 2) / ((pc[k] - pc[k].mean()) ** 2).sum()))
    R.table(["stat", "value"], [
        ["n marks", int(k.sum())],
        ["Pearson r (ctrl_pct vs tie_mm)", f"{rr:+.3f}"],
        ["slope, mm per percentile", f"{sl:+.2f} +/- {se_sl:.2f}"],
        ["fitted tie at p50, mm", f"{sl * 50 + ic:+.1f}"],
        ["fitted tie at p93 (mark 2210's siting), mm", f"{sl * 93 + ic:+.1f}"],
        ["fitted tie at p90 (mark 2036's siting), mm", f"{sl * 90 + ic:+.1f}"],
    ])

    print("\nOur point-cloud tie against the MnDNR table's own DEM error, at the same "
          "marks -- does the vendor Surface Z carry a bias our reconstruction does not?")
    ours = np.array([r["est"].tie_mm for r in recs])
    theirs = np.array([1000 * cp_err(r["cp"], cps) for r in recs])
    keep = np.isfinite(ours) & np.isfinite(theirs)
    R.table(["stat", "value"], [
        ["n marks", int(keep.sum())],
        ["median(ours - theirs), mm", f"{np.median(ours[keep] - theirs[keep]):+.1f}"],
        ["mean(ours - theirs), mm", f"{(ours[keep] - theirs[keep]).mean():+.1f}"],
        ["RMS(ours - theirs), mm", f"{np.sqrt(np.mean((ours[keep] - theirs[keep]) ** 2)):.1f}"],
        ["Pearson r", f"{np.corrcoef(ours[keep], theirs[keep])[0, 1]:+.3f}"],
    ])

    print("\nPlane fits to the MnDNR residual about the Elba reference point. The pooled "
          "fit is the one that gave the '2008 surface sits ~22 mm above its control' "
          "reading; the open-cover fit is the stratum an NVA mark belongs to:")
    plane_rows = []
    for label, m0 in (("all covers", np.ones(len(cps), bool)), ("L1O open only", cvr == "L1O")):
        for rad in (10, 20, 30, 50):
            m = m0 & (d < rad)
            if m.sum() < 6:
                plane_rows.append([label, rad, int(m.sum()), "--", "--", "--", "--"])
                continue
            M = np.c_[np.ones(m.sum()), (x[m] - REF[0]) / 1e3, (y[m] - REF[1]) / 1e3]
            c, *_ = np.linalg.lstsq(M, 1000 * err[m], rcond=None)
            r_ = 1000 * err[m] - M @ c
            cov = (r_ @ r_) / (m.sum() - 3) * np.linalg.inv(M.T @ M)
            se = np.sqrt(np.diag(cov))
            plane_rows.append([label, rad, int(m.sum()),
                               f"{c[0]:+.1f} +/- {se[0]:.1f}", f"{c[1]:+.2f} +/- {se[1]:.2f}",
                               f"{c[2]:+.2f} +/- {se[2]:.2f}", f"{np.sqrt((r_ ** 2).mean()):.0f}"])
    R.table(["stratum", "radius_km", "n", "intercept_mm", "dE_mm_per_km", "dN_mm_per_km",
             "resid_rms_mm"], plane_rows)

    print("\nPer-county mean of the MnDNR residual on OPEN cover -- gen1's delivered "
          "surface does not sit at one level across the project:")
    cty = np.array([_county_of(p_.point_id, cps) for p_ in cps])
    R.table(["county", "n", "mean_mm", "median_mm", "SE_mm", "km"],
            [[c, int(((cty == c) & (cvr == "L1O")).sum()),
              f"{1000 * err[(cty == c) & (cvr == 'L1O')].mean():+.1f}",
              f"{1000 * np.median(err[(cty == c) & (cvr == 'L1O')]):+.1f}",
              f"{1000 * err[(cty == c) & (cvr == 'L1O')].std(ddof=1) / np.sqrt(((cty == c) & (cvr == 'L1O')).sum()):.1f}",
              f"{d[(cty == c) & (cvr == 'L1O')].mean():.0f}"]
             for c in sorted(set(cty)) if ((cty == c) & (cvr == "L1O")).sum() > 2])

    print("\nReconciliation with the +22.7 mm anchor. Both are tie = surveyed - z_lidar, "
          "so both are POSITIVE when gen1 reads below the mark, and the geoid cancels "
          "out of the comparison (the anchor subtracts the shift it adds):")
    side = json.load(open(SIDECAR))
    R.input(SIDECAR, role="the +22.7 mm 2021-anchor constant, as shipped and since RETRACTED "
                          "-- read deliberately, as the thing this script compares against")
    print(f"    NOTE: the +22.7 mm anchor was RETRACTED on 2026-08-28. This comparison is "
          f"the argument that retired it; the adopted answer now lives in {ADOPTED}.")
    o = [r for r in recs if r["cp"].point_type == "L1O"]
    t = np.array([r["est"].tie_mm for r in o])
    sg = np.array([r["est"].sigma_mm for r in o])
    w = 1.0 / sg ** 2
    ivw = float((w * t).sum() / w.sum())
    rec_rows = [[f"2008 control, {r['cp'].point_id} ({r['dist']:.2f} km, open)",
                 f"{r['est'].tie_mm:+.1f} +/- {r['est'].sigma_mm:.1f}",
                 f"spread {1000 * r['spread']:.0f} mm at p{r['pct']:.0f}"] for r in o]
    rec_rows += [
        ["2008 control, open marks, inverse-variance weighted", f"{ivw:+.1f}",
         "weights are each mark's own radius sigma, as combine_ties does"],
        ["2008 control, open marks, plain mean", f"{t.mean():+.1f}",
         f"SE of the mean over {len(t)} marks {t.std(ddof=1) / np.sqrt(len(t)):.1f} mm"],
        ["2008 control, open marks, median", f"{np.median(t):+.1f}", "4 marks"],
    ]
    for ti in side["ties"]:
        rec_rows.append([f"2021 anchor, {ti['point_id']} (chain + geoid + lateral)",
                         f"{ti['tie_mm']:+.1f} +/- {ti['sigma_mm']:.1f}",
                         "the constant z_before_absolute.npy is built on"])
    rec_rows.append(["2021 anchor, combined (z_before_absolute.json)",
                     f"{side['datum_constant_mm']:+.1f} +/- {side['sigma_total_mm']:.1f}",
                     f"{side['unmodelled_bound_mm']:.1f} mm unmodelled bound held beside it"])
    diff = side["datum_constant_mm"] - np.median(t)
    rec_rows.append(["anchor - (2008 control, median of the open marks)",
                     f"{diff:+.1f}",
                     f"{diff / side['sigma_total_mm']:.2f} x the anchor's own sigma_total"])
    R.table(["quantity", "mm", "note"], rec_rows)

    if A.json:
        json.dump([dict(point=r["cp"].point_id, cover=r["cp"].point_type,
                        km=r["dist"], tile=r["tile"], tie_mm=r["est"].tie_mm,
                        tie_nolat_mm=r["raw"].tie_mm, sigma_mm=r["est"].sigma_mm,
                        spread_m=r["spread"], ctrl_pct=r["pct"], n_lines=r["n_lines"],
                        dnr_mm=1000 * cp_err(r["cp"], cps),
                        ladder=[dict(R=e.radius_m, n=e.n, z=e.z_lidar_m,
                                     fit_rms_mm=e.fit_rms_mm) for e in r["est"].curve])
                   for r in recs], open(A.json, "w"), indent=1)
        print(f"\nwrote {A.json}")

    R.done(headline="gen1's own 2008 control, measured on the point cloud with the "
                    "anchor's estimator")


_ERR = {}
_CTY = {}


def _county_of(point_id, cps):
    if not _CTY:
        with open(C._DATA / f"{CONTROL}.csv", newline="") as fh:
            for row in csv.DictReader(fh):
                _CTY.setdefault(row["point_id"], row["county"])
    return _CTY[point_id]


def cp_err(cp, cps):
    """The MnDNR table's own Error for this mark, m (a column of the bundled CSV)."""
    if not _ERR:
        with open(C._DATA / f"{CONTROL}.csv", newline="") as fh:
            for row in csv.DictReader(fh):
                _ERR[(row["point_id"], row["county"])] = float(row["dnr_error_m"])
        for (pid, _c), v in list(_ERR.items()):
            _ERR.setdefault(pid, v)
    return _ERR[cp.point_id]


if __name__ == "__main__":
    main()
