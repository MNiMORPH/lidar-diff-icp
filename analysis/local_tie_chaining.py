"""Is a gen1 flight-line tie a property of the LINE, or of the place it was fitted?

Elba's per-swath constants (``data/derived/elbaext/corrections_geoid.json``) were fitted
on a 4.5 x 4 km footprint and have been applied to surveyed marks up to 62.9 km along the
same lines. ``analysis/CROSS_LINE_FIT.md`` section 6 already showed the ACROSS-TRACK
coefficient is local to the overlap it was fitted on. This script asks the same question
of the CONSTANT, using ``lidar_diff_icp.localtie``:

  1. re-derive the 8.4 mm-per-link chain figure from the two saved Elba products;
  2. how much does the tie move with WINDOW SIZE at one place (both tie modes)?
  3. how much does it move ALONG TRACK, over the ~97 km of one flight-line pair that is
     on disk (both tie modes)?
  4. what does a locally-chained constant say at a real control mark, against the
     constant Elba's fit would have imported there?

Nothing here filters, drops or thresholds anything. Every window size, resolution, tie
mode and class filter is passed in explicitly; the ones I chose are declared ``MINE``.

Run:
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/local_tie_chaining.py \
        --section all
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from lidar_diff_icp import localtie as LT          # noqa: E402
from trust.provenance import Run                    # noqa: E402

warnings.filterwarnings("ignore", message=".*bottleneck.*")

ELBA_E, ELBA_N = 577825.0, 4884225.0               # corrections_geoid.json centroid
CACHE = "data/derived/localtie_cache"
CONTROL = "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv"
ELBAEXT = "data/derived/elbaext/corrections_geoid.json"
ELBA_OLD = "data/derived/elba/corrections.json"
ELBAEXT_OLD = "data/derived/elbaext/corrections.json"


def _tiles(pattern):
    return sorted(p for p in glob.glob(pattern) if "merged" not in p)


def _bbox(cache, path):
    return cache.header_bbox(path)


# ------------------------------------------------------------------ 1. the 8.4 mm/link

def section_perlink(R):
    """Re-derive the per-link chain error from the two Elba products, not from a memo."""
    a = json.load(open(ELBA_OLD))
    b = json.load(open(ELBAEXT_OLD))
    ca = {int(k): v[2] for k, v in a["per_swath_internal_alignment_dxdydz_m"].items()}
    cb = {int(k): v[2] for k, v in b["per_swath_internal_alignment_dxdydz_m"].items()}
    shared = sorted(set(ca) & set(cb))
    gauge = min(shared)                       # regauge both onto the lowest shared line
    rows = []
    per_link = []
    for s in shared:
        if s == gauge:
            continue
        d = ((ca[s] - ca[gauge]) - (cb[s] - cb[gauge])) * 1000.0
        nl = s - gauge                        # adjacent parallel lines: one link per step
        rows.append([str(s), str(nl), f"{d:+.1f}", f"{d / np.sqrt(nl):+.1f}"])
        per_link.append(d / np.sqrt(nl))
    R.column("line", "flight line, regauged onto line %d in BOTH products" % gauge)
    R.column("links", "steps from the gauge line along the adjacent-overlap chain")
    R.column("disagree_mm", "elba minus elbaext for that line's constant, mm")
    R.column("per_link_mm", "that disagreement divided by sqrt(links), mm")
    R.table(["line", "links", "disagree_mm", "per_link_mm"], rows)
    tot = np.array([float(r[2]) for r in rows])
    print(f"  RMS disagreement           = {np.sqrt(np.mean(tot ** 2)):.2f} mm")
    print(f"  RMS per link               = {np.sqrt(np.mean(np.square(per_link))):.2f} mm")
    print("  (the two products differ ONLY in the extent their constants were fitted on)")
    return float(np.sqrt(np.mean(np.square(per_link))))


# ------------------------------------------------------------------- 2. window ladder

LADDER_TILES = ["data/before/4342-29-64.laz", "data/before/4342-28-64.laz"]


def section_ladder(R, half_widths, res_m, exclude, shape):
    tiles = list(LADDER_TILES)
    cache = LT.TileCache(cache_dir=CACHE)
    R.column("pair", LT.LocalTie.columns()["pair"])
    R.column("tie_mode", "which coreg tie estimator: overlap_median or intercept")
    for k, v in LT.LocalTie.columns().items():
        if k != "pair":
            R.column(k, v)
    out = {}
    for a, b in [(135, 136), (136, 137), (137, 138)]:
        op = LT.nearest_overlap_point(tiles, a, b, easting=ELBA_E, northing=ELBA_N,
                                      res_m=res_m, exclude=exclude, cache=cache)
        print(f"\n  pair {a}-{b}: nearest overlap cell to the Elba centroid is "
              f"({op.easting:.0f}, {op.northing:.0f}), {op.distance_m:.0f} m away, "
              f"{op.n_overlap_cells} overlap cells in these two tiles")
        for tie in ("overlap_median", "intercept"):
            lad = LT.window_ladder(tiles, a, b, easting=op.easting, northing=op.northing,
                                   half_widths_m=half_widths, shape=shape, res_m=res_m,
                                   tie=tie, exclude=exclude, cache=cache)
            R.table(["pair", "tie_mode"] + [k for k in LT.LocalTie.columns() if k != "pair"],
                    [[r[0], tie] + r[1:] for r in lad.rows()])
            print(f"    {tie:15s} spread over the ladder = {lad.spread_mm:6.1f} mm, "
                  f"sd = {lad.sd_mm:5.1f} mm")
            out[(a, b, tie)] = lad
    cache.release()
    return out


# --------------------------------------------------------------- 3. along-track sweep

def section_alongtrack(R, pair, tile_glob, half_width, res_m, exclude, shape, centre):
    a, b = pair
    cache = LT.TileCache(cache_dir=CACHE)
    tiles = []
    for p in _tiles(tile_glob):
        try:
            cache.header_bbox(p)
        except Exception as e:                      # a tile being written by another job
            print(f"  SKIP {p}: {e}")
            continue
        tiles.append(p)
    rows = []
    ties = []
    for p in tiles:
        x0, y0, x1, y1 = cache.header_bbox(p)
        try:
            # two passes: the first finds the pair's overlap strip in this tile, the
            # second lands on the strip's own CENTRE at the tile's mid-northing. Centring
            # on the tile centre instead puts every window on the strip's western edge,
            # where the across-track term is sampled one-sidedly and the intercept tie
            # becomes an extrapolation -- measured, and reported in the write-up.
            first = LT.nearest_overlap_point([p], a, b, easting=(x0 + x1) / 2,
                                             northing=(y0 + y1) / 2, res_m=res_m,
                                             exclude=exclude, cache=cache)
            op = first if centre == "tile" else LT.nearest_overlap_point(
                [p], a, b, easting=first.median_easting, northing=(y0 + y1) / 2,
                res_m=res_m, exclude=exclude, cache=cache)
        except ValueError as e:
            print(f"  {os.path.basename(p)}: {e}")
            cache.release(p)
            continue
        row = [os.path.basename(p), f"{op.northing:.0f}",
               f"{(op.northing - ELBA_N) / 1000.0:+.1f}", f"{op.easting:.0f}"]
        got = {}
        for tie in ("overlap_median", "intercept"):
            t = LT.pair_tie_at([p], a, b, easting=op.easting, northing=op.northing,
                               half_width_m=half_width, shape=shape, res_m=res_m,
                               tie=tie, exclude=exclude, cache=cache)
            got[tie] = t
            row += [f"{t.dz_m * 1000:+.1f}", f"{t.dz_sigma_m * 1000:.1f}",
                    ("D" if t.degenerate else "") + ("X" if t.extrapolated else "")]
            ties.append((tie, op.northing, t))
        t_med = got["overlap_median"]
        row += [f"{t_med.dz_overlap_median_m * 1000:+.1f}",
                f"{got['intercept'].c_mm_per_tan:+.0f}", str(t_med.n_overlap_cells)]
        rows.append(row)
        print("  " + "  ".join(row), flush=True)
        cache.release(p)                            # shared laptop: one tile at a time
    hdr = ["tile", "northing", "d_elba_km", "easting", "dz_med_mm", "sig_med_mm",
           "flag_med", "dz_int_mm", "sig_int_mm", "flag_int", "dz_plain_mm", "c_mm_tan",
           "ovl_cells"]
    defs = {
        "tile": "the gen1 tile the window was cut from",
        "northing": "northing of the window centre (nearest overlap cell to tile centre), m",
        "d_elba_km": "signed along-track distance from the Elba centroid, km",
        "easting": "easting of the window centre, m",
        "dz_med_mm": f"tie for the pair (ref={a}, src={b}) in that window with "
                     f"tie='overlap_median' -- what to ADD to line {b} to bring it into "
                     f"line {a}'s frame there, mm",
        "sig_med_mm": "coreg's formal 1-sigma for it, mm",
        "flag_med": "D if the Nuth & Kaeaeb fit used zero cells there, X if the window "
                    "does not sample dtan = 0",
        "dz_int_mm": f"the same tie with tie='intercept' (LAD intercept at across-track "
                     f"position zero), mm",
        "sig_int_mm": "coreg's formal 1-sigma for it, mm",
        "flag_int": "same two flags for the intercept run",
        "dz_plain_mm": f"median of the UNSHIFTED overlap difference z({a}) - z({b}) in "
                       f"that window -- no horizontal solution applied, mm",
        "c_mm_tan": "across-track slope fitted with the intercept there, mm per unit tangent",
        "ovl_cells": "cells where both lines have terrain returns in that window",
    }
    for k, v in defs.items():
        R.column(k, v)
    R.table(hdr, rows)
    for tie, col in (("overlap_median", 4), ("intercept", 7)):
        v = np.array([float(r[col]) for r in rows])
        s = np.array([float(r[col + 1]) for r in rows])
        print(f"  {tie:15s} n={v.size}  mean {v.mean():+7.1f}  sd {v.std(ddof=1):6.1f}  "
              f"range {v.min():+7.1f} to {v.max():+7.1f}  "
              f"(mean formal sigma {s.mean():.1f} mm)")
    return rows, hdr


# ------------------------------------------------- 3c. is the scatter short or long range?

def section_repeat(R, pair, tile, half_widths, step_m, res_m, exclude, shape):
    """Tie the SAME pair at a row of windows a few hundred metres apart in ONE tile.

    Section 3 measures ~30 mm of scatter over 94 km, against a formal sigma near 1 mm.
    That scatter is only evidence of an along-track DRIFT if it is bigger than what the
    same estimator produces between windows a few hundred metres apart. This measures
    that directly, at two window sizes, so the two can be compared on the same footing.
    """
    a, b = pair
    cache = LT.TileCache(cache_dir=CACHE)
    x0, y0, x1, y1 = cache.header_bbox(tile)
    first = LT.nearest_overlap_point([tile], a, b, easting=(x0 + x1) / 2,
                                     northing=(y0 + y1) / 2, res_m=res_m,
                                     exclude=exclude, cache=cache)
    e = first.median_easting
    rows = []
    for hw in half_widths:
        ns = np.arange(y0 + hw + 1.0, y1 - hw - 1.0, step_m)
        vals = {"overlap_median": [], "intercept": []}
        for n in ns:
            for tie in vals:
                t = LT.pair_tie_at([tile], a, b, easting=e, northing=float(n),
                                   half_width_m=hw, shape=shape, res_m=res_m, tie=tie,
                                   exclude=exclude, cache=cache)
                vals[tie].append(t.dz_m * 1000.0)
        for tie, v in vals.items():
            v = np.array(v, float)
            rows.append([f"{hw:.0f}", tie, str(v.size), f"{v.mean():+.1f}",
                         f"{v.std(ddof=1):.1f}", f"{v.min():+.1f}", f"{v.max():+.1f}",
                         f"{(ns.max() - ns.min()) / 1000:.1f}"])
            print("  " + "  ".join(rows[-1]), flush=True)
    cache.release()
    defs = {
        "half_width_m": "window half-width about each sample location, m",
        "tie_mode": "which coreg tie estimator produced the number, overlap_median or intercept",
        "n_windows": f"windows tied along the (ref={a}, src={b}) overlap inside "
                     f"{os.path.basename(tile)}, count",
        "mean_mm": "mean of the tie over those windows, mm",
        "sd_mm": "standard deviation of the tie between those windows, mm",
        "min_mm": "smallest tie among those windows, mm",
        "max_mm": "largest tie among those windows, mm",
        "span_km": "along-track distance from the first window centre to the last, km",
    }
    for k, v in defs.items():
        R.column(k, v)
    R.table(["half_width_m", "tie_mode", "n_windows", "mean_mm", "sd_mm", "min_mm",
             "max_mm", "span_km"], rows)
    return rows


# ----------------------------------------------------------------- 4. marks and chains

def _marks_in(cache, tiles, emin, emax):
    rows = list(csv.DictReader(open(CONTROL)))
    boxes = [cache.header_bbox(p) for p in tiles]
    out = []
    for r in rows:
        e = float(r["easting"]); n = float(r["northing"])
        if not (emin <= e <= emax):
            continue
        for p, (x0, y0, x1, y1) in zip(tiles, boxes):
            if x0 <= e <= x1 and y0 <= n <= y1:
                out.append((r["point_id"], e, n, p))
                break
    return sorted(out, key=lambda t: abs(t[2] - ELBA_N))


def section_marks(R, half_width, res_m, exclude, shape, radius_m, target_line,
                  ladder_half_widths):
    cache = LT.TileCache(cache_dir=CACHE)
    tiles = [p for p in _tiles("data/before/4342-*.laz") + _tiles("data/before/5142-*.laz")]
    tiles = [p for p in tiles if 574000 <= (cache.header_bbox(p)[0]) <= 580500]
    marks = _marks_in(cache, tiles, 574000.0, 580500.0)
    imported = {int(k): v[2] for k, v in json.load(open(ELBAEXT))
                ["per_swath_internal_alignment_dxdydz_m"].items()}
    from lidar_diff_icp.groundtruth import chain as gt_chain
    rows = []
    boxes = {p: cache.header_bbox(p) for p in tiles}
    for pid, e, n, path in marks:
        # every tile of the same TILE ROW -- those whose northing range contains the
        # mark -- so a chain can walk east or west out of the mark's own tile. This is a
        # geometric selection (which tiles cover this latitude), not a tuned radius.
        near = [p for p in tiles if boxes[p][1] <= n <= boxes[p][3]]
        inv = LT.inventory_from_cache(near, exclude=exclude, cache=cache)
        cov = gt_chain.covering_lines(inv, e, n, radius_m, res=res_m)
        cov = {k: v for k, v in cov.items() if k in imported}
        if not cov:
            print(f"  {pid}: no line of the imported set within {radius_m} m")
            for p in near:
                cache.release(p)
            continue
        src = max(cov, key=cov.get)
        if target_line not in inv.lines:
            note = f"line {target_line} absent from this tile"
            rows.append([pid, f"{(n - ELBA_N) / 1000:+.1f}", str(src), str(cov[src]),
                         "-", "-", "-", "-", "-", note])
            print("  " + "  ".join(rows[-1]), flush=True)
            for p in near:
                cache.release(p)
            continue
        try:
            ch = LT.chain_local(near, easting=e, northing=n, source_line=src,
                                target_line=target_line, half_width_m=half_width,
                                shape=shape, res_m=res_m, tie="intercept",
                                exclude=exclude, cache=cache,
                                ladder_half_widths_m=ladder_half_widths)
        except ValueError as exc:
            rows.append([pid, f"{(n - ELBA_N) / 1000:+.1f}", str(src), str(cov[src]),
                         "-", "-", "-", "-", "-", str(exc)[:60]])
            print("  " + "  ".join(rows[-1]), flush=True)
            for p in near:
                cache.release(p)
            continue
        cmp = LT.compare_to_constants(ch, imported)
        rows.append([pid, f"{(n - ELBA_N) / 1000:+.1f}", str(src), str(cov[src]),
                     "-".join(str(x) for x in ch.nodes), f"{cmp.local_mm:+.1f}",
                     f"{ch.dz_sigma_formal_m * 1000:.1f}",
                     f"{ch.dz_sigma_window_m * 1000:.1f}",
                     f"{cmp.imported_mm:+.1f}", f"{cmp.difference_mm:+.1f}"])
        print("  " + "  ".join(rows[-1]), flush=True)
        for p in near:
            cache.release(p)
    hdr = ["mark", "d_elba_km", "line", "n_near", "path", "local_mm", "sig_formal_mm",
           "sig_window_mm", "imported_mm", "diff_mm"]
    defs = {
        "mark": "MnDNR 2008 control point id",
        "d_elba_km": "signed along-track distance from the Elba centroid, km",
        "line": f"flight line with the most terrain returns within {radius_m} m of the mark",
        "n_near": f"how many of its returns are within {radius_m} m",
        "path": f"lines walked from that line to line {target_line}; a single entry means "
                f"the mark is ALREADY on the target line and no transfer is needed",
        "local_mm": f"constant to ADD to the mark's line to put it into line "
                    f"{target_line}'s frame, measured AT THE MARK, mm",
        "sig_formal_mm": "coreg's per-link sigmas in quadrature, mm -- known optimistic",
        "sig_window_mm": "per-link window-ladder spreads in quadrature, mm",
        "imported_mm": f"the same constant implied by elbaext's fitted set "
                       f"(dz[line] - dz[{target_line}]), mm",
        "diff_mm": "local minus imported, mm -- what the extrapolation costs here",
    }
    for k, v in defs.items():
        R.column(k, v)
    R.table(hdr, rows)
    got = [(r[4], float(r[-1])) for r in rows
           if r[-1] not in ("-",) and not r[-1][0].isalpha()]
    for label, sel in (("all marks with a route", got),
                       ("of those, ALREADY on line %d (zero links)" % target_line,
                        [g for g in got if "-" not in g[0]]),
                       ("of those, needing at least one link", 
                        [g for g in got if "-" in g[0]])):
        d = [v for _, v in sel]
        if d:
            print(f"  {label:52s} n={len(d):2d}  mean {np.mean(d):+7.1f}  "
                  f"sd {np.std(d, ddof=1) if len(d) > 1 else float('nan'):6.1f}  "
                  f"RMS {np.sqrt(np.mean(np.square(d))):6.1f} mm")
    return rows


# ------------------------------------------------------------------- input declaration

def declare_inputs(R, A):
    """Every file this run will read, declared before the banner prints."""
    if A.section in ("all", "perlink"):
        R.input(ELBA_OLD, role="elba product: per-swath constants fitted on the "
                               "2.54 x 3.50 km elba footprint")
        R.input(ELBAEXT_OLD, role="elbaext product: the SAME swaths fitted on the "
                                  "4.45 x 4.05 km elbaext footprint")
    if A.section in ("all", "marks"):
        R.input(ELBAEXT, role="the per-swath constants actually applied to the control "
                              "marks today, and the swath_tie mode they were fitted with")
        R.input(CONTROL, role="the 1004 MnDNR 2008 ground control points -- gen1's own "
                              "control; used here only for their POSITIONS")
    seen = set()
    for p in (LADDER_TILES if A.section in ("all", "ladder", "repeat", "synthesis") else []):
        seen.add(p)
        R.input(p, role="gen1 tile over the Elba study area, for the window ladder and "
                        "the within-tile repeat")
    if A.section in ("all", "alongtrack", "marks", "synthesis"):
        for pat in ("data/before/*-*-63.laz", "data/before/*-*-64.laz"):
            for p in _tiles(pat):
                if p in seen:
                    continue
                seen.add(p)
                R.input(p, role="gen1 tile on the lines 133-138 corridor, an along-track "
                                "sample of their overlaps and a holder of control marks")


# ------------------------------------------------------------------------ 5. synthesis

def section_synthesis(R, along, repeat, half_width):
    """Is the along-track scatter bigger than the same estimator's own short-range scatter?

    Section 3 measures the tie at one window per tile over ~90 km; section 3c measures it
    at a row of windows a few hundred metres apart inside ONE tile, same pair, same
    estimator, same window size. If the long-range scatter is no bigger than the
    short-range scatter, the along-track variation is the estimator talking to itself and
    there is nothing to see. The variance ratio settles it.
    """
    from scipy import stats
    rows = []
    for tie, col in (("overlap_median", 4), ("intercept", 7)):
        long = np.array([float(r[col]) for r in along], float)
        # `repeat` rows carry summary statistics, not the individual windows: use the
        # sd and n it printed.
        rr = [r for r in repeat if r[0] == f"{half_width:.0f}" and r[1] == tie][0]
        s_short = float(rr[4]); n_short = int(rr[2])
        s_long = float(np.std(long, ddof=1)); n_long = long.size
        F = (s_long / s_short) ** 2
        p = float(stats.f.sf(F, n_long - 1, n_short - 1))
        # separation vs |difference|, over every pair of along-track samples
        y = np.array([float(r[1]) for r in along], float)
        iu = np.triu_indices(y.size, 1)
        sep = np.abs(y[iu[0]] - y[iu[1]]) / 1000.0
        dif = np.abs(long[iu[0]] - long[iu[1]])
        r_pearson = float(np.corrcoef(sep, dif)[0, 1])
        rows.append([tie, str(n_long), f"{s_long:.1f}", str(n_short), f"{s_short:.1f}",
                     f"{F:.2f}", f"{p:.4f}", str(sep.size), f"{r_pearson:+.3f}"])
        print("  " + "  ".join(rows[-1]))
    defs = {
        "tie_mode": "which coreg tie estimator both scatters were measured with",
        "n_long": "along-track samples, one window per tile over the corridor, count",
        "sd_long_mm": "standard deviation of the tie between them, mm",
        "n_short": "windows inside ONE tile, a few hundred metres apart, count",
        "sd_short_mm": "standard deviation of the tie between THOSE, mm",
        "F": "variance ratio, long-range over short-range, dimensionless",
        "p": "one-sided p for that ratio, F(n_long-1, n_short-1)",
        "n_pairs": "pairs of along-track samples used for the correlation, count",
        "r_sep_vs_diff": "Pearson r between the separation of two samples and the "
                         "absolute difference of their ties, dimensionless",
    }
    for k, v in defs.items():
        R.column(k, v)
    R.table(["tie_mode", "n_long", "sd_long_mm", "n_short", "sd_short_mm", "F", "p",
             "n_pairs", "r_sep_vs_diff"], rows)
    return rows


# ---------------------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all",
                    choices=["all", "perlink", "ladder", "alongtrack", "repeat",
                             "synthesis", "marks"])
    ap.add_argument("--res-m", type=float, default=2.0)
    ap.add_argument("--half-width-m", type=float, default=400.0)
    ap.add_argument("--ladder-m", type=float, nargs="+",
                    default=[100.0, 200.0, 400.0, 800.0, 1200.0])
    ap.add_argument("--shape", default="square", choices=["square", "disk"])
    ap.add_argument("--exclude", type=int, nargs="+", default=[5, 6, 9])
    ap.add_argument("--repeat-m", type=float, nargs="+", default=[400.0, 800.0])
    ap.add_argument("--repeat-step-m", type=float, default=300.0)
    ap.add_argument("--mark-radius-m", type=float, default=25.0)
    ap.add_argument("--target-line", type=int, default=137)
    ap.add_argument("--centre", nargs="+", default=["strip", "tile"],
                    choices=["strip", "tile"],
                    help="where in the pair's overlap the along-track windows sit: on the "
                         "STRIP's own centre, or at the overlap cell nearest the TILE centre")
    A = ap.parse_args()

    R = Run("is a gen1 flight-line vertical tie a property of the LINE, or of the place "
            "it was fitted?")
    R.param("res_m", A.res_m, src="repo",
            why="coreg.coregister_swaths / align_swaths default grid resolution (2.0 m)")
    R.param("exclude", tuple(A.exclude), src="repo",
            why="coreg.coregister_swaths terrain proxy ~isin(classification,(5,6,9)) -- "
                "the VENDOR classes, so no CSF run is needed")
    R.param("tie_modes", ("overlap_median", "intercept"), src="repo",
            why="the two estimators coreg.coregister_swaths supports; "
                "corrections_geoid.json records swath_tie='intercept' as the pipeline's")
    R.param("half_width_m", A.half_width_m, src="MINE",
            why="window half-width for the along-track sweep and the chains. There is NO "
                "library default and none is invented here; 400 m is one rung of the "
                "ladder in section 2, chosen so the window fits inside a single 2.5 km "
                "tile at every along-track sample. Section 2 measures what the choice "
                "costs; section 3 is re-runnable at any other rung with --half-width-m")
    R.param("ladder_half_widths_m", tuple(A.ladder_m), src="MINE",
            why="the window sizes reported in section 2. Chosen to span a factor of 12, "
                "from below one flight-line sidelap width to a whole tile; nothing is "
                "excluded on their basis, the whole ladder is printed")
    R.param("shape", A.shape, src="MINE",
            why="square window. A disk would give a ragged grid edge and fewer cells for "
                "the same half-width; localtie supports both and --shape switches it")
    R.param("mark_radius_m", A.mark_radius_m, src="MINE",
            why="radius within which a mark's flight line is identified by return count "
                "(groundtruth.chain.covering_lines). It selects WHICH LINE the mark is "
                "on, not which marks are used; every mark is reported")
    R.param("repeat_step_m", A.repeat_step_m, src="MINE",
            why="along-track spacing of the within-tile repeat windows in section 3c. "
                "300 m so that neighbouring 400 m-half-width windows overlap and "
                "neighbouring 800 m ones overlap heavily; it sets the SAMPLING of the "
                "short-range scatter, and excludes nothing")
    R.param("repeat_half_widths_m", tuple(A.repeat_m), src="MINE",
            why="the two window sizes section 3c repeats at, so that scatter which is "
                "estimator noise (shrinks with window) can be told from scatter which is "
                "local structure (does not)")
    R.param("centre", tuple(A.centre), src="MINE",
            why="where each along-track window sits across the overlap strip. 'strip' is "
                "the strip's own median easting -- the sidelap middle, where the "
                "intercept tie is an interpolation; 'tile' is the overlap cell nearest "
                "the tile centre, which on a N-S strip is its western EDGE. Both are run "
                "and both are reported; neither is used to exclude anything")
    R.param("target_line", A.target_line, src="MINE",
            why="the frame the marks are carried into. 137 is the middle Elba line and "
                "the one the brief names; --target-line changes it")
    declare_inputs(R, A)
    R.banner()

    if A.section in ("all", "perlink"):
        print("\n=== 1. per-link chain error, re-derived from the two Elba products ===")
        section_perlink(R)
    if A.section in ("all", "ladder"):
        print("\n=== 2. how the tie moves with WINDOW SIZE, at the Elba centroid ===")
        section_ladder(R, A.ladder_m, A.res_m, tuple(A.exclude), A.shape)
    along = {}
    repeat = None
    if A.section in ("all", "alongtrack", "synthesis"):
        for centre in A.centre:
            print(f"\n=== 3a. how the tie moves ALONG TRACK: pair 135-136, column-63 "
                  f"tiles, window centred on the overlap {centre} ===")
            along[(135, 136, centre)] = section_alongtrack(
                R, (135, 136), "data/before/*-*-63.laz", A.half_width_m,
                A.res_m, tuple(A.exclude), A.shape, centre)[0]
            print(f"\n=== 3b. the same for pair 136-137, column-64 tiles, {centre} ===")
            along[(136, 137, centre)] = section_alongtrack(
                R, (136, 137), "data/before/*-*-64.laz", A.half_width_m,
                A.res_m, tuple(A.exclude), A.shape, centre)[0]
    if A.section in ("all", "repeat", "synthesis"):
        print("\n=== 3c. the same pair tied at a row of windows inside ONE tile ===")
        repeat = section_repeat(R, (136, 137), "data/before/4342-29-64.laz", A.repeat_m,
                                A.repeat_step_m, A.res_m, tuple(A.exclude), A.shape)
    if A.section in ("all", "synthesis") and repeat is not None and along:
        print("\n=== 3d. long-range against short-range scatter, pair 136-137, "
              f"{A.half_width_m:.0f} m windows ===")
        section_synthesis(R, along[(136, 137, "strip")], repeat, A.half_width_m)
    if A.section in ("all", "marks"):
        print(f"\n=== 4. a local chain at each control mark, against elbaext's constants "
              f"(target line {A.target_line}) ===")
        section_marks(R, A.half_width_m, A.res_m, tuple(A.exclude), A.shape,
                      A.mark_radius_m, A.target_line, A.ladder_m)

    R.done()


if __name__ == "__main__":
    main()
