#!/usr/bin/env python3
"""Make the swath-to-swath tie extent-invariant, and measure what it buys.

``coreg.align_swaths`` ties each gen1 flight line to its neighbour by the MEDIAN height
difference over their overlap. ``analysis/SWATH_ACROSS_TRACK_TEST.md`` measured that this
difference is not flat across the sidelap: pooled ``c = +80.0 +- 5.7`` mm per unit tangent
(t = 14.1) over 87,737 cell-pairs, +108.0 +- 3.8 over all 506,705 in-grid overlap cells,
per flight-line pair rather than sensor-wide (Wald W = 314.2, df = 7, p = 6e-64), not
spatial and not canopy. The tie is therefore ``k + c * mean(dtan)`` -- an average over
whatever part of the sidelap the tile covers -- so **the tie depends on the tile extent**,
and two tiles get different constants for the same pair of flight lines.

This run implements the two changes that report proposes, WITHOUT changing any default:

1. **Extent-invariant tie.** ``coreg.coregister_swaths(..., tie="intercept")`` estimates
   the tie as the LAD (median-regression) intercept at ``tan(scan_ref) = tan(scan_src)``
   -- the middle of the sidelap, a fixed geometric position -- instead of the unweighted
   overlap median. LAD is used so the estimator reduces exactly to the shipped median when
   the across-track slope is zero, which is what makes the two numbers comparable.
2. **A zero line that is not an edge-cut swath.** ``pipeline.py:670`` pins the network to
   ``int(ps8.min())``, the lowest-numbered line, which at both Elba tiles is the one whose
   nadir track falls outside the tile (elba/135 at mean scan -11.79 deg, elbaext/133 at
   -12.64 deg, against interior lines near 0). Its constant and its across-track slope
   correlate at 0.989/0.991 there.

A shared roll term is NOT added: the across-track report rejects one at p = 6e-64 and the
repo's own ``boresight.estimate_boresight`` rejects it in its own numbers.

Everything is measured against the shipped pipeline on the same inputs:

* the ``overlap_median`` ties are re-derived from the cached CSF clouds and required to
  reproduce ``corrections.json`` to its own written precision before anything is compared;
* the DoD consequence is computed by re-reducing the per-cell gen1 ground from
  ``beam_offset_table.parquet`` with the new per-swath constants -- the pipeline's own
  estimator, ``ground_q = 0.50`` of the slope-normal residual per cell -- and that path is
  itself validated by reproducing an independently computed number
  (``analysis/STABLE_POINT_TILT_AUDIT.md`` section 5c, "per-swath constants" removed);
* the tilt and the stable-ground scatter are then read by running
  ``analysis/stable_point_tilt_audit.py --dod`` on the old and the new rasters, so the two
  are produced by identical code.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/swath_tie_intercept.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import laspy
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lidar_diff_icp import coreg, io
from lidar_diff_icp.swathdiff import _median_grid
from trust.provenance import Run

TILES = {"elba_fulldensity": "data/csf_cache/elba.las",
         "elbaext": "data/csf_cache/elbaext.las"}
JSON_PRECISION_M = 5e-5      # corrections.json writes the shifts rounded to 4 decimals (m)
REPRO_TOL_M = 1e-3           # see R.param("reproduction tolerance") -- declared, not silent


def load_pc(las_path, R):
    """The gen1 cloud exactly as ``pipeline.difference_dem`` builds it before aligning.

    ``pipeline.py`` reads the archived CSF cloud, converts the LAS 1.4 scan angle from its
    0.006 deg storage unit (``pipeline.py:660``), and hands ``coreg.align_swaths`` a
    ``PointCloud`` of ALL points -- ``align_swaths`` does its own terrain-class selection.
    Reproduced here line for line so the ties compared below are the pipeline's own.
    """
    f = laspy.read(R.input(las_path, role="gen1 2008 cloud with PDAL CSF ground "
                                          "classification, cached by pipeline.difference_dem "
                                          "(ground_source='csf'), raw coordinates"))
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    ps = np.asarray(f.point_source_id); cl = np.asarray(f.classification)
    sa = np.asarray(f.scan_angle).astype(float) * 0.006
    del f
    return io.PointCloud(x, y, z, ps, cl, np.zeros_like(z), sa, "EPSG:26915"), ps


def pair_mechanism(pc, a, b, res=2.0, exclude=(5, 6, 9)):
    """``(mean_dtan, c, k, med, cells)`` for one pair, on the co-registration's own grids.

    Rebuilds exactly what ``coreg.coregister_swaths`` builds -- same cells, same shift --
    and reports the across-track slope and the position the shipped median averages at, so
    the change in the tie can be read against the mechanism instead of only asserted.
    """
    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == a)
    mb = terr & (pc.point_source_id == b)
    x, y, z = pc.x, pc.y, pc.z
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    nx = int(np.ceil((x1 - x0) / res)); ny = int(np.ceil((y1 - y0) / res))
    za = _median_grid(x[ma], y[ma], z[ma], res, x0, y0, nx, ny)
    zb = _median_grid(x[mb], y[mb], z[mb], res, x0, y0, nx, ny)
    c0 = coreg.nuth_kaab(za, zb, res)
    sa = np.asarray(pc.scan_angle, float)
    ta = _median_grid(x[ma], y[ma], np.tan(np.radians(sa[ma])), res, x0, y0, nx, ny)
    tb = _median_grid(x[mb], y[mb], np.tan(np.radians(sa[mb])), res, x0, y0, nx, ny)
    dh = za - coreg._shift_grid(zb, c0.dx, c0.dy, res)
    dt = ta - coreg._shift_grid(tb, c0.dx, c0.dy, res)
    k, cc, n = coreg.across_track_tie(dh, dt)
    m = np.isfinite(dh) & np.isfinite(dt)
    return float(np.mean(dt[m])), cc, k, float(np.nanmedian(dh)), n


def edge_map(edges):
    return {(int(a), int(b)): (dx, dy, dz, n) for a, b, dx, dy, dz, n in edges}


def reref(d, s):
    return {k: v - d[s] for k, v in d.items()}


def f2(v):
    return f"{v:+.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", nargs="+", default=list(TILES))
    ap.add_argument("--dod-tile", default="data/derived/elba_fulldensity",
                    help="the tile whose DoD raster is rebuilt with the new constants")
    ap.add_argument("--dod", default="dod_cover_q2.npy")
    A = ap.parse_args()

    R = Run("does taking each swath-to-swath tie as the fitted intercept at across-track "
            "position zero, instead of the unweighted overlap median, remove the tile-extent "
            "dependence of the gen1 per-swath constants -- and what does it do to the DoD?")
    R.param("tie estimators compared", "overlap_median (shipped) vs intercept (new)", src="repo")
    R.param("res (pairwise co-registration grid)", 2.0, src="repo")
    R.param("exclude (non-terrain classes)", (5, 6, 9), src="repo")
    R.param("ground_q", 0.50, src="repo")
    R.param("minimum counts imposed", "none", src="repo")
    R.param("reproduction tolerance (vertical)", REPRO_TOL_M, src="MINE",
            why="the gate on whether the re-derived shipped ties ARE the shipped ties. 1 mm "
                "is the coarsest level at which a reproduction failure could still matter "
                "here: it is a sixth of the smallest tie change measured below and an "
                "eighth of the smallest tile-to-tile disagreement being explained. The "
                "achieved values are printed per tile and per axis; it excludes no data.")

    for c, d in [
        ("tile", "derived-product directory / cached cloud the numbers come from"),
        ("pair", "adjacent flight-line pair A-B (point_source_id), A < B"),
        ("cells", "overlap grid cells the pairwise tie used"),
        ("dz_med", "pairwise vertical tie from the shipped overlap median, mm "
                   "(+ = line B must be raised onto line A)"),
        ("dz_int", "the same tie from the across-track intercept at dtan = 0, mm"),
        ("delta", "dz_int - dz_med, mm: what the extent-dependence was worth on that link"),
        ("c_mm", "across-track slope fitted with the intercept, mm per unit tangent"),
        ("mean_dtan", "mean of tan(scan_A) - tan(scan_B) over the pair's overlap cells in "
                      "this tile: the across-track position the shipped tie averages at"),
        ("pred", "-c_mm * mean_dtan, mm: the shift the ramp predicts for this link, against "
                 "which delta is the measured value"),
        ("swath", "gen1 flight line (point_source_id)"),
        ("returns", "in-grid gen1 CSF ground returns of that line"),
        ("scan_mean", "mean scan angle of that line's in-grid returns, deg (signed)"),
        ("scan_sd", "standard deviation of that line's scan angle, deg"),
        ("scan_p5", "5th percentile of that line's scan angle, deg"),
        ("scan_p95", "95th percentile of that line's scan angle, deg"),
        ("asym", "|mean scan| of the line, deg: 0 means the tile samples it symmetrically "
                 "about nadir, large means the tile only ever sees it off to one side"),
        ("dz_old", "per-swath constant from the shipped tie, re-referenced as stated, mm"),
        ("dz_new", "per-swath constant from the intercept tie, re-referenced as stated, mm"),
        ("elba_old", "elba's constant for that swath, shipped tie, re-referenced to 135, mm"),
        ("ext_old", "elbaext's constant for the same swath, shipped tie, ref 135, mm"),
        ("dis_old", "elba_old - ext_old: the tile-to-tile disagreement, shipped tie, mm"),
        ("elba_new", "elba's constant, intercept tie, re-referenced to 135, mm"),
        ("ext_new", "elbaext's constant, intercept tie, ref 135, mm"),
        ("dis_new", "elba_new - ext_new: the same disagreement, intercept tie, mm"),
        ("removed", "dis_old - dis_new, mm, and as a percentage of dis_old"),
        ("line", "which flight line the level is expressed against -- the tile's own "
                 "ZERO LINE, the proposed COMMON LINE, or zero-mean -- and with which tie"),
        ("level", "point-count-weighted mean vertical shift of the gen1 cloud against that "
                  "line, mm: the level the tile's whole DoD inherits"),
        ("raster", "DoD raster file this run wrote"),
        ("md_mm", "median change in the DoD against the shipped raster, mm"),
        ("lo_mm", "most negative change in the DoD, mm"),
        ("hi_mm", "most positive change in the DoD, mm"),
    ]:
        R.column(c, d)

    # ---------------------------------------------------------------- 0. the ties
    sol, pcs = {}, {}
    for name in A.tiles:
        tile = f"data/derived/{name}"
        meta = json.load(open(R.input(f"{tile}/corrections.json",
                                      role="pipeline corrections actually shipped for this "
                                           "tile: grid geometry and the per-swath constants "
                                           "this run must reproduce before changing them")))
        pc, ps = load_pc(TILES[name], R)
        zero_line = int(ps.min())                    # pipeline.py:670
        out, zeromean = {}, {}
        for tie in ("overlap_median", "intercept"):
            corr, edges, _ = coreg.align_swaths(pc, ref=zero_line, tie=tie)
            out[tie] = (corr, edge_map(edges))
            # align_swaths' other origin, already supported: zero-mean over the swaths
            # present, i.e. no line is privileged. Solved rather than derived so the
            # level below is the function's own answer.
            zm, _, _ = coreg.align_swaths(pc, ref=None, tie=tie)
            zeromean[tie] = zm
        swaths = sorted(out["overlap_median"][0])
        mech = {(a, b): pair_mechanism(pc, a, b)
                for a, b in zip(swaths[:-1], swaths[1:])}
        del pc
        shipped = {int(k): v for k, v in
                   meta["per_swath_internal_alignment_dxdydz_m"].items()}
        # Gate against the estimator the tile ACTUALLY shipped, read from its own
        # corrections.json, not a hardcoded one. The intercept tie became the pipeline
        # default on 2026-08-26 at 17:16, so a fixed "overlap_median" here compares two
        # different estimators and fails by the size of the tie change itself (5.9 mm at
        # elba) -- the gate asking the wrong question rather than catching a fault.
        # Tiles whose corrections.json predates the change record no swath_tie; those
        # shipped the overlap median.
        shipped_tie = meta.get("swath_tie") or "overlap_median"
        if shipped_tie not in out:
            raise ValueError(f"{name}: corrections.json records swath_tie={shipped_tie!r}, "
                             f"which this script does not solve for ({sorted(out)})")
        got = out[shipped_tie][0]
        worst = [max(abs(got[s][k] - shipped[s][k]) for s in shipped) for k in range(3)]
        # The assertion is on dz, which is the axis this run changes and the axis every
        # number below is read from. The horizontal components are reported, not asserted:
        # elbaext's corrections.json predates the current cached cloud, so its dx/dy are
        # not expected to land on the same value and the size of that gap is printed.
        assert worst[2] <= REPRO_TOL_M, (
            f"{name}: the shipped VERTICAL tie ({shipped_tie}) was NOT reproduced "
            f"(worst {worst[2]:.3g} m against a {REPRO_TOL_M:g} m tolerance) -- "
            f"nothing below is comparable")
        sol[name] = dict(meta=meta, zero_line=zero_line, out=out, worst=worst, mech=mech,
                         zeromean=zeromean, shipped_tie=shipped_tie)

    R.banner()
    print()
    print("## 0. The shipped ties, re-derived from the cached clouds\n")
    print("  align_swaths(pc, ref=int(ps.min())) on the same cached CSF cloud reproduces the\n"
          "  per-swath constants in each corrections.json to the precision that file is\n"
          "  written to. Only after that is the intercept tie a like-for-like substitution\n"
          "  rather than a different pipeline.\n")
    for name in A.tiles:
        w = sol[name]["worst"]
        st = os.stat(f"data/derived/{name}/corrections.json").st_mtime
        sc = os.stat(TILES[name]).st_mtime
        print(f"  {name}: worst |re-derived - shipped| over all swaths -- "
              f"dx {1000*w[0]:.2f} mm, dy {1000*w[1]:.2f} mm, dz {1000*w[2]:.2f} mm "
              f"(file precision {1000*JSON_PRECISION_M:g} mm)")
        print(f"      corrections.json written {pd.Timestamp(st, unit='s'):%Y-%m-%d %H:%M}, "
              f"cached cloud written {pd.Timestamp(sc, unit='s'):%Y-%m-%d %H:%M}")

    # ------------------------------------------------------ 1. the pairwise ties
    print("\n## 1. Every pairwise tie, shipped against extent-invariant\n")
    print("  The pairwise tie is the only thing that changes. The horizontal solution is\n"
          "  identical by construction (the Nuth & Kaeaeb aspect fit removes the median of\n"
          "  dh before fitting, so it never saw the vertical constant) and the network\n"
          "  weight is deliberately left at the horizontal fit's cell count, so the two\n"
          "  networks differ in one thing only. `pred` is the mechanism's own prediction\n"
          "  for `delta`, computed from this tile's sampling of this overlap.\n")
    rows = []
    for name in A.tiles:
        em, ei = sol[name]["out"]["overlap_median"][1], sol[name]["out"]["intercept"][1]
        for k in sorted(sol[name]["mech"]):
            mdt, cc, kk, med, n = sol[name]["mech"][k]
            rows.append([name, f"{k[0]}-{k[1]}", n,
                         f2(em[k][2] * 1000), f2(ei[k][2] * 1000),
                         f2((ei[k][2] - em[k][2]) * 1000),
                         f"{cc * 1000:+.1f}", f"{mdt:+.4f}", f2(-cc * mdt * 1000)])
    R.table(["tile", "pair", "cells", "dz_med", "dz_int", "delta", "c_mm", "mean_dtan",
             "pred"], rows)

    # ------------------------------- 2. extent-invariance, in measurement
    print("\n## 2. The elba / elbaext disagreement about the SAME flight lines\n")
    print("  Both tiles carry swaths 135-138, so their constants can be put on a COMMON\n"
          "  LINE (swath 135) and compared. Under the shipped tie they disagree; an\n"
          "  extent-dependent tie is what predicts that disagreement, so an extent-invariant\n"
          "  one must remove it.\n")
    e_old = reref({s: v[2] * 1000 for s, v in
                   sol["elba_fulldensity"]["out"]["overlap_median"][0].items()}, 135)
    x_old = reref({s: v[2] * 1000 for s, v in
                   sol["elbaext"]["out"]["overlap_median"][0].items()}, 135)
    e_new = reref({s: v[2] * 1000 for s, v in
                   sol["elba_fulldensity"]["out"]["intercept"][0].items()}, 135)
    x_new = reref({s: v[2] * 1000 for s, v in
                   sol["elbaext"]["out"]["intercept"][0].items()}, 135)
    rows = []
    for s in (136, 137, 138):
        do, dn = e_old[s] - x_old[s], e_new[s] - x_new[s]
        rows.append([str(s), f2(e_old[s]), f2(x_old[s]), f2(do),
                     f2(e_new[s]), f2(x_new[s]), f2(dn),
                     f"{do - dn:+.2f} ({100 * (1 - abs(dn) / abs(do)):.0f}%)"])
    R.table(["swath", "elba_old", "ext_old", "dis_old", "elba_new", "ext_new", "dis_new",
             "removed"], rows)
    print("\n  Re-referencing to 135 is itself a choice of origin, so the disagreements are "
          "only determined\n  up to a common constant. The origin-free statement is their spread "
          "over the four\n  shared lines, and their RMS after the best common constant is "
          "removed:\n")
    R.column("subset", "which shared flight lines the summary is taken over")
    R.column("spread_old", "max - min of the tile-to-tile disagreement over those lines, "
                           "shipped tie, mm (independent of the common line)")
    R.column("spread_new", "the same spread with the intercept tie, mm")
    R.column("rms_old", "RMS of the disagreement after removing its mean over those lines, "
                        "shipped tie, mm (independent of the common line)")
    R.column("rms_new", "the same RMS with the intercept tie, mm")
    rows = []
    for lab, ss in (("135-138 (all shared)", (135, 136, 137, 138)),
                    ("135-137 (138 dropped)", (135, 136, 137))):
        vo = np.array([e_old[s] - x_old[s] for s in ss])
        vn = np.array([e_new[s] - x_new[s] for s in ss])
        rows.append([lab, f"{np.ptp(vo):.2f}", f"{np.ptp(vn):.2f}",
                     f"{np.sqrt(((vo - vo.mean()) ** 2).mean()):.2f}",
                     f"{np.sqrt(((vn - vn.mean()) ** 2).mean()):.2f}"])
    R.table(["subset", "spread_old", "spread_new", "rms_old", "rms_new"], rows)
    print("\n  Dropping 138 is NOT a filter on the result -- both rows are shown, and the "
          "whole\n  point of the second is that everything left after the fix is swath 138.")

    # ------------------------------------------------- 3. the common line
    print("\n## 3. The common line: how each tile samples the line it is re-expressed "
          "against\n")
    scan = {}
    rows = []
    for name in A.tiles:
        tile = f"data/derived/{name}"
        t = pq.read_table(R.input(f"{tile}/beam_offset_table.parquet",
                                  role="per-return gen1 CSF ground offsets to the gen2 "
                                       "surface, slope-normal mm, with the four registration "
                                       "terms and point_source_id / scan_angle per return"),
                          columns=["point_source_id", "scan_angle", "in_grid", "cell",
                                   "d_mm_corr", "dz_swath_mm"])
        g = t["in_grid"].to_numpy().astype(bool)
        df = pd.DataFrame({k: t[k].to_numpy()[g] for k in
                           ("point_source_id", "scan_angle", "cell", "d_mm_corr",
                            "dz_swath_mm")})
        del t
        scan[name] = df
        for s, q in df.groupby("point_source_id"):
            sa = q.scan_angle.to_numpy(float)
            rows.append([name, str(int(s)), len(q), f2(sa.mean()), f"{sa.std():.2f}",
                         f2(np.percentile(sa, 5)), f2(np.percentile(sa, 95)),
                         f"{abs(sa.mean()):.2f}"])
    R.table(["tile", "swath", "returns", "scan_mean", "scan_sd", "scan_p5", "scan_p95",
             "asym"], rows)

    both = sorted(set(scan["elba_fulldensity"].point_source_id.unique())
                  & set(scan["elbaext"].point_source_id.unique()))
    worst_asym = {int(s): max(abs(scan[n][scan[n].point_source_id == s].scan_angle.mean())
                              for n in A.tiles) for s in both}
    common = int(min(worst_asym, key=worst_asym.get))
    R.param("common line", common, src="MINE",
            why=f"chosen by a stated rule on measured quantities, not by taste: of the lines "
                f"BOTH tiles carry ({', '.join(str(s) for s in both)}), the one whose worst "
                f"|mean scan angle| across the two tiles is smallest -- the line both tiles "
                f"sample most symmetrically about nadir. Its effect is measured below and is "
                f"a LEVEL shift only (tests/test_coreg.py proves the choice changes no "
                f"swath-to-swath difference), so it excludes no data and moves no gradient.")
    print(f"\n  lines carried by both tiles: {both}")
    print("  worst |mean scan| across the two tiles: "
          + ", ".join(f"{s}: {v:.2f} deg" for s, v in sorted(worst_asym.items())))
    print(f"  -> common line = swath {common}; each tile's shipped ZERO LINE is "
          + ", ".join(f"{n}: {sol[n]['zero_line']}" for n in A.tiles))
    print("\n  what the choice is worth: the point-count-weighted mean shift of the gen1\n"
          "  cloud, which is the level the tile's whole DoD inherits.\n")
    rows = []
    for name in A.tiles:
        w = scan[name].groupby("point_source_id").size()
        lev = lambda d: sum(d[int(s)] * n for s, n in w.items()) / w.sum()
        for tie in ("overlap_median", "intercept"):
            dz = {s: v[2] * 1000 for s, v in sol[name]["out"][tie][0].items()}
            for g_ in (sol[name]["zero_line"], common):
                rows.append([name, f"{tie}, line {g_}", f2(lev(reref(dz, g_)))])
            rows.append([name, f"{tie}, zero-mean (ref=None)",
                         f2(lev({s: v[2] * 1000 for s, v in
                                 sol[name]["zeromean"][tie].items()}))])
    R.table(["tile", "line", "level"], rows)
    print("\n  The line is a LEVEL choice and nothing else -- it changes no swath-to-swath\n"
          "  difference (tests/test_coreg.py). But the level is the DoD's datum: gen1 has no\n"
          "  vertical tie to gen2 anywhere in the pipeline (only the geoid constant), so the\n"
          "  tile's whole DoD sits wherever the pinned line's raw z sits.")

    # ------------------------------------------- 4. what it does to the DoD
    print("\n## 4. What the new constants do to the DoD\n")
    T = A.dod_tile
    name = os.path.basename(T)
    meta = sol[name]["meta"]
    RES = float(meta["res_m"])
    dod = np.load(R.input(f"{T}/{A.dod}",
                          role="the DoD in use: gen2 at the cover-dependent percentile minus "
                               "gen1's registered median, m, + = elevation rose"))
    NY, NX = dod.shape
    zf = np.load(R.input(f"{T}/z_after.npy",
                         role="gen2 gridded ground, the DoD reference surface, m"))
    mfill = ~np.isfinite(zf)
    if mfill.any():
        zf = zf[tuple(distance_transform_edt(mfill, return_distances=False,
                                             return_indices=True))]
    gy, gx = np.gradient(zf, RES)
    nn = np.sqrt(gx.ravel() ** 2 + gy.ravel() ** 2 + 1.0)

    df = scan[name]
    cell = df.cell.to_numpy()
    ps = df.point_source_id.to_numpy()
    d0 = df.d_mm_corr.to_numpy(float)
    base = pd.Series(d0).groupby(cell).median()

    def dod_with(extra_mm):
        """The DoD after adding ``extra_mm`` (slope-normal mm, per return) to gen1.

        The per-cell gen1 ground is re-reduced with the pipeline's own estimator -- the
        ``ground_q = 0.50`` quantile of the slope-normal residual per cell,
        ``pipeline.py:603`` -- so a per-swath constant propagates through the median of a
        mixed-line cell exactly as the pipeline propagates it, not as an average. The
        slope-normal change is returned to the vertical by the cell's own |n|, which is
        how ``registration.swath_alignment_term`` divided by it in the first place.
        """
        s1 = pd.Series(d0 + extra_mm).groupby(cell).median()
        out = dod.ravel().copy() * 1000.0
        i = base.index.values
        out[i] = out[i] - (s1.values - base.values) * nn[i]
        return out.reshape(NY, NX) / 1000.0

    dz_old = {s: v[2] for s, v in sol[name]["out"]["overlap_median"][0].items()}
    dz_new = {s: v[2] for s, v in sol[name]["out"]["intercept"][0].items()}
    d_tie = {s: dz_new[s] - dz_old[s] for s in dz_old}
    d_all = {s: reref(dz_new, common)[s] - reref(dz_old, sol[name]["zero_line"])[s] for s in dz_old}

    print(f"  per-swath constant, mm, on the shipped zero line ({sol[name]['zero_line']}):\n")
    R.table(["swath", "dz_old", "dz_new", "delta"],
            [[str(s), f2(dz_old[s] * 1000), f2(dz_new[s] * 1000), f2(d_tie[s] * 1000)]
             for s in sorted(dz_old)])
    print(f"\n  and re-referenced to the proposed common line (swath {common}):\n")
    R.table(["swath", "dz_old", "dz_new", "delta"],
            [[str(s), f2(reref(dz_old, common)[s] * 1000), f2(reref(dz_new, common)[s] * 1000),
              f2(d_all[s] * 1000)] for s in sorted(dz_old)])

    nnc = nn[cell]
    cases = {
        # the validation case: subtract the whole shipped per-swath alignment term back out.
        # STABLE_POINT_TILT_AUDIT.md section 5c reached the same raster by a different route
        # (per-(cell,line) medians differenced against the per-cell median) and reports
        # mean -36.16, dE -24.25, dN -15.11 mm/km at 250 m blocks. Reproducing that is the
        # check that this re-reduction is the pipeline's estimator and not a look-alike.
        "noswath": -df.dz_swath_mm.to_numpy(float),
        "tie": np.array([1000.0 * d_tie[int(s)] for s in ps]) / nnc,
        # NOTE: this key becomes the raster filename suffix, so it still reads
        # "gauge" -- the retired word -- in dod_cover_q2_tie_gauge.npy, which is on
        # disk and cited twice in analysis/SWATH_TIE_INTERCEPT.md. Renaming the key
        # renames the product and orphans those references, so that is a separate
        # decision, not a side effect of fixing the vocabulary.
        "tie_gauge": np.array([1000.0 * d_all[int(s)] for s in ps]) / nnc,
    }
    rows = []
    for tag, extra in cases.items():
        arr = dod_with(extra)
        p = os.path.join(T, A.dod.replace(".npy", f"_{tag}.npy"))
        np.save(p, arr)
        ch = (arr - dod).ravel() * 1000.0
        f = np.isfinite(ch)
        rows.append([os.path.basename(p), f2(np.median(ch[f])), f2(ch[f].min()),
                     f2(ch[f].max())])
    R.table(["raster", "md_mm", "lo_mm", "hi_mm"], rows)
    print("\n  The tilt and the stable-ground scatter are read off these rasters by\n"
          "  analysis/stable_point_tilt_audit.py, so before and after come from identical\n"
          "  code. `noswath` is the validation case, not a proposal:\n")
    for tag in cases:
        print(f"    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python "
              f"analysis/stable_point_tilt_audit.py --tile {T} "
              f"--dod {A.dod.replace('.npy', f'_{tag}.npy')}")

    R.done(headline="intercept tie changes elba's constants by "
                    + ", ".join(f"{1000 * d_tie[s]:+.1f}" for s in sorted(d_tie))
                    + " mm; the elba/elbaext disagreement goes "
                    + ", ".join(f2(e_old[s] - x_old[s]) for s in (136, 137, 138)) + " -> "
                    + ", ".join(f2(e_new[s] - x_new[s]) for s in (136, 137, 138)) + " mm")


if __name__ == "__main__":
    main()
