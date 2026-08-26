"""Aggregate the per-mark, per-line gen1 ties on lines 133-138 and test them against
the pipeline's own per-swath constants.

Reads ``more_marks.csv`` / ``more_marks_perline.csv`` (from ``gen1_more_marks_tie.py``)
and ``line_tracks.json`` (from ``gen1_line_tracks.py``).  Computes nothing from the
point cloud; every number here is an aggregate of those.

The one cut in this file is the ``--ladder-cut`` sweep, and it is a MEASUREMENT of what
cutting does, not a cut applied to the answer: the headline numbers use every mark.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                       # noqa: E402

SCRATCH = os.environ.get("SCRATCH", ".")
LINES = (133, 134, 135, 136, 137, 138)
CORR = "data/derived/elbaext/corrections_geoid.json"


def se(v):
    v = np.asarray(v, float)
    return float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else float("nan")


def main():
    R = Run("Does widening the 2008-control catchment on gen1's own flight lines "
            "133-138 enlarge the set and tighten gen1's datum?")
    p_marks = R.input(f"{SCRATCH}/more_marks.csv",
                      role="per-mark ties + siting screen, vendor class-2 ground, "
                           "no geoid shift, no lateral shift")
    p_pl = R.input(f"{SCRATCH}/more_marks_perline.csv",
                   role="per-mark-per-line ties, restricted to that line's own returns")
    p_tr = R.input(f"{SCRATCH}/line_tracks.json",
                   role="whole-line straight ground-track fits; used for TARGETING only")
    p_c = R.input(CORR, role="the pipeline's own per-swath internal alignment, the thing "
                             "the marks are tested against")
    R.param("lines", list(LINES), src="andy", why="Elba's own flight lines")
    R.param("cover_classes", ["L1O", "L5U"], src="repo",
            why="ADDITIONAL_GROUND_CONTROL.md section 7.1 / GEN1_OWN_CONTROL_TIE.md section 4: "
                "the vegetated classes measure canopy, not datum")
    R.param("catchment_radius_m", 2000.0, src="MINE",
            why="widened from the previous 500 m purely to MEASURE where the catchment "
                "ends; it selects candidates to READ, and the line assignment that "
                "decides membership is read off point_source_id, not off this")
    R.param("ground_source", "vendor_class2", src="andy", why="CSF is ~460 s/tile")
    R.param("geoid_shift_m", 0.0, src="repo",
            why="2008 control is NAVD88(GEOID03) and so is the raw gen1 cloud -- same "
                "frame, nothing to convert")
    R.param("swath_shift_m", (0.0, 0.0, 0.0), src="repo",
            why="no gen2-derived term in a gen1-against-its-own-control comparison")

    marks = pd.read_csv(p_marks)
    pl = pd.read_csv(p_pl)
    tracks = json.load(open(p_tr))["whole"]
    corr = json.load(open(p_c))["per_swath_internal_alignment_dxdydz_m"]

    on = marks.n_lines_ours > 0
    R.mask("on_lines_133_138", on.to_numpy(), of=len(marks),
           defn="the mark's vendor class-2 ground returns inside the 7.5 m report radius "
                "include at least one point_source_id in 133-138")

    R.column("line", "point_source_id of the returns the tie was estimated from (unitless id)")
    R.column("n_marks", "count of open/urban 2008 control marks carrying that line's ground returns")
    R.column("mean_mm", "mean over those marks of tie = surveyed - z_lidar, mm; "
                        "positive = gen1 reads BELOW the mark")
    R.column("sd_mm", "sample sd of tie over those marks (ddof=1)")
    R.column("SE_mm", "standard error of the mean of the per-mark ties on that line, sd/sqrt(n_marks), mm")
    R.column("median_mm", "median of tie over those marks")
    R.column("km_min", "smallest distance from the Elba reference point, km")
    R.column("km_max", "largest distance from the Elba reference point, km")
    R.column("R_m", "candidate catchment: |cross-track| to the nearest 133-138 centreline")
    R.column("n_cand", "count of open/urban control marks inside R_m")
    R.column("n_on_lines", "count of those candidates whose ground returns actually include a point_source_id in 133-138")
    R.column("n_tiles", "count of distinct gen1 tiles the candidates fall in")
    R.column("n_new_tiles", "count of those tiles not already on disk before this run")
    R.column("point", "MnDNR 2008 control point identifier (unitless label)")
    R.column("cls", "MnDNR land-cover class label: L1O open terrain, L5U urban (unitless)")
    R.column("km", "distance from the Elba reference point (579705.72, 4883677.71), km")
    R.column("d_near_m", "|cross-track| to the nearest of the six centrelines")
    R.column("line_proxy", "line id assigned by nearest fitted centreline -- the OLD rule (unitless id)")
    R.column("line_psid", "line id(s) actually present in the ground returns at the mark -- the NEW rule (unitless id)")
    R.column("n_report", "count of ground returns inside the 7.5 m report radius")
    R.column("tie_mm", "surveyed - z_lidar at the report radius, that line's returns only")
    R.column("sigma_mm", "half the tie spread across the pipeline-scale radii (2.5-10 m)")
    R.column("ladder_mm", "full tie spread across the pipeline-scale radii")
    R.column("spread_mm", "local p05-p95 of ground z within 5 m -- section 7.1 siting")
    R.column("ctrl_pct", "percentile of the surveyed height in that local z distribution")
    R.column("dnr_mm", "MnDNR's own Control Z - Surface Z for this mark, mm")
    R.column("cut_mm", "the sweep value: marks with ladder_mm <= this are kept, mm")
    R.column("sigma_site_mm", "pooled sd of tie about its own line's mean (ANOVA within-line)")
    R.column("sd_all_mm", "sd of tie over all kept marks, ignoring line")
    R.column("quantity", "name of the statistic stated in this row (unitless label)")
    R.column("value", "the value of that statistic, with its unit printed in the cell")
    R.column("note", "how that value was obtained (unitless text)")
    R.column("dz_pipeline_mm", "the vertical component of corrections_geoid.json per_swath_internal_alignment_dxdydz_m for this line, converted to mm")
    R.column("mark_mean_mm", "the mark-derived mean tie for this line, mm")
    R.banner()

    # ---------------------------------------------------------------- 1. the catchment
    rows = []
    on_disk_before = set()
    for t in sorted(set(marks.tile)):
        on_disk_before.add(t)
    prev = json.load(open(f"{SCRATCH}/tiles_before.json")) if os.path.exists(
        f"{SCRATCH}/tiles_before.json") else None
    for Rm in (250, 500, 700, 1000, 1500, 2000):
        s = marks[marks.d_min <= Rm]
        tiles = sorted(set(s.tile))
        new = [t for t in tiles if prev is not None and t not in prev]
        rows.append([Rm, len(s), int((s.n_lines_ours > 0).sum()), len(tiles),
                     len(new) if prev is not None else -1])
    R.table(["R_m", "n_cand", "n_on_lines", "n_tiles", "n_new_tiles"], rows)

    # ------------------------------------------------- 2. the line assignment, corrected
    good = marks[on].copy()
    good["lines_set"] = [tuple(sorted(int(float(v)) for v in str(t).split(",")))
                         for t in good.lines_ours]
    good["line_proxy_wrong"] = [tuple([int(a)]) != b
                                for a, b in zip(good.line_nearest, good.lines_set)]
    print(f"\nline assignment: {len(good)} marks carry 133-138 ground; the nearest-centreline "
          f"proxy disagrees with point_source_id on {int(good.line_proxy_wrong.sum())} of them")
    if good.line_proxy_wrong.any():
        print(good.loc[good.line_proxy_wrong,
                       ["point", "d_min", "line_nearest", "lines_ours"]].to_string(index=False))
    # the mixed-line tie and the single-line tie must coincide wherever only one line is
    # present -- that identity is what says the previous run's 18 numbers are reproduced
    m1 = good.n_lines_ours == 1
    j = good[m1].merge(pl, left_on="point", right_on="point", suffixes=("", "_pl"))
    dmax = float(np.abs(j.tie_mixed_mm - j.tie_mm).max())
    print(f"mixed-line tie vs single-line tie, marks with one line: max |difference| "
          f"{dmax:.6f} mm over {len(j)} marks")
    print(f"marks carrying MORE THAN ONE of lines 133-138 in the report radius: "
          f"{int((good.n_lines_ours > 1).sum())}")
    off = marks[~on]
    print(f"candidates carrying NO 133-138 ground: {len(off)}; their point_source_ids: "
          f"{sorted(set(off.psid_at_mark.dropna()))}")

    # ------------------------------------------------------------------ 3. screen table
    srows = []
    for _, r in good.sort_values(["lines_ours", "km"]).iterrows():
        srows.append([r.point, r.cls, round(r.km, 2), round(r.d_min, 1), int(r.line_nearest),
                      r.lines_ours, int(r.n_report), round(r.tie_mixed_mm, 1),
                      round(r.sigma_mm, 1), round(r.ladder_mm, 1), round(r.spread_mm, 0),
                      round(r.ctrl_pct, 0), round(r.dnr_mm, 0)])
    R.table(["point", "cls", "km", "d_near_m", "line_proxy", "line_psid", "n_report",
             "tie_mm", "sigma_mm", "ladder_mm", "spread_mm", "ctrl_pct", "dnr_mm"], srows)

    # ------------------------------------------------------------------- 4. by line
    lrows = []
    for ln in LINES:
        v = pl[pl.line == ln].tie_mm.to_numpy()
        k = pl[pl.line == ln]
        lrows.append([ln, v.size, round(float(v.mean()), 1),
                      round(float(v.std(ddof=1)), 1) if v.size > 1 else float("nan"),
                      round(se(v), 1), round(float(np.median(v)), 1),
                      round(float(k.km.min()), 1), round(float(k.km.max()), 1)])
    R.table(["line", "n_marks", "mean_mm", "sd_mm", "SE_mm", "median_mm", "km_min", "km_max"],
            lrows)

    # --------------------------------------------------------------- 5. combined, 2 ways
    lm = np.array([pl[pl.line == ln].tie_mm.mean() for ln in LINES])
    allv = pl.tie_mm.to_numpy()
    # ANOVA over lines
    grand = allv.mean()
    ssb = sum(pl[pl.line == ln].tie_mm.size * (pl[pl.line == ln].tie_mm.mean() - grand) ** 2
              for ln in LINES)
    ssw = sum(((pl[pl.line == ln].tie_mm - pl[pl.line == ln].tie_mm.mean()) ** 2).sum()
              for ln in LINES)
    dfb, dfw = len(LINES) - 1, allv.size - len(LINES)
    F = (ssb / dfb) / (ssw / dfw)
    from scipy import stats
    pF = float(stats.f.sf(F, dfb, dfw))
    sigma_site = float(np.sqrt(ssw / dfw))
    crows = [
        ["mean of the six line means", f"{lm.mean():+.1f} mm",
         f"SE over the six lines {se(lm):.1f} mm"],
        ["mean over marks, marks independent", f"{allv.mean():+.1f} mm",
         f"SE {se(allv):.1f} mm, n = {allv.size}"],
        ["median over marks", f"{np.median(allv):+.1f} mm", f"n = {allv.size}"],
        ["ANOVA over the six lines", f"F = {F:.2f}",
         f"p = {pF:.4f}, df = ({dfb}, {dfw})"],
        ["sigma_site (within-line, all marks)", f"{sigma_site:.1f} mm",
         "pooled sd of a mark about its own line's mean"],
        ["sd over marks, ignoring line", f"{allv.std(ddof=1):.1f} mm", f"n = {allv.size}"],
    ]
    R.table(["quantity", "value", "note"], crows)

    # --------------------------------------------------- 6. does screening reduce sigma?
    crow = []
    for cut in (15, 20, 25, 30, 40, 50, 75, 100, 1e9):
        s = pl[pl.ladder_mm <= cut]
        if s.line.nunique() < 2 or len(s) <= s.line.nunique():
            crow.append([cut if cut < 1e8 else "no cut", len(s), s.line.nunique(),
                         float("nan"), float("nan"), float("nan")])
            continue
        g = grand = s.tie_mm.mean()
        w = sum(((s[s.line == ln].tie_mm - s[s.line == ln].tie_mm.mean()) ** 2).sum()
                for ln in s.line.unique())
        dw = len(s) - s.line.nunique()
        crow.append([cut if cut < 1e8 else "no cut", len(s), s.line.nunique(),
                     round(float(np.sqrt(w / dw)), 1), round(float(s.tie_mm.std(ddof=1)), 1),
                     round(float(s.tie_mm.mean()), 1)])
    R.column("n_kept", "count of mark-line ties surviving the cut")
    R.column("n_lines_kept", "count of distinct flight lines surviving the cut")
    R.column("mean_mm_kept", "mean of tie over the surviving mark-line ties, mm")
    R.table(["cut_mm", "n_kept", "n_lines_kept", "sigma_site_mm", "sd_all_mm", "mean_mm_kept"],
            crow)

    # ------------------------------------------------- 7. against the pipeline constants
    dz = {int(k): 1000.0 * v[2] for k, v in corr.items()}
    mk = {ln: float(pl[pl.line == ln].tie_mm.mean()) for ln in LINES}
    prows = [[ln, round(dz[ln], 1), round(mk[ln], 1)] for ln in LINES]
    R.table(["line", "dz_pipeline_mm", "mark_mean_mm"], prows)
    a = np.array([dz[ln] for ln in LINES]); b = np.array([mk[ln] for ln in LINES])
    r6, p6 = stats.pearsonr(a, b)
    idx = [ln for ln in LINES if ln != 135]
    a5 = np.array([dz[ln] - dz[135] for ln in idx])
    b5 = np.array([mk[ln] - mk[135] for ln in idx])
    r5, p5 = stats.pearsonr(a5, b5)
    idx3 = [ln for ln in LINES if ln != 133]
    a3 = np.array([dz[ln] - dz[133] for ln in idx3])
    b3 = np.array([mk[ln] - mk[133] for ln in idx3])
    r3, p3 = stats.pearsonr(a3, b3)
    rms5 = float(np.sqrt(np.mean((b5 - a5) ** 2)))
    rms6 = float(np.sqrt(np.mean(((b - b.mean()) - (a - a.mean())) ** 2)))
    R.table(["quantity", "value", "note"], [
        ["Pearson r, all six lines, raw", f"{r6:+.3f}",
         f"p = {p6:.3f}, n = 6; correlation is invariant to the reference line"],
        ["Pearson r, relative to line 135", f"{r5:+.3f}", f"p = {p5:.3f}, n = 5"],
        ["Pearson r, relative to line 133", f"{r3:+.3f}",
         f"p = {p3:.3f}, n = 5; 133 is the pipeline's own reference (dz = 0)"],
        ["RMS residual, relative to 135", f"{rms5:.1f} mm", "mark minus pipeline, n = 5"],
        ["RMS residual, both mean-removed", f"{rms6:.1f} mm", "n = 6"],
        ["spread of the pipeline dz", f"{a.max() - a.min():.1f} mm", "max - min over six lines"],
        ["spread of the mark line means", f"{b.max() - b.min():.1f} mm", "max - min over six lines"],
        ["SE of a single line mean, median", f"{np.nanmedian([se(pl[pl.line == ln].tie_mm) for ln in LINES]):.1f} mm",
         "over the five lines with more than one mark; line 136 has one mark and no SE"],
    ])

    # sign / reference variants, because the brief's r = -0.605 does not reproduce
    print("\nreference and sign variants of the pipeline test "
          "(mark line means vs per_swath_internal_alignment dz):")
    print(f"{'ref':>5} {'sign':>5} {'n':>3} {'r':>8} {'p':>8} {'rms_mm':>8}")
    for ref in list(LINES) + [None]:
        for sgn in (+1, -1):
            ii = [ln for ln in LINES if ln != ref]
            if ref is None:
                aa = np.array([dz[ln] for ln in LINES]); bb = sgn * np.array([mk[ln] for ln in LINES])
                aa = aa - aa.mean(); bb = bb - bb.mean()
            else:
                aa = np.array([dz[ln] - dz[ref] for ln in ii])
                bb = sgn * np.array([mk[ln] - mk[ref] for ln in ii])
            rr, pp = stats.pearsonr(aa, bb)
            print(f"{str(ref):>5} {sgn:>+5d} {len(aa):>3d} {rr:>+8.3f} {pp:>8.3f} "
                  f"{np.sqrt(np.mean((bb - aa) ** 2)):>8.1f}")

    # the 43-mark screen this widening started from, re-aggregated the same way, so the
    # "organised by flight line" claim and sigma_site can be compared like with like
    prev = pd.read_csv(f"{SCRATCH}/screen_results.csv")
    prev = prev[prev.status == "ok"]
    gl = prev.groupby("line").tie_mm
    grand2 = prev.tie_mm.mean()
    ssb2 = sum(g.size * (g.mean() - grand2) ** 2 for _, g in gl)
    ssw2 = sum(((g - g.mean()) ** 2).sum() for _, g in gl)
    dfb2, dfw2 = prev.line.nunique() - 1, len(prev) - prev.line.nunique()
    F2 = (ssb2 / dfb2) / (ssw2 / dfw2)
    print(f"\nthe 43-mark screen (lines 128-149, nearest-centreline assignment, as run "
          f"before this session):")
    print(f"  ANOVA over {prev.line.nunique()} lines: F = {F2:.2f}, "
          f"p = {stats.f.sf(F2, dfb2, dfw2):.2e}, df = ({dfb2}, {dfw2})")
    print(f"  sigma_site (within-line pooled sd) = {np.sqrt(ssw2 / dfw2):.1f} mm; "
          f"sd ignoring line = {prev.tie_mm.std(ddof=1):.1f} mm")
    print(f"  same statistics restricted to lines 133-138 in THAT table: ", end="")
    q = prev[prev.line.isin(LINES)]
    print(f"n = {len(q)}, mean = {q.tie_mm.mean():+.1f} mm, "
          f"sd = {q.tie_mm.std(ddof=1):.1f} mm")

    R.done(headline=f"18 marks, 6 lines, mean of line means {lm.mean():+.1f} +/- {se(lm):.1f} mm; "
                    f"widening the catchment from 500 m to 2000 m added 0 marks")


if __name__ == "__main__":
    main()
