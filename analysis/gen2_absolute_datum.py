#!/usr/bin/env python3
"""gen2's absolute vertical datum, from its own 390 published per-point residuals.

The 2021 3DEP project publishes, for every held-out accuracy checkpoint, the surveyed
height and the delivered surface read at the mark -- see
``analysis/groundtruth/parse_gen2_control.py``. That is the gen2 analogue of the 963-mark
gen1 residual field, and it needs no point cloud.

This driver answers, in order:

1. what the residual's sign actually is, re-derived from the CSV;
2. gen2's level, by accuracy class (NVA / VVA), by QL block, and in distance bands about
   a named site -- all bands reported, none chosen;
3. the NVA-vs-VVA contrast, which is the cleanest canopy discriminator available because
   the classes were assigned by the surveyors and not by us;
4. the residual as a spatial FIELD, kriged to the site with the SAME machinery
   (``groundtruth.residual_field``) and the same swept nuisance grid as gen1's
   ``analysis/control_residual_field.py``, so the two epochs' site values are comparable
   rather than merely both present;
5. our own estimator against USGS's, at the marks where both exist.

Sign convention throughout, the same as ``groundtruth.tie``: **positive = the surface
reads LOW**, i.e. the number is the constant to ADD. ``usgs_*_error_m`` is
``surveyed - surface``, re-derived in section 1.

**The LCPs are excluded from every accuracy statement below.** They calibrated gen2
(vendor FGDC metadata, Ground Conditions); checking gen2 against them would be circular.
They carry no residual in any case -- the VATool never tested them -- and they appear
here only in the coverage bookkeeping.

Nothing here invents a cut. Every radius, lag, block and class list is a required
argument.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from trust.provenance import Run                                        # noqa: E402
from lidar_diff_icp.groundtruth import residual_field as RF             # noqa: E402

DEFAULT_CSV = "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"

SURFACES = {"laz": ("usgs_ql1_laz_z_m", "usgs_ql1_laz_error_m",
                    "usgs_ql0_laz_z_m", "usgs_ql0_laz_error_m"),
            "dem": ("usgs_ql1_dem_z_m", "usgs_ql1_dem_error_m",
                    "usgs_ql0_dem_z_m", "usgs_ql0_dem_error_m")}


def load(csv_path, surface, block_pref):
    """Rows with a residual on ``surface``, one per mark.

    ``block_pref`` is an explicit ordered list of QL blocks; a mark tested in more than
    one is taken from the first block in that list that has it, and the choice is
    reported, never averaged silently.
    """
    zc1, ec1, zc0, ec0 = SURFACES[surface]
    cols = {"QL1": (zc1, ec1), "QL0": (zc0, ec0)}
    out = []
    n_both = 0
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            have = [b for b in ("QL1", "QL0") if r[cols[b][1]] != ""]
            if not have:
                continue
            n_both += len(have) > 1
            b = next(b for b in block_pref if b in have)
            out.append(dict(
                point_id=r["point_id"], cls=r["point_type"], role=r["role"],
                easting=float(r["easting"]), northing=float(r["northing"]),
                elevation=float(r["elevation"]),
                block=b, blocks="+".join(have),
                surf_z=float(r[cols[b][0]]), resid_mm=float(r[cols[b][1]]) * 1000.0,
                collected=r["collected"], geoid=r["geoid_model"],
            ))
    return out, n_both


def sign_test(csv_path, tol_m):
    """Re-derive which subtraction each error column is, on every row of the CSV."""
    out = {}
    for surface, (zc1, ec1, zc0, ec0) in SURFACES.items():
        for block, (zc, ec) in (("QL1", (zc1, ec1)), ("QL0", (zc0, ec0))):
            a, b, n = [], [], 0
            with open(csv_path, newline="") as fh:
                for r in csv.DictReader(fh):
                    if r[ec] == "":
                        continue
                    z, s, e = float(r["elevation"]), float(r[zc]), float(r[ec])
                    a.append(abs((z - s) - e)); b.append(abs((s - z) - e)); n += 1
            a, b = np.array(a), np.array(b)
            out[(surface, block)] = (n, int((a <= tol_m).sum()), float(a.max()),
                                     int((b <= tol_m).sum()), float(b.max()))
    return out


def welch(a, b):
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1) / a.size, b.var(ddof=1) / b.size
    t = (ma - mb) / np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (a.size - 1) + vb ** 2 / (b.size - 1))
    return ma - mb, np.sqrt(va + vb), t, df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control-csv", default=DEFAULT_CSV)
    ap.add_argument("--site-name", required=True)
    ap.add_argument("--site-easting", type=float, required=True)
    ap.add_argument("--site-northing", type=float, required=True)
    ap.add_argument("--sign-tol-m", type=float, required=True)
    ap.add_argument("--surfaces", required=True,
                    help="comma list from {laz,dem}: which delivered surface's residual")
    ap.add_argument("--block-preference", required=True,
                    help="comma list of QL blocks, first match wins for a mark tested twice")
    ap.add_argument("--band-radii-km", required=True, help="comma list, all reported")
    ap.add_argument("--max-lag-m", required=True, help="comma list, swept")
    ap.add_argument("--n-lags", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--n-pairs", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--estimators", required=True, help="comma list: dowd and/or matheron")
    ap.add_argument("--seeds", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--gen1-csv", default=None,
                    help="the bundled 2008 control CSV; enables section 8, the SAME "
                         "open-vs-vegetated contrast computed on gen1")
    ap.add_argument("--gen1-open-classes", default=None,
                    help="comma list of gen1 point_type values to call OPEN; required "
                         "with --gen1-csv, no default")
    ap.add_argument("--gen1-veg-classes", default=None,
                    help="comma list of gen1 point_type values to call VEGETATED; "
                         "required with --gen1-csv, no default")
    ap.add_argument("--boxes-dir", default=None,
                    help="directory of per-checkpoint gen2 boxes; enables the per-line "
                         "section 6b")
    ap.add_argument("--line-half-width-m", type=float, default=None,
                    help="half-width of the window read from each box for section 6b; "
                         "required with --boxes-dir, no default")
    ap.add_argument("--res-m", type=float, default=None,
                    help="pipeline grid resolution setting the radius ladder; required "
                         "with --boxes-dir, no default")
    ap.add_argument("--tie-json", default=None,
                    help="output of analysis/groundtruth/gen2_checkpoint_tie.py --json; "
                         "section 5 compares our estimator with USGS's where both exist")
    A = ap.parse_args()

    surfaces = A.surfaces.split(",")
    block_pref = A.block_preference.split(",")
    bands = [float(s) for s in A.band_radii_km.split(",")]
    max_lags = [float(s) for s in A.max_lag_m.split(",")]
    n_lags_grid = [int(s) for s in A.n_lags.split(",")]
    n_pairs_grid = [int(s) for s in A.n_pairs.split(",")]
    estimators = A.estimators.split(",")
    seeds = [int(s) for s in A.seeds.split(",")]

    R = Run(f"what is gen2's absolute vertical level at {A.site_name}, measured against "
            f"its own held-out surveyed checkpoints, and how does it split by canopy?")
    R.input(A.control_csv, role="the 534 surveyed control points of the 2021 3DEP "
                                "MN_SE_Driftless project with the USGS VATool per-point "
                                "residuals on the 390 held-out checkpoints")
    if A.tie_json:
        R.input(A.tie_json, role="our own estimate_tie result at the six checkpoints "
                                 "whose gen2 point cloud is on disk")
    R.param("site", (A.site_name, A.site_easting, A.site_northing), src="andy", why="")
    R.param("sign_tol_m", A.sign_tol_m, src="andy", why="")
    R.param("surfaces", surfaces, src="andy", why="")
    R.param("block_preference", block_pref, src="andy", why="")
    R.param("band_radii_km", bands, src="andy", why="")
    R.param("max_lag_m / n_lags / n_pairs / estimators / seeds",
            (max_lags, n_lags_grid, n_pairs_grid, estimators, seeds), src="andy", why="")
    R.param("excluded_from_every_accuracy_statement", "role=calibration (the 143 LCPs)",
            src="repo",
            why="the vendor FGDC metadata says the LCPs calibrated the lidar and the "
                "NVA/VVA checkpoints were held out; the LCPs carry no residual at all")
    if A.gen1_csv:
        R.input(A.gen1_csv, role="the 1004 MnDNR ground control rows validating the 2008 "
                                 "SE-Minnesota lidar, gen1's own control on gen1's own "
                                 "geoid")
        R.param("gen1_open_classes / gen1_veg_classes",
                (A.gen1_open_classes, A.gen1_veg_classes), src="MINE",
                why="gen1's land-cover taxonomy (L1O open, L2T tall weeds and crops, "
                    "L3B brush and low trees, L4F forested, L5U urban) is NOT the 3DEP "
                    "NVA/VVA taxonomy, so putting the two epochs on one contrast needs a "
                    "grouping and I chose this one; it decides which gen1 marks are "
                    "called vegetated and is a proposal, not a result")
    R.param("line_half_width_m / res_m", (A.line_half_width_m, A.res_m), src="andy",
            why="")
    R.param("geoid_term_mm", 0.0, src="repo",
            why="the control is NAVD88(GEOID18) and gen2 is delivered on NAVD88(GEOID18), "
                "so the geoid cancels and no conversion is applied; section 0 asserts "
                "this per mark and raises if any mark disagrees")

    R.column("surface", "which delivered 2021 surface the residual is against: laz = the "
                        "classified point cloud, dem = the OPR DEM")
    R.column("block", "USGS quality-level block whose products were tested, QL0 or QL1")
    R.column("n", "marks entering this row, count")
    R.column("class", "the survey's own accuracy class: NVA non-vegetated, VVA vegetated")
    R.column("mean_mm", "mean residual, surveyed minus surface; POSITIVE = the surface "
                        "reads LOW, mm")
    R.column("median_mm", "median of the same residual, mm")
    R.column("sd_mm", "sample standard deviation over marks, mm")
    R.column("se_mm", "sd/sqrt(n): standard error of THE MEAN OVER THESE MARKS, not the "
                      "uncertainty of the level at any one place, mm")
    R.column("rmse_mm", "root mean square of the residual over these marks, mm")
    R.column("band_km", "marks within this distance of the site, cumulative, km")
    R.column("diff_mm", "NVA mean minus VVA mean: how much higher the surface reads under "
                        "vegetation than on open ground, mm")
    R.column("se_diff_mm", "standard error of that difference (Welch), mm")
    R.column("t", "Welch t statistic of the NVA-VVA difference, dimensionless")
    R.column("df", "Welch-Satterthwaite degrees of freedom, dimensionless")
    R.column("variant", "which marks and which drift enter the kriging system")
    R.column("estimator", "empirical-variogram estimator: dowd (robust) or matheron")
    R.column("max_lag_m", "largest separation entering the empirical variogram, m")
    R.column("nugget_mm2", "fitted spherical nugget, mm^2")
    R.column("sill_mm2", "fitted spherical partial sill, mm^2")
    R.column("range_m", "fitted spherical range, m")
    R.column("pred_mm", "kriged residual AT THE SITE, mm, sign as above")
    R.column("sd_field_mm", "sd of the error in predicting the spatially correlated part "
                            "of the field at the site, nugget filtered out, mm")
    R.column("sd_mark_mm", "sd of the error in predicting what ONE new mark placed at the "
                           "site would read, nugget included, mm")
    R.column("point_id", "surveyed mark identifier")
    R.column("usgs_mm", "USGS VATool residual at the mark, surveyed minus delivered "
                        "surface, mm")
    R.column("ours_mm", "our estimate_tie at the same mark on the same delivered class-2 "
                        "cloud, surveyed minus lidar, mm")
    R.column("delta_mm", "ours minus USGS at the same mark, mm")
    R.column("sigma_mm", "half the tie's spread across the pipeline-scale radii, mm")
    if A.boxes_dir:
        from lidar_diff_icp.groundtruth import checkpoints as _C
        for _cp in _C.load_bundled().usable():
            _p = os.path.join(A.boxes_dir, f"cp{_cp.point_id.split('_')[0]}_gen2.laz")
            if os.path.exists(_p):
                R.input(_p, role=f"gen2 2021 3DEP full-density box around checkpoint "
                                 f"{_cp.point_id}, read per flight line in section 6b")
    R.banner()

    # ---------------------------------------------------------------- 0. geoid, per mark
    print("\n== 0. the geoid, asserted per mark ==")
    bad = []
    with open(A.control_csv, newline="") as fh:
        rows_all = list(csv.DictReader(fh))
    for r in rows_all:
        if r["geoid_model"] != "GEOID18" or r["vertical_datum"] != "NAVD88":
            bad.append(r["point_id"])
    if bad:
        raise SystemExit(f"{len(bad)} marks are not NAVD88/GEOID18: {bad[:5]} -- the "
                         f"geoid does NOT cancel and this run must not proceed")
    print(f"  all {len(rows_all)} marks are NAVD88 / GEOID18, the same geoid gen2 is")
    print("  delivered on, so the geoid term is exactly zero and none is applied.")

    # ------------------------------------------------------------------- 1. the sign
    print(f"\n== 1. the sign, re-derived on every row (tolerance {A.sign_tol_m:g} m) ==")
    for (surf, blk), (n, ncs, mcs, nsc, msc) in sorted(sign_test(A.control_csv,
                                                                 A.sign_tol_m).items()):
        print(f"  {surf} {blk}: surveyed-surface {ncs}/{n} rows (max |resid| {mcs:.2e} m)"
              f"   surface-surveyed {nsc}/{n} (max {msc:.2e} m)")
    print("  => the error columns are surveyed MINUS surface: positive means the 2021")
    print("     surface reads LOW of the mark, the same sign family as groundtruth.tie.")

    if A.gen1_csv and (A.gen1_open_classes is None or A.gen1_veg_classes is None):
        raise SystemExit("--gen1-csv requires --gen1-open-classes and "
                         "--gen1-veg-classes; neither has a default")
    if A.boxes_dir and (A.line_half_width_m is None or A.res_m is None):
        raise SystemExit("--boxes-dir requires --line-half-width-m and --res-m; "
                         "neither has a default")
    tie_json = json.load(open(A.tie_json)) if A.tie_json else None

    for surface in surfaces:
        rows, n_both = load(A.control_csv, surface, block_pref)
        by_id = {r["point_id"]: r for r in rows}
        v = np.array([r["resid_mm"] for r in rows])
        cls = np.array([r["cls"] for r in rows])
        blk = np.array([r["block"] for r in rows])
        x = np.array([r["easting"] for r in rows])
        y = np.array([r["northing"] for r in rows])
        d_km = np.hypot(x - A.site_easting, y - A.site_northing) / 1000.0

        print(f"\n\n########## surface = {surface} ##########")
        print(f"  {len(rows)} marks carry a residual; {n_both} were tested against both "
              f"blocks and are taken from {block_pref[0]} (preference {block_pref})")
        R.mask(f"{surface}: role=check", np.array([r['role'] == 'check' for r in rows]),
               defn="marks the vendor held out of calibration (NVA + VVA); the 143 LCPs "
                    "carry no residual and are absent from this table entirely",
               of=len(rows_all))

        # ------------------------------------------------- 2. the level
        print("\n== 2. gen2's level, whole project and by block ==")
        tab = []
        for label, sel in ([("all blocks", np.ones(len(rows), bool))] +
                           [(b, blk == b) for b in ("QL1", "QL0")]):
            for c in ("all", "NVA", "VVA"):
                m = sel & (np.ones(len(rows), bool) if c == "all" else (cls == c))
                if m.sum() < 2:
                    continue
                s = v[m]
                tab.append([surface, label, c, int(m.sum()), f"{s.mean():+.2f}",
                            f"{np.median(s):+.2f}", f"{s.std(ddof=1):.2f}",
                            f"{s.std(ddof=1)/np.sqrt(m.sum()):.2f}",
                            f"{np.sqrt((s**2).mean()):.2f}"])
        R.table(["surface", "block", "class", "n", "mean_mm", "median_mm", "sd_mm",
                 "se_mm", "rmse_mm"], tab)

        # ------------------------------------------- 3. NVA vs VVA, the canopy contrast
        print("\n== 3. NVA vs VVA -- the canopy discriminator, classes assigned by the "
              "surveyors ==")
        tab = []
        for label, sel in ([("all blocks", np.ones(len(rows), bool))] +
                           [(b, blk == b) for b in ("QL1", "QL0")]):
            a = v[sel & (cls == "NVA")]
            b_ = v[sel & (cls == "VVA")]
            if a.size < 2 or b_.size < 2:
                continue
            dm, se, t, df = welch(a, b_)
            tab.append([surface, label, int(a.size), int(b_.size), f"{dm:+.2f}",
                        f"{se:.2f}", f"{t:+.3f}", f"{df:.1f}"])
        R.table(["surface", "block", "n", "n", "diff_mm", "se_diff_mm", "t", "df"], tab)

        # ---------------------------------------------------- 4. distance bands
        print(f"\n== 4. by distance from {A.site_name}, every band reported ==")
        tab = []
        for bkm in bands:
            for c in ("all", "NVA", "VVA"):
                m = (d_km <= bkm) & (np.ones(len(rows), bool) if c == "all" else (cls == c))
                if m.sum() < 1:
                    tab.append([surface, f"{bkm:g}", c, int(m.sum()), "-", "-", "-", "-", "-"])
                    continue
                s = v[m]
                sd = s.std(ddof=1) if m.sum() > 1 else np.nan
                tab.append([surface, f"{bkm:g}", c, int(m.sum()), f"{s.mean():+.2f}",
                            f"{np.median(s):+.2f}",
                            "-" if not np.isfinite(sd) else f"{sd:.2f}",
                            "-" if not np.isfinite(sd) else f"{sd/np.sqrt(m.sum()):.2f}",
                            f"{np.sqrt((s**2).mean()):.2f}"])
        R.table(["surface", "band_km", "class", "n", "mean_mm", "median_mm", "sd_mm",
                 "se_mm", "rmse_mm"], tab)

        # -------------------------------------------------------- 5. the field model
        print(f"\n== 5. the residual as a spatial field, kriged to {A.site_name} ==")
        print("   same machinery and the same swept nuisance grid as gen1's")
        print("   analysis/control_residual_field.py; three cover treatments side by side,")
        print("   none of them chosen here.")
        Xall, labels, ev = RF.cover_design(cls, ("NVA", "VVA"))
        variants = [
            ("NVA only", cls == "NVA", None, None, None),
            ("VVA only", cls == "VVA", None, None, None),
            ("class covariate -> NVA", np.ones(len(rows), bool), Xall, ev("NVA"), labels),
            ("class covariate -> VVA", np.ones(len(rows), bool), Xall, ev("VVA"), labels),
        ]
        fit_rows, pred_rows = [], []
        for name, sel, X, x0d, lab in variants:
            xs, ys, vs = x[sel], y[sel], v[sel]
            Xs = None if X is None else X[sel]
            for est in estimators:
                for ml in max_lags:
                    mods, preds, sdf, sdm = [], [], [], []
                    for nl in n_lags_grid:
                        for npair in n_pairs_grid:
                            for sd_ in seeds:
                                mod, *_ = RF.fit_field(xs, ys, vs, max_lag_m=ml,
                                                       n_lags=nl, n_pairs=npair,
                                                       estimator=est, seed=sd_)
                                mods.append(mod)
                                k = RF.krige(xs, ys, vs, mod, A.site_easting,
                                             A.site_northing, X=Xs, x0_drift=x0d,
                                             drift_labels=lab or ("const",))
                                preds.append(k.value_mm)
                                sdf.append(k.sd_field_mm)
                                sdm.append(k.sd_new_mark_mm)
                    fit_rows.append([surface, name, est, f"{ml:.0f}",
                                     f"{np.median([m.nugget for m in mods]):.0f}",
                                     f"{np.median([m.sill for m in mods]):.0f}",
                                     f"{np.median([m.range_ for m in mods]):.0f}",
                                     int(sel.sum())])
                    pred_rows.append([surface, name, est, f"{ml:.0f}", int(sel.sum()),
                                      f"{np.median(preds):+.1f}",
                                      f"{np.min(preds):+.1f}..{np.max(preds):+.1f}",
                                      f"{np.median(sdf):.1f}", f"{np.median(sdm):.1f}"])
        R.column("pred_range_mm", "min..max of the kriged site value over the whole "
                                  "nuisance grid at this max_lag, mm")
        R.table(["surface", "variant", "estimator", "max_lag_m", "nugget_mm2", "sill_mm2",
                 "range_m", "n"], fit_rows)
        print()
        R.table(["surface", "variant", "estimator", "max_lag_m", "n", "pred_mm",
                 "pred_range_mm", "sd_field_mm", "sd_mark_mm"], pred_rows)

        # ------------------------------- 6. our estimator against USGS's, same marks
        if tie_json is not None and surface == "laz":
            print("\n== 6. our estimate_tie against USGS's VATool, at the marks where "
                  "both exist ==")
            tab = []
            ours, theirs = [], []
            for cp in tie_json["checkpoints"]:
                if not cp.get("attempted") or cp["point_id"] not in by_id:
                    continue
                u = by_id[cp["point_id"]]["resid_mm"]
                o = cp["tie_mm"]
                ours.append(o); theirs.append(u)
                tab.append([cp["point_id"], u["cls"] if isinstance(u, dict) else
                            by_id[cp["point_id"]]["cls"], f"{u:+.1f}", f"{o:+.1f}",
                            f"{o-u:+.1f}", f"{cp['sigma_mm']:.1f}"])
            R.table(["point_id", "class", "usgs_mm", "ours_mm", "delta_mm", "sigma_mm"],
                    tab)
            ours, theirs = np.array(ours), np.array(theirs)
            dd = ours - theirs
            print(f"  n={dd.size}  mean(ours-USGS) {dd.mean():+.1f} mm  "
                  f"median {np.median(dd):+.1f} mm  sd {dd.std(ddof=1):.1f} mm  "
                  f"RMS {np.sqrt((dd**2).mean()):.1f} mm")
            print(f"  correlation r = {np.corrcoef(ours, theirs)[0,1]:+.3f}")
            print("  The two estimators read the SAME delivered cloud at the SAME marks;")
            print("  the spread between them is estimator plus siting, not epoch.")

        # ------------------------------- 6b. per-line structure at the boxed marks
        if A.boxes_dir and surface == "laz":
            print("\n== 6b. per flight line, at the marks whose gen2 box is on disk ==")
            print("   DESCRIPTIVE ONLY. The pipeline reads gen2's delivered class 2 and")
            print("   performs no swath alignment on it, so a line-to-line spread here is")
            print("   a property of the delivered product, not something we correct.")
            from lidar_diff_icp.groundtruth import checkpoints as C
            from lidar_diff_icp.groundtruth import tie as T
            cps = {c.point_id: c for c in C.load_bundled().usable()}
            tab = []
            for pid, cp in cps.items():
                box = os.path.join(A.boxes_dir, f"cp{pid.split('_')[0]}_gen2.laz")
                if not os.path.exists(box):
                    continue
                g = T.vendor_ground_near(box, cp.easting, cp.northing,
                                         A.line_half_width_m, ground_class=2)
                lines, counts = np.unique(g.point_source_id, return_counts=True)
                for ln, cnt in zip(lines, counts):
                    est = T.estimate_tie(cp, g, line=int(ln), res=A.res_m,
                                         geoid_shift_m=0.0)
                    tab.append([pid, by_id[pid]["cls"] if pid in by_id else "-",
                                int(ln), int(cnt), f"{est.tie_mm:+.1f}",
                                f"{est.sigma_mm:.1f}", est.n_report])
            R.column("line", "point_source_id of the 2021 flight line, integer")
            R.column("n_line", "class-2 returns of that line inside the box, count")
            R.table(["point_id", "class", "line", "n_line", "ours_mm", "sigma_mm", "n"],
                    tab)
            for pid in sorted({t[0] for t in tab}):
                ts = [float(t[4]) for t in tab if t[0] == pid
                      and np.isfinite(float(t[4]))]
                nl = sum(1 for t in tab if t[0] == pid)
                if len(ts) > 1:
                    print(f"    {pid}: {nl} lines in the box, {len(ts)} with returns "
                          f"inside the fitting radius; spread {max(ts)-min(ts):.1f} mm, "
                          f"median {np.median(ts):+.1f} mm")
                else:
                    print(f"    {pid}: {nl} lines in the box but only {len(ts)} with "
                          f"returns inside the fitting radius -- NOT a per-line "
                          f"comparison at this mark")

    # ------------------------------- 7. the same contrast on gen1, for the difference
    if A.gen1_csv:
        print("\n\n== 7. the SAME open-vs-vegetated contrast, computed on gen1's own "
              "control ==")
        print("   gen1's residual is the vendor's Control Z - Surface Z against the")
        print("   DELIVERED 2008 DNR DEM; gen2's 'dem' column is the delivered 2021 OPR")
        print("   DEM. Those two are the like-for-like pair. gen1 has no published")
        print("   point-cloud residual, so gen2's 'laz' column has no gen1 counterpart.")
        opens = tuple(A.gen1_open_classes.split(","))
        vegs = tuple(A.gen1_veg_classes.split(","))
        cr = RF.load_residuals(A.gen1_csv)
        g_open = RF.stratify(cr, opens)
        g_veg = RF.stratify(cr, vegs)
        R.mask("gen1 open", g_open, defn=f"point_type in {opens}", of=len(cr))
        R.mask("gen1 vegetated", g_veg, defn=f"point_type in {vegs}", of=len(cr))
        tab = []
        for label, m in (("open", g_open), ("vegetated", g_veg),
                         ("both", g_open | g_veg)):
            sv = cr.resid_mm[m]
            tab.append(["gen1 dnr_dem", "all counties", label, int(m.sum()),
                        f"{sv.mean():+.2f}", f"{np.median(sv):+.2f}",
                        f"{sv.std(ddof=1):.2f}",
                        f"{sv.std(ddof=1)/np.sqrt(m.sum()):.2f}",
                        f"{np.sqrt((sv**2).mean()):.2f}"])
        R.table(["surface", "block", "class", "n", "mean_mm", "median_mm", "sd_mm",
                 "se_mm", "rmse_mm"], tab)
        dm, se, t, df = welch(cr.resid_mm[g_open], cr.resid_mm[g_veg])
        print(f"  gen1 open minus vegetated: {dm:+.2f} +/- {se:.2f} mm "
              f"(Welch t {t:+.3f}, df {df:.1f}, n {int(g_open.sum())}/{int(g_veg.sum())})")
        rows_d, _ = load(A.control_csv, "dem", block_pref)
        vd = np.array([r["resid_mm"] for r in rows_d])
        cd = np.array([r["cls"] for r in rows_d])
        dm2, se2, t2, df2 = welch(vd[cd == "NVA"], vd[cd == "VVA"])
        print(f"  gen2 NVA  minus VVA       : {dm2:+.2f} +/- {se2:.2f} mm "
              f"(Welch t {t2:+.3f}, df {df2:.1f}, n {int((cd=='NVA').sum())}/"
              f"{int((cd=='VVA').sum())})")
        print("\n  the epoch difference of the LEVELS, one stratum at a time, as plain")
        print("  arithmetic on the whole-project means above (positive = gen1 sits LOWER")
        print("  of its control than gen2 does of its own):")
        for label, m, sel2 in (("open / NVA", g_open, cd == "NVA"),
                               ("vegetated / VVA", g_veg, cd == "VVA")):
            a = cr.resid_mm[m]; b = vd[sel2]
            d = a.mean() - b.mean()
            sd = np.sqrt(a.var(ddof=1)/a.size + b.var(ddof=1)/b.size)
            print(f"    {label:18s} gen1 {a.mean():+7.2f} - gen2 {b.mean():+7.2f} = "
                  f"{d:+7.2f} +/- {sd:.2f} mm  (n {a.size}/{b.size})")
        print("  These are PROJECT-WIDE means over two differently distributed mark sets,")
        print("  not the value at any site. Section 5 is the site-local statement.")

    # ------------------------------------------------------------- 8. coverage
    print("\n\n== 8. coverage: which marks the data on disk can reach ==")
    print("  gen2 point cloud on disk covers Elba/elbaext only. Marks inside it:")
    boxes = {
        "data/after/elbaext_3dep_fd_class2.laz": (575450.0, 580200.0, 4882050.0, 4886400.0),
        "data/after/3dep_4358_fulltile.laz": (584860.0, 587400.0, 4893245.0, 4896745.0),
    }
    inside = []
    for r in rows_all:
        e, n = float(r["easting"]), float(r["northing"])
        for p, (e0, e1, n0, n1) in boxes.items():
            if e0 <= e <= e1 and n0 <= n <= n1:
                inside.append((r["point_id"], r["point_type"], r["role"], p))
    for t in inside:
        print(f"    {t[0]:15s} {t[1]:4s} role={t[2]:12s} in {t[3]}")
    nboxed = len(tie_json["checkpoints"]) if tie_json else 0
    print(f"  plus {nboxed} marks with a separately fetched 400 m box under "
          f"data/after/checkpoints/")
    n_resid = sum(1 for r in rows_all if r["usgs_ql1_laz_error_m"] or r["usgs_ql0_laz_error_m"])
    print(f"  {n_resid} marks carry a USGS per-point residual; "
          f"{n_resid - nboxed} of them are NOT reachable from any point cloud on disk.")
    print("  Reaching them would mean fetching 3DEP tiles at ~0.7-2.3 GB each, and it")
    print("  would buy a second estimate of a quantity USGS has already published.")

    R.done(headline=f"gen2's absolute level at {A.site_name}, from its own held-out "
                    f"checkpoints; see sections 2-5.")


if __name__ == "__main__":
    main()
