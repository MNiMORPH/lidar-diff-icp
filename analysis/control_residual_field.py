#!/usr/bin/env python
"""gen1's vertical offset AT ELBA from the 2008 control residual modelled as a FIELD.

Every prior estimate reports ``sd/sqrt(n)`` of a sample mean over some set of marks.
This one fits the spatial field of the vendor's own published residual
(``dnr_error_m`` = ``Control Z - Surface Z``) and kriges it to a named site, so the
uncertainty attached to the answer is a prediction variance AT THAT LOCATION.

Nothing here has a default. Every radius, lag count, pair count, bin count, block size
and cover stratum is a required command-line argument, and the ones that are nuisance
choices are swept, with every point of the sweep printed.

Usage is in ``analysis/CONTROL_RESIDUAL_FIELD.md``.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from lidar_diff_icp.groundtruth import residual_field as RF   # noqa: E402
from trust.provenance import Run                              # noqa: E402

# The three cover treatments the brief asks for, named. Not a filter I invented: they
# are the three hypotheses about what the vegetated classes are doing.
VARIANTS = (
    ("open", ("L1O",), "raw"),
    ("open+urban", ("L1O", "L5U"), "raw"),
    ("cover-covariate", RF.COVER_CLASSES, "ols_residual"),
)

SERIES_COLOR = {"open": "#2a78d6", "open+urban": "#eb6834", "cover-covariate": "#1baf7a"}


def parse_md_table(path, header_must_contain):
    """Pull a fixed-width table out of a committed markdown analysis document.

    The header row is located by the substrings in ``header_must_contain``; the point
    id is whatever is left of the last ``len(header)-1`` whitespace-separated fields, so
    ids with spaces survive. Returns (header, list-of-rows).
    """
    lines = open(path).read().splitlines()
    hits = [k for k, ln in enumerate(lines)
            if all(s in ln for s in header_must_contain) and ln.strip().startswith("point")]
    if len(hits) != 1:
        raise SystemExit(f"{path}: expected exactly one header matching "
                         f"{header_must_contain}, found {len(hits)}")
    h = hits[0]
    hdr = lines[h].split()
    rows = []
    for ln in lines[h + 2:]:
        if not ln.strip() or ln.strip().startswith("```"):
            break
        parts = ln.strip().split()
        nfix = len(hdr) - 1
        rows.append([" ".join(parts[:len(parts) - nfix])] + parts[len(parts) - nfix:])
    return hdr, rows


def fmt(v, nd=1):
    return "--" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site-name", required=True)
    ap.add_argument("--site-easting", type=float, required=True)
    ap.add_argument("--site-northing", type=float, required=True)
    ap.add_argument("--sign-tol-m", type=float, required=True,
                    help="arithmetic tolerance for calling the sign identity exact")
    ap.add_argument("--band-radii-km", required=True,
                    help="comma list of radii for the distance-band sample-mean table")
    ap.add_argument("--max-lag-m", required=True, help="comma list, swept")
    ap.add_argument("--n-lags", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--n-pairs", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--estimators", required=True, help="comma list: dowd and/or matheron")
    ap.add_argument("--seeds", required=True, help="comma list, swept (nuisance)")
    ap.add_argument("--block-m", required=True, help="comma list of blocked-CV block sides")
    ap.add_argument("--loo-verify-idx", required=True,
                    help="comma list of mark indices to brute-force check the LOO shortcut on")
    ap.add_argument("--vgm-csv", required=True, help="where to write the full empirical variogram sweep")
    ap.add_argument("--fig", required=True, help="output PNG")
    ap.add_argument("--tie-table-18", required=True)
    ap.add_argument("--tie-table-16", required=True)
    a = ap.parse_args()

    radii = [float(s) for s in a.band_radii_km.split(",")]
    max_lags = [float(s) for s in a.max_lag_m.split(",")]
    n_lags_grid = [int(s) for s in a.n_lags.split(",")]
    n_pairs_grid = [int(s) for s in a.n_pairs.split(",")]
    estimators = a.estimators.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    blocks_m = [float(s) for s in a.block_m.split(",")]
    verify_idx = [int(s) for s in a.loo_verify_idx.split(",")]

    R = Run("what is the 2008 control residual AT the named site, as a value of a fitted "
            "spatial field with a prediction variance there, rather than as the mean of "
            "whatever marks happened to be nearby?")
    R.input(str(RF.DEFAULT_CONTROL_CSV),
            role="MnDNR 2008 SE-MN validation control: surveyed Control Z, the delivered "
                 "2008 surface Z, and the vendor's own residual dnr_error_m; no lidar "
                 "processing of ours enters it")
    R.input(a.tie_table_18,
            role="another agent's committed table, printed by analysis/groundtruth/"
                 "gen1_more_marks_report.py this session: 18 open/urban marks on lines "
                 "133-138 with OUR tie_mm (vendor class-2 ground, no lateral shift) beside "
                 "the vendor dnr_mm; READ, not recomputed here")
    R.input(a.tie_table_16,
            role="the earlier committed table from analysis/groundtruth/gen1_own_control_"
                 "tie.py: 16 marks inside the 13 gen1 tiles on disk with OUR tie_mm (CSF "
                 "ground, elbaext lateral shift) beside the vendor dnr_mm; READ, not "
                 "recomputed here")
    R.param("site", (a.site_name, a.site_easting, a.site_northing), src="andy",
            why="")
    R.param("sign_tol_m", a.sign_tol_m, src="MINE",
            why="tolerance for calling the Control-minus-Surface identity exact; it "
                "excludes nothing from any estimate, it only sets when the identity is "
                "declared to hold, and the max residual is printed beside it")
    R.param("band_radii_km", radii, src="MINE",
            why="the radii at which the sample-mean-and-SE table is re-derived; they "
                "select nothing downstream, the whole field fit uses every mark, and the "
                "point of the table is that the mean MOVES with this parameter")
    R.param("max_lag_m", max_lags, src="MINE",
            why="swept, not chosen: the variogram range is known to track max_lag when a "
                "field has long-wavelength structure, so every value is fitted and every "
                "resulting site prediction is printed")
    R.param("n_lags / n_pairs / seeds", (n_lags_grid, n_pairs_grid, seeds), src="MINE",
            why="nuisance grid of the binning and random-pair sampling; results are "
                "reported as median and full min-max across this grid so no single "
                "setting is privileged")
    R.param("estimators", estimators, src="andy", why="")
    R.param("block_m", blocks_m, src="MINE",
            why="swept block sides for spatially blocked cross-validation; each is run "
                "and printed, and no mark is dropped at any of them")
    R.param("loo_verify_idx", verify_idx, src="MINE",
            why="marks at which the leave-one-out shortcut is checked against a full "
                "refit; it verifies arithmetic and excludes nothing")

    cr = RF.load_residuals()
    site = (a.site_easting, a.site_northing)

    # ---------------------------------------------------------------- column defs
    R.column("what", "name of the row's quantity or stratum")
    R.column("value", "the number, units named in the row")
    R.column("radius_km", "radius about the site within which marks are counted, km")
    R.column("n_all", "marks of ANY cover class inside the radius")
    R.column("mean_all_mm", "mean of dnr_error_m over those marks, mm")
    R.column("se_all_mm", "sd/sqrt(n) of that mean, mm -- an SE OF A SAMPLE MEAN, not a prediction sd")
    R.column("n_ou", "marks of cover class L1O or L5U inside the radius")
    R.column("mean_ou_mm", "mean of dnr_error_m over the L1O+L5U marks, mm")
    R.column("se_ou_mm", "sd/sqrt(n) of that mean, mm")
    R.column("cover", "MnDNR land-cover class from the CSV point_type column")
    R.column("meaning", "what that class code denotes in the validation reports")
    R.column("n", "number of de-duplicated marks in the row")
    R.column("mean_mm", "mean of dnr_error_m over the row's marks, mm")
    R.column("median_mm", "median of dnr_error_m over the row's marks, mm")
    R.column("sd_mm", "sample sd of dnr_error_m over the row's marks, mm")
    R.column("se_mm", "sd/sqrt(n) of the row's mean, mm")
    R.column("variant", "cover treatment: which marks enter and whether cover is a drift term")
    R.column("estimator", "empirical-variogram estimator: dowd (robust) or matheron (classical)")
    R.column("max_lag_m", "largest separation entering the empirical variogram, m")
    R.column("nugget_mm2", "fitted spherical nugget, mm^2 (median over the nuisance grid)")
    R.column("sill_mm2", "fitted spherical PARTIAL sill (correlated variance), mm^2, median over the nuisance grid")
    R.column("range_m", "fitted spherical range (correlation length), m, median over the nuisance grid")
    R.column("sqrt_tot_mm", "sqrt(nugget+partial sill), mm -- the modelled total sd of the field")
    R.column("pred_mm", "kriged residual at the site, mm, median over the nuisance grid; positive = the delivered 2008 surface reads LOW")
    R.column("pred_lo_mm", "minimum kriged residual at the site across the nuisance grid, mm")
    R.column("pred_hi_mm", "maximum kriged residual at the site across the nuisance grid, mm")
    R.column("sd_field_mm", "kriging sd of the CORRELATED field component at the site (nugget filtered out), mm, median over the nuisance grid")
    R.column("sd_mark_mm", "kriging sd for the residual A NEW CONTROL MARK at the site would show (nugget included), mm, median over the nuisance grid")
    R.column("lag_m", "centre of the empirical-variogram lag bin, m")
    R.column("gamma_mm2", "empirical semivariance in that bin, mm^2")
    R.column("pairs", "random pairs falling in that bin")
    R.column("rmse_krige_mm", "root-mean-square cross-validation error of the kriged field, mm")
    R.column("rmse_null_mm", "root-mean-square cross-validation error of one global constant (fold mean), mm")
    R.column("mae_krige_mm", "mean absolute cross-validation error of the kriged field, mm")
    R.column("mae_null_mm", "mean absolute cross-validation error of one global constant (fold median), mm")
    R.column("skill", "1 - rmse_krige/rmse_null; positive = the field model beats the constant")
    R.column("n_L1O", "of the row's marks, how many are open terrain (L1O)")
    R.column("rmse_L1O_mm", "root-mean-square CV error of the kriged field, over the L1O marks ONLY, mm -- the common mark set on which the three variants are comparable")
    R.column("null_L1O_mm", "root-mean-square CV error of the fold constant, over the same L1O marks, mm")
    R.column("skill_L1O", "1 - rmse_L1O/null_L1O")
    R.column("block_m", "side of the square spatial block held out as one CV fold, m")
    R.column("n_blocks", "number of occupied blocks, i.e. folds")
    R.column("county", "MnDNR county whose validation report the mark came from")
    R.column("F", "one-way ANOVA F statistic across the grouping")
    R.column("p", "ANOVA p value")
    R.column("sd_group_means_mm", "sd of the group means, mm")
    R.column("sd_within_mm", "pooled within-group sd, mm")
    R.column("set", "which committed per-mark table the row summarises")
    R.column("dmean_mm", "mean of (our tie_mm - vendor dnr_mm) over the set, mm")
    R.column("dmedian_mm", "median of (our tie_mm - vendor dnr_mm), mm")
    R.column("dsd_mm", "sample sd of (our tie_mm - vendor dnr_mm), mm")
    R.column("dse_mm", "sd/sqrt(n) of dmean_mm, mm")
    R.column("r", "Pearson correlation of our tie_mm against the vendor dnr_mm over the set")
    R.column("sel_mean_mm", "mean vendor residual of the SELECTED marks, mm")
    R.column("field_at_sel_mm", "mean of the leave-one-out kriged field at those same locations, mm")
    R.column("excess_mm", "sel_mean_mm - field_at_sel_mm, mm: how far the selected marks sit above their own local field")
    R.column("off_mm", "surface-to-surface offset (our tie minus vendor residual), mm, from the named comparison set")
    R.column("carried18_mm", "pred_mm + the 18-mark offset: the residual OUR reconstructed gen1 surface is predicted to have at the site, mm; positive = our surface reads LOW, i.e. the constant to ADD")
    R.column("sd18_mm", "sqrt(sd_field_mm^2 + the 18-mark offset SE^2), mm")
    R.column("carried16_mm", "pred_mm + the 16-mark offset, mm, same sign convention")
    R.column("sd16_mm", "sqrt(sd_field_mm^2 + the 16-mark offset SE^2), mm")
    R.column("max_abs_err_mm", "largest absolute disagreement between the LOO shortcut and a full refit, mm")
    R.column("max_abs_var_mm2", "largest absolute disagreement in LOO kriging variance, mm^2")

    R.mask("deduped", np.ones(len(cr), bool), of=cr.n_rows_in,
           defn="rows of the control CSV kept after dropping exact repeats of "
                "(easting, northing, elevation) -- the same physical mark printed in two "
                "counties' validation reports")
    for name, covers, _ in VARIANTS:
        m = RF.stratify(cr, covers)
        R.mask(name, m, of=len(cr),
               defn=f"marks whose point_type is in {covers}")

    R.banner()

    out = []

    # ------------------------------------------------------------------ §1
    print("\n## 1  What the file is, re-derived\n")
    n_rows, n_cms, mx_cms, n_smc, mx_smc = RF.check_sign_convention(tol_m=a.sign_tol_m)
    R.table(["what", "value"], [
        ["rows in CSV", n_rows],
        ["unique marks after (E,N,Z) de-duplication", len(cr)],
        ["rows dropped as repeats", cr.n_dup_rows],
        ["groups those repeats came from", cr.n_dup_groups],
        ["rows where elevation - dnr_surface_z_m == dnr_error_m", f"{n_cms}/{n_rows}"],
        ["max |residual| of that identity, m", f"{mx_cms:.3e}"],
        ["rows where dnr_surface_z_m - elevation == dnr_error_m", f"{n_smc}/{n_rows}"],
        ["max |residual| of the reverse order, m", f"{mx_smc:.6f}"],
    ])

    # ------------------------------------------------------------------ §2
    print(f"\n## 2  The sample mean moves with the radius  (site = {a.site_name})\n")
    d_site = np.hypot(cr.easting - site[0], cr.northing - site[1])
    ou = RF.stratify(cr, ("L1O", "L5U"))
    rows = []
    for rk in radii:
        m = d_site <= rk * 1000.0
        v = cr.resid_mm[m]
        vo = cr.resid_mm[m & ou]
        rows.append([fmt(rk, 0), v.size, fmt(v.mean()), fmt(v.std(ddof=1) / np.sqrt(v.size)),
                     vo.size, fmt(vo.mean()), fmt(vo.std(ddof=1) / np.sqrt(vo.size))])
    R.table(["radius_km", "n_all", "mean_all_mm", "se_all_mm", "n_ou", "mean_ou_mm", "se_ou_mm"], rows)

    # ------------------------------------------------------------------ §3
    print("\n## 3  By cover class, over all 963 marks\n")
    rows = []
    for c in RF.COVER_CLASSES:
        m = cr.cover == c
        v = cr.resid_mm[m]
        rows.append([c, RF.COVER_MEANING[c], v.size, fmt(v.mean()), fmt(np.median(v)),
                     fmt(v.std(ddof=1)), fmt(v.std(ddof=1) / np.sqrt(v.size))])
    R.table(["cover", "meaning", "n", "mean_mm", "median_mm", "sd_mm", "se_mm"], rows)

    # ------------------------------------------------------------------ §4
    print("\n## 4  County as a DIAGNOSTIC, not a parameter\n")
    from scipy import stats as _st
    rows = []
    for name, covers, _ in VARIANTS:
        m = RF.stratify(cr, covers)
        groups = [cr.resid_mm[m & (cr.county == c)] for c in np.unique(cr.county[m])]
        groups = [g for g in groups if g.size > 1]
        F, p = _st.f_oneway(*groups)
        gm = np.array([g.mean() for g in groups])
        ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
        dfw = sum(g.size - 1 for g in groups)
        rows.append([name, len(groups), int(m.sum()), f"{F:.2f}", f"{p:.3e}",
                     fmt(gm.std(ddof=1)), fmt(np.sqrt(ssw / dfw))])
    R.column("n_counties", "counties with more than one mark in the stratum")
    R.table(["variant", "n_counties", "n", "F", "p", "sd_group_means_mm", "sd_within_mm"], rows)

    # ------------------------------------------------------------------ §5-6
    print("\n## 5  The variogram sweep and the site prediction\n")
    X_all, drift_labels, drift_at = RF.cover_design(cr.cover, RF.COVER_CLASSES)
    vgm_rows_csv = [("variant,estimator,max_lag_m,n_lags,n_pairs,seed,lag_m,gamma_mm2,pairs")]
    table_rows = []
    fitted = {}          # (variant, estimator, max_lag) -> dict of arrays for the figure
    for name, covers, vgm_on in VARIANTS:
        m = RF.stratify(cr, covers)
        x, y = cr.easting[m], cr.northing[m]
        if vgm_on == "raw":
            X, x0 = None, None
            vv = cr.resid_mm[m]
            vg_target = vv
        else:
            X = X_all[m]
            x0 = drift_at("L1O")
            vv = cr.resid_mm[m]
            beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
            vg_target = vv - X @ beta
        for est in estimators:
            for ml in max_lags:
                nug, sil, rng_, pred, sdf, sdm = [], [], [], [], [], []
                keep_emp = None
                for nl in n_lags_grid:
                    for npair in n_pairs_grid:
                        for sd in seeds:
                            mod, cen, gam, cnt = RF.fit_field(
                                x, y, vg_target, max_lag_m=ml, n_lags=nl,
                                n_pairs=npair, estimator=est, seed=sd)
                            kr = RF.krige(x, y, vv, mod, site[0], site[1],
                                          X=X, x0_drift=x0, drift_labels=drift_labels)
                            nug.append(mod.nugget); sil.append(mod.sill); rng_.append(mod.range_)
                            pred.append(kr.value_mm); sdf.append(kr.sd_field_mm)
                            sdm.append(kr.sd_new_mark_mm)
                            for L, G, C in zip(cen, gam, cnt):
                                vgm_rows_csv.append(
                                    f"{name},{est},{ml:.0f},{nl},{npair},{sd},"
                                    f"{L:.1f},{G:.2f},{C}")
                            if keep_emp is None:
                                keep_emp = (cen, gam, cnt, mod)
                fitted[(name, est, ml)] = keep_emp
                table_rows.append([
                    name, est, fmt(ml, 0), fmt(np.median(nug), 0), fmt(np.median(sil), 0),
                    fmt(np.median(rng_), 0), fmt(np.sqrt(np.median(nug) + np.median(sil))),
                    fmt(np.median(pred)), fmt(min(pred)), fmt(max(pred)),
                    fmt(np.median(sdf)), fmt(np.median(sdm))])
    R.table(["variant", "estimator", "max_lag_m", "nugget_mm2", "sill_mm2", "range_m",
             "sqrt_tot_mm", "pred_mm", "pred_lo_mm", "pred_hi_mm", "sd_field_mm",
             "sd_mark_mm"], table_rows)
    with open(a.vgm_csv, "w") as f:
        f.write("\n".join(vgm_rows_csv) + "\n")
    print(f"\n  full empirical variogram sweep ({len(vgm_rows_csv)-1} bins) written to {a.vgm_csv}")

    print("\n## 5b  The empirical variogram itself, at the widest and narrowest swept "
          f"max_lag, n_lags={n_lags_grid[-1]}, n_pairs={n_pairs_grid[-1]}, seed={seeds[0]}, "
          f"estimator={estimators[0]}\n")
    emp_rows = []
    for name, _, _ in VARIANTS:
        for ml in (min(max_lags), max(max_lags)):
            cen, gam, cnt, mod = fitted[(name, estimators[0], ml)]
            for L, G, C in zip(cen, gam, cnt):
                emp_rows.append([name, fmt(ml, 0), fmt(L, 0), fmt(G, 0), int(C)])
    R.table(["variant", "max_lag_m", "lag_m", "gamma_mm2", "pairs"], emp_rows)

    print("\n## 5c  Under the cover-covariate variant, the site prediction for EACH cover "
          f"class (estimator={estimators[0]}, nuisance grid medians)\n")
    m = RF.stratify(cr, RF.COVER_CLASSES)
    x, y, vv = cr.easting[m], cr.northing[m], cr.resid_mm[m]
    X = X_all[m]
    beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
    tgt = vv - X @ beta
    cls_rows = []
    for ml in (min(max_lags), max(max_lags)):
        for c in RF.COVER_CLASSES:
            pv, sf = [], []
            for nl in n_lags_grid:
                for npair in n_pairs_grid:
                    for sd in seeds:
                        mod, *_ = RF.fit_field(x, y, tgt, max_lag_m=ml, n_lags=nl,
                                               n_pairs=npair, estimator=estimators[0], seed=sd)
                        kr = RF.krige(x, y, vv, mod, site[0], site[1], X=X,
                                      x0_drift=drift_at(c), drift_labels=drift_labels)
                        pv.append(kr.value_mm); sf.append(kr.sd_field_mm)
            cls_rows.append([c, RF.COVER_MEANING[c], fmt(ml, 0), fmt(np.median(pv)),
                             fmt(min(pv)), fmt(max(pv)), fmt(np.median(sf))])
    R.table(["cover", "meaning", "max_lag_m", "pred_mm", "pred_lo_mm", "pred_hi_mm",
             "sd_field_mm"], cls_rows)

    # ------------------------------------------------------------------ §6 LOO check
    print("\n## 6  The leave-one-out shortcut, checked against a full refit\n")
    vrows = []
    for name, covers, vgm_on in VARIANTS:
        m = RF.stratify(cr, covers)
        x, y, vv = cr.easting[m], cr.northing[m], cr.resid_mm[m]
        X = None if vgm_on == "raw" else X_all[m]
        tgt = vv
        if X is not None:
            beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
            tgt = vv - X @ beta
        mod, *_ = RF.fit_field(x, y, tgt, max_lag_m=max_lags[0], n_lags=n_lags_grid[0],
                               n_pairs=n_pairs_grid[0], estimator=estimators[0], seed=seeds[0])
        idx = [i for i in verify_idx if i < x.size]
        de, dv = RF.verify_loo_shortcut(x, y, vv, mod, idx, X=X)
        vrows.append([name, len(idx), f"{de:.2e}", f"{dv:.2e}"])
    R.table(["variant", "n", "max_abs_err_mm", "max_abs_var_mm2"], vrows)

    # ------------------------------------------------------------------ §7 LOO CV
    print("\n## 7  Leave-one-out cross-validation against the null of one global constant\n")
    cv_rows = []
    for name, covers, vgm_on in VARIANTS:
        m = RF.stratify(cr, covers)
        x, y, vv = cr.easting[m], cr.northing[m], cr.resid_mm[m]
        X = None if vgm_on == "raw" else X_all[m]
        tgt = vv
        if X is not None:
            beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
            tgt = vv - X @ beta
        em, ed = RF.constant_null_errors(vv, np.arange(vv.size))
        op = cr.cover[m] == "L1O"
        for est in estimators:
            for ml in max_lags:
                rm, ma, ro = [], [], []
                for nl in n_lags_grid:
                    for npair in n_pairs_grid:
                        for sd in seeds:
                            mod, *_ = RF.fit_field(x, y, tgt, max_lag_m=ml, n_lags=nl,
                                                   n_pairs=npair, estimator=est, seed=sd)
                            err, _ = RF.loo_errors(x, y, vv, mod, X=X)
                            rm.append(np.sqrt(np.mean(err ** 2)))
                            ma.append(np.mean(np.abs(err)))
                            ro.append(np.sqrt(np.mean(err[op] ** 2)))
                rn = np.sqrt(np.mean(em ** 2))
                rno = np.sqrt(np.mean(em[op] ** 2))
                cv_rows.append([name, est, fmt(ml, 0), vv.size,
                                fmt(np.median(rm)), fmt(rn),
                                fmt(np.median(ma)), fmt(np.mean(np.abs(ed))),
                                fmt(1 - np.median(rm) / rn, 3),
                                int(op.sum()), fmt(np.median(ro)), fmt(rno),
                                fmt(1 - np.median(ro) / rno, 3)])
    R.table(["variant", "estimator", "max_lag_m", "n", "rmse_krige_mm", "rmse_null_mm",
             "mae_krige_mm", "mae_null_mm", "skill",
             "n_L1O", "rmse_L1O_mm", "null_L1O_mm", "skill_L1O"], cv_rows)
    print("\n  The last four columns score all three variants on THE SAME marks (every")
    print("  L1O mark each variant contains), because rmse over different mark sets is not")
    print("  a comparison. The open variant's restricted and unrestricted columns are the")
    print("  same numbers by construction.\n")

    # ------------------------------------------------------------------ §8 block CV
    print("\n## 8  Spatially blocked cross-validation, variogram REFITTED inside every "
          f"training fold; nuisance settings held at n_lags={n_lags_grid[0]}, "
          f"n_pairs={n_pairs_grid[0]}, seed={seeds[0]}, estimator={estimators[0]}\n")
    bcv_rows = []
    for name, covers, vgm_on in VARIANTS:
        m = RF.stratify(cr, covers)
        x, y, vv = cr.easting[m], cr.northing[m], cr.resid_mm[m]
        X = None if vgm_on == "raw" else X_all[m]
        op = cr.cover[m] == "L1O"
        for bm in blocks_m:
            for ml in max_lags:
                err, bid, nb = RF.block_cv(
                    x, y, vv, block_m=bm, max_lag_m=ml, n_lags=n_lags_grid[0],
                    n_pairs=n_pairs_grid[0], estimator=estimators[0], seed=seeds[0],
                    X=X, refit_variogram=True, variogram_on=vgm_on)
                em, ed = RF.constant_null_errors(vv, bid)
                rk = np.sqrt(np.nanmean(err ** 2))
                rn = np.sqrt(np.nanmean(em ** 2))
                rko = np.sqrt(np.nanmean(err[op] ** 2))
                rno = np.sqrt(np.nanmean(em[op] ** 2))
                bcv_rows.append([name, fmt(bm, 0), nb, fmt(ml, 0), vv.size,
                                 fmt(rk), fmt(rn), fmt(np.nanmean(np.abs(err))),
                                 fmt(np.nanmean(np.abs(ed))), fmt(1 - rk / rn, 3),
                                 int(op.sum()), fmt(rko), fmt(rno),
                                 fmt(1 - rko / rno, 3)])
    R.table(["variant", "block_m", "n_blocks", "max_lag_m", "n", "rmse_krige_mm",
             "rmse_null_mm", "mae_krige_mm", "mae_null_mm", "skill",
             "n_L1O", "rmse_L1O_mm", "null_L1O_mm", "skill_L1O"], bcv_rows)

    # ------------------------------------------------------------------ §9 our surface
    print("\n## 9  Our reconstructed surface against the DELIVERED one, at the marks "
          "where both have been read\n")
    id_to_i = {p: i for i, p in enumerate(cr.point_id)}
    sets = {}
    for label, path, need in (
        ("18mk-lines133-138-class2-nolat", a.tie_table_18, ("cls", "dnr_mm", "tie_mm")),
        ("16mk-tiles-on-disk-CSF-lat", a.tie_table_16, ("cover", "dnr_mm", "tie_mm")),
    ):
        hdr, rows = parse_md_table(path, need)
        ti, di = hdr.index("tie_mm"), hdr.index("dnr_mm")
        ids = [r[0] for r in rows]
        miss = [p for p in ids if p not in id_to_i]
        if miss:
            raise SystemExit(f"{path}: point_id(s) not in the control CSV: {miss}")
        tie = np.array([float(r[ti]) for r in rows])
        dnr = np.array([float(r[di]) for r in rows])
        sets[label] = (ids, tie, dnr)

    off_rows = []
    for label, (ids, tie, dnr) in sets.items():
        d = tie - dnr
        off_rows.append([label, len(ids), fmt(d.mean()), fmt(np.median(d)), fmt(d.std(ddof=1)),
                         fmt(d.std(ddof=1) / np.sqrt(d.size)),
                         fmt(np.corrcoef(tie, dnr)[0, 1], 3)])
    R.table(["set", "n", "dmean_mm", "dmedian_mm", "dsd_mm", "dse_mm", "r"], off_rows)
    print("\n  (our tie - vendor residual) = (surveyed - z_ours) - (surveyed - z_vendor)")
    print("  = z_vendor - z_ours.  POSITIVE means our reconstructed surface sits BELOW the")
    print("  delivered 2008 surface, so this is the constant to ADD to a field prediction")
    print("  made on the vendor residual before it applies to our surface.\n")

    # ------------------------------------------------------------------ §10 selection bias
    print("\n## 10  Where the selected marks sit relative to their own local field\n")
    sb_rows = []
    for name, covers, vgm_on in VARIANTS:
        m = RF.stratify(cr, covers)
        x, y, vv = cr.easting[m], cr.northing[m], cr.resid_mm[m]
        X = None if vgm_on == "raw" else X_all[m]
        tgt = vv
        if X is not None:
            beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
            tgt = vv - X @ beta
        sub_idx = {p: k for k, p in enumerate(cr.point_id[m])}
        mod, *_ = RF.fit_field(x, y, tgt, max_lag_m=max_lags[0], n_lags=n_lags_grid[-1],
                               n_pairs=n_pairs_grid[-1], estimator=estimators[0], seed=seeds[0])
        err, _ = RF.loo_errors(x, y, vv, mod, X=X)
        loo_pred = err + vv
        for label, (ids, tie, dnr) in sets.items():
            k = [sub_idx[p] for p in ids if p in sub_idx]
            if not k:
                continue
            sb_rows.append([name, label, len(k), fmt(vv[k].mean()),
                            fmt(loo_pred[k].mean()),
                            fmt(vv[k].mean() - loo_pred[k].mean())])
    R.column("set2", "which committed per-mark table the selected marks come from")
    R.table(["variant", "set2", "n", "sel_mean_mm", "field_at_sel_mm", "excess_mm"],
            [[r[0], r[1], r[2], r[3], r[4], r[5]] for r in sb_rows])
    print("\n  The field prediction at each selected mark is LEAVE-ONE-OUT, so the mark")
    print("  itself does not enter the field it is compared against.\n")

    # population within each swept radius, for the same comparison the brief asks for
    print("\n## 10b  The selected marks against the plain local population mean\n")
    pop_rows = []
    for label, (ids, tie, dnr) in sets.items():
        k = [id_to_i[p] for p in ids]
        for rk in radii:
            mm = d_site <= rk * 1000.0
            pop_rows.append([label, len(k), fmt(rk, 0), fmt(cr.resid_mm[k].mean()),
                             fmt(cr.resid_mm[mm].mean()),
                             fmt(cr.resid_mm[k].mean() - cr.resid_mm[mm].mean())])
    R.column("pop_mean_mm", "mean vendor residual of ALL marks within radius_km of the site, mm")
    R.table(["set", "n", "radius_km", "sel_mean_mm", "pop_mean_mm", "excess_mm"], pop_rows)

    # ------------------------------------------------------------------ §11 carry
    print("\n## 11  The field prediction CARRIED to our own reconstructed surface\n")
    o18 = sets["18mk-lines133-138-class2-nolat"]
    o16 = sets["16mk-tiles-on-disk-CSF-lat"]
    d18 = o18[1] - o18[2]
    d16 = o16[1] - o16[2]
    m18, s18 = d18.mean(), d18.std(ddof=1) / np.sqrt(d18.size)
    m16, s16 = d16.mean(), d16.std(ddof=1) / np.sqrt(d16.size)
    carry_rows = []
    for r in table_rows:
        pv, sf = float(r[7]), float(r[10])
        carry_rows.append([r[0], r[1], r[2], fmt(pv), fmt(sf),
                           fmt(pv + m18), fmt(np.hypot(sf, s18)),
                           fmt(pv + m16), fmt(np.hypot(sf, s16))])
    R.table(["variant", "estimator", "max_lag_m", "pred_mm", "sd_field_mm",
             "carried18_mm", "sd18_mm", "carried16_mm", "sd16_mm"], carry_rows)
    print(f"\n  offsets carried: 18-mark set {m18:+.1f} +/- {s18:.1f} mm, "
          f"16-mark set {m16:+.1f} +/- {s16:.1f} mm (mean +/- sd/sqrt(n) of")
    print("  our tie minus the vendor residual at the marks where both have been read).")
    print("  The offset own uncertainty is an SE OF A MEAN over marks; only the field")
    print("  term is a prediction sd at the site, and the two are added in quadrature.\n")

    # ------------------------------------------------------------------ figure
    make_figure(a, cr, site, fitted, table_rows, max_lags, estimators, sets, id_to_i)
    print(f"\n  figure written to {a.fig}")

    R.done(headline=f"control residual field kriged to {a.site_name}")


def make_figure(a, cr, site, fitted, table_rows, max_lags, estimators, sets, id_to_i):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    SURF = "#fcfcfb"
    INK = "#0b0b0b"
    INK2 = "#52514e"
    GRID = "#e3e2de"

    fig = plt.figure(figsize=(13.0, 8.6), facecolor=SURF)
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.26)

    est = estimators[0]
    for k, (name, _, _) in enumerate(VARIANTS):
        ax = fig.add_subplot(gs[0, k])
        ax.set_facecolor(SURF)
        col = SERIES_COLOR[name]
        for j, ml in enumerate(max_lags):
            cen, gam, cnt, mod = fitted[(name, est, ml)]
            alpha = 0.30 + 0.65 * j / max(len(max_lags) - 1, 1)
            ax.plot(cen / 1000.0, gam / 1000.0, "o", ms=3.5, color=col, alpha=alpha,
                    mec="none", zorder=2)
            hh = np.linspace(0, ml, 300)
            hr = np.clip(hh / mod.range_, 0, 1)
            ax.plot(hh / 1000.0, (mod.nugget + mod.sill * (1.5 * hr - 0.5 * hr ** 3)) / 1000.0,
                    "-", lw=2.0, color=col, alpha=alpha, zorder=3)
        ax.set_title(name, color=INK, fontsize=11, loc="left", pad=6)
        ax.set_xlabel("lag  (km)", color=INK2, fontsize=9)
        if k == 0:
            ax.set_ylabel(r"semivariance  ($10^3$ mm$^2$)", color=INK2, fontsize=9)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8)

    # prediction vs max_lag
    ax = fig.add_subplot(gs[1, 0])
    ax.set_facecolor(SURF)
    for name, _, _ in VARIANTS:
        col = SERIES_COLOR[name]
        rows = [r for r in table_rows if r[0] == name and r[1] == est]
        ml = np.array([float(r[2]) for r in rows]) / 1000.0
        pv = np.array([float(r[7]) for r in rows])
        sdf = np.array([float(r[10]) for r in rows])
        ax.fill_between(ml, pv - sdf, pv + sdf, color=col, alpha=0.13, lw=0)
        ax.plot(ml, pv, "-o", color=col, lw=2.0, ms=5, label=name)
        ax.annotate(name, (ml[-1], pv[-1]), textcoords="offset points", xytext=(4, 3),
                    color=col, fontsize=8.5)
    ax.axhline(0, color=INK2, lw=0.8, ls="--")
    ax.set_xlabel("variogram max lag  (km)", color=INK2, fontsize=9)
    ax.set_ylabel(f"kriged residual at {a.site_name}  (mm)", color=INK2, fontsize=9)
    ax.set_title("prediction, +/- the field kriging sd", color=INK, fontsize=11, loc="left", pad=6)
    ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.legend(frameon=False, fontsize=8, loc="best", labelcolor=INK2)

    # map
    ax = fig.add_subplot(gs[1, 1:])
    ax.set_facecolor(SURF)
    lim = float(np.percentile(np.abs(cr.resid_mm), 98))
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    # Two hues either side of a neutral midpoint -- blue = surface reads HIGH,
    # orange = surface reads LOW. Never a rainbow, never a hue at the midpoint.
    cmap = LinearSegmentedColormap.from_list(
        "resid", ["#184f95", "#2a78d6", "#f0efec", "#eb6834", "#8c3410"])
    sc = ax.scatter(cr.easting / 1000.0, cr.northing / 1000.0, c=cr.resid_mm, s=13,
                    cmap=cmap, norm=norm, lw=0.2, edgecolor="#ffffff", zorder=2)
    ax.scatter([site[0] / 1000.0], [site[1] / 1000.0], marker="*", s=280, color="#0b0b0b",
               zorder=5)
    ax.annotate(a.site_name, (site[0] / 1000.0, site[1] / 1000.0),
                textcoords="offset points", xytext=(9, -12), color=INK, fontsize=10)
    for label, (ids, tie, dnr) in sets.items():
        k = [id_to_i[p] for p in ids]
        ax.scatter(cr.easting[k] / 1000.0, cr.northing[k] / 1000.0, s=64,
                   facecolor="none", edgecolor="#0b0b0b", lw=1.1, zorder=4)
    ax.set_aspect("equal")
    ax.set_xlabel("EPSG:26915 easting  (km)", color=INK2, fontsize=9)
    ax.set_ylabel("northing  (km)", color=INK2, fontsize=9)
    ax.set_title("vendor residual Control Z - Surface Z; ringed = marks our own tie has "
                 "also been read at", color=INK, fontsize=10.5, loc="left", pad=6)
    ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    cb = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.035)
    cb.set_label("mm  (positive = delivered surface reads LOW)", color=INK2, fontsize=8.5)
    cb.ax.tick_params(colors=INK2, labelsize=8)
    cb.outline.set_edgecolor(GRID)

    fig.savefig(a.fig, dpi=170, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
