#!/usr/bin/env python3
"""Head to head: our along-track drift spline vs DeLong's IDW correction surface.

Andy, 2026-09-11: "What is better: the DeLong-style corrected original product or our
custom product?"

THEY ARE NOT RIVALS FOR THE SAME ERROR, and half the answer needs no run at all.
`pipeline.difference_dem` applies the correction surface as ``C[iyp, ixp]`` -- a function of
GRID POSITION ONLY. Two flight lines sharing a cell therefore receive the IDENTICAL
correction, so it cancels EXACTLY in a paired same-cell difference: **a spatial correction
surface is structurally incapable of changing a per-line step.** `fit_along_track_drift` is
per-swath in `gps_time`, so it can. The overlap test found ~17.3 mm of the 65.2 mm apparent
line structure to be genuinely per-line; only the drift can touch that part. What remains to
be measured is the other ~73%: the spatial field, which both methods CAN address.

WHY OUT-OF-SAMPLE IS MANDATORY. Both are FITTED to the stable-ground residual, so both
improve it in-sample by construction, and the more flexible one wins every time. A flexible
IDW on the target has burned us before: an earlier convex-prior IDW correction looked clean,
dropped false convex deposition 12% -> 1%, and was OVERFIT -- its out-of-sample skill was
0.14 and the extra ~35 mm it removed was terrain-correlated noise, not warp. In-sample
numbers here would repeat that mistake.

THE BLOCK SIZE DECIDES THE ANSWER, SO IT IS SWEPT, NOT CHOSEN. Held-out blocks must be
larger than the correlation length or an IDW simply interpolates across the gap and scores
its own neighbours. Small blocks flatter the IDW; large blocks starve it. Picking one would
BE the result, so every size is reported and the trend is the finding. The drift spline is
refit on the same training points each fold, so both methods pay the same price.

Base for every route: d_mm + geoid + lateral + swath (the terms both methods share).
  ours    = base + along-track drift spline, refit per fold
  delong  = base + IDW correction surface, refit per fold
Reported on HELD-OUT cells only.

    ./lidar-icp/bin/python ground_control/drift_vs_correction_surface.py --tile elba_fulldensity
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from lidar_diff_icp import coreg
from trust.provenance import Run

COLS = ["d_mm", "dz_geoid_mm", "dz_lateral_mm", "dz_swath_mm", "dz_drift_mm",
        "gps_time", "point_source_id", "cell", "slope", "in_grid"]


def nmad(v):
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", required=True)
    ap.add_argument("--block-m", type=float, nargs="+",
                    default=[100.0, 200.0, 400.0, 800.0, 1600.0],
                    help="held-out block sizes to SWEEP. Not a choice: the block size "
                         "decides which method wins, so the trend is the result")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260911)
    ap.add_argument("--slope-max-deg", type=float, default=3.0)
    ap.add_argument("--dz-max-mm", type=float, default=700.0)
    ap.add_argument("--radius-m", type=float, default=400.0)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    d = f"data/derived/{a.tile}"
    tab, fpp = os.path.join(d, "beam_offset_table.parquet"), os.path.join(d, "floodplain_mask.npy")
    corr = json.load(open(os.path.join(d, "corrections.json")))
    res = float(corr["res_m"]); X0, Y0 = corr["bounds"][0], corr["bounds"][1]

    R = Run("Out of sample, which removes more of the SPATIAL stable-ground residual: our "
            "per-swath along-track drift spline, or DeLong's IDW correction surface?")
    R.input(tab, role="per-return offsets with every correction term separable, so any "
                      "combination of routes can be rebuilt from the SAME returns")
    R.input(fpp, role="elevation-cut floodplain mask; stands in for DeLong's 100 m stream buffer")
    R.input(os.path.join(d, "corrections.json"), role="grid geometry (res, origin) for the "
                                                      "correction surface")
    R.param("block_m", tuple(a.block_m), src="MINE",
            why="SWEPT, never chosen: held-out blocks smaller than the correlation length "
                "let an IDW score its own neighbours, larger ones starve it, so a single "
                "value would BE the answer. Excludes no data at any size -- every cell is "
                "held out exactly once per size")
    R.param("folds", a.folds, src="MINE",
            why="k-fold over whole blocks. Every cell is predicted exactly once from a model "
                "that never saw its block; nothing is discarded")
    R.param("seed", a.seed, src="MINE", why="fixed so the fold assignment reproduces")
    R.param("radius_m", a.radius_m, src="repo",
            why="coreg.correction_surface's own default and DeLong's 400 m; the honest reach "
                "of the interpolation, beyond which it injects error")
    R.param("slope_max_deg", a.slope_max_deg, src="MINE",
            why="DeLong s4.4 'areas with slope >3 deg', verified verbatim in the extracted "
                "text. MINE because a published constant has no category here and adopting "
                "it was my choice; the value is checkable in the paper")
    R.param("dz_max_mm", a.dz_max_mm, src="MINE",
            why="DeLong's companion cut, 'DoD values >0.7 m in flat areas'. Applied to the "
                "RAW d_mm so the identical returns serve every route")
    R.param("tile", a.tile, src="MINE", why="the tile with a built table; not chosen by its answer")
    R.column("route", "base = geoid+lateral+swath (shared); ours = base+drift spline; "
                      "delong = base+IDW correction surface")
    R.column("oos_sd", "sd of the held-out stable cell residual, ddof=1, mm")
    R.column("oos_nmad", "1.4826*MAD of the same, mm -- robust twin of oos_sd")
    R.column("extrap%", "percentage of held-out points lying OUTSIDE their own swath's "
                        "training gps_time span. A spline cannot extrapolate, so this "
                        "attributes a large-block blow-up rather than leaving it mysterious")
    R.column("skill", "1 - var(route)/var(base) on held-out cells. Positive = better than "
                      "doing nothing. NEGATIVE MEANS IT INJECTS ERROR")

    t = pd.read_parquet(tab, columns=COLS)
    n_all = len(t)
    fpm = np.load(fpp); ny, nx = fpm.shape
    iy, ix = np.divmod(t["cell"].to_numpy(), nx)
    ok = (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx)
    infp = np.zeros(n_all, bool); infp[ok] = fpm[iy[ok], ix[ok]].astype(bool)
    stable = (t["in_grid"].to_numpy().astype(bool) & ok & ~infp
              & (t["slope"].to_numpy() <= a.slope_max_deg)
              & (np.abs(t["d_mm"].to_numpy()) <= a.dz_max_mm))
    R.mask("stable", stable,
           defn=f"in grid, slope <= {a.slope_max_deg} deg, |raw dz| <= {a.dz_max_mm} mm, not "
                f"floodplain -- DeLong's definition with their stream buffer substituted",
           of=n_all)

    s = t[stable].copy()
    del t
    base_pt = (s["d_mm"] + s["dz_geoid_mm"] + s["dz_lateral_mm"] + s["dz_swath_mm"]).to_numpy()
    gt = s["gps_time"].to_numpy(); ps = s["point_source_id"].to_numpy()
    cell = s["cell"].to_numpy(); ciy, cix = np.divmod(cell, nx)
    print(f"  stable returns {len(s):,} of {n_all:,} ({100*len(s)/n_all:.1f}%), "
          f"grid {ny}x{nx} at {res:g} m")

    # ---- IN-SAMPLE SANITY CHECK, not a result. Both methods are fitted, so both MUST
    # improve in sample. A negative in-sample skill means the harness has a sign or
    # convention error, not that the method is bad. This caught a flipped dz once.
    allm = np.ones(len(s), bool)
    dr_in, _ = coreg.fit_along_track_drift(gt, -base_pt / 1000.0, allm, ps)
    Zb0 = np.full((ny, nx), np.nan)
    cm0 = pd.Series(base_pt).groupby(cell).mean()
    Zb0[np.divmod(cm0.index.to_numpy(), nx)] = cm0.to_numpy()
    C0 = coreg.correction_surface(np.zeros_like(Zb0), Zb0 / 1000.0, res, X0, Y0,
                                  radius=a.radius_m, slope_thresh_deg=90.0,
                                  dz_thresh=1e9)["C"]
    c0 = C0[ciy, cix] * 1000.0
    ins = {"base": base_pt, "ours": base_pt + dr_in * 1000.0,
           "delong": base_pt + np.where(np.isfinite(c0), c0, 0.0)}
    vb0 = np.nanvar(pd.Series(ins["base"]).groupby(cell).mean().to_numpy())
    print("\n  IN-SAMPLE sanity (must be positive for both, else the harness is wrong):")
    for lab in ("ours", "delong"):
        vv = pd.Series(ins[lab]).groupby(cell).mean().to_numpy()
        print(f"    {lab:>8} in-sample skill {1.0 - np.nanvar(vv)/vb0:+.3f}")

    rows = []
    print(f"\n{'block_m':>8}{'route':>9}{'oos_sd':>10}{'oos_nmad':>10}{'skill':>9}"
          f"{'extrap%':>9}{'cells':>10}")
    for L in a.block_m:
        per = max(1, int(round(L / res)))
        blk = (ciy // per).astype(np.int64) * 100000 + (cix // per)
        ub = np.unique(blk)
        rng = np.random.default_rng(a.seed)
        fold_of = dict(zip(ub, rng.integers(0, a.folds, ub.size)))
        fold = np.array([fold_of[b] for b in blk])

        n_extrap, n_te = [0], [0]
        pred_ours = np.full(len(s), np.nan)
        pred_del = np.full(len(s), np.nan)
        for k in range(a.folds):
            tr, te = fold != k, fold == k
            if te.sum() == 0 or tr.sum() < 1000:
                continue
            # --- ours: refit the per-swath drift spline on TRAINING points only
            drift, _ = coreg.fit_along_track_drift(gt, -base_pt / 1000.0, tr, ps)
            pred_ours[te] = base_pt[te] + drift[te] * 1000.0
            # a spline cannot extrapolate: count held-out points outside their own
            # swath's TRAINING gps_time span, so a blow-up is attributable not mysterious
            for sw in np.unique(ps[te]):
                mtr, mte = tr & (ps == sw), te & (ps == sw)
                if mtr.sum() == 0:
                    n_extrap[0] += mte.sum(); continue
                lo, hi = gt[mtr].min(), gt[mtr].max()
                n_extrap[0] += int(((gt[mte] < lo) | (gt[mte] > hi)).sum())
            n_te[0] += int(te.sum())
            # --- delong: IDW correction surface from TRAINING cells only
            Zb = np.full((ny, nx), np.nan)
            cm = pd.Series(base_pt[tr]).groupby(blk[tr] * 0 + cell[tr]).mean()
            Zb[np.divmod(cm.index.to_numpy(), nx)] = cm.to_numpy()
            # SIGN: correction_surface forms dz = z_ref - z_src and its C is ADDED to
            # z_src (pipeline.py does `zc += Cpt`). base_pt carries gen1's sense, so
            # z_src = base and z_ref = 0 gives dz = -base, the true gen2-gen1 residual.
            # Passing -base instead flips C and DOUBLES the error -- the in-sample column
            # below exists to catch exactly that, and did.
            C = coreg.correction_surface(np.zeros_like(Zb), Zb / 1000.0, res, X0, Y0,
                                         radius=a.radius_m,
                                         slope_thresh_deg=90.0, dz_thresh=1e9)["C"]
            cpt = C[ciy[te], cix[te]] * 1000.0
            pred_del[te] = base_pt[te] + np.where(np.isfinite(cpt), cpt, 0.0)

        # aggregate to CELLS before scoring: a per-return score would count the same
        # ground thousands of times and is not an independent sample
        cb = pd.Series(base_pt).groupby(cell).mean().to_numpy()
        co = pd.Series(pred_ours).groupby(cell).mean().to_numpy()
        cd = pd.Series(pred_del).groupby(cell).mean().to_numpy()
        good = np.isfinite(cb) & np.isfinite(co) & np.isfinite(cd)
        vb = np.nanvar(cb[good])
        for lab, v in (("base", cb[good]), ("ours", co[good]), ("delong", cd[good])):
            sk = float("nan") if lab == "base" else 1.0 - np.nanvar(v) / vb
            ex = "" if lab != "ours" else f"{100.0*n_extrap[0]/max(n_te[0],1):>9.1f}"
            print(f"{L:>8.0f}{lab:>9}{np.nanstd(v, ddof=1):>10.2f}{nmad(v):>10.2f}"
                  f"{'' if lab=='base' else f'{sk:>9.3f}'}{ex if lab=='ours' else '':>9}"
                  f"{good.sum():>10,}")
            rows.append(dict(block_m=L, route=lab, oos_sd=float(np.nanstd(v, ddof=1)),
                             oos_nmad=nmad(v), skill=None if lab == "base" else float(sk),
                             extrap_frac=(n_extrap[0]/max(n_te[0],1)) if lab == "ours" else None,
                             n_cells=int(good.sum())))
        print()

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(dict(tile=a.tile, rows=rows), fh, indent=1); fh.write("\n")
        print(f"    wrote {a.out}")
    best = {r["block_m"]: r for r in rows if r["route"] == "ours"}
    R.done(headline="; ".join(
        f"{int(L)} m: ours skill {best[L]['skill']:+.3f} vs delong "
        f"{[r for r in rows if r['block_m']==L and r['route']=='delong'][0]['skill']:+.3f}"
        for L in sorted(best)))


if __name__ == "__main__":
    main()
