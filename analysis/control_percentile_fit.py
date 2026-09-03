#!/usr/bin/env python3
"""WHERE the surveyed ground sits inside gen2's own return column, as a function of lowveg.

The control analogue of the ridgeline q2 fit. Every percentile-picking decision on this
project -- ``q2 = 0.5 - 0.19*cover`` and its descendants -- has inferred the right quantile
from gen1-vs-gen2 agreement on divide cells, with no external truth. Here the answer is
measured directly: the histogram supplies the SHAPE, the surveyed elevation supplies the
POSITION, and true ground lands at a definite rank inside the observed distribution.

    q_ctrl = fraction of the mark's near-ground returns lying BELOW the surveyed ground

Nothing is re-specified. The mark population, the lowveg metric and the band edges come
from control_lowveg_offset.load()/lowveg(), so this fit and the -290 mm/unit fit rest on
exactly the same marks. Bins are the same uniform 0.06 in lowveg.

THE ONE COORDINATE DECISION, AND WHY IT IS NOT A CHOICE. The histogram's zero is OUR
order-2 least-squares surface through the box's class-2 returns, and heights are
``(z - S)/|n|`` with ``|n| = sqrt(1 + gx^2 + gy^2)`` (cover_at_control_marks._sn_hist). To
locate the surveyed ground inside that histogram it must be expressed in the SAME frame:

    h_survey = (surveyed_Z - surface_coef[0]) / |n|          [m, slope-normal]

The USGS published residual cannot be used for this. It is referenced to the vendor's own
DEM, not to our surface, and the two diverge in vegetation -- measured at these very marks
as -138 mm (ours) against -76 mm (published). Mixing them would put the surveyed height in
one frame and the returns in another. The consequence is stated rather than hidden: q_ctrl
inherits whatever bias our least-squares surface carries under canopy.

    ./lidar-icp/bin/python analysis/control_percentile_fit.py
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_lowveg_offset import STRUCT, load          # same marks, same lowveg, same edges

BAND_LO, BAND_HI = 0.15, 2.00       # CONTROL_LOWVEG_OFFSET.md's exact definition
BIN_W = 0.06                        # the mm fit's uniform bin width
SET = "gen2_2021_control"

_ap = argparse.ArgumentParser()
# The stored boxes carry TWO windows. The near-ground one stops at +2.00 m, which is a
# convention, not a physical edge: it was chosen to exclude tree crowns. Whether that choice
# moves the answer is testable here rather than arguable, because the tall window covers
# -2..+45 m and so has no upper truncation at all. Its bins are 0.25 m against 0.02 m, so it
# trades vertical resolution for completeness -- both numbers are reported, neither is
# DEFAULT = full, and NO clipping of any kind (Andy, 2026-09-03). The near window's +2.00 m
# top excludes tree crowns, which is a decision about what counts as "the return column"
# taken before knowing whether it matters. Measured, it barely does: on the same marks the
# slope moves -0.849 -> -0.922 and the intercept +0.4385 -> +0.4561, one to two SE. Since it
# is arbitrary and nearly free, the untruncated column is the honest default. It also drops
# bottom-clipped marks from 13 to 3, because it reaches a metre lower.
#
# FUTURE (noted, not done): clip around the PERCEIVED GROUND SURFACE rather than at fixed
# heights. Both windows are fixed distances from our fitted surface, so in tall canopy the
# denominator is dominated by returns that have nothing to do with where the ground is. A
# window defined relative to the ground return population itself -- rather than -1/+2 or
# -2/+45 m -- would make the rank mean the same thing in a field and under a closed canopy.
# That is a change to what q_ctrl MEASURES, not a filter, so it needs its own comparison
# against these numbers before it could replace them.
_ap.add_argument("--window", choices=["near", "full"], default="full",
                 help="near = ng_all, -1..+2 m at 0.02 m; full = can_all, -2..+45 m at 0.25 m")
_ap.add_argument("--out", default="data/derived/control_q_ctrl_fit.json",
                 help="where to write the relation, so a consumer READS it instead of "
                      "carrying a typed coefficient")
_ap.add_argument("--drop-bottom-clipped", action="store_true",
                 help="drop marks whose surveyed ground lies BELOW the window, whose rank "
                      "therefore clips at 0 and can only flatten the slope")
ARGS = _ap.parse_args()
EDGE_KEY = "ng_edges" if ARGS.window == "near" else "can_edges"
HIST_KEY = "ng_all" if ARGS.window == "near" else "can_all"


def q_at_surveyed(point_id, surveyed_z):
    """Rank of the surveyed ground within this mark's near-ground return column."""
    f = os.path.join(STRUCT, f"{SET}__{point_id}.npz")
    if not os.path.exists(f):
        return np.nan, np.nan
    z = np.load(f)
    coef = z["surface_coef"]
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    h_s = (float(surveyed_z) - float(coef[0])) / nn
    e = z[EDGE_KEY]
    h = z[HIST_KEY].astype(float)
    t = h.sum()
    if t <= 0 or not np.isfinite(h_s):
        return np.nan, np.nan
    # Linear interpolation of the empirical CDF at h_s. Outside the stored window the rank
    # is 0 or 1 and is reported as such -- those marks are counted, not dropped.
    # The rank is not observed exactly; it is INFERRED from a finite sample. Treat the
    # returns as n draws from the underlying height distribution and the count below the
    # surveyed ground as binomial. With the Jeffreys prior Beta(1/2, 1/2) -- the reference
    # prior for a binomial proportion -- the posterior is
    #
    #     p | k, n  ~  Beta(k + 1/2,  n - k + 1/2)
    #
    # and the point estimate is its mean, (k + 1/2) / (n + 1). This needs NO special case
    # for a mark with nothing below the surveyed ground: k = 0 is one value of k, and the
    # posterior handles it exactly as it handles any other. k is kept fractional (the
    # straddling bin contributes its share), which a Beta admits, rather than rounded to
    # force an integer.
    #
    # THE ASSUMPTION: the n returns are treated as independent. They are not. Two things
    # follow, and they are not the same thing, so keep them apart:
    #
    # IDENTICALLY DISTRIBUTED -- assume the vegetation structure is statistically uniform
    #   within the cell / the 7.5 m disc (Andy, 2026-09-03). Then every return, whichever
    #   pulse it came from, is a draw from the SAME height distribution. That is what the
    #   posterior MEAN needs, so the rank estimate is approximately unbiased under it. It
    #   is an assumption about the site, and it fails where a cell straddles a stand edge.
    #
    # INDEPENDENT -- uniformity does not give this, and only the VARIANCE depends on it.
    #   Returns from one pulse are ordered along a single ray. But the effect is bounded
    #   and small, measured on the gen2 cloud these marks come from:
    #
    #     returns 182,923,322   pulses 102,885,312   mean returns per pulse r = 1.7779
    #     number_of_returns:  1: 30.0%   2: 30.9%   3: 26.5%   4: 10.3%   5: 2.0%
    #
    #   If returns within a pulse carried NO independent information -- the worst case --
    #   the effective sample size is n/r and the posterior SD is understated by
    #   sqrt(r) = 1.3334x. That is a bound, and a loose one: a pulse's returns sit at
    #   DIFFERENT heights (canopy, then ground), so they are informative about different
    #   parts of the distribution rather than duplicates. The truth lies between 1.00x and
    #   1.33x on the SD, and nothing on the mean.
    #
    #   Not bounded by that number: spatial correlation BETWEEN pulses within the disc.
    #   Uniformity is the assumption doing the work there, and it is the one to revisit
    #   before trusting these SEs at face value.
    c = np.concatenate([[0.0], np.cumsum(h) / t])
    k = float(np.interp(h_s, e, c)) * t
    a = k + 0.5
    b = t - k + 0.5
    q = a / (a + b)
    sd = np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1.0)))
    return q, h_s * 1000.0, sd


m = load(BAND_LO, BAND_HI)
res = [q_at_surveyed(r.point_id, r.elevation) for r in m.itertuples()]
m["q_ctrl"] = [a for a, _, _ in res]
m["h_survey_mm"] = [b for _, b, _ in res]
m["q_sd"] = [c for _, _, c in res]        # per-mark posterior SD of the rank
m = m.dropna(subset=["q_ctrl"]).copy()
_EPS_RANK = 1e-6      # below this a rank came from the zero-event rule, not from returns

_e0 = np.load(os.path.join(STRUCT, f"{SET}__{m.point_id.iloc[0]}.npz"))[EDGE_KEY]
print(f"window: {ARGS.window}  ({_e0[0]:+.2f} to {_e0[-1]:+.2f} m, "
      f"{_e0[1]-_e0[0]:.3f} m bins, {len(_e0)-1} bins)")
print(f"marks: {len(m):,} with lowveg, a published residual and a locatable surveyed height")
n_lo = int((m.q_ctrl <= 1e-6).sum()); n_hi = int((m.q_ctrl >= 1 - 1e-6).sum())
print(f"  rank from Beta(k+1/2, n-k+1/2) posterior mean; no mark is special-cased")
print(f"  marks with NO return below the surveyed ground: "
      f"{int((m.q_ctrl * (m.q_ctrl.notna()) < 1e-3).sum())} at rank < 1e-3   "
      f"per-mark posterior SD: median {m.q_sd.median():.5f}  max {m.q_sd.max():.5f}")
if ARGS.drop_bottom_clipped and n_lo:
    m = m[m.q_ctrl > 0].copy()
    print(f"  DROPPED the {n_lo} bottom-clipped marks as asked; {len(m):,} remain")
print(f"  q_ctrl  median {m.q_ctrl.median():.4f}   IQR {m.q_ctrl.quantile(.25):.4f}"
      f"-{m.q_ctrl.quantile(.75):.4f}")
print(f"  h_survey median {m.h_survey_mm.median():+.1f} mm   "
      f"(negative = surveyed ground sits BELOW our fitted surface)")

E = np.arange(0, m.lowveg.max() + BIN_W, BIN_W)
m["bin"] = np.digitize(m.lowveg, E) - 1
print(f"\n{'bin':>12s} {'n':>5s} {'mean lowveg':>12s} {'q_ctrl':>8s} {'SE':>7s} "
      f"{'h_survey mm':>12s}")
X, Y, S, NN = [], [], [], []
for b, g in m.groupby("bin"):
    if len(g) < 1:
        continue
    # Total variance of the bin mean = scatter between marks + the marks' own posterior
    # variance, which is what makes this a propagation rather than a substitution.
    se = (np.sqrt(g.q_ctrl.var(ddof=1) / len(g) + (g.q_sd ** 2).sum() / len(g) ** 2)
          if len(g) > 1 else np.nan)
    print(f"  {E[b]:.2f}-{E[b]+BIN_W:<5.2f} {len(g):5d} {g.lowveg.mean():12.4f} "
          f"{g.q_ctrl.mean():8.4f} {se if np.isfinite(se) else float('nan'):7.4f} "
          f"{g.h_survey_mm.mean():12.1f}")
    if len(g) > 1:
        X.append(g.lowveg.mean()); Y.append(g.q_ctrl.mean()); S.append(se); NN.append(len(g))
X, Y, S, NN = map(np.asarray, (X, Y, S, NN))


def wls(x, y, w):
    A = np.vstack([np.ones_like(x), x]).T
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ (A.T @ W @ y)
    r = y - A @ p
    dof = max(len(x) - 2, 1)
    s2 = float((w * r ** 2).sum() / dof)
    err = np.sqrt(np.diag(cov) * (s2 if w.sum() else 1.0))
    return p, err, float((w * r ** 2).sum() / dof)


print(f"\nlinear q_ctrl = a + b*lowveg, on {len(X)} bins with n>1:")
for lab, w in (("DESIGN-weighted 1/SE^2", 1.0 / S ** 2), ("ABUNDANCE-weighted by n", NN * 1.0)):
    p, e, chi2 = wls(X, Y, w)
    print(f"  {lab:24s} intercept {p[0]:+.4f} +/- {e[0]:.4f}   "
          f"slope {p[1]:+.4f} +/- {e[1]:.4f}   chi2/dof {chi2:.2f}")
pm = np.polyfit(m.lowveg, m.q_ctrl, 1)
print(f"  per-mark, unweighted     intercept {pm[1]:+.4f}                slope {pm[0]:+.4f}"
      f"   (n = {len(m):,})")
print("\nFor comparison, the ridgeline route on the same quantity (Elba, lowveg, free "
      "intercept): intercept +0.4853  slope -0.0200")


# ---------------------------------------------------------------------------------------
# The relation as DATA, with the population that produced it attached. A consumer that
# applies these two numbers to a differently-built return column would be using a rank
# measured in one population to index another -- the same class of error as reading a
# percentile off class-2 near-ground returns when it was fitted on the full column. So the
# spec travels with the coefficients and dod_cover_corrected.py refuses on a mismatch.
p_des, e_des, chi2_des = wls(X, Y, 1.0 / S ** 2)
_rel = {
    "relation": "q_ctrl = a + b * lowveg",
    "intercept": float(p_des[0]), "intercept_se": float(e_des[0]),
    "slope": float(p_des[1]), "slope_se": float(e_des[1]),
    "chi2_per_dof": float(chi2_des),
    "weighting": "DESIGN 1/SE^2 on uniform lowveg bins",
    "rank_estimator": "Beta(k+1/2, n-k+1/2) posterior mean, Jeffreys prior",
    "rank_estimator_note": (
        "the rank is inferred, not observed: k returns below the surveyed ground out of n "
        "is binomial, and the Jeffreys prior gives p|k,n ~ Beta(k+1/2, n-k+1/2). A mark "
        "with nothing below is k=0 and needs no rule of its own. Bin SEs carry the marks' "
        "posterior variances as well as the scatter between them. ASSUMES the n returns "
        "are independent, which lidar returns are not. Under the assumption that structure "
        "is statistically uniform within the cell, the returns are at least identically "
        "distributed, so the MEAN is approximately unbiased; only the variance suffers. "
        "Measured on this cloud, r = 1.7779 returns per pulse, so the SD is understated by "
        "at most sqrt(r) = 1.3334x, and by less than that because a pulse's returns sit at "
        "different heights rather than repeating one."),
    "fitted_on": {"set": SET, "marks": int(len(m)), "bins": int(len(X)),
                  "bin_width_lowveg": BIN_W,
                  "lowveg_max_observed": float(m.lowveg.max()),
                  "bottom_clipped_dropped": bool(ARGS.drop_bottom_clipped)},
    # What q_ctrl is a rank WITHIN. A consumer must reproduce this population.
    "percentile_population": {
        "window_lo_m": float(_e0[0]), "window_hi_m": float(_e0[-1]),
        "bin_m": float(_e0[1] - _e0[0]), "n_bins": int(len(_e0) - 1),
        "classes": "ALL returns (no classification filter)",
        "height": "slope-normal to the local ground surface, (z - S)/|n|",
    },
    # What lowveg is, so the covariate is built the same way too.
    "covariate": {
        "name": "lowveg", "band_lo_m": BAND_LO, "band_hi_m": BAND_HI,
        "denominator_lo_m": -1.0, "denominator_hi_m": 2.0, "bin_m": 0.02,
        "classes": "ALL returns (no classification filter)",
        "test": "bin CENTRE, so with 0.02 m bins this is 'fraction above 0.14 m'",
        "source": "analysis/CONTROL_LOWVEG_OFFSET.md",
    },
    # The one difference a tile CANNOT reproduce, stated rather than left to be discovered.
    "known_scale_difference": (
        "at the marks the reference surface is an order-2 least-squares fit to class-2 "
        "returns within 7.5 m of the point; on a tile it is the gridded z_after plane of a "
        "5 m cell. Mark-scale and grid-scale ground estimates are not the same surface, and "
        "this relation was fitted against the former."),
}
json.dump(_rel, open(ARGS.out, "w"), indent=2)
print(f"\nwrote {ARGS.out}  (intercept {_rel['intercept']:+.4f}, slope {_rel['slope']:+.4f})")
