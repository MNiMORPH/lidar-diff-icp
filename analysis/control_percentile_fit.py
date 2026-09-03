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
# A mark whose surveyed ground sits below every observed return has an EMPIRICAL rank of
# exactly 0, which is not a measurement of "no mass below" -- it is the resolution limit of a
# finite sample. Zero observed in n trials still bounds the true fraction, and the bound
# depends on n. Ignoring that treats a 1,725-return column and a 9,079-return column as
# equally certain, which they are not.
#
# The rules below all answer "what fraction lies below, having seen none in n?", and differ
# only in how much prior mass they place. None is invented here; they are the standard
# zero-event estimators. `ground` additionally uses the LIKELIHOOD OF A GROUND RETURN: only
# ground returns can fall below the true surface, so the informative sample size is the
# mark's own class-2 count, not its total, and the resulting fraction is rescaled back onto
# the full column by the class-2 share.
_ap.add_argument("--below-rule", default="ground",
                 choices=["zero", "jeffreys", "laplace", "rule3", "ground"],
                 help="rank for marks whose surveyed ground lies below every return. "
                      "zero = the empirical 0; "
                      "jeffreys = 0.5/(n+1); laplace = 1/(n+2); rule3 = 3/n (a 95%% upper "
                      "bound, deliberately generous); ground = Jeffreys on the class-2 "
                      "count, rescaled by the class-2 share (DEFAULT). Measured: on "
                      "these marks every rule agrees to 4 decimals, because a column of "
                      "1,725-9,079 returns makes 'none below' a tight constraint -- the "
                      "substituted ranks are 5.5e-05..1.7e-03 against a smallest measured "
                      "bin mean of 0.0222.")
_ap.add_argument("--out", default="data/derived/control_q_ctrl_fit.json",
                 help="where to write the relation, so a consumer READS it instead of "
                      "carrying a typed coefficient")
_ap.add_argument("--drop-bottom-clipped", action="store_true",
                 help="drop marks whose surveyed ground lies BELOW the window, whose rank "
                      "therefore clips at 0 and can only flatten the slope")
ARGS = _ap.parse_args()
EDGE_KEY = "ng_edges" if ARGS.window == "near" else "can_edges"
HIST_KEY = "ng_all" if ARGS.window == "near" else "can_all"


def _below_zero_rank(n_total, n_class2):
    """Plausible rank when NO return was seen below the surveyed ground."""
    if ARGS.below_rule == "zero":
        return 0.0
    if ARGS.below_rule == "jeffreys":
        return 0.5 / (n_total + 1.0)
    if ARGS.below_rule == "laplace":
        return 1.0 / (n_total + 2.0)
    if ARGS.below_rule == "rule3":
        return 3.0 / max(n_total, 1.0)
    # ground: only a ground return can land below the true surface, so condition on those.
    if not n_class2 or not np.isfinite(n_class2) or n_class2 <= 0:
        return 0.5 / (n_total + 1.0)
    frac_below_given_ground = 0.5 / (n_class2 + 1.0)
    return frac_below_given_ground * (n_class2 / max(n_total, 1.0))


def q_at_surveyed(point_id, surveyed_z, n_class2=np.nan):
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
    c = np.concatenate([[0.0], np.cumsum(h) / t])
    q = float(np.interp(h_s, e, c))
    if q <= 0.0:
        q = _below_zero_rank(float(t), float(n_class2))
    return q, h_s * 1000.0


m = load(BAND_LO, BAND_HI)
_nc2 = dict(zip(m.point_id, m.get("n_struct_class2", pd.Series(np.nan, index=m.index))))
res = [q_at_surveyed(r.point_id, r.elevation, _nc2.get(r.point_id, np.nan))
       for r in m.itertuples()]
m["q_ctrl"] = [a for a, _ in res]
m["h_survey_mm"] = [b for _, b in res]
m = m.dropna(subset=["q_ctrl"]).copy()
_EPS_RANK = 1e-6      # below this a rank came from the zero-event rule, not from returns

_e0 = np.load(os.path.join(STRUCT, f"{SET}__{m.point_id.iloc[0]}.npz"))[EDGE_KEY]
print(f"window: {ARGS.window}  ({_e0[0]:+.2f} to {_e0[-1]:+.2f} m, "
      f"{_e0[1]-_e0[0]:.3f} m bins, {len(_e0)-1} bins)")
print(f"marks: {len(m):,} with lowveg, a published residual and a locatable surveyed height")
n_lo0 = int((m.q_ctrl <= 0).sum()) if ARGS.below_rule == "zero" else int(
    (m.q_ctrl > 0).sum() - (m.q_ctrl > _EPS_RANK).sum())
n_lo = int((m.q_ctrl <= 0).sum()); n_hi = int((m.q_ctrl >= 1).sum())
print(f"  surveyed ground below EVERY return: {n_lo0}   above every return: {n_hi}   "
      f"(--below-rule {ARGS.below_rule})")
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
    se = g.q_ctrl.std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.nan
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
    "below_rule": ARGS.below_rule,
    "below_rule_note": (
        "rank assigned where NO return lies below the surveyed ground. 'ground' takes "
        "Jeffreys on the mark's class-2 count -- only a ground return can fall below the "
        "true surface -- rescaled onto the full column by the class-2 share. On these 3 of "
        "389 marks every rule (zero/jeffreys/laplace/rule3/ground) gives the same fit to 4 "
        "decimals; the estimator matters for sparser columns, not for these."),
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
