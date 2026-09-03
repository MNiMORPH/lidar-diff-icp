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
    c = np.concatenate([[0.0], np.cumsum(h) / t])
    return float(np.interp(h_s, e, c)), h_s * 1000.0


m = load(BAND_LO, BAND_HI)
res = [q_at_surveyed(r.point_id, r.elevation) for r in m.itertuples()]
m["q_ctrl"] = [a for a, _ in res]
m["h_survey_mm"] = [b for _, b in res]
m = m.dropna(subset=["q_ctrl"]).copy()

_e0 = np.load(os.path.join(STRUCT, f"{SET}__{m.point_id.iloc[0]}.npz"))[EDGE_KEY]
print(f"window: {ARGS.window}  ({_e0[0]:+.2f} to {_e0[-1]:+.2f} m, "
      f"{_e0[1]-_e0[0]:.3f} m bins, {len(_e0)-1} bins)")
print(f"marks: {len(m):,} with lowveg, a published residual and a locatable surveyed height")
n_lo = int((m.q_ctrl <= 0).sum()); n_hi = int((m.q_ctrl >= 1).sum())
print(f"  surveyed ground BELOW the whole column: {n_lo}   ABOVE it: {n_hi}")
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
