"""The ground percentile as a function of the cell's own class-2 spread.

WHY THIS EXISTS. The pipeline took the per-cell MEDIAN of the ground-class returns,
``ground_q = 0.50``. Measured against 519 surveyed gen2 control marks, that is right where
the ground class is clean and wrong where it is not -- and the spread of the class-2 returns
says which is which, with no cover product, no windows and no external layer. Held out on
10 km spatially blocked folds it removed a +8.1 mm median bias and cut RMS 124.5 -> 104.6 mm.

The full derivation, the shape, and the limits are in
``analysis/GROUND_Q_FROM_CLASS2_SPREAD.md``; the curve is produced by
``analysis/calibrate_ground_q.py``.

WHAT THE CURVE IS. An isotonic (monotone non-increasing) regression of

    rank of the surveyed ground within a mark's class-2 returns

on ``log(class-2 standard deviation in mm)``. Monotone because more contamination cannot
mean a HIGHER ground rank; isotonic rather than a fitted form because the shape -- flat while
the class is no wider than bare-ground noise, falling once it is wider -- should come from
the data rather than from a functional family or a threshold.

PER EPOCH, ALWAYS. A curve is valid only for the epoch it was calibrated on. 2008 was flown
leaf-off in November and 2021 at green-up in May, and the deliveries used different
classifiers. The loader REFUSES to hand back a curve whose recorded epoch does not match the
one asked for, because a silent mismatch here would bias every elevation on the tile.
"""
from __future__ import annotations

import os

import numpy as np

#: Where calibrate_ground_q.py writes its curves.
CURVE_DIR = os.path.join("data", "derived")


def curve_path(epoch):
    """Canonical path for an epoch's curve, e.g. 'gen2_2021_control'."""
    return os.path.join(CURVE_DIR, f"ground_q_vs_class2sd_{epoch}.npz")


def load_curve(path_or_epoch, *, expect_epoch=None):
    """Load a calibration curve, refusing a mismatch rather than silently applying it.

    Accepts a path or an epoch name. ``expect_epoch`` is checked against the epoch recorded
    IN the file, not against the filename, so a renamed file cannot smuggle the wrong
    calibration through.
    """
    p = path_or_epoch if str(path_or_epoch).endswith(".npz") else curve_path(path_or_epoch)
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"no ground-q curve at {p}. It is produced by\n"
            f"    ./lidar-icp/bin/python analysis/calibrate_ground_q.py --set <epoch>\n"
            f"from the surveyed control marks. There is no default curve and no fallback: "
            f"applying an uncalibrated percentile would put an unmeasured bias into every "
            f"elevation on the tile.")
    z = np.load(p, allow_pickle=True)
    got = str(z["set"]) if "set" in z else None
    if expect_epoch is not None and got != expect_epoch:
        raise ValueError(
            f"{p} was calibrated on {got!r}, but {expect_epoch!r} was asked for. A curve is "
            f"valid only for its own epoch -- 2008 is leaf-off November, 2021 is green-up "
            f"May, and the classifiers differ. Calibrate the epoch you mean.")
    return {"log_sd_mm": np.asarray(z["log_sd_mm"], float),
            "q": np.asarray(z["q"], float),
            "epoch": got,
            "n_marks": int(z["n_marks"]) if "n_marks" in z else None,
            "path": p,
            "provenance": {k: str(z[k]) for k in
                           ("fitted_on", "response", "covariate", "shape", "cv",
                            "known_limits") if k in z}}


def q_from_spread(sd_mm, curve, *, min_count=None, count=None):
    """Ground percentile per cell, from that cell's class-2 spread in MILLIMETRES.

    Outside the calibrated range the curve is held at its end values rather than
    extrapolated: an isotonic fit has no meaningful behaviour beyond its data, and a linear
    continuation would run q out of [0, 1] -- the failure that dogged every cover-relation
    version of this correction.

    Returns NaN where the spread is unusable, and where ``count`` is below ``min_count`` if
    both are given. NaN means "this method declines to estimate here", which the caller must
    handle; it is deliberately not the 0.50 default, because falling back silently is how an
    uncorrected cell would end up looking corrected.
    """
    sd = np.asarray(sd_mm, float)
    q = np.full(sd.shape, np.nan)
    ok = np.isfinite(sd) & (sd > 0)
    if count is not None and min_count is not None:
        ok &= np.asarray(count, float) >= min_count
    if ok.any():
        q[ok] = np.interp(np.log(sd[ok]), curve["log_sd_mm"], curve["q"])
    return q


def describe(curve, q=None):
    """One block of text naming the curve and what it did -- for a run's own output."""
    out = [f"ground_q from the class-2 spread: {curve['path']}",
           f"  epoch {curve['epoch']}, calibrated on {curve['n_marks']} surveyed marks"]
    for k in ("shape", "cv", "known_limits"):
        if k in curve["provenance"]:
            out.append(f"  {k}: {curve['provenance'][k]}")
    if q is not None:
        qq = np.asarray(q, float)
        fin = np.isfinite(qq)
        if fin.any():
            out.append(f"  q applied: median {np.nanmedian(qq):.3f}  "
                       f"p10 {np.nanpercentile(qq, 10):.3f}  min {np.nanmin(qq):.3f}; "
                       f"{int(np.sum(qq < 0.45)):,} cells below 0.45, "
                       f"{int((~fin).sum()):,} cells declined")
    return "\n".join(out)
