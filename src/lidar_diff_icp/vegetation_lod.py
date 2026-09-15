"""Inflate the LoD where gen1 could not see the ground through vegetation.

Andy, 2026-09-15: "Let's do this to the LoD for the floodplain because we cannot measure
vegetation. This will still give us bank erosion, which is good."

WHY AN ERROR TERM AND NOT A CORRECTION. gen1 cannot characterise its own vegetation: 16
returns per 5 m cell in the -1..+2 m band against gen2's 224, 97.8% single-return, and every
gen1 structure proxy is dead against the percentile that would fix it (|corr| <= 0.20; only
gen2's spread carries signal, +0.476, and gen2 is May 2021 while the bias is November 2008).
So the floodplain surface in 2008 is NOT MEASURED, and a 78 mm LoD over it is false
precision. This says so in the LoD instead of implying a correction exists.

THE MAGNITUDE IS gen1's OWN, PER CELL: ``lift = median - p05`` of that cell's slope-normal
offsets. It is how far gen1's median sits above its own lowest returns -- an estimate, from
gen1 alone, of how much vegetation could be holding the surface up. Added in QUADRATURE, so
it is an extra independent uncertainty rather than a replacement.

WHY A SLOPE CUT IS STILL NEEDED, having tried to avoid one. The lift does NOT separate
vegetation from terrain: measured at elbaext it is 149 mm in the vegetated patches and
150 mm on the valley banks, a ratio of 0.99. Slope-normal removes the PLANE, but a bank's
within-cell curvature and roughness spread the returns just as much as grass does. So the
lift is a good MAGNITUDE and a useless DISCRIMINATOR, and where to apply it has to come
from terrain. Applying it everywhere would inflate the banks and destroy exactly the bank
erosion this is meant to preserve.

``slope_max_deg`` HAS NO DEFAULT. It decides how much of the valley stops being measurable,
and that is a scientific choice. Measured at elbaext, of 45,954 detections inside the valley
mask:

    cut     cells inflated (% of valley)    detections spared on steeper ground
    1 deg   18,983  ( 7.3%)                 39,851  (86.7%)
    2 deg   36,243  (13.9%)                 35,447  (77.1%)
    3 deg   49,287  (18.9%)                 32,433  (70.6%)
    5 deg   70,098  (26.9%)                 28,374  (61.7%)

THE FLOODPLAIN MASK IS AN ELEVATION CUT and is known to be imperfect: at elbaext it is
everything below the 230 m valley top, so it climbs the valley sides (slope p10 1.4, p50
12.8, p90 28.3 deg). Andy, 2026-09-15: "if you are using your current floodplain cut (which
goes partway up the hills), this is okay for now with a note to correct later." THIS IS THAT
NOTE. A hypsometric definition would be better; the slope cut here is doing work the
floodplain definition should be doing itself.
"""
from __future__ import annotations

import numpy as np

__all__ = ["gen1_lift", "inflate_lod"]


def gen1_lift(cell, d_mm, shape, *, min_returns=8, q_low=0.05):
    """Per-cell ``median - p05`` of gen1's slope-normal offsets, in mm.

    ``q_low`` is the low tail taken to stand for ground. Cells with fewer than
    ``min_returns`` get NaN rather than a number from too few samples -- at 16 returns per
    cell a p05 is already thin, and below 8 it is meaningless.
    """
    import pandas as pd

    d = np.asarray(d_mm, float)
    ok = np.isfinite(d)
    s = pd.Series(d[ok]).groupby(np.asarray(cell)[ok])
    n = s.size()
    lift = (s.quantile(0.50) - s.quantile(q_low))[n >= min_returns]
    out = np.full(int(np.prod(shape)), np.nan)
    out[lift.index.to_numpy()] = lift.to_numpy()
    return out.reshape(shape)


def inflate_lod(lod, lift_mm, floodplain, slope_deg, *, slope_max_deg, verbose=True):
    """Add the vegetation lift in quadrature to the LoD, on FLAT floodplain only.

    Returns ``(lod_new, applied_mask, info)``. Cells outside the mask are untouched, so
    bank erosion on the valley sides keeps the LoD it had.

    A cell inside the mask with no lift (too few gen1 returns) is left ALONE rather than
    given a guessed inflation: it is already a cell we know little about, and inventing an
    error bar for it would be the same false precision in the other direction. Their count
    is returned so the gap is visible.
    """
    if slope_max_deg is None:
        raise ValueError(
            "slope_max_deg has no default. It decides how much of the valley stops being "
            "measurable -- at elbaext 1 deg inflates 7.3% of the valley and spares 86.7% "
            "of its detections, 5 deg inflates 26.9% and spares 61.7%. State it.")
    fp = np.asarray(floodplain, bool)
    flat = fp & (np.asarray(slope_deg, float) <= float(slope_max_deg))
    have = flat & np.isfinite(lift_mm) & np.isfinite(lod)
    out = np.array(lod, float, copy=True)
    out[have] = np.hypot(out[have], np.asarray(lift_mm, float)[have] / 1000.0)
    info = {"slope_max_deg": float(slope_max_deg),
            "cells_in_floodplain": int(fp.sum()),
            "cells_flat_floodplain": int(flat.sum()),
            "cells_inflated": int(have.sum()),
            "cells_flat_without_lift": int((flat & ~have & np.isfinite(lod)).sum()),
            "median_lod_before_mm": float(1000 * np.nanmedian(lod[have])) if have.any() else None,
            "median_lod_after_mm": float(1000 * np.nanmedian(out[have])) if have.any() else None,
            "median_lift_mm": float(np.nanmedian(np.asarray(lift_mm)[have])) if have.any() else None}
    if verbose and have.any():
        print(f"  vegetation LoD: {info['cells_inflated']:,} flat-floodplain cells "
              f"(slope <= {slope_max_deg:g} deg) inflated "
              f"{info['median_lod_before_mm']:.0f} -> {info['median_lod_after_mm']:.0f} mm "
              f"by a median lift of {info['median_lift_mm']:.0f} mm; "
              f"{info['cells_flat_without_lift']:,} flat cells had too few gen1 returns "
              f"for a lift and were LEFT ALONE", flush=True)
    return out, have, info
