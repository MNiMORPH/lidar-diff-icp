"""Widen the LoD on flat floodplain, where gen1's vegetation state is unrecoverable.

Andy, 2026-09-16: "adjust the LoD on the flat floodplain" -- after we established that
correcting is the less principled move. The reasoning, all measured at elbaext:

WHY NOT CORRECT. A correction `DoD + 0.294 * gen1_spread` is available and fitted, but:

  * it WORKS where the bias is small and FAILS where it is large. On flat floodplain,
    cells whose gen1 returns reach gen2 (48,574) read DoD p50 +1.3 mm and correct cleanly
    to +62.0. Cells where NO gen1 return reaches gen2 (5,097) read -216.2 mm, and the
    correction removes only 76 of those 216 mm -- 35% -- leaving -140.2. Their spread is
    259 mm against the others' 207: the covariate barely responds, because when no pulse
    reaches ground the distribution NARROWS onto the vegetation top.
  * on the 90% that do penetrate, applying it converts "no detectable change" (+1.3 mm)
    into "60 mm of deposition" (+62.0). That conclusion rests entirely on the fit's
    intercept b = +55.0 mm, which cannot be separated from a residual level offset
    (gen2's own absolute level at Elba: -6.6 mm, field sd 31.1 mm).

So the un-applied correction becomes the ERROR TERM instead. We decline to move the DoD and
say instead how far it could be wrong.

THE MAGNITUDE, PER CELL: `m_mm_per_mm * spread`, the correction we are choosing not to
apply, added in QUADRATURE. `m` is READ, never assumed -- it is the fitted slope of
DoD against gen1's own p90-p10 spread on flat floodplain, -0.294 at elbaext (whole
floodplain -0.176; banks dilute it because their spread is terrain).

THE PENETRATION-FAILURE CELLS GET NO DERIVED VALUE, deliberately. Where no gen1 return
reaches gen2 the surface was never measured and NOTHING in gen1 bounds the error: its own
spread is not informative there, as the 35% figure above shows. Their floor is a stated
parameter with no default. Measured, they need ~450 mm to silence 90% of their detections
against ~154 mm for the rest -- but that number comes from the DoD we are judging, so it is
offered as evidence, not adopted as a default.

WHY gen1's SPREAD AND NOT gen2's. gen2 shows NO relation between spread and DoD -- m =
+0.006, r = +0.070 with its full unclassified near-ground column (251 returns/cell,
vegetation included) against gen1's m = -0.294, r = -0.353 (17 returns/cell). If this were
terrain roughness or gridding, the better-sampled epoch would show it too. It does not, so
the effect is gen1-specific. Read as season: gen1 flew November 2008 with floodplain sedge
standing as tall dead biomass, gen2 May 2021 with that flattened and regrowth short.

WHY A SLOPE CUT INSIDE THE ELEVATION CUT. `floodplain_mask` is a plain elevation cut, so it
contains the banks (slope p90 11.2 deg at elbaext). Applying this there would suppress bank
erosion, which is the signal Andy most wants kept. `slope_max_deg` has NO default.
"""
from __future__ import annotations

import numpy as np

__all__ = ["gen1_spread", "penetration_failure", "inflate_lod"]


def gen1_spread(cell, d_mm, shape, *, min_returns=10, lo=0.10, hi=0.90):
    """Per-cell ``p90 - p10`` of gen1's slope-normal near-ground offsets, in mm.

    Cells with fewer than ``min_returns`` get NaN: at ~17 returns per 5 m cell a p10/p90
    pair is already thin, and below 10 it is noise.
    """
    import pandas as pd

    d = np.asarray(d_mm, float)
    ok = np.isfinite(d)
    s = pd.Series(d[ok]).groupby(np.asarray(cell)[ok])
    n = s.size()
    v = (s.quantile(hi) - s.quantile(lo))[n >= min_returns]
    out = np.full(int(np.prod(shape)), np.nan)
    out[v.index.to_numpy()] = v.to_numpy()
    return out.reshape(shape)


def penetration_failure(cell, d_mm, shape, *, min_returns=10):
    """Cells where NO gen1 return reaches gen2's surface: ``min(d) > 0``.

    A yes/no physical test -- no threshold, no percentile, nothing fitted. One penetrating
    return clears a cell, so it flags only cells with genuinely no ground signal. Pooling
    proves these are unrecoverable: in the largest clusters 0.0% of returns reach gen2 even
    pooled over 774 cells and 11,091 shots.
    """
    import pandas as pd

    d = np.asarray(d_mm, float)
    ok = np.isfinite(d)
    s = pd.Series(d[ok]).groupby(np.asarray(cell)[ok])
    n = s.size()
    v = s.min()[n >= min_returns]
    out = np.zeros(int(np.prod(shape)), bool)
    out[v.index.to_numpy()] = v.to_numpy() > 0
    return out.reshape(shape)


def inflate_lod(lod, spread_mm, floodplain, slope_deg, *, slope_max_deg,
                m_mm_per_mm, failure=None, failure_floor_mm=None, verbose=True):
    """Add the declined vegetation correction to the LoD, in quadrature, on flat floodplain.

    Returns ``(lod_new, applied, info)``. Cells outside the mask are untouched, so bank
    erosion keeps the LoD it had.

    A flat-floodplain cell with no spread (too few gen1 returns) is LEFT ALONE rather than
    given a guessed inflation; its count is reported so the gap is visible.
    """
    if slope_max_deg is None:
        raise ValueError(
            "slope_max_deg has no default: it decides how much of the valley stops being "
            "measurable, and the floodplain mask is an elevation cut that includes the "
            "banks. At elbaext 1 deg inflates 7.3% of the valley, 5 deg inflates 26.9%.")
    if m_mm_per_mm is None:
        raise ValueError(
            "m_mm_per_mm has no default: it is the FITTED slope of DoD against gen1 spread "
            "on THIS site's flat floodplain, fitted over EVERY cell there with enough "
            "gen1 returns to measure spread. At elbaext (<=2 deg): -0.483, n=53,671 of "
            "58,767. It is a per-site measurement, not a constant. NOTE -0.294, used "
            "until 2026-09-16, was this fit restricted to the 16,092 cells (30%) that the "
            "near-ground cube built for the return-structure work happens to cover -- a "
            "file footprint, not a chosen population, and it undersized the LoD by ~60%.")
    fp = np.asarray(floodplain, bool)
    flat = fp & (np.asarray(slope_deg, float) <= float(slope_max_deg))
    sp = np.asarray(spread_mm, float)
    add = np.abs(float(m_mm_per_mm)) * sp / 1000.0          # metres

    if failure is not None:
        fail = flat & np.asarray(failure, bool)
        if failure_floor_mm is None:
            raise ValueError(
                f"{int(fail.sum()):,} flat-floodplain cells have NO gen1 return reaching "
                f"gen2, so gen1 never measured that ground and its spread does not bound "
                f"the error -- the correction recovers only 35% of their offset. "
                f"failure_floor_mm has no default and must be stated. Evidence, not a "
                f"recommendation: they need ~450 mm to silence 90% of their detections, "
                f"against ~154 mm for cells that do penetrate; that number comes from the "
                f"DoD being judged, so adopting it silently would be circular.")
        add = np.where(fail, np.maximum(add, float(failure_floor_mm) / 1000.0), add)
    else:
        fail = np.zeros_like(flat)

    have = flat & np.isfinite(add) & np.isfinite(lod)
    out = np.array(lod, float, copy=True)
    out[have] = np.hypot(out[have], add[have])
    info = {"slope_max_deg": float(slope_max_deg), "m_mm_per_mm": float(m_mm_per_mm),
            "failure_floor_mm": failure_floor_mm,
            "cells_flat_floodplain": int(flat.sum()), "cells_inflated": int(have.sum()),
            "cells_penetration_failure": int((fail & have).sum()),
            "cells_flat_without_spread": int((flat & ~have & np.isfinite(lod)).sum()),
            "median_lod_before_mm": float(1000 * np.nanmedian(lod[have])) if have.any() else None,
            "median_lod_after_mm": float(1000 * np.nanmedian(out[have])) if have.any() else None}
    if verbose and have.any():
        print(f"  vegetation LoD: {info['cells_inflated']:,} flat-floodplain cells "
              f"(slope <= {slope_max_deg:g} deg) {info['median_lod_before_mm']:.0f} -> "
              f"{info['median_lod_after_mm']:.0f} mm via |m|={abs(m_mm_per_mm):.3f} x spread; "
              f"{info['cells_penetration_failure']:,} penetration-failure cells floored at "
              f"{failure_floor_mm} mm; {info['cells_flat_without_spread']:,} left alone "
              f"(too few gen1 returns)", flush=True)
    return out, have, info
