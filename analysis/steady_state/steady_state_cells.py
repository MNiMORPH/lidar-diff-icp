"""Geomorphically-motivated steady-state stable cells for DoD / datum checking.

Concept
-------
For linear hillslope diffusion  dz/dt = K * grad^2(z), a cell with ZERO
planform Laplacian curvature (grad^2 z = 0, i.e. locally planar) has zero flux
divergence: it receives as much sediment as it delivers.  In a landscape near
topographic steady state such a cell's surface elevation should not change over
time, PROVIDED it is not mass-wasting.  Restricting to slope < max_slope (default
15 deg, conservative) excludes mass-wasting cells.

These planar cells are therefore an INDEPENDENT, geomorphically-motivated set of
"expected-zero-change" cells.  The gen2 - gen1 elevation difference (DoD) over
them should center on ~0 mm if
    (a) the landscape is near steady state, and
    (b) the vertical datum tie between epochs is correct.
It is a check on the DoD / datum that does NOT rely on hard (built) surfaces.

CAVEAT (assumption, not fact): topographic steady state is ASSUMED.  A landscape
that is net-aggrading or net-incising will move planar cells coherently and this
test cannot separate that real signal from a datum offset.  On soil-mantled
hillslopes over a ~13 yr interval the expected steady-state departure is tiny
(see the diffusive-signal budget below), so the assumption is defensible here,
but it remains an assumption.

This module is data-agnostic: pass in DEM / curvature / slope / cover-mask arrays
and thresholds; it returns masks, extracted elevations, and PDF statistics.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# robust statistics helpers
# ---------------------------------------------------------------------------
def nmad(x):
    """Normalized median absolute deviation (robust std estimator)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def diff_stats(diff):
    """Robust + classical summary of a difference sample (metres in -> mm out)."""
    diff = np.asarray(diff, dtype=float)
    diff = diff[np.isfinite(diff)]
    if diff.size == 0:
        return dict(n=0, median_mm=np.nan, nmad_mm=np.nan, mean_mm=np.nan,
                    std_mm=np.nan, iqr_mm=np.nan, q25_mm=np.nan, q75_mm=np.nan)
    q25, q75 = np.percentile(diff, [25, 75])
    return dict(
        n=int(diff.size),
        median_mm=float(np.median(diff) * 1e3),
        nmad_mm=float(nmad(diff) * 1e3),
        mean_mm=float(np.mean(diff) * 1e3),
        std_mm=float(np.std(diff) * 1e3),
        iqr_mm=float((q75 - q25) * 1e3),
        q25_mm=float(q25 * 1e3),
        q75_mm=float(q75 * 1e3),
    )


# ---------------------------------------------------------------------------
# eps_curv selection (principled)
# ---------------------------------------------------------------------------
def eps_curv_from_quantile(curv, base_mask, central_frac=0.30):
    """Symmetric |curvature| band that keeps the central `central_frac` of the
    curvature distribution over `base_mask`.

    Isolates the genuinely planar cells: |grad^2 z| below the
    `central_frac`-quantile of |curvature|.  Returns the eps value (1/m).
    """
    c = np.abs(curv[base_mask])
    c = c[np.isfinite(c)]
    return float(np.percentile(c, central_frac * 100.0))


def eps_curv_from_diffusion(K, dt, signal_budget_m=0.005):
    """|curvature| below which the linear-diffusion elevation change over `dt`
    stays under `signal_budget_m`:  K * eps * dt < budget  ->  eps = budget/(K*dt).

    K in m^2/yr, dt in yr, budget in m; returns eps in 1/m.
    """
    return float(signal_budget_m / (K * dt))


# ---------------------------------------------------------------------------
# core selection
# ---------------------------------------------------------------------------
def steady_state_mask(curv, slope_deg, dod, cover_mask=None,
                      eps_curv=None, max_slope=15.0):
    """Boolean mask of geomorphic steady-state cells.

    A cell qualifies when ALL hold:
        |grad^2 z| < eps_curv        (locally planar -> zero flux divergence)
        slope < max_slope [deg]      (not mass-wasting)
        cover_mask is True           (e.g. core forest), if supplied
        curv, slope, dod all finite  (usable data)

    Parameters
    ----------
    curv, slope_deg, dod : 2-D arrays, same shape
    cover_mask : bool array or None
    eps_curv : float (1/m).  Required.
    max_slope : float, degrees

    Returns
    -------
    mask : bool array
    """
    if eps_curv is None:
        raise ValueError("eps_curv must be supplied (see eps_curv_from_* helpers)")
    curv = np.asarray(curv, float)
    slope_deg = np.asarray(slope_deg, float)
    dod = np.asarray(dod, float)
    finite = np.isfinite(curv) & np.isfinite(slope_deg) & np.isfinite(dod)
    mask = finite & (np.abs(curv) < eps_curv) & (slope_deg < max_slope)
    if cover_mask is not None:
        mask &= np.asarray(cover_mask, bool)
    return mask


def extract_diff(dod, mask):
    """gen2 - gen1 elevation difference (the DoD) over `mask`, finite values only.

    The DoD *is* the gen2 - gen1 difference; gen1 elevation at a cell is
    z_after - dod.  For the steady-state test we only need the difference, which
    is exactly the DoD sampled on the selected cells.
    """
    return dod[mask][np.isfinite(dod[mask])]


def extract_elevations(z_after, dod, mask):
    """Return (gen1, gen2, diff) elevation samples over `mask`.

    gen2 = z_after ; gen1 = z_after - dod ; diff = gen2 - gen1 = dod.
    """
    m = mask & np.isfinite(z_after) & np.isfinite(dod)
    gen2 = z_after[m]
    gen1 = z_after[m] - dod[m]
    return gen1, gen2, gen2 - gen1
