"""The stable reference-cell population used to calibrate gen1 against gen2.

Everything that estimates *how* gen1's returns map onto gen2's -- the near-ground rank
gap, the matching quantile, the cover curve -- is calibrated on cells that should not
have changed between epochs. Getting that population right is the whole experiment, so it
is defined ONCE, here, instead of being re-derived per script.

**The governing constraint: the stability test must not look at the vertical offset.**
The quantity being calibrated (a few mm in the open, up to ~130 mm under dense canopy) is
the same size as the LoD (median ~105 mm on these cells). So excluding cells where
``|DoD| > LoD`` -- the standard change detector -- would strip out precisely the
high-cover cells carrying the signal, and the calibration would then "prove" that
vegetation has no effect. That is circular, and it is available here only as an explicit
sensitivity (``gross_change_mm`` set low), never as the default.

Every default criterion is therefore either geometric or vegetation-structural:

* **divide cells, low concavity** (``ridge_mask``, ``|curv_laplacian| <= curv_max``) --
  hilltops shed rather than collect, so they are the geomorphic no-change population.
* **gentle slope** (``slope_max``) -- excludes cells where mass wasting is plausible and
  where intra-cell relief dominates the return column. On these tiles it does not move
  any answer; it trims the noise.
* **no building returns in either epoch** (``n_bldg``) -- structures come and go, and
  their returns are not ground.
* **not in a known retreat zone** (``blufftop_margin_mask``, where the tile has one).
* **no gross change** (``|DoD| > gross_change_mm``, default 500 mm) -- quarries, road
  construction, fills. This DOES look at the offset, but at 25x its scale and ~4x the
  LoD, so it cannot remove a vegetation-driven cell. It bites 0.65-0.67% at Elba.
* **not clear-cut between epochs** (``frac_gen1 - frac_gen2 > clearcut_drop``) -- canopy
  present in 2008 and gone in 2021 means the ground itself was worked over. The test is
  ONE-SIDED on purpose: gen1 is leaf-off November and gen2 leaf-on May, so a deciduous
  stand legitimately gains canopy fraction between epochs. Excluding on ``|delta frac|``
  would strip deciduous stands for their phenology and bias the cover strata.
"""
from __future__ import annotations

import os

import numpy as np

__all__ = ["reference_cells"]


def _opt(d, name):
    p = os.path.join(d, name)
    return np.load(p) if os.path.exists(p) else None


def reference_cells(tile_dir, *, cells=None, curv_max=0.015, slope_max=12.0,
                    gross_change_mm=500.0, clearcut_drop=0.30, require_ridge=True,
                    exclude_valley=True, valley_top_m=None,
                    use_floodplain_mask=True):
    """Boolean mask of stable reference cells, plus a report of what each cut removed.

    ``cells`` is an optional array of flat cell indices (e.g. the near-ground cube's
    ``cells``); the mask is returned over those, otherwise over the whole raveled grid.
    Any criterion whose input the tile does not carry is skipped and recorded as such.

    Returns ``(mask, report)`` where ``report`` maps criterion name -> number of cells it
    removed from the running mask, in application order, with ``"start"`` and ``"kept"``.
    """
    sl = np.load(os.path.join(tile_dir, "slope.npy")).ravel()
    idx = np.arange(sl.size) if cells is None else np.asarray(cells)
    m = np.ones(idx.size, bool)
    rep = {"start": int(m.sum())}

    def cut(name, ok):
        nonlocal m
        before = int(m.sum())
        m = m & ok
        rep[name] = before - int(m.sum())

    ridge = _opt(tile_dir, "ridge_mask.npy")
    if require_ridge and ridge is not None:
        cut("not a divide cell", ridge.astype(bool).ravel()[idx])
    curv = _opt(tile_dir, "curv_laplacian.npy")
    if curv is not None:
        c = np.abs(curv.ravel()[idx])
        cut(f"|curv| > {curv_max:g}", np.isfinite(c) & (c <= curv_max))
    cut(f"slope >= {slope_max:g} deg", np.isfinite(sl[idx]) & (sl[idx] < slope_max))

    g1 = _opt(tile_dir, "gen1_canopy_frac.npz")
    g2 = _opt(tile_dir, "gen2_canopy_frac.npz")
    if g1 is not None and g2 is not None:
        cut("building returns", (g1["n_bldg"].ravel()[idx] == 0)
            & (g2["n_bldg"].ravel()[idx] == 0))
        cut(f"clear-cut (frac drop > {clearcut_drop:g})",
            (g1["frac"].ravel()[idx] - g2["frac"].ravel()[idx]) <= clearcut_drop)

    bluff = _opt(tile_dir, "blufftop_margin_mask.npy")
    if bluff is not None:
        cut("blufftop retreat margin", ~bluff.astype(bool).ravel()[idx])

    # Valley floor. Ridge-ness, curvature and slope do NOT exclude flat valley-bottom
    # ground, so terraces and floodplain enter the "divide" population and dominate it:
    # a 19% valley limb produced an easting gradient of -84.8 mm/km against +3.5 on the
    # upland (analysis/STABLE_POINT_TILT_AUDIT.md), and a 27%-floodplain flat-slope bin
    # produced a spurious +7.9 mm rise. Excluded by default for anything divide-based.
    # Two cuts, neither with an invented number: the crude floodplain mask, and the
    # ANTIMODE of this tile's own elevation histogram, computed here, which separates the
    # upland plateau from the valley terraces.
    if exclude_valley:
        fld = _opt(tile_dir, "floodplain_mask.npy") if use_floodplain_mask else None
        if fld is not None:
            cut("floodplain mask", ~fld.astype(bool).ravel()[idx])
        elif use_floodplain_mask:
            # A MISSING input used to skip this cut in silence, which made two tiles'
            # populations differ by 39,038 cells with nothing to notice it.
            rep["floodplain mask MISSING -- cut NOT applied"] = 0
        zf = _opt(tile_dir, "z_after.npy")
        if zf is not None and valley_top_m is not None:
            # An EXPLICIT valley top, so tiles of the same landscape share one threshold.
            # The per-tile antimode below does not: elba computes 228.9 m and elbaext
            # 237.1 m on overlapping ground, a difference of 31,242 cells at elbaext.
            # 230 m is the established Elba value (analysis/steady_state/
            # run_steady_state_strata.py VALLEY_TOP, and ALLFOREST_BLUFFLAND.md).
            z = zf.ravel()[idx]
            cut(f"below valley top {float(valley_top_m):.1f} m (explicit)",
                ~(np.isfinite(z) & (z < float(valley_top_m))))
        elif zf is not None:
            z = zf.ravel()[idx]
            zc = z[m & np.isfinite(z)]
            if zc.size > 100:
                h, e = np.histogram(zc, bins=60)
                lo, hi = int(np.argmax(h[:30])), 30 + int(np.argmax(h[30:]))
                if hi > lo + 2:                      # bimodal: cut at the antimode
                    anti = lo + int(np.argmin(h[lo:hi]))
                    zthr = 0.5 * (e[anti] + e[anti + 1])
                    cut(f"below elevation antimode {zthr:.1f} m",
                        ~(np.isfinite(z) & (z < zthr)))

    dod = _opt(tile_dir, "dod.npy")
    if dod is not None and gross_change_mm is not None:
        d = dod.ravel()[idx] * 1000.0
        cut(f"|DoD| > {gross_change_mm:g} mm", np.isfinite(d)
            & (np.abs(d) <= gross_change_mm))

    rep["kept"] = int(m.sum())
    return m, rep
