"""Terrain masks from a gridded reference surface: slope, curvature, valley, stable.

DEFINED ONCE, HERE. It used to be computed inline inside ``difference_dem``, which meant
there was nothing to import, so five other scripts re-derived it by hand -- and when the
valley cut moved from TPI to elevation on 2026-09-04 only the pipeline's copy was fixed.
``analysis/stable_mask_repro.py`` existed for no other reason than to reproduce this block
by hand; its name was the diagnosis.

TWO DIFFERENT QUESTIONS, TWO DIFFERENT POPULATIONS. This module answers "which ground can
calibrate the INSTRUMENT" -- registration, stable_sigma, the LoD. It wants MANY cells
spanning the slope range, because the lateral (Nuth-Kaeae) shift is estimated from how dh
varies with slope and aspect, so a mask that removes slope removes its signal.
:func:`lidar_diff_icp.refcells.reference_cells` answers a different question -- "which
ground provably did NOT CHANGE" -- for calibrating corrections, and is far stricter
(divides, buildings, clear-cut, gross change). They are not substitutes.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt as _edt, gaussian_filter

from . import coreg

__all__ = ["terrain_masks", "resolve_valley_top"]


def resolve_valley_top(valley_top_m, tile_dir=None):
    """Return ``(elevation_or_None, source_label)``. THE CALLER ALWAYS SAYS WHICH.

    ``valley_top_m`` is an elevation in metres, ``"registry"`` (the established, cited value
    for this tile), or ``"histogram"`` (computed from the LANDSCAPE's pooled elevations, so
    tiles sharing ground share one cut). There is no default: a stated value and a computed
    one are different claims about the same tile, and a chain that silently substitutes one
    for the other hides which was used -- which is how elba once rebuilt on a computed
    226.9 m instead of its established 230.0 m with nothing saying so.
    """
    from .refcells import VALLEY_TOP_M, valley_top_for_landscape
    import os

    if valley_top_m is None:
        raise ValueError(
            "a valley top is required: an elevation in metres, 'registry' for this tile's "
            "established value, or 'histogram' to compute it from the landscape. It will "
            "not be chosen for you.")
    if valley_top_m == "registry":
        if tile_dir is None:
            raise ValueError("'registry' needs tile_dir to know which tile to look up")
        name = os.path.basename(str(tile_dir).rstrip("/"))
        v = VALLEY_TOP_M.get(name)
        if v is None:
            raise ValueError(
                f"no established valley top for {name!r}; add one to refcells.VALLEY_TOP_M "
                f"with its source, or ask for 'histogram'")
        return float(v), "registry"
    if valley_top_m == "histogram":
        v, members = valley_top_for_landscape(tile_dir or ".")
        return (None if v is None else float(v)), f"histogram over {len(members)} tile(s)"
    return float(valley_top_m), "stated"


def terrain_masks(z_grid, res, *, valley_top_m, tile_dir=None, curv_max=0.005,
                  verbose=True):
    """Slope, curvature, valley floor and the STABLE mask, from a gridded ground surface.

    Returns a dict with ``slope_deg``, ``laplacian``, ``floodplain``, ``upland``,
    ``stable``, ``valley_top_m``, ``valley_top_source`` and ``report`` (cells removed by
    each cut, in order), so a caller can record what it used rather than restate it.

    ``stable`` is LOW CURVATURE ABOVE THE VALLEY TOP -- not ridgetops, and with no slope
    restriction. The previous mask, ``(slope<3 & upland) | (5<slope<35 & upland & lap<0)``,
    left whitewater at a median slope of 3.0 deg, nearly flat and close to useless for the
    lateral fit; this keeps p50 ~10 deg and p90 ~28 deg at both tiles it was measured on.

    ``curv_max`` applies to the Laplacian of a 25 m-smoothed surface, which is NOT
    ``curv_laplacian.npy``: refcells' 0.015 sits past p99 here and would cut nothing.
    """
    z = np.asarray(z_grid, float)
    zf = z.copy()
    nanm = np.isnan(zf)
    if nanm.any():
        zf = zf[tuple(_edt(nanm, return_distances=False, return_indices=True))]

    vt, src = resolve_valley_top(valley_top_m, tile_dir)
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(zf, 2.0), res)[0])
    zsm = gaussian_filter(zf, 50 / res / 2)
    lap = (np.gradient(np.gradient(zsm, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(zsm, res, axis=1), res, axis=1))

    report = {"start": int(z.size)}
    if vt is None:
        floodplain = np.zeros(z.shape, bool)
        report["valley cut"] = 0
        if verbose:
            print("  NO valley cut: this landscape's elevation histogram has no minimum "
                  "above its dominant mode, so floodplain cells are NOT excluded from the "
                  "stable set. stable_sigma and the LoD include any valley change.",
                  flush=True)
    else:
        floodplain = np.isfinite(z) & (z < vt)
        report[f"below valley top {vt:.1f} m ({src})"] = int(floodplain.sum())
        if verbose:
            print(f"  valley cut at {vt:.1f} m ({src}): {int(floodplain.sum()):,} cells "
                  f"excluded from the stable set ({100 * floodplain.mean():.1f}% of the "
                  f"grid)", flush=True)
    upland = ~floodplain
    stable = upland & (np.abs(lap) <= curv_max) & np.isfinite(z)
    report[f"|laplacian| > {curv_max:g}"] = int((upland & ~stable & np.isfinite(z)).sum())
    report["kept"] = int(stable.sum())
    if verbose:
        print(f"  stable: |laplacian| <= {curv_max:g} and above the valley top -> "
              f"{int(stable.sum()):,} cells ({100 * stable.mean():.1f}% of the grid); "
              f"slope p50 {np.median(sdeg[stable]):.1f} deg p90 "
              f"{np.percentile(sdeg[stable], 90):.1f}", flush=True)
    return {"slope_deg": sdeg, "laplacian": lap, "filled": zf,
            "floodplain": floodplain, "upland": upland, "stable": stable,
            "valley_top_m": vt, "valley_top_source": src, "report": report}
