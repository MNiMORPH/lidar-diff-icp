"""Apply the pipeline's cross-epoch registration to PER-RETURN offsets.

The DoD pipeline registers gen1 to gen2 in two steps (see :mod:`pipeline`): a constant
lateral shift ``(dx, dy)`` from the order-0 Nuth & Kaeaeb tie, then the geoid-difference
datum added to z. The per-return analysis chain (``gen1_csf_angles`` -> the per-beam
table) instead measures raw gen1 LAS coordinates against the gen2 surface, so its
``d_mm`` carries BOTH uncorrected. This module supplies the two correction terms so a
registration-corrected per-return offset can be built without re-running anything.

Why it matters, and why the lateral term is the dangerous one: the geoid datum is a
constant (plus a sub-mm/km tilt), so it moves an offset distribution bodily and cannot
create structure. The lateral shift does not -- projected onto a slope it becomes

    dz = -(gx*dx + gy*dy)      i.e.  |shift| * tan(slope) * cos(aspect - shift azimuth)

which is a tan(slope) signature of exactly the form a slope-dependent instrument error
would produce, and at Elba-scale shifts (~0.75 m) it reaches hundreds of mm on steep
ground -- larger than the effects being measured. Any f(slope) fitted to uncorrected
per-return offsets is therefore confounded with residual misregistration, and the two
cannot be told apart after the fact. That is the whole reason these columns exist.

NOT included: the per-swath internal alignment (``per_swath_internal_alignment_dxdydz_m``
in the same corrections file), which is per-flight-line and at this site is LARGER than
the cross-epoch shift. Correcting the cross-epoch term alone is an improvement in datum
but not a guarantee of an improvement in slope-dependence; treat corrected and uncorrected
offsets as two readings to compare, not as wrong and right.

Sign convention: the returned terms are ADDED to a slope-normal offset defined as
``d = (z_gen1 - plane_gen2(x, y)) / |n|``, positive = gen1 above gen2.
"""
from __future__ import annotations

import json
import os

import numpy as np

__all__ = ["surface_gradients", "read_cross_epoch_datum", "geoid_term", "lateral_term",
           "corrected_offset"]

_CORRECTION_FILES = ("corrections_geoid.json", "corrections.json")


def surface_gradients(z, res):
    """``(gx, gy, nnorm)`` of a gen2 elevation grid, matching the per-return producer.

    NaNs are filled from the nearest finite cell before differencing (as
    ``gen1_save_angles_slope`` does) so the gradient is defined everywhere the grid is.
    ``nnorm = sqrt(gx^2 + gy^2 + 1)`` converts a vertical difference to a slope-normal one.
    """
    from scipy.ndimage import distance_transform_edt
    zf = np.asarray(z, float).copy()
    miss = ~np.isfinite(zf)
    if miss.any():
        zf = zf[tuple(distance_transform_edt(miss, return_distances=False, return_indices=True))]
    gy, gx = np.gradient(zf, float(res))
    return gx, gy, np.sqrt(gx ** 2 + gy ** 2 + 1.0)


def read_cross_epoch_datum(tile_dir):
    """The tile's ``cross_epoch_datum`` block (geoid const/tilt/centroid + lateral shift).

    Raises FileNotFoundError if the tile carries no corrections file with that block --
    silence would mean shipping an "uncorrected" column that is really an unread file.
    """
    for fn in _CORRECTION_FILES:
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            j = json.load(open(p))
            if "cross_epoch_datum" in j:
                return j["cross_epoch_datum"]
    raise FileNotFoundError(
        f"no cross_epoch_datum in {tile_dir} (looked for {', '.join(_CORRECTION_FILES)})")


def geoid_term(x, y, nnorm_pt, datum):
    """Slope-normal mm to add for the geoid datum: ``(const + tilt.(r - centroid)) / |n|``."""
    gc = float(datum["const_m"])
    gb = float(datum.get("tilt_b_m_per_km", 0.0)); gcc = float(datum.get("tilt_c_m_per_km", 0.0))
    cx, cy = datum.get("centroid", (0.0, 0.0))
    dz = gc + gb * (np.asarray(x, float) - cx) / 1000.0 + gcc * (np.asarray(y, float) - cy) / 1000.0
    return 1000.0 * dz / np.asarray(nnorm_pt, float)


def lateral_term(gx_pt, gy_pt, nnorm_pt, datum):
    """Slope-normal mm to add for the lateral shift: ``-(gx*dx + gy*dy) / |n|``.

    Zero on flat ground by construction, growing as ``|shift| * tan(slope)`` and changing
    sign with aspect -- which is why it can masquerade as, or cancel, a slope-dependent
    instrument error.
    """
    dx, dy = (float(v) for v in datum["horizontal_shift_m"])
    return -1000.0 * (np.asarray(gx_pt, float) * dx + np.asarray(gy_pt, float) * dy) \
        / np.asarray(nnorm_pt, float)


def corrected_offset(d_mm, x, y, gx_pt, gy_pt, nnorm_pt, datum):
    """``(d_corrected_mm, geoid_mm, lateral_mm)`` -- the terms are returned alongside the
    sum so a caller can store them separately and undo either one."""
    g = geoid_term(x, y, nnorm_pt, datum)
    lat = lateral_term(gx_pt, gy_pt, nnorm_pt, datum)
    return np.asarray(d_mm, float) + g + lat, g, lat
