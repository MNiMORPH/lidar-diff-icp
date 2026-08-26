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

The corrections come in two families, and BOTH matter:

* CROSS-EPOCH (gen1 vs gen2): the geoid datum and the constant lateral tie.
* INTERNAL (gen1 vs itself): the per-swath alignment to the lowest-numbered flight line,
  and the per-swath along-track GNSS drift. Without these the point returns carry raw
  swath-to-swath disagreement -- the same ground measured twice from two flight lines
  differs by up to 1.4 m laterally and ~44 mm vertically here -- so the cloud is not even
  internally consistent, and any statistic pooled across swaths mixes that in.

The internal terms are also what makes two tiles comparable: the alignment is relative to
each tile's OWN lowest swath, so tiles built from different swath sets sit on different
gauges (at this site the elba and elbaext lateral ties differ by 197 mm, of which all but
25 mm is gauge). Applying the alignment puts them on a common frame.

Terms compose additively: every shift is small enough that the gen2 surface is linear over
it, so the total is the sum of the individual terms and each can be stored and undone
separately. NOT included, and deliberately so: the boresight-RESIDUAL roll. The vendor's TerraMatch
boresight is already applied in the delivered 2008 data; our opt-in term searches for a
RESIDUAL roll on top of it, and when that search was run it found nothing resolvable --
hence ``boresight_roll_mm_per_deg: None`` in the corrections files. That is a settled
negative result, not an outstanding correction to add later.

Sign convention: the returned terms are ADDED to a slope-normal offset defined as
``d = (z_gen1 - plane_gen2(x, y)) / |n|``, positive = gen1 above gen2.
"""
from __future__ import annotations

import json
import os

import numpy as np

__all__ = ["surface_gradients", "read_cross_epoch_datum", "read_swath_alignment",
           "read_drift_curves", "geoid_term", "lateral_term", "swath_alignment_term",
           "along_track_drift_term", "corrected_offset", "registration_terms"]

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


def read_corrections(tile_dir):
    """The tile's corrections sidecar, by the project-wide precedence.

    ``corrections_geoid.json`` wins over ``corrections.json`` where both exist: a tile
    rebuilt onto the geoid-only datum writes the former and LEAVES the older
    ``reference_plane`` product in place beside it (elbaext carries both). Reading
    ``corrections.json`` by name therefore silently picks the obsolete datum on exactly
    the tiles that have been brought up to date -- which is the bug this exists to stop.

    Precedence is by FILENAME, matching :func:`read_cross_epoch_datum` and the rest of the
    codebase. That is a convention, not a check on content; if a tile ever writes a
    non-geoid datum to the ``_geoid`` name this picks the wrong one. Callers that care
    should assert on ``cross_epoch_datum["method"]``.
    """
    for fn in _CORRECTION_FILES:
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            return json.load(open(p))
    raise FileNotFoundError(
        f"no corrections file in {tile_dir} (looked for {', '.join(_CORRECTION_FILES)})")


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


def _load(tile_dir, key):
    for fn in _CORRECTION_FILES:
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            j = json.load(open(p))
            if key in j:
                return j[key]
    raise FileNotFoundError(f"no {key} in {tile_dir} (looked for {', '.join(_CORRECTION_FILES)})")


def read_swath_alignment(tile_dir):
    """``{point_source_id: (dx, dy, dz)}`` -- the internal alignment of each flight line
    to the lowest-numbered one (which is the reference and maps to zeros)."""
    return {int(k): tuple(float(c) for c in v)
            for k, v in _load(tile_dir, "per_swath_internal_alignment_dxdydz_m").items()}


def read_drift_curves(tile_dir):
    """``{point_source_id: (gps_time, drift_m)}`` -- the per-swath along-track GNSS drift
    curve, to be interpolated in gps_time within each swath."""
    out = {}
    for k, v in _load(tile_dir, "along_track_drift_gpsTime_to_m").items():
        t = np.asarray(v["gps_time"], float); d = np.asarray(v["drift_m"], float)
        order = np.argsort(t)
        out[int(k)] = (t[order], d[order])
    return out


def swath_alignment_term(psid, gx_pt, gy_pt, nnorm_pt, align):
    """Slope-normal mm to add for the per-swath internal alignment.

    Each flight line gets its own ``(dx, dy, dz)``: the vertical part enters directly and
    the lateral part through the surface gradient, exactly as the cross-epoch shift does.
    An unmapped point_source_id raises rather than silently contributing zero -- a swath
    quietly left uncorrected is the failure this whole module exists to prevent.
    """
    psid = np.asarray(psid)
    gx_pt = np.asarray(gx_pt, float); gy_pt = np.asarray(gy_pt, float)
    nnorm_pt = np.asarray(nnorm_pt, float)
    missing = set(np.unique(psid).tolist()) - set(align)
    if missing:
        raise KeyError(f"no alignment for point_source_id {sorted(missing)}; "
                       f"known: {sorted(align)}")
    dx = np.empty(psid.shape, float); dy = np.empty(psid.shape, float); dz = np.empty(psid.shape, float)
    for s, (ax, ay, az) in align.items():
        m = psid == s
        if m.any():
            dx[m] = ax; dy[m] = ay; dz[m] = az
    return 1000.0 * (dz - (gx_pt * dx + gy_pt * dy)) / nnorm_pt


def along_track_drift_term(psid, gps_time, nnorm_pt, curves):
    """Slope-normal mm to add for the per-swath along-track GNSS drift.

    The drift is a vertical, time-varying, per-flight-line term: each return is
    interpolated on its own swath's curve (clamped to the curve ends outside its span).
    """
    psid = np.asarray(psid); gps_time = np.asarray(gps_time, float)
    nnorm_pt = np.asarray(nnorm_pt, float)
    missing = set(np.unique(psid).tolist()) - set(curves)
    if missing:
        raise KeyError(f"no drift curve for point_source_id {sorted(missing)}; "
                       f"known: {sorted(curves)}")
    dz = np.zeros(psid.shape, float)
    for s, (t, d) in curves.items():
        m = psid == s
        if m.any():
            dz[m] = np.interp(gps_time[m], t, d)
    return 1000.0 * dz / nnorm_pt


def registration_terms(d_mm, x, y, gps_time, psid, gx_pt, gy_pt, nnorm_pt, tile_dir):
    """Every registration term for a per-return offset, plus their sum.

    Returns a dict with ``geoid``, ``lateral``, ``swath``, ``drift`` (each slope-normal mm)
    and ``d_corr`` = ``d_mm`` + all four -- the offset as the DoD pipeline would measure it.
    Terms are kept separate so any one can be inspected, excluded, or undone.
    """
    datum = read_cross_epoch_datum(tile_dir)
    geoid = geoid_term(x, y, nnorm_pt, datum)
    lateral = lateral_term(gx_pt, gy_pt, nnorm_pt, datum)
    swath = swath_alignment_term(psid, gx_pt, gy_pt, nnorm_pt, read_swath_alignment(tile_dir))
    drift = along_track_drift_term(psid, gps_time, nnorm_pt, read_drift_curves(tile_dir))
    return {"geoid": geoid, "lateral": lateral, "swath": swath, "drift": drift,
            "d_corr": np.asarray(d_mm, float) + geoid + lateral + swath + drift}
