"""Stable flat-reference-surface detection for co-registration datum estimation.

Roofs fail as vertical control on rural tiles: they are small and pitched, so any
horizontal misregistration turns roof slope into an apparent vertical offset, and the
2008 sensor barely resolves them (validated on Elba -- see analysis). The reliable
references are **flat, hard, stable** surfaces -- parking lots, courts, running
tracks, cemeteries, concrete pads, and generally low-slope low-roughness ground:

* **flat** -> no pitch, so no slope x misregistration confound; a robust median of
  (before - after) IS the vertical datum, no plane fitting;
* **hard** (low within-cell roughness) -> densely and accurately sensed in *both*
  epochs, unlike small pitched roofs;
* **stable** -> pavement/concrete/mown-lawn do not erode over the survey interval.

This module PICKS THEM OUT -- geometrically (the workhorse, especially rural) plus
optional OSM labels -- and estimates the vertical datum offset between two epochs on
them. On the Elba pilot the geometric detector yields ~470 hard-flat cells across the
tile and a datum of gen1 - gen2 = -55 mm (SE ~3 mm), consistent across surface
hardness (pavement/compacted/soil) -- far more reliable than roofs. The datum tilt is
better constrained by the independent geoid-model difference (~4 mm/km here) than by
the references themselves (clustered, vulnerable to a resurfaced lot), so this returns
the CONSTANT; keep any linear term to the geoid.
"""
from __future__ import annotations

import numpy as np


def _cell_stats(x, y, z, bounds, res):
    """Per-cell median, robust roughness (IQR->sigma), and count on a `res` grid."""
    import pandas as pd
    X0, Y0, X1, Y1 = bounds
    nx = int((X1 - X0) / res); ny = int((Y1 - Y0) / res)
    ok = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)
    c = ((y[ok] - Y0) / res).astype(int) * nx + ((x[ok] - X0) / res).astype(int)
    g = pd.DataFrame({"c": c, "z": z[ok]}).groupby("c")["z"]
    q = g.quantile([0.25, 0.5, 0.75]).unstack()
    return q[0.5], 0.7413 * (q[0.75] - q[0.25]), g.size(), nx, ny


def flat_hard_cells(bx, by, bz, ax, ay, az, bounds, *, res=2.0, max_slope_deg=4.0,
                    max_rough_m=0.04, min_before=4, min_after=6):
    """Detect flat, hard, both-epochs stable reference cells and their raw offset.

    ``bx,by,bz`` / ``ax,ay,az``: before/after point coordinates (any returns -- on a
    hard flat surface every return is the surface). ``bounds``=(X0,Y0,X1,Y1) CRS units.
    A cell is a reference if: slope < ``max_slope_deg`` (flat), after-epoch within-cell
    roughness < ``max_rough_m`` (hard), and it has >= ``min_before``/``min_after``
    returns (measured in both epochs). Returns a dict of arrays for the kept cells:
    ``x, y, before_z, after_z, roughness, offset`` where ``offset = before_z - after_z``
    (m; the raw datum, gen1 relative to gen2).
    """
    from scipy.ndimage import distance_transform_edt as edt
    import pandas as pd
    X0, Y0, X1, Y1 = bounds
    am, ar, an, nx, ny = _cell_stats(ax, ay, az, bounds, res)
    bm, br, bn, _, _ = _cell_stats(bx, by, bz, bounds, res)
    Z = np.full(nx * ny, np.nan); Z[am.index.values] = am.values; Z = Z.reshape(ny, nx)
    nanm = ~np.isfinite(Z); Zf = Z.copy()
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    gy, gx = np.gradient(Zf, res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy))).ravel()
    common = np.intersect1d(am.index.values, bm.index.values)
    d = pd.DataFrame({"c": common, "am": am.reindex(common).values,
                      "bm": bm.reindex(common).values, "ar": ar.reindex(common).values,
                      "an": an.reindex(common).values, "bn": bn.reindex(common).values,
                      "slope": slope[common]})
    keep = d[(d.slope < max_slope_deg) & (d.ar < max_rough_m) &
             (d.bn >= min_before) & (d.an >= min_after)]
    return dict(x=X0 + (keep.c.values % nx + 0.5) * res,
                y=Y0 + (keep.c.values // nx + 0.5) * res,
                before_z=keep.bm.values, after_z=keep.am.values,
                roughness=keep.ar.values, offset=keep.bm.values - keep.am.values)


def datum_offset(cells, *, keep_frac=1.0):
    """Robust vertical datum from stable flat-hard reference cells.

    Returns a dict: ``dz_before`` (the shift to ADD to the before-epoch so stable
    surfaces match the after-epoch, = -median(offset), m), ``raw`` (median before-after,
    m), ``nmad`` (m), ``se`` (standard error of the median, m), ``n``, and
    ``by_hardness`` (offset by roughness band -- a consistency check: a datum that does
    not depend on surface hardness is real, not a surface artifact). ``keep_frac`` < 1
    trims the most extreme cells (robustness to a resurfaced lot).
    """
    o = np.asarray(cells["offset"], float); r = np.asarray(cells["roughness"], float)
    o = o[np.isfinite(o)]
    if keep_frac < 1.0:
        lo, hi = np.quantile(o, [(1 - keep_frac) / 2, 1 - (1 - keep_frac) / 2])
        o = o[(o >= lo) & (o <= hi)]
    med = float(np.median(o)); nmad = float(1.4826 * np.median(np.abs(o - med)))
    by = {}
    for a, b, lbl in [(0, .025, "pavement<2.5cm"), (.025, .04, "compacted"), (.04, .06, "soil")]:
        s = cells["offset"][(r >= a) & (r < b)]
        if s.size >= 8:
            by[lbl] = (round(1000 * float(np.median(s)), 1), int(s.size))
    return dict(dz_before=-med, raw=med, nmad=nmad, se=nmad / np.sqrt(max(len(o), 1)),
                n=len(o), by_hardness=by)


def osm_flat_references(bbox_latlon, *, to_epsg=26915, timeout=100):
    """Optional: fetch OSM flat-reference footprints (parking, pitches, tracks,
    cemeteries) for ``bbox_latlon``=(lat0,lon0,lat1,lon1). Needs network + requests +
    pyproj. Returns a list of dicts ``{kind, name, poly (list of [x,y] in to_epsg)}``.
    Sparse on rural tiles (use :func:`flat_hard_cells` as the workhorse); rich in towns.
    """
    import requests
    from pyproj import Transformer
    lat0, lon0, lat1, lon1 = bbox_latlon
    bb = f"{lat0},{lon0},{lat1},{lon1}"
    sel = ["way[amenity=parking]", "relation[amenity=parking]", "way[leisure=pitch]",
           "way[leisure=track]", "way[landuse=cemetery]", "relation[landuse=cemetery]"]
    q = "[out:json][timeout:60];(" + "".join(f"{s}({bb});" for s in sel) + ");out geom tags;"
    h = {"User-Agent": "lidar-diff-icp/0.1 (research)"}
    r = requests.get("https://overpass-api.de/api/interpreter", params={"data": q},
                     headers=h, timeout=timeout)
    els = r.json().get("elements", []) if r.text.strip().startswith("{") else []
    t = Transformer.from_crs(4326, to_epsg, always_xy=True)
    out = []
    for e in els:
        tg = e.get("tags", {}); g = e.get("geometry") or []
        if len(g) < 4:
            continue
        kind = ("parking" if tg.get("amenity") == "parking" else
                "cemetery" if tg.get("landuse") == "cemetery" else
                f"pitch:{tg.get('sport', '?')}" if tg.get("leisure") == "pitch" else
                tg.get("leisure", "other"))
        out.append({"kind": kind, "name": tg.get("name", ""),
                    "poly": [list(t.transform(p["lon"], p["lat"])) for p in g]})
    return out
