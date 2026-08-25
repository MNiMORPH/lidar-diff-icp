"""Deterministic cross-epoch vertical datum from the geoid-model difference.

Two lidar epochs of the same ground can carry different NAVD88 realizations because
they were reduced with different geoid models -- e.g. 2008 MN lidar on GEOID03 vs 2021
3DEP on GEOID18. The orthometric-height datum shift to ADD to the before-epoch (gen1)
so it matches the after-epoch (gen2) is exactly the geoid-undulation difference
``N_before - N_after``: a smooth, nearly planar field over a tile (~sub-mm from a
plane on the Elba pilot). :func:`geoid_difference` GATHERS that field from the PROJ
geoid grids -- no hard-coded constants -- and returns it as ``const + tilt``, the
``geoid_datum`` tuple the pipeline adds to gen1 after the lateral (x,y) registration.

This is the principled datum: get x,y right (Nuth-Kaeaeb lateral shift), then apply the
required geodetic z offset (this). No arbitrary plane is fitted to "stable" surfaces --
messy residuals are left for later analysis, not baked into the datum.

:func:`osm_flat_references` is an optional standalone utility (OSM flat footprints);
it is not part of the datum path.
"""
from __future__ import annotations

import numpy as np


def geoid_difference(bounds, crs, *, before_geoid="us_noaa_geoid03_conus.tif",
                     after_geoid="us_noaa_g2018u0.tif", n=7,
                     proj_data="/usr/share/proj"):
    """Geoid-model datum shift ``(const_m, b, c)`` to ADD to the before-epoch (gen1).

    ``before_geoid`` / ``after_geoid`` are PROJ geoid grid names for the two epochs'
    NAVD88 realizations (defaults: GEOID03 for 2008 gen1, GEOID18 for 2021 gen2). The
    shift is ``N_before - N_after`` (geoid-undulation difference), sampled on an ``n x n``
    grid over ``bounds`` and fit as::

        shift(E,N) = a + b*(E - cx)/1000 + c*(N - cy)/1000       [m; tilt in m/km]

    about the ``bounds`` centroid ``(cx, cy)`` -- the same centroid the pipeline uses --
    so the returned ``(a, b, c)`` is a drop-in ``geoid_datum``. On the Elba pilot this is
    ~(+0.0673, +0.00078, -0.00057), reproducing the independently derived +67 mm.

    ``bounds``=(X0,Y0,X1,Y1) in ``crs`` (EPSG int or CRS). pyproj's data dir is set to
    ``proj_data`` explicitly so this works when PROJ_DATA is unset for GDAL/rasterio --
    pyproj and GDAL do not share a data dir, so this does not disturb rasterio. PROJ grid
    NETWORK access is enabled by default (set ``PROJ_NETWORK=OFF`` to forbid it): the geoid
    grids ship with no PROJ install, and fetching them is what makes the datum reproducible
    on a new machine rather than something to hard-code per tile.
    """
    import os
    import pyproj
    pyproj.datadir.set_data_dir(proj_data)     # pyproj-only; leaves GDAL/rasterio untouched
    # The GEOID03/GEOID18 grids are large and are NOT part of a normal PROJ install, so on a
    # machine without them this dies with "could not find required grid(s)". PROJ can fetch
    # and cache them itself, so network access is enabled by DEFAULT here -- an explicit
    # PROJ_NETWORK=OFF is still honoured for anyone who wants a strictly offline run.
    if os.environ.get("PROJ_NETWORK", "").upper() != "OFF":
        os.environ.setdefault("PROJ_NETWORK", "ON")
        pyproj.network.set_network_enabled(True)
    from pyproj import Transformer
    X0, Y0, X1, Y1 = bounds
    cx = 0.5 * (X0 + X1); cy = 0.5 * (Y0 + Y1)
    xs = np.linspace(X0, X1, n); ys = np.linspace(Y0, Y1, n)
    XX, YY = (a.ravel() for a in np.meshgrid(xs, ys))
    lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(XX, YY)

    def undulation(grid):
        tr = Transformer.from_pipeline(f"+proj=vgridshift +grids={grid} +multiplier=1")
        return np.asarray(tr.transform(lon, lat, np.zeros_like(lon))[2])

    diff = undulation(before_geoid) - undulation(after_geoid)   # N_gen1 - N_gen2 (m)
    A = np.c_[np.ones_like(XX), (XX - cx) / 1000.0, (YY - cy) / 1000.0]
    (a, b, c), *_ = np.linalg.lstsq(A, diff, rcond=None)
    return float(a), float(b), float(c)


def osm_flat_references(bbox_latlon, *, to_epsg=26915, timeout=100):
    """Optional: fetch OSM flat-reference footprints (parking, pitches, tracks,
    cemeteries) for ``bbox_latlon``=(lat0,lon0,lat1,lon1). Needs network + requests +
    pyproj. Returns a list of dicts ``{kind, name, poly (list of [x,y] in to_epsg)}``.
    Standalone utility, not part of the datum path.
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
