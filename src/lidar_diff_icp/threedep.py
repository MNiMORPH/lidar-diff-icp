"""Automatic discovery of the covering USGS 3DEP (Second-Generation) EPT project
for a tile, with a mandatory coverage check.

The Second-Generation reference is a USGS 3DEP lidar project served as an Entwine
Point Tile (EPT) store on the public ``usgs-lidar-public`` S3 bucket. The hobu/
usgs-lidar repo publishes a GeoJSON index of every public EPT project (name, URL,
point count, and a lon/lat boundary polygon). We use it to (1) find which
project(s) cover a point, (2) rank them so the recent gen2 acquisition is chosen
over an incidental gen1-era reprocessing (e.g. an Arrowhead 2011 collection also
lives on the bucket), and (3) verify the project boundary fully contains the tile
bbox BEFORE downloading -- the check whose absence let a silent tile-count
truncation masquerade as a project-boundary coverage gap.
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

BOUNDARIES_URL = (
    "https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson"
)
# statewide/national mosaics: valid, but vintage varies within them, so they are a
# poor clean-vintage gen2 reference -- rank them last, never auto-pick over a
# specific dated project.
_MOSAIC = re.compile(r"FullState|_National|Mosaic", re.I)


def _load_boundaries(cache: str | Path | None = None) -> dict:
    """Fetch (and optionally cache) the EPT project boundary index."""
    cache = Path(cache) if cache else None
    if cache and cache.exists():
        return json.loads(cache.read_text())
    with urllib.request.urlopen(BOUNDARIES_URL, timeout=180) as r:
        gj = json.load(r)
    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(gj))
    return gj


def _years(name: str) -> list[int]:
    return sorted(int(y) for y in re.findall(r"(?:19|20)\d{2}", name))


def find_projects(lon: float, lat: float, *, cache: str | Path | None = None
                  ) -> list[dict]:
    """3DEP EPT projects whose boundary contains (lon, lat), most recent first.

    Each entry: ``{name, url, count, years, latest, is_mosaic, geom}`` (``geom`` a
    shapely polygon in lon/lat). Ranked by latest name-year descending with
    statewide mosaics pushed last -- so ``[0]`` is the best gen2 candidate, but the
    caller should still confirm the vintage (name-year parsing is heuristic).
    """
    from shapely.geometry import shape, Point
    p = Point(lon, lat)
    out = []
    for f in _load_boundaries(cache)["features"]:
        g = shape(f["geometry"])
        if not g.contains(p):
            continue
        pr = f["properties"]
        name = pr.get("name", "")
        yrs = _years(name)
        out.append(dict(name=name, url=pr.get("url"), count=pr.get("count"),
                        years=yrs, latest=(yrs[-1] if yrs else 0),
                        is_mosaic=bool(_MOSAIC.search(name)), geom=g))
    out.sort(key=lambda d: (not d["is_mosaic"], d["latest"]), reverse=True)
    return out


def bbox_covered(geom, bbox_lonlat: tuple[float, float, float, float]) -> bool:
    """Does a project boundary fully contain a lon/lat bbox (minlon,minlat,maxlon,
    maxlat)? Conservative: the whole rectangle must be inside."""
    from shapely.geometry import box
    return geom.contains(box(*bbox_lonlat))


def resolve_reference(lon: float, lat: float,
                      bbox_lonlat: tuple[float, float, float, float] | None = None,
                      *, cache: str | Path | None = None) -> dict:
    """Pick the gen2 3DEP project for a point and (if a bbox is given) require it
    to fully cover the bbox. Returns the chosen project dict (with a ``covers``
    flag when a bbox was checked). Raises ``LookupError`` if none cover the point,
    or if a bbox is given and no covering project fully contains it.
    """
    cands = find_projects(lon, lat, cache=cache)
    if not cands:
        raise LookupError(f"no 3DEP EPT project covers ({lon}, {lat})")
    if bbox_lonlat is None:
        return cands[0]
    for c in cands:                                   # first (most-recent) full cover
        if bbox_covered(c["geom"], bbox_lonlat):
            return {**c, "covers": True}
    raise LookupError(
        f"a 3DEP project covers ({lon}, {lat}) but none fully covers the bbox "
        f"{bbox_lonlat}; candidates: "
        + ", ".join(f"{c['name']}({c['latest']})" for c in cands))
