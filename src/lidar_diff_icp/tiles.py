"""Discovery and retrieval of MnGeo First-Generation (2008-2012) lidar LAZ tiles.

MnGeo stores the statewide First-Generation lidar per county at
``.../elevation/lidar/county/<county>/``, where ``<county>`` is the county name
lowercased with all spaces and punctuation removed (``lesueur``, ``stlouis``,
``lakeofthewoods``). Each county has a ``<county>_tile_list.txt`` and a ``laz/``
directory. Tiles are named ``AAAA-BB-CC`` on a grid skewed to 1:24k-quad
(lat/lon) edges rather than UTM axes, and the LAZ files carry **no embedded CRS**.
Rather than reverse-engineer the naming, we read each tile's LAS public-header
bounding box directly; the server supports HTTP range requests, so a tile's bbox
costs one ~512-byte read instead of downloading the whole ~20 MB file.

Pick the county from a coordinate with :func:`county_for_lonlat` (the county name
is resolved from the point via the FCC Census Area API and verified against the
live county listing). Everything else takes a ``county`` argument.
"""

from __future__ import annotations

import json
import re
import struct
import urllib.request
from pathlib import Path

COUNTY_ROOT = (
    "https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/"
)
DEFAULT_COUNTY = "winona"          # the Elba pilot's county (back-compatible default)

# Byte offsets into the LAS public header block (LAS 1.1-1.4, little-endian).
# Max/Min X, Max/Min Y are six f64 starting at offset 179 (spec order:
# maxx, minx, maxy, miny, maxz, minz).
_BBOX_OFFSET = 179


def _county_urls(county: str) -> tuple[str, str, str]:
    base = f"{COUNTY_ROOT}{county}/"
    return base, base + "laz/", f"{base}{county}_tile_list.txt"


def _get_json(url: str):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def list_counties() -> set[str]:
    """The set of valid county directory names on the MnGeo lidar server.

    This is the authority for the naming convention -- resolve any county name to
    one of these rather than guessing punctuation/spacing.
    """
    with urllib.request.urlopen(COUNTY_ROOT, timeout=60) as r:
        html = r.read().decode(errors="replace")
    return set(re.findall(r'href="([a-z]+)/"', html))


def normalize_county(name: str) -> str:
    """County display name -> MnGeo directory name: lowercase, drop the word
    'county', strip everything but a-z. 'St. Louis County' -> 'stlouis'."""
    return re.sub(r"[^a-z]", "", name.lower().replace("county", ""))


def county_for_lonlat(lon: float, lat: float, *, verify: bool = True) -> str:
    """Resolve the MnGeo county directory name containing a (lon, lat) point.

    Uses the FCC Census Area API (point -> county), normalizes to the directory
    convention, and (default) verifies the name exists in the live county listing.
    Raises ``ValueError`` if the point is not in Minnesota, ``LookupError`` if the
    normalized name is not a valid county directory.

    Note: a point on a county line (e.g. a river that is the boundary) resolves to
    ONE county; the tile covering it may be listed under the neighbor. If
    :func:`find_tile` then finds no containing tile, pass the neighboring county
    explicitly.
    """
    d = _get_json(
        f"https://geo.fcc.gov/api/census/area?lat={lat}&lon={lon}&format=json"
    )
    res = d.get("results") or []
    if not res:
        raise LookupError(f"FCC returned no area for ({lon}, {lat})")
    state = res[0].get("state_name")
    if state != "Minnesota":
        raise ValueError(
            f"({lon}, {lat}) is in {state}, not Minnesota -- this workflow is MN-only")
    county = normalize_county(res[0]["county_name"])
    if verify and county not in list_counties():
        raise LookupError(
            f"resolved county '{county}' ({res[0]['county_name']}) is not a MnGeo "
            f"lidar county directory")
    return county


def list_tiles(county: str = DEFAULT_COUNTY) -> list[str]:
    """Return the ``AAAA-BB-CC`` tile names available for a county."""
    _, _, tile_list_url = _county_urls(county)
    with urllib.request.urlopen(tile_list_url, timeout=60) as r:
        raw = r.read().decode().split()
    return sorted(
        t for t in raw
        if t.count("-") == 2 and all(p.isdigit() for p in t.split("-"))
    )


def header_bbox(name: str, county: str = DEFAULT_COUNTY
                ) -> tuple[float, float, float, float]:
    """(minx, miny, maxx, maxy) for a tile, via a single HTTP range read."""
    _, laz_url, _ = _county_urls(county)
    req = urllib.request.Request(
        f"{laz_url}{name}.laz", headers={"Range": "bytes=0-511"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        b = r.read()
    if b[:4] != b"LASF":
        raise ValueError(f"{name}: not a LAS/LAZ file")
    maxx, minx, maxy, miny = struct.unpack_from("<4d", b, _BBOX_OFFSET)
    return minx, miny, maxx, maxy


def _cache_path(cache: str | Path | None, county: str) -> Path | None:
    """Per-county cache path so indices for different counties never collide."""
    if cache is None:
        return None
    cache = Path(cache)
    return cache.with_suffix(f".{county}{cache.suffix or '.json'}")


def build_index(county: str = DEFAULT_COUNTY,
                cache: str | Path | None = None) -> dict[str, list[float]]:
    """Read every tile's bbox once for a county and (optionally) cache to JSON.

    Returns a mapping ``name -> [minx, miny, maxx, maxy]``. The cache is keyed by
    county so reruns read headers only once per county.
    """
    cpath = _cache_path(cache, county)
    if cpath and cpath.exists():
        return json.loads(cpath.read_text())
    index: dict[str, list[float]] = {}
    for name in list_tiles(county):
        try:
            index[name] = list(header_bbox(name, county))
        except Exception as exc:  # transient network / missing tile
            print(f"  skip {name}: {exc}")
    if cpath:
        cpath.parent.mkdir(parents=True, exist_ok=True)
        cpath.write_text(json.dumps(index))
    return index


LIDAR_ROOT = COUNTY_ROOT.split("county/")[0]            # .../elevation/lidar/
STATE_TILE_INDEX = LIDAR_ROOT + "tile_index/indx_q006kpy4.gdb.zip"
DEFAULT_TILE_INDEX_CACHE = "data/mn_tile_centroids.csv"


def _find_ogr2ogr():
    import glob, os, shutil
    found = shutil.which("ogr2ogr")
    if found:
        return found
    for pat in ("~/anaconda3/bin/ogr2ogr", "~/miniconda3/bin/ogr2ogr",
                "/opt/conda/bin/ogr2ogr"):
        for g in glob.glob(os.path.expanduser(pat)):
            return g
    raise FileNotFoundError("ogr2ogr (GDAL) not found -- needed to build the tile index")


def centroid_index(cache: str | Path = DEFAULT_TILE_INDEX_CACHE):
    """Statewide tile name -> UTM 15N centroid, from MnGeo's authoritative lidar
    tile index (``indx_q006kpy4.gdb``, whose ``DNR_QQQ_ID`` equals our AAAA-BB-CC
    tile name). Returns ``(names, xy)`` with ``xy`` an N x 2 float array. Built
    once via GDAL from the remote geodatabase (one range-read, not a per-tile
    header scan) and cached to a small CSV -- the state index is our friend, and
    this is gentle on the server.
    """
    import numpy as np
    cache = Path(cache)
    if not cache.exists():
        import subprocess
        cache.parent.mkdir(parents=True, exist_ok=True)
        url = f"/vsizip//vsicurl/{STATE_TILE_INDEX}/indx_q006kpy4.gdb"
        subprocess.run([_find_ogr2ogr(), "-f", "CSV", str(cache), url,
                        "indx_q006kpt4", "-select", "DNR_QQQ_ID,x,y"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    names, xs, ys = [], [], []
    with open(cache) as fh:
        next(fh)                                          # header line
        for line in fh:
            p = line.rstrip("\n").split(",")
            if len(p) >= 3 and p[1] and p[2]:
                names.append(p[0]); xs.append(float(p[1])); ys.append(float(p[2]))
    return names, np.column_stack([xs, ys])


def find_tile(easting: float, northing: float, county: str | None = None,
              cache: str | Path = DEFAULT_TILE_INDEX_CACHE) -> str:
    """Return the ``AAAA-BB-CC`` tile whose grid cell contains a UTM 15N (E, N)
    point, via the statewide centroid index -- the nearest tile centroid IS the
    containing cell of this regular grid. One cached index, correct across county
    lines, and NO per-tile header hammering (which self-throttles on big metro
    counties, so throttled 404s were misread as missing tiles). ``county`` is
    accepted but ignored (kept for backward compatibility).

    Returns the DNR_QQQ_ID tile name, which IS the LAZ filename in non-metro
    counties. **Metro exception:** the Metro-2011 counties (Ramsey, Hennepin, ...)
    store the same tile under a suffixed name (e.g. ``4342-03-32_b_a.laz``) in
    their ``laz/`` directory, given by the ``las_tile_name`` field of the county's
    own ``elevation_data.gdb`` ``tile_index`` layer (which also carries the
    footprint polygons). For those, resolve the filename from the county gdb.
    """
    import numpy as np
    names, xy = centroid_index(cache)
    d2 = (xy[:, 0] - easting) ** 2 + (xy[:, 1] - northing) ** 2
    return names[int(d2.argmin())]


def download_tile(name: str, out_dir: str | Path,
                  county: str = DEFAULT_COUNTY, *, tries: int = 4) -> Path:
    """Download a full tile LAZ from a county to ``out_dir``; return the path.
    Retries transient failures (throttling/5xx) with backoff so a rate-limited
    response is not mistaken for a missing tile."""
    import time
    import urllib.error
    _, laz_url, _ = _county_urls(county)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{name}.laz"
    url = f"{laz_url}{name}.laz"
    for k in range(tries):
        try:
            urllib.request.urlretrieve(url, dest)
            return dest
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and k < tries - 1:
                time.sleep(3 * (k + 1)); continue
            raise
        except Exception:
            if k < tries - 1:
                time.sleep(3 * (k + 1)); continue
            raise
    return dest
