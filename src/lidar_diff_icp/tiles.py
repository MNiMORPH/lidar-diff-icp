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


def find_tile(easting: float, northing: float, county: str = DEFAULT_COUNTY,
              cache: str | Path | None = None) -> str:
    """Return the county tile whose bbox contains a UTM 15N (E, N) point.

    Raises ``LookupError`` if no tile in ``county`` contains the point (which can
    happen for a point on a county boundary -- try the neighboring county). If
    more than one bbox contains it (bboxes overlap slightly due to the grid skew),
    the tightest containing tile is returned.
    """
    index = build_index(county, cache)
    hits = [
        (name, bb) for name, bb in index.items()
        if bb[0] <= easting <= bb[2] and bb[1] <= northing <= bb[3]
    ]
    if not hits:
        raise LookupError(
            f"no tile in county '{county}' contains ({easting}, {northing}); "
            f"if this point is on a county line, pass the neighboring county")
    # tightest = smallest area, in case of overlapping bboxes at an edge
    name, _ = min(hits, key=lambda kv: (kv[1][2] - kv[1][0]) * (kv[1][3] - kv[1][1]))
    return name


def download_tile(name: str, out_dir: str | Path,
                  county: str = DEFAULT_COUNTY) -> Path:
    """Download a full tile LAZ from a county to ``out_dir``; return the path."""
    _, laz_url, _ = _county_urls(county)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{name}.laz"
    urllib.request.urlretrieve(f"{laz_url}{name}.laz", dest)
    return dest
