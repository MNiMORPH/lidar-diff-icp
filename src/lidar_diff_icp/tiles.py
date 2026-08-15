"""Discovery and retrieval of MnGeo Winona County 2008 lidar LAZ tiles.

The tiles are named ``AAAA-BB-CC`` on a grid skewed to 1:24k-quad (lat/lon)
edges rather than UTM axes, and the LAZ files carry **no embedded CRS**. Rather
than reverse-engineer the naming, we read each tile's LAS public-header bounding
box directly. The MnGeo server supports HTTP range requests, so a tile's bbox
costs one ~512-byte read instead of downloading the whole ~20 MB file.
"""

from __future__ import annotations

import json
import struct
import urllib.request
from pathlib import Path

BASE_URL = (
    "https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/winona/"
)
LAZ_URL = BASE_URL + "laz/"
TILE_LIST_URL = BASE_URL + "winona_tile_list.txt"

# Byte offsets into the LAS public header block (LAS 1.1-1.4, little-endian).
# Max/Min X, Max/Min Y are six f64 starting at offset 179 (spec order:
# maxx, minx, maxy, miny, maxz, minz).
_BBOX_OFFSET = 179


def list_tiles() -> list[str]:
    """Return the ``AAAA-BB-CC`` tile names available for Winona County."""
    with urllib.request.urlopen(TILE_LIST_URL, timeout=60) as r:
        raw = r.read().decode().split()
    return sorted(
        t for t in raw
        if t.count("-") == 2 and all(p.isdigit() for p in t.split("-"))
    )


def header_bbox(name: str) -> tuple[float, float, float, float]:
    """(minx, miny, maxx, maxy) for a tile, via a single HTTP range read."""
    req = urllib.request.Request(
        f"{LAZ_URL}{name}.laz", headers={"Range": "bytes=0-511"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        b = r.read()
    if b[:4] != b"LASF":
        raise ValueError(f"{name}: not a LAS/LAZ file")
    maxx, minx, maxy, miny = struct.unpack_from("<4d", b, _BBOX_OFFSET)
    return minx, miny, maxx, maxy


def build_index(cache: str | Path | None = None) -> dict[str, list[float]]:
    """Read every tile's bbox once and (optionally) cache to JSON.

    Returns a mapping ``name -> [minx, miny, maxx, maxy]``. Reuses the cache
    file if present so the ~220 header reads happen only once.
    """
    cache = Path(cache) if cache else None
    if cache and cache.exists():
        return json.loads(cache.read_text())
    index: dict[str, list[float]] = {}
    for name in list_tiles():
        try:
            index[name] = list(header_bbox(name))
        except Exception as exc:  # transient network / missing tile
            print(f"  skip {name}: {exc}")
    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(index))
    return index


def find_tile(easting: float, northing: float,
              cache: str | Path | None = None) -> str:
    """Return the tile whose bounding box contains a UTM 15N (E, N) point.

    Raises ``LookupError`` if no tile contains the point. If more than one bbox
    contains it (bboxes overlap slightly due to the grid skew), the tightest
    containing tile is returned.
    """
    index = build_index(cache)
    hits = [
        (name, bb) for name, bb in index.items()
        if bb[0] <= easting <= bb[2] and bb[1] <= northing <= bb[3]
    ]
    if not hits:
        raise LookupError(f"no tile contains ({easting}, {northing})")
    # tightest = smallest area, in case of overlapping bboxes at an edge
    name, _ = min(hits, key=lambda kv: (kv[1][2] - kv[1][0]) * (kv[1][3] - kv[1][1]))
    return name


def download_tile(name: str, out_dir: str | Path) -> Path:
    """Download a full tile LAZ to ``out_dir`` and return the local path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{name}.laz"
    urllib.request.urlretrieve(f"{LAZ_URL}{name}.laz", dest)
    return dest
