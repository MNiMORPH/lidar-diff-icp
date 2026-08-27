"""gen1's datum at a site from marks on the SITE'S OWN flight lines.

Why this is a different estimator from the residual field
---------------------------------------------------------
gen1's residual is organised per flight line (``analysis/FRAME_2026-08-26-PM.md``:
one-way ANOVA F = 8.63, p < 0.001, line means spanning -146.6 to +113.5 mm).  A kriged
field over all marks carries no line term, so it pools across lines and pulls a site's
answer toward the regional mean.  This module does the opposite: it uses only marks the
site's OWN lines illuminated, and takes the LINE as the unit of replication.

The pass, not the point_source_id, is the line
----------------------------------------------
``point_source_id`` is reused across this acquisition.  Measured on Elba's own psids from
``ground_control/data/gen1_line_tracks.json``: two passes sharing a psid sit **10.4 to
83.3 km apart along track** (median 52.1) and **67 to 1059 m apart across track** -- the
latter being about one 961 m line spacing.  They are different physical lines, not one
line with a hole where tiles are missing.

**Where that distinction acts is the SEARCH, not the aggregation.** Only one pass per psid
comes near any given site (at Elba the next-nearest pass of the same psid is tens of km
away), so once the search is restricted to the site's own passes, the psid read off the
returns identifies the pass uniquely and ``gen1_datum.combine_datum`` can be used
UNCHANGED.  Searching by psid instead reaches into another survey block and labels marks
tens of kilometres away as the site's line.

:func:`estimate` therefore takes ``scope="pass"`` or ``scope="psid"`` and the two differ in
**one** thing: which tracks the search walks.  The tile resolution, the ground source, the
tie estimator and the combination are byte-identical between them, so a difference between
the two answers is the search and nothing else.

Assignment still comes from the returns
---------------------------------------
A track only TARGETS a mark.  Which line illuminated it is settled inside
``gen1_datum.measure_site`` by ``assign_line_from_returns``.  This module never assigns by
distance to a track.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))

import lines as L  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G  # noqa: E402

SCOPES = ("pass", "psid")

#: Half of the measured gen1 line spacing (942/987/932/988/956 m, mean 961) -- the
#: vendor's own class-12 overlap seam, at which bare-earth ground is cut.  From
#: ``analysis/GEN1_DATUM_MORE_MARKS.md`` section 1.  A MEASUREMENT of this acquisition,
#: not a chosen search radius: past it, ground belongs to the neighbouring line.
SEAM_HALF_SPACING_M = 481.0


@dataclass(frozen=True)
class Scope:
    """Which tracks the search walked, and why those."""

    scope: str
    psids: tuple
    track_keys: tuple
    n_tracks: int
    dropped_track_keys: tuple
    note: str


def site_scope(trackset: L.TrackSet, *, psids, easting, northing, scope: str) -> Scope:
    """The tracks to search along, for one site.

    ``scope="pass"``  -- for each psid, ONLY the pass whose track runs nearest the site.
    ``scope="psid"``  -- every pass of every psid, which is what code keying on
                         point_source_id alone effectively does.
    """
    if scope not in SCOPES:
        raise ValueError(f"scope={scope!r} must be one of {SCOPES}")
    psids = tuple(int(p) for p in psids)
    keep, dropped = [], []
    for p in psids:
        cands = trackset.by_psid(p)
        if not cands:
            continue
        if scope == "psid":
            keep.extend(q.key for q in cands)
            continue
        best = min(cands, key=lambda q: _track_distance(easting, northing, q.vertices))
        keep.append(best.key)
        dropped.extend(q.key for q in cands if q.key != best.key)
    note = ("only the pass of each psid running nearest the site; the others are "
            "different physical lines that reuse the id"
            if scope == "pass" else
            "EVERY pass of each psid, including ones in other survey blocks -- the "
            "behaviour of code that keys on point_source_id alone")
    return Scope(scope=scope, psids=psids, track_keys=tuple(keep), n_tracks=len(keep),
                 dropped_track_keys=tuple(dropped), note=note)


def _track_distance(px, py, vertices) -> float:
    v = np.asarray(vertices, float)
    a, b = v[:-1], v[1:]
    ab = b - a
    denom = np.maximum((ab ** 2).sum(1), 1e-9)
    t = np.clip(((px - a[:, 0]) * ab[:, 0] + (py - a[:, 1]) * ab[:, 1]) / denom, 0.0, 1.0)
    c = a + t[:, None] * ab
    return float(np.hypot(px - c[:, 0], py - c[:, 1]).min())


def discover(control, trackset: L.TrackSet, sc: Scope, *, half_width_m: float):
    """Marks within ``half_width_m`` of any track in the scope.  A SEARCH, not an assignment."""
    tracks = {k: v for k, v in trackset.as_search_tracks().items() if k in set(sc.track_keys)}
    return G.discover_near_lines(control, tracks, half_width_m)


def estimate(trackset, *, psids, easting, northing, scope, half_width_m, covers,
             tile_dirs, res, on_missing="skip", control=None):
    """Discover -> resolve -> measure -> combine, for one scope.

    Returns ``(Scope, sites, measurements, skipped, Gen1DatumEstimate)``.  Everything
    downstream of :func:`discover` is ``gen1_datum``'s own code, unchanged.
    """
    control = G.load_control() if control is None else control
    sc = site_scope(trackset, psids=psids, easting=easting, northing=northing, scope=scope)
    sites = discover(control, trackset, sc, half_width_m=half_width_m)
    if covers is not None:
        keep = set(covers)
        sites = [s for s in sites if s.mark.cover_class in keep]
    resolution = G.resolve_tiles(sites, tile_dirs)
    meas, skipped = G.measure_sites(sites, resolution, on_missing=on_missing, res=res)
    est = G.combine_datum(meas, mode="per_line",
                          notes=(f"search scope: {sc.scope} -- {sc.note}",
                                 f"tracks searched: {', '.join(sc.track_keys)}"))
    return sc, sites, meas, skipped, est


def marks_on_scope_psids(meas, psids):
    """Measurements whose RETURNS put them on one of the site's psids."""
    want = {int(p) for p in psids}
    return [m for m in meas if m.line_id is not None and int(m.line_id) in want]
