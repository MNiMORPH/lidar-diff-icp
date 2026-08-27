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

SCOPES = ("pass", "psid", "track")

#: Half of the measured gen1 line spacing (942/987/932/988/956 m, mean 961) -- the
#: vendor's own class-12 overlap seam, at which bare-earth ground is cut.  From
#: ``analysis/GEN1_DATUM_MORE_MARKS.md`` section 1.  A MEASUREMENT of this acquisition,
#: not a chosen search radius: past it, ground belongs to the neighbouring line.
SEAM_HALF_SPACING_M = 481.0


def collinearity_sigma(A, B):
    """How many prediction-sd away pass B sits from pass A's extrapolated track.

    A near-N-S line at heading 179.3 deg drifts ~1.1 km in easting over 94 km of track,
    so a raw easting separation of ~800 m between two passes 56 km apart is what ONE
    continuous line looks like -- it is NOT evidence of two lines.  This extrapolates the
    longer pass's fitted track to the shorter one's position and scales the miss by the
    extrapolation's own prediction sd, which is the only way the comparison has meaning
    tens of kilometres beyond the data.
    """
    if len(A.vertices) < len(B.vertices):
        A, B = B, A
    va = np.asarray(A.vertices, float); vb = np.asarray(B.vertices, float)
    xn, ye = va[:, 1], va[:, 0]
    if np.ptp(xn) < 100 or xn.size < 3:
        return float("inf")
    w = np.polyfit(xn, ye, 1)
    r = ye - np.polyval(w, xn)
    s = float(np.sqrt((r ** 2).sum() / max(xn.size - 2, 1)))
    nb = float(np.median(vb[:, 1]))
    Sxx = float(((xn - xn.mean()) ** 2).sum())
    pred_sd = s * np.sqrt(1.0 / xn.size + (nb - xn.mean()) ** 2 / Sxx)
    if pred_sd <= 0:
        return float("inf")
    return abs(np.polyval(w, nb) - np.median(vb[:, 0])) / pred_sd


def collinear_groups(trackset, psid, *, sigma):
    """Passes of one psid, merged into physical LINES by collinearity.

    ``sigma`` is the caller's: it is how far, in units of the extrapolation's own
    uncertainty, two passes may sit apart and still be called one line.  There is no
    default -- at Elba's six psids the verdicts run 0.1 to 21.6 sigma, so where the line
    is drawn changes the grouping.
    """
    ps = list(trackset.by_psid(psid))
    parent = {p.key: p.key for p in ps}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k

    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            if collinearity_sigma(ps[i], ps[j]) < sigma:
                a, b = find(ps[i].key), find(ps[j].key)
                if a != b:
                    parent[a] = b
    out = {}
    for p in ps:
        out.setdefault(find(p.key), []).append(p.key)
    return list(out.values())


@dataclass(frozen=True)
class Scope:
    """Which tracks the search walked, and why those."""

    scope: str
    psids: tuple
    track_keys: tuple
    n_tracks: int
    dropped_track_keys: tuple
    note: str


def site_scope(trackset: L.TrackSet, *, psids, easting, northing, scope: str,
               collinear_sigma: float | None = None) -> Scope:
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
        if scope == "pass":
            keep.append(best.key)
            dropped.extend(q.key for q in cands if q.key != best.key)
            continue
        if collinear_sigma is None:
            raise ValueError("scope='track' requires collinear_sigma; where the line is "
                             "drawn changes the grouping and has no default")
        grp = next(g for g in collinear_groups(trackset, p, sigma=collinear_sigma)
                   if best.key in g)
        keep.extend(grp)
        dropped.extend(q.key for q in cands if q.key not in grp)
    note = ("only the pass of each psid running nearest the site -- the most "
            "conservative choice, and it DISCARDS same-line marks whenever a gap was a "
            "tile hole rather than a separate flight"
            if scope == "pass" else
            "every pass COLLINEAR with the site's nearest one, i.e. the physical flight "
            "line, merged across gaps that are tile holes"
            if scope == "track" else
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
             tile_dirs, res, on_missing="skip", control=None, collinear_sigma=None):
    """Discover -> resolve -> measure -> combine, for one scope.

    Returns ``(Scope, sites, measurements, skipped, Gen1DatumEstimate)``.  Everything
    downstream of :func:`discover` is ``gen1_datum``'s own code, unchanged.
    """
    control = G.load_control() if control is None else control
    sc = site_scope(trackset, psids=psids, easting=easting, northing=northing,
                    scope=scope, collinear_sigma=collinear_sigma)
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


# ------------------------------------- the catchment-free estimator (preferred)

def marks_in_tiles(control, tile_dirs, *, covers=None):
    """Every control mark that falls inside a gen1 tile ON DISK.

    This replaces the catchment search entirely.  The catchment was only ever a compute
    bound on candidates -- ``assign_line_from_returns`` does the assigning, and it can
    only reject -- so bounding by "is the tile here" is both simpler and strictly more
    complete than bounding by distance to a fitted track.  It also removes the confound
    a radius introduces: at Elba, widening 481 m -> 2000 m added four marks that were ALL
    urban, so the radius was silently shifting the COVER mix.
    """
    import glob
    import laspy

    paths = []
    for d in tile_dirs:
        paths.extend(p for p in glob.glob(f"{d}/*.laz") if "merged" not in p)
    boxes = {}
    for p in sorted(paths):
        with laspy.open(p) as f:
            h = f.header
            boxes[p] = (h.mins[0], h.mins[1], h.maxs[0], h.maxs[1])
    out = []
    for m in control:
        if covers is not None and m.cover_class not in set(covers):
            continue
        e, n = m.checkpoint.easting, m.checkpoint.northing
        for p, b in boxes.items():
            if b[0] <= e <= b[2] and b[1] <= n <= b[3]:
                out.append((G.MarkSite(m, 0.0, f"inside tile {Path(p).name}"), p))
                break
    return out


def on_site_line(trackset, mark_easting, mark_northing, psid, *, site_easting,
                 site_northing, collinear_sigma):
    """Is this mark on the SAME PHYSICAL LINE as the site's pass of ``psid``?

    Only matters for a psid carrying more than one pass -- 16 of 41 here.  The mark is
    attributed to the pass of that psid whose track runs nearest IT, and that pass is
    then tested for collinearity with the pass nearest the SITE.

    ``gps_time`` cannot do this job and must not be substituted: measured on Elba's
    psids, the correlation between the gps_time gap and the collinearity sigma is
    **-0.32** -- the wrong sign.  Adjacent lines flown back-to-back are close in time and
    different (138.0/138.1: 135 s apart, 21.6 sigma), while one line interrupted by
    missing tiles has a long gap and is the same (133.0/133.1: 682 s apart, 0.2 sigma).
    """
    cands = trackset.by_psid(int(psid))
    if len(cands) <= 1:
        return True, "single pass: no ambiguity"
    site_pass = min(cands, key=lambda q: _track_distance(site_easting, site_northing,
                                                         q.vertices))
    mark_pass = min(cands, key=lambda q: _track_distance(mark_easting, mark_northing,
                                                         q.vertices))
    if mark_pass.key == site_pass.key:
        return True, f"same pass {mark_pass.key}"
    sg = collinearity_sigma(site_pass, mark_pass)
    ok = sg < collinear_sigma
    return ok, f"{mark_pass.key} vs site {site_pass.key}: {sg:.1f} sigma"


def estimate_by_returns(trackset, *, psids, easting, northing, covers, tile_dirs, res,
                        collinear_sigma, control=None, on_missing="skip"):
    """The catchment-free estimate: every mark in a tile, assigned by its own returns.

    Returns ``(measurements, kept, rejected, Gen1DatumEstimate)``.
    """
    control = G.load_control() if control is None else control
    cand = marks_in_tiles(control, tile_dirs, covers=covers)
    sites = [s for s, _ in cand]
    resolution = G.resolve_tiles(sites, tile_dirs)
    meas, skipped = G.measure_sites(sites, resolution, on_missing=on_missing, res=res)
    want = {int(p) for p in psids}
    kept, rejected = [], []
    for m in meas:
        if m.line_id is None:
            rejected.append((m.point_id, None, "no ground returns: no line"))
            continue
        if int(m.line_id) not in want:
            rejected.append((m.point_id, m.line_id, "returns place it on another line"))
            continue
        ok, why = on_site_line(trackset, m.site.mark.easting, m.site.mark.northing,
                               m.line_id, site_easting=easting, site_northing=northing,
                               collinear_sigma=collinear_sigma)
        (kept if ok else rejected).append(
            (m.point_id, m.line_id, why) if not ok else m)
    est = G.combine_datum([k for k in kept], mode="per_line",
                          notes=("candidates: every control mark inside a tile on disk; "
                                 "NO catchment radius",
                                 "assignment: assign_line_from_returns",
                                 f"reused psids disambiguated by collinearity at "
                                 f"{collinear_sigma} sigma"))
    return meas, kept, rejected, est
