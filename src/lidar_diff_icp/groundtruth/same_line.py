"""Attribute a control mark to a physical FLIGHT LINE, from the returns themselves.

Promoted from ``ground_control/same_line.py`` on 2026-09-09, per its INTEGRATION.md: the
site-local estimator, which the statewide datum needs everywhere rather than at Elba only.

What came across is the RETURNS route and the collinearity that supports it. What stayed
behind is the CATCHMENT route -- ``Scope``, ``site_scope``, ``discover``,
``SEAM_HALF_SPACING_M`` and the catchment ``estimate`` -- which is superseded for measuring
a datum, because it conflates "found near a track" with "belongs to that line" and cost
~9 mm at Elba. It is not deleted: five drivers still call it, and three of those exist to
COMPARE the two routes, so they need the superseded one by definition.

INTEGRATION.md also asked for ``estimate_by_returns`` to be renamed ``estimate`` here. It
is NOT, deliberately: while the catchment ``estimate`` survives in ``ground_control``, two
functions named ``estimate`` with opposite meanings is the same trap as the project's two
``datum.py``. Rename it when the catchment route retires, not before.

A psid is not a flight line -- see :mod:`lidar_diff_icp.groundtruth.lines`. Passes of one
psid are merged into physical lines by COLLINEARITY, never by ``gps_time``: measured on
Elba's psids the correlation between the time gap and the collinearity sigma is -0.32,
the wrong sign.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lidar_diff_icp.groundtruth import lines as L      # noqa: F401  (used in type text)
from lidar_diff_icp.groundtruth import gen1_datum as G

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
    ps = list(trackset.by_line(trackset.the_project, psid))
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


def _track_distance(px, py, vertices) -> float:
    v = np.asarray(vertices, float)
    a, b = v[:-1], v[1:]
    ab = b - a
    denom = np.maximum((ab ** 2).sum(1), 1e-9)
    t = np.clip(((px - a[:, 0]) * ab[:, 0] + (py - a[:, 1]) * ab[:, 1]) / denom, 0.0, 1.0)
    c = a + t[:, None] * ab
    return float(np.hypot(px - c[:, 0], py - c[:, 1]).min())


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
    cands = trackset.by_line(trackset.the_project, int(psid))
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
