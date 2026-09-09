"""Flight-line ground tracks for the 2008 gen1 acquisition -- derived, and COMMITTED.

Why this module exists rather than a scratchpad script
-----------------------------------------------------
``analysis/investigations/gen1_same_line_not_reproduced/gen1_line_tracks.py`` derives the same thing but writes
``line_tracks.json`` to ``$SCRATCH``.  That file is gone, and with it the ability to
reproduce or extend the gen1 same-line datum -- see
``ground_control/GEN1_SAME_LINE_NOT_REPRODUCED.md``.  The tracks here are written into
the repository and committed, because the statewide goal needs this estimator at every
site and a number whose inputs evaporate cannot be checked or re-run.

The method is that module's, unchanged, and its four parameters are inherited from it
rather than re-chosen -- see :data:`INHERITED_PARAMS`.

Two things are needed from the flight lines and must not be confused
--------------------------------------------------------------------
* **Targeting** -- where does a line run, so which marks could it have seen?  That is what
  this module builds.
* **Assignment** -- which line actually illuminated a given mark?  That is NOT read off a
  track.  It comes from the ``point_source_id`` of the ground returns at the mark, via
  ``gen1_datum.assign_line_from_returns``.  A track fitted at one latitude walks off by
  hundreds of metres over tens of kilometres, so centreline distance mislabels marks
  exactly as the search widens.

``point_source_id`` is reused across missions
---------------------------------------------
In this acquisition one psid can span days (the source module measured psid 151 across
337,057 s, ~3.9 days).  A psid is therefore NOT a flight line: it is split into **passes**
at ``gps_time`` gaps longer than ``gap_s``, and each pass carries its own track.  Treating
a psid as one line would merge two passes flown days apart into a single "line" and
average their offsets -- which is the very structure the per-line datum exists to keep
apart.

A gap is not by itself proof of a second pass: it is equally a stretch of line under tiles
we do not hold.  Every pass therefore reports its gap length so the caller can judge.

Memory
------
This is a shared working laptop.  Tiles are streamed with ``laspy.chunk_iterator`` and
only the near-nadir strip is retained, so peak RSS stays far below reading a tile whole.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

#: Inherited verbatim from ``analysis/investigations/gen1_same_line_not_reproduced/gen1_line_tracks.py`` so the tracks are
#: the same object that module built.  They are that module's author's choices, not new
#: ones; each is passed explicitly and printed by the driver.
INHERITED_PARAMS = {
    "stride": 7,          # subsample stride; thins the near-nadir strip only
    "nadir_deg": 2.0,     # |scan_angle_rank| kept, degrees -- defines "near nadir"
    "bin_s": 1.0,         # gps_time bin width, seconds
    "gap_s": 120.0,       # split a psid into passes at gps_time gaps longer than this
    "min_bin_points": 20,  # a gps_time bin needs this many returns to yield a centroid
}


@dataclass(frozen=True)
class Pass:
    """One flight-line pass: a psid over one continuous stretch of mission time."""

    project: str              # the ACQUISITION this line belongs to, e.g. lidar_semn2008
    psid: int
    pass_index: int
    t0: float
    t1: float
    n_bins: int
    span_km: float
    heading_deg: float
    axis: str                 # "N" = fitted x on y (north-south run), "E" = the reverse
    vertices: list            # [(x, y), ...] along-track centroids, time-ordered
    resid_med_m: float        # median |perpendicular| residual of centroids to the fit
    resid_p95_m: float
    resid_max_m: float

    @staticmethod
    def psid_of_key(key: str) -> int:
        """The psid back out of a ``project:psid.pass_index`` key.

        Provided because the key format CHANGED on 2026-09-09 and the old
        ``key.split(".")[0]`` now returns ``project:psid`` -- which parses as neither an
        int nor an error in every caller. Use this rather than splitting by hand.
        """
        return int(str(key).split(":")[-1].split(".")[0])

    @property
    def key(self) -> str:
        """``project:psid.pass_index``. The PROJECT is part of the key because psid is a
        PER-PROJECT line number, not a global one: gen1 is four acquisitions across the six
        pilot sites, and cook uses psids 8-11 while battlecreek uses 1012-1014. Keyed by
        psid alone, a statewide TrackSet would eventually merge two different flight lines
        that happen to share a number -- silently, which is how this project's costliest
        errors have all presented."""
        return f"{self.project}:{self.psid}.{self.pass_index}"

    def endpoints(self):
        """Two endpoints of the fitted straight track, for a segment search."""
        return (tuple(self.vertices[0]), tuple(self.vertices[-1]))


@dataclass
class TrackSet:
    """Every pass found, plus how it was made."""

    passes: list
    params: dict
    tiles: list
    n_tiles_read: int
    n_returns_kept: int

    def by_line(self, project, psid):
        """Passes of one flight line. BOTH arguments are required: see ``Pass.key``."""
        return [p for p in self.passes
                if p.project == str(project) and p.psid == int(psid)]

    @property
    def projects(self) -> set:
        return {p.project for p in self.passes}

    @property
    def the_project(self) -> str:
        """The single acquisition this set covers, or raise.

        Callers that hold a bare psid need one. Raising when the set spans several is the
        point: a psid does not identify a flight line across acquisitions, so the caller
        must say which survey it means rather than get whichever matched first.
        """
        ps = self.projects
        if len(ps) != 1:
            raise ValueError(
                f"this TrackSet covers {sorted(ps)}; a bare psid does not identify a line "
                f"across acquisitions. Name the project explicitly.")
        return next(iter(ps))

    def as_search_tracks(self, psids=None, *, project=None) -> dict:
        """``{pass_key: [(x, y), ...]}`` in the shape ``discover_near_lines`` wants.

        Keyed by ``project:psid.pass_index``, never by psid alone, so neither two passes of
        one psid nor two projects' identically-numbered lines are merged. Filtering by
        ``psids`` REQUIRES ``project`` when the set spans more than one, because a bare
        psid does not identify a line across acquisitions.
        """
        if psids is not None and project is None and len(self.projects) > 1:
            raise ValueError(
                f"this TrackSet spans {sorted(self.projects)}, so a psid alone does not "
                f"identify a flight line. Pass project= as well.")
        want = None if psids is None else {int(s) for s in psids}
        out = {}
        for p in self.passes:
            if project is not None and p.project != str(project):
                continue
            if want is None or p.psid in want:
                out[p.key] = [tuple(v) for v in p.vertices]
        return out

    def to_dict(self):
        return dict(passes=[asdict(p) for p in self.passes], params=self.params,
                    tiles=self.tiles, n_tiles_read=self.n_tiles_read,
                    n_returns_kept=self.n_returns_kept)


def save_tracks(ts: TrackSet, path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(ts.to_dict(), indent=1))
    return p


def load_tracks(path) -> TrackSet:
    """Load a track set. REFUSES a file whose passes do not record their acquisition.

    Files written before 2026-09-09 carry no ``project``. Stamping one on load would be a
    guess, and guessing which survey a flight line belongs to is exactly the error the key
    exists to prevent -- so such a file must be migrated explicitly, with the evidence
    recorded. See scripts/migrate_line_tracks.py.
    """
    d = json.loads(Path(path).read_text())
    missing = [q for q in d["passes"] if "project" not in q]
    if missing:
        raise ValueError(
            f"{path}: {len(missing)} of {len(d['passes'])} passes record no `project`. "
            f"psid is a PER-PROJECT line number, so this file cannot be keyed safely. "
            f"Migrate it explicitly (scripts/migrate_line_tracks.py) rather than assuming "
            f"an acquisition.")
    return TrackSet(passes=[Pass(**q) for q in d["passes"]], params=d["params"],
                    tiles=d["tiles"], n_tiles_read=d["n_tiles_read"],
                    n_returns_kept=d["n_returns_kept"])


def _collect(tile_paths, *, stride, nadir_deg, chunk_size, progress=None):
    """Stream the tiles, keeping only the strided near-nadir strip."""
    import laspy

    acc: dict[int, list] = {}
    n_kept = 0
    used = []
    for i, path in enumerate(map(str, tile_paths)):
        with laspy.open(path) as fh:
            for chunk in fh.chunk_iterator(chunk_size):
                sa = np.abs(np.asarray(chunk.scan_angle_rank).astype(float))
                keep = sa <= nadir_deg
                if not keep.any():
                    continue
                idx = np.flatnonzero(keep)[::stride]
                if idx.size == 0:
                    continue
                x = np.asarray(chunk.x)[idx]
                y = np.asarray(chunk.y)[idx]
                t = np.asarray(chunk.gps_time)[idx]
                s = np.asarray(chunk.point_source_id)[idx].astype(int)
                for ln in np.unique(s):
                    m = s == ln
                    acc.setdefault(int(ln), []).append(
                        np.column_stack([t[m], x[m], y[m]]))
                n_kept += idx.size
        used.append(path)
        if progress:
            progress(i + 1, len(list(tile_paths)), path, n_kept)
    return {k: np.vstack(v) for k, v in acc.items()}, n_kept, used


def _fit_pass(t, x, y, *, bin_s, min_bin_points, psid, pass_index, project) -> Pass | None:
    """Bin one pass in gps_time, fit a straight track, return it."""
    k = np.round(t / bin_s).astype(np.int64)
    order = np.argsort(k, kind="stable")
    k, x, y, t = k[order], x[order], y[order], t[order]
    uk, idx, cnt = np.unique(k, return_index=True, return_counts=True)
    cx = np.add.reduceat(x, idx) / cnt
    cy = np.add.reduceat(y, idx) / cnt
    ct = np.add.reduceat(t, idx) / cnt
    good = cnt >= min_bin_points
    cx, cy, ct = cx[good], cy[good], ct[good]
    if ct.size < 2:
        return None
    # fit the longer axis on the shorter so the regression is well conditioned
    if np.ptp(cy) >= np.ptp(cx):
        w = np.polyfit(cy, cx, 1)
        r = np.abs(cx - np.polyval(w, cy)) / math.hypot(1.0, w[0])
        axis = "N"
    else:
        w = np.polyfit(cx, cy, 1)
        r = np.abs(cy - np.polyval(w, cx)) / math.hypot(1.0, w[0])
        axis = "E"
    span = math.hypot(np.ptp(cx), np.ptp(cy)) / 1000.0
    head = math.degrees(math.atan2(cx[-1] - cx[0], cy[-1] - cy[0])) % 360.0
    o = np.argsort(ct)
    return Pass(project=str(project), psid=int(psid), pass_index=int(pass_index),
                t0=float(ct.min()), t1=float(ct.max()), n_bins=int(ct.size),
                span_km=float(span), heading_deg=float(head), axis=axis,
                vertices=[[float(a), float(b)] for a, b in zip(cx[o], cy[o])],
                resid_med_m=float(np.median(r)), resid_p95_m=float(np.percentile(r, 95)),
                resid_max_m=float(r.max()))


def derive_tracks(tile_paths, *, project, stride, nadir_deg, bin_s, gap_s, min_bin_points,
                  chunk_size, progress=None) -> TrackSet:
    """Derive one track per flight-line PASS from the tiles on disk.

    Every parameter is required.  ``INHERITED_PARAMS`` holds the values
    ``gen1_line_tracks.py`` used; pass those to reproduce it.
    """
    tile_paths = [str(p) for p in tile_paths]
    raw, n_kept, used = _collect(tile_paths, stride=stride, nadir_deg=nadir_deg,
                                 chunk_size=chunk_size, progress=progress)
    passes = []
    for psid in sorted(raw):
        a = raw[psid]
        a = a[np.argsort(a[:, 0], kind="stable")]
        t, x, y = a[:, 0], a[:, 1], a[:, 2]
        # split into passes at mission-time gaps: a psid is not a flight line
        cut = np.flatnonzero(np.diff(t) > gap_s) + 1
        for j, (lo, hi) in enumerate(zip(np.r_[0, cut], np.r_[cut, t.size])):
            p = _fit_pass(t[lo:hi], x[lo:hi], y[lo:hi], bin_s=bin_s,
                          min_bin_points=min_bin_points, psid=psid, pass_index=j,
                          project=project)
            if p is not None:
                passes.append(p)
    return TrackSet(passes=passes, tiles=used, n_tiles_read=len(used),
                    n_returns_kept=int(n_kept),
                    params=dict(stride=stride, nadir_deg=nadir_deg, bin_s=bin_s,
                                gap_s=gap_s, min_bin_points=min_bin_points,
                                chunk_size=chunk_size))
