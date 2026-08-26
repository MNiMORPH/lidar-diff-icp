"""gen1's absolute vertical datum at ANY Minnesota site, from its OWN 2008 control.

What this is for
----------------
The 2008 MN DNR lidar (gen1) floats in z. The cheapest, cleanest reference to pin it
against is the control the 2008 acquisition was itself validated on: 1 004 rows of
surveyed checkpoints transcribed from the eight MnGeo county validation reports and
bundled as ``data/mn_dnr_2008_control_semn.csv``. That control is on **the same vertical
datum and the same geoid model as the raw gen1 cloud** (NAVD88, GEOID03), so a
comparison against it carries **no geoid conversion and no cross-epoch term at all**.

``analysis/GEN1_OWN_CONTROL_TIE.md`` did this once, by hand, at Elba. This module is the
same measurement made reusable and pointed anywhere: discover the marks near a site (or
near a set of flight lines), say which tiles are needed and which are already on disk,
read each mark with the committed estimator in :mod:`~lidar_diff_icp.groundtruth.tie`,
and combine.

The four things it refuses to do
--------------------------------
1. **It never converts a geoid.** :func:`assert_no_geoid_conversion` compares the marks'
   recorded geoid model against the lidar's and RAISES on a mismatch. The 2008 case is
   an equality, not a cancellation that has to be argued -- but it is asserted in code
   on every run rather than assumed in a comment.
2. **It never brings in a gen2-derived term.** The pipeline's Nuth & Kaeaeb lateral
   shift is measured *between epochs*; it has no place in gen1-against-its-own-control.
   ``lateral_shift_m`` therefore defaults to ``None`` and must be passed deliberately.
3. **It never downloads.** :func:`resolve_tiles` reports which tiles hold the marks and
   which of those are on disk. Fetching is the caller's act, one tile at a time.
4. **It invents no cut.** There is no minimum ``n``, no maximum slope, no spread
   threshold, no distance cut with a default. Every screen statistic is *returned*;
   which marks to keep is the caller's decision and it is recorded when they make it.

Flight lines come from the RETURNS
----------------------------------
Each mark is assigned to a flight line by the ``point_source_id`` of the ground returns
at the mark (:func:`assign_line_from_returns`), not by distance to a fitted centreline.
Line spacing in this acquisition is ~1 km and the nadir tracks were fitted at one
latitude, so a centreline assignment mislabels marks as the search widens -- and the
scatter this module has to model is organised *by line*, so a mislabel goes straight
into the error bar. The returns cannot be wrong about which sortie hit the ground.

The error bar, and what it is the error of
------------------------------------------
Per-mark ties are **not independent**: marks under one flight line share that swath's
unknown constant. A one-way ANOVA over the line groups is computed on every run
(:attr:`Gen1DatumEstimate.anova_F`) and the flight line, not the mark, is the unit of
replication:

    value = mean over lines of ( mean over that line's marks )
    SE    = sd(line means) / sqrt(number of lines)

so ``se_of`` reads *"SE of the mean over flight lines of the within-line mean tie"* --
which is what :attr:`Gen1DatumEstimate.se_of` literally contains, printed with the
number. The per-mark ("independence assumed") SE is returned beside it, together with
the design effect, because the ratio is the size of the mistake that assumption makes.

Two combination modes
---------------------
``mode="per_line"``
    Each line's swath constant is unknown and is treated as an independent draw. This is
    the mode for a fresh site where no ``corrections.json`` exists yet.
``mode="common_datum"``
    The per-swath constants solved by :func:`lidar_diff_icp.coreg.align_swaths` (stored
    in a tile's ``corrections.json`` / ``corrections_geoid.json``) are applied to the
    returns before the estimate, putting every mark in ONE frame. The scalar is then the
    datum of that frame -- and the **per-line residuals are returned**, because they are
    an out-of-sample test of an overlap-derived swath network against ground truth,
    which the network has no internal redundancy to provide for itself.

Sign convention, unchanged from :mod:`~lidar_diff_icp.groundtruth.tie`: a tie is
``surveyed - z_lidar``, i.e. the constant to **ADD to gen1**. Positive means gen1 reads
low.
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .checkpoints import Checkpoint, read_checkpoint_csv
from .provenance import Param
from .tie import GroundReturns, TieEstimate, estimate_tie, radius_ladder, vendor_ground_near

_DATA = Path(__file__).with_name("data")

#: The vertical datum and geoid model of the delivered 2008 SE-MN lidar. Asserted, never
#: converted. Source: the dataset metadata page ``lidar_semn2008.html``, "Vertical datum:
#: NAVD88 (Geoid03)"; recorded per row in the bundled control CSV's ``verified`` column,
#: which also states that the validation reports themselves print no datum.
GEN1_VERTICAL_DATUM = "NAVD88"
GEN1_GEOID_MODEL = "GEOID03"
_SRC_GEN1_DATUM = ("lidar_semn2008.html 'Vertical datum: NAVD88 (Geoid03)'; the bundled "
                   "control CSV's verified column carries the same statement per row and "
                   "flags it as a DATASET-level assertion")

#: The bundled 2008 control set.
DEFAULT_CONTROL = "mn_dnr_2008_control_semn"


class DatumMismatchError(ValueError):
    """The control and the lidar are not on the same vertical datum / geoid model.

    Raised instead of converting. This module exists because the 2008 control needs no
    conversion; the moment it would need one, the caller is doing a different experiment
    and must say so with :func:`~lidar_diff_icp.groundtruth.tie.geoid_shift_for`.
    """


# ----------------------------------------------------------------- the control set

@dataclass(frozen=True)
class ControlMark:
    """One physical surveyed mark, with every report that publishes it.

    A mark on a county line is printed in **both** counties' validation reports, so the
    transcription has more rows than marks. The rows are merged here on exact equality of
    ``(easting, northing, elevation)`` -- no tolerance, no rounding -- and every id and
    report that contributed is kept in ``aliases`` and ``reports``, so nothing about
    where the mark came from is lost by the merge.
    """

    checkpoint: Checkpoint
    cover_class: str                 # L1O open / L2T weeds+crops / L3B brush / L4F forest / L5U urban
    counties: tuple                  # every county report the mark appears in
    aliases: tuple                   # every point_id spelling seen for it
    reports: tuple                   # every source string seen for it
    dnr_surface_z_m: float | None    # the report's own delivered-surface elevation
    dnr_error_m: float | None        # the report's own Control Z - Surface Z

    @property
    def point_id(self) -> str:
        return self.checkpoint.point_id

    @property
    def easting(self) -> float:
        return self.checkpoint.easting

    @property
    def northing(self) -> float:
        return self.checkpoint.northing


@dataclass
class ControlSet:
    """Every mark in a bundled control transcription, with the merge reported.

    ``n_rows`` is what the CSV holds; ``len(self)`` is how many physical marks that is.
    The difference is listed row by row in ``merges`` -- a duplicate is never dropped
    silently, and the count that changed is printed with the reason it changed.
    """

    marks: list
    origin: str
    n_rows: int
    merges: list = field(default_factory=list)   # [(kept_id, [merged_ids], [counties])]

    def __len__(self) -> int:
        return len(self.marks)

    def __iter__(self):
        return iter(self.marks)

    def __getitem__(self, key):
        if isinstance(key, str):
            for m in self.marks:
                if m.point_id == key or key in m.aliases:
                    return m
            raise KeyError(key)
        return self.marks[key]

    @property
    def merge_note(self) -> str:
        return (f"{self.n_rows} transcribed rows -> {len(self.marks)} physical marks; "
                f"{self.n_rows - len(self.marks)} rows merged as the same mark published "
                f"in more than one county report (exact match on easting, northing and "
                f"elevation)")

    def xy(self) -> np.ndarray:
        return np.array([[m.easting, m.northing] for m in self.marks], float)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_control(name: str = DEFAULT_CONTROL) -> ControlSet:
    """Load a bundled control transcription, merging rows that are the same mark.

    Two rows are the same mark when their easting, northing **and** elevation are
    exactly equal as read from the file. That is an identity test, not a tolerance: two
    distinct monuments cannot share a coordinate to the millimetre. Rows that share a
    ``point_id`` but sit at different positions raise, rather than one silently winning.
    """
    path = _DATA / f"{name}.csv"
    if not path.exists():
        raise FileNotFoundError(f"no bundled control set {name!r} in {_DATA}")
    cset = read_checkpoint_csv(path)
    with open(path, newline="") as fh:
        raw = list(csv.DictReader(fh))
    if len(raw) != len(cset.points):
        raise AssertionError(f"{path}: {len(raw)} raw rows vs {len(cset.points)} parsed")

    seen_id, order, groups = {}, [], {}
    for cp, row in zip(cset.points, raw):
        key = (cp.easting, cp.northing, cp.elevation)
        prev = seen_id.setdefault(cp.point_id, key)
        if prev != key:
            raise ValueError(
                f"{path}: point_id {cp.point_id!r} appears at two different positions "
                f"{prev} and {key}. Two marks cannot share an id; fix the transcription "
                "rather than letting one of them win silently.")
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((cp, row))

    marks, merges = [], []
    for key in order:
        members = groups[key]
        head = members[0][0]
        ids = tuple(dict.fromkeys(cp.point_id for cp, _ in members))
        counties = tuple(dict.fromkeys(r.get("county", "") for _, r in members if r.get("county")))
        reports = tuple(dict.fromkeys(r.get("source", "") for _, r in members if r.get("source")))
        first = members[0][1]
        marks.append(ControlMark(
            checkpoint=head, cover_class=head.point_type, counties=counties,
            aliases=ids, reports=reports,
            dnr_surface_z_m=_f(first.get("dnr_surface_z_m")),
            dnr_error_m=_f(first.get("dnr_error_m"))))
        if len(members) > 1:
            merges.append((head.point_id, list(ids[1:]), list(counties)))
    return ControlSet(marks=marks, origin=str(path.resolve()),
                      n_rows=len(raw), merges=merges)


# --------------------------------------------------------- the no-conversion assertion

def _as_checkpoint(obj) -> Checkpoint:
    """Unwrap a MarkSite / ControlMark / Checkpoint down to the Checkpoint."""
    for _ in range(4):
        if isinstance(obj, Checkpoint):
            return obj
        nxt = getattr(obj, "checkpoint", None)
        if nxt is None:
            break
        obj = nxt
    raise TypeError(f"expected a Checkpoint, ControlMark or MarkSite, got "
                    f"{type(obj).__name__}")


def assert_no_geoid_conversion(marks, *, lidar_vertical_datum: str = GEN1_VERTICAL_DATUM,
                               lidar_geoid_model: str = GEN1_GEOID_MODEL) -> str:
    """Assert control and lidar share a vertical datum and geoid model; return the claim.

    The returned string is the sentence a run should print: it names both sides and the
    number of marks checked, so "no geoid conversion" appears in the output as a checked
    statement rather than as an author's assurance. Raises :class:`DatumMismatchError`
    naming the offending marks otherwise.
    """
    marks = list(marks)
    if not marks:
        raise ValueError("assert_no_geoid_conversion needs at least one mark")
    want = (lidar_vertical_datum.strip().upper(),
            lidar_geoid_model.strip().upper().replace(" ", ""))
    bad = []
    for m in marks:
        cp = _as_checkpoint(m)
        cp.require_datum()
        got = (cp.vertical_datum.strip().upper(),
               cp.geoid_model.strip().upper().replace(" ", ""))
        if got != want:
            bad.append((cp.point_id, got))
    if bad:
        raise DatumMismatchError(
            f"{len(bad)} mark(s) are not on the lidar's datum "
            f"{lidar_vertical_datum}/{lidar_geoid_model}: "
            + ", ".join(f"{p} on {a}/{b}" for p, (a, b) in bad[:5])
            + (" ..." if len(bad) > 5 else "")
            + ". This module does not convert geoids -- a conversion is a different "
              "experiment and must be requested explicitly with tie.geoid_shift_for.")
    return (f"no geoid conversion: all {len(marks)} marks and the lidar are on "
            f"{lidar_vertical_datum}({lidar_geoid_model}); checked, not assumed")


# ------------------------------------------------------------------------- discovery

@dataclass(frozen=True)
class MarkSite:
    """A control mark placed relative to whatever the search was about."""

    mark: ControlMark
    distance_m: float
    relative_to: str            # "site (E, N)" or "flight line <id>"
    nearest_feature: str = ""   # which line, when the search was over lines

    @property
    def checkpoint(self) -> Checkpoint:
        return self.mark.checkpoint

    @property
    def point_id(self) -> str:
        return self.mark.point_id


def discover_near_point(control: ControlSet, easting: float, northing: float,
                        radius_m: float) -> list:
    """Marks within ``radius_m`` of a site, nearest first.

    ``radius_m`` has no default: how far to reach for control is a decision about the
    site (how far a datum may be assumed constant), not a property of the code.
    """
    radius_m = float(radius_m)
    out = []
    for m in control:
        d = math.hypot(m.easting - easting, m.northing - northing)
        if d <= radius_m:
            out.append(MarkSite(m, d, f"site ({easting:.1f}, {northing:.1f})"))
    return sorted(out, key=lambda s: s.distance_m)


def _point_segment_distance(px, py, x0, y0, x1, y1):
    vx, vy = x1 - x0, y1 - y0
    L2 = vx * vx + vy * vy
    if L2 == 0.0:
        return math.hypot(px - x0, py - y0)
    t = max(0.0, min(1.0, ((px - x0) * vx + (py - y0) * vy) / L2))
    return math.hypot(px - (x0 + t * vx), py - (y0 + t * vy))


def discover_near_lines(control: ControlSet, lines: dict, half_width_m: float) -> list:
    """Marks within ``half_width_m`` of any of a set of flight-line tracks.

    ``lines`` maps a line id to a track: either two endpoints ``((x0, y0), (x1, y1))`` or
    a sequence of vertices. ``half_width_m`` has no default -- it is the swath half-width
    of the acquisition being studied, which the caller measures.

    **This is a search, not an assignment.** Which line actually hit a mark is settled by
    :func:`assign_line_from_returns` from the returns themselves; a track fitted at one
    latitude walks off by hundreds of metres over tens of kilometres, which is precisely
    the mislabelling this module's line statistics must not inherit.
    """
    half_width_m = float(half_width_m)
    out = []
    for m in control:
        best, best_line = math.inf, None
        for lid, track in lines.items():
            pts = list(track)
            for a, b in zip(pts[:-1], pts[1:]):
                d = _point_segment_distance(m.easting, m.northing, a[0], a[1], b[0], b[1])
                if d < best:
                    best, best_line = d, lid
        if best <= half_width_m:
            out.append(MarkSite(m, best, f"nearest of {len(lines)} flight-line tracks",
                                nearest_feature=str(best_line)))
    return sorted(out, key=lambda s: s.distance_m)


# -------------------------------------------------------------------- tile resolution

@dataclass(frozen=True)
class TileNeed:
    """One tile, the marks it holds, and whether it is already on disk."""

    tile: str
    path: str | None
    marks: tuple

    @property
    def on_disk(self) -> bool:
        return self.path is not None


@dataclass
class TileResolution:
    """Which tile each mark falls in, split into on-disk and to-be-fetched.

    Nothing here fetches anything. ``to_fetch`` is a list for a human to work through one
    tile at a time; :func:`lidar_diff_icp.tiles.download_tile` is the caller's to run.
    """

    per_mark: dict                    # point_id -> tile name
    needs: list                       # [TileNeed]
    search_dirs: tuple
    index_cache: str

    @property
    def on_disk(self) -> list:
        return [n for n in self.needs if n.on_disk]

    @property
    def to_fetch(self) -> list:
        return [n for n in self.needs if not n.on_disk]

    def path_for(self, point_id: str) -> str | None:
        t = self.per_mark.get(point_id)
        for n in self.needs:
            if n.tile == t:
                return n.path
        return None

    @staticmethod
    def table_columns() -> dict:
        return {
            "tile": "MnGeo AAAA-BB-CC tile name from tiles.find_tile (statewide centroid index)",
            "on_disk": "True when a file of that name was found in one of search_dirs",
            "n_marks": "control marks falling in that tile, count",
            "marks": "their point ids",
            "path": "the file found, or '' when the tile must be fetched",
        }

    def table_rows(self) -> list:
        return [[n.tile, n.on_disk, len(n.marks), ", ".join(n.marks[:4])
                 + (" ..." if len(n.marks) > 4 else ""), n.path or ""]
                for n in self.needs]


def resolve_tiles(sites, search_dirs, *, cache=None, suffixes=(".laz", ".las")
                  ) -> TileResolution:
    """Name the tile each mark falls in and say which of those are already on disk.

    ``search_dirs``  directories to look in for ``<tile><suffix>``. No download is ever
                     attempted -- a missing tile is reported, not fetched, because on a
                     shared machine the fetch is a decision with a cost.
    ``cache``        the statewide centroid index CSV passed to
                     :func:`lidar_diff_icp.tiles.find_tile`; ``None`` uses that
                     function's own repo default. If the cache does not exist the lookup
                     will try to build it, which needs the network -- pass a built cache
                     to keep this offline.
    """
    from ..tiles import find_tile

    dirs = tuple(str(d) for d in (search_dirs if not isinstance(search_dirs, (str, Path))
                                  else [search_dirs]))
    kw = {} if cache is None else {"cache": cache}
    per_mark, by_tile = {}, {}
    for s in sites:
        t = find_tile(s.mark.easting, s.mark.northing, **kw)
        per_mark[s.point_id] = t
        by_tile.setdefault(t, []).append(s.point_id)

    needs = []
    for t, ids in by_tile.items():
        found = None
        for d in dirs:
            for suf in suffixes:
                p = os.path.join(d, f"{t}{suf}")
                if os.path.exists(p):
                    found = p
                    break
            if found:
                break
        needs.append(TileNeed(tile=t, path=found, marks=tuple(ids)))
    needs.sort(key=lambda n: (not n.on_disk, n.tile))
    from ..tiles import DEFAULT_TILE_INDEX_CACHE
    return TileResolution(per_mark=per_mark, needs=needs, search_dirs=dirs,
                          index_cache=str(cache or DEFAULT_TILE_INDEX_CACHE))
