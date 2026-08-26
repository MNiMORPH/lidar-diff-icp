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


# ------------------------------------------------------- flight lines FROM THE RETURNS

@dataclass(frozen=True)
class LineAssignment:
    """Which flight line(s) hit the ground at a mark, counted from ``point_source_id``.

    ``dominant`` is the line contributing the most ground returns inside ``radius_m``.
    ``mixed`` is True whenever more than one line is present -- a fact, not a rejection:
    a mixed mark's tie folds the two lines' relative offset into itself, and the caller
    decides what to do about that with the per-line ties in
    :attr:`MarkMeasurement.per_line_tie_mm`.
    """

    counts: dict                 # line id -> ground returns inside radius_m
    dominant: int | None
    dominant_fraction: float
    radius_m: float
    source: str = "point_source_id of the ground returns at the mark"

    @property
    def n_lines(self) -> int:
        return len(self.counts)

    @property
    def mixed(self) -> bool:
        return self.n_lines > 1

    @property
    def n(self) -> int:
        return int(sum(self.counts.values()))

    @staticmethod
    def table_columns() -> dict:
        return {
            "line": "point_source_id contributing the most ground returns within radius_m",
            "line_frac": "that line's share of the ground returns within radius_m, 0-1",
            "n_lines": "distinct point_source_id values within radius_m, count",
            "line_counts": "every point_source_id within radius_m with its return count",
        }


def assign_line_from_returns(ground: GroundReturns, easting: float, northing: float,
                             radius_m: float) -> LineAssignment:
    """Assign a mark to a flight line by counting ground returns per ``point_source_id``.

    This replaces assigning by distance to a fitted nadir track. The tracks in this
    acquisition are ~1 km apart and were fitted at one latitude; a 1-3 deg heading error
    displaces a track by several hundred metres over 20 km, so a centreline assignment
    mislabels marks exactly as the search widens. The returns are the acquisition's own
    record of which sortie illuminated the ground and cannot be wrong about it.
    """
    radius_m = float(radius_m)
    r = np.hypot(np.asarray(ground.x) - easting, np.asarray(ground.y) - northing)
    sel = r <= radius_m
    ids = np.asarray(ground.point_source_id)[sel]
    if ids.size == 0:
        return LineAssignment(counts={}, dominant=None, dominant_fraction=float("nan"),
                              radius_m=radius_m)
    vals, cnt = np.unique(ids, return_counts=True)
    counts = {int(v): int(c) for v, c in zip(vals, cnt)}
    k = int(np.argmax(cnt))
    return LineAssignment(counts=counts, dominant=int(vals[k]),
                          dominant_fraction=float(cnt[k] / cnt.sum()), radius_m=radius_m)


# ---------------------------------------------------------------------- siting screen

@dataclass(frozen=True)
class SitingScreen:
    """How the ground behaves AT a mark, measured before any tie is selected on.

    Every field is a property of the point cloud at the mark and is independent of the
    tie's value, so screening on it cannot select the answer. **No threshold lives here.**
    The one screen this project has actually tested -- ``radius_spread_mm`` -- did *not*
    reduce site-to-site scatter (``analysis/GEN1_DATUM_MODULE.md`` §5), so the cut is not
    only the caller's to choose, it is a cut they should first check is buying anything.
    """

    n: int
    slope_deg: float
    relief_mm: float
    fit_rms_mm: float
    radius_spread_mm: float
    radius_m: float

    @staticmethod
    def table_columns() -> dict:
        return {
            "n": "ground returns within the report radius, count",
            "slope_deg": "|grad S| of the local order-2 surface at the mark, degrees",
            "relief_mm": "p95 - p05 of return elevations within the report radius, mm",
            "fit_rms_mm": "RMS of z_i - S(x_i, y_i) for that surface, mm",
            "radius_spread_mm": ("max - min of the tie across the pipeline-scale radii "
                                 "(res/2 to 2*res), mm -- how much the answer depends on "
                                 "an unstated window"),
        }


def _screen_from(est: TieEstimate) -> SitingScreen:
    rep = next((r for r in est.curve if r.radius_m == est.report_radius_m), None)
    if rep is None or not rep.ok:
        return SitingScreen(n=(rep.n if rep else 0), slope_deg=float("nan"),
                            relief_mm=float("nan"), fit_rms_mm=float("nan"),
                            radius_spread_mm=est.radius_spread_mm,
                            radius_m=est.report_radius_m)
    return SitingScreen(n=rep.n, slope_deg=rep.slope_deg, relief_mm=rep.relief_mm,
                        fit_rms_mm=rep.fit_rms_mm,
                        radius_spread_mm=est.radius_spread_mm,
                        radius_m=est.report_radius_m)


# ------------------------------------------------------------------ swath constants

def swath_constants_from_corrections(path) -> tuple:
    """Per-swath ``(dx, dy, dz)`` from a tile's ``corrections.json``, and its provenance.

    Reads ``per_swath_internal_alignment_dxdydz_m``, written by
    :func:`lidar_diff_icp.coreg.align_swaths` -- gen1-internal swath-to-swath alignment,
    with no cross-epoch content. Returns ``(constants, source)`` where ``constants`` maps
    an int line id to a metre triple and ``source`` is a sentence naming the file and the
    tie method recorded in it, for the run banner.
    """
    path = str(path)
    with open(path) as fh:
        d = json.load(fh)
    key = "per_swath_internal_alignment_dxdydz_m"
    if key not in d:
        raise KeyError(f"{path} has no {key!r}; keys are {sorted(d)}")
    const = {int(k): tuple(float(v) for v in val) for k, val in d[key].items()}
    src = (f"{os.path.abspath(path)} {key}, swath_tie={d.get('swath_tie', '?')!r}, "
           f"ground_source={d.get('ground_source', '?')!r}, res_m={d.get('res_m', '?')}")
    return const, src


# ------------------------------------------------------------------- one measurement

@dataclass
class MarkMeasurement:
    """One mark measured: the tie, the line it came from, and how it is sited."""

    site: MarkSite
    tie: TieEstimate
    line: LineAssignment
    screen: SitingScreen
    swath_shift_m: tuple
    swath_constant_source: str
    per_line_tie_mm: dict = field(default_factory=dict)
    params: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def point_id(self) -> str:
        return self.site.point_id

    @property
    def tie_mm(self) -> float:
        return self.tie.tie_mm

    @property
    def sigma_mm(self) -> float:
        return self.tie.sigma_mm

    @property
    def line_id(self):
        return self.line.dominant

    @staticmethod
    def table_columns() -> dict:
        cols = {
            "point": "control mark id (first spelling in the bundled CSV)",
            "cover": ("MnDNR land-cover class of the mark: L1O open, L2T tall weeds/crops, "
                      "L3B brush/low trees, L4F forested, L5U urban"),
            "km": "distance from whatever the search was relative to, km",
            "tile": "the gen1 tile the mark falls in",
            "tie_mm": ("surveyed elevation minus gen1 ground at the mark, mm; the constant "
                       "to ADD to gen1. No geoid term (both sides GEOID03)"),
            "sigma_mm": "half the tie's radius spread over the pipeline-scale radii, mm",
        }
        cols.update(LineAssignment.table_columns())
        cols.update(SitingScreen.table_columns())
        return cols

    def table_row(self, columns) -> list:
        v = {
            "point": self.point_id, "cover": self.site.mark.cover_class,
            "km": f"{self.site.distance_m / 1000:.2f}",
            "tile": getattr(self, "tile", ""),
            "tie_mm": f"{self.tie_mm:+.1f}", "sigma_mm": f"{self.sigma_mm:.1f}",
            "line": self.line.dominant, "line_frac": f"{self.line.dominant_fraction:.2f}",
            "n_lines": self.line.n_lines,
            "line_counts": ", ".join(f"{k}:{v}" for k, v in sorted(self.line.counts.items())),
            "n": self.screen.n, "slope_deg": f"{self.screen.slope_deg:.1f}",
            "relief_mm": f"{self.screen.relief_mm:.0f}",
            "fit_rms_mm": f"{self.screen.fit_rms_mm:.0f}",
            "radius_spread_mm": f"{self.screen.radius_spread_mm:.1f}",
        }
        return [v.get(c, "") for c in columns]


def default_ground_loader(tile_path, easting, northing, half_width_m) -> GroundReturns:
    """Vendor class-2 ground in a square window about a mark.

    Vendor ground, not CSF: CSF is ~460 s per tile and the 2026-08-26 control run
    measured the whole verdict's ground-source dependence at 6.5 mm median absolute over
    16 marks. A caller who wants CSF passes their own loader (and their own crop
    half-width, which for CSF is a real parameter because the cloth sees the crop edge).
    """
    return vendor_ground_near(tile_path, easting, northing, half_width_m)


class TileGroundCache:
    """A ``ground_loader`` that decompresses each tile ONCE and crops in memory.

    Several marks usually share a tile, and re-reading a 20 MB LAZ per mark is the whole
    cost of a run on a shared laptop. This holds **one** tile's ground at a time and
    drops it when the tile changes, so peak memory is one tile's class-2 returns rather
    than the whole search. It changes nothing about the answer -- the crop it hands back
    is the same square :func:`default_ground_loader` would build -- and
    :func:`measure_sites` groups its work by tile so the cache hits.
    """

    def __init__(self, ground_class: int = 2):
        self.ground_class = int(ground_class)
        self._path = None
        self._all = None
        self.reads = 0

    def __call__(self, tile_path, easting, northing, half_width_m) -> GroundReturns:
        tile_path = str(tile_path)
        if tile_path != self._path:
            self._all = None            # free the previous tile before reading the next
            self._all = vendor_ground_near(tile_path, 0.0, 0.0, float("inf"),
                                           ground_class=self.ground_class)
            self._path = tile_path
            self.reads += 1
        g = self._all
        m = ((np.abs(g.x - easting) <= half_width_m)
             & (np.abs(g.y - northing) <= half_width_m))
        return GroundReturns(g.x[m], g.y[m], g.z[m], g.scan_angle[m],
                             g.point_source_id[m], g.source,
                             f"{g.origin} cropped +/-{half_width_m} m about "
                             f"({easting:.1f}, {northing:.1f})", int(m.sum()))


def measure_site(site: MarkSite, tile_path, *, res: float = 5.0,
                 ground_loader=default_ground_loader, crop_half_width_m=None,
                 swath_constants=None, swath_constants_source: str = "",
                 lateral_shift_m=None,
                 lidar_vertical_datum: str = GEN1_VERTICAL_DATUM,
                 lidar_geoid_model: str = GEN1_GEOID_MODEL,
                 per_line_ties: bool = True) -> MarkMeasurement:
    """Measure one mark: read the ground, assign its line, estimate the tie, screen it.

    ``res``               the pipeline grid resolution that sets the radius ladder
                          (repo default 5.0 m, ``corrections.json res_m``).
    ``crop_half_width_m`` window half-width handed to the loader. ``None`` DERIVES it as
                          the largest radius on the ladder (``5*res``) -- the smallest
                          square that contains every fitting window, so for a vendor-class
                          loader the answer cannot depend on it. It is not a chosen
                          number; for a CSF loader, where it *does* matter, pass one.
    ``swath_constants``   ``{line: (dx, dy, dz)}`` from
                          :func:`swath_constants_from_corrections`. When given, the
                          mark's OWN returns-assigned line's constant is applied before
                          estimating, putting the mark in the swath network's frame
                          (``mode="common_datum"``). A mark whose line has no constant is
                          measured unshifted and flagged; it is not dropped here.
    ``lateral_shift_m``   OPT-IN gen2-derived Nuth & Kaeaeb shift ``(dx, dy)``. Default
                          ``None`` = not applied, because a cross-epoch term does not
                          belong in gen1-against-its-own-control. Passing it is a
                          deliberate, recorded act.

    The geoid is asserted equal, never converted: ``geoid_shift_m`` is hard 0.0 and
    :func:`assert_no_geoid_conversion` runs first.
    """
    res = float(res)
    claim = assert_no_geoid_conversion([site.mark],
                                       lidar_vertical_datum=lidar_vertical_datum,
                                       lidar_geoid_model=lidar_geoid_model)
    ladder = radius_ladder(res)
    half = float(max(ladder)) if crop_half_width_m is None else float(crop_half_width_m)
    if half < max(ladder):
        raise ValueError(
            f"crop_half_width_m={half} is smaller than the largest fitting radius "
            f"{max(ladder)} m; the outer rungs of the radius ladder would be truncated "
            "by the crop and the ladder would silently stop meaning what it says.")
    cp = site.checkpoint
    ground = ground_loader(tile_path, cp.easting, cp.northing, half)

    report_radius = 1.5 * res
    line = assign_line_from_returns(ground, cp.easting, cp.northing, report_radius)

    dx, dy, dz = 0.0, 0.0, 0.0
    const_src, notes = "", []
    if swath_constants is not None:
        if line.dominant in swath_constants:
            dx, dy, dz = swath_constants[line.dominant]
            const_src = swath_constants_source
        else:
            notes.append(f"no swath constant for line {line.dominant}: measured unshifted "
                         "and NOT in the common frame")
    if lateral_shift_m is not None:
        dx += float(lateral_shift_m[0])
        dy += float(lateral_shift_m[1])
        notes.append(f"gen2-derived lateral shift {tuple(lateral_shift_m)} m applied by "
                     "explicit request; it is a cross-epoch term")
    if line.mixed:
        notes.append(f"{line.n_lines} flight lines at the mark "
                     f"({', '.join(f'{k}:{v}' for k, v in sorted(line.counts.items()))}) "
                     "-- the tie mixes their relative offsets")

    est = estimate_tie(cp, ground, line=None, res=res, swath_shift_m=(dx, dy, dz),
                       geoid_shift_m=0.0)

    per_line = {}
    if per_line_ties and line.mixed:
        for lid in sorted(line.counts):
            sdx, sdy, sdz = ((swath_constants or {}).get(lid, (0.0, 0.0, 0.0))
                             if swath_constants is not None else (0.0, 0.0, 0.0))
            if lateral_shift_m is not None:
                sdx += float(lateral_shift_m[0]); sdy += float(lateral_shift_m[1])
            try:
                e = estimate_tie(cp, ground, line=lid, res=res,
                                 swath_shift_m=(sdx, sdy, sdz), geoid_shift_m=0.0)
                per_line[lid] = e.tie_mm
            except Exception as exc:                                # pragma: no cover
                per_line[lid] = float("nan")
                notes.append(f"per-line tie for {lid} failed: {exc}")

    params = list(est.params) + [
        Param("geoid_shift_m", 0.0, "repo", _SRC_GEN1_DATUM + "; " + claim),
        Param("crop_half_width_m", half, "repo" if crop_half_width_m is None else "andy",
              "DERIVED as max(radius_ladder) = 5*res, the smallest square containing "
              "every fitting window, so a vendor-class read cannot depend on it"
              if crop_half_width_m is None else "supplied by the caller"),
        Param("line_assignment", line.dominant, "repo",
              f"dominant point_source_id of the {line.n} ground returns within "
              f"{report_radius} m of the mark -- from the RETURNS, not from distance to a "
              "fitted centreline"),
        Param("swath_constants", ("applied" if const_src else "none"),
              "repo" if const_src else "repo",
              const_src or "mode=per_line: each line's swath constant is left unknown"),
        Param("lateral_shift_m", tuple(lateral_shift_m) if lateral_shift_m else None,
              "andy" if lateral_shift_m else "repo",
              "gen2-derived Nuth & Kaeaeb shift, applied on explicit request"
              if lateral_shift_m else
              "NOT applied: a cross-epoch term has no place in gen1 vs its own control"),
    ]
    out = MarkMeasurement(site=site, tie=est, line=line, screen=_screen_from(est),
                          swath_shift_m=(dx, dy, dz), swath_constant_source=const_src,
                          per_line_tie_mm=per_line, params=params, notes=notes)
    out.tile = os.path.splitext(os.path.basename(str(tile_path)))[0]
    return out


def measure_sites(sites, resolution: TileResolution, *, on_missing: str = "skip",
                  **kw) -> tuple:
    """Measure every site whose tile is on disk. Returns ``(measurements, skipped)``.

    ``on_missing="skip"`` records the mark in ``skipped`` with its tile name; nothing is
    fetched. ``on_missing="raise"`` refuses to produce a partial answer.

    Work is grouped by tile so a :class:`TileGroundCache` (the default loader here)
    decompresses each tile once and holds only one at a time; the measurements come back
    in the order the sites were given, not in tile order.
    """
    sites = list(sites)
    kw.setdefault("ground_loader", TileGroundCache())
    order = {id(s): k for k, s in enumerate(sites)}
    by_tile = {}
    skipped = []
    for s in sites:
        p = resolution.path_for(s.point_id)
        if p is None:
            if on_missing == "raise":
                raise FileNotFoundError(
                    f"{s.point_id}: tile {resolution.per_mark.get(s.point_id)} is not on "
                    "disk. This module never downloads; fetch it deliberately, one tile "
                    "at a time, with tiles.download_tile.")
            skipped.append((s.point_id, resolution.per_mark.get(s.point_id), "tile not on disk"))
            continue
        by_tile.setdefault(p, []).append(s)

    done = []
    for p, group in by_tile.items():                # one tile decompressed at a time
        for s in group:
            try:
                done.append((order[id(s)], measure_site(s, p, **kw)))
            except Exception as exc:
                skipped.append((s.point_id, resolution.per_mark.get(s.point_id),
                                f"{type(exc).__name__}: {exc}"))
    return [m for _, m in sorted(done, key=lambda t: t[0])], skipped
