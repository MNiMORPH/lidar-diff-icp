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
