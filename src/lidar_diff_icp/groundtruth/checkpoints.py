"""Surveyed checkpoints, with their datum stated rather than assumed.

A tie between a lidar epoch and a surveyed mark is only as good as the agreement of
their vertical datums. gen1 is NAVD88 on **GEOID03**; the 2021 3DEP checkpoints are
NAVD88 on **GEOID18**; the two differ by ~67 mm at Elba. A checkpoint whose geoid model
is not recorded cannot be used for a tie at all, because the error it would introduce is
exactly the size of the thing being measured. So :class:`Checkpoint` carries the datum
as data, and :meth:`CheckpointSet.usable` **raises** rather than guessing.

Two sources:

* :func:`read_3dep_va_shapefile` -- the authoritative USGS 3DEP vertical-accuracy point
  shapefile, where it is on disk.
* :func:`load_bundled` -- a checked-in CSV transcription for offline work
  (``data/mn_se_driftless_2021_ql1_near_elba.csv``; see that directory's README).

Neither reader converts units. The 3DEP contractor shapefile mislabels its elevation
field as US Feet when the values are metres, so ``source_ele_units`` has no default:
a caller must say what the numbers are.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from pathlib import Path

_DATA = Path(__file__).with_name("data")

#: Strings that mean "nobody recorded the geoid model". Matched case-insensitively
#: after stripping. Anything here is a refusal, not a default.
UNKNOWN_TOKENS = ("", "unknown", "none", "null", "n/a", "na", "?", "nan")


class UnknownDatumError(ValueError):
    """A checkpoint's vertical datum or geoid model is not recorded.

    Raised instead of assuming a model. A silent geoid mismatch is a ~70 mm error at
    Elba -- the same size as the offset a tie is trying to measure -- so it must not be
    possible to reach a tie number without having stated the model.
    """


@dataclass(frozen=True)
class Checkpoint:
    """One surveyed mark.

    ``point_type`` follows the 3DEP vertical-accuracy convention: ``"NVA"`` = non-
    vegetated vertical accuracy (open, hard ground), ``"VVA"`` = vegetated vertical
    accuracy (under canopy, where the published spread is 2-3x larger). It is kept as
    reported and never used to filter -- a VVA mark is a noisier tie, not an invalid one,
    and the diagnostics are what should say so.
    """

    point_id: str
    point_type: str
    easting: float
    northing: float
    elevation: float
    elevation_units: str
    horizontal_crs: str
    vertical_datum: str
    geoid_model: str
    project_id: str = ""
    collected: str = ""
    source: str = ""
    verified: str = ""

    @property
    def datum_known(self) -> bool:
        """True when both the vertical datum and the geoid model are recorded."""
        return not (self.vertical_datum.strip().lower() in UNKNOWN_TOKENS
                    or self.geoid_model.strip().lower() in UNKNOWN_TOKENS)

    def require_datum(self) -> "Checkpoint":
        """Return self, or raise :class:`UnknownDatumError` if the datum is not recorded."""
        if not self.datum_known:
            raise UnknownDatumError(
                f"checkpoint {self.point_id!r}: vertical_datum="
                f"{self.vertical_datum!r} geoid_model={self.geoid_model!r}. "
                "A tie cannot be computed without the geoid model -- the GEOID03/GEOID18 "
                "difference alone is ~67 mm at Elba. Record the model or drop the point.")
        return self

    @property
    def elevation_m(self) -> float:
        """Elevation in metres. Raises if the recorded units are not metres.

        No conversion is performed: the one unit problem this data actually has is a
        *mislabel* (see the data README), and a converter would launder it.
        """
        u = self.elevation_units.strip().lower()
        if u in ("m", "metre", "metres", "meter", "meters"):
            return float(self.elevation)
        raise ValueError(
            f"checkpoint {self.point_id!r}: elevation_units={self.elevation_units!r}. "
            "This reader does not convert units; state the elevation in metres, having "
            "checked what the source actually holds.")


@dataclass
class CheckpointSet:
    """A set of checkpoints plus where it came from."""

    points: list[Checkpoint]
    origin: str = ""
    fields: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.points)

    def __iter__(self):
        return iter(self.points)

    def __getitem__(self, key):
        if isinstance(key, str):
            for p in self.points:
                if p.point_id == key:
                    return p
            raise KeyError(key)
        return self.points[key]

    @property
    def ids(self) -> list[str]:
        return [p.point_id for p in self.points]

    def usable(self) -> list[Checkpoint]:
        """Every point, each having passed :meth:`Checkpoint.require_datum`.

        Raises on the first point with an unrecorded datum. It does **not** silently
        drop it: a control point that cannot be used is a reported result.
        """
        return [p.require_datum() for p in self.points]

    def within(self, bounds) -> "CheckpointSet":
        """Subset falling inside ``bounds`` = (x0, y0, x1, y1), same CRS as the points."""
        x0, y0, x1, y1 = bounds
        keep = [p for p in self.points
                if x0 <= p.easting <= x1 and y0 <= p.northing <= y1]
        return CheckpointSet(keep, origin=f"{self.origin} within {bounds}", fields=self.fields)


# ----------------------------------------------------------------------- readers

def list_bundled() -> list[str]:
    """Names of the checked-in checkpoint CSVs (without the ``.csv``)."""
    return sorted(p.stem for p in _DATA.glob("*.csv"))


def load_bundled(name: str = "mn_se_driftless_2021_ql1_near_elba") -> CheckpointSet:
    """Load a checked-in CSV transcription of a surveyed checkpoint set.

    The CSV carries the datum per row (``horizontal_crs``, ``vertical_datum``,
    ``geoid_model``, ``elevation_units``) so that a subset extracted from it cannot lose
    the metadata. Missing columns raise here rather than defaulting.
    """
    path = _DATA / f"{name}.csv"
    if not path.exists():
        raise FileNotFoundError(f"no bundled checkpoint set {name!r}; have {list_bundled()}")
    return read_checkpoint_csv(path)


def read_checkpoint_csv(path) -> CheckpointSet:
    """Read a checkpoint CSV whose header names the fields of :class:`Checkpoint`."""
    required = ("point_id", "point_type", "easting", "northing", "elevation",
                "elevation_units", "horizontal_crs", "vertical_datum", "geoid_model")
    pts = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        missing = [c for c in required if c not in (rd.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing required column(s) {missing}; "
                             f"have {rd.fieldnames}")
        for row in rd:
            pts.append(Checkpoint(
                point_id=row["point_id"], point_type=row["point_type"],
                easting=float(row["easting"]), northing=float(row["northing"]),
                elevation=float(row["elevation"]),
                elevation_units=row["elevation_units"],
                horizontal_crs=row["horizontal_crs"],
                vertical_datum=row["vertical_datum"],
                geoid_model=row["geoid_model"],
                project_id=row.get("project_id", ""),
                collected=row.get("collected", ""),
                source=row.get("source", ""),
                verified=row.get("verified", "")))
    return CheckpointSet(pts, origin=str(Path(path).resolve()))


#: Field names looked for in a 3DEP vertical-accuracy shapefile, in order. The USGS
#: contractor deliveries are not schema-stable, so each role lists the spellings seen or
#: documented rather than one hard-coded name. If none matches, the reader raises and
#: prints the fields the file actually has -- it never guesses.
VA_FIELDS = {
    "point_id": ("unique_ind", "unique_id", "point_id", "pointid", "id", "name"),
    "point_type": ("nva_vva", "va_type", "point_type", "type", "class", "category"),
    "elevation": ("source_ele", "source_elev", "survey_ele", "elev", "elevation", "z"),
    "geoid_model": ("geoid", "geoid_mode", "geoid_model", "geoidmodel", "vert_geoid"),
    "vertical_datum": ("vert_datum", "vdatum", "vertical_d", "v_datum"),
    "project_id": ("prj_id", "project_id", "proj_id", "workunit"),
    "collected": ("collect_da", "collected", "survey_dat", "date"),
}


def _pick(fields, candidates):
    low = {f.lower(): f for f in fields}
    for c in candidates:
        if c in low:
            return low[c]
    return None


def read_3dep_va_shapefile(path, *, source_ele_units: str,
                           horizontal_crs: str | None = None,
                           vertical_datum: str | None = None,
                           field_map: dict | None = None) -> CheckpointSet:
    """Read a USGS 3DEP vertical-accuracy checkpoint shapefile.

    ``source_ele_units`` has **no default and must be passed**. The contractor
    shapefile for MN_SE_Driftless_2021_B21 labels its elevation field "US Feet" while
    the values are metres (verified against USGS EPQS; see the data README). A default
    here would let that mislabel propagate into a tie silently, which is precisely the
    failure this module exists to prevent.

    The geoid model is read from the file. If the file has no geoid field and none is
    supplied through ``field_map``, this raises :class:`UnknownDatumError` listing the
    fields present -- it does not fall back to GEOID18.

    ``horizontal_crs`` defaults to the ``.prj`` sidecar if pyproj can read it.
    Needs ``pyshp`` (``import shapefile``).
    """
    import shapefile  # pyshp

    path = str(path)
    r = shapefile.Reader(path)
    names = [f[0] for f in r.fields[1:]]           # drop the DeletionFlag pseudo-field
    fm = dict(field_map or {})
    for role, cands in VA_FIELDS.items():
        fm.setdefault(role, _pick(names, cands))

    if not fm.get("elevation"):
        raise ValueError(f"{path}: no elevation field found among {names}; "
                         "pass field_map={'elevation': '<name>'}")
    if not fm.get("geoid_model"):
        raise UnknownDatumError(
            f"{path}: no geoid-model field among {names}. The geoid model is not "
            "optional -- GEOID03 vs GEOID18 is ~67 mm at Elba. Pass "
            "field_map={'geoid_model': '<name>'} once you have confirmed which field "
            "holds it, or use a source that records it.")

    if horizontal_crs is None:
        prj = os.path.splitext(path)[0] + ".prj"
        if os.path.exists(prj):
            try:
                from pyproj import CRS
                c = CRS.from_wkt(open(prj).read())
                horizontal_crs = c.to_string() if c.to_epsg() is None else f"EPSG:{c.to_epsg()}"
            except Exception:
                horizontal_crs = f"(unparsed .prj: {prj})"
        else:
            horizontal_crs = ""

    idx = {k: (names.index(v) if v else None) for k, v in fm.items()}
    pts = []
    for sr in r.iterShapeRecords():
        rec = sr.record
        x, y = sr.shape.points[0][:2]

        def g(role, default=""):
            i = idx.get(role)
            return default if i is None else str(rec[i]).strip()

        pts.append(Checkpoint(
            point_id=g("point_id") or f"{x:.1f}_{y:.1f}",
            point_type=g("point_type"),
            easting=float(x), northing=float(y),
            elevation=float(rec[idx["elevation"]]),
            elevation_units=source_ele_units,
            horizontal_crs=horizontal_crs or "",
            vertical_datum=(vertical_datum if vertical_datum is not None
                            else g("vertical_datum") or "NAVD88"),
            geoid_model=g("geoid_model"),
            project_id=g("project_id"), collected=g("collected"),
            source=f"{os.path.basename(path)} field {fm['elevation']}",
            verified="read from shapefile"))
    return CheckpointSet(pts, origin=os.path.abspath(path), fields=fm)
