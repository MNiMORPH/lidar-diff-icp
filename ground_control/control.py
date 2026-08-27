"""Epoch-agnostic access to the bundled control tables.

`residual_field.py` already carries the machinery this subsystem needs -- variogram
fitting, kriging with a prediction variance, LOO and spatially blocked CV -- and every
one of those entry points (`fit_field`, `krige`, `krige_many`, `loo_errors`, `block_cv`,
`cover_design`) takes plain ``(x, y, v)`` arrays and is already epoch-agnostic.

Only its *edges* are not.  `load_residuals`, `check_sign_convention` and `stratify` are
written against gen1's schema (``point_type``, ``dnr_error_m``); gen2's table has a
different one (``role``, and FOUR residual columns).  Making `residual_field` itself
epoch-agnostic would mean editing `src/`, which this session may not do -- so this module
adapts both tables INTO `residual_field.ControlResiduals` and hands them to the existing,
untouched machinery.  That is a deliberate adapter, not a fork: no estimator is
reimplemented here.

If this subsystem is ever promoted into `src/lidar_diff_icp/`, the right move is to move
these loaders into `residual_field` and delete this module.  See INTEGRATION.md.

Sign convention, verified rather than inherited -- :func:`verify_sign_convention`
re-derives it on every row of both tables:

    tie = surveyed - z_lidar        POSITIVE = the surface reads LOW (add the constant)

gen1's ``dnr_error_m`` is ``Control Z - Surface Z`` and gen2's four ``*_error_m`` columns
are ``Z - <surface>z``; both are that same subtraction, so no flip separates them.

**gen2 has four surfaces, not one.**  The delivered DEM and the delivered point cloud, at
each of two quality levels, are four different answers to "what does the surface read
here".  ``surface=`` is therefore REQUIRED for gen2 and has no default: picking one
silently would hide a choice that changes the number.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_PKG = Path(__file__).resolve().parent.parent / "src" / "lidar_diff_icp" / "groundtruth"
DATA_DIR = _PKG / "data"

GEN1_CSV = DATA_DIR / "mn_dnr_2008_control_semn.csv"
GEN2_CSV = DATA_DIR / "mn_se_driftless_2021_control.csv"

EPOCHS = ("gen1", "gen2")

#: gen2's four delivered surfaces.  ``<quality-level>_<product>``.
GEN2_SURFACES = ("ql1_dem", "ql1_laz", "ql0_dem", "ql0_laz")

#: The horizontal CRS each table's coordinates are published in.  These are carried and
#: asserted, never silently mixed: PROJ returns a null transform between them here, but
#: the NAD83(1986) -> NAD83(2011) realization difference is not modelled by that null
#: transform and is therefore NOT accounted for by treating them as equal.
EPOCH_CRS = {"gen1": "EPSG:26915", "gen2": "EPSG:6344"}

#: gen1's cover taxonomy is MnDNR's five classes; gen2's is USGS's binary accuracy class.
#: They are DIFFERENT TAXONOMIES.  Nothing in this module maps one onto the other; any
#: such mapping is a scientific choice and belongs to the caller, stated as one.
GEN1_COVER_CLASSES = ("L1O", "L2T", "L3B", "L4F", "L5U", "other")
GEN2_COVER_CLASSES = ("NVA", "VVA", "LCP")


def _f(v):
    v = (v or "").strip()
    return float(v) if v else None


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ------------------------------------------------------------------ sign convention

@dataclass(frozen=True)
class SignCheck:
    """Which subtraction the table's residual column actually is, re-derived per row."""

    epoch: str
    column: str
    n_rows_checked: int
    n_exact_surveyed_minus_lidar: int
    n_exact_lidar_minus_surveyed: int
    worst_miss_surveyed_minus_lidar_m: float
    worst_miss_lidar_minus_surveyed_m: float

    @property
    def is_surveyed_minus_lidar(self) -> bool:
        return (self.n_rows_checked > 0
                and self.n_exact_surveyed_minus_lidar == self.n_rows_checked)


def _sign_check(epoch, rows, z_col, err_col, column_label, tol_m):
    okf = okr = n = 0
    wf = wr = 0.0
    for r in rows:
        z, s, e = _f(r["elevation"]), _f(r[z_col]), _f(r[err_col])
        if None in (z, s, e):
            continue
        n += 1
        df, dr = abs((z - s) - e), abs((s - z) - e)
        okf += df < tol_m
        okr += dr < tol_m
        wf, wr = max(wf, df), max(wr, dr)
    return SignCheck(epoch, column_label, n, okf, okr, wf, wr)


def verify_sign_convention(epoch: str, *, tol_m: float, surface: str | None = None):
    """Re-derive the sign convention on every row.  ``tol_m`` is required, not assumed."""
    if epoch == "gen1":
        return [_sign_check("gen1", _rows(GEN1_CSV), "dnr_surface_z_m", "dnr_error_m",
                            "dnr_error_m", tol_m)]
    if epoch != "gen2":
        raise ValueError(f"epoch={epoch!r} must be one of {EPOCHS}")
    rows = _rows(GEN2_CSV)
    surfaces = GEN2_SURFACES if surface is None else (surface,)
    return [_sign_check("gen2", rows, f"usgs_{s}_z_m", f"usgs_{s}_error_m",
                        f"usgs_{s}_error_m", tol_m) for s in surfaces]


# ------------------------------------------------------------------------- loading

@dataclass(frozen=True)
class ControlLoad:
    """A loaded control set, plus everything about HOW it was selected.

    ``residuals`` is a :class:`residual_field.ControlResiduals`, so every estimator in
    that module accepts it unchanged.  The rest of the fields exist so a run cannot
    report a number without also being able to report which rows produced it.
    """

    residuals: object          # residual_field.ControlResiduals
    epoch: str
    crs: str
    surface: str | None
    n_rows_in: int
    n_with_residual: int
    n_marks_out: int
    n_dup_rows: int
    n_dup_groups: int
    roles_present: dict
    covers_present: dict
    dropped_no_residual: int


def load_control(epoch: str, *, surface: str | None = None,
                 roles: tuple | None = None) -> ControlLoad:
    """Load one epoch's control into ``residual_field.ControlResiduals``.

    Parameters
    ----------
    epoch : {"gen1", "gen2"}
    surface : str, required for gen2
        Which delivered surface the residual describes; one of :data:`GEN2_SURFACES`.
        There is no default -- gen2 publishes four and they are four different answers.
    roles : tuple, optional
        Restrict to these ``role`` values (gen2 only; gen1's table has no role column).
        Left as ``None`` this does NOT silently filter: it keeps every row that carries a
        residual, and reports in ``roles_present`` which roles survived.  gen2's 143 LCPs
        carry no residual in any of the four columns, so they drop out as a FACT of the
        table rather than as a choice made here -- see REPORT.md.

    De-duplication is on the exact ``(easting, northing, elevation)`` triple, no
    tolerance, keeping the first occurrence -- a tolerance would be an invented parameter.
    A mark on a county line is printed in both counties' validation reports and must not
    enter a variogram twice.
    """
    import sys
    sys.path.insert(0, str(_PKG.parent.parent))
    from lidar_diff_icp.groundtruth.residual_field import ControlResiduals

    if epoch not in EPOCHS:
        raise ValueError(f"epoch={epoch!r} must be one of {EPOCHS}")
    if epoch == "gen2" and surface is None:
        raise ValueError(
            "epoch='gen2' requires surface=; it publishes four residual columns "
            f"({', '.join(GEN2_SURFACES)}) and they are four different answers. "
            "Picking one here would hide a choice that changes the number.")
    if epoch == "gen1" and surface is not None:
        raise ValueError("epoch='gen1' publishes ONE surface; surface= is meaningless")
    if epoch == "gen2" and surface not in GEN2_SURFACES:
        raise ValueError(f"surface={surface!r} must be one of {GEN2_SURFACES}")

    path = GEN1_CSV if epoch == "gen1" else GEN2_CSV
    rows = _rows(path)
    n_rows_in = len(rows)
    err_col = "dnr_error_m" if epoch == "gen1" else f"usgs_{surface}_error_m"

    keep_rows = [r for r in rows if _f(r[err_col]) is not None]
    n_with_residual = len(keep_rows)
    if roles is not None:
        if epoch == "gen1":
            raise ValueError("gen1's table has no 'role' column; roles= is meaningless")
        keep_rows = [r for r in keep_rows if r["role"] in roles]

    seen, keep, dup = {}, [], 0
    for r in keep_rows:
        k = (r["easting"], r["northing"], r["elevation"])
        if k in seen:
            seen[k] += 1
            dup += 1
            continue
        seen[k] = 1
        keep.append(r)
    n_dup_groups = sum(1 for v in seen.values() if v > 1)

    county = ([r.get("county", "") for r in keep] if epoch == "gen1"
              else [r.get("va_blocks", "") for r in keep])

    cr = ControlResiduals(
        point_id=np.array([r["point_id"] for r in keep]),
        county=np.array(county),
        cover=np.array([r["point_type"] for r in keep]),
        easting=np.array([float(r["easting"]) for r in keep]),
        northing=np.array([float(r["northing"]) for r in keep]),
        resid_mm=np.array([float(r[err_col]) * 1000.0 for r in keep]),
        n_rows_in=n_rows_in,
        n_dup_rows=dup,
        n_dup_groups=n_dup_groups,
    )

    def _count(key):
        out = {}
        for r in keep:
            out[r.get(key, "")] = out.get(r.get(key, ""), 0) + 1
        return dict(sorted(out.items()))

    return ControlLoad(
        residuals=cr, epoch=epoch, crs=EPOCH_CRS[epoch], surface=surface,
        n_rows_in=n_rows_in, n_with_residual=n_with_residual, n_marks_out=len(keep),
        n_dup_rows=dup, n_dup_groups=n_dup_groups,
        roles_present=_count("role") if epoch == "gen2" else {},
        covers_present=_count("point_type"),
        dropped_no_residual=n_rows_in - n_with_residual,
    )


def cover_from_point_id_prefix(point_id):
    """The BROKEN way to read gen1's cover, kept so its breakage can be demonstrated.

    21 gen1 marks carry a ``point_id`` beginning ``L10`` with a DIGIT ZERO rather than
    ``L1O`` with the letter O.  Parsing cover from the prefix silently drops them, giving
    209 open marks instead of the CSV's own 230.  Never use this to select marks; it
    exists so ``tests/`` can assert the difference and keep the trap documented.
    """
    return str(point_id)[:3]
