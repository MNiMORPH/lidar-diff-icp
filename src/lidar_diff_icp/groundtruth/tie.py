"""The tie estimator: lidar ground elevation at a surveyed mark, with a radius curve.

Why this is not a one-liner
---------------------------
The obvious estimator -- fit a plane to the ground returns within R of the mark and read
it at the mark -- is **strongly radius-dependent** on real checkpoints. Measured on gen1
at 3DEP checkpoint 2210 (line 128, 201 class-2 returns within 10 m), lidar minus
surveyed:

===========  ========  ========  ========  ========  ========
R (m)             2         5        10        15        20
plane          -109      -200      -589      -977     -1169
plane RMS         39       132       543       776       824
===========  ========  ========  ========  ========  ========

(mm.) A number that moves by a metre depending on an unstated radius is not a
measurement.

**The cause is curvature, not noise.** The mark sits on a local topographic high -- a
road crown/shoulder, with the ground falling away within 10 m; the surveyed 349.288 m
sits at the p95 of gen1 returns inside 5 m. A least-squares *plane* has no curvature
term, so over a convex patch of curvature k its value at the centre is pulled low by
~k*R^2/4: quadratic in R, which is exactly the shape of the row above (-109, -200, -589,
-1169 against R^2 = 4, 25, 100, 400).

**The fix is the one this codebase already made for the same reason.**
:func:`lidar_diff_icp.pipeline.difference_dem` carries a ``ground="poly2"`` estimator
whose docstring says it plainly: a windowed 2nd-order polynomial read at its constant
term is "curvature-UNBIASED, unlike the per-cell median (which carries the cell's
curvature) or a plane (which has no curvature term)". Using order 2 here, on the same
data:

===========  ========  ========  ========  ========  ========
R (m)           2.5         5       7.5        10        15
order 2         -75       -85       -86       -81      -284
===========  ========  ========  ========  ========  ========

Flat to ~11 mm across a 4x range of radius, then breaking down at 15 m where the patch
leaves a quadratic's reach (and says so: the fit RMS goes 19 -> 43 -> 101 -> 183 -> 385 mm).

What the estimator is
---------------------
The project's slope-normal ground read, generalised from a cell centre to an arbitrary
point. In :func:`~lidar_diff_icp.pipeline.difference_dem` with ``ground="slope_normal"``
the ground of a cell is

    Zreg(cell centre) + quantile_q( z_i - [Zreg + dE*dZreg/dE + dN*dZreg/dN] )

i.e. *a smooth local reference surface, read at the target point, plus a quantile of the
vertical residuals to it*. Here the smooth reference surface is fitted locally, at order
``surface_order`` (2 by default, for the reason above), and the target point is the
checkpoint instead of a cell centre:

    z_lidar(mark) = S(mark) + quantile_q( z_i - S(x_i, y_i) )

with ``q = 0.50``, the pipeline's ``ground_percentile``. Setting ``surface_order=1``
reproduces the plane pathology, which is how the regression test bites.

Radius is a first-class output
------------------------------
:func:`estimate_tie` never returns a single number. It returns the whole
:class:`RadiusEstimate` ladder, and the headline uncertainty is *how far the estimate
moves across the radii the pipeline itself works at* (cell half-width ``res/2`` to
``2*res``), not the standard error of a fit -- because on this data the radius term
dominates and the standard error is optimistic by an order of magnitude.

Whether a given checkpoint is usable is a question about tolerance, and the tolerance
belongs to the caller: :meth:`TieEstimate.verdict` requires an explicit
``tolerance_mm`` and states what it was. Nothing in this module drops a point.
"""
from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass, field

import numpy as np

from .provenance import Param

#: Where the repo's own values come from, quoted in the Param records.
_SRC_RES = "pipeline.difference_dem(res=5.0); data/derived/elbaext/corrections.json res_m"
_SRC_Q = ("pipeline.difference_dem(ground_q=0.50); data/derived/elbaext/corrections.json "
          "ground_percentile = 0.5")
_SRC_ORDER = ("pipeline._poly2_ground is order 2 and its docstring gives the reason: a plane "
              "'has no curvature term'; the checkpoints sit on local highs")
_SRC_WINDOW = ("pipeline._poly2_ground reads its constant term over a 3x3 cell window, i.e. a "
               "half-width of 1.5*res")
_SRC_GROUND = "data/derived/elbaext/corrections.json ground_source = csf"


def scan_angle_deg(las) -> np.ndarray:
    """Scan angle in DEGREES from a laspy object, both LAS point-format families.

    Point formats <= 5 carry ``scan_angle_rank`` in integer degrees (0 = nadir); 6+
    carry ``scan_angle`` in 0.006-degree units. This is the same two-branch read as
    ``analysis/ridgelines/gen1_save_angles_slope.py``. It RAISES when neither dimension
    is present rather than returning zeros: PDAL rewrites a point-format-1 crop as
    format 7, and a silent zero here would print "all beams at nadir" for a swath edge.
    """
    dims = set(las.point_format.dimension_names)
    if "scan_angle" in dims:
        return np.asarray(las.scan_angle).astype(float) * 0.006
    if "scan_angle_rank" in dims:
        return np.asarray(las.scan_angle_rank).astype(float)
    raise ValueError(f"no scan angle dimension in {sorted(dims)}")


@dataclass
class GroundReturns:
    """Ground-classified returns near one point, plus how they were classified."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    scan_angle: np.ndarray
    point_source_id: np.ndarray
    source: str                 # "csf" | "vendor_class2" | ...
    origin: str = ""            # file(s) the returns came from
    n_input: int = 0            # returns in the crop before ground selection

    def __len__(self) -> int:
        return int(self.x.size)

    def for_line(self, line: int | None) -> "GroundReturns":
        """Subset to one flight line (``None`` = keep all)."""
        if line is None:
            return self
        m = self.point_source_id == line
        return GroundReturns(self.x[m], self.y[m], self.z[m], self.scan_angle[m],
                             self.point_source_id[m], self.source,
                             f"{self.origin} line {line}", self.n_input)

    def shifted(self, dx: float, dy: float, dz: float) -> "GroundReturns":
        """Apply a rigid 3-D shift (a swath alignment, or a chain-accumulated one)."""
        return GroundReturns(self.x + dx, self.y + dy, self.z + dz, self.scan_angle,
                             self.point_source_id, self.source,
                             f"{self.origin} shifted by ({dx:+.4f},{dy:+.4f},{dz:+.4f}) m",
                             self.n_input)


@dataclass
class RadiusEstimate:
    """The ground read, and its diagnostics, at one radius."""

    radius_m: float
    n: int
    z_lidar_m: float
    fit_rms_mm: float           # RMS of z_i - S(x_i,y_i); how well the surface describes the patch
    median_resid_mm: float      # the quantile term itself, the part beyond the LS surface
    slope_deg: float            # |grad S| at the mark
    relief_mm: float            # p95 - p05 of the returns' z inside the radius
    scan_p50: float
    scan_p05: float
    scan_p95: float
    n_lines: int
    ok: bool = True
    note: str = ""


@dataclass
class TieEstimate:
    """Lidar ground elevation at a checkpoint, its radius curve, and the tie.

    ``tie_mm`` is the constant to **ADD to gen1** (already in the reference-swath frame
    and already geoid-shifted to the checkpoint's geoid model) to place it on the
    surveyed datum: ``tie = surveyed - z_lidar_corrected``, read at ``report_radius_m``.

    ``tie_median_mm`` is the same quantity taken as the median over the pipeline-scale
    radii. It is a robustness companion, not a competing answer: where the two differ by
    more than ``sigma_mm`` the radius curve is wobbling and the curve, not either number,
    is the result.
    """

    point_id: str
    point_type: str
    line: int | None
    curve: list[RadiusEstimate]
    report_radius_m: float
    z_lidar_raw_m: float          # at report_radius, before any datum term
    swath_shift_m: tuple          # (dx, dy, dz) applied before estimating
    geoid_shift_m: float          # added to gen1 z (N_gen1 - N_checkpoint)
    surveyed_m: float
    tie_mm: float
    tie_median_mm: float          # median of the tie over the pipeline-scale radii
    radius_spread_mm: float       # max-min of tie over the pipeline-scale radii
    radius_spread_all_mm: float   # ... and over the whole ladder
    fit_se_mm: float              # fit_rms/sqrt(n) at report_radius -- OPTIMISTIC, see docs
    sigma_mm: float               # reported uncertainty = half the pipeline-scale spread
    params: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def n_report(self) -> int:
        for r in self.curve:
            if r.radius_m == self.report_radius_m:
                return r.n
        return 0

    def verdict(self, tolerance_mm: float, *, tolerance_source: str) -> tuple[bool, str]:
        """Usable or not, against a tolerance the CALLER supplies and names.

        There is no default. A checkpoint is "unusable" only relative to how much
        radius sensitivity you are willing to carry, and that is a decision about the
        science, not a property of the code. The returned string states the tolerance
        and its source so it appears in any output that prints the verdict.
        """
        why = (f"radius spread {self.radius_spread_mm:.1f} mm over "
               f"{self._pipeline_radii()} m vs tolerance {tolerance_mm:.1f} mm "
               f"({tolerance_source})")
        if not any(r.ok for r in self.curve):
            return False, "no radius produced an estimate: " + "; ".join(
                r.note for r in self.curve if r.note)
        return (self.radius_spread_mm <= tolerance_mm), why

    def _pipeline_radii(self):
        return [r.radius_m for r in self.curve if r.ok and r.radius_m <= self.report_radius_m
                * 4.0 / 3.0]

    def table_rows(self):
        """Rows for the radius curve, matching :meth:`table_columns`."""
        rows = []
        for r in self.curve:
            if not r.ok:
                rows.append([f"{r.radius_m:.1f}", r.n, "-", "-", "-", "-", "-", "-", r.note])
                continue
            tie = 1000.0 * (self.surveyed_m - (r.z_lidar_m + self.geoid_shift_m))
            rows.append([f"{r.radius_m:.1f}", r.n, f"{tie:+.1f}", f"{r.fit_rms_mm:.0f}",
                         f"{r.median_resid_mm:+.1f}", f"{r.slope_deg:.1f}",
                         f"{r.relief_mm/1000:.2f}", r.n_lines,
                         f"{r.scan_p05:.0f}/{r.scan_p50:.0f}/{r.scan_p95:.0f}"])
        return rows

    @staticmethod
    def table_columns() -> dict:
        """Column name -> definition (with units), for ``trust.provenance.Run.column``."""
        return {
            "R_m": "radius of the fitting window about the checkpoint, m",
            "n": "ground returns inside that radius, count",
            "tie_mm": ("surveyed elevation minus lidar ground at the mark, after the swath "
                       "and geoid terms; the constant to ADD to gen1, mm"),
            "fit_rms_mm": "RMS of z_i - S(x_i,y_i) for the local surface fit S, mm",
            "medres_mm": "median of those residuals -- the slope-normal quantile term, mm",
            "slope_deg": "|grad S| at the checkpoint, degrees",
            "relief_m": "p95 - p05 of return elevations inside the radius, m",
            "n_lines": "distinct flight lines contributing returns inside the radius, count",
            "scan_deg": "|scan angle| p05/p50/p95 of those returns, degrees (0 = nadir)",
        }


# --------------------------------------------------------------------- estimator core

def _design(u, v, order):
    cols = [np.ones_like(u), u, v]
    if order >= 2:
        cols += [u * u, v * v, u * v]
    return np.column_stack(cols)


def ground_elevation_at(x, y, z, px, py, radius, *, surface_order=2, quantile=0.50):
    """Slope-normal ground elevation at ``(px, py)`` from returns within ``radius``.

    Fits ``S`` of order ``surface_order`` to the returns in the window by least squares,
    reads it at the target point, and adds the ``quantile`` of the vertical residuals --
    the same two-step form the pipeline's ``ground="slope_normal"`` estimator uses at a
    cell centre (smooth local reference surface + residual quantile).

    Returns ``(z_hat, info)``; ``z_hat`` is NaN and ``info["note"]`` says why when the
    window holds too few points to determine the surface (3 coefficients at order 1,
    6 at order 2) or the normal matrix is singular.
    """
    r = np.hypot(x - px, y - py)
    m = r <= radius
    n = int(m.sum())
    k = 3 if surface_order == 1 else 6
    info = {"n": n, "note": "", "fit_rms_mm": np.nan, "median_resid_mm": np.nan,
            "slope_deg": np.nan}
    if n < k:
        info["note"] = f"n={n} < {k} coefficients at order {surface_order}"
        return np.nan, info
    u = x[m] - px
    v = y[m] - py
    zz = z[m]
    A = _design(u, v, surface_order)
    try:
        coef, *_ = np.linalg.lstsq(A, zz, rcond=None)
    except np.linalg.LinAlgError:                                  # pragma: no cover
        info["note"] = "least-squares fit failed"
        return np.nan, info
    if not np.isfinite(coef).all():                                # pragma: no cover
        info["note"] = "singular design (collinear returns)"
        return np.nan, info
    resid = zz - A @ coef
    info["fit_rms_mm"] = 1000.0 * float(np.sqrt(np.mean(resid ** 2)))
    q = float(np.quantile(resid, quantile))
    info["median_resid_mm"] = 1000.0 * q
    info["slope_deg"] = float(np.degrees(np.arctan(np.hypot(coef[1], coef[2]))))
    return float(coef[0] + q), info


def radius_ladder(res: float) -> tuple:
    """The radii swept by default, all multiples of the pipeline's grid resolution.

    ``res/2`` is the cell half-width the slope-normal estimator reads at; ``1.5*res`` is
    the half-width of ``_poly2_ground``'s 3x3 window; the rest run out to ``5*res`` so
    the breakdown of the local fit is visible rather than cropped.
    """
    return (0.5 * res, 1.0 * res, 1.5 * res, 2.0 * res, 3.0 * res, 4.0 * res, 5.0 * res)


# ------------------------------------------------------------------------ ground crops

def vendor_ground_near(tile_path, easting, northing, half_width, *, ground_class=2):
    """Vendor-classified ground returns in a square window about a point.

    ``ground_class=2`` is ASPRS bare earth, the class
    :func:`lidar_diff_icp.pipeline.read_after_ground` uses in ``mode="class2"``.
    """
    import laspy
    f = laspy.read(str(tile_path))
    x = np.asarray(f.x); y = np.asarray(f.y)
    m = (np.abs(x - easting) <= half_width) & (np.abs(y - northing) <= half_width)
    n_in = int(m.sum())
    cl = np.asarray(f.classification)[m]
    g = cl == ground_class
    sa = scan_angle_deg(f)[m][g]
    return GroundReturns(x[m][g], y[m][g], np.asarray(f.z)[m][g], sa,
                         np.asarray(f.point_source_id)[m][g],
                         source=f"vendor_class{ground_class}",
                         origin=str(tile_path), n_input=n_in)


def csf_ground_near(tile_path, easting, northing, half_width, *, pdal=None, cache_dir=None,
                    **csf_kwargs):
    """CSF-classified ground returns in a square window about a point.

    This is the pipeline's own ground source (``ground_source="csf"`` in
    ``corrections.json``) via :func:`lidar_diff_icp.ground.classify_ground_csf`, so the
    tie is measured on the same bare earth the product is built from. The window is
    cropped first: CSF on a 600 m box of 2008 MN lidar is ~350 k points and ~4 s, against
    minutes for a whole tile.

    ``half_width`` must be large enough that the cloth is not dominated by the crop edge.
    The caller chooses it and it is reported; there is no hidden default here.
    """
    import laspy
    from ..ground import classify_ground_csf

    tmp = tempfile.mkdtemp(prefix="gt_csf_")
    try:
        f = laspy.read(str(tile_path))
        x = np.asarray(f.x); y = np.asarray(f.y)
        m = (np.abs(x - easting) <= half_width) & (np.abs(y - northing) <= half_width)
        n_in = int(m.sum())
        if n_in == 0:
            return GroundReturns(*[np.empty(0) for _ in range(3)], np.empty(0),
                                 np.empty(0, int), "csf", str(tile_path), 0)
        crop = os.path.join(tmp, "crop.las")
        sub = laspy.LasData(f.header)
        sub.points = f.points[m]
        sub.write(crop)
        del f, sub, x, y
        out = os.path.join(tmp, "ground.las")
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            out = os.path.join(cache_dir,
                               f"{os.path.basename(str(tile_path))}_"
                               f"{easting:.0f}_{northing:.0f}_{half_width:.0f}_csf.las")
        if not (cache_dir and os.path.exists(out)):
            classify_ground_csf(crop, out, pdal=pdal, **csf_kwargs)
        g = laspy.read(out)
        return GroundReturns(np.asarray(g.x), np.asarray(g.y), np.asarray(g.z),
                             scan_angle_deg(g),
                             np.asarray(g.point_source_id), source="csf",
                             origin=str(tile_path), n_input=n_in)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------------- the estimate

def geoid_shift_for(checkpoint, *, box_m=500.0, before_geoid="us_noaa_geoid03_conus.tif",
                    crs="EPSG:26915"):
    """Geoid shift to ADD to gen1 z so it lands on the checkpoint's geoid model.

    Delegates to :func:`lidar_diff_icp.references.geoid_difference`, which gathers
    ``N_before - N_after`` from the PROJ grids -- computed, never hard-coded. The
    checkpoint's own ``geoid_model`` selects the after-grid, so a checkpoint on a
    different model automatically gets a different shift. Raises for a model this
    function has no grid name for, rather than substituting one.
    """
    from ..references import geoid_difference

    checkpoint.require_datum()
    grids = {"GEOID18": "us_noaa_g2018u0.tif", "GEOID12B": "us_noaa_g2012bu0.tif",
             "GEOID12A": "us_noaa_g2012au0.tif", "GEOID09": "us_noaa_geoid09_conus.tif",
             "GEOID03": "us_noaa_geoid03_conus.tif"}
    key = checkpoint.geoid_model.strip().upper().replace(" ", "")
    if key not in grids:
        raise ValueError(
            f"checkpoint {checkpoint.point_id!r}: no PROJ grid known for geoid model "
            f"{checkpoint.geoid_model!r}; known: {sorted(grids)}. Add the grid name "
            "explicitly rather than letting the tie use the wrong model.")
    e, n = checkpoint.easting, checkpoint.northing
    a, b, c = geoid_difference((e - box_m, n - box_m, e + box_m, n + box_m), crs,
                               before_geoid=before_geoid, after_geoid=grids[key])
    return float(a)          # evaluated at the box centre = the checkpoint


def estimate_tie(checkpoint, ground: GroundReturns, *, line=None, res=5.0,
                 radii=None, report_radius=None, surface_order=2, quantile=0.50,
                 swath_shift_m=(0.0, 0.0, 0.0), geoid_shift_m=None,
                 geoid_before="us_noaa_geoid03_conus.tif", crs="EPSG:26915") -> TieEstimate:
    """Tie one checkpoint to a lidar cloud, over a ladder of radii.

    ``ground``      ground-classified returns near the mark (see :func:`csf_ground_near`
                    / :func:`vendor_ground_near`).
    ``line``        restrict to one ``point_source_id``. A checkpoint under a single
                    line's nadir gives an unmixed tie; mixing lines folds their relative
                    offsets into the answer.
    ``swath_shift_m`` (dx, dy, dz) applied to the returns first -- the chain-accumulated
                    alignment that brings ``line`` into the reference swath's frame.
    ``geoid_shift_m`` added to gen1 z. ``None`` computes it from the PROJ grids for this
                    checkpoint's geoid model (:func:`geoid_shift_for`).

    Everything that could bias the answer is returned in ``.params`` with its origin, so
    a caller can declare it into a ``trust.provenance.Run`` without retyping it.
    """
    cp = checkpoint.require_datum()
    surveyed = cp.elevation_m
    res = float(res)
    radii = tuple(radius_ladder(res)) if radii is None else tuple(float(r) for r in radii)
    report_radius = 1.5 * res if report_radius is None else float(report_radius)
    if report_radius not in radii:
        radii = tuple(sorted(set(radii + (report_radius,))))

    g = ground.for_line(line).shifted(*swath_shift_m)
    if geoid_shift_m is None:
        geoid_shift_m = geoid_shift_for(cp, before_geoid=geoid_before, crs=crs)
    geoid_shift_m = float(geoid_shift_m)

    curve = []
    for R in radii:
        zh, info = ground_elevation_at(g.x, g.y, g.z, cp.easting, cp.northing, R,
                                       surface_order=surface_order, quantile=quantile)
        sel = np.hypot(g.x - cp.easting, g.y - cp.northing) <= R
        zs = g.z[sel]
        sa = np.abs(g.scan_angle[sel])
        relief = 1000.0 * float(np.percentile(zs, 95) - np.percentile(zs, 5)) if zs.size else np.nan
        p05, p50, p95 = ((float(np.percentile(sa, 5)), float(np.percentile(sa, 50)),
                          float(np.percentile(sa, 95))) if sa.size else (np.nan,) * 3)
        curve.append(RadiusEstimate(
            radius_m=R, n=info["n"], z_lidar_m=zh, fit_rms_mm=info["fit_rms_mm"],
            median_resid_mm=info["median_resid_mm"], slope_deg=info["slope_deg"],
            relief_mm=relief, scan_p50=p50, scan_p05=p05, scan_p95=p95,
            n_lines=int(np.unique(g.point_source_id[sel]).size),
            ok=bool(np.isfinite(zh)), note=info["note"]))

    def tie_at(r):
        return 1000.0 * (surveyed - (r.z_lidar_m + geoid_shift_m))

    good = [r for r in curve if r.ok]
    # Pipeline-scale radii: the window sizes the slope-normal / poly2 estimators
    # actually work at, res/2 (cell half-width) through 2*res. The spread ACROSS those
    # is the honest uncertainty; the tail radii are kept in the curve to show the
    # local fit breaking down, not to widen the error bar.
    scale = [r for r in good if 0.5 * res <= r.radius_m <= 2.0 * res]
    ties_scale = [tie_at(r) for r in scale]
    ties_all = [tie_at(r) for r in good]
    spread = (max(ties_scale) - min(ties_scale)) if len(ties_scale) > 1 else np.nan
    spread_all = (max(ties_all) - min(ties_all)) if len(ties_all) > 1 else np.nan

    tie_median = float(np.median(ties_scale)) if ties_scale else np.nan
    rep = next((r for r in curve if r.radius_m == report_radius), None)
    z_raw = rep.z_lidar_m if (rep and rep.ok) else np.nan
    tie = 1000.0 * (surveyed - (z_raw + geoid_shift_m))
    se = (rep.fit_rms_mm / np.sqrt(rep.n) if (rep and rep.ok and rep.n) else np.nan)

    params = [
        Param("res_m", res, "repo", _SRC_RES),
        Param("ground_quantile", quantile, "repo", _SRC_Q),
        Param("surface_order", surface_order, "repo", _SRC_ORDER),
        Param("report_radius_m", report_radius, "repo", _SRC_WINDOW),
        Param("radius_ladder_m", list(radii), "repo",
              "multiples of res: res/2 (cell half-width) .. 5*res, so the local fit's "
              "breakdown is shown rather than cropped"),
        Param("ground_source", g.source, "repo", _SRC_GROUND),
        Param("flight_line", line, "repo",
              "the single line covering this mark; mixing lines would fold their "
              "relative offsets into the tie"),
        Param("swath_shift_m", tuple(round(float(v), 4) for v in swath_shift_m), "repo",
              "chain-accumulated coreg.align_swaths correction into the reference frame"),
        Param("geoid_shift_m", round(geoid_shift_m, 5), "repo",
              f"references.geoid_difference, {geoid_before} -> {cp.geoid_model}, "
              "computed from the PROJ grids at the checkpoint"),
    ]
    notes = []
    if rep is not None and rep.n_lines > 1:
        notes.append(f"{rep.n_lines} flight lines inside the report radius -- the tie mixes them")
    if cp.point_type.upper() == "VVA":
        notes.append("VVA checkpoint (under vegetation): the published 3DEP spread for this "
                     "class is 27 cm at the 95th percentile, against 3.5 cm RMSE for NVA")
    return TieEstimate(
        point_id=cp.point_id, point_type=cp.point_type, line=line, curve=curve,
        report_radius_m=report_radius, z_lidar_raw_m=z_raw,
        swath_shift_m=tuple(float(v) for v in swath_shift_m), geoid_shift_m=geoid_shift_m,
        surveyed_m=surveyed, tie_mm=tie, tie_median_mm=tie_median, radius_spread_mm=spread,
        radius_spread_all_mm=spread_all, fit_se_mm=se,
        sigma_mm=(0.5 * spread if np.isfinite(spread) else np.nan),
        params=params, notes=notes)
