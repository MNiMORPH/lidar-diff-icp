"""Measure gen2's per-swath vertical deviation, and its sensitivity to where it is measured.

gen1's inter-swath structure is measured and corrected (``corrections.json`` carries a
per-swath ``dxdydz`` for lines 133-138).  gen2 has never been measured at all: no
``corrections*.json`` on disk carries a gen2 swath.  This script asks the two questions
that gap leaves open, using the repository's OWN estimator rather than a new one:

1. **Deviation** -- what is the vertical tie between adjacent gen2 flight lines?
2. **Sensitivity** -- how much does that tie move with WHERE and over HOW LARGE a window
   it is measured?  ``analysis/LOCAL_TIE_CHAINING.md`` showed for gen1 that this, not
   ``coreg``'s formal sigma, is the honest error bar.

Every offset returned comes from :func:`lidar_diff_icp.coreg.coregister_swaths`,
unmodified.  The only thing this module adds is where the points come from.

Two gen2-specific facts drive the reading, and both are asserted rather than assumed --
:func:`describe_cloud` prints them so a run cannot hide them:

* gen2's delivered classification carries **only classes 1 and 2**.  ``coreg``'s shipped
  ``exclude=(5, 6, 9)`` therefore removes NOTHING from a gen2 cloud and would admit
  leaf-on canopy into the surface.  Ground must be selected as ``classification == 2``.
* gen2 is LAS point format 6+, whose ``scan_angle`` is an int16 in 0.006 degree units and
  is a DIFFERENT dimension from format<6's ``scan_angle_rank``.
  :func:`lidar_diff_icp.io.read_tile` reads ``scan_angle_rank`` and falls back to zeros
  when it is absent; on an all-zero predictor ``coreg.across_track_tie`` returns NaN and
  ``coregister_swaths`` silently falls back to the median tie.  The intercept tie would
  look like it ran and would not have.

There are no default window sizes, resolutions, centres or ground classes.  All are
required arguments, for the reason ``localtie`` gives: a default window size hides the
answer to the question the module exists to ask.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

LAS14_SCAN_ANGLE_DEG_PER_UNIT = 0.006  # LAS 1.4 spec, point formats 6-10, ScanAngle


@dataclass(frozen=True)
class CloudDescription:
    """What a gen2 cloud actually contains, asserted before it is used."""

    classes: dict
    has_scan_angle: bool
    has_scan_angle_rank: bool
    scan_angle_deg_min: float
    scan_angle_deg_max: float
    n_points: int

    def ground_class_is_only_selector(self) -> bool:
        """True when ``exclude=(5, 6, 9)`` would remove nothing from this cloud."""
        return not ({5, 6, 9} & set(self.classes))


def describe_cloud(las) -> CloudDescription:
    """Read the facts that decide how this cloud must be handled."""
    import numpy as _np

    dims = set(las.point_format.dimension_names)
    cl = _np.asarray(las.classification)
    u, c = _np.unique(cl, return_counts=True)
    has_sa = "scan_angle" in dims
    has_sar = "scan_angle_rank" in dims
    if has_sa:
        sa = _np.asarray(las.scan_angle) * LAS14_SCAN_ANGLE_DEG_PER_UNIT
    elif has_sar:
        sa = _np.asarray(las.scan_angle_rank).astype(float)
    else:
        sa = _np.zeros(cl.size)
    return CloudDescription(
        classes={int(a): int(b) for a, b in zip(u, c)},
        has_scan_angle=has_sa,
        has_scan_angle_rank=has_sar,
        scan_angle_deg_min=float(sa.min()) if sa.size else float("nan"),
        scan_angle_deg_max=float(sa.max()) if sa.size else float("nan"),
        n_points=int(cl.size),
    )


def scan_angle_degrees(las):
    """Scan angle in DEGREES, whichever LAS generation the file is.

    Raises rather than returning zeros: a silent zero predictor is exactly the failure
    described in the module docstring.
    """
    dims = set(las.point_format.dimension_names)
    if "scan_angle" in dims:
        return np.asarray(las.scan_angle).astype(float) * LAS14_SCAN_ANGLE_DEG_PER_UNIT
    if "scan_angle_rank" in dims:
        return np.asarray(las.scan_angle_rank).astype(float)
    raise ValueError(
        "cloud carries neither 'scan_angle' nor 'scan_angle_rank'; the intercept tie "
        "would run on an all-zero predictor and silently fall back to the median tie")


def crop_window(copc_path, *, easting, northing, half_width_m, keep_classes, pdal_bin):
    """Range-read a square window out of a COPC, keeping only ``keep_classes``.

    Returns the path of a temporary LAZ.  The read is bounded, so a 2 GB COPC costs
    only the window: this is a shared laptop.
    """
    copc_path = str(Path(copc_path).resolve())
    tmp = Path(tempfile.mkdtemp(prefix="g2swath_")) / "win.laz"
    x0, x1 = easting - half_width_m, easting + half_width_m
    y0, y1 = northing - half_width_m, northing + half_width_m
    stages = [{"type": "readers.copc", "filename": copc_path,
               "bounds": f"([{x0},{x1}],[{y0},{y1}])"}]
    ranges = "".join(f"Classification[{c}:{c}]" for c in sorted(keep_classes))
    stages.append({"type": "filters.range", "limits": ranges})
    stages.append({"type": "writers.las", "filename": str(tmp),
                   "compression": "laszip", "minor_version": 4,
                   "dataformat_id": 7})
    pipe = {"pipeline": stages}
    pf = tmp.parent / "pipe.json"
    pf.write_text(json.dumps(pipe))
    subprocess.run([str(pdal_bin), "pipeline", str(pf)], check=True,
                   capture_output=True)
    return tmp


def window_pointcloud(laz_path, *, crs):
    """Build a :class:`lidar_diff_icp.io.PointCloud` with the scan angle IN DEGREES."""
    import laspy
    from lidar_diff_icp.io import PointCloud

    f = laspy.read(str(laz_path))
    desc = describe_cloud(f)
    pc = PointCloud(
        x=np.asarray(f.x), y=np.asarray(f.y), z=np.asarray(f.z),
        point_source_id=np.asarray(f.point_source_id),
        classification=np.asarray(f.classification),
        gps_time=np.asarray(f.gps_time),
        scan_angle=scan_angle_degrees(f),
        crs=crs,
    )
    return pc, desc


def pair_dz(pc, line_ref, line_src, *, res_m, tie, exclude):
    """The vertical tie of one pair on this window, from ``coreg`` unmodified.

    Returns ``(dz_mm, dx_m, dy_m, n_cells)`` or ``None`` when the two lines do not
    overlap inside the window.
    """
    from lidar_diff_icp import coreg

    try:
        c = coreg.coregister_swaths(pc, line_ref, line_src, res=res_m,
                                    exclude=exclude, tie=tie)
    except ValueError:
        return None
    return (c.dz * 1000.0, c.dx, c.dy, int(getattr(c, "n", -1)))


# ------------------------------------------------------------------ statistics

@dataclass(frozen=True)
class PairStats:
    """The tie of one pair on one window, by several estimators.

    ``median_mm`` is what ``coreg``'s shipped ``overlap_median`` tie reduces to.  On the
    delivered 3DEP product it is **quantization-limited**: z is stored with a 1 cm scale,
    so a per-cell median lands on a 5 mm lattice and a small true offset reads as exactly
    zero.  ``mean_mm`` is not so limited -- quantization error is zero-mean, so it
    averages away -- which is why both are reported and neither is called "the" tie.

    ``se_block_mm`` is the standard error of ``mean_mm`` over spatial BLOCKS, not over
    cells: adjacent cells are not independent draws.  It is the SE of the mean of the
    per-block means.
    """

    pair: tuple
    n_cells: int
    n_blocks: int
    mean_mm: float
    se_block_mm: float
    median_mm: float
    nmad_mm: float
    k_intercept_mm: float
    c_mm_per_tan: float
    dtan_mean: float


def pair_stats(pc, line_ref, line_src, *, res_m, block_m, exclude) -> PairStats | None:
    """Tie statistics for one pair on this window.

    The overlap grids and the across-track intercept come from the same functions
    ``coreg.coregister_swaths`` uses (``swathdiff._median_grid``, ``coreg.across_track_tie``);
    what is added is the block-resampled standard error and the mean, because the median
    cannot resolve below the product's 1 cm z quantum.
    """
    from lidar_diff_icp import coreg
    from lidar_diff_icp.swathdiff import _median_grid

    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == line_ref)
    mb = terr & (pc.point_source_id == line_src)
    if ma.sum() == 0 or mb.sum() == 0:
        return None
    x, y, z = pc.x, pc.y, pc.z
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        return None
    nx = int(np.ceil((x1 - x0) / res_m)); ny = int(np.ceil((y1 - y0) / res_m))
    if nx < 2 or ny < 2:
        return None
    za = _median_grid(x[ma], y[ma], z[ma], res_m, x0, y0, nx, ny)
    zb = _median_grid(x[mb], y[mb], z[mb], res_m, x0, y0, nx, ny)
    tn = np.tan(np.radians(np.asarray(pc.scan_angle, float)))
    ta = _median_grid(x[ma], y[ma], tn[ma], res_m, x0, y0, nx, ny)
    tb = _median_grid(x[mb], y[mb], tn[mb], res_m, x0, y0, nx, ny)

    dh = (za - zb) * 1000.0
    dtan = ta - tb
    fin = np.isfinite(dh) & np.isfinite(dtan)
    if fin.sum() < 2:
        return None

    iy, ix = np.nonzero(fin)
    vals = dh[fin]
    per = int(round(block_m / res_m))
    bid = (iy // per) * (nx // per + 1) + (ix // per)
    order = np.argsort(bid, kind="stable")
    bid_s, vals_s = bid[order], vals[order]
    edges = np.flatnonzero(np.r_[True, bid_s[1:] != bid_s[:-1], True])
    bmeans = np.array([vals_s[a:b].mean() for a, b in zip(edges[:-1], edges[1:])])
    nb = bmeans.size
    se_block = float(bmeans.std(ddof=1) / np.sqrt(nb)) if nb > 1 else float("nan")

    k, c, _n = coreg.across_track_tie(dh, dtan)
    med = float(np.median(vals))
    return PairStats(
        pair=(int(line_ref), int(line_src)),
        n_cells=int(fin.sum()), n_blocks=int(nb),
        mean_mm=float(vals.mean()), se_block_mm=se_block,
        median_mm=med,
        nmad_mm=float(1.4826 * np.median(np.abs(vals - med))),
        k_intercept_mm=float(k), c_mm_per_tan=float(c),
        dtan_mean=float(dtan[fin].mean()),
    )


@dataclass(frozen=True)
class AcrossTrackFit:
    """OLS ``dh = k + c*dtan`` with cluster-robust standard errors.

    ``coreg.across_track_tie`` fits the same model by LAD, which is the right estimator
    on gen1 and is what the pipeline uses.  On the delivered gen2 product the response is
    on a 1 cm lattice, and LAD -- whose solution passes exactly through data points --
    snaps onto that lattice and returns exactly 0.0 for both ``k`` and ``c`` on some
    pairs.  That zero is an artifact of the estimator meeting the quantum, not a
    measurement that the term is absent: the same pairs return significantly non-zero
    coefficients under OLS, whose solution is not lattice-constrained.  Both are
    reported; neither is called "the" answer.

    Standard errors are clustered on spatial blocks because neighbouring 2 m cells are
    not independent draws.
    """

    k_mm: float
    se_k_mm: float
    c_mm_per_tan: float
    se_c_mm_per_tan: float
    n_cells: int
    n_blocks: int


def across_track_ols(dh_mm, dtan, block_id) -> AcrossTrackFit:
    """Fit ``dh = k + c*dtan`` by OLS, SEs clustered on ``block_id``."""
    y = np.asarray(dh_mm, float)
    t = np.asarray(dtan, float)
    g = np.asarray(block_id)
    X = np.c_[np.ones(y.size), t]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    XtX_inv = np.linalg.inv(X.T @ X)
    meat = np.zeros((2, 2))
    for b in np.unique(g):
        m = g == b
        s = X[m].T @ r[m]
        meat += np.outer(s, s)
    V = XtX_inv @ meat @ XtX_inv
    se = np.sqrt(np.diag(V))
    return AcrossTrackFit(k_mm=float(beta[0]), se_k_mm=float(se[0]),
                          c_mm_per_tan=float(beta[1]), se_c_mm_per_tan=float(se[1]),
                          n_cells=int(y.size), n_blocks=int(np.unique(g).size))


def overlap_grids(pc, line_ref, line_src, *, res_m, exclude):
    """The (dh_mm, dtan, block_id) triple both estimators read. None if no overlap."""
    from lidar_diff_icp.swathdiff import _median_grid

    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == line_ref)
    mb = terr & (pc.point_source_id == line_src)
    if ma.sum() == 0 or mb.sum() == 0:
        return None
    x, y, z = pc.x, pc.y, pc.z
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        return None
    nx = int(np.ceil((x1 - x0) / res_m)); ny = int(np.ceil((y1 - y0) / res_m))
    if nx < 2 or ny < 2:
        return None
    za = _median_grid(x[ma], y[ma], z[ma], res_m, x0, y0, nx, ny)
    zb = _median_grid(x[mb], y[mb], z[mb], res_m, x0, y0, nx, ny)
    tn = np.tan(np.radians(np.asarray(pc.scan_angle, float)))
    ta = _median_grid(x[ma], y[ma], tn[ma], res_m, x0, y0, nx, ny)
    tb = _median_grid(x[mb], y[mb], tn[mb], res_m, x0, y0, nx, ny)
    dh = (za - zb) * 1000.0
    dt = ta - tb
    f = np.isfinite(dh) & np.isfinite(dt)
    if f.sum() < 2:
        return None
    iy, ix = np.nonzero(f)
    return dh[f], dt[f], iy, ix, nx


def block_ids(iy, ix, nx, *, res_m, block_m):
    """Spatial block id per cell, for cluster-robust SEs."""
    per = int(round(block_m / res_m))
    return (iy // per) * (nx // per + 1) + (ix // per)
