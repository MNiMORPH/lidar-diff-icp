"""Predict a lidar epoch-difference offset jointly from terrain SLOPE and CANOPY COVER,
keeping the covariance between those two predictors visible instead of assumed away.

Why this module exists: on dissected terrain the steep ground is preferentially forested
(at Elba, ``corr(slope, cover) ~ +0.57``), so a model fitted on one predictor silently
absorbs the other. Every routine here is built so that confounding stays in view:

* :func:`predictor_covariance` reports the correlation and the variance-inflation factor.
* :func:`median_surface` bins the response on a (slope, cover) grid and returns the CELL
  COUNTS with it, so combinations the terrain never supplies are reported as unsupported
  rather than quietly extrapolated through.
* :func:`matched_band_effects` holds one predictor inside a narrow band and fits the
  other -- the direct control, free of any functional-form assumption.
* :func:`partial_correlations` gives each predictor's association net of the other, in
  both Pearson and Spearman (rank) form, the latter robust to the heavy per-cell tails.
* :func:`fit_offset_model` fits additive and interaction forms and reports both, because
  on this data the interaction is not a refinement: slope and cover push the offset in
  OPPOSITE directions, and an additive model averages that structure away.

Two response scales matter and must not be confused. Per-RETURN or per-CELL values carry
scatter roughly two orders of magnitude larger than the systematic term, so a per-cell fit
returns R^2 ~ 1e-3 even when the systematic structure is real and reproducible; fitting
the binned MEDIAN surface (with ``weights=counts``) is what answers "predict the median
offset". :func:`fit_offset_model` is agnostic -- pass whichever you mean, and read the
returned ``r2`` against that scale.

Sign convention is inherited from the caller's response: with ``d_mm`` as gen1-minus-gen2
slope-normal offset, negative means the older epoch reads the ground LOW.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["OffsetModel", "predictor_covariance", "median_surface", "matched_band_effects",
           "partial_correlations", "fit_offset_model", "cell_reduce"]


@dataclass
class OffsetModel:
    """A fitted offset model, in PHYSICAL units (slope in degrees, cover as a fraction).

    ``coeffs`` is ``(intercept, slope, cover)`` for an additive fit and
    ``(intercept, slope, cover, slope*cover)`` when ``interaction`` is True.
    """
    coeffs: np.ndarray
    interaction: bool
    r2: float
    rms: float
    n: int
    weighted: bool = False
    meta: dict = field(default_factory=dict)

    def predict(self, slope, cover):
        """Predicted offset at ``slope`` (deg) and ``cover`` (fraction)."""
        slope = np.asarray(slope, float); cover = np.asarray(cover, float)
        b = self.coeffs
        out = b[0] + b[1] * slope + b[2] * cover
        if self.interaction:
            out = out + b[3] * slope * cover
        return out

    def d_dcover(self, slope):
        """Sensitivity to cover (per full cover unit) at a given slope."""
        return self.coeffs[2] + (self.coeffs[3] * np.asarray(slope, float) if self.interaction else 0.0)

    def d_dslope(self, cover):
        """Sensitivity to slope (per degree) at a given cover."""
        return self.coeffs[1] + (self.coeffs[3] * np.asarray(cover, float) if self.interaction else 0.0)

    def __str__(self):
        b = self.coeffs
        s = f"d = {b[0]:+.1f} {b[1]:+.3f}*slope {b[2]:+.1f}*cover"
        if self.interaction:
            s += f" {b[3]:+.3f}*slope*cover"
        return s + f"   (R2 {self.r2:.3f}, RMS {self.rms:.1f}, n {self.n})"


def cell_reduce(cell, value, min_n=3, reducer=np.median):
    """Reduce a per-return ``value`` to one number per ``cell`` (default: the median).

    Returns ``(cells, reduced, counts)`` for cells with at least ``min_n`` returns.
    Per-cell medians are the natural response here: they difference out per-return noise
    while preserving the per-cell systematic term.
    """
    cell = np.asarray(cell); value = np.asarray(value, float)
    order = np.argsort(cell, kind="stable")
    c_sorted = cell[order]; v_sorted = value[order]
    cells, start, counts = np.unique(c_sorted, return_index=True, return_counts=True)
    keep = counts >= min_n
    out = np.array([reducer(v_sorted[s:s + n]) for s, n in zip(start[keep], counts[keep])])
    return cells[keep], out, counts[keep]


def predictor_covariance(slope, cover):
    """``(corr, vif)`` for the two predictors. ``vif = 1`` means orthogonal."""
    r = float(np.corrcoef(np.asarray(slope, float), np.asarray(cover, float))[0, 1])
    return r, float(1.0 / (1.0 - r ** 2))


def _lstsq(X, y, w=None):
    if w is None:
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return b
    rw = np.sqrt(np.asarray(w, float))[:, None]
    b, *_ = np.linalg.lstsq(X * rw, y * rw[:, 0], rcond=None)
    return b


def _design(slope, cover, interaction):
    cols = [np.ones(len(slope)), slope, cover]
    if interaction:
        cols.append(slope * cover)
    return np.column_stack(cols)


def fit_offset_model(slope, cover, value, *, weights=None, interaction=True):
    """Least-squares fit of ``value`` on slope and cover, in physical units.

    Pass ``weights`` (e.g. the cell counts behind each binned median) to fit a median
    surface; leave it None for a per-cell fit. Returns an :class:`OffsetModel`.
    """
    slope = np.asarray(slope, float); cover = np.asarray(cover, float)
    value = np.asarray(value, float)
    X = _design(slope, cover, interaction)
    b = _lstsq(X, value, weights)
    pred = X @ b
    if weights is None:
        w = np.ones_like(value)
    else:
        w = np.asarray(weights, float)
    mean = np.sum(w * value) / np.sum(w)
    ss_res = np.sum(w * (value - pred) ** 2); ss_tot = np.sum(w * (value - mean) ** 2)
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    rms = float(np.sqrt(ss_res / np.sum(w)))
    return OffsetModel(coeffs=b, interaction=interaction, r2=r2, rms=rms, n=int(value.size),
                       weighted=weights is not None)


def median_surface(slope, cover, value, slope_edges, cover_edges, min_cells=30):
    """Bin ``value`` on the (slope, cover) grid; return ``(grid, counts)``.

    ``grid`` holds the median of ``value`` per box and is NaN wherever fewer than
    ``min_cells`` samples fall in the box -- those are combinations the terrain does not
    supply, and NaN is the honest answer there rather than an extrapolated number.
    """
    slope = np.asarray(slope, float); cover = np.asarray(cover, float)
    value = np.asarray(value, float)
    slope_edges = np.asarray(slope_edges, float); cover_edges = np.asarray(cover_edges, float)
    ns, nc = len(slope_edges) - 1, len(cover_edges) - 1
    si = np.clip(np.digitize(slope, slope_edges) - 1, 0, ns - 1)
    ci = np.clip(np.digitize(cover, cover_edges) - 1, 0, nc - 1)
    grid = np.full((ns, nc), np.nan); counts = np.zeros((ns, nc), int)
    for i in range(ns):
        mi = si == i
        for j in range(nc):
            m = mi & (ci == j)
            counts[i, j] = int(m.sum())
            if counts[i, j] >= min_cells:
                grid[i, j] = float(np.median(value[m]))
    return grid, counts


def surface_centres(slope_edges, cover_edges):
    """Box-centre (slope, cover) meshes matching :func:`median_surface` output."""
    s = 0.5 * (np.asarray(slope_edges, float)[:-1] + np.asarray(slope_edges, float)[1:])
    c = 0.5 * (np.asarray(cover_edges, float)[:-1] + np.asarray(cover_edges, float)[1:])
    return np.meshgrid(s, c, indexing="ij")


def matched_band_effects(held, varied, value, edges, min_n=200):
    """Fit ``value`` on ``varied`` inside each band of ``held``: the covariance control.

    Holding one predictor inside a narrow band and regressing on the other removes the
    confounding by construction, with no functional form imposed on the held variable.
    Returns a list of ``(lo, hi, n, d value/d varied, varied_min, varied_max)``; the
    gradient is NaN for bands with fewer than ``min_n`` samples.
    """
    held = np.asarray(held, float); varied = np.asarray(varied, float)
    value = np.asarray(value, float); edges = np.asarray(edges, float)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (held >= lo) & (held < hi)
        n = int(m.sum())
        if n < min_n:
            rows.append((float(lo), float(hi), n, float("nan"), float("nan"), float("nan")))
            continue
        b = _lstsq(np.column_stack([np.ones(n), varied[m]]), value[m])
        rows.append((float(lo), float(hi), n, float(b[1]),
                     float(varied[m].min()), float(varied[m].max())))
    return rows


def _resid(y, x):
    X = np.column_stack([np.ones(len(x)), x])
    return y - X @ _lstsq(X, y)


def partial_correlations(slope, cover, value):
    """Each predictor's correlation with ``value`` NET of the other.

    Returns ``{"cover|slope": (pearson, spearman), "slope|cover": (pearson, spearman)}``.
    The Spearman form ranks first, which matters because per-cell offsets are heavy-tailed.
    """
    from scipy.stats import rankdata
    slope = np.asarray(slope, float); cover = np.asarray(cover, float)
    value = np.asarray(value, float)
    out = {}
    for name, x1, x2 in (("cover|slope", cover, slope), ("slope|cover", slope, cover)):
        rp = float(np.corrcoef(_resid(value, x2), _resid(x1, x2))[0, 1])
        rs = float(np.corrcoef(_resid(rankdata(value), rankdata(x2)),
                               _resid(rankdata(x1), rankdata(x2)))[0, 1])
        out[name] = (rp, rs)
    return out
