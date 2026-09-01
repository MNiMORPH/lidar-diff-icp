"""q2(cover): the gen2 percentile whose ground elevation matches gen1's median.

WHAT IT IS. gen1 and gen2 do not see the same ground under vegetation, so differencing
their medians differences two different surfaces. This finds, per canopy-cover bin, the
percentile of gen2's near-ground column that lands on gen1's median, and fits how that
percentile declines with cover:

    q2(c) = 0.50 + b * c

**It is a matching relation, not a claim about which epoch is right.** It brings the two
epochs onto a common surface for differencing.

THE RECIPE IS FIXED AND DOCUMENTED, and the defaults here ARE it -- see
``analysis/ridgelines/Q2_COVER_RELATION.md``:

  * ``q2(0) = 0.50`` is IMPOSED, not fitted: at zero cover both epochs see the true ground.
    gen1's own percentile ``q1`` is likewise fixed at 0.50.
  * uniform 0.05 cover bins, so no bin gains leverage from where an edge was drawn.
  * weight by the cells behind each point, and keep EVERY bin -- a minimum count would cut
    exactly the sparse high-cover regime the relation is about.
  * per cell, >= 5 gen1 returns and >= 10 gen2 class-2 returns.
  * linear in cover. The doc's variant table shows the linear term is stable across binning
    and weighting while the power exponent is not: "the curvature is an artefact of binning
    and weighting choices."

THE FIT IS PER SITE. It depends on each pair's phenology, so there is no site-invariant
slope to carry between tiles. Measured with the floodplain excluded: elba_fulldensity
-0.1247, elbaext -0.1871 -- and the tiles differ because their gen2 undergrowth differs,
which the bare-ground bins show directly (q2* 0.572 vs 0.469).

THE VALLEY IS EXCLUDED BY DEFAULT, via ``refcells.reference_cells``, which applies
``floodplain_mask.npy``. Floodplain aggradation is possible, so valley cells are not valid
reference ground. ``exclude_valley=False`` restores the pre-2026-08-26-16:49 population and
is needed to reconstruct the originally published -0.1922.

A MATCHING PERCENTILE NEED NOT EXIST. If gen1's median lies outside gen2's whole
near-ground column for a bin, that bin is UNMATCHABLE and is reported, not silently
dropped and not crashed on. At elba with the valley excluded, bin 0.70-0.75 (n=15) is
unmatchable: gen1's median sits below gen2's entire column.
"""
from __future__ import annotations

import os
import numpy as np
from scipy.optimize import brentq

from .refcells import reference_cells

# The documented recipe. Changing any of these changes the relation, not just its precision.
BIN_WIDTH = 0.05
MIN_GEN1_RETURNS = 5
MIN_GEN2_RETURNS = 10
Q1 = 0.50            # gen1's own percentile: the median
Q2_AT_ZERO = 0.50    # imposed, not fitted


def ragged_sorted(cell, val, ncell):
    """Flat array of per-cell sorted values, with per-cell offsets and counts."""
    o = np.lexsort((val, cell))
    cs = cell[o]; vs = val[o]
    n = np.bincount(cs, minlength=ncell)
    off = np.r_[0, np.cumsum(n)[:-1]]
    return vs, off, n


def ragged_quantile(vs, off, n, q, sel):
    """Linear-interpolated quantile (numpy 'linear' convention) for the selected cells."""
    nn = n[sel].astype(float)
    pos = q * (nn - 1.0)
    lo = np.floor(pos).astype(np.int64); hi = np.minimum(lo + 1, nn.astype(np.int64) - 1)
    f = pos - lo
    b = off[sel]
    return vs[b + lo] * (1 - f) + vs[b + hi] * f


def hist_quantile(C, ntot, q, zlo, dz):
    """Quantile of a histogram row, interpolated WITHIN the bin, in mm."""
    r = q * ntot
    k = (C >= r[:, None]).argmax(1)
    below = np.where(k > 0, C[np.arange(C.shape[0]), np.maximum(k - 1, 0)], 0.0)
    inbin = C[np.arange(C.shape[0]), k] - below
    f = np.where(inbin > 0, (r - below) / np.maximum(inbin, 1e-9), 0.0)
    return (zlo + (k + np.clip(f, 0, 1)) * dz) * 1000.0


def q2_at(cover, slope, q2_zero=Q2_AT_ZERO):
    """The gen2 percentile to sample at a given canopy cover.

    Clipped to (0, 1): a percentile outside that is not a percentile, and the relation is
    linear so it will leave the range if extrapolated far enough.
    """
    return np.clip(q2_zero + slope * np.asarray(cover, float), 1e-4, 1 - 1e-4)


def solve_q2_per_bin(cover, gen1_median_mm, C, ntot, zlo, dz, bin_width=BIN_WIDTH):
    """Per cover bin, the gen2 percentile matching gen1's median.

    Returns (bins, unmatchable). Each bin carries its centre, cell count, solved q2*, and
    the mm-per-0.01-rank sensitivity. `unmatchable` lists bins where no percentile works,
    with the side gen1's median fell on -- these are a property of the data, not an error.
    """
    edges = np.arange(0.0, float(np.max(cover)) + bin_width, bin_width)
    bins, unmatchable = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (cover > lo - 1e-9) & (cover <= hi)
        if not m.any():
            continue
        f = lambda q: float(np.median(gen1_median_mm[m] - hist_quantile(C[m], ntot[m], q, zlo, dz)))
        fa, fb = f(1e-4), f(1 - 1e-4)
        if fa * fb > 0:
            unmatchable.append(dict(lo=float(lo), hi=float(hi), n=int(m.sum()),
                                    side="above" if fa > 0 else "below",
                                    resid_lo_mm=fa, resid_hi_mm=fb))
            continue
        q = brentq(f, 1e-4, 1 - 1e-4, xtol=1e-6)
        sens = (f(max(q - 0.05, 0.01)) - f(min(q + 0.05, 0.99))) / 10.0
        bins.append(dict(lo=float(lo), hi=float(hi), n=int(m.sum()),
                         cover=float(np.mean(cover[m])), q2=float(q),
                         mm_per_centirank=float(sens)))
    return bins, unmatchable


def fit_slope(bins, q2_zero=Q2_AT_ZERO):
    """Weighted least squares of q2 on cover, THROUGH the imposed anchor q2(0) = q2_zero.

    Weights are cell counts, per the documented recipe. Returns (slope, se).
    """
    c = np.array([b["cover"] for b in bins], float)
    q = np.array([b["q2"] for b in bins], float)
    w = np.array([b["n"] for b in bins], float)
    if c.size < 2:
        return float("nan"), float("nan")
    slope = float((w * c * (q - q2_zero)).sum() / (w * c * c).sum())
    resid = q - (q2_zero + slope * c)
    dof = max(c.size - 1, 1)
    se = float(np.sqrt((w * resid ** 2).sum() / dof / (w * c * c).sum()))
    return slope, se


def free_intercept(bins):
    """Unconstrained linear fit, as a CHECK on the imposed anchor -- never as the product.

    Returns (a, b, se_a, se_b). If `a` departs from 0.50 the anchor is doing work the data
    do not support, which is worth knowing even though the anchor is theoretical.
    """
    c = np.array([b["cover"] for b in bins], float)
    q = np.array([b["q2"] for b in bins], float)
    w = np.array([b["n"] for b in bins], float)
    if c.size < 3:
        return (float("nan"),) * 4
    A = np.c_[np.ones(c.size), c]
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    beta = cov @ (A.T @ W @ q)
    r = q - A @ beta
    s2 = float((w * r ** 2).sum() / (c.size - 2))
    err = np.sqrt(np.diag(cov * s2))
    return float(beta[0]), float(beta[1]), float(err[0]), float(err[1])


def fit_tile(tile_dir, *, exclude_valley=True, valley_top_m=None,
             use_floodplain_mask=True, bin_width=BIN_WIDTH,
             min_gen1=MIN_GEN1_RETURNS, min_gen2=MIN_GEN2_RETURNS, slope_max=90.0):
    """Fit q2(cover) for one tile, from that tile's own derived products.

    Needs, in `tile_dir`: nearground_cells_sn.npz, z_after.npy, canopy_cover_pfs.npy,
    beam_offset_table.parquet, nearground_gen2_class_split.npz.

    Returns a dict with the slope and its SE, the per-bin table, the unmatchable bins, the
    free-intercept check, and the cell accounting -- so a caller cannot quote the slope
    without the population it came from.
    """
    import pyarrow.parquet as pq

    D = str(tile_dir).rstrip("/")
    need = ["nearground_cells_sn.npz", "z_after.npy", "canopy_cover_pfs.npy",
            "beam_offset_table.parquet", "nearground_gen2_class_split.npz"]
    missing = [f for f in need if not os.path.exists(os.path.join(D, f))]
    if missing:
        raise FileNotFoundError(
            f"{D} is missing {missing}. Each has a --tile-parameterized producer; "
            f"nearground_gen2_class_split.npz comes from "
            f"analysis/ridgelines/nearground_class_split.py.")

    cube = np.load(f"{D}/nearground_cells_sn.npz")
    cells = cube["cells"]; dz = float(cube["dz"]); zlo = float(cube["zlo"])
    N = np.load(f"{D}/z_after.npy").size
    cover_all = np.load(f"{D}/canopy_cover_pfs.npy").ravel()[cells]

    t = pq.read_table(f"{D}/beam_offset_table.parquet",
                      columns=["cell", "d_mm_corr", "in_grid"])
    g = t["in_grid"].to_numpy().astype(bool)
    vs, off, n1 = ragged_sorted(t["cell"].to_numpy()[g],
                                t["d_mm_corr"].to_numpy()[g].astype(float), N)

    Hg = np.load(f"{D}/nearground_gen2_class_split.npz")["Hg"]
    Cg = np.cumsum(Hg, 1).astype(float); ng = Cg[:, -1]

    stable, cuts = reference_cells(D, cells=cells, slope_max=slope_max,
                                   exclude_valley=exclude_valley,
                                   valley_top_m=valley_top_m,
                                   use_floodplain_mask=use_floodplain_mask)
    ok = (stable & (n1[cells] >= min_gen1) & (ng >= min_gen2) & np.isfinite(cover_all))
    sel = cells[ok]
    g1 = ragged_quantile(vs, off, n1, Q1, sel)

    bins, unmatchable = solve_q2_per_bin(cover_all[ok], g1, Cg[ok], ng[ok], zlo, dz,
                                         bin_width=bin_width)
    slope, se = fit_slope(bins)
    a, b, sa, sb = free_intercept(bins)
    return dict(
        tile=os.path.basename(D), slope=slope, slope_se=se,
        relation=f"q2 = {Q2_AT_ZERO} + {slope:.4f} * cover",
        q2_at_zero_imposed=Q2_AT_ZERO, q1=Q1,
        bins=bins, unmatchable=unmatchable,
        free_intercept=dict(a=a, b=b, se_a=sa, se_b=sb,
                            sigma_from_imposed=(a - Q2_AT_ZERO) / sa if sa == sa else float("nan")),
        cuts=cuts,
        cells=dict(stable=int(stable.sum()), used=int(ok.sum()),
                   exclude_valley=bool(exclude_valley), valley_top_m=valley_top_m,
                   use_floodplain_mask=bool(use_floodplain_mask),
                   min_gen1=min_gen1, min_gen2=min_gen2, bin_width=bin_width),
    )
