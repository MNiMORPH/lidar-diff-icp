"""One way to bin, summarise and WEIGHT spatially correlated per-return data.

Every analysis here reduces millions of returns to a handful of binned medians and then
fits something to them, so the binning and the error bars decide the answer. Doing that
ad hoc per script produced two recurring mistakes, both of which this module exists to
prevent:

1. **Errors that count returns.** Returns inside one woodlot, field or swath are not
   independent, so ``SE = 1.2533 * NMAD / sqrt(n_returns)`` is too small -- by 4-5x on the
   Elba tiles. Fits weighted that way let a huge, homogeneous open-ground bin outvote the
   forest bins by ~20x more than its real information content. :func:`binned_stats` reports
   BOTH the naive return-based SE and a cluster-robust one built from spatial block
   medians, and prefers the latter.

2. **Bins chosen for round numbers.** A skewed covariate (canopy cover here is 64% near
   zero) puts almost all the leverage in bins holding almost no data. :func:`quantile_edges`
   spaces bins by quantile so each carries comparable weight, while always spanning the FULL
   observed range -- sparse extremes are kept and given honest uncertainty rather than
   truncated away, because in this data the extremes carry the largest effect.

The weighting rule follows from those two: bins with more INDEPENDENT information weigh
more, where independence is counted in spatial blocks rather than returns.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["nmad", "block_ids", "quantile_edges", "binned_stats", "BinStats"]


def nmad(a):
    """Normalised median absolute deviation -- a median-consistent robust sigma."""
    a = np.asarray(a, float)
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def block_ids(cell, nx, res, block_m):
    """Spatial block index per observation, from a C-order ``cell = iy*nx + ix``.

    Blocks are the unit of independence for the cluster-robust error: make them larger than
    the correlation length of whatever is being measured (canopy patches, swath overlap
    stripes), not larger than the features being compared.
    """
    cell = np.asarray(cell, np.int64)
    step = max(1, int(round(float(block_m) / float(res))))
    return (cell // nx // step) * (2**32) + (cell % nx // step)


@dataclass
class BinStats:
    """Binned robust summary. ``se`` is cluster-robust where blocks were supplied."""
    x: np.ndarray            # mean covariate value in the bin
    y: np.ndarray            # median response
    n: np.ndarray            # returns in the bin
    n_block: np.ndarray      # independent spatial blocks in the bin
    se_return: np.ndarray    # naive SE, counting returns (too small when clustered)
    se_block: np.ndarray     # cluster-robust SE, counting blocks
    se: np.ndarray           # the one to use: se_block where available, else se_return
    lo: np.ndarray           # bin edges
    hi: np.ndarray

    @property
    def weights(self):
        """Inverse-variance weights on the preferred SE: more independent data, more weight."""
        return 1.0 / np.maximum(self.se, 1e-12) ** 2

    def __len__(self):
        return self.x.size


def quantile_edges(x, nbins, *, first_edge=None, span_all=True):
    """Bin edges at equal quantiles of ``x``, spanning the full observed range.

    ``first_edge`` splits off a leading spike (e.g. the mass of near-zero canopy cover) into
    its own bin before quantiling the remainder, so that spike does not swallow the bins.
    With ``span_all`` the last edge is nudged past ``max(x)`` so no observation is dropped --
    the sparse tail is kept and carries its own (large) uncertainty.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("no finite values to bin")
    if first_edge is None:
        edges = np.quantile(x, np.linspace(0, 1, int(nbins) + 1))
    else:
        upper = x[x >= first_edge]
        q = (np.quantile(upper, np.linspace(0, 1, int(nbins) + 1)) if upper.size
             else np.array([first_edge]))
        edges = np.concatenate([[x.min(), float(first_edge)], q])
    edges = np.unique(edges)
    if span_all:
        edges[0] = min(edges[0], x.min())
        edges[-1] = np.nextafter(max(edges[-1], x.max()), np.inf)
    return edges


def binned_stats(x, y, edges, *, block=None, min_n=1, min_block=2):
    """Robust binned medians of ``y`` vs ``x`` with cluster-robust standard errors.

    ``block`` gives each observation a spatial block id (see :func:`block_ids`); the
    cluster-robust SE is then the robust spread of the per-block medians divided by
    sqrt(number of blocks). Bins with fewer than ``min_block`` blocks fall back to the
    return-based SE and are still reported -- a bin is dropped only for having fewer than
    ``min_n`` observations, never for being inconveniently sparse.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    edges = np.asarray(edges, float)
    blk = None if block is None else np.asarray(block)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if blk is not None:
        blk = blk[ok]
    cols = {k: [] for k in ("x", "y", "n", "n_block", "se_return", "se_block", "se", "lo", "hi")}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() < min_n:
            continue
        v = y[m]
        s_ret = 1.2533 * nmad(v) / np.sqrt(m.sum())
        nb, s_blk = 0, np.nan
        if blk is not None:
            b = blk[m]
            uniq = np.unique(b)
            if uniq.size >= min_block:
                bm = np.array([np.median(v[b == u]) for u in uniq])
                nb = uniq.size
                s_blk = 1.2533 * nmad(bm) / np.sqrt(nb)
                if not np.isfinite(s_blk) or s_blk <= 0:      # all block medians identical
                    s_blk = s_ret
        cols["x"].append(float(x[m].mean())); cols["y"].append(float(np.median(v)))
        cols["n"].append(int(m.sum())); cols["n_block"].append(int(nb))
        cols["se_return"].append(float(s_ret)); cols["se_block"].append(float(s_blk))
        cols["se"].append(float(s_blk if np.isfinite(s_blk) else s_ret))
        cols["lo"].append(float(lo)); cols["hi"].append(float(hi))
    return BinStats(**{k: np.asarray(v) for k, v in cols.items()})
