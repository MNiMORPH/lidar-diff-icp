"""Shared pieces of the gen1 return-shape fingerprint: the pooled histogram and the shape
features, copied from analysis/rcg/return_shape.py so that fingerprint_fit.py (which fits
the model on Elba) and fingerprint_apply.py (which applies it to any tile) compute the SAME
quantities from the same code. Nothing here is new: every constant and every line of
shape_feats is return_shape.py's, including its choices (20 mm bins over +/-1.5 m; 3x3
pooling by scipy uniform_filter with mode="nearest", times 9, rounded; the half-sample
mode read off the histogram; +/-50 mm peak window; +/-100 mm tail thresholds; the e-folding
length as the mean excess above mode+100 mm; a minimum of 30 pooled returns per cell).

ONE INCIDENTAL CHOICE OF return_shape.py IS REPRODUCED ON PURPOSE: its run() joins the
per-cell median intensity (a median over the cell's OWN in-range returns) into the feature
frame before the all-finite filter, so a cell with no return of its own -- even with >= 30
pooled returns from its neighbours -- is dropped. Intensity is not in the model, but the
population is: with it, Elba pass B has 12,135 cells / 926 positives (the record's run);
without it 12,158 / 932. So the domain here also requires OWN_RETURNS >= 1
(REQUIRE_OWN_RETURN), and gen1_pooled_hist returns the own-cell count for that.

The only departures from return_shape.py are mechanical, and verified to give the same
numbers on Elba (fingerprint_apply.py asserts it):
  - the 3x3 pooling runs in bin-chunks so the float32 copy never holds the whole
    (cells x bins) array twice (the filter kernel is size 1 along the bin axis, so chunking
    over bins cannot change a value);
  - the histogram is accumulated directly in float32 (exact for counts < 2^24; the
    original accumulated int32 and cast to float32 before the filter).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import uniform_filter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from hsm import hsm_from_hist  # noqa: E402

# return_shape.py: LO, HI, DZ = -1500.0, 1500.0, 20.0
LO, HI, DZ = -1500.0, 1500.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
# return_shape.py run(): keep = F.sum(1) >= 30  (pooled returns per cell)
MIN_POOLED_RETURNS = 30
# return_shape.py run(): the intensity extra (median over the cell's OWN in-range returns) is NaN
# where the cell has none, and the all-finite filter then drops the cell -- see the module docstring
REQUIRE_OWN_RETURN = True
# return_shape.py run(): shape_cols = [c for c in feats.columns if c not in ("n", "intensity", "c2_share")]
FEATURES = ["lo_w", "hi_w", "mode_med", "bowley", "peak", "tail_L", "tail_mass", "low_mass"]
FEATURE_DEFS = {
    "lo_w": "p50 - p10 of the pooled gen1 return heights, mm (lower width)",
    "hi_w": "p90 - p50, mm (upper width)",
    "mode_med": "half-sample mode (hsm_from_hist, middle-of-ties) - p50, mm",
    "bowley": "(p75 + p25 - 2 p50) / max(p75 - p25, 20 mm) (Bowley skew, dimensionless)",
    "peak": "share of pooled returns within +/-50 mm of the mode",
    "tail_L": "mean of (h - mode - 100 mm) over returns above mode + 100 mm (e-folding length of the upper tail), mm",
    "tail_mass": "share of pooled returns above mode + 100 mm",
    "low_mass": "share of pooled returns below mode - 100 mm",
}
# quantiles are read as the first bin centre whose cumulative share reaches p (return_shape.py q())


def grid_from_corrections(tile):
    """Grid geometry from <tile>/corrections.json: bounds [X0, Y0, X1, Y1] and res_m."""
    with open(os.path.join(tile, "corrections.json")) as f:
        c = json.load(f)
    X0, Y0, X1, Y1 = (float(v) for v in c["bounds"]); res = float(c["res_m"])
    ny = int(round((Y1 - Y0) / res)); nx = int(round((X1 - X0) / res))
    return X0, Y0, res, ny, nx


def gen1_pooled_hist(parquet_path, ny, nx, chunk_bins=25):
    """3x3-pooled per-cell histograms of gen1 CSF ground returns about the slope-normal
    gen2 surface (column d_mm_corr), exactly as return_shape.py builds H1p.

    Returns (H1p, own, n_rows_in_grid, n_rows_in_range): H1p is float32 (ny*nx, NB); own is
    the UNPOOLED in-range return count per cell (ny*nx,), for REQUIRE_OWN_RETURN."""
    t = pq.read_table(parquet_path, columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
    cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); del t
    ok = np.isfinite(d) & (d > LO) & (d < HI)
    H1 = np.zeros((ny * nx, NB), np.float32)
    np.add.at(H1, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1)
    own = np.bincount(cell[ok], minlength=ny * nx)
    H1 = H1.reshape(ny, nx, NB); H1p = np.empty_like(H1)
    for b0 in range(0, NB, chunk_bins):
        b1 = min(NB, b0 + chunk_bins)
        H1p[:, :, b0:b1] = np.rint(uniform_filter(H1[:, :, b0:b1], size=(3, 3, 1), mode="nearest") * 9)
    return H1p.reshape(-1, NB), own, int(g.sum()), int(ok.sum())


def shape_feats(F):
    """return_shape.py shape_feats, verbatim. F: (cells, NB) pooled counts."""
    tot = F.sum(1); P = F / np.maximum(tot, 1)[:, None]; C = np.cumsum(P, 1)
    q = lambda p: ctr[np.argmax(C >= p, axis=1)]  # noqa: E731
    p10, p25, p50, p75, p90 = (q(p) for p in (0.1, 0.25, 0.5, 0.75, 0.9))
    mode = np.array([hsm_from_hist(F[i], ctr) if tot[i] > 0 else np.nan for i in range(len(F))])
    rel = ctr[None, :] - mode[:, None]
    peak = (P * ((rel >= -50) & (rel <= 50))).sum(1); tail_mass = (P * (rel > 100)).sum(1); low_mass = (P * (rel < -100)).sum(1)
    w = P * (rel > 100); tail_L = (w * (rel - 100)).sum(1) / np.maximum(w.sum(1), 1e-9)
    return pd.DataFrame(dict(lo_w=p50 - p10, hi_w=p90 - p50, mode_med=mode - p50, bowley=(p75 + p25 - 2 * p50) / np.maximum(p75 - p25, DZ), peak=peak, tail_L=tail_L, tail_mass=tail_mass, low_mass=low_mass, n=tot)), mode, p50


def blocks(idx, nx):
    """return_shape.py blocks(): 5 folds over 100-cell (500 m at 5 m) blocks, rng(0)."""
    blk = ((idx // nx) // 100) * 1000 + ((idx % nx) // 100); ub, binv = np.unique(blk, return_inverse=True)
    return np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]


def score_from_model(model, feats):
    """LDA decision function from the saved JSON model: ((x - mean) / std) @ coef + intercept."""
    X = feats[model["feature_names"]].to_numpy(float)
    Z = (X - np.asarray(model["mean"])) / np.asarray(model["std"])
    return Z @ np.asarray(model["coef"]) + float(model["intercept"])
