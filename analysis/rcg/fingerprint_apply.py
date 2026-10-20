#!/usr/bin/env python3
"""Apply the saved gen1 return-shape fingerprint (data/derived/rcg/fingerprint_lda.json,
fitted on Elba by fingerprint_fit.py) to any tile directory and write per-cell scores.

INPUTS, in <tile_dir>:
  beam_offset_table.parquet   gen1 CSF ground returns; columns cell, d_mm_corr, in_grid
  corrections.json            grid geometry (bounds, res_m)
  floodplain_mask.npy         the domain -- or floodplain_mask_clubb.npy if present, which
                              is then used INSTEAD and announced in the banner
  slope.npy                   degrees; the domain keeps slope <= 2 (return_shape.py's flat)

DOMAIN. floodplain & slope <= 2 deg, then the model's own population rules, copied from
return_shape.py via fingerprint_common: >= 30 pooled returns in the 3x3 block, >= 1 in-range
return of the cell's own, all eight features finite. Everything else is NaN. NOTE the
training domain on Elba additionally required (i) a finite DoD and (ii) no gen2
full-density return above 2 m in the cell (open floodplain, from the 2021 cloud); neither is
available from the inputs above, so the scores here cover wooded and DoD-void floodplain
cells too -- an extrapolation of the model beyond its training population, which the banner
counts and the record (FINGERPRINT_MODEL.md) flags.

OUTPUTS, in <tile_dir>:
  fingerprint_score.npy  (ny, nx) float64: the LDA decision function; NaN off-domain
  fingerprint_p.npy      (ny, nx) float64: 1 / (1 + exp(-score)); NaN off-domain

CHECK on the training tile: with --check-insample <npz> (default: the file fingerprint_fit.py
wrote, used automatically when the tile is the training tile), the scores at the training
cells are compared with the refit's in-sample scores and the max |difference| is asserted
below 1e-9 and printed.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/rcg/fingerprint_apply.py --tile data/derived/elbaext_delong_star_nosurf
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.getcwd())
from fingerprint_common import (DZ, FEATURE_DEFS, FEATURES, HI, LO, MIN_POOLED_RETURNS, NB, REQUIRE_OWN_RETURN,  # noqa: E402
                                gen1_pooled_hist, grid_from_corrections, score_from_model, shape_feats)
from trust.provenance import Run  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True, help="tile directory holding beam_offset_table.parquet, corrections.json, floodplain_mask.npy, slope.npy")
ap.add_argument("--model", default="data/derived/rcg/fingerprint_lda.json")
ap.add_argument("--check-insample", default=None, help="npz with idx, score from fingerprint_fit.py; default: auto on the training tile")
ap.add_argument("--slope-max", type=float, default=2.0, help="degrees; return_shape.py's flat domain (repo default, do not change for comparability)")
A = ap.parse_args(); T = A.tile.rstrip("/")

R = Run(f"Where on tile {os.path.basename(T)} does the gen1 return shape carry the Elba-fitted fingerprint of the standing herbaceous sward?")
with open(R.input(A.model, role="the fitted LDA: feature names, standardisation mean/std, coef, intercept; fitted by fingerprint_fit.py")) as f: M = json.load(f)
X0, Y0, res, ny, nx = grid_from_corrections(T); R.input(f"{T}/corrections.json", role="grid geometry: bounds and res_m")
clubb = os.path.exists(f"{T}/floodplain_mask_clubb.npy"); fp_path = f"{T}/floodplain_mask_clubb.npy" if clubb else f"{T}/floodplain_mask.npy"
fp = np.load(R.input(fp_path, role=("Clubb (lsdtt-valley-metrics) floodplain mask, USED AS THE DOMAIN because it is present" if clubb else "pipeline floodplain mask, the domain"))).astype(bool)
sl = np.load(R.input(f"{T}/slope.npy", role="slope, degrees; domain keeps <= slope-max"))
pq_path = R.input(f"{T}/beam_offset_table.parquet", role="gen1 CSF ground returns: cell, d_mm_corr (mm above the slope-normal gen2 surface, corrected), in_grid")
assert fp.shape == (ny, nx) == sl.shape, (fp.shape, sl.shape, ny, nx)
R.param("slope_max_deg", A.slope_max, src="repo"); R.param("min_pooled_returns", MIN_POOLED_RETURNS, src="repo"); R.param("require_own_return", REQUIRE_OWN_RETURN, src="repo")
R.param("bins", f"{DZ:.0f} mm over [{LO:.0f}, {HI:.0f}] mm, {NB} bins; 3x3 pooled", src="repo"); R.param("features", FEATURES, src="repo")
R.param("model_training", f"{M['training_tile']}; n {M['n_cells']}; AUC CV {M['auc_cv_500m_block']:.3f}, resub {M['auc_resubstitution']:.3f}; git {M['git_commit'][:12]}", src="repo")
dom = fp & (sl <= A.slope_max)
R.mask("floodplain", fp, defn=("floodplain_mask_clubb.npy == True" if clubb else "floodplain_mask.npy == True"))
R.mask("domain", dom, defn=f"floodplain & slope <= {A.slope_max:g} deg (return_shape.py's flat, minus its finite-DoD and gen2-open conditions, which need files this tile need not have)")
for k, v in FEATURE_DEFS.items(): R.column(k, v)
R.column("score", "LDA decision function on the standardised shape features, dimensionless; > 0 leans towards the imagery class"); R.column("p", "1 / (1 + exp(-score)), the LDA's class-1 posterior under its fitted priors (dimensionless)")
R.column("what", "row label of the summary quantity"); R.column("value", "the summary quantity, count or dimensionless score as the row label says")
R.banner()

H1p, own, n_grid, n_range = gen1_pooled_hist(pq_path, ny, nx)
print(f"gen1 rows in_grid {n_grid:,}; in [{LO:.0f},{HI:.0f}] mm {n_range:,}", flush=True)
idx = np.flatnonzero(dom.ravel()); F = H1p[idx]; keep = F.sum(1) >= MIN_POOLED_RETURNS; n_dom = idx.size; idx, F = idx[keep], F[keep]; del H1p
feats, mode, p50 = shape_feats(F); del F
if REQUIRE_OWN_RETURN: feats.loc[own[idx] < 1, FEATURES[0]] = np.nan
fin = np.all(np.isfinite(feats[FEATURES].to_numpy(float)), 1); feats, idx = feats[fin], idx[fin]
score = score_from_model(M, feats); p = 1 / (1 + np.exp(-score))
S_ = np.full(ny * nx, np.nan); S_[idx] = score; P_ = np.full(ny * nx, np.nan); P_[idx] = p
np.save(f"{T}/fingerprint_score.npy", S_.reshape(ny, nx)); np.save(f"{T}/fingerprint_p.npy", P_.reshape(ny, nx))
rows = [["domain cells", n_dom], [f"with >= {MIN_POOLED_RETURNS} pooled returns", int(keep.sum())], ["scored (own return, finite features)", int(idx.size)],
        ["score median", f"{np.median(score):.4f}"], ["score p10 / p90", f"{np.percentile(score, 10):.4f} / {np.percentile(score, 90):.4f}"],
        ["cells with p > 0.5", int((p > 0.5).sum())], ["share with p > 0.5", f"{(p > 0.5).mean():.4f}"]]
R.table(["what", "value"], rows)

chk = A.check_insample
if chk is None and os.path.abspath(T) == os.path.abspath(M["training_tile"]) and os.path.exists("data/derived/rcg/fingerprint_fit_insample.npz"): chk = "data/derived/rcg/fingerprint_fit_insample.npz"
if chk:
    z = np.load(R.input(chk, role="fingerprint_fit.py's in-sample training cells (idx) and their refit scores (score)"))
    tr_idx, tr_score = z["idx"], z["score"]; here = S_[tr_idx]; missing = int(np.isnan(here).sum())
    assert missing == 0, f"{missing} training cells are not scored by the apply domain"
    dmax = float(np.max(np.abs(here - tr_score)))
    print(f"\nIN-SAMPLE CHECK on the training tile: {tr_idx.size:,} training cells; max |apply - refit| = {dmax:.3e}")
    assert dmax < 1e-9, dmax
    R.notes.append(f"in-sample check: {tr_idx.size} cells, max abs diff {dmax:.3e}")
print(f"\nwrote {T}/fingerprint_score.npy and {T}/fingerprint_p.npy ({ny} x {nx})")
R.done(headline=f"{idx.size} cells scored; {int((p > 0.5).sum())} with p > 0.5")
