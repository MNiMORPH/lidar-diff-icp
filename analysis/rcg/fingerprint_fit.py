#!/usr/bin/env python3
"""Refit the gen1 return-shape fingerprint on the Elba tile and SAVE it, so that the detector
return_shape.py only ever evaluated in-script can be applied to another tile.

WHAT IS FITTED (return_shape.py, gen1 pass B, reproduced here line for line):
  features  the eight translation-invariant shape features of the 3x3-pooled gen1 CSF
            return-height histogram about the cell's own half-sample mode
            (fingerprint_common.FEATURES; definitions in FEATURE_DEFS). Return count and
            intensity are NOT in the model (return_shape.py's "shape only" set).
  label     the 2008 imagery class of the survey's own year: NE-like probability from
            Landsat-5 2008-10-29 (rcg/nelike_p_ls_20081029.tif, map_ne_like.py), sampled
            to the 5 m grid by nearest 30 m pixel, in class at p > 0.5.
  domain    open flat floodplain, NE patch EXCLUDED (return_shape.py dom_b):
              flat  = floodplain_mask & slope <= 2 deg & finite DoD
              open5 = no gen2 full-density return above 2 m in the cell
                      (rcg/gen2_cellstats.npz, above2m == 0; lidar_fingerprint.py)
              NE    = the largest connected patch of dod_floodplain_vegcorr < -LoD on the
                      floodplain (8-connected), excluded because the lidar helped define it
            then cells with >= 30 pooled returns, finite features and finite logit(p),
            AND at least one in-range return of the cell's own (return_shape.py's
            intensity extra, a median over the cell's own returns, is NaN otherwise and its
            all-finite filter drops the cell; fingerprint_common.REQUIRE_OWN_RETURN).
  model     sklearn LinearDiscriminantAnalysis (default solver "svd"), fitted on the
            STANDARDISED features (z = (x - mean) / std, std = population std as X.std(0),
            the same scaling return_shape.py used to report its weights). LDA is affine
            equivariant, so this gives the same decision function as the raw-feature fit
            return_shape.py performed; the script checks that and prints the max |diff|.
  scoring   AUC of the decision function, 5-fold CV by 500 m block (return_shape.py
            blocks(), rng(0)); and the resubstitution AUC of the fit on all cells.
  reference pass A (NE patch INCLUDED, return_shape.py dom_a) is fitted and scored the
            same way and recorded in the JSON for comparison, but is NOT the saved model.

Choices copied from return_shape.py and NOT changed here: 20 mm bins over +/-1.5 m; the
3x3 pooling; the 30-return minimum; the +/-50 mm peak window and +/-100 mm tail thresholds;
p > 0.5 for the class; slope <= 2 deg; 100-cell blocks; the rng(0) fold permutation.

Output: data/derived/rcg/fingerprint_lda.json (feature names, mean, std, coef, intercept,
the raw-scale coefficients too, the domain definition, the AUCs, cell counts, git commit).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/rcg/fingerprint_fit.py
"""
import json
import os
import subprocess
import sys

import numpy as np
import rioxarray
import sklearn
from scipy.ndimage import label
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.getcwd())
from fingerprint_common import (DZ, FEATURE_DEFS, FEATURES, HI, LO, MIN_POOLED_RETURNS, NB, REQUIRE_OWN_RETURN, blocks,  # noqa: E402
                                gen1_pooled_hist, grid_from_corrections, score_from_model, shape_feats)
from trust.provenance import Run  # noqa: E402

S = "data/derived/elbaext_delong_star_nosurf"; R_ = f"{S}/rcg"; OUT = "data/derived/rcg/fingerprint_lda.json"
R = Run("What are the gen1 return-shape LDA coefficients that return_shape.py fitted on Elba, saved so the detector can be applied to another tile?")
X0, Y0, res, ny, nx = grid_from_corrections(S); R.input(f"{S}/corrections.json", role="grid geometry: bounds and res_m (5 m, EPSG:26915)")
dod = np.load(R.input(f"{S}/dod.npy", role="DoD gen2 - gen1, m; only its finiteness enters (flat domain)")) * 1e3
new = np.load(R.input(f"{S}/dod_floodplain_vegcorr.npy", role="vegetation-corrected floodplain DoD, m; defines the NE patch (< -LoD)")) * 1e3
fp = np.load(R.input(f"{S}/floodplain_mask.npy", role="pipeline floodplain (valley-floor) mask, bool")).astype(bool)
sl = np.load(R.input(f"{S}/slope.npy", role="gen2 slope, degrees")); lod = np.load(R.input(f"{S}/lod.npy", role="per-cell LoD, m; NE patch threshold")) * 1e3
assert dod.shape == (ny, nx), (dod.shape, ny, nx)
pq_path = R.input(f"{S}/beam_offset_table.parquet", role="gen1 CSF ground returns: cell index, d_mm_corr = height above the slope-normal gen2 surface after the star/lateral corrections, in_grid")
zc = np.load(R.input(f"{R_}/gen2_cellstats.npz", role="gen2 full-density per-cell stats (lidar_fingerprint.py cache); above2m defines OPEN cells"))
tif = R.input(f"{R_}/nelike_p_ls_20081029.tif", role="imagery class: p(NE-like) from Landsat-5 2008-10-29, 30 m (map_ne_like.py) -- the TRAINING LABEL")
R.param("bins", f"{DZ:.0f} mm over [{LO:.0f}, {HI:.0f}] mm, {NB} bins", src="repo"); R.param("pooling", "3x3 cells, uniform_filter mode=nearest, x9, rounded", src="repo")
R.param("min_pooled_returns", MIN_POOLED_RETURNS, src="repo"); R.param("require_own_return", REQUIRE_OWN_RETURN, src="repo", why="return_shape.py's intensity extra drops cells with no return of their own (see fingerprint_common docstring)"); R.param("class_threshold_p", 0.5, src="repo"); R.param("slope_max_deg", 2, src="repo")
R.param("open_cell", "gen2 above2m == 0", src="repo"); R.param("cv", "5 folds by 100-cell (500 m) block, rng(0)", src="repo")
R.param("features", FEATURES, src="repo"); R.param("lda_solver", "svd (sklearn default)", src="repo"); R.param("standardise", "z = (x - mean) / X.std(0) before the LDA; affine-equivalent to return_shape.py's raw fit (checked below)", src="repo")

flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod)
lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
open5 = np.zeros(ny * nx, bool); open5[zc["cells"][zc["above2m"] == 0]] = True; open5 = open5.reshape(ny, nx)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):  # return_shape.py sample(): nearest 30 m pixel at each 5 m cell centre
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(tif); logit = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4)))  # noqa: E731
R.mask("flat", flat, defn="floodplain_mask & slope <= 2 deg & finite dod"); R.mask("open5", open5, defn="gen2 full-density above2m == 0 (lidar_fingerprint.py cache)")
R.mask("NE", NE, defn="largest 8-connected patch of dod_floodplain_vegcorr < -lod on the floodplain (EXCLUDED from training)")
R.mask("dom_b (training)", flat & open5 & ~NE, defn="flat & open5 & ~NE (return_shape.py dom_b)"); R.mask("dom_a (reference)", flat & open5, defn="flat & open5 (return_shape.py dom_a)")
R.mask("class_2008", flat & open5 & ~NE & (p08 > 0.5), of=int((flat & open5 & ~NE).sum()), defn="p(NE-like, Landsat-5 2008-10-29) > 0.5 within dom_b")
for k, v in FEATURE_DEFS.items(): R.column(k, v)
R.column("pass_", "B = NE patch excluded (the saved model); A = NE patch included (reference only)")
R.column("cells", "count of domain cells with >= 30 pooled returns and finite features/label"); R.column("positives", "count of those cells in the 2008 imagery class")
R.column("auc_cv", "ROC AUC of the LDA decision function, 5-fold CV by 500 m block (dimensionless)"); R.column("auc_resub", "ROC AUC of the LDA fitted on all cells, scored on the same cells (resubstitution)")
R.column("weight", "LDA coefficient on the standardised feature (per 1 population-sd of the feature; sign + = towards the class)")
R.banner()

H1p, own, n_grid, n_range = gen1_pooled_hist(pq_path, ny, nx)
print(f"gen1 rows in_grid {n_grid:,}; in [{LO:.0f},{HI:.0f}] mm {n_range:,}; pooled histogram {H1p.shape}", flush=True)
R.notes.append(f"gen1 rows in_grid {n_grid}; within histogram range {n_range}")

def fit(domain, name):
    idx = np.flatnonzero(domain.ravel()); F = H1p[idx]; keep = F.sum(1) >= MIN_POOLED_RETURNS; idx, F = idx[keep], F[keep]
    feats, mode, p50 = shape_feats(F)
    if REQUIRE_OWN_RETURN: feats.loc[own[idx] < 1, FEATURES[0]] = np.nan  # dropped by the all-finite filter below, as in return_shape.py
    y = (p08 > 0.5).ravel()[idx].astype(int); s_ = logit(p08).ravel()[idx]; fin = np.all(np.isfinite(feats.to_numpy(float)), 1) & np.isfinite(s_)
    feats, y, idx = feats[fin], y[fin], idx[fin]; fold = blocks(idx, nx)
    X = feats[FEATURES].to_numpy(float); mean = X.mean(0); std = X.std(0); Z = (X - mean) / std
    pred = np.zeros(y.size)
    for k in range(5):
        tr = fold != k; m = LinearDiscriminantAnalysis().fit(Z[tr], y[tr]); pred[~tr] = m.decision_function(Z[~tr])
    pred_raw = np.zeros(y.size)  # return_shape.py's own form, raw features: the CV AUC must agree
    for k in range(5):
        tr = fold != k; m = LinearDiscriminantAnalysis().fit(X[tr], y[tr]); pred_raw[~tr] = m.decision_function(X[~tr])
    auc_cv = roc_auc_score(y, pred); auc_cv_raw = roc_auc_score(y, pred_raw)
    lda = LinearDiscriminantAnalysis().fit(Z, y); lda_raw = LinearDiscriminantAnalysis().fit(X, y)
    score = lda.decision_function(Z); score_raw = lda_raw.decision_function(X); auc_resub = roc_auc_score(y, score)
    model = dict(feature_names=FEATURES, mean=mean.tolist(), std=std.tolist(), coef=lda.coef_[0].tolist(), intercept=float(lda.intercept_[0]),
                 coef_raw=lda_raw.coef_[0].tolist(), intercept_raw=float(lda_raw.intercept_[0]))
    dmax = float(np.max(np.abs(score_from_model(model, feats) - score_raw)))
    print(f"\n[{name}] cells {idx.size:,}; positives {y.sum():,}; folds sizes {np.bincount(fold).tolist()}")
    print(f"  AUC block-CV {auc_cv:.4f} (raw-feature form {auc_cv_raw:.4f}); AUC resubstitution {auc_resub:.4f}")
    print(f"  max |score(standardised model) - score(raw LDA)| = {dmax:.3e}  (affine equivariance check)")
    R.table(["weight"] + FEATURES, [[name] + [f"{c:+.3f}" for c in lda.coef_[0]]])
    return dict(pass_=name, cells=int(idx.size), positives=int(y.sum()), auc_cv=float(auc_cv), auc_cv_raw_form=float(auc_cv_raw), auc_resub=float(auc_resub),
                affine_check_max_abs_diff=dmax, fold_sizes=np.bincount(fold).tolist(), model=model, idx=idx, score=score)

B = fit(flat & open5 & ~NE, "B"); A = fit(flat & open5, "A")
R.table(["pass_", "cells", "positives", "auc_cv", "auc_resub"], [[p["pass_"], p["cells"], p["positives"], f"{p['auc_cv']:.4f}", f"{p['auc_resub']:.4f}"] for p in (B, A)])
assert B["affine_check_max_abs_diff"] < 1e-6, B["affine_check_max_abs_diff"]

def git(*a):
    try: return subprocess.check_output(["git", *a], text=True).strip()
    except Exception as e: return f"unavailable ({e})"
out = dict(
    model="gen1 return-shape fingerprint: LDA on standardised shape features; score = ((x - mean) / std) @ coef + intercept  (equivalently x @ coef_raw + intercept_raw); p = 1 / (1 + exp(-score))",
    fitted_from="analysis/rcg/fingerprint_fit.py, reproducing analysis/rcg/return_shape.py gen1 pass B (NE patch excluded)",
    git_commit=git("rev-parse", "HEAD"), git_dirty=bool(git("status", "--porcelain", "--", "analysis/rcg", "analysis/hsm.py", "trust")), sklearn=sklearn.__version__,
    training_tile=S, label=dict(raster=tif, defn="p(NE-like) from Landsat-5 2008-10-29 two-class ML (map_ne_like.py), nearest 30 m pixel at the 5 m cell centre; in class at p > 0.5", caveat="an imagery spectral class, not a species map"),
    training_domain=dict(defn="flat & open5 & ~NE, then >= 30 pooled returns, >= 1 own in-range return (REQUIRE_OWN_RETURN), finite features, finite logit(p)",
                         flat="floodplain_mask & slope <= 2 deg & finite dod", open5="gen2 full-density cloud has no return above 2 m in the cell (rcg/gen2_cellstats.npz above2m == 0)",
                         NE="largest 8-connected patch of dod_floodplain_vegcorr < -lod on the floodplain, EXCLUDED",
                         counts=dict(flat=int(flat.sum()), open5=int(open5.sum()), NE=int(NE.sum()), dom_b=int((flat & open5 & ~NE).sum()))),
    histogram=dict(source="beam_offset_table.parquet columns cell, d_mm_corr, in_grid", lo_mm=LO, hi_mm=HI, bin_mm=DZ, n_bins=NB, pooling="3x3 cells, scipy uniform_filter mode=nearest, x9, rint", min_pooled_returns=MIN_POOLED_RETURNS, require_own_return=REQUIRE_OWN_RETURN),
    feature_names=FEATURES, feature_defs=FEATURE_DEFS, **{k: B["model"][k] for k in ("mean", "std", "coef", "intercept", "coef_raw", "intercept_raw")},
    n_cells=B["cells"], n_positives=B["positives"], auc_resubstitution=B["auc_resub"], auc_cv_500m_block=B["auc_cv"], cv=dict(folds=5, block_cells=100, block_m=500, seed=0, fold_sizes=B["fold_sizes"]),
    affine_check_max_abs_diff=B["affine_check_max_abs_diff"],
    reference_pass_A_NE_included=dict(n_cells=A["cells"], n_positives=A["positives"], auc_resubstitution=A["auc_resub"], auc_cv_500m_block=A["auc_cv"], coef=A["model"]["coef"], intercept=A["model"]["intercept"], mean=A["model"]["mean"], std=A["model"]["std"]),
    on_record="RCG_PHENOLOGY.md section S: shape-only AUC 0.71 (patch excluded), 0.82 (patch included)",
)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f: json.dump(out, f, indent=1)
np.savez_compressed("data/derived/rcg/fingerprint_fit_insample.npz", idx=B["idx"], score=B["score"])
print(f"\nmodel: {OUT}; in-sample scores (for fingerprint_apply.py's check): data/derived/rcg/fingerprint_fit_insample.npz")
R.done(headline=f"pass B: n {B['cells']}, AUC CV {B['auc_cv']:.3f}, resub {B['auc_resub']:.3f}; pass A: AUC CV {A['auc_cv']:.3f}")
