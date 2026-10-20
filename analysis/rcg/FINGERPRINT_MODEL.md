# The gen1 return-shape fingerprint, saved and applicable (2026-10-20)

**Model file:** `data/derived/rcg/fingerprint_lda.json` (tracked).
**Fit:** `analysis/rcg/fingerprint_fit.py`. **Apply:** `analysis/rcg/fingerprint_apply.py --tile <dir>`.
**Shared feature code:** `analysis/rcg/fingerprint_common.py` (copied from `return_shape.py`).
Context: `analysis/RCG_PHENOLOGY.md` section S, which fitted this detector in-script on Elba and
reported block-CV AUC 0.71 (NE patch excluded) and 0.82 (included) without saving it.

## What the model is

A linear discriminant (sklearn `LinearDiscriminantAnalysis`, default `svd` solver) on eight
translation-invariant features of the SHAPE of the per-cell distribution of gen1 (November
2008) CSF ground returns about the cell's own mode. Nothing cross-epoch enters the features:
the heights are `d_mm_corr` (mm above the slope-normal gen2 surface after the star and
lateral corrections), but every feature is a difference, share or length measured relative to
the cell's own mode or median, so the gen2 surface only supplies the frame. The detector runs
on the 2008 survey by itself.

The score is the LDA decision function,

    score = ((x - mean) / std) @ coef + intercept        (= x @ coef_raw + intercept_raw)
    p     = 1 / (1 + exp(-score))

with `mean`, `std` (population std, `X.std(0)`), `coef`, `intercept`, `coef_raw`,
`intercept_raw` all in the JSON. `p` is the LDA's class-1 posterior under the TRAINING prior
(926 of 12,135 cells = 7.6% positive), so `p > 0.5` is a conservative cut; what the AUC
validates is the RANKING by `score`, not a threshold. No threshold has been chosen or tuned.

## The features, exactly

Histogram: gen1 returns with `in_grid`, finite `d_mm_corr` in (-1500, 1500) mm, binned at
20 mm (150 bins, centres `ctr`), per 5 m cell; then 3x3-pooled with
`scipy.ndimage.uniform_filter(size=(3,3,1), mode="nearest") * 9`, rounded (`np.rint`). A cell
is scored only if its pooled histogram holds >= 30 returns AND the cell has >= 1 in-range
return of its own (see "one incidental choice" below). Quantiles `p10 ... p90` are read as
the first bin centre whose cumulative share reaches the level; `mode` is the half-sample mode
read off the histogram (`analysis/hsm.py` `hsm_from_hist`, ties broken at the middle window);
`rel = ctr - mode`; `P` = bin share.

| feature | definition |
|---|---|
| `lo_w` | p50 - p10, mm |
| `hi_w` | p90 - p50, mm |
| `mode_med` | mode - p50, mm |
| `bowley` | (p75 + p25 - 2 p50) / max(p75 - p25, 20 mm) |
| `peak` | share of returns with -50 <= rel <= 50 mm |
| `tail_L` | mean of (rel - 100) over returns with rel > 100 mm (upper-tail e-folding length), mm |
| `tail_mass` | share of returns with rel > 100 mm |
| `low_mass` | share of returns with rel < -100 mm |

Return count `n` and intensity are NOT in the model (`return_shape.py`'s "shape only" set).

## The training population (Elba, `data/derived/elbaext_delong_star_nosurf`)

Label: the 2008 imagery class of the survey's own year, `rcg/nelike_p_ls_20081029.tif`
(p(NE-like) from a two-class maximum-likelihood on Landsat-5 2008-10-29, `map_ne_like.py`),
sampled at each 5 m cell centre by nearest 30 m pixel; in class at p > 0.5.

Domain (`return_shape.py` pass B, "the independent fingerprint"), with the counts the fit
printed:

    flat   = floodplain_mask & slope <= 2 deg & finite dod          54,548 cells
    open5  = no gen2 full-density return above 2 m in the cell      13,579
             (rcg/gen2_cellstats.npz above2m == 0)
    NE     = largest 8-connected patch of dod_floodplain_vegcorr     2,355   EXCLUDED
             < -lod on the floodplain (the lidar helped define it)
    dom_b  = flat & open5 & ~NE                                     12,215
    then >= 30 pooled returns, >= 1 own return, finite features and finite logit(p):
                                                                     12,135 cells, 926 positives

One incidental choice reproduced on purpose: `return_shape.py` joins the per-cell median
intensity (a median over the cell's OWN in-range returns) into the feature frame before its
all-finite filter, so a cell with no return of its own is dropped even when its 3x3 pool has
>= 30. Without it the population is 12,158 / 932 and the CV AUC 0.7054; with it the numbers
below match `return_shape.py`'s own run to the digit (re-run at HEAD 2026-10-20: 12,135 cells,
926 in class, AUC 0.708 / 0.818). `fingerprint_common.REQUIRE_OWN_RETURN` carries it.

## The numbers the fit printed (2026-10-20, git a9c8011, sklearn 1.4.1.post1)

    pass B (saved model):  cells 12,135; positives 926
      AUC, 5-fold CV by 500 m block:  0.7077   (record: 0.71)
      AUC, resubstitution:            0.8100
      fold sizes [5771, 849, 2993, 1209, 1313]
    pass A (NE patch included, reference only, in the JSON under reference_pass_A_NE_included):
      cells 13,499; positives 2,168; AUC CV 0.8177 (record: 0.82); resubstitution 0.8731

    standardised weights (coef, per 1 population sd; + = towards the class):
      lo_w -0.197  hi_w -0.143  mode_med -0.112  bowley +0.055  peak +0.486
      tail_L -0.325  tail_mass +1.876  low_mass +0.966        intercept -3.328
    mean  [101.89, 129.64, -8.77, 0.0447, 0.5353, 79.88, 0.1227, 0.0640]
    std   [ 81.35, 112.99, 44.65, 0.2281, 0.1770, 88.44, 0.1339, 0.0883]
    raw-scale coef [-0.002426, -0.001266, -0.002500, 0.2411, 2.748, -0.003675, 14.01, 10.95]; intercept_raw -6.547

The weights equal those `return_shape.py` prints (tail_mass +1.88, low_mass +0.97, peak +0.49,
tail_L -0.33, lo_w -0.20, hi_w -0.14, mode_med -0.11, bowley +0.05). The LDA was fitted on
standardised features; LDA is affine equivariant, and the fit checks its decision function
against the raw-feature fit `return_shape.py` performed: max |diff| 1.0e-12.

## Applying it

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/rcg/fingerprint_apply.py --tile <tile_dir>

Needs in `<tile_dir>`: `beam_offset_table.parquet` (only `cell`, `d_mm_corr`, `in_grid` are
read), `corrections.json` (bounds, res_m), `floodplain_mask.npy` -- or
`floodplain_mask_clubb.npy`, used instead when present and announced in the banner -- and
`slope.npy`. Writes `fingerprint_score.npy` and `fingerprint_p.npy` on the tile grid, NaN
outside floodplain & slope <= 2 deg, below 30 pooled returns, or without an own return.
Elba, 2026-10-20: 58,767 domain cells, 54,541 scored, 3,804 with p > 0.5; peak RSS 1.6 GB.

**Verification on the training tile.** Applying the saved JSON reproduces the refit's
in-sample scores at all 12,135 training cells: max |apply - refit| = 8.9e-16 (asserted
< 1e-9 by the script).

**Held-out check, applied scores on Elba** (bare numbers; the NE patch was excluded from
training and 87% of its flat cells are in the 2008 imagery class):

    population                                            n     score median   share p>0.5   imagery class share
    NE patch, flat (held out)                           1,690      -0.910         0.163            0.873
    training domain, in 2008 class                        926      -1.203         0.199            1.000
    training domain, out of class                      11,209      -4.018         0.026            0.000
    flat floodplain NOT open (gen2 > 2 m): extrapolation 40,796    -2.244         0.074            0.052

## Caveats

1. **The label is an imagery spectral class, not a species.** "NE-like" means October-2008
   Landsat reflectance like the NE patch's (a bright, senesced, dense sward); the Wisconsin
   reports list its confusions (grazed pasture, cut hay). The detector learns to find that
   class in the return shape; it has never seen a field-verified reed-canary-grass polygon.
2. **It is a 2008 (leaf-off, standing senesced sward) capability.** Section S: the May-2021
   survey carries no shape signal at all. Apply it to gen1 only.
3. **The apply domain is wider than the training domain.** Training used open flat floodplain
   (no gen2 return above 2 m, finite DoD). `fingerprint_apply.py` scores all floodplain cells
   with slope <= 2 deg because the open screen needs the gen2 full-density cloud. On Elba,
   40,796 of the 54,541 scored cells are that extrapolation (wooded or shrubby floodplain),
   where the return shape has other causes. Mask with a canopy layer before interpreting.
4. **Presence, not cover.** Section S: block-CV R^2 of logit(p) on the shape features is 0.03.
5. **p is a posterior under the Elba training prior (7.6% positive).** Use the score's rank;
   do not read p > 0.5 as a calibrated presence threshold on another tile.
6. Fitted on one tile, one survey vendor configuration (SE-MN 2008, ~0.5 ground returns/m^2).
   Another gen1 block (different density, scan geometry, season) may shift every feature's
   scale; the standardisation is Elba's.
