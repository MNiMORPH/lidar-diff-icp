# FRAME — the current state. Read this first; the dated FRAME_*.md files are superseded.

Written 2026-09-06. Verify every structural claim against git and the files before relying
on it: this is a map of ROLES and NEXT ACTIONS, not a results log.

    ./lidar-icp/bin/python -m pytest -q        expect 431 passed
    git rev-list --count origin/main..HEAD     expect ~105, all UNPUSHED
    lidar-diff-workflow --tile data/derived/elba --check

## The pipeline is six steps

    completeness -> base -> slope -> ridge_mask -> convexity -> curvature

`base` is `pipeline.difference_dem` via `run_all_sites`, and it produces the shipped
product: `dod.npy`, `lod.npy`, `z_after.npy`, `change.npy`, the three GeoTIFFs,
`corrections.json`, `regions.json`. Everything else is downstream.

**It takes no cloud arguments** — `base` gets them from the Site record in `sites.py`.

`lidar-diff-workflow --tile <dir> --run` executes it, stopping at the first failure, and
refuses before running anything if a selected step lacks an argument. `--only`, `--force`,
`--dry-run`.

## Alongside, reached over to, not part of it

    --with vegetation_correction     10 steps: class2_spread, q2_fit, dod_cover, lod_cover,
                                     and their feeders pfs_cover, gen1_angles, beam_table,
                                     nearground, nearground_split, canopy_struct

The correction is MEASURED AND NOT ADOPTED: on open ground the calibrated curve was WORSE
than the median (RMS 52.5 vs 49.1 mm, held out on 227 NVA marks). Its outputs are
`dod_cover_q2.npy` / `lod_cover_q2.npy`, never `dod.npy` / `lod.npy`.

The split was measured, not judged: a step is alongside when its products are read only by
`q2cover.py`, the correction's own library. `slope`/`ridge_mask`/`convexity`/`curvature`
look similar but are read by `refcells.py`, which 20+ scripts use, so they are pipeline.

`analysis/tools/` holds nine hand-run tools, declared in `workflow.TOOLS` with each
script's own docstring line.

## The swath alignment — see analysis/SWATH_ALIGNMENT_METHOD.md

Eight steps, and step 7 is the one that is easy to miss: the free network deliberately
leaves the absolute level free, and **the ground-control datum is what fixes it**, exactly
(44.60 mm of zero-line spread -> < 1e-9).

Since 2026-09-05 the network is weighted by **1/variance, not by cell count**. That changed
a shipped product: Battle Creek's swath 1102 moved +13.7 -> +27.1 mm, because its only edge
was a 36-cell sliver with a 43:1 lever arm returning -3.4640 m, which under `w = n` entered
at weight zero and left the swath at an `lstsq` minimum-norm value.

`across_track_tie` (the intercept tie) is OURS: the regression is USGS/ASPRS standard
(Sampath et al. 2016) but they constrain the intercept to zero and report the slope. See
`analysis/ACROSS_TRACK_TIE_LITERATURE.md`. **Open**: our single covariate estimates a
COMMON-MODE roll only; regressing on `tan θ_ref` and `tan θ_src` separately would test
whether `k` carries a differential bias.

## Sites

All six rebuilt on current code. Five are byte-identical across the reorganization
(dod/lod/z_after), which is the evidence that moving 55 files changed no science.

    site         stable_sigma  median LoD   gen2 completeness
    elba            0.055 m      0.113 m       0.978
    whitewater      0.055        0.105         0.971
    carlton         0.041        0.078         0.989
    cook            0.088        0.160         0.936
    mnrv            0.056        0.103         UNMEASURABLE (bbox spans two 3DEP projects)
    battlecreek     0.044        0.090         0.998

## Open, and needing Andy

1. **Battle Creek's valley top.** Its histogram cut removes 72.4% of the grid (283.6 m) in
   a built environment where graded lots set the modal elevation. Its `stable_sigma` is
   computed on that set. He must state an elevation or a fraction ceiling.
2. **A 1-cell instability at Battle Creek, cause NOT established.** Two runs gave 4,841 vs
   4,842 stable cells at the same valley top with identical `z_after`, moving the DoD by
   <= 0.883 mm through the drift fit. Re-running twice more is byte-identical, so the
   pipeline IS deterministic given identical inputs — something in the input differed and I
   did not find it. Sub-mm, but unexplained.
3. **The deletion candidates** — `analysis/DELETION_CANDIDATES.md`, 26 scripts, nothing
   deleted. The eleven dated 2026-09-01 or later are FALSE POSITIVES and are named there.
4. **105 unpushed commits.**

## Guardrails

* **PLOT WHAT YOU MAKE**, same turn, unprompted.
* **A claim that would change course must be verified in the same turn it is made.**
* **When Andy points at prior work, LIST THE FILES.** A repeated request is evidence the
  premise is wrong.
* **No invented thresholds.** Say it in one sentence and let him decide.
* **A regression test must be shown to FAIL without the fix.** Two vacuous tests were caught
  this way on 2026-09-05; both breaks passed before the scene was fixed.
* **Numbers in prose are pasted from command output.** Four claims went out unverified on
  2026-09-05-06 — mnrv's 0.67, two correlations, and a `swath_tie` defect that was my own
  code. All four were one command away.
* **Evidence of use hides from static scans.** Three times: scripts that COMPUTE a quantity
  rather than loading it; documents naming a file without its path; work recorded only in a
  commit message. Inventory by function, import AND artifact, and read the result.
* Shared laptop: one heavy job at a time, watch `free -h`.
* Commit granularly; never push/tag/release without explicit current-message authorization.
