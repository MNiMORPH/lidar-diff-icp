# Next session: the prepared plan

Written 2026-09-06 at the end of a long session. `analysis/FRAME.md` is the state; this is
the queue. Verify both against git and the files before acting — that is the standing rule
and it caught four wrong claims of mine in the last two days.

## 1. Finish the last 7 moves — READY, blocked on a one-line fix per test

`scripts/reorganize_analysis.py --plan` wants 7 files in 4 clusters. Six go to
`analysis/lib/`, and three clusters were REVERTED by the gate on the last attempt. The
cause is known and small: each test reaches its subject with a `sys.path` insert naming the
subject's OLD directory.

    tests/test_crossline_fit.py:21   ".." , "analysis"               -> needs "analysis", "lib"
    tests/test_csf_tiled.py:19       "..", "analysis", "slope_bias"  -> needs "analysis", "lib"
    tests/test_forest_metrics.py:11  "..", "analysis"                -> needs "analysis", "lib"

Do it in this order, or the gate will keep reverting: **fix the three inserts and the move
in the SAME commit**, since neither half passes alone. Then

    ./lidar-icp/bin/python scripts/reorganize_analysis.py --apply

The fourth cluster, `analysis/groundtruth/gen1_datum_at_site.py`, is depth-blocked in a
non-`sys.path` way and `groundtruth/` is already a coherent home. Leaving it is defensible;
moving it does nothing for the pipeline.

## 2. The decisions waiting for Andy

* **The deletion list** — `analysis/DELETION_CANDIDATES.md`, 26 scripts, nothing deleted.
  The 11 dated 2026-09-01 or later are FALSE POSITIVES and are named individually there.
  The 14 dated 2026-08-2x are the real candidates: refuted slope-bias and beam-geometry
  work, the same class as the 33 penetration scripts already removed.
* **105 unpushed commits.** A review precedes any push and the push needs explicit
  authorization in the message that asks for it.
* **The two-covariate tie test** (below) — worth doing before this method is written up.

## 3. The open science, in priority order

**The intercept tie may carry a differential-roll bias.** From the literature search
(`analysis/ACROSS_TRACK_TIE_LITERATURE.md`): our regressor `dtan = tan θ_ref − tan θ_src`
reduces the two-line problem to a COMMON-MODE term, valid only if both lines share the same
roll error. Our own p ≈ 6e-64 heterogeneity of per-pair `c` is evidence against that. The
test is cheap and uses caches we hold: regress `dh` on `tan θ_ref` and `tan θ_src` as TWO
covariates, on elba and battlecreek, and check the conditioning before trusting the split —
within one pair they may be near-collinear, which is why the field solves it as a block.

**Battle Creek's 1-cell instability (#50).** Two runs, 4,841 vs 4,842 stable cells, same
valley top, identical `z_after`, DoD moving ≤ 0.883 mm through the drift fit. The pipeline
is deterministic given identical inputs (verified). Hypothesis, NOT established:
`terrain_masks` works from `Z21` while `valley_top_for_landscape` reads `z_after.npy` FROM
DISK — the previous run's product — which would make `base` impure on every tile using
`valley_top="histogram"` (all but elba). Andy's read is that it is likely the tests and not
worth much worry.

**Untouched all session**: #11 statewide gen1 datum, #16 grass-lift, #19 elbaext's
`dod_cover_q2`, #25 the Whitewater correction.

## 4. What NOT to redo

* **battlecreek's valley top is not a number to improve.** It has no river valley, so the
  histogram cut's concept does not apply there. Read its stable set as "low-curvature ground
  above an arbitrary elevation". Recorded in `sites.py`.
* **The vegetation correction is measured and NOT adopted.** Worse than the median on open
  ground, RMS 52.5 vs 49.1 mm. It is alongside the pipeline, reached with
  `--with vegetation_correction`.
* **`overlap_median` stays the library tie default.** `intercept` is the better estimand but
  fails unbounded on poor geometry; the pipeline opts in explicitly.

## 5. The standing fragility, so the next session does not trip on it

The reorganization rule is only as good as `META_DOCS`. Writing `DELETION_CANDIDATES.md`
made it the most specific citer of 15 scripts, and the rule wanted to move them into
`investigations/deletion_candidates/` — a folder named after a report about them, brought
into existence by writing the report. **Every project-level document written from now on
must be added to `META_DOCS`**, or it will start claiming scripts.
