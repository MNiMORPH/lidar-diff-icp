# `penetration.npy`: RETIRED 2026-09-05

Andy's call: *"Let's let penetration go and remove it."* The layer, its producer and its
library module are gone from the live pipeline as of commits 978ece3 and e0f962e. The code
is recoverable from git history; this file records what was removed and what was not.

## Why it went

It was not helping. The workflow graph settled it: the only step requiring
`penetration.npy` was `strata_core`, itself optional AND blocked. The base chain never
touched it, none of the five vegetation-correction steps touched it, and `convexity` ran
`--without penetration` deliberately. `canopy.py` already carried the verdict — *"do not use
`penetration` as a canopy measure; `canopy_cover_pfs` is the cover measure, computed
identically on every tile."*

## Removed

    src/lidar_diff_icp/canopy.py     the whole module (ground_penetration was the layer;
                                     leafon_slope_flag and inflate_lod were reachable only
                                     from their own tests, the flag having been retired
                                     2026-09-02)
    scripts/make_penetration.py      the producer
    tests/test_canopy.py             its tests
    Step("penetration")              from the graph
    Step("strata_core")              its last consumer -- already blocked, and penetration
                                     DEFINED its two classes, so it could not outlive it
    run_all_sites                    the computation, the save, --no-penetration

## NOT removed, and the distinction matters

**The seven `data/derived/*/penetration.npy` files (17 MB) are still on disk.**
`data/derived` is gitignored, so unlike the code above they were **never in git** — deleting
them is not reversible by checkout, only by re-running the producer recovered from c616e68.

**37 analysis scripts still load the layer by path.** They are the historical beam-angle and
strata investigation, not the settled method. They will fail if the files go.

Both are open questions for Andy, deliberately left rather than decided by a cleanup pass.

## What the layer actually measured, for the record

The quantity was **geometry-confounded** and a poor canopy proxy. As recorded in the repo,
ground-return fraction correlates **-0.84 with scan angle**
(`analysis/ridgelines/gen1_intensity_fit.py:6`, which names it as the source of the 0.6
bimodal split, and `gen1_save_angles_slope.py:90`).

CORRECTED 2026-09-05: an earlier version of this file added "-0.91" for overlap density and
"-0.33" for canopy, and said both were stated in canopy.py and AUDIT_findings.md. Neither
number appears anywhere in this repository, and neither file makes that statement. They came
from a session memory note and were written here as if cited. Only the -0.84 is supported.

## The defect, measured 2026-09-05

`data/derived/elba_fulldensity/penetration.npy` predates the producer and carries a
zero-fill where the function returns NaN:

| tile | cells | exactly 0.0 | NaN |
|---|---|---|---|
| elba | 355,600 | 206 | 677 |
| elba_fulldensity | 355,600 | 883 | **0** |

The two layers are computed from the same cloud over the same grid, and this is verified,
not assumed:

* where both are finite — 354,923 cells — they are **identical, max abs difference 0.000e+00**
* the 677 cells that are NaN in `elba` are **all exactly 0.0** in `elba_fulldensity`
* 677 + 206 = 883, exactly the zero count in `elba_fulldensity`

So `elba_fulldensity`'s 883 zeros are 206 genuine (cells with returns but no ground returns,
where 0.0 is the right answer) plus **677 cells that have NO gen2 returns at all and should
be NaN**.

**Consequence.** Every consumer cuts forest at `pen < 0.25`. Those 677 cells — 0.190% of the
grid — therefore read as *maximally closed canopy* rather than as unmeasured. An absent
measurement is being used as evidence of dense forest, which is the same error as
`nan_to_num` on an unmeasurable variance: the most extreme available value assigned exactly
where there is no information.

`elba_fulldensity` is the tile roughly twenty analysis scripts hardcode by path.

## The zero-fill is now moot

The fix was one command against a layer nothing live consumes. With the layer retired it is
not worth running: the 677 cells matter only to the historical scripts above, and only if
those are ever revived.
