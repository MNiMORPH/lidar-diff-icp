# `penetration.npy`: it has a producer, and one shipped copy still carries a defect

Resolves the standing question "retire it, or give it a producer". **Neither retire nor
build: it already has one**, and it is far too widely consumed to retire.

    producer   scripts/make_penetration.py          (declared as Step("penetration"))
    consumers  34 tracked files
    tiles      battlecreek, carlton, cook, elba, elba_fulldensity, mnrv, whitewater

The caveat that matters is not about the producer but about the QUANTITY: penetration is
**geometry-confounded** and is a poor canopy proxy. `canopy.py` states the verdict
directly — *"do not use `penetration` as a canopy measure; `canopy_cover_pfs` is the cover
measure, computed identically on every tile"* — and `AUDIT_findings.md` puts the layer in
Tier 1 as a gen2-derived variable, warning against binning any gen1 quantity against it.

The measurement behind that, as recorded in the repo: ground-return fraction correlates
**-0.84 with scan angle** (`analysis/ridgelines/gen1_intensity_fit.py:6`, which names it as
the source of the 0.6 bimodal split, and `gen1_save_angles_slope.py:90`).

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

## The fix, and why it is not applied here

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/make_penetration.py \
        --tile elba_fulldensity --after data/after/3dep2021_fulldensity.laz

The producer writes NaN by construction, so this yields exactly `elba`'s array.

It is NOT run unprompted because it is a data product, not code. But the blast radius is
SMALLER than first written here, and the workflow graph is what settles it: the only step
that requires `penetration.npy` is `strata_core`, which is optional AND blocked
(`canopy_struct.npz` has no producer). The base chain does not touch it, none of the five
vegetation-correction steps touch it, and `convexity` runs `--without penetration`
deliberately.

So regenerating it would invalidate the HISTORICAL beam-angle and strata analyses that load
the layer by path — `geology_forest_split`, `curvature_diffusion`, `run_steady_state_strata`
and the rest — and nothing the current DoD or the current correction depends on. That makes
it low priority rather than a decision blocking anything.

The magnitude is small (0.19% of cells) and the direction is known (677 cells leave the
forest class). Whether that is worth a re-run is exactly the judgement to be made out loud.
