# Open items

Each entry states the problem, what has been measured, and what remains.
Numbers here are pasted from a run, with the script that produced them.

---

## gen1 may not resolve ground under dense grass

**The problem.** gen1's per-cell ground elevation is the MEDIAN of its CSF-classified
ground returns. Under a dense sward, returns that stop in the grass enter that class and
the median floats above true ground. gen2 has roughly 12x the ground-return density
(5.78 vs 0.49 pts/m^2), so far more of its shots reach the real surface through the same
grass. The DoD is `gen2 - gen1`, so the result reads as **false erosion**.

This is a SAMPLING argument, not a leaf-state one. It is independent of, and opposite in
sign to, the leaf-on canopy bias that lifts gen2 and reads as false deposition, which the
DoD already corrects via `q2(cover)` (`analysis/ridgelines/q2_cover_fit.py`). Two epochs,
two mechanisms, opposite signs.

**What is measured** (`analysis/gen1_grass_lift.py`, commit `6e41db1`; gen1-INTERNAL,
`p50 - p10` of each cell's own gen1 ground returns from `gen1_csf_angles.npz`):

```
  group                        ncell   p50-p10   slope   cover    DoD mm
  floodplain, cover<0.15       38894      96.0     5.4   0.010     +18.8
  upland,     cover<0.15       99555      60.2     3.8   0.000      -5.2

  floodplain EROSION cells  n=12,789  lift 220.0 mm  slope 10.0  cover 0.170  DoD -188.6
  floodplain, not eroding  n=114,724  lift 106.7 mm  slope 13.8  cover 0.230  DoD  +28.5
```

Open floodplain against open upland: `96.0` vs `60.2 mm`, at near-zero canopy in both and
similar low slope -- meadow sward against cropland. The cells detected as erosion carry
twice the lift while being FLATTER and LESS canopied, so neither intra-cell relief nor
leaf-on explains them; both would push the other way. The lift, `220 mm`, is about the
size of the apparent erosion, `-189 mm`.

Independent support: Google Earth shows the eroding floodplain as open meadow in
September 2008 and a dense sapling thicket in April 2021.

Per-cell lift grid: `data/derived/elba_fulldensity/gen1_lift_p50_p10.npy`.

**Two ways forward, and the choice is open.**

1. *Correct gen1.* Replace the per-cell median with a lower statistic. **`p10` is not
   proposed** -- it was used to measure the spread and is an unjustified estimator.
   A defensible replacement needs either calibration against surveyed 2008 control on
   grass, or a floor derived from return count and footprint. Whatever is chosen must be
   applied the same way to both epochs, or its effect on gen2 measured, or it manufactures
   a difference.

2. *Blank the affected cells.* Mask cells whose gen1 lift exceeds a threshold and report
   them as no-data rather than as change. Cheaper, honest, and loses coverage exactly
   where the grass is thickest. The threshold is still a choice that has to be justified.

**Caution before any tile-wide version.** The lift also grows with intra-cell relief. The
flat floodplain is a clean control precisely because relief is low there; a global
correction needs the slope dependence separated first.

**Still to do.** Confirm the high-lift cells coincide with the 2008 open meadow; choose
between correction and masking; re-run the floodplain DoD against whichever is chosen.

---

## Generalization: what stops a second tile running end to end

**State as of 2026-09-01, after the two producers were parameterized.** elbaext now has
`floodplain_mask.npy` and `crest_mask.npy`; it is still missing `penetration.npy`,
`canopy_struct.npz`, `core_forest.npy` and `core_open.npy`.

### 1. DONE: both producers are parameterized by tile

`convexity_dod_landcover.py` (commit `94a3833`) and `strata_core.py` (commit `3f051d3`).
Grid origin and resolution are read from each tile's own `corrections.json`; every output
path is built from `--tile`. Both reproduce elba **byte-identically** (`floodplain_mask`,
`crest_mask`, `kappa_L10/20/30`, `core_forest`, `core_open`; `ridgecrest_pixels.npz` matches
in column order and in every column but `dod_m`, which differs because
`elba_refdatum/dod_geoid.npy` was regenerated after that table was written).

**The blocking gate is now clear.** `reference_cells` runs on both tiles: 28,533 cells at
elba_fulldensity, 78,270 at elbaext.

Two things this turned up, both recorded and neither acted on:

- **`curvature_diffusion.py` had the same defect; fixed** in commit `86beb3b`. `--dod`,
  `--gen1-date` and `--gen2-date` are now explicit (a wrong dt rescales K directly), the
  silent `sys.exit(0)` states what it lacks, and the in-place augmentation of
  `ridgecrest_pixels.npz` with `curv_xx`/`curv_yy`/`curv_laplacian` is documented as an
  ORDERING requirement: run it AFTER the convexity producer or the three columns vanish.

### The DoD grids are NOT on one frame -- OPEN

Checked because the elbaext run needed a `--dod` and I picked one. The pairing I used is
wrong, and not in the way expected: `elbaext/dod_geoid.npy` is fine; the ELBA grid it was
compared against is the odd one out. What each shipped grid actually applied, read from its
own `corrections*.json`:

    grid                             method            const   tilt mm/km      swath_tie
    elba_refdatum/dod_geoid.npy      geoid, HARDCODED  67.000  +0.610 -0.730   None
    elba_fulldensity/dod.npy         geoid, auto       67.281  +0.778 -0.568   intercept
    elbaext/dod_geoid.npy            geoid, elba plane 66.701  +0.778 -0.568   intercept
    elbaext/dod.npy                  reference_plane  -84.930  +6.940 +3.540   None

`elbaext/dod_geoid.npy` was built by `elbaext_geoid_regrid.py` to share ONE frame with elba:
it fits the geoid plane on ELBA's bounds and re-expresses it about elbaext's centroid, which
is why the tilt matches `elba_fulldensity/dod.npy` exactly and the constant differs by the
0.580 mm the re-centring implies. Both carry `swath_tie: intercept`.

`elba_refdatum/dod_geoid.npy` is the outlier: `run_elba_dod.py`'s own docstring says it
SUPERSEDES it, for passing a hardcoded geoid triple. Recomputing the planes now,
`geoid_difference` gives elba (+0.778, -0.568) and elbaext (+0.353, -1.031) mm/km, so the
hardcoded (+0.610, -0.730) matches NEITHER tile's own-bounds fit.

**Consequence:** the two crest tables produced on 2026-09-01 are not comparable across tiles.
The same-method pair is `elba_fulldensity/dod.npy` with `elbaext/dod_geoid.npy`; better still,
once #19 lands, both tiles' `dod_cover_q2`. The elba crest run was left on the superseded grid
deliberately, because reproducing it byte-identically was the check. **Andy's call which grid
the crest suite should read from here.**

### 2. STILL OPEN: retire `penetration.npy`, or give it a producer

**How it was produced -- established 2026-09-01, by reproducing it.** No tracked code writes
the file, but `src/lidar_diff_icp/canopy.py:19` `ground_penetration()` does compute it, and it
is the producer. Recomputing that function's arithmetic on
`data/after/3dep2021_fulldensity.laz` (182,923,322 points) over elba's bounds, res 5.0,
`ground_class=2`, `noise_class=7`, accumulating the bincounts in chunks so the computation is
the same one rather than an approximation of it:

    cells finite in BOTH                     354,923
    of those, cells that DIFFER                    0
    cells the recomputation calls NaN            677   <- the whole disagreement

So the file is that function's output, and it is regenerable. Timing agrees: the file is
dated 2026-08-21 13:13 and `canopy.py` was committed at 13:58 the same day (`dc3c667`), i.e.
it was run interactively before the function was committed. `scripts/run_all_sites.py:135`
calls `ground_penetration` but persists only the derived `leafon_flag.npy`, which is why no
write of `penetration.npy` appears in the source.

**The one deviation is a no-data-as-zero substitution.** The function returns NaN where a cell
has no returns at all. The stored file holds `0.0` in all 677 such cells -- a single distinct
value, so it was filled, not measured. Under the strata cut used across the repo
(forest `pen < 0.25`), **all 677 are classified forest**: cells with no gen2 returns read as
maximally closed canopy. They do NOT reach `core_forest` (0 of 677 survive the purity and
cluster filters), so the contamination is confined to the raw `forest0` mask and anything
built directly on `pen < 0.25`.

This is the same defect class as the ones fixed this session, in the data rather than the
code. Regenerating the file from `canopy.ground_penetration` would fix it by construction,
because the NaN is already what the function returns.

**`canopy_struct.npz` is a second orphan of the same kind** — no producer either, read by six
files, present only for elba. `strata_core.py` needs BOTH, so it is blocked twice over.

What is still gated on the decision: the forest/open crest split (step 4 of the convexity
producer) and `core_forest`/`core_open`. Nothing else.

**Recommendation: retire it.** `analysis/ridgelines/AUDIT_findings.md` flags `penetration.npy`
as a gen2-derived variable contaminating gen1-internal conclusions; it is referenced by 34
tracked files. `canopy_cover_pfs.npy`, `forest_pfs.npy` and `open_pfs.npy` already exist for
BOTH tiles and are the cover measure meant to replace it. The cost is not free: swapping the
crest split to the PFS layers CHANGES elba's step-4 answer, so it needs the two cover
definitions measured against each other on elba first, not substituted. Andy's call.

### 3. Then produce elbaext's cover-corrected products

`dod_cover_corrected.py` and `lod_cover_q2.py` are ALREADY tile-parameterized; they have
simply never been run there. Two commands. The q2 slope must be elbaext's OWN fit, because
the relation is per-site (it depends on each pair's phenology and undergrowth).

### 4. Loose end, not blocking

`nearground_cells_sn.npz`, `gen1_canopy_frac.npz` and `ridge_mask.npy` showed no producer in
the scan, but the scan cannot see filenames built from flags (`nearground_cells.py --out
..._sn.npz` is invisible to it). All exist for BOTH tiles, so they do not block; the
producers just are not identified. `curv_laplacian.npy` IS identified —
`curvature_diffusion.py` writes it, along with `curv_xx` and `curv_yy`.

### The rule these now follow

A missing input REFUSES; it does not run differently. Running without an optional layer is
allowed and often right, but must be stated -- `use_floodplain_mask=False`,
`--without penetration,core_forest`. Silent zero-fill turned an ABSENT layer into a
MEASURED EMPTY one, which is how a table came to look like a finding of "no forest", and how
two tiles' reference populations came to differ by 39,038 cells unnoticed.

---

## The canopy k is now read, not typed -- and its form is wrong

**Closed 2026-09-01** (`c99c141`, `97798e4`). `dod_cover_attribution.py` carried
`k = 48.4 if TILE == "elbaext" else 49.6`; neither number existed anywhere else in the repo,
both were read off a run and retyped, and the `else` branch handed every other tile Elba's
value silently. `cover_offset_reference.py` now writes its coefficients as JSON and the
consumer reads them, prints the population behind the fit, and refuses when it is absent.

Re-running the calibration on the current registration:

    tile               typed     read now
    elba_fulldensity    49.6       48.73
    elbaext             48.4       48.62

The typed pair differed by 1.2 mm between tiles; measured, they differ by 0.11 -- the
tile-to-tile difference that pair implied was an artefact of when each was typed.

**Still open, and larger.** Neither tile's own model selection picks the LINEAR form that the
additive canopy model `pred = k * cover` requires:

    elba_fulldensity   selects quadratic      AIC 108.7  vs linear 372.6
    elbaext            selects optical depth  AIC  86.8  vs linear 152.5

The run now says so whenever the selected form is not linear, instead of using a coefficient
the data does not favour without comment. Whether the attribution should move to the selected
form -- and it is a DIFFERENT form on each tile, which is itself worth explaining before
adopting either -- is not decided here.

## Can the WHOLE pipeline run on a new region? Audited 2026-09-01

**Core DoD: yes, already.** `scripts/run_all_sites.py` carries a `SITES` dict and runs
`difference_dem` + `detect_change_standard` per site, persisting every output.

**Optional pieces: all parameterized now, but not wired together.** Every producer takes a
tile: `trace_ridgelines.py` (positional `tile_dir`, `--out ridge_mask.npy` -- this identifies
the "unknown producer" of `ridge_mask.npy` listed below as a loose end),
`forest_metrics_pfs.py` (`tile_dir`, `after_laz`), `beam_offset_table.py` and
`gen1_save_angles_slope.py` (positional), `convexity_dod_landcover.py`, `strata_core.py`,
`curvature_diffusion.py`, `cover_offset_reference.py`, `dod_cover_attribution.py`,
`q2_cover_fit.py`, `dod_cover_corrected.py`, `lod_cover_q2.py`, `scripts/make_penetration.py`.

### DONE 2026-09-01 (1), (2), (3) and the packaging question

* **(1) `08e7513`** -- the q2 slope is read from the tile's own `q2_cover_fit.json`, which
  `q2_cover_fit.py` now writes, and refuses when absent. Three values were in play:
  `-0.1922` typed in the code, `-0.1835` in the shipped product, `-0.1792` from today's
  refit -- and the JSON records its inputs' mtimes, so the consumer warns that the refit
  itself read a `beam_offset_table` older than `corrections.json`.
* **(2) canopy_struct.npz was NOT an orphan.** `analysis/ridgelines/canopy_struct.py`
  produces it; my write-form search missed it because the path is a module constant. It was
  hardcoded to elba. Parameterized; all six fields byte-identical on elba.
* **(3) `97c0001`** -- `src/lidar_diff_icp/workflow.py` declares the 15-step graph and
  derives the order. `--check` reports MISSING / STALE / OK per step; `--plan` prints the
  commands in order. It runs nothing.
* **Packaging (`977352c`)** -- four runtime deps the LIBRARY imports were undeclared
  (`pyarrow`, `requests`, `pyshp`, `shapely`), so a fresh `pip install` produced a package
  that failed on import. `testpaths = ["tests"]` turns a bare `pytest` from a 25-minute
  abort into 8 s. `lidar-diff-workflow` is a console script.

### The chain is STALE on both tiles, and that is the finding

`lidar-diff-workflow --check` on elba_fulldensity reports `beam_table` STALE against
`corrections.json`, and everything downstream of it inherits that. So the adopted
`dod_cover_q2` rests on registration that has since been superseded. Regenerating it is a
decision, not a chore, and is NOT taken here.

### What still stops a new region

1. ~~**`dod_cover_corrected.py --slope` defaults to `-0.1922`**~~ FIXED, see (1) above.
   Original note kept for the reasoning: -- Elba's own q2 slope, as a
   DEFAULT. Run it on another region without noticing and it applies ELBA's cover correction
   and writes a `dod_cover_q2` that looks finished. This is worse than the `k = 49.6` case
   that was fixed today, because the number IS the correction rather than an attribution
   term. It must be required, or read from that tile's own `q2cover.fit_tile`. `--gen2`
   likewise defaults to Elba's cloud.
2. **`canopy_struct.npz` has no producer** -- blocks `strata_core.py` on any new tile, the
   same way `penetration.npy` did until today.
3. **No driver runs the optional chain in order.** The dependencies are real --
   ridge_mask -> convexity -> floodplain/crest; PFS cover -> beam_offset_table -> q2 fit ->
   dod_cover_corrected -> lod_cover_q2 -- and nothing encodes them, so the order lives only
   in whoever remembers it.
4. **The shared vertical frame is a per-region decision, not a rule.** elba came from
   `run_elba_dod.py` (geoid auto-computed on its own bounds); elbaext from
   `elbaext_geoid_regrid.py`, which deliberately fits the plane on ELBA's bounds and
   re-expresses it, because "both auto-computed" is the same METHOD but not the same FRAME.
   A third region needs that choice made explicitly; there is no generic answer.
5. **gen1's absolute datum statewide** is unresolved (task #11).

### And there is no reproduce-elba test

Every reproduction this session was verified by hand, artifact by artifact. Nothing runs the
chain and asserts it still lands on elba's shipped products, so the "does it still reproduce
elba" half of the question has no automated answer. An ordered driver (3) would give that
test something to invoke.

## First end-to-end run on a new site: Carlton, 2026-09-02

Andy's instruction: run the chain on Carlton, EXCLUDING the forest-elevation adjustment,
which was fitted to gen2 undergrowth at Elba and is phenology-specific. So `q2_fit`,
`dod_cover`, `lod_cover` and `cover_calibration` were out of scope, and with them
`gen1_angles` and `beam_table`, which in this graph feed only those -- and which could not
have run regardless, since `data/csf_cache/carlton.las` does not exist.

**Five steps, all succeeded, every script unmodified on a tile it had never seen:**

    slope        698 x 484 at 5 m, 12,455 gap cells filled; slope median 6.43 deg
    ridge_mask   99,097 ridge cells, 49% on highs (threshold 200, the script's default)
    convexity    17,769 crest cells of 99,097 divide cells; 36,269 divide cells floodplain
    curvature    crest Laplacian median -0.0112 1/m
    pfs_cover    337,310 of 337,832 cells carry a cover value

Three optional pieces declined and each NAMED what it lacked rather than substituting: the
ridge tracer's furrow/forest QC, the convexity producer's step 4, and the curvature
producer's diffusion part.

### What the run found that reasoning had not

1. **`slope.npy` was not a base input** (`61cb160`). I had declared it as one; nothing
   produced it but a script hardcoded to elba, which is why carlton held every other base
   product and none of the chain could start. Now `scripts/make_slope.py`, bit-identical to
   the shipped elba file.
2. **`curvature` reported STALE the moment it succeeded** (`3523732`) -- on carlton AND on
   elba. It augments `ridgecrest_pixels.npz` in place, so declaring that file a REQUIREMENT
   made its own write invalidate it. `Step.mutates` now models in-place augmentation.

Both were defects in the graph declared the day before, and neither would have surfaced
without running it.

### Two measurements to carry forward, NOT interpreted

**CORRECTED.** I first reported "carlton has no cross-epoch datum at all". That was wrong,
and it was a reading error: I called `.get("cross_epoch_datum", {}).get("method")`, got
`None`, and read absence of that KEY as absence of a datum. Carlton carries one, under a
different key -- `cross_epoch_tie_order2_coef`, the retired order-2 PARABOLA.

**Only one derived tile is on the current datum method** (`analysis/datum_method_audit.py`):

    elba_fulldensity    geoid: geoid_difference     <- the only one
    elba_refdatum       geoid: reference_plane      (also retired)
    elbaext             geoid: reference_plane      (also retired)
    elba, final, carlton, carlton_density, cook, mnrv, ne,
    whitewater, battlecreek                         order-2 PARABOLA (retired)

Retiring a method from the code did not retire it from the products. Measured on carlton,
its stored coefficients evaluated against the geoid plane over its own grid:

    parabola (applied)   median  +93.82   ptp  51.77 mm
    geoid plane          median  +78.10   ptp   5.17 mm
    parabola - geoid     median  +15.91   range -19.89 .. +35.48 mm

The applied datum imposes ~52 mm of spatial structure where the physical geoid difference
has 5 mm, against carlton's `stable_1sigma` of 44.4 mm. That is precisely the
fitted-surface-absorbs-real-signal failure the parabola was retired for, still in force.

**The crest DoD is opposite in sign to Elba's at steep slopes:**

    carlton  all crests n=16,927  median -3.6 mm   15-90 deg: -27.7 mm
    elba     all crests n=12,459  median +2.4 mm   15-90 deg:  +8.7 mm

The two are NOT comparable as they stand, and the reason is now specific: they are on
DIFFERENT DATUM METHODS -- carlton on the retired parabola, elba_fulldensity on the geoid --
and carlton's parabola alone carries 52 mm of spatial structure. Add to that different
acquisitions, and no cover correction on carlton by instruction. Do not read this as a site
difference.

**The cover thresholds classify little at Carlton, and that is EXPECTED, not a problem**
(Andy, 2026-09-02): they are Elba-specific, and nothing at Carlton splits or corrects on
them -- the cover grid is a measured layer here, nothing more.

    forest (cover >= 0.5)    13,575    4.0%
    open   (cover <= 0.1)   108,082   32.0%
    NEITHER                 215,653   63.9%

Recorded as the measurement it is. My first note called for "attention to the thresholds",
which was drift toward the cover-based adjustment that was explicitly excluded from this run.

## Six sites re-run onto the geoid datum, 2026-09-02

Andy's instruction after the Carlton run exposed that only `elba_fulldensity` was on the
current datum. Snapshot first: `data/derived/_parabola_era_snapshot/` (383 MB, all six,
with a README recording why the old products cannot simply be corrected).

**Why re-running was the only option.** The parabola fits `dx`, `dy` AND `dz`, and the
horizontal part was applied by RESAMPLING, so it is not invertible on the grid. Measured on
carlton: horizontal field median `|dxy|` 103 mm, max 304 mm, which on that tile's own slopes
implies a vertical error of median 10.4 / p90 48.4 / max 219.2 mm -- against a
`stable_1sigma` of 44.4 mm. Undoing it costs about the whole detection limit.

**No new flags were needed.** `difference_dem` no longer implements the parabola and
`run_all_sites.py` passes none of the datum arguments, so the re-run took the current
defaults: `tie="reference"`, `geoid_datum=None` (auto per site), `swath_tie="intercept"`,
`along_track_drift=True`. Postcondition checked with `analysis/datum_method_audit.py`: all
six now record `geoid: geoid_difference` with `swath_tie: intercept`.

**What changed (new minus snapshot, on cells finite in both):**

    site           n both    median    NMAD      p1      p99   sigma old  sigma new
    battlecreek    20,375     -10.5    34.3  -293.8    297.5      0.0322     0.0305
    carlton       305,799      -1.3    46.8  -271.7    561.2      0.0444     0.0391
    cook          287,900     -23.5    75.4  -598.4    375.9      0.0705     0.0828
    mnrv          242,935     -42.8    58.8  -757.1    640.7      0.0606     0.0652
    whitewater    282,959     -28.6    64.8  -321.6    406.0      0.0804     0.0863
    elba          339,829       1.2    52.3  -469.5    312.5      0.0599     0.0550

**CORRECTION -- that table is NOT the datum change alone.** The old products also predate
`ed8ab82` (2026-08-20), which made the MEDIAN the default ground estimator. Read from the
products themselves: carlton's snapshot records `ground_percentile=0.1`, the new run records
`0.5`. So the re-run changed TWO things at once, and the difference above is their combined
effect:

    1. cross-epoch datum:  parabola -> geoid
    2. ground estimator:   p10 -> median (ground_q 0.1 -> 0.5), BOTH epochs

The second is much the larger. Carlton's gen2 ground surface alone moved by

    z_after new - old:  n=325,377   median +244.69 mm   NMAD 258.43 mm   max|d| 6157.7 mm

which is why the DoD moved only tens of millimetres: two large shifts, one per epoch, that
largely cancel in the difference. Do not attribute the DoD table to the datum. Separating
them would need a re-run holding `ground_q=0.1`, which is not worth doing -- the median is
the adopted estimator and the geoid the adopted datum, so the new products are the right
ones on both counts.

Stable-ground sigma fell at three sites and rose at three. INTERPRETATION, flagged as such:
a rise is what removing a FITTED surface should produce, because the parabola was fitted to
minimise residuals on stable ground and the geoid cannot flatter itself that way. A lower
sigma under the parabola was never evidence it was right. The three sites where sigma
IMPROVED are not explained by that argument and are not explained here.

**Still on a retired method, and NOT re-run** (not in the `SITES` registry): `carlton_density`
and `final` and `ne` on the parabola; `elba_refdatum` and `elbaext` on `reference_plane`.
`elbaext` matters most -- it is the elba-overlapping tile and shares elba's frame by
construction.

## Canopy cover is OPT-IN, not part of a standard run

**Decided 2026-09-02 (Andy).** "In future pipelines it should be opt-in, depending on
whether a cover correction is needed."

The trigger: I ran `pfs_cover` at all six sites because I pictured a uniform layer set,
not because anything needed it, and then reported the resulting 29-column beam tables as a
CAVEAT -- as though the absent `canopy_cover` / `pfs_forest` / `pfs_open` columns were a gap
to close. They are the optional-column mechanism working as designed.

Checked, not assumed: `slope`, `ridge_mask`, `convexity` and `curvature` never read
`canopy_cover_pfs`; `gen1_angles` and `beam_table` carry it as a COLUMN only and both ran
clean without it at four sites. Its only real consumers are `q2_fit`, `dod_cover` and
`cover_calibration` -- the cover-adjustment family -- where it is definitional.

It is not free: a large tile needs an `untwine` COPC first, and cook's took ~14 GB of
scratch for a 0.98 GB output.

**So: build cover when a cover correction is wanted at that site, and not otherwise.** A
tile without it is complete, not deficient.

# Decisions

Closed questions, with the reasoning and the conditions that would reopen them.

## gen2's horizontal accuracy is not checked, and does not need to be

**Decided 2026-08-31 (Andy).** Closes next-action 3 of
`ground_control/HANDOFF_FROM_GROUND_CONTROL.md`, which listed gen2's unchecked horizontal
accuracy as an open item on the grounds that step 4 registers gen1 *to* gen2, so a lateral
error in gen2 propagates into gen1 invisibly and becomes a vertical error on a slope.

**Why it is closed, and it is not "3DEP is probably fine".** The DoD is invariant to it.
Registration applies a Nuth-Kaab lateral shift that moves gen1 onto gen2, so a UNIFORM
horizontal error in gen2 is absorbed: both epochs end up in the same frame and the
difference of two clouds that agree with each other does not depend on where that shared
frame sits. The quantity we compute cannot see the error.

**Where it does not cancel** is the tie to control, because the marks carry independently
surveyed horizontal coordinates, so a lateral error samples gen2 at the wrong place. The
geometry makes that cheap, because control is sited flat by design:

```
  slope at the 389 gen2 control marks (deg): p50=2.6  p90=6.5  p99=15.5

  vertical error produced by a lateral error, at those marks (mm)
     lateral      p50      p90      p99
       10 cm      4.5     11.4     27.8
       20 cm      9.0     22.9     55.6
       50 cm     22.6     57.2    139.0

  for scale: the gen1 Elba datum SE is +/- 23.38 mm
```

A 20 cm lateral error costs 9 mm at the median mark. It also averages down over 389 marks,
because a fixed lateral shift produces vertical errors whose sign follows local aspect, so
with mixed aspects it behaves as noise rather than as an offset.

**Two conditions that would reopen it.**

1. *Non-uniform lateral error.* The cancellation assumes one shift fits the whole tile. A
   per-swath or along-track lateral drift would not be absorbed, and this acquisition has
   already required a per-swath VERTICAL tie, so a lateral analogue is not far-fetched. It
   would appear as edge-of-swath artefacts, not as a whole-tile offset.
2. *Comparing gen2 to anything it has not been registered to.* The flatness that makes the
   error cheap holds at the marks, NOT on the tile: Elba slope is p50 9.2 deg, p90 28.6, so
   a 20 cm lateral error there is 32.5 mm at the median and 108.9 mm at p90. Those numbers
   are harmless only while the error cancels.

Numbers from `analysis/control_lowveg_offset.py` mark slopes and
`data/derived/elba_fulldensity/slope.npy`.
