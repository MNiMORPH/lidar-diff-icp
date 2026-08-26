# Pinning gen1's vertical datum to surveyed ground control

**Module:** `src/lidar_diff_icp/groundtruth/`
**Runnable examples:** `analysis/groundtruth/elba_absolute_tie.py` (the gen1 ties),
`analysis/groundtruth/gen2_checkpoint_tie.py` (gen2 against the same marks, no chain),
`analysis/groundtruth/elba_datum_constant.py` (one constant, its budget, the product),
`analysis/groundtruth/reference_swath_bias.py` (is the reference swath biasing anything?)
**Tests:** `tests/test_groundtruth_{checkpoints,tie,chain,datum}.py`
**Report:** `analysis/ABSOLUTE_BASIS_ELBA.md`

---

## 1. The problem

The pipeline registers gen1 (2008 MN DNR) to gen2 laterally, applies the
GEOID03 → GEOID18 datum shift, and aligns gen1's flight lines to one another. What it
never does is pin the *group* in z. `coreg.align_swaths` says so in its own docstring:

> Either way the group's absolute offset from another epoch must be tied separately.

So gen1's elevations at Elba float. This module supplies the missing tie, from surveyed
3DEP QA checkpoints, and it is built to be reused — the statewide goal is a deterministic
per-swath correction of the 2008 lidar, of which Elba is the pilot.

The obstacle is geometric, and no amount of data fixes it directly. gen1's lines run
**north–south** with a measured swath half-width of ~710–730 m; every checkpoint near
Elba is displaced **east–west** by 3.1–15.6 km. Lines 135–138, the ones over Elba, cover
no checkpoint and cannot be made to. But adjacent lines **overlap** — 1.31–1.50 km² per
pair, measured — so a tie taken on a line that *does* cover a checkpoint propagates,
link by link, to a line that covers the site.

---

## 2. Sign conventions

| quantity | meaning |
|---|---|
| `tie_mm` | the constant to **ADD to gen1** — already in the reference swath's frame, already geoid-shifted to the checkpoint's geoid model — to place it on the surveyed datum. `tie = surveyed − z_lidar_corrected`. |
| `Link.dz_m` | add to the **src** line to align it onto the **ref** line (`coreg.coregister_swaths(pc, swath_ref=a, swath_src=b)`). Walking a path outward from a reference: `offset[b] = offset[a] + dz`. |
| `ChainSolution.dz_total_m` | add to the far line's z to bring it into the reference line's frame. |
| `geoid_shift_m` | add to gen1's z: `N_gen1 − N_checkpoint`, from `references.geoid_difference`, computed from the PROJ grids at each mark. |

Note this is the *opposite* sign family from `dod.npy` (= gen2 − gen1) and the *same*
family as the `d_mm` per-beam table only in that both are "gen1 relative to a reference";
see `analysis/ridgelines/FRAME_2026-08-26.md` §1 before mixing them.

---

## 3. Part 1 — checkpoint ingestion (`checkpoints.py`)

Each `Checkpoint` carries its own datum as data: horizontal CRS, vertical datum, geoid
model, elevation units, and the 3DEP accuracy class (`NVA` = open ground, `VVA` = under
vegetation). Two readers:

* `read_3dep_va_shapefile` — the authoritative USGS vertical-accuracy point shapefile.
* `load_bundled` — a checked-in CSV of the six marks near Elba, for offline work
  (`groundtruth/data/`, transcribed from `analysis/ridgelines/ABSOLUTE_ELEVATION_REFS.md`
  §1a).

Two refusals, both deliberate:

* **Unknown geoid model → `UnknownDatumError`.** GEOID03 vs GEOID18 is ~67–74 mm at
  Elba, the same size as the offset a tie is trying to measure. Assuming a model would
  silently manufacture the answer.
* **Units are never converted.** The contractor `.dbf` labels `source_ele` "US Feet"
  while the values are metres (verified against USGS EPQS). A converter would launder
  that mislabel, so `read_3dep_va_shapefile` has **no default** for `source_ele_units`.

A checkpoint that cannot be used raises. It is never dropped from the set, so it still
appears in any listing — an unusable control point is a reported result.

---

## 4. Part 2 — the tie estimator (`tie.py`), and the radius problem

### The pathology

Fit a plane to the ground returns within *R* of the mark, read it at the mark, and the
answer walks by a metre with *R*. Measured on gen1 at checkpoint 2210 (CSF ground, no
datum terms), lidar − surveyed in mm:

| R (m) | 2.5 | 5.0 | 7.5 | 10.0 | 15.0 | 20.0 | 25.0 |
|---|---|---|---|---|---|---|---|
| **order 1 (plane)** | −130 | −190 | −261 | −385 | −994 | −1297 | −1379 |
| its fit RMS | 45 | 133 | 329 | 553 | 777 | 827 | 825 |
| **order 2 (the fix)** | **−75** | **−85** | **−86** | **−81** | −284 | −684 | −964 |
| its fit RMS | 19 | 42 | 102 | 196 | 392 | 492 | 527 |

### The diagnosis: curvature, not noise

These marks are sited for *survey* convenience — they sit on local topographic highs. At
2210 the surveyed 349.288 m sits at the **p95** of gen1 returns within 5 m, with the
ground falling away inside 10 m: a road crown or shoulder. A least-squares plane has no
curvature term, so over a patch of curvature *k* its value at the centre is displaced by
≈ *kR²*/2 — quadratic in *R*, exactly the shape of the first row. (The synthetic dome in
`test_plane_fit_is_radius_dependent_on_a_local_high` checks that law analytically: over a
disc the plane's constant term is the *mean* of *z*, and mean(*r*²) = *R*²/2.)

### The fix, which this codebase had already made

`pipeline._poly2_ground` exists for the same reason and says so:

> Curvature-UNBIASED, unlike the per-cell median (which carries the cell's curvature) or
> a plane (which has no curvature term).

Reading the constant term of a **2nd-order** local surface removes the radius walk: 11 mm
across a 4× range of window (R = 2.5 → 10 m), against 255 mm for the plane. Nothing is
tuned per checkpoint; order 2 is the order the pipeline already reads its ground at, and
`surface_order` is an argument, so `surface_order=1` reproduces the pathology on demand.

### What the estimator actually is

The project's slope-normal ground read, generalised from a cell centre to an arbitrary
point. With `ground="slope_normal"` the pipeline's ground of a cell is

```
Zreg(cell centre) + quantile_q( z_i − [Zreg + dE·∂Zreg/∂E + dN·∂Zreg/∂N] )
```

— *a smooth local reference surface read at the target point, plus a quantile of the
vertical residuals to it*. Here the smooth surface is fitted locally at `surface_order`
and the target point is the checkpoint:

```
z_lidar(mark) = S(mark) + quantile_q( z_i − S(x_i, y_i) )
```

with `q = 0.50` — `ground_percentile` from `data/derived/elbaext/corrections_geoid.json`.
The ground returns themselves come from the pipeline's own ground source
(`ground_source = csf`), run on a crop rather than a tile; on this data CSF and the vendor
class-2 ground give ties within 1–8 mm of one another, so the choice does not carry the
answer.

### Radius is an output, not a setting

`estimate_tie` never returns a bare number. It returns the whole `RadiusEstimate` ladder
— radii at `res/2, res, 1.5·res, 2·res, 3·res, 4·res, 5·res`, all multiples of the
pipeline's 5 m grid — each with its return count, fit RMS, median residual, local slope,
relief, contributing flight lines and scan-angle spread. The headline is read at
`1.5·res = 7.5 m`, the half-width of `_poly2_ground`'s 3×3 window; `tie_median_mm` gives
the median over the pipeline-scale radii as a robustness companion. **The uncertainty is
half the spread across those radii**, not the standard error of the fit — the fit SE is
optimistic by an order of magnitude here because returns a metre apart are not
independent.

### Usable or not is the caller's call

`TieEstimate.verdict(tolerance_mm, tolerance_source=...)` has **no default tolerance** and
prints the source it was handed. Whether a radius spread of 105 mm disqualifies a mark is
a decision about the science, not a property of the code. The example script passes
161 mm — gen1's own published vertical RMSEz for Winona County (InPort 68818) — on the
grounds that a tie whose radius sensitivity exceeds the dataset's own accuracy cannot
inform anything.

---

## 4a. Part 2b – combining ties (`datum.py`)

Two ties that agree to 7.5 mm, both with sigmas larger than that, are a **consistency
check**. Averaging them and quoting `sigma/sqrt(2)` would turn the check into a claim,
because the two ties share one extrapolated lateral shift, one alignment estimator, one
ground source and one un-applied drift term.

`combine_ties` therefore takes `BudgetTerm`s that each carry a `kind`: `"random"` averages
down with the marks, `"common"` enters the total at full size no matter how many marks are
added, and `"unmodelled"` is a bound reported beside the total and never folded into it. A
term with no `source` is refused, so an unattributed error number cannot enter a budget.
The uncertainty is a **table**, not a single plus-or-minus.

## 5. Part 3 — the swath chain (`chain.py`)

### Search order is minimum work, and it is code, not a comment

1. **Along-swath, zero links.** `plan_path` first asks whether any line that actually
   puts returns on the mark *is* a target line. If so the tie applies directly, no
   cross-swath transfer happens, and nothing else is searched. At Elba the answer is no —
   but the test runs every time, because at another site it will be yes and the chain will
   be free.
2. **Otherwise the shortest chain**, breadth-first, minimising **link count** — not
   distance, not tile count. Each link contributes its own alignment error and a chain has
   no redundancy to absorb it.

`covering_lines` answers step 1 from the returns themselves (are there terrain returns
within *R* of the mark?), not from a fitted nadir track, so it cannot be fooled by a line
that ends before it gets there.

### Each link costs its overlap, not its tiles

`solve_link` loads only the two lines' points inside their shared bounding box, plus two
cells of margin for the Nuth & Kääb slope/aspect gradients, and hands them to the repo's
own `coreg.coregister_swaths`. Only the tiles the *chosen path* names are ever read: the
path is planned from the line geometry first, then data is pulled for exactly those links.
That ordering is the difference between four tiles and thirty-three.

Overlap areas are measured on a **global** cell grid anchored at (0, 0), so a pair whose
overlap straddles two tiles is counted once rather than twice.

Everything runs on the **vendor** classification (`~isin(class, (5,6,9))`, exactly
`coregister_swaths`' selection), so the chain needs no CSF pass at all.

### Error along a chain

A chain's misclosure is identically zero, so its internal residuals cannot see an
accumulated error, and the quadrature sum of per-link σ (sub-millimetre here) is not to be
believed. The only real check is a **second, independent path**: `compare_paths` reports
the disagreement between routes and puts the formal σ beside it for contrast.

---

## 6. What the example run produced

`analysis/groundtruth/elba_absolute_tie.py`, seven gen1 tiles already on disk, nothing
downloaded at run time. Every checkpoint reported separately, before any combination:

| checkpoint | type | line | links | n | tie (mm) | σ (mm) | median over radii | chain (mm) | geoid (mm) |
|---|---|---|---|---|---|---|---|---|---|
| 2210 | NVA | 128 | 5 (west) | 124 | **+21.3** | 12.4 | +29.1 | −36.9 | +74.0 |
| 3056 | VVA | 128 | 5 (west) | 121 | −103.2 | 52.3 | −138.1 | −36.9 | +73.8 |
| 2024 | NVA | 129 | 4 (west) | 200 | +156.6 | 54.5 | +139.3 | +53.3 | +70.2 |
| 2036 | NVA | 144 | 6 (east) | 115 | **+28.9** | 27.0 | −4.6 | −12.1 | +67.8 |
| 2099 | NVA | — | — | — | *not attempted* | | | | |
| 3089 | VVA | — | — | — | *not attempted* | | | | |

2099 and 3089 need tiles `4342-26-61` and `4358-26-02`, which are not on disk. They are
reported, not dropped.

Three independent disagreements, each of which means something different:

| pair | what a disagreement means | measured |
|---|---|---|
| 2210 vs 3056 | two marks on **one line**: the ESTIMATOR | +124.5 mm |
| 2210 vs 2024 | two links of one chain: the LINK 128–129 | −135.2 mm |
| 2210 vs 2036 | **west chain vs east chain**: the CHAIN | **−7.5 mm** |

Two strata, both fixed by the data and the method before any tie was computed:

| stratum | n | spread | mean |
|---|---|---|---|
| all control points | 4 | 259.8 mm | +25.9 mm |
| on the corridor band (\|ΔN\| ≤ 1 km) | 3 | 132.0 mm | −17.7 mm |
| on band **and** NVA (open ground) | 2 | **7.5 mm** | **+25.1 mm** |

**The result.** Two fully independent chains — five links west through lines
128-129-130-131-132-133, six links east through 144-143-142-141-140-139-138 — land within
**7.5 mm** of each other: gen1 at Elba reads **≈25 mm low** against surveyed
NAVD88(GEOID18), i.e. **+25 mm** must be added. Each mark's own radius uncertainty
(±12 and ±27 mm) is larger than their disagreement, so 7.5 mm should be read as "the two
chains do not contradict each other", not as the accuracy of the tie.

**What does not agree, and why that is worth knowing:**

* **3056 (VVA) vs 2210, 125 mm apart on the same line, 120 m apart on the ground.** Its
  radius curve never plateaus (−202 → −173 → −103 → −97 → +14 → +77 → +102) and its median
  residual peaks at +50 mm at R = 10 m — a one-sided lift consistent with vegetation in the
  ground class. 3DEP's own published VVA spread for this block is 27 cm at the 95th
  percentile against 3.5 cm RMSE for NVA. **This is an estimator/terrain limit at a
  vegetated mark, not a chain error**, and it is exactly what the two-marks-on-one-line
  comparison exists to separate.
* **2024, 135 mm from 2210.** It sits **3.24 km north** of the corridor band the links were
  solved in. The pipeline's per-swath along-track drift term, which this module cannot
  apply (it is fitted against the gen2 grid, which does not exist at the checkpoints),
  measures **11–29 mm/km** on the elbaext swaths — so ~50 mm of the 135 mm is unmodelled
  drift. **The remaining ~85 mm is unexplained.** Its radius curve also rises monotonically
  (+61 → +122 → +157 → +170), so the mark itself is not on a stable patch.

---

## 7. What this does not do

* **No along-track drift.** `pipeline.fit_along_track_drift` regresses against the gen2
  grid, and gen2 was deliberately not fetched at the checkpoints. Every tie is therefore
  valid at the along-track position of the mark, and carries ~16 mm/km of unmodelled drift
  in transferring that to Elba. This is the single largest known gap.
* **gen2 at the marks: now done, and it moved the budget.**
  `analysis/groundtruth/gen2_checkpoint_tie.py` reads gen2 against all six marks with no
  chain, no geoid and no lateral term, from six full-density 400 m boxes. gen2 is not
  systematically low (median +17.2 mm over four NVA marks), but those four span **103 mm**
  and give **40.8 mm RMS** on the cleanest case the method has. That scatter is a property
  of the marks, it applies equally to gen1, and it is the largest single term in the datum
  budget. See `analysis/ABSOLUTE_BASIS_ELBA.md` §3.
* **The datum constant itself, and its budget:** `analysis/ABSOLUTE_BASIS_ELBA.md`.
  The two ties below combine to **+22.7 mm** to be added to gen1, with sigma_total
  **39.7 mm** and a 42.6 mm unmodelled bound held outside it, via
  :mod:`lidar_diff_icp.groundtruth.datum`, which keeps common-mode terms at full size.
* **The lateral term is extrapolated.** The elbaext Nuth & Kääb shift (−0.750, −0.189 m)
  is applied to put gen1 in the checkpoints' horizontal frame, but it was measured against
  gen2 at Elba, 7–16 km away. Its effect on each tie is reported (−7 to +43 mm), so the
  sensitivity is visible; its validity out at line 128 is an assumption.
* **Two checkpoints untried.** 2099 and 3089 need two more tiles.
* **Two independent paths, not three.** A third would turn a "they agree" into a
  distribution. 2099 (line 130, same western chain, 9 km north) would add one for two more
  tiles, but it would inherit the same along-track problem as 2024.

---

## 8. Reproducing

```bash
env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
    analysis/groundtruth/elba_absolute_tie.py            # west + east, csf ground
env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
    analysis/groundtruth/elba_absolute_tie.py --no-east --ground vendor
```

Needs, all local: `data/before/4342-29-6{1,2,3,4}.laz` (west corridor),
`data/before/4358-29-0{1,2,3}.laz` (east corridor), `data/before/4342-28-61.laz`
(checkpoint 2024), `data/derived/elbaext/corrections_geoid.json`, and PDAL with
`filters.csf` (the conda `lidar-icp` env). The first run builds a terrain-point cache
under `data/derived/groundtruth/`; later runs reuse it. `references.geoid_difference`
fetches the PROJ geoid grids on first use unless `PROJ_NETWORK=OFF`.
