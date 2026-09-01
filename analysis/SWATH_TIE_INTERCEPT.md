# Making the swath-to-swath tie extent-invariant, and what it is worth

**Date:** 2026-08-26
**Code:** `coreg.across_track_tie` + `tie=` on `coregister_swaths`/`align_swaths`
(commit `89b6234`), regression test (`203c9c4`), measurement script
`analysis/swath_tie_intercept.py` (`4456295`).
**Runs:**

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/swath_tie_intercept.py
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/stable_point_tilt_audit.py --dod dod_cover_q2.npy
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/stable_point_tilt_audit.py --dod dod_cover_q2_tie.npy

Every input, parameter, mask and column is declared through `trust/provenance.py`.
**Nothing in `coreg.py` or `pipeline.py` changes behaviour**: the new tie is an opt-in
argument that defaults to the shipped estimator, and elba's re-derived constants still
reproduce `corrections.json`.

## Vocabulary note, 2026-09-01

This report predates the rename and says "gauge" throughout. Two different things wear
that word, and they are now named separately:

* **zero line** -- the flight line defined as zero when ONE tile's swath network is
  solved. Per-tile, arbitrary, sets only the level the tile inherits, cancelled exactly by
  an absolute datum. Recorded as `zero_line` in each `corrections.json`
  (`scripts/backfill_zero_line.py`). Where this report says "the shipped gauge" or
  "gauged on", read *zero line*.
* **common line** -- a flight line present in BOTH tiles, used to re-express each against
  one reference so their constants can be compared. Where this report says "re-referenced
  to 135" or "gauge 137", read *common line*.

Neither changes any swath-to-swath difference, and neither sets goodness of fit: every
line is used in the network solve regardless.

## What was implemented

`analysis/SWATH_ACROSS_TRACK_TEST.md` established that the between-line height difference
is not flat across the sidelap, that the coefficient is per flight-line pair rather than
sensor-wide, and that the tie `align_swaths` takes is therefore `k + c·mean(dtan)` — an
average over whatever part of the sidelap the tile happens to cover. Two changes follow
from that, and both are implemented here. A shared roll term is **not** added; that report
rejects one at p = 6×10⁻⁶⁴ and so do the repo's own `estimate_boresight` numbers.

1. **The tie is taken at across-track position zero.** `coreg.across_track_tie` fits
   `dh = k + c·dtan` over the pair's overlap cells and returns `k`, the value at
   `tan θ_ref = tan θ_src` — the middle of the sidelap, a position defined by the flight
   geometry rather than by the tile boundary. The fit is **LAD (median regression)**, not
   OLS, so that at `c = 0` it reduces *exactly* to the median the shipped tie already
   takes. That nesting is what makes the old and new numbers comparable rather than two
   different estimators.
2. **The gauge.** `align_swaths` gained no new gauge — `ref=` and `ref=None` already
   exist — but the choice is now measured rather than inherited, and its docstring says
   what it costs. See §4.

### Extent-invariance by construction, and shown to bite

`tests/test_coreg.py::test_swath_tie_intercept_is_extent_invariant` builds two N–S lines
flown there-and-back (so `tan θ_A + tan θ_B` is the fixed spacing/height ratio, the
geometry measured at Elba), each carrying its own across-track error, and samples them
over a narrow easting window and a wide one. The shipped tie moves between the two extents
by `c·(mean dtan₁ − mean dtan₂)` — asserted against that closed form, not merely against
"they differ". The intercept tie returns the same number on both, and the injected
constant. **Shown to fail without the fix:** with `across_track_tie`'s result discarded so
the intercept path falls back to the median, the two extents disagree by **21.4 mm** and
the last two assertions fail. The full test suite passes with the fix in place.

The estimand argument is the same one: `median(dh)` estimates `k + c·E[dtan | extent]`,
which is a function of the extent; the intercept estimates `k`, which is not.

### The shipped ties are reproduced before anything is compared to them

| tile | worst \|re-derived − shipped\| dx | dy | **dz** | file precision |
|---|---|---|---|---|
| elba_fulldensity | 0.02 mm | 0.04 mm | **0.05 mm** | 0.05 mm |
| elbaext | 0.22 mm | 5.35 mm | **0.10 mm** | 0.05 mm |

`align_swaths(pc, ref=int(ps.min()))` on each tile's cached CSF cloud returns
`corrections.json`'s own constants. The gate is on `dz`, the axis this work changes, at a
**1 mm** tolerance I chose and declared (`R.param(..., src="MINE")`) — a sixth of the
smallest tie change measured below. elbaext's 5.35 mm horizontal gap on swath 138 is
reported, not hidden: its `corrections.json` was written 2026-08-22 and the cached cloud
it is being re-derived from was rewritten 2026-08-25, so they are not the same run. Its
*vertical* solution still lands within 0.10 mm.

## 1. The new ties against the old

`c_mm` is the across-track slope fitted alongside the intercept; `mean_dtan` is where in
the sidelap this tile's coverage sits; `pred = −c·mean_dtan` is what the mechanism says
`delta` should be.

| tile | pair | cells | dz_med | dz_int | **delta** | c (mm/unit tan) | mean_dtan | pred |
|---|---|---|---|---|---|---|---|---|
| elba | 135-136 | 221,308 | −23.90 | −18.03 | **+5.87** | +106.9 | −0.0607 | +6.49 |
| elba | 136-137 | 309,425 | −8.56 | −8.67 | −0.11 | +157.7 | +0.0093 | −1.47 |
| elba | 137-138 | 307,230 | −11.29 | −11.94 | −0.65 | +45.3 | +0.0065 | −0.29 |
| elbaext | 133-134 | 379,429 | +18.88 | +21.99 | **+3.11** | +187.7 | −0.0164 | +3.08 |
| elbaext | 134-135 | 390,485 | −16.17 | −15.80 | +0.37 | +51.8 | +0.0025 | −0.13 |
| elbaext | 135-136 | 414,285 | −15.88 | −15.98 | −0.10 | +117.0 | −0.0051 | +0.60 |
| elbaext | 136-137 | 385,725 | −6.92 | −8.60 | −1.68 | +150.2 | +0.0110 | −1.65 |
| elbaext | 137-138 | 401,072 | −3.54 | −4.24 | −0.70 | +41.7 | +0.0112 | −0.47 |

Three things worth stating.

**The change is concentrated exactly where the mechanism says it should be.** The only two
links whose tie moves by more than 2 mm are the two edge-cut links — elba's 135-136 and
elbaext's 133-134 — and they are the two whose `mean_dtan` is largest in magnitude
(−0.0607, −0.0164 against ≤0.0112 elsewhere). Every interior link moves by under 1.7 mm.
That is the extent-dependence, isolated.

**`c` measured here bears on the across-track report's one anomaly.** These coefficients
come from the co-registration's own 2 m vendor-class grids (221k–414k cells) rather than
the 5 m CSF reference cells (7k–15k), so exact agreement is not expected and is not the
claim. Six of the eight pairs have the same sign and the same order as that report's
Table 1, differing by 0.2–3.7 of its own standard errors. The two exceptions are the
informative ones: for both 137-138 pairs the report's raw fit gave **−43.8** and **−34.0**,
and only after adding a per-pair along-track control did they become **+34.5** and
**+36.7**. On the full 2 m overlap they come out **+45.3** and **+41.7** — the
along-track-controlled values, not the raw ones. The report's reading of those two pairs
as along-track-contaminated is supported here on an independent population.

**`pred` is not a perfect predictor of `delta`, and that is expected.** It matches on the
two links that matter (+6.49 against +5.87; +3.08 against +3.11) but over-predicts the
small interior ones. `pred` uses the unweighted mean `dtan`; `dz_med` is a *median* over
cells that are not uniform in `dtan`. Where the effect is ~0.5 mm that difference shows.

## 2. The elba / elbaext disagreement: 8.0, 9.8, 17.4 mm

Re-referenced to swath 135, as in the across-track report:

| swath | elba old | ext old | **dis old** | elba new | ext new | **dis new** | removed |
|---|---|---|---|---|---|---|---|
| 136 | −23.90 | −15.88 | **−8.02** | −18.03 | −15.98 | **−2.04** | −5.97 (75%) |
| 137 | −32.46 | −22.80 | **−9.66** | −26.70 | −24.59 | **−2.11** | −7.55 (78%) |
| 138 | −43.75 | −26.34 | **−17.41** | −38.63 | −28.82 | **−9.81** | −7.60 (44%) |

(The earlier report's 8.0 / 9.8 / 17.4 came from `corrections.json`, which stores the
shifts rounded to 0.1 mm; the unrounded network solve gives 8.02 / 9.66 / 17.41.)

Re-referencing is itself a gauge, so those numbers are determined only up to a common
constant. The gauge-free statements:

| shared lines | spread old | spread new | RMS old (mean removed) | RMS new |
|---|---|---|---|---|
| 135–138 (all) | 17.41 | **9.81** | 6.18 | **3.75** |
| 135–137 (138 dropped) | 9.66 | **2.11** | 4.22 | **0.98** |

Both rows are shown; the second is not a filter but the point — **everything left after
the fix is swath 138.** Excluding it, the two tiles agree about the remaining three
shared flight lines to **0.98 mm RMS**, from 4.22 mm.

**Swath 138's −8.9 mm survives.** The across-track report predicted −8.5 mm of the −17.4
from extent-dependent sampling and left −8.9 mm outstanding. Measured: the fix removes
−7.60 mm and leaves **−9.81 mm**. The outstanding residual is real, is essentially the size
the report predicted, and is not a tie-averaging artefact. It is now the *only* thing
separating the two tiles' swath solutions.

## 3. The tilt

**Read this against `analysis/STABLE_POINT_TILT_AUDIT.md`, which landed while this work
was running and changed the answer.** The `dE = −14.19 ± 5.15, dN = −16.70 ± 3.65 mm/km,
n = 24,287` figures the task names were measured on the pre-`7701383` reference-cell
population. Commit `7701383` then excluded the valley floor from `reference_cells` by
default, because 19% of that population sat on valley terraces and carried the whole
easting term. The comparison below is therefore run **twice**, and in each case the before
and the after are produced by the same code on the same population — the only thing that
differs is the raster.

The tilt is re-derived by `analysis/stable_point_tilt_audit.py` (committed at `199cef3`,
before this work), not inline. It reproduces `−14.19 / −16.70 / n = 24,287 / intercept
−10.30` exactly on the old population. (One correction to the note that recorded it: those
SEs of 5.15 and 3.65 are not what the block bootstrap returns — 500 replicates at 50 m
give **1.81 / 3.37 / 2.21**, and the cluster-robust sandwich gives 3.24 / 2.10. The point
estimates are exact; the quoted uncertainties are not reproduced.)

### On the current (valley-excluded) population, n = 19,683, 50 m blocks

| raster | mean (mm) | **dE (mm/km)** | **dN (mm/km)** | NMAD (mm) |
|---|---|---|---|---|
| `dod_cover_q2.npy` (shipped) | −4.74 | **+3.51 ± 2.61** | **−17.90 ± 1.76** | **47.4** |
| `dod_cover_q2_tie.npy` | −10.01 | **+2.51 ± 2.61** | **−17.84 ± 1.76** | **47.5** |
| `..._tie_gauge.npy` | −36.71 | +2.51 ± 2.61 | −17.84 ± 1.76 | 47.5 |

### On the pre-`7701383` population, n = 24,287, 50 m blocks

| raster | mean (mm) | dE (mm/km) | dN (mm/km) | NMAD (mm) |
|---|---|---|---|---|
| `dod_cover_q2.npy` (shipped) | −10.30 | −14.19 ± 3.24 | −16.70 ± 2.10 | 51.5 |
| `dod_cover_q2_tie.npy` | −15.61 | −14.77 ± 3.22 | −16.64 ± 2.10 | 51.4 |

**The extent-invariant tie does not move the tilt.** `dE` moves by −1.00 mm/km on the
current population and by +0.58 mm/km (i.e. the *wrong way*) on the old one; both are
0.38 σ and 0.18 σ of their own cluster-robust errors. `dN` moves by 0.06 mm/km — nothing. Under L1 the same: dE +4.85 → +4.00,
dN −18.97 → −18.89.

**This is a negative result and it is the expected one.** The hypothesis was that `dE` is
residual across-swath structure, since the per-swath constants run ≈ −22 mm/km in easting.
But the constants become *flatter* under the new tie, not steeper: their west-to-east
spread goes from 43.75 mm to 38.63 mm across the ~2 km of flight-line spacing, so if
anything they remove *less* easting gradient than before. That is why `dE` gets marginally
more negative on the old population. The tilt audit had already reached the same
conclusion by a different route — `dE` lives in the valley limb (−84.8 mm/km there against
+3.5 on the upland), is absent inside the two flight lines that have easting leverage, and
is not significant at any block size ≥ 250 m — and this measurement is consistent with it.
**The per-swath constants are not the easting tilt, before or after the fix.**

**Stable-ground scatter: flat.** NMAD 47.4 → 47.5 mm (current population), 51.5 → 51.4 mm
(old). The tie is a per-swath constant; it moves levels, not scatter, and the measurement
says so. It neither improves nor damages the DoD's precision.

### How the DoD was rebuilt, and how that path was checked

The DoD was not recomputed by re-running `difference_dem`. The per-cell gen1 ground was
re-reduced from `beam_offset_table.parquet` with the pipeline's own estimator — the
`ground_q = 0.50` quantile of the slope-normal residual per cell (`pipeline.py:603`) — with
the new per-swath constant added to each return, so a changed constant propagates through
the median of a mixed-line cell exactly as the pipeline propagates it.

That path is **validated, not assumed.** Subtracting the whole shipped per-swath term back
out and re-fitting reproduces `STABLE_POINT_TILT_AUDIT.md` §5c, which reached the same
quantity by a different route (per-(cell,line) medians differenced against the per-cell
median), to the printed precision **on both populations**:

| population | route | mean | dE | dN |
|---|---|---|---|---|
| pre-`7701383` | audit §5c | −36.16 | −24.25 | −15.11 |
| pre-`7701383` | this run's `_noswath` raster | **−36.16** | **−24.25** | **−15.11** |
| current | audit §5c | −29.22 | −8.41 | −14.92 |
| current | this run's `_noswath` raster | **−29.22** | **−8.41** | **−14.92** |

**One limitation, stated because it would otherwise be a trap.** Sections 5, 5a and 5c of
`stable_point_tilt_audit.py` decompose the raster by flight line using `d_mm_corr` from the
parquet, which still carries the *old* per-swath term. On a rebuilt raster those sections
mix a new per-cell DoD with old per-line offsets and are **not** interpretable. Only §2
(the plane, the scatter, the L1 fit, the block-median fit) is used above.

## 4. The gauge

Measured per swath over the in-grid gen1 CSF ground returns of each tile:

| tile | swath | returns | mean scan | sd | p5 | p95 | \|mean\| |
|---|---|---|---|---|---|---|---|
| elba | **135** | 498,847 | **−11.79** | 1.92 | −15.00 | −9.00 | **11.79** |
| elba | 136 | 2,424,563 | +1.77 | 7.68 | −10.00 | +14.00 | 1.77 |
| elba | 137 | 2,544,119 | +0.58 | 8.52 | −13.00 | +14.00 | 0.58 |
| elba | 138 | 1,303,950 | −7.17 | 4.86 | −15.00 | +0.00 | 7.17 |
| elbaext | **133** | 558,218 | **−12.14** | 2.05 | −15.00 | −9.00 | **12.14** |
| elbaext | 134 | 2,692,161 | +1.70 | 7.60 | −10.00 | +14.00 | 1.70 |
| elbaext | 135 | 3,042,616 | +0.27 | 8.62 | −13.00 | +14.00 | 0.27 |
| elbaext | 136 | 3,124,488 | +0.13 | 8.68 | −14.00 | +14.00 | 0.13 |
| elbaext | 137 | 2,966,122 | +0.55 | 8.52 | −13.00 | +14.00 | 0.55 |
| elbaext | 138 | 1,544,958 | −7.29 | 4.85 | −15.00 | +0.00 | 7.29 |

Both tiles pin their network to a line they see only from one side, over a 6° one-sided
window, while the interior lines are sampled across roughly ±13–14°. Of the lines
**both** tiles carry, the one whose worst |mean scan| across the two tiles is smallest is
**swath 137** (0.58 / 0.55°). That rule — stated, computed, flagged `MINE` in the
provenance banner — is how the candidate was picked.

**But the gauge is not a free improvement, and I am not recommending it be flipped as part
of this change.** `tests/test_coreg.py::test_align_swaths_reference_is_a_gauge_not_an_observation`
already proves the gauge changes no swath-to-swath difference. It changes exactly one
thing: **the level the whole tile's DoD sits at**, because nothing in the pipeline ties
gen1 to gen2 vertically — the geoid term is a constant computed from PROJ grids, not a
fit. Point-count-weighted mean shift of the gen1 cloud:

| tile | tie | gauge 135/133 (shipped) | gauge 137 | zero-mean (`ref=None`) |
|---|---|---|---|---|
| elba | overlap_median | −29.18 | +3.28 | −4.15 |
| elba | **intercept** | **−23.92** | **+2.77** | **−3.08** |
| elbaext | overlap_median | −5.62 | +14.48 | +0.27 |
| elbaext | **intercept** | **−3.03** | **+15.37** | **+0.75** |

Moving elba from gauge 135 to gauge 137 raises the gen1 cloud by 32.5 mm, which **lowers
the DoD by 32.5 mm everywhere**: the stable open divide mean goes from −4.74 to −36.71 mm.
That is a large, consequential move, and the evidence points both ways:

* **Against.** Stable upland divides would then read 37 mm of erosion over 13 years. The
  divides are the most precise constraint the project has on the epoch difference.
* **For, weakly.** `analysis/ABSOLUTE_BASIS_ELBA.md` puts gen1 **+22.7 ± 39.7 mm low**
  against two surveyed NAVD88(GEOID18) marks, i.e. it asks for gen1 to be *raised* — the
  same direction. On gauge 137 those marks would read gen1 as ~10 mm high instead of
  ~23 mm low, a smaller residual. But σ = 39.7 mm cannot choose between the two, and the
  tie itself is measured against the gauged cloud, so it would have to be re-derived.
* **Zero-mean is not the answer for cross-tile work.** It looks attractive above (the two
  tiles' levels land 3.8 mm apart) but that is coincidence: the two tiles carry different
  swath sets, so "zero-mean" means a different thing in each, and it does **not** put a
  shared flight line at a shared level. Only a *named line present in both tiles* does
  that. On gauge 137 with the intercept tie, elba and elbaext apply
  constants to swath 136 that differ by **0.07 mm** (from 1.64 mm), to swath 135 by
  **2.11 mm** (from 9.66 mm), to 137 by zero by construction — and to swath 138 by
  **7.70 mm**, essentially unchanged from 7.75 mm. The residual is 138 and only 138.

So the gauge is a level decision that belongs with the absolute-datum work, not a
by-product of a tie fix. `dod_cover_q2_tie_gauge.npy` is on disk so the option can be
inspected; it is not proposed for adoption here.

## 5. The integration I propose

**Adopt now — the tie.** In `pipeline.difference_dem`, change line 670 from

    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))

to

    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()), tie="intercept")

and record it in `corrections.json` as a new key (`swath_tie: "intercept"`), so an older
product is identifiable rather than silently different. Nothing else changes: the
horizontal solution, the network, the weights, the gauge and every other term are
untouched. Blast radius, measured rather than estimated:

* Every tile's `corrections.json` and every `dod*.npy`, `lod*.npy`, `z_before*.npy` derived
  from it must be rebuilt, or explicitly kept and labelled as `swath_tie: "overlap_median"`.
  At elba the DoD moves by −5.76 mm (median; range −5.87 to 0.00). Existing products remain
  reproducible: `tie="overlap_median"` is the default and reproduces them exactly.
* `beam_offset_table.parquet`'s `dz_swath_mm` column is read from `corrections.json`, so
  the two tiles' beam tables must be regenerated in step with it or they will mix ties.
  Every analysis that uses `d_mm_corr` inherits that.
* `analysis/ABSOLUTE_BASIS_ELBA.md`'s +22.7 mm constant is transported through swath
  chains and would move; at elba the level moves by +5.26 mm (−29.18 → −23.92).
* The gain: the tile-to-tile disagreement about the same flight lines falls from 6.18 to
  3.75 mm RMS over all four shared lines, and from 4.22 to 0.98 mm once swath 138 is set
  aside. Nothing measurable is lost: the tilt does not move (≤0.4 σ) and the stable-ground
  NMAD does not move (0.1 mm).

**Do not bundle — the gauge.** Keep `ref=int(ps8.min())` for now. Raise the gauge as its
own decision, argued against the absolute datum, and if it is taken, take it as a *named
line present in every tile of the project* rather than `min()` or zero-mean. The
`align_swaths` docstring now states what the gauge costs so the next reader does not have
to rediscover it.

**Still open.** Swath 138. Its −9.81 mm is now the whole of the remaining tile-to-tile
disagreement, it was predicted to be left over, and it is not an extent artefact. Two
things about 138 are on record and worth putting together: it is the other one-sidedly
sampled line (mean scan −7.2°, p95 = 0.00 — the tile never sees its far side), and both
tiles' 137-138 link is the one with the along-track contamination the across-track report
isolated. That is where to look next.
