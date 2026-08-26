# The cross line: what it settles, and the thing it settles against the model

**Date:** 2026-08-26
**Code:** `analysis/crossline_fit.py` (the fit), `analysis/slope_bias/csf_tiled.py` (the CSF pass)
**Tests:** `tests/test_crossline_fit.py`, `tests/test_csf_tiled.py` — 12 pass
**Runs:**

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/csf_tiled.py \
        --src data/before/4342-28-64.laz --out data/csf_cache/4342-28-64.las --nx 2 --ny 2
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/crossline_fit.py

**Read, not redone:** `analysis/SWATH_DEGENERACY_BREAKING.md` (which found the line and ranked
this first), `analysis/SWATH_ACROSS_TRACK_TEST.md` (the model and the estimator),
`analysis/SWATH_TIE_INTERCEPT.md` (the sibling change to `coreg.py`),
`analysis/FRAME_2026-08-26-PM.md`. Every table below is emitted by the script that computed it
through `trust/provenance.py`; the banner carries each input's size, mtime and digest, each
parameter's source, and the two parameters flagged `MINE`. `coreg.py` and `pipeline.py` are
imported, not modified.

---

## The headline, in three sentences

**The CSF pass worked and the degeneracy is broken.** Line 10010 now carries **640,397** ground
returns at **0.45 pts/m²**, inside the 0.46–0.52 of the N-S lines it crosses, with the scan
angle intact through the round trip; on its three pairs `SE(q)/SE(p)` is **1.02–1.12** against
4.95–7.11 on the N-S pairs of the same tile, so each line's own coefficient is estimable.

**The coefficients come out, and then fail their first out-of-sample test.** Solved from the
cross pairs alone, `c_136 = +162.0 ± 26.0`, `c_137 = −7.1 ± 8.4`, `c_138 = −36.1 ± 27.8`,
`c_10010 = −73.6 ± 17.8` mm per unit tangent. Those predict `(c_136+c_137)/2 = +77.4` and
`(c_137+c_138)/2 = −21.6`, against **+135.8 ± 11.9** and **+52.8 ± 9.8** measured on the *same
tile's* own N-S overlaps: residuals **−58.3 (−3.22σ)** and **−74.3 (−4.27σ)**.

**Chasing that miss is the actual finding: `c_s` is not a property of the flight line.** Refit
the N-S pairs inside the cross line's own northing span — where the prediction was made — and
136-137 moves from +135.8 to **+61.6 ± 20.6**, which agrees with the +77.4 prediction to
**+0.64σ**. Without any band at all, the along-track gradient of the coefficient is
`+83.2 ± 27.1` mm per unit tangent per km on 135-136 (t = +3.1), and in six equipopulated
northing bands the same pair runs **+141, +115, +23, +31, +429, +274**. The across-track term is
real, but it is **local to the piece of overlap it is measured on**, not a line constant — so a
per-line `c_s` correction is not the object to build.

---

## 1. The CSF pass: it worked

`analysis/slope_bias/csf_tiled.py`, 2×2 tiles with the standard 150 m halo, calling
`ground.classify_ground_csf` at its repo defaults (rigidness 1, everything else PDAL's own).
Nothing retuned. **6,941,881 class-2 ground points** written to `data/csf_cache/4342-28-64.las`.

| psid | n_raw | vendor frac2 | vendor frac12 | n_csf | n_csf_last | csf_rate | pts/m² |
|---|---|---|---|---|---|---|---|
| 135 | 549,598 | 0.190 | 0.773 | 444,915 | 444,915 | 0.810 | 0.46 |
| 136 | 2,632,245 | 0.650 | 0.208 | 2,143,062 | 2,143,059 | 0.814 | 0.52 |
| 137 | 2,841,248 | 0.570 | 0.292 | 2,382,403 | 2,382,400 | 0.839 | 0.50 |
| 138 | 1,531,974 | 0.546 | 0.296 | 1,331,104 | 1,331,093 | 0.869 | 0.52 |
| **10010** | **747,107** | **0.000** | **0.998** | **640,397** | **640,397** | **0.857** | **0.45** |

The line with **zero** vendor bare-earth returns comes out of CSF with a ground rate of 0.857,
which is *inside* the 0.810–0.869 of the four vendor-classified lines, and a ground density of
0.45 pts/m² against their 0.46–0.52. **The blocker is removed, and the classifier does not
treat the cross line differently from the lines it crosses.**

**The scan angle survived**, which this project has silently lost twice: after the PF1 → PF7
promotion line 10010 reads **−17.00 to +2.00°, 91.0% nonzero**. `csf_tiled.py` asserts it before
shipping the cache, and `crossline_fit.scan_angle_deg` refuses an all-zero field rather than
returning it; `tests/test_crossline_fit.py` builds a zero-angle LAS and requires the refusal.

**One defect found and fixed on the way.** The tiled-CSF core rule widened the outer edges with
`core |= (gx >= cx1)`, a whole half-plane, so every tile in the last column claimed the
`x.max()` point. Measured on the real elbaext cloud with its shipped 3×3 tiling: the old rule
double-claims **8 of 17,354,919** points, the new rule claims every point exactly once. Real,
and small enough that `data/csf_cache/elbaext.las` does not need rebuilding. (Separately
measured so the two are not confused: **7.23%** of the raw `elbaext_gen1_merged.laz` returns and
**7.02%** of `4342-28-64.laz` are exact duplicates in `(x,y,z,gps_time,psid)` in the vendor
delivery itself. That is upstream of everything here.)

---

## 2. The one substitution, measured rather than argued

The eight N-S pairs are reduced with `d_mm_corr` — the slope-normal residual of a gen1 CSF
ground return to the **gen2** surface, after the pipeline's registration terms. **No gen2 covers
the cross-line tile**, which sits 2.6–3.5 km north of the elbaext grid, so `d_mm_corr` cannot be
formed there. This run therefore uses the identical estimator with the reference surface built
from the **pooled gen1 ground of the same tile**: per 5 m cell the median of all lines' ground
returns, gap-filled, smoothed by the pipeline's own `sn_smooth_cells = 1.2`, read as each line's
median residual to that common tilted plane. Flagged `MINE` in the banner. In a between-line
difference a common reference cancels — but that is an argument, so it is measured on elbaext's
N-S pairs, where both estimators can be run:

| pair | pipeline `d_mm_corr` | pooled-gen1 reference | difference | difference / √(SE²+SE²) |
|---|---|---|---|---|
| 133-134 | +215.1 ± 19.4 | +206.0 ± 32.8 | −9.1 | 0.24 |
| 134-135 | +56.8 ± 6.7 | +68.9 ± 9.3 | +12.1 | 1.06 |
| 135-136 | +143.1 ± 9.9 | +147.0 ± 20.6 | +3.9 | 0.17 |
| 136-137 | +134.0 ± 9.3 | +110.8 ± 17.5 | −23.2 | 1.17 |
| 137-138 | +48.7 ± 39.1 | +124.7 ± 20.6 | +76.0 | 1.72 |

**No pair's difference reaches twice that ratio**, and four of five stay under 1.2. Two
honesties about the ratio: the two estimators read the same CSF cloud through different
reductions, so they are *positively correlated* and a quadrature sum of their standard errors
overstates the error on their difference — the ratio is therefore optimistic, not conservative,
and it is a descriptive scale rather than a test. What the table actually licenses is a level:
**the substitution is validated at roughly ±25 mm per unit tangent**, no better. On 137-138 —
the pair both this document and the across-track report keep finding to be the awkward one — the
two estimators are 76 mm apart. Nothing below leans on the substitution being tighter than that,
and the two decisive results (§5, §6) compare estimates that share no data at all.

**A staleness warning that belongs here, because it changes numbers other documents quote.**
`SWATH_ACROSS_TRACK_TEST.md` §2 reports elbaext 135-136 / 136-137 / 137-138 as
+133.4 / +183.7 / **−34.0** and elba as +125.2 / +169.6 / **−43.8**. Re-run this session through
that script's own code, the same pairs are +143.1 / +134.0 / **+48.7** and
+141.4 / +146.7 / **+83.2**. The cause is verified, not guessed: commit `7701383` ("Exclude the
valley floor from divide reference cells by default") landed **after** that document's commit
`f76a767` (`git log 7701383..f76a767` is empty, i.e. the document is an ancestor), and it added
a floodplain mask and an elevation-antimode cut to `reference_cells` — 20,335 and 10,033 cells
at elba. **The 137-138 coefficient changes sign under that population change.** Any downstream
use of that table's numbers should be re-derived. (`data/derived/elba_fulldensity/dod.npy` was
also being rewritten by concurrent work during these runs; the banner records the digest each
run actually read, and the `|DoD| > 500 mm` cut it feeds removes 43 cells at elba, so the effect
is the ≤0.9 mm drift visible between successive runs of this script.)

---

## 3. The geometry, on real classified ground

`se_ratio = SE(q)/SE(p)` from the design alone; `q = (c_A − c_B)/2` is the combination the N-S
chain cannot see.

| pair | kind | cells | area | corr(tanA,tanB) | sd(sum) | sd(dif) | **SE(q)/SE(p)** |
|---|---|---|---|---|---|---|---|
| 135-136 | N-S/N-S | 38,104 | 95.3 ha | −0.954 | 0.0112 | 0.0722 | 6.45 |
| 136-137 | N-S/N-S | 61,482 | 153.7 ha | −0.921 | 0.0220 | 0.1088 | 4.95 |
| **136-10010** | **CROSS** | 20,766 | **51.9 ha** | **−0.039** | 0.1176 | 0.1224 | **1.04** |
| 137-138 | N-S/N-S | 67,894 | 169.7 ha | −0.961 | 0.0161 | 0.1145 | 7.11 |
| **137-10010** | **CROSS** | 42,339 | **105.8 ha** | **−0.134** | 0.1694 | 0.1896 | **1.12** |
| **138-10010** | **CROSS** | 22,868 | **57.2 ha** | **−0.024** | 0.1217 | 0.1246 | **1.02** |

`SWATH_DEGENERACY_BREAKING.md` §2 predicted 1.06 / 1.14 / 1.05 from the vendor classes and the
raw cloud; on our CSF ground the same numbers are **1.04 / 1.12 / 1.02**. **The geometric claim
reproduces.**

**And here is the fact that governs how the result must be read.** Within a flight line the
ground coordinate and the scan angle are the same variable. Measured, not assumed — the
correlation of each line's per-cell tangent with cell easting and northing:

| pair | r(tanA, E) | r(tanA, N) | r(tanB, E) | r(tanB, N) |
|---|---|---|---|---|
| 136-10010 | **+1.00** | −0.01 | +0.00 | **+1.00** |
| 137-10010 | **−1.00** | −0.00 | +0.12 | **+0.99** |
| 138-10010 | **+1.00** | +0.00 | −0.02 | **+0.98** |

On a cross pair, `c_A` **is** the easting gradient of the between-line difference field and
`c_B` **is** its northing gradient. That is exactly why the pair separates them — and it is also
why the separation cannot distinguish a per-line beam term from any other spatially structured
between-line difference. It is one fact with both consequences, and §6 is where it bites.

*(A note on what this rules out as a control: linear position controls would be collinear with
the two coefficients being estimated and would delete the signal rather than control it. The
sensitivity in §8 therefore uses quadratic terms `E², N², EN`, which are not.)*

---

## 4. Each line's own coefficient

Two estimators, and they must both be reported because they disagree.

**Cross pairs alone** (`cross_only_solve`: one intercept per pair, coefficients shared through
the cross line, 85,973 cell-pairs over 648 blocks, r² = 0.0197). Nothing an N-S overlap measures
enters, so this is the estimate whose predictions are testable:

| line | c (mm per unit tangent) | SE |
|---|---|---|
| 136 | **+162.0** | 26.0 |
| 137 | **−7.1** | 8.4 |
| 138 | **−36.1** | 27.8 |
| **10010** | **−73.6** | 17.8 |

**Joint within-cell solve of all five lines at once** (`g(cell,line) = μ_cell + a_line +
c_line·tan`, 418,916 (cell, line) rows on 194,793 multiply-covered cells, 2,276 blocks, within-
cell r² = 0.0472, design condition number 16.9). Cell fixed effects absorb terrain, cover and
any spatial field exactly as differencing within a cell does — it *is* differencing within a
cell, generalised from two lines to five — so this uses every constraint the tile offers,
including the N-S overlaps:

| line | c | SE |
|---|---|---|
| 135 | +120.4 | 58.0 |
| 136 | **+252.4** | 16.6 |
| 137 | −25.3 | 7.4 |
| 138 | **+82.3** | 18.3 |
| 10010 | −56.9 | 18.6 |

They differ by **118.4** mm per unit tangent on line 138 (−36.1 ± 27.8 against +82.3 ± 18.3)
and by **90.4** on line 136 — in both cases more than three times either standard error. (No
σ is quoted for the gap: the joint solve contains the cross-only data, so the two are correlated
and their difference has no clean error.) That gap is not a bug in either estimator: it is the
redundancy tension of §5, showing up as the price of adding the N-S constraints.

**The cross line's own coefficient is not self-consistent either.** Fitted pair by pair it is
**−193.6 ± 30.9** (from 136), **−31.2 ± 20.2** (from 137) and **−38.6 ± 30.4** (from 138):
χ² = 20.75 on 2 degrees of freedom, p = 3.1×10⁻⁵. Two of the three agree well; the estimate from
the 136 overlap is the outlier.

---

## 5. Redundancy — the first real test the model has had

The distinction that makes this a test rather than a tautology: the joint solve *uses* this
tile's N-S overlaps, so its same-tile rows are in-sample residuals of an overdetermined fit
(5 coefficients against 9 slope-informative combinations: 3 N-S pair sums plus 2 apiece from
the 3 cross pairs), **not** predictions. Only the cross-only
rows, and the elba/elbaext rows, are out of sample.

| source of prediction | pair | predicted | observed | resid | σ | what it is |
|---|---|---|---|---|---|---|
| **cross-only, same tile** | **136-137** | **+77.4** | **+135.8** | **−58.3** | **−3.22** | **OUT OF SAMPLE** |
| **cross-only, same tile** | **137-138** | **−21.6** | **+52.8** | **−74.3** | **−4.27** | **OUT OF SAMPLE** |
| joint, same tile | 135-136 | +186.4 | +171.3 | +15.0 | +0.38 | in-sample residual |
| joint, same tile | 136-137 | +113.5 | +135.8 | −22.2 | −1.47 | in-sample residual |
| joint, same tile | 137-138 | +28.5 | +52.8 | −24.3 | −1.82 | in-sample residual |
| joint → elba | 135-136 | +186.4 | +141.4 | +44.9 | +1.36 | out of sample, 2.6–3.5 km S |
| joint → elba | 136-137 | +113.5 | +146.7 | −33.1 | −2.49 | out of sample, 2.6–3.5 km S |
| joint → elba | 137-138 | +28.5 | +83.2 | −54.7 | −1.83 | out of sample, 2.6–3.5 km S |
| joint → elbaext | 135-136 | +186.4 | +143.1 | +43.2 | +1.46 | out of sample, 2.6–3.5 km S |
| joint → elbaext | 136-137 | +113.5 | +134.0 | −20.4 | −1.55 | out of sample, 2.6–3.5 km S |
| joint → elbaext | 137-138 | +28.5 | +48.7 | −20.2 | −0.50 | out of sample, 2.6–3.5 km S |

**Read the two bold rows first.** They are the same tile, the same cells, the same estimator,
the same session — every escape route the elba/elbaext comparison leaves open is closed — and
they miss by −3.22σ and −4.27σ. **On its first genuine test, the per-line model
`err_s(θ) = a_s + c_s·tan θ` does not reproduce what the N-S overlaps measure.**

The six southern rows are milder (−2.49σ worst, four of six under 1.6σ), which is worth
recording but should not be leaned on: they compare a northern coefficient set against southern
overlaps, so a mechanism that varies with position is being *averaged*, not tested.

---

## 6. Following the miss: `c_s` is not a property of the flight line

The cross line covers only the northern **820 m** of the tile's 3,500 m, so the cross-only
prediction is made *there* while the N-S pair sums above are measured over the whole tile. If
the coefficient varies along track, restricting the N-S fits to that same band should move them
toward the prediction. The band is not chosen: it is the exact northing span of the cross line's
own overlap cells, computed from them.

| pair | c, whole tile | **c, cross line's band** | SE | cross-only prediction | resid | σ |
|---|---|---|---|---|---|---|
| 135-136 | +171.3 | **+250.5** | 44.6 | – (135 not on the cross line) | – | – |
| **136-137** | +135.8 | **+61.6** | 20.6 | **+77.4** | **+15.8** | **+0.64** |
| 137-138 | +52.8 | **+32.2** | 13.6 | −21.6 | −53.8 | −2.71 |

**136-137 is resolved.** Its −3.22σ miss becomes **+0.64σ** the moment the two quantities are
about the same ground. 137-138 halves (−4.27σ → −2.71σ) but does not close. And the coefficient
itself moves by a factor of two or more on every pair: +135.8 → +61.6, +52.8 → +32.2,
+171.3 → +250.5.

Without any band at all — fit `D = k + κ·Y + c·dtan + m·dtan·Y` with `Y` the northing in km, so
`m = d(c)/d(along-track distance)` — and beside it the same two-parameter coefficient in six
equipopulated northing bands, south to north (six is `swath_across_track_test.py`'s own binning
of its raw medians, not a new choice):

| pair | kind | c at mid | **dc/dY** | SE | t | c in six northing bands |
|---|---|---|---|---|---|---|
| 135-136 | N-S/N-S | +166.3 | **+83.2** | 27.1 | **+3.1** | +141 +115 +23 +31 **+429** +274 |
| 136-137 | N-S/N-S | +134.7 | +13.3 | 9.4 | +1.4 | +8 +184 +172 **+352** −36 +158 |
| 137-138 | N-S/N-S | +51.8 | −7.7 | 7.7 | −1.0 | +39 +57 +161 −59 +65 +35 |

**The banded columns are the finding, and they need no model.** Over 3.5 km of one flight-line
pair, on one tile, the across-track coefficient ranges +23 to +429 (135-136) and −36 to +352
(136-137) — a spread five to ten times the ±9.8 to ±27.9 mm standard error the whole-pair fit reports.
The smooth linear gradient is significant on only one pair (+83.2 ± 27.1, t = +3.1); the rest of
the variation is not a ramp, it is structure.

**So the object being measured is not what the model calls it.** A between-line difference that
varies with across-track position is real — §3 of the across-track report established that on
eight pairs and it stands — but the coefficient describing it is a property of *the piece of
overlap it was fitted on*, not of the flight line. That is the honest reading of a −3.2σ and a
−4.3σ out-of-sample miss followed by a +0.64σ agreement as soon as the region is matched.

*(The cross pairs' own `dc/dY` — −304.3 ± 40.9 and −801.1 ± 77.0 — must NOT be read as an
along-track gradient. On a cross pair `tan_10010` is itself essentially northing (§3), so
`dtan·Y` is a quadratic in northing, not an interaction with position. It is reported in the run
for completeness and is not interpretable here.)*

---

## 7. The even/odd separation: possible in principle, not supported by this one line

`SWATH_ACROSS_TRACK_TEST.md` §0 records that inside a there-and-back N-S pair an even error
`c₂·tan²θ` collapses onto the odd predictor, because `tan²θ_A − tan²θ_B = (tan θ_A + tan θ_B)
(tan θ_A − tan θ_B)` and the first factor is pinned. A cross line's tangents are not pinned, so
the two forms are separable for the first time. Fitting both per line in the joint solve:

| line | c (odd) | SE | c₂ (even) | SE | t |
|---|---|---|---|---|---|
| 135 | −2865.0 | 632.1 | −7757.6 | 1445.0 | −5.4 |
| 136 | +110.7 | 52.6 | +486.3 | 161.3 | +3.0 |
| 137 | −29.6 | 7.6 | −43.9 | 49.7 | −0.9 |
| 138 | +63.1 | 83.4 | −16.2 | 237.0 | −0.1 |
| 10010 | −472.4 | 58.6 | −1564.5 | 231.3 | −6.8 |

**This is a negative and it is a conditioning one, not a physical one.** Within-cell r² rises
only 0.0472 → 0.0515 for five extra parameters, while the design condition number rises from
**16.9 to 192.4**, and line 135 — which this tile sees over a 0.0722 tangent span, one-sidedly —
returns −2865 ± 632 with an even term of −7758. Those are not measurements. Line 135 reaches
this tile only through its overlap with 136, whose tangent difference has sd **0.0722** against
0.1088 and 0.1145 on the other two N-S pairs, and line 10010's scan angles run **−17.00 to
+2.00°** — its nadir and its whole `+` side are in the tile to the north. **One cross line, on a
tile where both of those lines are sampled one-sidedly, does not support the even/odd split.**
It would need the cross line's own nadir track, which lies north in `4342-27-64`.

---

## 8. Sensitivity

Every coefficient re-derived under each variant; entries are the shift from the headline run, in
mm per unit tangent.

| variant | Δc₁₃₆ | Δc₁₃₇ | Δc₁₃₈ | Δc₁₀₀₁₀ |
|---|---|---|---|---|
| pairwise, `min_cell_line = 3` | −0.6 | +2.1 | +0.5 | −1.2 |
| pairwise, `min_cell_line = 5` | −7.9 | +3.4 | +5.7 | +3.2 |
| pairwise, quadratic position controls (E², N², EN) | +5.7 | −0.0 | −11.1 | −14.0 |
| joint, `min_cell_line = 3` | −4.6 | +5.9 | −2.6 | −5.2 |
| joint, `min_cell_line = 5` | −13.5 | +8.6 | −2.5 | −3.3 |

`min_cell_line = 1` is flagged `MINE`: it is the definitional floor — a median needs one value,
so it imposes no cut — and it is `swath_across_track_test.py`'s own default. Raising it to 3 or
5 moves every coefficient by at most 13.5 mm per unit tangent, and the quadratic position
controls by at most 14.0, against standard errors of 7.4–27.8
and a band-to-band spread of hundreds. **None of the parameters I chose is load-bearing for
anything above.** In particular they do not touch the −3.22σ and −4.27σ misses.

---

## 9. What this does NOT determine

**The global cross-track tilt of the gen1 mosaic is untouched, and the problem is not closed.**
`SWATH_DEGENERACY_BREAKING.md` §1 proves the null direction of the overlap network is
`e(x) = g·(x − x₀)`, a tilt identical for both lines at a shared ground point, so it cancels in
**every** between-line difference — adjacent, non-adjacent or crossing — to **0.000e+00 mm**.
Every coefficient in this document is a between-line difference, and the joint solve's cell
fixed effects remove exactly the same thing. **A cross line is blind to the tilt by
construction, and this run does nothing to change that.** The tilt requires an absolute external
reference; the instrument for it is ground control, where the two surveyed anchors 15.49 km
apart already give `g = −0.485 ± 1.916 mm/km` (that document §3, not re-derived here).

Also not determined: **line 133's coefficient in any useful sense.** Propagating the joint
solve down the elba+elbaext pair sums gives `c_134 = −6.8 ± 59.6` and `c_133 = +437.0 ± 71.1` —
each step doubles the previous error and adds twice the sum's, so the chain is numerically
useless after two links even before §6 removes the premise that the propagation rests on.

And not determined: **whether more cross lines exist.** `point_source_id` is not in the LAS
header, so each probe costs one tile download. 10010's own nadir track and entire `+` side lie
in `4342-27-64`, and it continues east into `4358-28-01` over lines 139–140 — the cheapest way
to test §7's conditioning problem and §6's locality at once.

---

## 10. The integration I propose — and it is *not* a per-line correction

**Nothing here should change `coreg.py` or `pipeline.py`, and specifically: do not add a
per-line across-track term.** That was the point of getting the cross line, and the cross line
is what argues against it. A `c_s` fitted on one overlap does not reproduce a neighbouring
overlap of the same two lines on the same tile (−3.22σ, −4.27σ), and the same coefficient runs
+23 to +429 across 3.5 km of a single pair. Applying such a term would import the local
structure of whichever overlap happened to be fitted into every cell of the swath. That is a
larger error than the one being corrected.

**Adopt instead what is already proposed, and for a reason this run strengthens.**
`analysis/SWATH_TIE_INTERCEPT.md` proposes `coreg.align_swaths(pc, ref=int(ps8.min()),
tie="intercept")`, estimating the tie at `dtan = 0` rather than as the overlap mean. **That
proposal survives everything found here, and this document makes its case stronger rather than
weaker:** the intercept tie does not assert that `c` is an instrument constant, it only removes
the tie's dependence on *where in the sidelap the tile happens to look*. If `c` is local — which
is what §6 measures — then the overlap-mean tie is contaminated by local structure in exactly
the way the intercept tie is not. I propose no change to that recommendation and nothing that
conflicts with it: **adopt `tie="intercept"`, record `swath_tie` in `corrections.json`, and keep
the gauge decision separate**, as written there.

**Two things worth doing that this run makes concrete, in order:**

1. **Fetch `4342-27-64` and `4358-28-01`** (~25 MB each). The first gives line 10010 its own
   nadir track and its `+` side, which is the single change that would fix §7's conditioning
   (135 and 10010 are both one-sided here) and let the even/odd split be decided. The second
   carries the cross line east over lines 139–142 and tests §6's locality on lines this tile
   never sees. Both are cheap and both are now well-posed questions rather than exploratory ones.
2. **Re-derive `SWATH_ACROSS_TRACK_TEST.md` §2** on the current `reference_cells`. Its
   137-138 coefficient changes sign under commit `7701383` (§2 above), and that pair is the one
   every anomaly in this line of work has converged on — the outstanding −8.9 mm of tile
   disagreement, the along-track contamination, the residual −2.71σ in §6, and the one pair
   where the two estimators of §2 differ by 76 mm.

**What I am not proposing, explicitly.** No new tie type, no per-line `c_s`, no shared roll
(rejected at p = 6×10⁻⁶⁴ by the across-track report and by `boresight.estimate_boresight`'s own
numbers), and no use of the cross line as a *vertical* constraint beyond what `align_swaths`
already does with any overlapping line. The cross line's value is diagnostic, and what it
diagnosed is that the model it was fetched to identify is the wrong model.
