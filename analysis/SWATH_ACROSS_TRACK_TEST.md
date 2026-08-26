# Is the gen1 per-swath vertical offset a constant, or a function of across-track position?

**Date:** 2026-08-26
**Script:** `analysis/swath_across_track_test.py`
**Run:** `env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/swath_across_track_test.py`
**Inputs:** `beam_offset_table.parquet`, `corrections.json`, `z_after.npy` and
`canopy_cover_pfs.npy` for `elba_fulldensity` and `elbaext`. Every input, parameter, mask and
column is declared through `trust/provenance.py`; the full banner prints with each run.
`coreg.py` and `pipeline.py` are imported, not modified.

## The question, and the form in which it can be falsified

`coreg.align_swaths` gives each gen1 flight line one constant `(dx, dy, dz)`. If the real
error instead varies with **across-track position**, that constant is an average over
whatever part of each overlap the tile happened to cover, and three otherwise unexplained
facts follow at once: the elba and elbaext tiles disagree by 8.0, 9.8 and 17.4 mm about the
same four swaths; both tiles pin their whole network to an edge-cut swath sampled at mean
scan −11.8° and −12.6°; and the two tiles see swath 135 at mean scan −11.78° and +0.65°.

Regressing the registered offset `d_mm_corr` on `scan_angle` within a swath cannot decide
it, and this run measures why rather than asserting it: within a flight line the ground
easting is `x = x_track + h·tan(scan)`, with r² = 0.95–0.997 for eight of the ten lines and
0.655 for the one edge-cut line (§10 of the run). Scan angle *is* across-track position, so
any error field that varies across the tile projects onto it.

**The overlap is the decisive test.** Where two lines cover the same cell they measure the
same ground at the same epoch, so terrain, land cover, vegetation, the gen2 reference and
any spatially varying field cancel *exactly* in the between-line difference `D`, while an
across-track error does not. Three hypotheses, separated by form:

1. **constant offset** – `D` is flat in scan angle.
2. **one shared coefficient** – `D = k_pair + c·(tan θ_A − tan θ_B)` with the *same* `c` on
   every pair and both tiles. That is what a mounting roll is: a pointing error `δ` at range
   `R = h/cos θ` displaces the point by `R·δ` normal to the beam, whose vertical part is
   `R·δ·sin θ = h·δ·tan θ`, so `c = h·δ` is a sensor constant.
3. **per-pair coefficients** – different `c` per pair; not an instrument constant.

**Reduction.** The per-cell, per-line ground estimate is the **median of `d_mm_corr`** – the
pipeline's own estimator (`ground_estimator = "slope_normal"`, `ground_q = 0.50`,
`pipeline.py:603` takes the `ground_q` quantile of the slope-normal residual per cell),
evaluated one flight line at a time. Nothing new is invented. `min_cell_line` defaults to 1,
the definitional floor; §9 of the run shows 3 and 5 move `c` by 0.0 and 0.3 mm per unit
tangent. Errors are cluster-robust on 50 m blocks, the repo default.

## 0. The flight geometry the test rests on, and what it can and cannot resolve

Measured, not assumed. Within every pair `tan θ_A + tan θ_B` is **constant**: its standard
deviation is 0.017–0.024 against a mean of 0.359–0.385 in magnitude. Both lines therefore see
the shared ground on the same body-fixed side, which happens only if consecutive lines fly
opposite directions. The fitted `h` in §10 confirms it independently, alternating sign line to
line (−2424, +2564, −2626, +2473 m at elba) because `scan_angle` is body-fixed and the
aircraft turns around.

Two independent routes give the flying height: median |h| from the fit is **2562 m**, and the
mean |sum_tan| of 0.369 is `S/h`, which with a nadir-track spacing of 786–991 m (mean 916)
gives 2482 m. A full swath is then `2h·tan17° = 1566 m` against that 916 m spacing, a **41%
sidelap 650 m wide** – the overlap this test uses. For swaths 136–138 the elba and elbaext
fits put the same line's nadir track within 3–8 m of each other.

Two consequences, both structural:

1. **Only the antisymmetric combination is identifiable.** With `tan θ_A + tan θ_B` fixed at
   `S/h`, separate coefficients on `tan θ_A` and `tan θ_B` are collinear by construction.
2. **An EVEN scan error mimics an odd one.** `tan²θ_A − tan²θ_B = (tan θ_A + tan θ_B)(tan θ_A
   − tan θ_B)`, so a symmetric error looks like a roll with a coefficient proportional to
   `S/h` – which **flips sign** between pairs that share the `+` side and pairs that share
   the `−` side. Sign alternation, not magnitude, is what tells the two apart.

## 1. Is it a constant? No.

| tile | pair | cells | ret A | ret B | 50 m blocks | mean scan A | mean scan B | dtan range | dtan sd |
|---|---|---|---|---|---|---|---|---|---|
| elba | 135-136 | 7,394 | 92,276 | 104,290 | 317 | −11.90 | −8.43 | −0.205…+0.073 | 0.0675 |
| elba | 136-137 | 9,817 | 132,886 | 130,191 | 440 | +10.79 | +10.58 | −0.218…+0.236 | 0.1062 |
| elba | 137-138 | 10,411 | 133,007 | 138,504 | 431 | −9.62 | −10.72 | −0.217…+0.251 | 0.1117 |
| elbaext | 133-134 | 7,218 | 91,548 | 100,992 | 327 | −12.71 | −7.92 | −0.218…+0.072 | 0.0637 |
| elbaext | 134-135 | 14,751 | 195,310 | 195,082 | 553 | +10.41 | +11.33 | −0.217…+0.199 | 0.1002 |
| elbaext | 135-136 | 14,582 | 193,232 | 196,568 | 567 | −9.96 | −10.44 | −0.212…+0.234 | 0.1049 |
| elbaext | 136-137 | 11,534 | 154,280 | 151,327 | 524 | +10.84 | +10.62 | −0.218…+0.236 | 0.1040 |
| elbaext | 137-138 | 12,030 | 153,246 | 159,505 | 501 | −9.55 | −10.74 | −0.217…+0.253 | 0.1129 |

`c` in mm per unit tangent (multiply by π/180 for mm/deg at small angle), cluster-robust SE:

| tile | pair | k (mm) | c | SE | mm/deg | t | c, terrain-grad controlled | c, along-track controlled |
|---|---|---|---|---|---|---|---|---|
| elba | 135-136 | −4.48 | **+125.2** | 18.0 | +2.185 | +6.9 | +131.8 | +115.1 |
| elba | 136-137 | −3.37 | **+169.6** | 8.0 | +2.959 | +21.3 | +168.4 | +167.3 |
| elba | 137-138 | +17.18 | **−43.8** | 21.2 | −0.764 | −2.1 | −43.5 | +34.5 |
| elbaext | 133-134 | −22.89 | **+190.5** | 17.7 | +3.325 | +10.8 | +194.5 | +192.7 |
| elbaext | 134-135 | +5.48 | **+62.6** | 6.6 | +1.093 | +9.5 | +62.8 | +57.2 |
| elbaext | 135-136 | +16.98 | **+133.4** | 9.3 | +2.328 | +14.3 | +134.8 | +126.8 |
| elbaext | 136-137 | +4.64 | **+183.7** | 9.1 | +3.206 | +20.1 | +185.0 | +172.1 |
| elbaext | 137-138 | +6.35 | **−34.0** | 18.0 | −0.593 | −1.9 | −33.6 | +36.7 |

**The between-line difference is not flat in across-track position.** It is not a fitting
artefact. The raw binned medians, six equipopulated `dtan` bins per pair with no model at
all, rise **strictly monotonically** on four of the eight overlaps – elbaext 135-136 runs
−0.7 → +40.8 mm, elbaext 136-137 −21.8 → +26.4 mm, elba 136-137 −27.5 → +19.7 mm, elbaext
133-134 −54.4 → −22.6 mm – and on two more apart from a single adjacent inversion well
inside the error bars, which run 1.5–7.5 mm. Pooled over both tiles with pair fixed effects,
adding one across-track term takes r² from 0.0584 to 0.0748, `c = +80.0 ± 5.7` (t = 14.1).

**The two 137-138 pairs are the only non-monotone ones, and they are informative.** They are
also the only pairs whose `dtan` correlates appreciably with *along*-track position
(r = −0.26, −0.28, against ≤0.14 elsewhere), their binned medians are peaked, and adding a
per-pair along-track control moves them from −43.8 → +34.5 and −34.0 → +36.7. With that
control every pair on both tiles has a **positive** across-track coefficient, +34.5 to
+192.7, and the pooled `c` rises to +100.6 ± 3.9 with r² 0.159.

## 2. Does one shared coefficient fit? No – and not for a spatial reason.

| model (both tiles pooled, pair fixed effects, 87,737 cell-pairs) | c | SE | t | r² |
|---|---|---|---|---|
| constant only, no scan term | – | – | – | 0.0584 |
| shared `c·(tan θ_A − tan θ_B)` | +80.0 | 5.7 | +14.1 | 0.0748 |
| shared `c·(θ_A − θ_B)`, degrees | +1.4 | 0.1 | +14.1 | 0.0747 |
| shared `c₂·(tan²θ_A − tan²θ_B)` (even) | +92.0 | 15.8 | +5.8 | 0.0614 |
| shared c, + terrain-gradient controls | +80.3 | 5.7 | +14.1 | 0.0754 |
| shared c, + per-pair along-track controls | +100.6 | 3.9 | +25.6 | 0.1590 |
| shared c, ALL in-grid overlap cells (506,705) | +108.0 | 3.8 | +28.6 | 0.0346 |

Cluster-robust Wald test of `H₀: all eight c_pair equal`: **W = 314.2, df = 7, p = 6×10⁻⁶⁴**;
with the along-track controls, W = 290.8, p = 6×10⁻⁵⁹. Letting each pair have its own
coefficient lifts r² from 0.0748 to 0.0953 for seven extra parameters. **A single shared
coefficient is rejected decisively.**

The even model is not the answer either. `c/sum_tan` – constant if the error were symmetric
in scan angle – runs −348, +448, +122, −522, +163, −370, +483, +95, and the shared even
model fits worse than the shared odd one (r² 0.0614 against 0.0748).

**But per-pair does not mean spatial.** Three independent checks say the heterogeneity is a
property of the flight lines:

* **It repeats across tiles.** Three pairs are estimated separately by both tiles, on
  different extents: +125.2 vs +133.4, +169.6 vs +183.7, −43.8 vs −34.0. Differences are
  0.40, 1.17 and 0.35 standard errors. *Caveat:* elba's footprint lies mostly inside
  elbaext's, so this is repeatability under a **change of extent** – which is exactly the
  failure mode at issue – not an independent replication on new ground.
* **It survives the terrain-gradient control.** A residual lateral misregistration *between*
  the two lines is the one error that does not cancel by being shared ground; it lands on
  `D` as `−(gx·Δδx + gy·Δδy)`. Adding the gen2 surface gradient from
  `registration.surface_gradients` – the same array `lateral_term` uses – moves `c` by at
  most 6.6 mm per unit tangent on any pair, and the pooled `c` from +80.0 to +80.3.
* **It is not the canopy.** Path length through vegetation grows as `1/cos θ`, so a
  penetration effect would give a coefficient near zero on open ground and growing with
  cover, the *same way in every pair*. The opposite is observed. In open cells (cover < 0.05)
  the coefficient is +167.0, +168.3, −48.2, +207.5, +58.2, +164.7, +205.6, −11.9 – larger in
  magnitude than the pooled value on five of the eight pairs, within 8% of it on two more,
  and smaller only on elbaext 137-138. It is never small where the pooled value is large,
  and under dense cover it is *smaller*, not larger (elbaext 135-136: +164.7 open against
  +60.0 dense).

So the across-track term is real, is carried at full strength by bare open ground, and
differs from flight line to flight line. It is a **per-line** across-track error – the same
object as a per-line roll, since a per-line error linear in cross-track ground coordinate and
a per-line roll are not distinguishable within a 650 m overlap strip.

**Per-swath coefficients cannot be recovered from overlaps alone.** Under a per-line model
`err_s = c_s·tan θ_s`, the fixed `tan θ_A + tan θ_B` collapses the pair fit to
`c_pair = (c_A + c_B)/2`. A chain of *n* lines gives *n*−1 such sums, so the `c_s` are
determined only up to adding an alternating vector `(+v, −v, +v, …)`. Breaking that needs a
loop – a non-adjacent overlap – or an external reference; neither tile has one.

**Magnitude.** With h = 2562 m the pair coefficients are equivalent roll errors of −3.5,
−2.7, and +5.0 to +15.3 arcseconds. Between two lines the difference sweeps 15–48 mm across
a single overlap. At the representative `c_pair` of 130 mm per unit tangent, and if the two
lines of a pair share the term equally, one line's own error runs ±40 mm between nadir and
±17°. That is the size of the cover-driven offsets being calibrated (−5 to −60 mm per cell)
and three times the 12.4 mm repeatability of the dz estimator reported in
`analysis/MISSION_TIME_DRIFT.md` §4.

## 3. tan against linear in angle: the data cannot separate them

Over the observed |scan| ≤ 17° the two predictors differ by at most 3.0% in shape and
correlate at **0.999992**. The pooled residual sums of squares are 3.22592×10⁸ (tan) against
3.22618×10⁸ (linear) – a **0.008%** difference. **This is a negative result: the functional
form is not identifiable from these data.** Prefer `tan` on mechanism, not on fit.

## 4. The elba / elbaext disagreement: two of the three, explained

`align_swaths` ties each line to its neighbour over their overlap, so if the difference
varies across that overlap the tie is that difference **averaged over the part of the overlap
the tile covers**. Two extents then get ties differing by `c·(mean dtan₁ − mean dtan₂)`, and
the disagreement in a swath's `dz` accumulates along the chain from the gauge swath 135.

| swath | dz elba | dz elbaext | disagreement | mean dtan elba | mean dtan ext | c_pair | step | cumulative | cum, along-track c | cum, shared c |
|---|---|---|---|---|---|---|---|---|---|---|
| 136 | −23.9 | −15.9 | **−8.0** | −0.0626 | +0.0031 | +129.3 | −8.50 | **−8.50** | −7.95 | −5.26 |
| 137 | −32.5 | −22.7 | **−9.8** | −0.0015 | −0.0009 | +176.6 | −0.11 | **−8.61** | −8.06 | −5.31 |
| 138 | −43.7 | −26.3 | **−17.4** | +0.0102 | +0.0136 | −38.9 | +0.13 | **−8.48** | −8.18 | −5.58 |

Means are over all in-grid overlap cells. (`coregister_swaths` uses the whole overlap at 2 m
on the vendor terrain classes rather than 5 m on CSF ground, so this is a close proxy, not
the identical population.)

The mechanism accounts for swath 136 to **0.5 mm** and swath 137 to **1.2 mm**, and it
names the cause. elba's western boundary (x₀ = 577,492.8) lies ~350 m east of swath 135's
fitted nadir track (577,143), so elba sees that line only from |scan| 8–16° and samples the
135-136 overlap at mean dtan −0.063; elbaext, whose western edge is 1,893 m further west,
samples the same overlap at +0.003. Almost the whole disagreement enters at that first link
and is then carried down the chain. Swath 138 is **not** explained:
−17.4 mm observed against −8.5 mm predicted, leaving −8.9 mm outstanding. A single shared
coefficient would explain only −5.3 of the −8.0 mm and none of the rest.

## 5. The edge swaths: how much of +6.77 and +4.27 mm/deg is identifiable

Per-return regression of `d_mm_corr` on `scan_angle`, same reference cells. The `before`
column reproduces the tabulated residual slopes at the precision they were quoted, so this
is the same object:

| tile | swath | returns | mean scan | scan sd | slope mm/deg | SE | t | corr(const, slope) | SE of fit at scan=0 | SE of fit at own mean | slope after shared c |
|---|---|---|---|---|---|---|---|---|---|---|---|
| elba | **135** | 92,322 | −11.78 | 1.81 | **+6.768** | 2.110 | +3.2 | **+0.989** | **25.3** | 3.8 | +5.308 |
| elba | 136 | 398,151 | +1.40 | 7.88 | +0.121 | 0.205 | +0.6 | −0.325 | 1.8 | 1.7 | −1.294 |
| elba | 137 | 423,019 | +0.41 | 8.33 | −1.402 | 0.247 | −5.7 | −0.296 | 2.1 | 2.0 | −2.818 |
| elba | 138 | 245,638 | −6.69 | 5.04 | +3.261 | 0.713 | +4.6 | +0.824 | 6.6 | 3.8 | +1.835 |
| elbaext | **133** | 91,554 | −12.64 | 1.94 | **+4.265** | 1.720 | +2.5 | **+0.991** | **22.3** | 3.1 | +2.798 |
| elbaext | 134 | 544,474 | +2.35 | 7.17 | −1.569 | 0.197 | −8.0 | −0.254 | 1.4 | 1.4 | −2.983 |
| elbaext | 135 | 633,137 | +0.65 | 8.66 | −1.671 | 0.157 | −10.7 | −0.133 | 1.3 | 1.3 | −3.088 |
| elbaext | 136 | 546,888 | −0.63 | 8.87 | −0.017 | 0.146 | −0.1 | −0.018 | 1.3 | 1.3 | −1.435 |
| elbaext | 137 | 488,855 | +0.52 | 8.34 | −1.746 | 0.223 | −7.8 | −0.272 | 2.0 | 1.9 | −3.163 |
| elbaext | 138 | 278,347 | −6.79 | 5.02 | +1.634 | 0.734 | +2.2 | +0.877 | 7.1 | 3.6 | +0.207 |

**What is identifiable** for elba/135: the end-to-end change of the fitted line across its own
observed range, 6.768 × 8° = **54 ± 17 mm** (3.2σ). Its neighbour-pair coefficient from the
overlap, `(c₁₃₅ + c₁₃₆)/2 = +2.185 ± 0.314 mm/deg`, is measured cleanly and independently of
any spatial field. For elbaext/133 the same pair quantity is `+3.325 ± 0.309 mm/deg`.

**What is not identifiable** is the split into a constant plus a slope. The constant and slope
estimates correlate at **+0.989** and **+0.991**. The fitted level is pinned to 3.8 mm and
3.1 mm *where the swath was observed* and only to **25.3 mm and 22.3 mm** extrapolated to
scan = 0, where a per-swath constant places it – a 6.7× and 7.2× inflation. Equivalently, any
constant shift of ±25.3 mm trades against a slope change of ∓2.15 mm/deg, which is **32%** of
the +6.77; for elbaext/133, ±1.76 mm/deg, **41%** of the +4.27. Both tiles pin their whole
alignment network to exactly this swath.

Subtracting a single shared roll does not clean these up: it shifts every slope by about
−1.4 mm/deg and leaves the spread intact (+5.31 to −3.16 mm/deg). Only a per-line term could,
and §2 shows per-line terms are not uniquely recoverable from overlaps alone.

One further sign the within-swath regressions are contaminated: the overlap gives
`(c_A + c_B)/2` cleanly, and reading the same quantity off the gen2-referenced within-swath
slopes disagrees by −4.09 to +1.69 mm/deg, negative on five of the eight pairs. The
gen2-referenced regression systematically *under*-reads the across-track term.

## 6. How this stands against the earlier boresight exclusion

`analysis/ridgelines/GLENNIE_SCANANGLE_SWATH_TEST.md` excluded a boresight / scan-mirror
artefact as the driver of the steep-slope gen1-low, on three grounds: the sign was wrong
(residual improved toward the edge), the effect was nadir-worst, and the edge-nadir delta was
incoherent between swaths.

**That conclusion is not contradicted, and this run reinforces its central point.** Two things
must be said plainly:

1. **The earlier test could not have seen a roll, by construction.** It binned on
   `|scan angle|`, unsigned, and its per-swath table compares `nadir<5` against `edge>15` in
   absolute value. A roll error is *odd* in scan angle: it reads high on one side and low on
   the other, and averages to zero under `|scan|`. It also worked on `d_mm`, before
   registration, against gen2 rather than between lines. So it tested a different quantity in
   a frame blind to this one. Its verdict on a symmetric, edge-worst artefact stands; it never
   bore on a signed across-track term in the inter-swath offset.
2. **This run reaches the same verdict about a shared instrument roll, with a sharper
   measurement.** A single sensor-constant roll is rejected at p = 6×10⁻⁶⁴. The earlier test's
   third ground – "between-swath edge behavior is incoherent" – is precisely what is found
   here: the coefficients are per-line, +34.5 to +192.7 mm per unit tangent, and homogeneity
   is rejected at p = 6×10⁻⁵⁹. The repo's own `boresight.estimate_boresight` says the same
   thing in its own output (§8 of the run): pooled b = +2.19 to +2.27 mm/deg with a bootstrap
   SE of 0.02–0.04 but a **between-pair** standard deviation of 0.71–1.06 mm/deg, 25–45×
   larger. That is why `boresight_roll_mm_per_deg: None` sits in the elba corrections file
   (elbaext's carries no boresight key at all), and it should stay there.

The disagreement with the earlier document is narrow and worth stating: it read the
across-track structure as "a mild co-registration / flight-line offset, not a scan-angle
ramp." The overlap test says it *is* a scan-angle ramp, just a per-line one rather than a
sensor-wide one, and it is large enough to matter.

## 7. Verdict

* **The per-swath offset is not a constant.** In the overlap, where terrain, cover and any
  spatial field cancel exactly, the between-line difference sweeps 15–48 mm with across-track
  position. Pooled, `c = +80.0 ± 5.7` mm per unit tangent (t = 14.1), rising to +100.6 ± 3.9
  once per-pair along-track structure is controlled, and +108.0 ± 3.8 on all 506,705 in-grid
  overlap cells.
* **One shared coefficient does not fit.** W = 314.2, df = 7, p = 6×10⁻⁶⁴. Nor does a shared
  even (symmetric) term. The coefficients are per flight-line pair.
* **Per-pair here does not mean spatial.** It repeats across tile extents to <1.2σ, survives
  the terrain-gradient control to 6.6 mm per unit tangent, and is carried at full strength by
  bare open ground where no canopy exists. It is a per-line across-track term – equivalently,
  a per-line roll of −3.5 to +15.3 arcseconds at the measured 2562 m flying height.
* **tan against linear is a negative result.** The two forms differ by 0.008% in fit; the data
  cannot choose. Take `tan` on mechanism.
* **Two of the three tile disagreements are accounted for.** Swaths 136 and 137 are predicted
  to 0.5 and 1.2 mm from the tiles' own across-track sampling of the 135-136 overlap. Swath
  138 is not: −8.5 mm predicted against −17.4 mm observed, leaving −8.9 mm outstanding.
* **The edge swaths cannot support the constant/slope split they are being asked for.** Their
  constant and slope correlate at 0.99, and their level is known to 3–4 mm where observed and
  22–25 mm where the model places it. Both tiles gauge their whole network on such a swath.

## 8. What I would propose, and what I have not done

I have changed nothing in `coreg.py` or `pipeline.py`. Three things follow, in order of
confidence:

1. **Stop gauging the network on an edge-cut swath.** `pipeline.difference_dem` calls
   `align_swaths(ref=int(ps8.min()))` at `pipeline.py:670` – the lowest-numbered line, which
   at both Elba tiles is the one whose nadir track falls outside the tile (elba/135 at
   577,143 against x₀ = 577,492.8; elbaext/133 at 575,340 against x₀ = 575,600), so it is
   sampled over an 8–10° one-sided scan range. Gauging on the line with the *widest two-sided*
   scan coverage, or on the zero-mean gauge `align_swaths` already supports, costs nothing and
   removes the worst-conditioned reference.
2. **Fit the pairwise tie at a stated across-track position, not as an unweighted mean.**
   The tie's extent-dependence is entirely the mean `dtan` the tile happens to sample.
   Estimating `k` at `dtan = 0` – the fitted intercept, which this script already produces –
   makes the tie an extent-invariant quantity. It would have removed the −8.50 mm step at the
   135-136 link, which is what generates the −8.0 mm disagreement at swath 136 and carries
   down the chain to 137.
3. **Do not add a shared roll term.** It is rejected here, `estimate_boresight` rejects it in
   the repo's own numbers, and applying `c_shared` leaves the within-swath slopes as spread as
   it found them. A per-line across-track term is what the data support, and it is **not
   uniquely recoverable** from a chain of overlaps – it needs a loop, a cross-line, or an
   external reference. Building one from these tiles alone would be fitting the alternating
   gauge, not the instrument.

The open item is swath 138's residual −8.9 mm, which this mechanism does not reach.
