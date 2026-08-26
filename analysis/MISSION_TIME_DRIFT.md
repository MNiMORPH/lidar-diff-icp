# Is the gen1 per-swath offset one continuous drift in mission time?

**Date:** 2026-08-26
**Scripts:** `analysis/mission_time_drift.py` (the test), `analysis/mission_time_drift_fit.py` (the model)
**Run:** `env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/mission_time_drift.py`
**Inputs:** the seven gen1 LAZ files the products were built from, and each site's
`data/derived/<site>/corrections.json`. Both scripts declare every input, parameter and
column through `trust/provenance.py`; the full banner is printed with each run.

## The hypothesis, and the form in which it can be falsified

The pipeline models two things per flight line: a constant `(dx, dy, dz)` from
`coreg.align_swaths`, and a within-swath curve in `gps_time` from
`coreg.fit_along_track_drift`. Andy's proposal is that this split is artificial – that
there is one continuous position error **F(t)** over mission time, and the per-swath
constants are that function evaluated during each line's 45-second window.

The supporting observation is real: ordered by time, elba's dz runs −43.7, −32.5, −23.9,
0.0 mm and elbaext's runs −23.6, −20.0, −13.2, +2.7, +18.8 mm, both monotonic, while
consecutive lines are separated by 1,000 to 12,300 seconds.

A correlation with `gps_time` cannot decide it. Within a tile, mission time, flight-line
ordinal and across-track position are collinear: lines are flown in order and laid down
side by side, so "dz rises with `gps_time`" is the same statement as "dz rises with line
number." The gap structure is what separates them. Consecutive-line gaps in this dataset
run from 360 s to 70,974 s, a factor of 197. A continuous F(t) predicts the step between
two lines scales with the **elapsed time** between them; a per-line quantity predicts it
scales with the **number of lines**. That is the test.

## 1. The per-site slopes

Reference swath excluded throughout: `pipeline.difference_dem` calls
`align_swaths(ref=int(ps8.min()))`, which pins that swath to `(0, 0, 0)`. It is a gauge
choice, not an observation. Correlations with it included are printed in the script's
`r_incl_ref` column; they are systematically weaker for dz (elba 0.841 rather than 0.977),
because in four of the seven sites the pinned swath is the last in time and its forced
zero sits at one end of the range.

| site | swaths fitted | dz r | dz mm/h | dx r | dx m/h | dy r | dy m/h |
|---|---|---|---|---|---|---|---|
| elba | 3 | +0.977 | **+8.1** | −0.989 | **−0.304** | −0.833 | −0.094 |
| elbaext | 5 | +0.819 | **+12.0** | −0.973 | **−0.376** | −0.463 | −0.042 |
| whitewater | 3 | −0.951 | **−74.7** | −0.932 | **−0.585** | −0.732 | −0.230 |
| mnrv | 6 | +0.338 | +12.1 | −0.610 | −0.261 | +0.527 | +0.109 |
| cook | 3 | −0.854 | **−62.9** | +0.998 | **+0.995** | −0.995 | **−1.291** |
| carlton | 4 | +0.999 | **+46.9** | +0.460 | +0.028 | +1.000 | **+0.146** |
| battlecreek | 4 | −0.255 | −0.3 | −0.752 | −0.004 | −0.764 | −0.003 |

**Do dx and dy drift too?** Yes, and far more strongly than dz in absolute terms. The
lateral offsets reach 1.39 m (elbaext swath 138) and 0.90 m (cook swath 11), and they are
more nearly monotonic than dz is. Furthermore they are **across-track**. Every site's
flight heading, fitted per swath from x(t) and y(t), shows lines alternating there and
back along one axis: elba, elbaext, whitewater and mnrv fly north–south, so the offset
lands almost entirely on dx; cook flies east–west, so it lands on dy as well; carlton
flies on a −40° axis and splits between the two. Rotating (dx, dy) into each site's own
track frame confirms it – the across-track component carries essentially the whole signal
(elba: across −0.306 m/h against along −0.089 m/h).

**Is the slope consistent across sites?** No. dz spans +46.9 to −74.7 mm/h and changes
sign; the lateral rate spans −0.585 to +0.995 m/h and changes sign. There is no shared
instrument property here. The sharpest case is internal to one acquisition: elba, elbaext
and whitewater all sit on **GPS day of week 2**, with disjoint windows (whitewater
226,000–229,475 s; elba/elbaext 237,350–249,711 s), so if they share a GPS week they are
one flight day. Their dz rates are +8.1, +12.0 and **−74.7** mm/h. One F(t) covering that
day would have to reverse sign in the 7,900-second gap where no returns exist. The GPS
week is not recorded in these LAS 1.1 files (`global_encoding` bit 0 = 0, week seconds),
so the same-day reading is consistent with the data but not established by it.

## 2. The falsifiable prediction, and where it fails

Site-by-site steps in dz between consecutive lines, against what each site's own fitted
rate predicts for that gap:

| site | pair | gap (s) | Δdz (mm) | implied mm/h | F(t) predicts (mm) |
|---|---|---|---|---|---|
| elba | 138→137 | 6,290 | +11.2 | 6.4 | +14.1 |
| elba | 137→136 | 1,998 | +8.6 | 15.5 | +4.5 |
| elba | 136→135 | 1,022 | +23.9 | 84.2 | +2.3 |
| elbaext | 138→137 | 6,297 | +3.6 | 2.1 | +21.1 |
| elbaext | 136→135 | 1,029 | +15.9 | 55.6 | +3.4 |
| elbaext | 134→133 | 1,060 | −18.8 | −63.8 | +3.5 |
| battlecreek | 1013→1014 | 689 | **−36.7** | −191.8 | −0.1 |
| battlecreek | 1014→1101 | **70,974** | **+10.0** | 0.5 | −6.8 |

Three results, each independent of the others.

**(1) The steps do not scale with elapsed time.** At elba the longest gap (6,290 s) gives
the smallest step (+11.2 mm) and the shortest gap (1,022 s) gives the largest (+23.9 mm) –
the reverse of the prediction, and a factor of 13 discrepancy in implied rate across three
pairs of the same site. Expressed as a coefficient of variation over each site's steps,
per-line is the tighter description of dz at five of the seven sites (elba 0.46 against
0.98; whitewater 0.12 against 0.44; cook 7.75 against 53.50).

**(2) Battle Creek's 19.7-hour mission break produces the smallest step on the site.**
Lines 1012–1014 were flown within 25 minutes of each other on GPS day 0 and differ by
+30.1 and −36.7 mm. Line 1101 was flown 70,974 s later on GPS day 1 and sits +10.0 mm from
1014. A continuous drift of any rate that explains 36.7 mm in 689 s cannot then produce
10.0 mm in 70,974 s. This is the one case in the dataset where the two candidate
explanations are cleanly separated by the design of the flight, and it goes against
mission time.

**(3) A two-parameter model in line ordinal beats a two-parameter model in elapsed time.**
Same parameter count, same data, R² on dz:

| site | R² linear in gps_time | R² linear in line ordinal |
|---|---|---|
| elba | 0.955 | **0.994** |
| elbaext | 0.671 | **0.930** |
| whitewater | 0.904 | **1.000** |
| mnrv | 0.114 | **0.190** |
| cook | 0.730 | **0.865** |
| carlton | **0.997** | 0.934 |
| battlecreek | 0.065 | **0.250** |

Ordinal wins at six of seven sites for dz. Carlton is the single exception and it is a
real one, not noise: among carlton's three steps that do not involve the pinned reference
swath, the gaps alternate 1,622 / 416 / 1,673 s and the dz steps follow them (+22.4 / +3.0
/ +23.2 mm), which is what a continuous drift looks like. The fourth step, 87→86 over
360 s, is +16.4 mm and does not follow, but 86 is the gauge swath and is excluded from the
fit. Carlton is worth revisiting on its own; it is not enough to carry the other six.

## 3. The model, and what it costs

`analysis/mission_time_drift_fit.py` fits **F(t) = a + b·t** per axis to the saved
per-swath constants. It fits to the saved constants rather than re-running
`align_swaths`, because five of the seven sites were aligned on a CSF-classified cloud
whose cache (`data/csf_cache/<site>.las`) no longer exists, and `coregister_swaths`
excludes classes (5, 6, 9) which CSF overwrites. A re-run on the vendor classification
would be a different method, and the two models would no longer be compared the same way.
Fitting to the constants keeps the method identical and makes the residual directly the
per-swath vertical bias, in mm, that the substitution would introduce.

**Parameter count.** Per axis, the current model has *n* swaths minus one for the gauge:
3 (elba, whitewater), 4 (carlton, battlecreek), 5 (elbaext), 6 (mnrv). Linear F(t) has
one free parameter after the gauge, at every site. The reduction is genuine, and it is a
handful of numbers.

**Residual.** RMS misfit of F(t) to the saved constants, and of the same-cost ordinal
model:

| site | swaths | dz RMS, F(t) | dz RMS, ordinal | dz worst residual, F(t) | site stable 1σ |
|---|---|---|---|---|---|
| elba | 3 | 1.7 mm | 0.6 mm | 2.3 mm | 59.9 mm |
| elbaext | 5 | 9.0 mm | 4.2 mm | 11.8 mm | 52.4 mm |
| whitewater | 3 | 8.2 mm | 0.1 mm | 11.0 mm | 80.4 mm |
| mnrv | 6 | 16.7 mm | 16.0 mm | 30.4 mm | 60.6 mm |
| cook | 3 | 8.4 mm | 6.0 mm | 11.7 mm | 70.5 mm |
| carlton | 4 | 0.9 mm | 4.4 mm | 1.2 mm | 44.4 mm |
| battlecreek | 4 | 13.0 mm | 11.7 mm | 18.4 mm | 32.2 mm |

With three observations a two-parameter line has one residual degree of freedom, so the
small numbers at elba, whitewater and cook are close to automatic; elbaext and mnrv carry
the weight. The bar these have to clear is the estimator's own repeatability, measured
below at 12.4 mm RMS on dz. F(t)'s residual reaches or exceeds that bar at mnrv (16.7 mm)
and battlecreek (13.0 mm) and comes close at elbaext (9.0 mm, worst single swath 11.8 mm).
The ordinal model is the lower of the two at six of the seven sites.

## 4. The elba/elbaext anomaly

The two Elba tiles align the same four flight lines from two different gen1 extents (a
2.5 × 3.5 km single tile against a 4.45 × 4.05 km merge). Re-referenced to swath 135 they
disagree, in mm:

| axis | sw135 | sw136 | sw137 | sw138 | RMS |
|---|---|---|---|---|---|
| dx | 0.0 | +2.2 | −7.6 | +3.3 | **4.9** |
| dy | 0.0 | −24.1 | −15.2 | +6.8 | **16.9** |
| dz | 0.0 | +8.0 | +9.8 | +17.4 | **12.4** |

**The mission-time explanation of this anomaly is refuted quantitatively.** Under the
hypothesis, each tile's constant averages F(t) over that tile's along-track segment, so
the disagreement is a claim about the tiles' `gps_time` midpoints. Those midpoints differ
by **+3.98, −2.87, +3.10 and −3.39 seconds** for swaths 135–138. To produce +17.4 mm from
−3.39 s requires **5.1 mm/s**, against elba's own fitted rate of **0.0022 mm/s** – short
by a factor of **2,300**. The tiles' segments differ by seconds, not by the hours the
explanation needs.

What the anomaly is instead: this estimator's repeatability under a change of extent.
`coregister_swaths` runs Nuth & Kääb on whatever overlap the two tiles share, and the
larger merge samples different terrain. dx, the axis carrying the metre-scale offsets,
repeats to 4.9 mm; dz repeats to 12.4 mm; dy is worst at 16.9 mm.

**Does F(t) make the tiles agree?** Partly, and not for the reason proposed. Fitting each
tile separately and differencing the fitted curves at the shared swaths gives 6.3 mm RMS
for F(t) against 12.4 mm as saved – but a two-parameter model has less freedom to
disagree than five constants do, so some of that reduction is arithmetic rather than
evidence. The informative comparison is against the equally-cheap ordinal model, which
gives **1.8 mm**. Smoothing removes the anomaly; mission time is not what does the
smoothing.

## 5. What the gaps cannot constrain

Returns exist in 45–60 second windows separated by hours. The fraction of each
inter-swath interval that carries any data at all runs from **11.5%** (carlton 89→88, a
416 s gap) down to **0.02%** (battlecreek 1014→1101, the 19.7-hour break). Between the
windows F(t) is set entirely by its assumed form. Three consequences:

1. **The functional form is doing the work, and must be stated.** Linear is the only form
   the gaps can support without inventing structure. Anything more flexible interpolates
   through intervals where 99% of the span is unobserved.
2. **The absolute level is degenerate.** The whole curve floats vertically; only
   differences between swaths are constrained. That free constant is exactly the quantity
   the ground-control chain work is separately trying to pin, and it should stay free
   here.
3. **One curve cannot carry both timescales.** The within-swath along-track drift already
   fitted by `coreg.fit_along_track_drift` has a peak-to-peak of 7.1 to 148.0 mm inside a
   single 45-second window, which is **524 to 10,638 mm/hour**. The between-swath rates
   fitted above are 0.3 to 74.7 mm/hour. The two differ by two to three orders of
   magnitude. A single linear F(t) therefore does not replace both terms of the current
   model; it replaces the per-swath constants only, and the within-swath curve survives
   unchanged.

## 6. Verdict

**The hypothesis holds nowhere in the form proposed, at one site in a weaker form, and
its central prediction fails everywhere it can be tested.**

- dz correlates with `gps_time` at five of the seven sites, |r| from 0.82 to 0.999. That
  much of Andy's observation is confirmed and is not a two-tile accident.
- The correlation does not survive the gap test. Steps do not scale with elapsed time
  (all sites but carlton); a 19.7-hour mission break produces the smallest step on its
  site (battlecreek); and a same-cost model in line ordinal fits dz better at six of
  seven sites.
- The rate is not a shared instrument property. It spans +46.9 to −74.7 mm/h and reverses
  sign, including within what may be a single flight day (elba +8.1 against whitewater
  −74.7).
- The elba/elbaext disagreement is not a mission-time effect. It needs 5.1 mm/s where the
  site supplies 0.0022 mm/s, and it is fully accounted for as the estimator's own
  extent-dependent repeatability.
- Carlton is the one site whose dz steps do track their gaps (R² 0.997 in time against
  0.934 in ordinal). It is a single site and it should be checked on its own rather than
  read as partial support.

**What survives, and what I would propose next.** The per-swath offsets are strongly
organised, just not by clock time. They are monotone in flight-line index and they lie
across-track. At Elba the lateral step per line is of order 0.3 m (elba 0.320, 0.276,
0.457 m; elbaext 0.254, 0.079, 0.322, 0.266, 0.468 m), and the dx values themselves
reproduce between the two tiles to 4.9 mm RMS. Line index and across-track position are not separable within a
single tile, so the next question is which of the two it is – and that is answerable, by
comparing tiles whose line spacing or line numbering differ.

I have not changed `coreg.py` or `pipeline.py`. If the ordinal/across-track structure
holds up under that test, the integration to propose is a per-axis linear term in
across-track line position replacing the free per-swath constants, which would cut the
alignment from *n*−1 parameters per axis to one and would carry a physical claim about
the flight geometry. Replacing them with a linear function of `gps_time` is not supported
by this evidence, and at elbaext, mnrv and battlecreek it would leave a per-swath bias at
or above the estimator's own noise floor.
