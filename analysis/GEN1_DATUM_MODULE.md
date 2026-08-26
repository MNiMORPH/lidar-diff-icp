# `gen1_datum` – gen1's absolute vertical datum at any Minnesota site

**Date:** 2026-08-26
**Module:** `src/lidar_diff_icp/groundtruth/gen1_datum.py`
**Tests:** `tests/test_groundtruth_gen1_datum.py` (38)
**Driver:** `analysis/groundtruth/gen1_datum_at_site.py`
**Manual:** `docs/groundtruth.md` §5a
**Supersedes as method, not as result:** the throwaway `screen_marks.py` in the session
scratchpad, which measured 43 marks with the line taken from distance to a fitted
centreline
**Written in parallel with:** `analysis/GEN1_DATUM_MORE_MARKS.md` and
`analysis/groundtruth/gen1_more_marks_{tie,report}.py`, which do the same measurement as
scripts inside lines 133-138. The two agree where they overlap (§4, §7); folding those
scripts onto this module is an open item, not done here.

Every number below is pasted from a command run while writing this file. Where a number
from another document could not be re-derived it is marked `UNVERIFIED` and its source
named. The parameters I chose are §6 and each carries its measured effect.

---

## 1. Bottom line

1. **The line assignment was the thing worth fixing.** Taken from the `point_source_id`
   of the ground returns instead of distance to a centreline, the assignment changes for
   **9 of 31** marks around Elba, and the change is entirely in the far field: **9 of 9
   agree within 10 km, 8 of 12 at 10–15 km, 5 of 10 at 15–25 km**. On those same 31 ties
   the per-line datum estimate moves from **−44.8 mm to −7.8 mm** – a **37.0 mm** swing
   produced by labelling alone.

2. **The bundled control has 963 marks, not 1 004.** 41 rows are the same physical mark
   printed in two counties' validation reports; 39 groups, one of them a triple. The
   named case `L2T-6126 Wabasha RTK` is one of them.

3. **gen1 at Elba, measured three ways, with what each is:**

   | what | value | SE | n marks | k lines |
   |---|---|---|---|---|
   | all covers, 20 km, `per_line` | **−55.5 mm** | 18.3 | 56 | 27 |
   | all covers, 20 km, `common_datum` (elbaext frame) | **+44.4 mm** | 30.5 | 14 | 6 |
   | L1O+L5U only, 20 km, `common_datum` (elbaext frame) | **+53.7 mm** | 27.2 | 7 | 4 |
   | all covers, 20 km, `common_datum` (elba frame) | **+72.3 mm** | 33.5 | 7 | 4 |

   Sign: the constant to **ADD to gen1**, positive = gen1 reads low. Every SE is the SE
   of *the mean over flight lines of the within-line mean tie*. The `per_line` and
   `common_datum` numbers are **not** the same quantity: `per_line` averages over 27
   lines whose constants are unknown, `common_datum` reports the datum of one named
   swath frame using only the marks that frame can hold.

4. **The scatter is organised by flight line, confirmed on the returns-based grouping.**
   Over the 56 marks within 20 km: **F = 3.37, p = 0.000952** (df 26, 29), ICC 0.537.
   Pooling the marks as independent would understate the SE by **1.42×**.

5. **Negative result: the radius-spread screen does not reduce the site-to-site scatter.**
   Over those 56 marks the sd of the per-mark ties runs **105.6 / 95.4 / 92.6 / 92.9 /
   96.5 mm** at cuts of 15 / 25 / 50 / 100 mm and no cut. There is no cut in the module
   and there should not be one until something is shown to buy anything.

6. **The swath network's per-line residuals against ground truth are large.** In the
   elbaext frame, RMS **68.3 mm** over six lines, worst **+125.6 mm** on line 135 (one
   mark). With three of the six lines carrying one mark each, this is not yet a test the
   network fails – it is a test that cannot be run at this n.

---

## 2. The API

```python
from lidar_diff_icp.groundtruth import gen1_datum as G

control = G.load_control()                       # 1004 rows -> 963 marks, merges reported
G.assert_no_geoid_conversion(control)            # raises unless both sides are GEOID03

sites = G.discover_near_point(control, easting, northing, radius_m)     # no default radius
sites = G.discover_near_lines(control, {137: ((x0, y0), (x1, y1))}, half_width_m)

res   = G.resolve_tiles(sites, ["data/before"], cache="data/mn_tile_centroids.csv")
res.on_disk, res.to_fetch                        # NOTHING is downloaded

const, src = G.swath_constants_from_corrections("…/corrections.json")   # optional
meas, skipped = G.measure_sites(sites, res, swath_constants=const,
                                swath_constants_source=src)
est = G.combine_datum(meas, mode="common_datum", swath_constants_source=src)
print(est.summary()); est.to_dict()
```

Returned objects, and what each carries:

| object | the numbers on it |
|---|---|
| `ControlSet` | `.marks`, `.n_rows`, `.merges`, `.merge_note` |
| `MarkSite` | the mark and its distance to whatever the search was about |
| `TileResolution` | `.per_mark`, `.on_disk`, `.to_fetch`, `.table_rows()` |
| `LineAssignment` | `.counts` (every `point_source_id` with its count), `.dominant`, `.dominant_fraction`, `.n_lines`, `.mixed` |
| `SitingScreen` | `n`, `slope_deg`, `relief_mm`, `fit_rms_mm`, `radius_spread_mm` |
| `MarkMeasurement` | `.tie`(the full `TieEstimate` with its radius ladder), `.line`, `.screen`, `.swath_shift_m`, `.per_line_tie_mm`, `.params`, `.notes` |
| `Gen1DatumEstimate` | `.value_mm`, `.se_mm`, **`.se_of`**, `.groups`, `.anova_F/p/df`, `.icc`, `.mean_over_marks_mm`, `.se_over_marks_mm`, `.design_effect`, `.line_residual_mm`, `.excluded` |

`table_columns()` on `TileResolution`, `MarkMeasurement`, `SitingScreen` and
`Gen1DatumEstimate` returns *name → definition with units*, which the driver feeds
straight into `trust.provenance.Run.column`. A column cannot be printed that the code has
not defined.

---

## 3. Design decisions, and why

**(a) The flight line comes from the returns.** `assign_line_from_returns` counts
`point_source_id` inside the estimator's report radius. Line spacing here is ~1 km and
the nadir tracks were fitted at one latitude, so a 1–3° heading error walks a track
several hundred metres over 20 km. The returns are the acquisition's own record of which
sortie lit the ground and cannot be wrong about it. §4 measures what this bought.

**(b) The flight line is the unit of replication.** Marks under one line share that
swath's unknown constant, so `value = mean over lines of (mean over that line's marks)`
and `SE = sd(line means)/sqrt(k)`. `se_of` carries that sentence as data, so the number
cannot be printed without it. The per-mark SE, the ANOVA, the ICC and the design effect
come back beside it, because the ratio is the size of the mistake the independence
assumption makes.

**(c) The mode is checked, not trusted.** `combine_datum(mode="per_line")` raises if the
measurements were made with swath constants applied, and `mode="common_datum"` excludes
(with a reason, and a count) any mark whose line has no constant. A per-line average must
not be labellable as a common-frame one.

**(d) The geoid is asserted, never converted.** `assert_no_geoid_conversion` compares each
mark's recorded model against the lidar's and raises on a mismatch; it *returns* the
sentence a run should print, so "no geoid conversion" appears in the output as a checked
statement rather than an author's assurance. It runs inside `measure_site`, so it is not
possible to reach a tie without it.

**(e) The gen2-derived lateral shift is opt-in.** `lateral_shift_m=None` by default. gen2
is not in a gen1-against-its-own-control comparison, so its Nuth & Kääb shift has no place
in one. Passing it is recorded as `src="andy"` in the params and noted on the measurement.

**(f) Nothing is downloaded and nothing is cut.** `resolve_tiles` reports; a test asserts
`tiles.download_tile` is never called. Every screen statistic is returned and there is no
threshold, no minimum `n` and no distance default anywhere in the module.

**(g) The duplicates are merged in code, and the CSV is left as a faithful
transcription.** A mark on a county line genuinely *is* published in two reports, and
deleting one of the rows would lose which reports carry it. `load_control` merges on
exact equality of `(easting, northing, elevation)` – an identity test with no tolerance –
keeps every id spelling and source report on the merged mark, and prints the count that
changed with the reason it changed. Rows that share a `point_id` but sit at different
positions raise rather than one winning silently (there are none in this file, checked).

---

## 4. What the returns-based assignment bought – measured

31 marks near Elba appear both in this module's 20 km run and in the scratchpad candidate
list that carried a centreline assignment. Same ties on both sides; only the grouping
differs.

```
marks in both the returns run and the centreline candidate list: 31
returns-assigned line EQUALS centreline-assigned line: 22 of 31 (71%)

disagreements:
            point_id      km_x  line  line_centreline    d_line_m  line_counts
L1O-6113 Wabasha RTK 14.697406   124              128 3676.087060  {'124': 96}
              L1O154 14.812691   153              149 3933.471196 {'153': 100}
              L5U153 14.880466   153              149 3995.145949  {'153': 80}
L5U-6112 Wabasha RTK 14.885250   124              128 3802.771782  {'124': 93}
               L1O69 15.988234   151              149 2331.164155  {'151': 96}
               L1O64 16.648176  1542              149 4490.083586 {'1542': 85}
              L5U121 17.249191   150              149  511.543313  {'150': 92}
L5U-1022 Olmsted VRS 17.672919   127              128 1380.009600  {'127': 91}
               L1O63 17.895325   156              149 6693.740778 {'156': 108}

agreement by distance band:
           n  agree
(0, 5]     2      2
(5, 10]    7      7
(10, 15]  12      8
(15, 25]  10      5
```

Note the mechanism, visible in the `d_line_m` column: the candidate list only held the
Elba network's own tracks, so a mark 6.7 km from every one of them was still handed the
nearest, and every mark whose returns say 150/151/153/156/1542 was labelled 149.

The effect on the answer, both groupings applied to the same 31 ties:

```
returns-based      k=19  F= 3.761  p=0.01174  df=(18, 12)  ICC=0.632  sd(line means)=  94.4  mean of line means=  -44.8  SE over lines=21.7
centreline-based   k=13  F= 5.023  p=0.001137  df=(12, 18)  ICC=0.643  sd(line means)=  66.7  mean of line means=   -7.8  SE over lines=18.5
```

**Read this carefully, because the ANOVA points the wrong way.** The centreline grouping
gives the *higher* F. That is not evidence it is the better grouping: collapsing 19 real
lines into 13 labels merges marks from lines that were never flown together, and an F
statistic cannot tell a real grouping from a tidy one. What matters is the 37.0 mm move
in the estimate, and that the 9 relabelled marks are, on the acquisition's own record,
under different aircraft passes than the label said.

### Where the centreline proxy actually fails, and where it does not

The mechanism is the **track set**, not the geometry. The candidate list held only
15 tracks – `[128, 131, 132, 133, 134, 135, 136, 137, 138, 140, 142, 143, 144, 145, 149]`
– so a mark 6.7 km from every one of them was still handed the nearest. **All 9
disagreements have a returns-line outside that set** (124, 127, 150, 151, 153, 156, 1542),
verified:

```
track set the candidate list knew: [128, 131, 132, 133, 134, 135, 136, 137, 138, 140, 142, 143, 144, 145, 149]
disagreements whose RETURNS-line is outside that track set: 9 of 9
```

The complementary half, measured here: the vendor's class-2 ground at a mark is almost
never shared between lines.

```
marks whose class-2 ground inside the report radius carries EXACTLY ONE point_source_id: 55 of 56
the one exception: [('L2T-6114 Wabasha RTK', {'124': 4, '125': 87})]
```

So within a complete local track set a centreline partition and the returns agree by
construction, and the returns-based assignment is not a correction to a sloppy proxy at
short range – it is what makes the proxy safe to abandon when the track set is partial,
which is every site the statewide workflow has not already fitted tracks for.
`analysis/GEN1_DATUM_MORE_MARKS.md` §1 and §3 give the mechanism behind that single-line
coverage (an explicit overlap class, cutting the ground class at a seam on the
perpendicular bisector) with the swath and spacing measurements behind it. **That is that
document's result and I have not re-derived it**; the 55-of-56 line above is mine.

---

## 5. Negative result: the siting screen does not help

From the 20 km `per_line` run, sd of the per-mark ties by radius-spread cut:

```
  marks measured                            56
  sd of the per-mark ties, mm             96.5
  sd of the line means, mm                94.9

  marks with radius_spread <= 15.0 mm       11   sd 105.6   mean  -48.8
  marks with radius_spread <= 25.0 mm       21   sd  95.4   mean  -80.3
  marks with radius_spread <= 50.0 mm       39   sd  92.6   mean  -45.2
  marks with radius_spread <= 100.0 mm      52   sd  92.9   mean  -51.5
  marks with radius_spread <= inf mm        56   sd  96.5   mean  -52.3
```

The tightest cut has the **largest** scatter. Whatever limits a single mark, it is not
radius instability, and the screen is not a filter worth applying. It stays in the module
as **reported statistics with no threshold**, which is also what
`ADDITIONAL_GROUND_CONTROL.md` §7.1 asked for – it proposed the criterion and said
explicitly that 17 measurements were not enough to set the cut. Measured on 56, the cut
is not there to be set.

This reproduces, on a different mark set and a different line assignment, the
`UNVERIFIED` figure carried in `FRAME_2026-08-26-PM.md` ("σ_site 83–97 mm at every ladder
cut"). My own numbers are 92.6–105.6 mm; the same conclusion, different values, and I
have not tried to reconcile the two sets.

---

## 6. Parameters I chose, each with its measured effect

Everything else is a repo default with a named source (`res_m = 5.0`,
`ground_quantile = 0.50`, `surface_order = 2`, `report_radius = 1.5·res`, the radius
ladder, `ground_source = vendor class 2`) or a caller argument with no default
(`radius_m`, `half_width_m`, `mode`, `--cover`, `tolerance_mm`).

| parameter | value | measured effect |
|---|---|---|
| **merge key** | exact `(easting, northing, elevation)` | 1 004 rows → **963 marks**, 41 merged in 39 groups. One merged mark falls in the 20 km measured set (`L2T-6126 Wabasha RTK`, −201.3 mm, line 130). Keeping the duplicate: line-first estimate **−55.55 → −55.79 mm (0.24 mm)**, mark-pooled mean **−52.31 → −54.92 mm (2.61 mm)**. The line-first average is nearly immune to a within-line duplicate; a pooled mean is not. |
| **`crop_half_width_m` DERIVED**, not chosen | `max(radius_ladder) = 5·res = 25 m` | the smallest square containing every fitting window. Against the 300 m used by the 2026-08-26 run, on five real gen1 marks: **`+0.000e+00 mm` on all five** (L1O101, L3B99, L5U171, L3B146, L4F145). A unit test asserts the same synthetically. This is the one parameter the earlier work flagged `MINE`; it is now derived and its effect is exactly zero for a vendor-class read. |
| **line-first averaging in `common_datum` too** | mean of line means, not mean of marks | elbaext frame, 14 marks: **+44.4 mm (line-first) vs +31.5 mm (marks pooled)**, a 12.9 mm difference, and the SE goes 22.1 → 30.5 mm. I chose line-first in both modes because the residual scatter is still organised by line after the constants are applied; the pooled number is returned as `mean_over_marks_mm` so the choice is visible. |
| **dominant line for a mark lit by two** | the `point_source_id` with the most returns | 1 of the 56 marks is mixed (`L2T-6114 Wabasha RTK`, 4 returns on 124 against 87 on 125). Assigning it to the minority line instead: **−55.55 → −55.32 mm**. Dropping it: **−55.54 mm**. Both per-line ties are returned in `per_line_tie_mm` and the mark is flagged, never dropped. |
| **ICC estimator** | one-way random effects with the unequal-size `n0` | reported for information beside F and p; nothing in the estimate depends on it. |

---

## 7. Verification

### Tests, and the three proven to bite

38 tests; the whole suite is **231 passed**. Three are regression tests in the strict
sense – each was broken in the committed source, the failure observed, and the source
restored with `git checkout --` (safe: it was committed first).

**(1) The line assignment comes from the returns.** Break: `measure_site` trusts
`site.nearest_feature` (the search's centreline answer) when it is set.

```
E       AssertionError: assert 136 == 137
E        +  where 136 = MarkMeasurement(...).line_id
1 failed, 37 deselected
```

**(2) The clustered SE does not shrink when marks are added to one line.** Break:
`combine_datum` pools the marks and takes `sd/sqrt(n)`.

```
E       assert 58.87840577551128 == 38.544964466378545 ± 1.0e-06
1 failed, 37 deselected
```

**(3) Duplicate rows are merged.** Break: `load_control` keys the groups per row, so
nothing ever merges.

```
E       AssertionError: assert 1004 == 963
E       AssertionError: assert 0 == 39
2 failed, 36 deselected
```

After each restore, `git status --short` on the module returned 0 lines, and the full
suite is back at 231 passed.

### What was run, on what

Three real runs at the Elba reference point (E 579 705.72, N 4 883 677.71), vendor class-2
ground, no CSF, no download, one at a time:

```
--radius-km 20 --mode per_line                                      56 marks / 27 lines
--radius-km 20 --mode common_datum --corrections …/elbaext/corrections_geoid.json
--radius-km 20 --mode common_datum --corrections …/elba_fulldensity/corrections.json
```

The banner of the first records: `1004 transcribed rows -> 963 physical marks; 41 rows
merged`; `no geoid conversion: all 98 marks and the lidar are on NAVD88(GEOID03);
checked, not assumed`; `27 of 65 tiles are on disk; 38 would have to be fetched and NONE
was`. Ledgers are in `.trust/runs/`.

### What I could NOT verify

* **`F = 8.63` over 43 marks – resolved, and it is not over 43 marks.** A one-way ANOVA on
  the scratchpad `screen_results.csv` using all 43 `ok` rows and all 15 line groups gives
  **`F=5.9404 p=3.233e-05 df=(14, 28)`**. Dropping the five singleton lines reproduces the
  published figure exactly:

  ```
  singleton lines dropped: k=10 n=38 F=8.6291 p=4.61e-06 df=(9, 28)
  ```

  (`analysis/GEN1_DATUM_MORE_MARKS.md` §5 reached the same reading independently; I
  re-derived the number rather than relaying it.) So `FRAME_2026-08-26-PM.md`'s "ANOVA
  over all 43 screened marks: F = 8.63" is **38 marks on 10 lines**, and the discarded
  singletons are what raises F from 5.94 to 8.63. `combine_datum` keeps singleton groups:
  a line with one mark contributes to the between-group sum of squares and nothing to the
  within-group one, which is what a singleton genuinely tells you.
* **`+53.6 ± 13.0 mm` on 18 open/urban marks.** My nearest equivalent – open and urban
  only, 20 km, elbaext frame – is **+53.7 ± 27.2 mm on 7 marks over 4 lines**. The value
  agrees to 0.1 mm and the SE does not, because it is the SE of a different statistic (SE
  over 4 lines, not over 18 marks) and because the returns-based assignment leaves only 7
  marks genuinely under lines 133–138 where the centreline assignment supplied 18. I have
  not reproduced the 18-mark set and it should not be assumed to be recoverable.
* **The 2008 control's datum.** The module asserts NAVD88(GEOID03) on both sides, but the
  claim itself comes from `lidar_semn2008.html` at dataset level; the validation reports
  print no datum and no geoid. The bundled CSV says so per row in its `verified` column.
  If that is wrong, every number here moves together.
* **Whether the swath constants transfer beyond the tile they were fitted on.** The
  per-line residuals (elbaext frame: 133 −105.5, 134 +21.1, 135 +125.6, 136 −16.1, 137
  −17.5, 138 −7.5 mm; RMS 68.3) mix network error with mark error, and three of the six
  lines carry one mark. This is the test `common_datum` mode exists to make possible; it
  is not a test that can be passed or failed at this n.
* **The absolute level itself.** `per_line` at −55.5 mm and `common_datum` at +44.4 mm are
  not in conflict – they are different quantities – but nothing here decides gen1's level
  at Elba, and this document does not claim to.

---

## 8. What is still open

* **The `per_line` and `common_datum` numbers have not been reconciled.** Doing it
  properly means asking what the 21 lines outside the elbaext frame are worth, which is a
  chain question (`chain.py`) and brings back the per-link error the 2008 route exists to
  avoid.
* **`discover_near_lines` has no measured user.** It is tested synthetically and unused by
  the driver; the statewide workflow will need it, and it needs real tracks to be exercised
  against.
* **38 tiles within 20 km of Elba hold control and are not on disk.** They are listed by
  every run. Fetching them is a decision with a cost and belongs to Andy, one tile at a
  time.
* **Line `1542`.** One mark (`L1O64`) sits under a `point_source_id` of 1542, four digits
  where every other line here is three. It is reported as its own group and nothing was
  done about it.
