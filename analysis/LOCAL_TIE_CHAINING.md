# The swath tie is not a constant of the flight line: measuring it where it is used

**Date:** 2026-08-26
**Code:** `src/lidar_diff_icp/localtie.py` (the module), `analysis/local_tie_chaining.py` (the run)
**Tests:** `tests/test_localtie.py` – 23 pass; the regression test is proven to bite three ways (§8)
**Run:**

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/local_tie_chaining.py --section all

Every number below is pasted from that command's output, ledger
`.trust/runs/20260826T195303-1380438.json`. The run is deterministic: re-running it after a
change that touched only column *labels* reproduced every digit of every table.

**Read, not redone:** `analysis/CROSS_LINE_FIT.md` (which established the same thing about the
across-track coefficient and is the direct ancestor of this work),
`analysis/ABSOLUTE_BASIS_ELBA.md` §2 (the 8.4 mm-per-link figure, re-derived here rather than
quoted), `analysis/SWATH_TIE_INTERCEPT.md` (the intercept tie). `coreg.py`, `pipeline.py` and
everything under `src/lidar_diff_icp/groundtruth/` are **imported, not modified**.

---

## The headline, in three sentences

**The vertical tie between two gen1 flight lines is local, the same way the across-track
coefficient is local.** On pair 136-137, tied in a 400 m-half-width window once per tile over
90 km of track, the tie runs **-57.1 to +28.8 mm, sd 28.3 mm** – while the same estimator, at
the same window size, tied at nine windows a few hundred metres apart inside **one** tile,
scatters by only **6.0 mm**. The variance ratio is **F = 22.18 on (5, 8) d.f., p = 0.0002**.

**`coreg`'s own uncertainty on that tie is 0.4 mm.** It is not wrong – it is the standard error
of a median over ~80 000 cells – it is simply not an estimate of the quantity anyone cares
about. The honest error bar is how far the number moves when the window moves or changes size:
on that pair, **28.3 mm between places and 29.5 mm between window sizes**, against the 0.4 mm
`coreg` prints beside it.

**At a real control mark this costs tens of millimetres.** Carrying a mark on line 136 into line
137's frame, elbaext's fitted constants supply **+8.6 mm everywhere**; measured at the mark, it
is **+12.4 mm** 3.8 km from Elba and **+78.9** and **+71.1 mm** at 28.6 and 29.6 km. Over the
twelve marks that need at least one link, local minus imported is **mean +20.8, sd 36.6,
RMS 40.7 mm**.

---

## 1. What was built

`src/lidar_diff_icp/localtie.py`. It does not reimplement the tie: every vertical offset it
returns comes out of `coreg.coregister_swaths` (or `coreg.align_swaths` for the network case),
unmodified, in one of its two tie modes. The only thing added is **where the points come from**.

| entry point | what it answers |
|---|---|
| `window_cloud(tiles, easting=, northing=, half_width_m=, shape=)` | the gen1 returns about a location, **with their scan angles** |
| `local_pair_tie(window, a, b, res_m=, tie=, exclude=)` → `LocalTie` | the tie of one pair **on that window** |
| `pair_tie_at(tiles, a, b, easting=, northing=, half_width_m=, …)` | crop and tie in one call |
| `window_ladder(…, half_widths_m=[…])` → `WindowLadder` | how the tie moves with window size; `.spread_mm`, `.sd_mm` |
| `local_network(tiles, lines, …, ref_line=)` → `LocalNetwork` | every line in the window solved into one line's frame |
| `nearest_overlap_point(tiles, a, b, easting=, northing=, …)` → `OverlapPoint` | the closest place to a mark where a pair *can* be tied, and how far away that is |
| `chain_local(tiles, easting=, northing=, source_line=, target_line=, …)` → `LocalChain` | **the headline API** |
| `compare_to_constants(chain, imported)` → `TieComparison` | local against a constant fitted somewhere else |
| `plan_path_local`, `inventory_from_cache` | route planning, delegated to `groundtruth.chain` |
| `TileCache` | per-tile arrays, in memory and on disk; `release()` because this is a shared laptop |

The question the brief asks, in code:

```python
ch = localtie.chain_local(tiles, easting=E, northing=N, source_line=L, target_line=137,
                          half_width_m=400.0, shape="square", res_m=2.0, tie="intercept",
                          exclude=(5, 6, 9), cache=cache,
                          ladder_half_widths_m=[100., 200., 400., 800., 1200.])
ch.dz_total_mm            # add this to the mark's gen1 elevation
ch.dz_sigma_window_m      # the error to quote (per-link ladder spreads, in quadrature)
ch.dz_sigma_formal_m      # coreg's own sigma – measured below to be ~30x too small
ch.max_solve_distance_m   # how far from the mark the farthest link had to be solved
```

**There is no default for `half_width_m`, `res_m`, `tie`, `shape` or `exclude`.** All five are
required keyword arguments. A default window size would have hidden the answer to the question
the module exists to ask; `tests/test_localtie.py::test_there_is_no_default_window_size_or_tie_mode`
holds that open by requiring `TypeError` when one is omitted.

---

## 2. Design decisions, and why

**Reuse `groundtruth/chain.py` for the graph and the route; do not reuse its point cloud.**
`chain.plan_path` already does the right thing – the along-swath (zero-link) case first, then
breadth-first over the *measured* overlap graph, minimising link count because each link adds
error – and it is called, not copied. `inventory_from_cache` fills `chain.SwathInventory` from
arrays this module has already read, so `chain.overlap_graph` and `chain.plan_path` run without
decompressing anything twice.

**The one deliberate fork, and its reason.** `chain._pair_cloud` builds its `PointCloud` with
`scan_angle=np.zeros(...)`. `tie="intercept"` reads that field; on an all-zero predictor
`coreg.across_track_tie` returns NaN and `coregister_swaths` silently falls back to the median
tie. A chain built on `chain.py`'s cache therefore *cannot* use the pipeline's own tie mode and
would not say so. `localtie` carries its own reader, which keeps every class and the scan angle,
read through `groundtruth.tie.scan_angle_deg` – which raises rather than returning zeros when
the dimension is absent. This is a fork of the *cache format*, not of the estimator.

**A degenerate window is reported, not hidden.** `coreg.nuth_kaab` abandons its fit when fewer
than 100 grid cells clear its 3° slope floor, and in that branch returns **`dz = 0.0` exactly**
with `n = 0` – not NaN, not the overlap median. On a whole tile that branch never fires; on a
small window over flat ground it fires easily, and a silent 0.0 mm tie is the worst failure
available to a module whose job is small windows. Every `LocalTie` therefore carries
`degenerate` and, beside it, `dz_overlap_median_m` – an independent read of the same overlap by
`swathdiff.swath_difference`. Nothing is dropped on that basis.

**The window is centred where a link *can* be solved, and the displacement is an output.** A
chain from a mark's line to a target line moves sideways across the flight direction, so the far
links cannot be solved at the mark. `nearest_overlap_point` finds the closest cell of the pair's
overlap and `LocalChain` reports `solve_distance_m` per link – the distance over which that
link's tie is being assumed constant. That is precisely the quantity that went unstated when
Elba's constants were applied 62.9 km away.

**`nearest_overlap_point` also reports the strip's own centre, because the nearest cell is an
edge.** On a north-south sidelap the cell nearest a point to its west sits on the strip's
**western edge**, and an edge window samples the across-track term one-sidedly. Both centrings
are run in §5 and both are reported.

---

## 3. The 8.4 mm per link, re-derived rather than quoted

`elba` and `elbaext` align the same flight lines from different extents and nothing else differs,
so their disagreement *is* the extent-dependence. Regauged onto line 135, which both contain:

| line | links | disagree_mm | per_link_mm |
|---|---|---|---|
| 136 | 1 | -8.0 | -8.0 |
| 137 | 2 | -9.8 | -6.9 |
| 138 | 3 | -17.4 | -10.0 |

    RMS disagreement           = 12.42 mm
    RMS per link               = 8.42 mm

**8.42 mm per link, verified.** `analysis/ABSOLUTE_BASIS_ELBA.md` §2 reports 12.4 mm RMS and
8.4 mm per link; this is the same quantity computed from the two saved products in this session.

---

## 4. Window size is a first-class output, and it dominates

Three pairs, tied at the Elba centroid, over a factor-of-12 ladder of window half-widths. `sig`
is `coreg`'s own 1σ; `spread` is max - min over the ladder.

| pair | tie mode | 100 m | 200 m | 400 m | 800 m | 1200 m | spread | `coreg` σ range |
|---|---|---|---|---|---|---|---|---|
| 135-136 | overlap_median | +7.2 | -6.5 | -19.2 | -25.1 | -24.6 | **32.2** | 0.4–4.2 |
| 135-136 | intercept | +10.1 | -21.0 | -19.3 | -22.3 | -22.5 | **32.6** | 0.4–4.2 |
| 136-137 | overlap_median | -35.3 | -27.0 | -9.6 | -5.8 | -9.1 | **29.5** | 0.4–1.6 |
| 136-137 | intercept | +1.3 | -12.5 | +3.9 | -1.5 | -7.1 | **16.4** | 0.4–1.6 |
| 137-138 | overlap_median | -6.9 | -46.7 | -33.0 | -25.2 | -21.1 | **39.8** | 0.5–3.8 |
| 137-138 | intercept | +71.7 | -66.8 | -17.8 | -22.7 | -17.8 | **138.5** | 0.5–3.8 |

(mm; the sign convention is `LocalTie`'s: add `dz` to the **src** line to reach the **ref**
line's frame.)

**Two things come out of this table.**

**(a) The formal σ is not the uncertainty of interest.** It runs 0.4 to 4.2 mm while the number
itself moves 16 to 140 mm across the ladder. This is the same shape of result as
`groundtruth/tie.py`'s radius pathology, and the same answer: report the ladder, and take the
spread over it as the headline uncertainty.

**(b) The intercept tie is unstable on small windows, and the module says why.** At half-widths
100 and 200 m, **every pair is flagged `extrap`** – the sampled across-track coordinate
`dtan = tan(scan_ref) - tan(scan_src)` does not reach zero, so the intercept is read outside the
data. On 137-138 at 100 m the window spans `dtan` +0.162 to +0.253 and the intercept is
**+71.7 mm** against -17.8 mm at 1200 m. That is not a tie, it is a lever arm. Where `dtan = 0`
*is* sampled (400 m and up), the intercept is the better-behaved of the two on 136-137
(spread 16.4 against 29.5) and is no worse on the others.

**This does not contradict `SWATH_TIE_INTERCEPT.md`; it bounds it.** The intercept tie removes
the dependence on *which part of the sidelap* an extent covers, and on a tile-sized extent it
does exactly that. It cannot do so on a window that does not contain the sidelap centre, and
`LocalTie.extrapolated` is how a caller finds out.

---

## 5. The tie along track: 90 km of one pair

One 400 m window per tile, centred on the pair's overlap strip. `d_elba_km` is signed
along-track distance from the elbaext centroid.

**Pair 136-137, column-64 tiles** (the clean case: 73 000–95 000 overlap cells per window, no
`extrap` flag anywhere):

| tile | d_elba_km | dz_med_mm | σ_med | dz_int_mm | dz_plain_mm | c (mm/tan) | ovl cells |
|---|---|---|---|---|---|---|---|
| 4342-21-64 | +28.0 | **-57.1** | 0.4 | -59.0 | -60.0 | +197 | 95 056 |
| 4342-28-64 | +3.7 | -16.3 | 0.6 | -19.2 | -5.0 | +247 | 80 177 |
| 4342-29-64 | +0.3 | -4.6 | 0.4 | +0.8 | +0.0 | +191 | 83 958 |
| 4342-30-64 | -3.2 | -9.8 | 0.6 | -12.5 | -20.0 | +295 | 84 747 |
| 5142-14-64 | -58.7 | **+28.8** | 0.3 | +27.6 | +25.0 | +96 | 75 896 |
| 5142-15-64 | -62.2 | +4.3 | 0.4 | +3.3 | +10.0 | +63 | 73 393 |

    overlap_median  n=6  mean    -9.1  sd   28.3  range   -57.1 to   +28.8  (mean formal sigma 0.4 mm)
    intercept       n=6  mean    -9.8  sd   29.0  range   -59.0 to   +27.6  (mean formal sigma 0.4 mm)

**Pair 135-136, column-63 tiles** – reported with a warning attached. What these tiles hold is
the *western sliver* of the 135-136 sidelap, not the sidelap: **28 422 to 37 582** overlap cells
per window against 73 393 to 95 056 for 136-137, and `extrap` is flagged in 8 of the 9 windows.
(A separate `point_source_id` scan of the tiles on disk agrees: line 136 contributes 250 711 to
379 115 returns to a column-63 tile, out of ~7 M.) The numbers are real measurements of an edge:

    overlap_median  n=9  mean    -1.9  sd   30.7  range   -44.4 to   +61.7  (mean formal sigma 1.0 mm)
    intercept       n=9  mean   -31.2  sd   36.7  range   -77.7 to   +29.9  (mean formal sigma 1.0 mm)

**None of the §6 conclusion rests on this pair.**

### The control that makes §5 mean something

A scatter of 28 mm over 90 km is only evidence of *position dependence* if it exceeds what the
same estimator produces between windows a few hundred metres apart. Same pair, same tie modes,
same window sizes, nine and seven windows inside **one** tile:

| half_width_m | tie_mode | n_windows | mean_mm | sd_mm | min_mm | max_mm | span_km |
|---|---|---|---|---|---|---|---|
| 400 | overlap_median | 9 | -12.0 | **6.0** | -20.9 | -4.4 | 2.4 |
| 400 | intercept | 9 | -13.7 | **10.7** | -26.0 | +1.1 | 2.4 |
| 800 | overlap_median | 7 | -9.3 | **2.5** | -12.7 | -5.8 | 1.8 |
| 800 | intercept | 7 | -7.9 | **5.1** | -17.0 | -1.9 | 1.8 |

Short-range scatter **shrinks with window size** (6.0 → 2.5 mm), which is what estimator noise
does. The long-range scatter does not go away.

| tie_mode | n_long | sd_long_mm | n_short | sd_short_mm | F | p | n_pairs | r(sep, abs diff) |
|---|---|---|---|---|---|---|---|---|
| overlap_median | 6 | 28.3 | 9 | 6.0 | **22.18** | **0.0002** | 15 | +0.564 |
| intercept | 6 | 29.0 | 9 | 10.7 | **7.34** | **0.0074** | 15 | +0.484 |

**The tie depends on where along the line it is measured, and the dependence is far larger than
the estimator's own repeatability.** The correlation between the separation of two samples and
the absolute difference of their ties is +0.564 over all 15 pairs – positive, as a spatially
structured quantity requires, but on 15 non-independent pairs it is a description, not a test.
The F ratio is the test.

**Sensitivity to the window centring** (strip centre against the cell nearest the tile centre;
both are printed in the run):

| pair | sd, strip-centred | sd, tile-centred |
|---|---|---|
| 136-137 | 28.3 | 26.3 |
| 135-136 | 30.7 | 31.2 |

The choice moves the answer by ~2 mm against a 28 mm effect. It is not load-bearing.

---

## 6. What it costs at a real control mark

Each of the 23 MnDNR-2008 control marks inside the lines 133-138 corridor, on the tiles on disk.
`local_mm` is what `chain_local` measures **at the mark**; `imported_mm` is the same constant
implied by `data/derived/elbaext/corrections_geoid.json` (`dz[line] - dz[137]`, so the gauge
cancels); `sig_window_mm` is the per-link window-ladder spread in quadrature.

|                  mark |  d_elba_km | line |                path | local_mm | σ_formal | σ_window | imported_mm | diff_mm |
|---|---|---|---|---|---|---|---|---|
|                 L2T51 |      -0.7 |  134 |     134-135-136-137 |    +43.4 |      1.0 |    121.9 |       +40.4 |    +3.0 |
|                L1O101 |      -1.5 |  137 |                 137 |     +0.0 |      0.0 |      0.0 |        +0.0 |    +0.0 |
|                L5U171 |      -2.1 |  137 |                 137 |     +0.0 |      0.0 |      0.0 |        +0.0 |    +0.0 |
|                L3B143 |      -2.4 |  133 | 133-134-135-136-137 |    -18.4 |      1.2 |    121.7 |       +18.4 |   -36.8 |
|                 L2T54 |      +3.8 |  136 |             136-137 |    +12.4 |      0.6 |     28.1 |        +8.6 |    +3.8 |
|                L2T98 |      +6.1 |  134 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
|                L5U172 |      +7.9 |  134 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
|  L1O-6123 Wabasha RTK |     +14.3 |  133 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
|  L2T-6100 Wabasha RTK |     +14.4 |  135 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
|  L1O-6196 Wabasha RTK |     +24.5 |  134 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
|  L3B-6190 Wabasha RTK |     +27.4 |  135 |         135-136-137 |    +18.3 |      1.1 |     58.1 |       +24.6 |    -6.3 |
|  L1O-6189 Wabasha RTK |     +27.7 |  135 |         135-136-137 |    +17.9 |      1.1 |     78.3 |       +24.6 |    -6.7 |
|  L1O-6182 Wabasha RTK |     +27.7 |  137 |                 137 |     +0.0 |      0.0 |      0.0 |        +0.0 |    +0.0 |
|  L2T-6183 Wabasha RTK |     +28.6 |  136 |             136-137 |    +78.9 |      0.7 |     21.6 |        +8.6 | **+70.3** |
|  L3B-6181 Wabasha RTK |     +29.5 |  134 |     134-135-136-137 |    +81.9 |      1.0 |     39.1 |       +40.4 |   +41.5 |
|  L5U-6188 Wabasha RTK |     +29.6 |  136 |             136-137 |    +71.1 |      0.7 |     23.9 |        +8.6 | **+62.5** |
|  L5U-6187 Wabasha RTK |     +30.8 |  135 |                   – |        – |        – |        – |           – | line 137 absent from this tile |
| L5U-2122 Fillmore VRS |     -60.4 |  137 |                 137 |     +0.0 |      0.0 |      0.0 |        +0.0 |    +0.0 |
| L4F-2123 Fillmore VRS |     -60.9 |  137 |                 137 |     +0.0 |      0.0 |      0.0 |        +0.0 |    +0.0 |
| L5U-2119 Fillmore VRS |     -61.5 |  138 |             138-137 |    +44.2 |      0.7 |     23.9 |        -4.2 |   +48.4 |
| L5U-2120 Fillmore VRS |     -61.7 |  138 |             138-137 |    +52.5 |      0.6 |     16.7 |        -4.2 |   +56.7 |
| L5U-2121 Fillmore VRS |     -62.5 |  138 |             138-137 |    +35.6 |      0.5 |     10.3 |        -4.2 |   +39.8 |
| L1O-2124 Fillmore VRS |     -62.6 |  133 |       133-10011-137 |     -8.4 |      0.6 |    117.6 |       +18.4 |   -26.8 |

    all marks with a route                               n=17  mean   +14.7  sd   31.9  RMS   34.2 mm
    of those, ALREADY on line 137 (zero links)           n= 5  mean    +0.0  sd    0.0  RMS    0.0 mm
    of those, needing at least one link                  n=12  mean   +20.8  sd   36.6  RMS   40.7 mm

**Four things to read here.**

**The zero-link rows are a check, not a result.** Five marks sit under line 137 itself;
`plan_path`'s along-swath-first rule returns a zero-link chain, both constants are 0 by
construction, and they agree exactly. That the machinery gets the trivial case trivially right
is worth one line and no more – and it is why the summary is stratified by link count rather
than pooled.

**The single-link rows carry the argument.** The imported constant for a 136 → 137 transfer is
`+8.6 mm` at every mark, because it is one number. Measured locally it is `+12.4` at 3.8 km and
`+78.9` / `+71.1` at 28.6 / 29.6 km. For 138 → 137 the imported constant is `-4.2 mm`
everywhere; measured at the three Fillmore marks 61–63 km south it is `+44.2`, `+52.5`, `+35.6`.
Note that this is **one** near mark against five far ones – the *pattern* rests on n = 1 at the
near end, and the evidence that the tie moves with distance is §5, not this column.

**The long chains say so themselves.** The four- and five-link paths at Elba carry
`σ_window` of 121.9 and 121.7 mm. That is the module reporting that a four-link chain of
locally-measured ties is not worth having, which is the same conclusion `CROSS_LINE_FIT.md` §9
reached about propagating coefficients down a chain ("numerically useless after two links").

**The route planner used a cross line without being told to.** `L1O-2124` was routed
`133-10011-137`: two links through flight line **10011**, a cross line in the 5142-14/15 tiles,
instead of four links along 134-135-136. That is `chain.plan_path` minimising link count over
the *measured* overlap graph, and it is the first time this project has chained through a cross
line. Its `σ_window` is 117.6 mm, so it is not yet a useful route – but it was found.

---

## 7. Every parameter I chose, and what it does

None of these has a library default; all are arguments of `analysis/local_tie_chaining.py` and
all are printed by the provenance banner with `src="MINE"`.

| parameter | value | measured effect |
|---|---|---|
| `half_width_m` | 400 | **large.** §4: the tie moves 16.4–138.5 mm across half-widths 100–1200 m. 400 m is one rung of that ladder, chosen so a window fits inside a 2.5 km tile at every along-track sample. Re-runnable at any rung with `--half-width-m` |
| `ladder_half_widths_m` | 100, 200, 400, 800, 1200 | it *is* §4. Nothing is excluded on its basis; the whole ladder prints, including the two rungs where the intercept extrapolates |
| `shape` | square | not measured. A disk gives a ragged grid edge and fewer cells per half-width; `--shape disk` switches it and the library supports both |
| `centre` | strip **and** tile | ±2 mm on a 28 mm effect (§5 sensitivity table). Both are run |
| `repeat_step_m` | 300 | sets the sampling of the short-range control in §5, not a cut. Neighbouring 400 m windows overlap at this spacing |
| `repeat_half_widths_m` | 400, 800 | the comparison that shows short-range scatter shrinking with window (6.0 → 2.5 mm) while long-range scatter does not |
| `mark_radius_m` | 25 | **it changes the answer for 3 of 12 marks.** Re-running `chain.covering_lines` at 10/25/50/100 m: `L4F-2123` reads line 138 at 10 m and 137 at 25 m and above; `L5U-2120` and `L5U-2121` read 138 except at 50 m. All three are Fillmore marks lying in the 137/138 sidelap, where "the mark's line" is genuinely ambiguous. It selects **which line** a mark is attributed to, never which marks are used – every mark is printed |
| `target_line` | 137 | the frame the marks are carried into; `--target-line` changes it |

`res_m = 2.0` and `exclude = (5, 6, 9)` are declared `repo`, not `MINE`: they are
`coreg.coregister_swaths`'s own values, and the class filter is the **vendor** classification, so
no CSF run happens anywhere in this work.

---

## 8. What was verified, and how

- **8.42 mm per link** – recomputed from `data/derived/elba/corrections.json` and
  `data/derived/elbaext/corrections.json` (§3), not quoted.
- **The module reproduces `coreg`'s own intercept.** `LocalTie.k_check_m` is computed by this
  module's mirror of `coregister_swaths`'s intercept branch;
  `test_the_across_track_diagnostic_reproduces_coreg_own_intercept` pins it to the `dz` that
  function returns, so the mirror cannot drift unnoticed.
- **The window is the size asked for**, in both shapes, and **carries the scan angle in
  degrees** – the field `chain.py` zeroes.
- **The tie follows a spatially varying offset**: a synthetic line whose offset ramps 0 → 100 mm
  from south to north is read as 15 mm in the south window and 85 mm in the north one.
- **The ladder grows when the offset varies across the window** and not when it is constant.
- **Regression test, proven to bite.** `test_a_flat_window_is_flagged_degenerate_and_the_median_read_survives`
  builds a flat overlap with a known 40 mm offset, where `nuth_kaab` abandons its fit and returns
  `dz = 0.0` with `n = 0`. Three breaks, each reverted, each confirmed against
  `git diff --stat` and a final byte-identical restore:

  | break | result |
  |---|---|
  | `degenerate=(int(c.n) == 0)` → `degenerate=False` | 1 failed (the flag) |
  | `dz_overlap_median_m=sd.median_offset` → `= c.dz` | 1 failed (the independent read) |
  | window mask → `np.ones(x.size, bool)` | 8 failed (the window itself) |

  Restored: `23 passed`.
- **The synthetic writer flies adjacent lines in opposite directions**, because otherwise `dtan`
  never reaches zero and `tie="intercept"` is an extrapolation – which is how the first draft of
  the test read a **2.24 m** tie for a 40 mm offset. That bug in my own fixture is what led to
  the `extrapolated` flag, and thence to the explanation of §4(b).

---

## 9. What I could NOT verify

**Whether a local constant is more ACCURATE than Elba's.** This work shows the two disagree by
tens of millimetres and that the disagreement is not estimator noise. It does **not** show which
is closer to the truth, because nothing here has an external vertical reference – every number
is a gen1-to-gen1 between-line difference. The experiment that would settle it is: read each
surveyed mark's lidar elevation with `groundtruth.tie.estimate_tie`, and compare the
between-line scatter of (lidar - surveyed) when the imported constants are applied against when
the local ones are. **I did not run it, deliberately.** That estimator is being re-run right now
by the agent that owns `src/lidar_diff_icp/groundtruth/` and `analysis/GEN1_DATUM_MORE_MARKS.md`,
and a second, separately-driven copy of the same mark ties would manufacture exactly the
method-artifact difference that makes the comparison worse than useless. It should be done once,
in that pipeline, with `localtie.chain_local` supplying the constants.

**Whether the along-track variation is navigation or terrain.** The within-tile control (§5)
bounds what the estimator does between windows a few hundred metres apart on the *same* terrain.
It does not bound what the estimator does when the terrain, land cover and canopy change between
tiles 30 km apart. A slowly-varying, cover-driven bias in the per-cell median ground surface
would produce a qualitatively similar signature. Distinguishing them needs a covariate the run
does not carry.

*(The reading I find natural, marked as hypothesis and not tested here: a between-line offset
that drifts along track is what a per-swath GPS/trajectory drift looks like, which is the
mechanism this project's MN-wide goal is aimed at. Nothing above depends on that being right.)*

**Pair 135-136 outside the sidelap.** In column-63 tiles that pair is sampled on the western
sliver of its overlap only (§5), with `extrap` set in 8 of 9 windows. Its numbers are honest
measurements of an edge and should not be read as that pair's tie.

**Six of the 23 marks could not be routed to line 137**, because the column-64 tile of their
tile row is not on disk. They are printed with the reason rather than dropped.

**The absolute level is untouched.** As with `align_swaths`, every quantity here is a difference
between lines; the group's absolute offset from any external datum is unchanged and still
requires ground control.

---

## 10. What I propose

1. **Do not treat a fitted per-swath constant as transferable along track.** `corrections*.json`
   should carry the extent it was fitted on (it already carries `bounds`) and any use of it
   outside that extent should quote a transfer uncertainty. The measured scale is the §5 sd,
   ~28 mm per pair, not the 0.4 mm `coreg` reports.
2. **When a mark is not on the study lines, measure its transfer at the mark**, with
   `chain_local(..., ladder_half_widths_m=...)`, and quote `dz_sigma_window_m`. Reject chains of
   more than two links: at Elba the four- and five-link paths carry 122 mm.
3. **Prefer marks that lie under a target line.** Five of the 23 need no transfer at all, and
   `plan_path`'s along-swath-first rule finds them for free.
4. **Feed §9's experiment into the mark-tie pipeline, not into a second one.** That is the one
   thing that would turn "the constants disagree" into "this constant is better".
5. **Nothing here asks for a change to `coreg.py` or `pipeline.py`.** The `tie="intercept"`
   recommendation of `SWATH_TIE_INTERCEPT.md` survives §4 with one addition worth carrying into
   any future use: on a window that does not sample `dtan = 0`, the intercept is a lever arm, and
   `LocalTie.extrapolated` is how to find out before believing it.
