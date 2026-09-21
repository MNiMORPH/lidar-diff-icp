# Which local gen1 tile breaks the heading/flight-order confound, and what it would cost

**Date:** 2026-09-22
**Status:** RECONNAISSANCE ONLY. Nothing was built, no CSF was run, no tile was fetched.
Every number below is pasted from a command run for this note; the command and its input
file are named beside it.

---

## 0. The confound, restated from the data

At `elbaext` the six gen1 passes alternate heading in perfect step with flight order.
Measured here, not recalled — in-grid returns from
`data/derived/elbaext/beam_offset_table.parquet` (columns `gps_time`,
`point_source_id`, `cell`, `in_grid`; northing recovered as
`Y0 + (cell//NX + 0.5)*RES`, the indexing stated at
`analysis/modules/vegetation_correction/beam_offset_table.py:35,74`), heading from OLS of
easting and northing on `gps_time`:

```
### elbaext  --  gen1 data/before/elbaext_gen1_merged.laz  --  lidar_semn2008 (GEOID03, collected 2008)
    grid 810x890 @ 5.0 m, bounds [575600.0, 4882200.0, 580050.0, 4886250.0], EPSG:26915; gps_time_type=0
    13,928,563 in-grid gen1 CSF-ground returns in data/derived/elbaext/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1    138    1544958        0.0   52.3       GPSwk Tue 17:55   359.0    N    4045     825   579630
      2    137    2966122     6297.7   51.9       GPSwk Tue 19:40   179.4    S    4045    1570   579032
      3    136    3124488     8288.5   53.8       GPSwk Tue 20:14   358.3    N    4045    1555   578060
      4    135    3042616     9319.0   51.8       GPSwk Tue 20:31   180.1    S    4045    1485   577110
      5    134    2692161    11246.3   51.6       GPSwk Tue 21:03   358.8    N    4045    1295   576250
      6    133     558218    12306.2   50.0       GPSwk Tue 21:20   179.4    S    4045     345   575775
```

Three things this table establishes that matter for everything below:

1. **One sortie.** All six passes fall inside 3 h 25 min of a single Tuesday afternoon
   and evening (pass starts at GPS week seconds 237,350 to 249,657, the last pass
   ending at 249,707; the LAS header carries
   `global_encoding.gps_time_type = 0`, GPS seconds of week, so the day of week and time
   of day are recoverable but the week is not).
2. **Heading is a function of pass number.** N, S, N, S, N, S with no exception.
3. **Spatial order follows too** — `E_mid` decreases monotonically, 579,630 to 575,775 m.

So heading, pass parity, time-within-sortie and across-track position all move together.
A heading-dependent instrument term and a time-varying GNSS term make the same prediction.

---

## 1. Method, and the three parameters I chose

Two readers were used, both light:

- **Registered sites** — the canonical `beam_offset_table.parquet` in each tile directory,
  four columns only (`gps_time`, `point_source_id`, `cell`, `in_grid`). Nothing else was
  read. Quantities are *in-grid*: a duration is time-in-tile, not the length of the whole
  flight line.
- **Raw LAZ** — `laspy.open(...).chunk_iterator(2_000_000)` keeping every 25th point
  (`STRIDE = 25`), one file at a time. 74 files, no cloud held whole. The whole sweep of
  `data/before/*.laz` took 1 m 37 s wall.

A **pass** is a maximal run of returns sharing one `point_source_id` with no internal
`gps_time` gap longer than **`GAP_S = 600 s`**. A **session** is a cluster of pass
mid-times with no neighbour gap over **`SESSION_S = 3600 s`**. Passes under
**`MIN_PTS = 50,000`** points are dropped as turn fragments.

**All three are MY CHOICE and are proposals, not established parameters.** What each is
worth:

- `GAP_S = 600 s`: a pass crosses a 3.5 km DNR quarter-quarter tile in ~45 s at ~80 m/s,
  so a 600 s hole *inside one tile* is a second visit. This is **not** the cross-tile
  test — memory `psid-is-not-a-flight-line` is explicit that across tiles a `gps_time`
  gap is usually one continuous line under tiles we do not hold. Used only within a tile.
- `MIN_PTS = 50,000`: over all 74 files the seven smallest passes carry 5,425–9,375 points
  and span 156–604 m along track (turns clipped at a tile corner); the next-smallest real
  pass carries 171,675 points over 868 m. The cut sits inside an 18x gap, so nothing in the
  data is near it. Without it, turn fragments masquerade as cross-tie lines — they did, in
  `4342-20-62` and `4342-20-63`, until the cut was applied.
- `SESSION_S = 3600 s`: unambiguous for the multi-day tiles (separations of 19.7 h to
  3.9 d). **Ambiguous at the 1–3 h scale** — it splits elbaext's own single sortie in two
  at its 6,298 s gap between passes 1 and 2, which is almost certainly the aircraft working
  an adjacent block rather than a new flight. Read "session" as "flight" only where the
  separation is many hours.

Three alias measures are reported per tile. The dominant along-track axis is the
point-count-weighted circular mean of headings mod 180; each pass is `+`, `-`, or a
cross-tie (more than 45 deg off the axis).

| measure | meaning | 1.000 means |
|---|---|---|
| `parity_alias` | fraction of `+/-` passes whose sign the best odd/even rule predicts | elbaext's exact confound: the pass NUMBER tells you the heading |
| `\|r(sgn,t)\|` | \|Pearson r\| of sign against pass mid `gps_time` | all of one heading early, all of the other late — heading aliased with TIME |
| `\|r(sgn,x)\|` | the same against across-track position | heading tells you WHERE the pass is |

A tile helps only if all three are well below 1, **or** it carries a cross-tie (no
alternation can produce one).

---

## 2. The headline: the parity metric is not the whole question

Ranking by `parity_alias` alone is misleading, and the survey showed why twice:

- `4342-18-11` (before_mnrv) scores a good `parity_alias = 0.667` — but its layout is
  `+++---`: three N passes, then three S passes 27,478 s later. Heading is now aliased with
  *time* instead (`|r(sgn,t)| = 0.993`). Same disease, different clothes.
- `4342-28-56` scores the worst possible `parity_alias = 1.000` — and yet it is one of the
  best tiles on disk, because its four passes are N,S on a Monday afternoon and N,S again
  on the Wednesday morning 39.6 h later. Both headings in both flights: `|r(sgn,t)| = 0.011`.

**The design that actually separates the two hypotheses is heading CROSSED WITH FLIGHT** —
both headings flown in each of two flights far apart in time. Then the heading effect is
the within-flight contrast and the flight (GNSS session) effect is the within-heading
contrast, and neither borrows from the other. Here is every local file that has it, or
that carries a cross-tie line:

```
tile                     dir                 npass nsess sessions_with_BOTH sep_between_them_s xtie      layout_by_session     tile_pts
4358-29-06.laz           before                  5     2                  2             335675    0                 -+-|+-    8,608,547
4342-28-56.laz           before                  4     2                  2             142605    0                  +-|+-    6,689,720
4342-21-08.laz           before_mnrv             6     2                  2              12304    0                +-|+--+   11,453,188
4342-16-08.laz           before_mnrv             6     2                  2              12296    0                +-|+--+   10,929,320
4342-28-61.laz           before                  4     2                  2              10167    0                  +-|+-   10,406,516
4342-29-61.laz           before                  4     2                  2              10163    0                  +-|+-   10,549,493
5158-01-05.laz           before                  7     4                  1                  0    2             X|+-+|-X|+    8,522,865
4342-28-64.laz           before                  5     2                  1                  0    1                 +X|-+-    8,302,172
4358-28-01.laz           before                  5     1                  1                  0    1                  -+-+X    8,674,160
5142-01-60.laz           before                  5     2                  1                  0    1                 +-+-|X    8,441,104
5142-14-64.laz           before                  5     3                  1                  0    1                +|-+-|X    8,254,787
5142-15-63.laz           before                  5     2                  1                  0    1                 +-+-|X    7,774,272
5142-15-64.laz           before                  5     3                  1                  0    1                +|-+-|X    8,187,305

(tiles with fewer than two crossed sessions AND no cross-tie are omitted; 60 of 73 files)
```

Reading the layout column: `-+-|+-` means a first flight with S, N, S and a second flight
with N, S. `X` is a cross-tie.

---

## 3. The single best candidate: `data/before/4358-29-06.laz`

Winona county, **`lidar_semn2008` — the same acquisition as elbaext** (county resolved with
`lidar_diff_icp.tiles.county_for_lonlat` from the tile centroid in
`data/mn_tile_centroids.csv`; project from `lidar_diff_icp.acquisitions.for_county`).
8,608,547 points. Five passes, all full-length, both headings flown in **both** flights,
the flights **3 d 21.2 h apart** (mean pass time per flight, 335,675 s):

```
data/before/4358-29-06.laz   8,608,547 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
     153   2560600      64668.808    Sun 17:57:48   45.1   180.8    S     3504   594356   4884693
     152   2869325      66573.557    Sun 18:29:33   47.6   357.5    N     3509   593452   4884680
    1542    605075      67046.183    Sun 18:37:26   42.3   179.5    S     3472   594849   4884700
     151   1273175     401498.860    Thu 15:31:38   42.5   358.3    N     3480   592888   4884671
    1512   1300375     402044.823    Thu 15:40:44   44.1   179.5    S     3474   592898   4884672
```

- **Flight A**, Sunday 17:57–18:38: S (153), N (152), S (1542).
- **Flight B**, Thursday 15:31–15:41: N (151), S (1512).
- Same-heading separation **337,376 s = 3.905 d** (S 153 on Sunday against S 1512 on
  Thursday); opposite-heading separation **470 s** inside flight A (152 against 1542) and
  **546 s** inside flight B.
- `parity_alias = 1.000`, but that is the metric's blind spot here: the alternation happens
  to continue across the flight boundary. The measures that matter are
  `|r(sgn,t)| = 0.167` and a crossed 2x2.
- The two flights are at different times of day (evening vs mid-afternoon), which also
  decorrelates anything diurnal.
- GPS week seconds give day-of-week and time-of-day but not the week, so "3.905 d" assumes
  both flights are in the same GPS week. If they are not, they are *further* apart, not
  closer — the conclusion that they are separate flights holds either way.

**Runner-up, same acquisition: `data/before/4342-28-56.laz`** (Wabasha county,
`lidar_semn2008`, 6,689,720 points), 39.6 h apart and balanced 2 N / 2 S:

```
data/before/4342-28-56.laz   6,689,720 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
     116   2660675     145359.693    Mon 16:22:39   45.0   358.2    N     3498   558779   4887751
     115   2003450     146508.355    Mon 16:41:48   44.5   180.5    S     3486   557996   4887744
    1142    147500     287601.735    Wed 07:53:21   47.3   359.5    N     3469   557519   4887740
     117   1878050     289476.167    Wed 08:24:36   43.9   180.1    S     3478   559489   4887757
```

`|r(sgn,t)| = 0.011` — as close to orthogonal as heading and time get on disk. Its weakness
is pass 1142, which contributes only 147,500 points over a 95 m ACROSS-track span (a clipped
edge of the tile).

**Third, and the cheapest by a wide margin: `data/before/4342-28-64.laz`.** This is the
tile immediately NORTH of elbaext, flown on the *same sortie* by the *same four lines*
(138, 137, 136, 135 at Tue 17:56, 19:40, 20:14, 20:30 — compare elbaext's table in Sec 0),
plus **a due-west cross line, psid 10010, flown at Tue 18:12:00, between elbaext's pass 1
and pass 2**:

```
data/before/4342-28-64.laz   8,302,172 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
     138   1531950     237403.583    Tue 17:56:43   44.9   359.0    N     3474   579583   4887969
   10010    747225     238320.017    Tue 18:12:00   24.4   270.0    W     1959   578980   4889299
     137   2841125     243608.301    Tue 19:40:08   44.6   180.9    S     3508   578983   4887962
     136   2632100     245692.223    Tue 20:14:52   46.0   358.6    N     3485   578110   4887952
     135    549775     246629.321    Tue 20:30:29   42.9   179.7    S     3472   577624   4887946
```

A 270 deg heading cannot be produced by any odd/even rule, and it sits 916 s into the same
sortie whose GNSS drift is the competing explanation. Six further semn2008 tiles carry the
same family of cross lines (psid 10008, 10009, 10011, 10012): `5158-01-05` (two of them,
in two different flights), `4358-28-01`, `5142-01-60`, `5142-14-64`, `5142-15-63`,
`5142-15-64`.

### 3.1 Acquisition attribution for every tile named above

Resolved the way the package itself resolves it: centroid from
`data/mn_tile_centroids.csv` -> lon/lat -> `lidar_diff_icp.tiles.county_for_lonlat`
(FCC Census Area lookup, `verify=False` so no county listing is fetched) ->
`lidar_diff_icp.acquisitions.for_county`. No geoid is ever defaulted.

```
4358-29-06  centroid (593769.2, 4884684.2) EPSG:26915 -> (-91.82828, 44.10934)  county=winona  project=lidar_semn2008  geoid=GEOID03
4342-28-56  centroid (558725.0, 4887750.1) EPSG:26915 -> (-92.26580, 44.14059)  county=wabasha  project=lidar_semn2008  geoid=GEOID03
4342-28-64  centroid (578721.9, 4887958.8) EPSG:26915 -> (-92.01579, 44.14059)  county=winona  project=lidar_semn2008  geoid=GEOID03
5158-01-05  centroid (591460.6, 4870765.4) EPSG:26915 -> (-91.85953, 43.98434)  county=winona  project=lidar_semn2008  geoid=GEOID03
4342-28-61  centroid (571223.1, 4887873.4) EPSG:26915 -> (-92.10954, 44.14059)  county=wabasha  project=lidar_semn2008  geoid=GEOID03
4342-29-61  centroid (571260.7, 4884402.4) EPSG:26915 -> (-92.10954, 44.10934)  county=wabasha  project=lidar_semn2008  geoid=GEOID03
4358-28-07  centroid (596219.3, 4888191.3) EPSG:26915 -> (-91.79703, 44.14059)  county=winona  project=lidar_semn2008  geoid=GEOID03
4358-32-06  centroid (593917.4, 4874271.4) EPSG:26915 -> (-91.82828, 44.01559)  county=winona  project=lidar_semn2008  geoid=GEOID03
```

**Every candidate is `lidar_semn2008`, GEOID03 — the same survey as elbaext.** Two tiles
elsewhere in `data/before` (`4358-26-05`, `4358-27-07`) have centroids that the lookup
places in Wisconsin, across the Mississippi; their MN portions are almost certainly Winona
but the package's own authority declines to say so, and I leave it UNRESOLVED rather than
assert it. `4358-26-05` is in any case only 2,190,404 points and half turn fragments.

---

## 4. What already exists, and therefore what it would cost

### 4.1 Analysis inputs present today

Presence check over every directory under `data/derived/` (`os.path.exists`):

```
tile_dir                                   stable  z_after_diffe  curv_laplacia  floodplain_ma  beam_offset_t  gen1_csf_angl        z_after            dod          slope  csf_cache
battlecreek                                    --             --            YES            YES            YES            YES            YES            YES            YES  YES
carlton                                        --             --            YES            YES            YES            YES            YES            YES            YES  YES
cook                                           --             --            YES            YES            YES            YES            YES            YES            YES  YES
elba                                           --             --            YES            YES            YES            YES            YES            YES            YES  YES
elba_fulldensity                               --             --            YES            YES            YES            YES            YES            YES            YES  --
elba_refdatum                                  --             --             --             --             --             --            YES            YES             --  --
elbaext                                        --             --            YES            YES            YES            YES            YES            YES            YES  YES
elbaext_delong                                YES            YES            YES            YES            YES            YES            YES            YES            YES  --
elbaext_delong_boresight                      YES            YES             --             --             --             --            YES            YES             --  --
elbaext_delong_fpsupport                      YES            YES             --             --             --             --            YES            YES             --  --
elbaext_independent                           YES            YES            YES            YES            YES            YES            YES            YES            YES  --
mnrv                                           --             --            YES            YES            YES            YES            YES            YES            YES  YES
whitewater                                     --             --            YES            YES            YES            YES            YES            YES            YES  YES
```

Note for anyone reaching for `stable.npy` / `z_after_differenced.npy`: they exist **only**
under the `elbaext_delong*` and `elbaext_independent` product directories, not under the
six registered sites. The `elbaext` directory itself carries `stable_geoid.npy`, which is
a different product.

**No candidate tile named in Sec 3 has a tile directory at all.** The one thing any of them
already has is a CSF ground cache:

```
data/csf_cache/4342-28-64.las            pts=      6941881 fmt=7 x[577455,579992] y[4886209,4889709]
data/csf_cache/elbaext.las               pts=     15581591 fmt=7 x[575450,580043] y[4882050,4886400]
```

and `4342-28-64.las` is class-2-ground-only with the cross line intact — counted directly:

```
classification counts: {2: 6941881}
GROUND (class 2) by point_source_id: {135: 444915, 136: 2143062, 137: 2382403, 138: 1331104, 10010: 640397}
```

### 4.2 The cost that dominates everything: there is no gen2

Header extents of every gen2 cloud on disk (`laspy`, `data/after*/*.la[sz]`):

```
data/after/3dep2021_fulldensity.laz                  pts=   182923322 x[577493,580035] y[4882738,4886238]
data/after/3dep_4358_fulldensity.laz                 pts=   148050625 x[584860,587400] y[4893245,4896745]
data/after/elbaext_3dep_fulldensity.laz              pts=   415080034 x[575450,580200] y[4882050,4886400]
data/after_battlecreek/battlecreek_3dep.laz          pts=     8437357 x[498750,499367] y[4975136,4976004]
data/after_battlecreek/gen2_4tile.laz                pts=    41184507 x[498135,499365] y[4975136,4976876]
data/after_carlton/carlton_3dep.laz                  pts=   199332626 x[547805,550225] y[5163676,5167168]
data/after_mnrv/mnrv_3dep2021.laz                    pts=   460784964 x[420190,422720] y[4903565,4907065]
data/after_ne/ne_3dep_fulldensity.laz                pts=   240897623 x[709531,711986] y[5323589,5327144]
```

Every gen2 extent is exactly a registered site's analysis grid (plus the elbaext and
battlecreek blocks). **None covers any candidate in Sec 3.** The best case is
`4342-28-64` (x 577,455–579,992, y 4,886,209–4,889,709) against the elbaext gen2
(y up to 4,886,400) — a 191 m sliver, and the cross line 10010 sits at N_mid 4,889,299,
2.9–3.5 km north of it, so the overlap is zero where it matters.

So the `beam_offset_table` quantity `d_mm_corr` — the per-return slope-normal residual of a
gen1 CSF ground return to the gen2 surface — **cannot be formed on any candidate tile**
without fetching a 3DEP cloud of the order of 10^8 points. That is outside this brief.

### 4.3 The route that needs no download

This has been faced before in this repo and solved, and the record is on disk.
`analysis/crossline_fit.py` used exactly `4342-28-64` and its line 10010 to break a
*different* degeneracy (the per-line across-track coefficient), and its docstring states
the substitution verbatim:

> No gen2 exists over the cross-line tile, which sits 2.6-3.5 km NORTH of the elbaext grid,
> so `d_mm_corr` cannot be formed there. This run therefore uses the identical estimator
> with the reference surface built from the POOLED gen1 ground of the same tile instead of
> gen2 [...] In a between-line difference a COMMON reference cancels [...] but that is an
> argument, so Sec 2 MEASURES it: both estimators are run on elbaext's N-S pairs and the
> coefficients compared.

`analysis/CROSS_LINE_FIT.md` reports that pass worked: line 10010 carries **640,397** ground
returns at **0.45 pts/m²**, inside the 0.46–0.52 of the N-S lines it crosses, and on its
three pairs `SE(q)/SE(p)` is **1.02–1.12** against 4.95–7.11 on the N-S pairs of the same
tile. `analysis/SWATH_DEGENERACY_BREAKING.md` is the survey that found it.

The gen1-internal machinery is already in the package and takes no gen2:
`lidar_diff_icp.swathdiff` (per-cell median bare-earth per swath, differenced on the shared
grid) and `coreg.align_swaths` (free-network least squares over swath pairs). Both key on
`point_source_id`, which every candidate LAZ carries.

### 4.4 Cost, per candidate

| candidate | acquisition | gen1 on disk | CSF cache | gen2 | derived products | what is missing |
|---|---|---|---|---|---|---|
| `4358-29-06` | `lidar_semn2008` (same as elbaext) | YES, 8,608,547 pts | **no** | **no, and none nearby** | none | a CSF run; then a gen1-internal estimator only |
| `4342-28-56` | `lidar_semn2008` | YES, 6,689,720 pts | **no** | **no** | none | same |
| `4342-28-64` | `lidar_semn2008`, same sortie as elbaext | YES, 8,302,172 pts | **YES**, 6,941,881 class-2 pts, 10010 intact | 191 m sliver only, useless | none | nothing, for the gen1-internal route — `analysis/crossline_fit.py` already runs on it |
| `5158-01-05` | `lidar_semn2008` | YES, 8,522,865 pts | **no** | **no** | none | a CSF run |
| `4342-28-61`, `4342-29-61` | `lidar_semn2008` | YES, 10.4 / 10.5 M pts | **no** | **no** | none | a CSF run; and their two flights are only 2.8 h apart, inside the ambiguous window of Sec 1 |
| battlecreek 4-tile block | `lidar_metro2011`, **different survey** | YES, 4 quarter-tiles | site tile only (`battlecreek.las`) | **YES**, `gen2_4tile.laz`, 41,184,507 pts, x[498135,499365] y[4975136,4976876] | full set for the site quarter only | a CSF run over the block; the merged gen1 is unusable (see 4.5) |

I have **no measured runtime** for `csf_tiled.py` on a single quarter-quad tile in this
session. Project memory records ~460 s for the elba CSF; treat that as the order, not a
measurement.

### 4.5 One defect worth knowing about

`data/before_battlecreek/gen1_4tile.laz` (24,642,050 points, written by laspy 2.5.3) has
**lost `point_source_id` and `gps_time` entirely**. Read directly:

```
pts 24642050 fmt 1 gps_type 0 gen laspy 2.5.3 sysid OTHER
n 24642050 psids [0] gps_time range 0.0 0.0
```

Any per-pass work on the battlecreek block must go to the four quarter-tile LAZs, which
retain both. (`data/before/elbaext_gen1_merged.laz`, also laspy-written, did keep them.
The untracked `scripts/merge_gen1_tiles.py` copies every dimension generically and would
not have this problem — stated as a fact about the tooling, not as a recommendation.)

---

## 5. The six registered sites, in-grid

Read from each site's `beam_offset_table.parquet` (four columns). The `clock` column
respects each gen1 header's `global_encoding.gps_time_type`: 0 is GPS week seconds
(day-of-week and time-of-day recoverable, week not), 1 is Adjusted Standard GPS Time,
converted to UTC with the GPS epoch 1980-01-06 and the leap-second offset for that
survey's year. The carlton conversion lands on **2012-11-02**, which matches the date this
project already holds for that tile — an independent check that the conversion is right.

```
### elba  --  gen1 data/before/4342-29-64.laz  --  lidar_semn2008 (GEOID03, collected 2008)
    grid 700x508 @ 5.0 m, bounds [577492.8, 4882737.6, 580032.8, 4886237.6], EPSG:26915; gps_time_type=0
    6,771,479 in-grid gen1 CSF-ground returns in data/derived/elba/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1    138    1303950        0.0   44.7       GPSwk Tue 17:55   358.8    N    3475     815   579623
      2    137    2544119     6290.6   45.0       GPSwk Tue 19:40   180.1    S    3490    1565   579033
      3    136    2424563     8288.0   46.3       GPSwk Tue 20:14   358.4    N    3485    1330   578160
      4    135     498847     9312.2   43.3       GPSwk Tue 20:31   179.6    S    3470     335   577663

### whitewater  --  gen1 data/before/4358-26-03.laz  --  lidar_semn2008 (GEOID03, collected 2008)
    grid 700x508 @ 5.0 m, bounds [584860.0, 4893245.0, 587400.0, 4896745.0], EPSG:26915; gps_time_type=0
    5,869,244 in-grid gen1 CSF-ground returns in data/derived/whitewater/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1    146     704373        0.0   44.8       GPSwk Tue 14:46     1.3    N    3475     630   587082
      2    145    2029353      599.5   44.5       GPSwk Tue 14:56   177.4    S    3480    1540   586582
      3    144    2290705     2813.4   44.1       GPSwk Tue 15:33   358.9    N    3485    1520   585628
      4    143     844813     3431.2   43.5       GPSwk Tue 15:43   179.1    S    3450     580   585152

### mnrv  --  gen1 data/before_mnrv/4342-23-01.laz  --  lidar_swmn2010 (GEOID03, collected 2010)
    grid 700x506 @ 5.0 m, bounds [420190.0, 4903565.0, 422720.0, 4907065.0], EPSG:26915; gps_time_type=0
    8,413,627 in-grid gen1 CSF-ground returns in data/derived/mnrv/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1   6191     131054        0.0   43.1       GPSwk Fri 15:10   180.6    S    3470     255   422590
      2   6201    1110057      460.5   43.0       GPSwk Fri 15:18     0.3    N    3475     785   422325
      3   6211    2179182     1889.0   46.1       GPSwk Fri 15:41   181.5    S    3485    1525   421955
      4   6221    2054611     2360.2   43.5       GPSwk Fri 15:49     1.3    N    3485    1525   421425
      5   6231    2044886     3792.9   46.2       GPSwk Fri 16:13   179.2    S    3480    1270   420828
      6   6241     892758     4427.8   41.7       GPSwk Fri 16:24     1.9    N    3475     660   420522
      7   6251       1079     5902.9    1.9       GPSwk Fri 16:48   179.9    S     145      10   420198

### cook  --  gen1 data/before_ne/1158-31-59.laz  --  lidar_arrowhead2011 (GEOID09, collected 2011)
    grid 711x491 @ 5.0 m, bounds [709531.0, 5323589.0, 711986.0, 5327144.0], EPSG:26915; gps_time_type=1
    4,244,559 in-grid gen1 CSF-ground returns in data/derived/cook/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1      8     600931        0.0   30.4  2011-05-27 15:00 UTC    86.1    E     990    2355   710711
      2      9    1264559      564.7   32.2  2011-05-27 15:09 UTC   263.2    W    1635    2380   710738
      3     10    1397106     1806.0   29.5  2011-05-27 15:30 UTC    93.9    E    1590    2380   710778
      4     11     981963     2476.5   32.3  2011-05-27 15:41 UTC   265.9    W    1000    2360   710804

### carlton  --  gen1 data/before_carlton/2742-12-53.laz  --  lidar_duluth2012 (GEOID09, collected 2012)
    grid 698x484 @ 5.0 m, bounds [547805.0, 5163676.0, 550225.0, 5167166.0], EPSG:26915; gps_time_type=1
    8,152,421 in-grid gen1 CSF-ground returns in data/derived/carlton/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1     90     161842        0.0   12.5  2012-11-02 18:19 UTC   318.4    N     655     625   548140
      2     89    1536130     1608.3   40.1  2012-11-02 18:46 UTC   135.4    S    2240    2105   548865
      3     88    3331929     2017.9   55.4  2012-11-02 18:53 UTC   326.1    N    3465    2415   549015
      4     87    2547411     3694.1   49.7  2012-11-02 19:21 UTC   140.5    S    3090    2415   549015
      5     86     575109     4065.3   25.7  2012-11-02 19:27 UTC   317.1    N    1485    1385   549515

### battlecreek  --  gen1 data/before_battlecreek/4342-03-32_b_a.laz  --  lidar_metro2011 (GEOID09, collected 2011)
    grid 174x123 @ 5.0 m, bounds [498750.0, 4975136.0, 499365.0, 4976006.0], EPSG:26915; gps_time_type=0
    4,158,192 in-grid gen1 CSF-ground returns in data/derived/battlecreek/beam_offset_table.parquet
    ord   psid  n_in_grid     t-t0_s  dur_s                 clock     hdg  dir  N_span  E_span    E_mid
      1   1012     696388        0.0   16.5       GPSwk Sun 09:45   177.8    S     865     295   499215
      2   1013    1498375      768.8   16.1       GPSwk Sun 09:58   358.0    N     865     585   499070
      3   1014    1375964     1457.3   16.9       GPSwk Sun 10:10   179.1    S     865     550   499028
      4   1101     587376    72432.0   15.8       GPSwk Mon 05:53   359.4    N     865     265   498885
```

Every one of the six is confounded the elbaext way: heading alternates with pass order in a
single flight. Three remarks:

- **cook** flies east-west (headings 86.1, 263.2, 93.9, 265.9), not north-south. Its
  heading axis is orthogonal to elbaext's, which is informative in its own right — but it
  is `lidar_arrowhead2011`, a different aircraft and GNSS session, and it alternates just
  as perfectly within itself.
- **carlton** flies diagonally (headings 318.4, 135.4, 326.1, 140.5, 317.1) with markedly
  unequal in-tile durations (12.5 to 55.4 s) and 3 NW against 2 SE — but all five passes
  are in one 68-minute flight, so heading is still a function of pass number.
- **battlecreek** is the only registered site whose analysis grid sees two flights
  (Sun 09:45–10:10 and Mon 05:53; in-grid first-return-to-first-return 72,432.0 s =
  20.13 h, and 71,663.2 s = 19.91 h between the two NORTH passes 1013 and 1101). Its grid
  clips only ONE pass
  from the Monday flight, so at the site scale heading is still aliased. At the four-tile
  BLOCK scale it is not — see Sec 6.

---

## 6. The battlecreek block: the one place both epochs already cover a crossed design

The four quarter-tiles, read from the raw LAZ:

```
data/before_battlecreek/4342-02-32_c_d.laz   6,084,356 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
    1012    970175      35137.591    Sun 09:45:37   16.2   179.9    S      868   499223   4976438
    1013   2130600      35933.073    Sun 09:58:53   16.3     2.1    N      885   499093   4976438
    1014   2077025      36594.336    Sun 10:09:54   17.0   179.9    S      868   499030   4976438
    1101    906575     107596.229    Mon 05:53:16   16.1     1.1    N      872   498877   4976438

data/before_battlecreek/4342-02-32_d_c.laz   5,981,571 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
    1014    631900      36593.956    Sun 10:09:53   16.9   179.2    S      870   498640   4976438
    1101   1898050     107596.489    Mon 05:53:16   16.0   358.6    N      878   498502   4976438
    1102   2116000     108289.846    Mon 06:04:49   15.3   179.3    S      875   498441   4976438
    1103   1163950     109000.130    Mon 06:16:40   16.3     0.9    N      872   498300   4976438
    1104    171675     109679.141    Mon 06:27:59   16.0   179.2    S      868   498181   4976438

data/before_battlecreek/4342-03-32_a_b.laz   6,367,451 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
    1014    766900      36608.050    Sun 10:10:08   16.4   180.7    S      870   498622   4975570
    1101   2058175     107583.248    Mon 05:53:03   15.8     0.6    N      872   498508   4975570
    1102   2142725     108302.506    Mon 06:05:02   15.5   178.5    S      882   498442   4975570
    1103   1161350     108986.470    Mon 06:16:26   16.5   359.1    N      872   498303   4975570
    1104    238325     109692.419    Mon 06:28:12   16.2   179.0    S      868   498189   4975570

data/before_battlecreek/4342-03-32_b_a.laz   6,208,672 pts   (gps_time_type=0 -> GPS seconds of week)
    psid    ~n_pts     gps_time_s           GPSwk  dur_s     hdg  dir  along_m    E_mid     N_mid
    1012    973875      35150.973    Sun 09:45:50   16.5   179.0    S      872   499218   4975570
    1013   2164775      35919.735    Sun 09:58:39   16.1     0.4    N      871   499073   4975570
    1014   2048100      36608.290    Sun 10:10:08   16.9   178.0    S      885   499024   4975570
    1101   1021850     107582.948    Mon 05:53:02   15.8   357.7    N      876   498884   4975570
```

Pooled over the block: flight A on Sunday 09:45–10:10 flies S (1012), N (1013), S (1014);
flight B on Monday 05:53–06:28 flies N (1101), S (1102), N (1103), S (1104). **Both
headings in both flights, 20.2 h apart** (mean pass time per flight, 72,753.0 s).
Extents:

```
data/before_battlecreek/4342-02-32_c_d.laz       pts=   6084356 x[498750,499367] y[4976004,4976872]
data/before_battlecreek/4342-02-32_d_c.laz       pts=   5981571 x[498134,498750] y[4976004,4976872]
data/before_battlecreek/4342-03-32_a_b.laz       pts=   6367451 x[498133,498750] y[4975136,4976004]
data/before_battlecreek/4342-03-32_b_a.laz       pts=   6208672 x[498750,499367] y[4975136,4976004]
data/after_battlecreek/gen2_4tile.laz            pts=  41184507 x[498135,499365] y[4975136,4976876]
```

`gen2_4tile.laz` covers the whole block. So this is the only crossed heading x flight
design on disk for which **both epochs are already present** and no fetch is needed.

Its cost: it is `lidar_metro2011` on GEOID09, **a different survey from elbaext**
(different aircraft, different GNSS sessions, different year). It is an independent
replicate of the question, not a measurement that can be pooled with elbaext as if the
same instrument produced it.

---

## 7. Full ranking, least confounded first

73 files (the 74th, `gen1_4tile.laz`, has no usable pass structure — Sec 4.5). `WORST` is
the largest of the three alias measures; `xtie` counts cross-tie passes, which break the
alias regardless of what the other columns say. `max_same_dt_s` is the longest interval
between two passes of the SAME heading — the lever for separating a time-varying error.

```
rank tile                     dir                  np  n+  n- xtie     time_seq parity_alias |r(sgn,t)| |r(sgn,x)|  WORST max_same_dt_s    span_s
   1 4342-21-05.laz           before_mnrv           7   3   4    0      -+--+-+        0.571      0.223      0.294  0.571          3479      4073
   2 4342-20-09.laz           before_mnrv           6   5   1    0       +++++-        0.667      0.637      0.636  0.667         96264     97381
   3 4342-21-08.laz           before_mnrv           6   3   3    0       +-+--+        0.667      0.019      0.123  0.667         15861     15861
   4 4342-16-08.laz           before_mnrv           6   3   3    0       +-+--+        0.667      0.016      0.137  0.667         15848     15848
   5 4342-24-05.laz           before_mnrv           8   4   4    0     -++--+-+        0.750      0.032      0.208  0.750          3481      3813
   6 4342-16-14.laz           before_mnrv           6   5   1    0       +++++-        0.667      0.758      0.638  0.758          9461     14052
   7 4342-21-14.laz           before_mnrv           6   5   1    0       +++++-        0.667      0.772      0.634  0.772          9459     14470
   8 4358-28-07.laz           before                4   1   3    0         --+-        0.750      0.576      0.778  0.778        183117    183117
   9 4342-25-14.laz           before_mnrv           6   5   1    0       +++++-        0.667      0.782      0.635  0.782          9464     14808
  10 4358-32-06.laz           before                5   3   2    0        -+-++        0.800      0.625      0.553  0.800        334928    336576
  11 4358-27-07.laz           before                4   1   3    0         --+-        0.750      0.576      0.809  0.809        183105    183105
  12 4342-22-12.laz           before_mnrv           6   4   2    0       +++-+-        0.833      0.552      0.799  0.833          6607      7362
  13 4342-17-12.laz           before_mnrv           6   4   2    0       +++-+-        0.833      0.495      0.801  0.833          6601      6943
  14 4358-26-05.laz           before                4   2   2    0         +--+        0.500      0.475      0.932  0.932        322758    322758
  15 4342-18-11.laz           before_mnrv           6   3   3    0       +++---        0.667      0.993      0.923  0.993          4866     37060
  16 4342-22-11.laz           before_mnrv           6   3   3    0       +++---        0.667      0.993      0.920  0.993          4872     37405
  17 4342-24-11.laz           before_mnrv           6   3   3    0       +++---        0.667      0.993      0.923  0.993          4870     37573
  18 5158-01-05.laz           before                7   3   2    2      X+-+-X+        1.000      0.285      0.553  1.000        322610    333336
  19 5142-14-64.laz           before                5   2   2    1        +-+-X        1.000      0.636      0.364  1.000          8297     20693
  20 5142-15-64.laz           before                5   2   2    1        +-+-X        1.000      0.642      0.366  1.000          8296     20737
  21 4342-28-64.laz           before                5   2   2    1        +X-+-        1.000      0.497      0.358  1.000          8289      9225
  22 5142-01-60.laz           before                5   2   2    1        +-+-X        1.000      0.440      0.342  1.000          3333     16831
  23 5142-15-63.laz           before                5   2   2    1        +-+-X        1.000      0.656      0.383  1.000          2988     12471
  24 4358-28-01.laz           before                5   2   2    1        -+-+X        1.000      0.571      0.349  1.000          2976      5865
  25 4358-29-06.laz           before                5   2   3    0        -+-+-        1.000      0.167      0.537  1.000        337376    337376
  26 4342-28-56.laz           before                4   2   2    0         +-+-        1.000      0.011      0.395  1.000        142968    144116
  27 4342-02-32_d_c.laz       before_battlecreek    5   2   3    0        -+-+-        1.000      0.398      0.060  1.000         73085     73085
  28 4342-03-32_a_b.laz       before_battlecreek    5   2   3    0        -+-+-        1.000      0.398      0.040  1.000         73084     73084
  29 4342-03-32_b_a.laz       before_battlecreek    4   2   2    0         -+-+        1.000      0.578      0.598  1.000         71663     72432
  30 4342-02-32_c_d.laz       before_battlecreek    4   2   2    0         -+-+        1.000      0.578      0.568  1.000         71663     72459
  31 4342-20-07.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.434      0.256  1.000         11264     13019
  32 4342-23-07.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.409      0.273  1.000         11262     12760
  33 elbaext_gen1_merged.laz  before                6   3   3    0       +-+-+-        1.000      0.346      0.235  1.000         11246     12305
  34 4342-18-06.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.428      0.244  1.000         10677     11534
  35 4342-28-61.laz           before                4   2   2    0         +-+-        1.000      0.112      0.400  1.000         10218     11309
  36 4342-29-61.laz           before                4   2   2    0         +-+-        1.000      0.120      0.407  1.000         10217     11393
  37 4342-21-64.laz           before                4   2   2    0         +-+-        1.000      0.429      0.360  1.000          8291      8619
  38 4342-30-64.laz           before                4   2   2    0         +-+-        1.000      0.515      0.365  1.000          8289      9398
  39 4342-29-64.laz           before                4   2   2    0         +-+-        1.000      0.506      0.359  1.000          8289      9311
  40 4342-16-04.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.225      0.277  1.000          5317      6340
  41 4342-20-02.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.266      0.258  1.000          4369      5083
  42 2742-12-53.laz           before_carlton        5   3   2    0        +-+-+        1.000      0.209      0.007  1.000          4072      4072
  43 4342-23-01.laz           before_mnrv           6   3   3    0       -+-+-+        1.000      0.162      0.246  1.000          3967      4427
  44 4342-20-62.laz           before                4   2   2    0         -+-+        1.000      0.624      0.335  1.000          3377      6050
  45 5142-01-62.laz           before                4   2   2    0         -+-+        1.000      0.425      0.332  1.000          3371      4942
  46 4342-29-62.laz           before                4   2   2    0         -+-+        1.000      0.495      0.338  1.000          3371      5280
  47 5142-02-62.laz           before                4   2   2    0         -+-+        1.000      0.406      0.333  1.000          3370      4857
  48 4342-30-62.laz           before                4   2   2    0         -+-+        1.000      0.479      0.340  1.000          3370      5195
  49 5142-04-62.laz           before                4   2   2    0         -+-+        1.000      0.367      0.327  1.000          3369      4687
  50 5142-11-62.laz           before                4   2   2    0         -+-+        1.000      0.213      0.335  1.000          3366      4092
  51 4342-27-59.laz           before                4   2   2    0         -+-+        1.000      0.564      0.357  1.000          3279      5518
  52 4342-20-63.laz           before                4   2   2    0         +-+-        1.000      0.049      0.368  1.000          3064      3106
  53 4342-28-63.laz           before                4   2   2    0         +-+-        1.000      0.305      0.372  1.000          2986      3922
  54 4342-30-63.laz           before                4   2   2    0         +-+-        1.000      0.354      0.381  1.000          2986      4096
  55 4342-29-63.laz           before                4   2   2    0         +-+-        1.000      0.330      0.372  1.000          2986      4009
  56 4342-27-63.laz           before                4   2   2    0         +-+-        1.000      0.280      0.379  1.000          2986      3835
  57 4342-25-63.laz           before                4   2   2    0         +-+-        1.000      0.228      0.382  1.000          2985      3660
  58 4342-22-63.laz           before                4   2   2    0         +-+-        1.000      0.146      0.378  1.000          2983      3397
  59 4342-21-63.laz           before                4   2   2    0         +-+-        1.000      0.117      0.375  1.000          2982      3309
  60 4358-23-01.laz           before                4   2   2    0         -+-+        1.000      0.643      0.356  1.000          2980      5417
  61 4358-29-01.laz           before                4   2   2    0         -+-+        1.000      0.554      0.343  1.000          2976      4899
  62 5158-01-01.laz           before                4   2   2    0         -+-+        1.000      0.481      0.354  1.000          2972      4550
  63 5158-04-01.laz           before                4   2   2    0         -+-+        1.000      0.419      0.353  1.000          2967      4288
  64 5158-10-01.laz           before                4   2   2    0         -+-+        1.000      0.274      0.347  1.000          2956      3765
  65 4358-29-02.laz           before                4   2   2    0         -+-+        1.000      0.558      0.399  1.000          2899      4822
  66 4358-30-02.laz           before                4   2   2    0         -+-+        1.000      0.540      0.407  1.000          2899      4736
  67 5158-01-02.laz           before                4   2   2    0         -+-+        1.000      0.484      0.394  1.000          2899      4477
  68 4358-26-03.laz           before                4   2   2    0         +-+-        1.000      0.211      0.322  1.000          2831      3431
  69 4358-29-03.laz           before                4   2   2    0         +-+-        1.000      0.294      0.325  1.000          2830      3693
  70 4358-30-03.laz           before                4   2   2    0         +-+-        1.000      0.321      0.314  1.000          2830      3780
  71 4358-32-03.laz           before                4   2   2    0         +-+-        1.000      0.371      0.306  1.000          2829      3954
  72 1158-31-59.laz           before_ne             4   2   2    0         +-+-        1.000      0.316      0.380  1.000          1912      2477
  73 4342-26-08.laz           before_mnrv           4   2   2    0         -+-+        1.000      0.185      0.360  1.000          1843      2169
```

Caveats on reading this table, stated because they bit during the survey:

1. `parity_alias = 1.000` does **not** mean useless — `4358-29-06` and `4342-28-56` both
   score 1.000 and are the two best tiles on disk (Sec 2). Use Sec 2's crossed-design table
   as the primary screen and this one as supporting detail.
2. A low `parity_alias` with a high `|r(sgn,t)|` is not progress: `4342-18-11`,
   `4342-22-11` and `4342-24-11` all pair 0.667 with 0.993.
3. Everything in `before_mnrv` is `lidar_swmn2010` and everything in `before_carlton` /
   `before_ne` / `before_battlecreek` is its own survey. Only `data/before/*` is elbaext's
   acquisition.

---

## 8. The plain answer

**No local tile solves the problem as a drop-in replacement for elbaext, because none of
the good ones has a gen2.** That is the binding constraint, and it is not close: every
3DEP cloud on disk stops at a registered site's boundary.

What the local data *does* contain, stated straight:

1. **The confound can be broken within elbaext's own acquisition, on data already on
   disk** — `4358-29-06` and `4342-28-56` both fly both headings in two flights 3.9 d and
   39.6 h apart, and `4342-28-64` carries a 270 deg cross line flown 916 s into elbaext's
   own sortie. All three are `lidar_semn2008`.
2. **Only a gen1-internal estimator can use them without a download.** Swath-to-swath
   differencing (`swathdiff`, `align_swaths`) and the pooled-gen1 reference substitution
   already validated in `analysis/CROSS_LINE_FIT.md` need no gen2. Anything phrased in
   terms of `d_mm_corr` against the gen2 surface does.
3. **`4342-28-64` costs nothing at all to start**: its CSF ground cache is built
   (6,941,881 class-2 points, line 10010 carrying 640,397 of them) and
   `analysis/crossline_fit.py` already runs on it.
4. **The battlecreek 4-tile block is the only crossed design with both epochs present** —
   but it is `lidar_metro2011`, so it answers the question for a different survey.
5. **What would be needed for a like-for-like elbaext replacement**: a 3DEP fetch over
   `4358-29-06` (Winona, LAS header extent x[592494,595045] y[4882931,4886438]) or
   `4342-28-56` (Wabasha, x[557460,559991] y[4886004,4889497]), then the full tile build. That is a download of the
   order of 10^8 points plus a CSF pass and a pipeline run — explicitly out of scope here,
   and recorded only so the cost is visible.

---

## 9. Reproducing this

Scripts used are in the session scratchpad, not committed (they are one-offs). No code in
the repo was modified and no data was written; this note is the only file added. The two
readers are eight lines each and their method is restated in Sec 1.
`analysis/investigations/swath_degeneracy_breaking/degeneracy_flightline_inventory.py`
is the repo's prior inventory of the same clouds — it asked whether a PERPENDICULAR tie
line exists (it does, psid 10010) and covered the 13 tiles then on disk. This survey covers
the 74 gen1 files on disk today and asks the different question of whether heading is
separable from flight order and flight session.
