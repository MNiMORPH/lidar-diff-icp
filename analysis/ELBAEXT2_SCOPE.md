# ELBAEXT2 scope — reaching a surveyed absolute datum from Elba

**Question.** Our vertical datum floats: there is no absolute tie at Elba. gen1
(2008 MN DNR) is corrected per flight line, so if a gen1 flight line passes over
both Elba *and* a surveyed 3DEP QA checkpoint, tying that line at the checkpoint
constrains its offset at Elba. How large an extension buys that tie?

**Answer, in one line.** It is not an extent problem at all. Flight lines 133–138
run **north–south**, and every checkpoint is displaced **east–west** by 3.1–15.6 km —
cross-track distance no extension can shrink. But the lines **chain**, and the chain
from Elba to the western checkpoints is **five links long and needs exactly two extra
gen1 tiles, both of which are now on disk**. Candidate extents A/B/C are all far
larger than the problem requires.

Compiled 2026-08-26. Every number below is labelled with its source. Nothing is
committed to git; probe tiles live in the session scratchpad.

---

## Bottom line (read this first)

1. **Lines 135–138 pass over none of the six checkpoints, and cannot.** They are
   N–S lines with a measured swath half-width of ~710–730 m. The closest any of the
   four comes to any checkpoint is **3 091 m** (PSID 138 → 3089_VVA) — 4.2× the
   half-width. Extending the study area north keeps you on the same lines but does
   not move them east or west, so **no extent, of any size, puts 135–138 over a
   checkpoint.** This is measured, not extrapolated (see §2).

2. **The checkpoints *are* directly overflown — by other lines.** Probing the four
   tiles that contain checkpoints shows each sits essentially under a single line's
   nadir track: **2210 and 3056 on PSID 128** (147 m and 131 m off-track),
   **2024 on PSID 129** (32 m — effectively nadir), **2036 on PSID 144** (227 m).
   Each has **173–201 class-2 gen1 ground returns within 10 m**, all near-nadir
   (|scan_angle_rank| median 0–5°), and only one line covers each — a clean,
   unmixed single-line tie.

3. **The chain to Elba is short, and the western half is already paid for.**
   Lines are spaced a uniform ~960 m and every adjacent pair overlaps by
   **1.31–1.48 km²** (measured, 2 m cells) — ample for Nuth & Kääb. The path is
   `128 – 129 – 130 – 131 – 132 – 133`, and 133 is elbaext's reference swath, already
   aligned to 134–138. **Five links from the checkpoint line to Elba's swath network.**
   All the overlaps needed live in tiles `4342-29-61`, `4342-29-62`, `4342-29-63`,
   `4342-29-64` — one tile row, all four now local.

4. **An independent eastern closure check costs two more tiles.** 2036 sits on line
   144; `138 → … → 144` is six links through tiles `4358-29-01` and `4358-29-02`
   (not yet fetched) plus `4358-29-03` (fetched). West tie and east tie should agree
   at Elba — that disagreement is the honest error bar on the whole chain, and there
   is no other way to get one (a chain has zero redundancy by construction; see
   `align_swaths` docstring, `src/lidar_diff_icp/coreg.py:444`).

5. **No CSF is needed for the chain.** `coregister_swaths`
   (`src/lidar_diff_icp/coreg.py:530`) selects points by `~isin(classification,
   (5,6,9))` — the *vendor* classification. The corridor chain can be solved on the
   raw tiles with `coreg.align_swaths`, skipping the multi-hour CSF that dominates
   the cost of extents A/B/C.

6. **gen2 is barely needed.** The gen1 chain is gen2-free. gen2 is wanted only as an
   independent check that the 3DEP surface itself hits the checkpoint — a 200 m box
   is ~0.8 M points, ~4 MB. Six boxes ≈ 24 MB. gen2 over a whole corridor or over
   extents A/B/C is 5.9–21.7 GB and would not fit (§6).

7. **Caution — the tie is available but not yet a number.** A naive local fit at the
   checkpoints is radius-sensitive and would mislead (§4). Extracting a defensible
   offset needs the project's slope-normal ground estimator and the GEOID03→GEOID18
   conversion. Do not treat the raw values in §4 as the answer.

---

## 1. What was fetched (be gentle with the server)

| what | source | cost |
|---|---|---|
| tile name → containing tile, 6 checkpoints | `tiles.find_tile` on the **cached** statewide centroid CSV `data/mn_tile_centroids.csv` (55 296 tiles, dated 2026-08-18) | **0 requests** |
| tile enumeration for extents A/B/C | same cached CSV | **0 requests** |
| bbox of 7 candidate tiles | `tiles.header_bbox` — HTTP **range read, 512 bytes each**, 2 s apart | 7 × 512 B |
| 4 gen1 tiles, downloaded | `tiles.download_tile`, 5 s apart | **105.9 MB** |
| 3DEP EPT project per checkpoint | `threedep.find_projects` on the **cached** `data/ept_boundaries.json` | **0 requests** |
| gen2 point data | — | **nothing fetched** |

Tiles downloaded, to the session scratchpad (not `data/`, not committed):

| tile | county dir | size | points | bbox (UTM 15N) |
|---|---|---|---|---|
| `4342-29-62` | winona | 23.6 MB | 7 485 709 | E 572 492–575 032, N 4 882 681–4 886 180 |
| `4342-29-61` | **olmsted** | 32.8 MB | 10 549 493 | E 569 992–572 530, N 4 882 654–4 886 151 |
| `4342-28-61` | **wabasha** | 29.6 MB | 10 406 516 | E 569 955–572 492, N 4 886 125–4 889 623 |
| `4358-29-03` | winona | 19.9 MB | 7 455 152 | E 584 994–587 539, N 4 882 830–4 886 334 |

Note the county directories differ (olmsted, wabasha, winona) — `find_tile` is
correct across county lines but `download_tile` still needs the right county dir.

The two remaining checkpoints, 2099 (`4342-26-61`) and 3089 (`4358-26-02`), were
**not** probed. Their covering lines are stated below as *unverified inference*.

---

## 2. Flight-line geometry — the crux

**Method.** For each `point_source_id`, the nadir ground track is fitted from
returns with `|scan_angle_rank| <= 1` (near-vertical beams sit under the aircraft),
binned by `gps_time` and fitted linearly in x and y. This is the heading fit of
`analysis/ridgelines/gen1_save_angles_slope.py` made robust to tile clipping — a
plain centroid track is biased by up to ~1 km when a swath is cut by the tile edge,
which is why the nadir version is used here. Swath half-width is the 99.5th
percentile of |cross-track| over **all** returns of that line.

**Self-consistency check:** lines 128, 129 and 130 appear in two independent tiles
(`4342-29-61` and `4342-28-61`, 3.5 km apart in northing). Their fitted track
eastings agree to **7–20 m**, and half-widths to **1–3 m**. The fit is sound.

**Fitted nadir tracks** (easting where the track crosses N 4 884 126; source: the
four probe tiles and `data/derived/elbaext/beam_offset_table.parquet`):

| PSID | track E | half-width | | PSID | track E | half-width |
|---|---|---|---|---|---|---|
| 128 | 570 345 | 717 m | | 136 | 578 054 | 720 m |
| 129 | 571 311 | 696 m | | 137 | 579 046 | 725 m |
| 130 | 572 275 | 698 m | | 138 | 579 989 | 730 m |
| 131 | 573 286 | 742 m | | *(139–143 unprobed)* | | |
| 132 | 574 222 | 720 m | | 144 | 585 756 | 698 m |
| *133* | *~575 18X (see note)* | — | | 145 | 586 764 | 717 m |
| 134 | 576 129 | 707 m | | | | |
| 135 | 577 119 | 707 m | | | | |

Spacing is uniformly **935–1 011 m** (mean ~960 m), so adjacent swaths overlap by
~460 m. Headings alternate 179–183° / 357–359° — a plain N–S boustrophedon. Speed
68–84 m/s.

*Note on 133:* its nadir band (~E 575 18X, interpolated from 132 and 134) falls in
the strip E 575 032–575 600, which lies **between** tile `4342-29-62` and the elbaext
grid's west edge. It is present in the already-local `4342-29-63.laz` but was cropped
out of the elbaext product. Any corridor build must include that strip.

**Cross-track distance, line to checkpoint** (m; `*` = inside that swath):

| checkpoint | 128 | 129 | 130 | 131 | 132 | 134 | 135 | 136 | 137 | 138 | 144 | 145 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2210_NVA | **−147\*** | −819 | 1 783 | −2 791 | 3 729 | 5 636 | −6 626 | 7 560 | −8 550 | 9 495 | 15 263 | −16 271 |
| 3056_VVA | **−131\*** | −836 | 1 800 | −2 805 | 3 744 | 5 653 | −6 644 | 7 575 | −8 564 | 9 511 | 15 281 | −16 288 |
| 2024_NVA | −957 | **−32\*** | 983 | −1 899 | 2 899 | 4 839 | −5 846 | 6 724 | −7 685 | 8 676 | 14 486 | −15 480 |
| 2036_NVA | −15 637 | 14 672 | −13 707 | 12 691 | −11 760 | −9 854 | 8 864 | −7 929 | 6 936 | −5 995 | **−227\*** | −780 |
| 2099_NVA | −1 765 | 739 | **191\*** | −952 | 2 057 | 4 051 | −5 086 | 5 872 | −6 784 | 7 853 | 13 732 | −14 704 |
| 3089_VVA | −12 707 | 11 676 | −10 749 | 10 008 | −8 890 | −6 888 | 5 848 | −5 076 | 4 171 | −3 091 | 2 799 | −3 766 |

**Reading the table.** For lines **135, 136, 137, 138** — the four in the Elba pilot —
every entry is ≥ 3 091 m against a half-width of ~730 m. They reach nothing, and
extending the extent cannot change a cross-track distance. 134 and 133 are no better.
The same holds structurally: at 960 m line spacing, a checkpoint 6 km cross-track is
six flight lines away, not a bigger box away.

**Two entries need a caveat.** 2210/3056/2024/2036 are *measured* — the covering line
was found in the tile that contains the checkpoint. **2099 on line 130 (+191 m) is
extrapolated** ~9 km north of the data and assumes line 130 runs that far; the nadir
fit residual is 16–23 m over 3.5 km, so the extrapolation is probably good to
~100 m against a 698 m half-width, but *whether line 130 exists at N 4 893 677 is
unverified and needs tile `4342-26-61`.* **3089 has no covering line among those
probed** (nearest 144 at 2 799 m); its line lies further east and is unknown without
`4358-26-02`.

**On along-track extent.** The `gps_time` gaps between successive lines at Elba are
1 031–1 926 s (with breaks of 6 298 s after 138 and 8 929 s after 129), which at
~80 m/s is 82–154 km of flight path between passes. That is far more than a
back-and-forth over our 4 km window, so the lines are long — but their true ends
**cannot be bounded from our data** and it does not matter, because the limit is
cross-track, not along-track.

**One survey block, one night.** `global_encoding = 0` (LAS 1.1), so `gps_time` is
GPS **seconds-of-week**. Lines 146→127 run monotonically from 225 866 s to 266 376 s
— GPS week-day 2, 14:44 UTC through week-day 3, 01:58 UTC — with easting decreasing
monotonically as time increases. The whole set is **one continuous progression flown
east-to-west in a single night's sortie sequence**, not separate campaigns. That is
what makes chaining physically reasonable rather than hopeful.

---

## 3. The chain — measured overlaps

`align_swaths` solves a free-network least-squares over pairwise Nuth & Kääb
observations. Measured overlap area per adjacent pair (2 m cells, classes 5/6/9
excluded, i.e. exactly the `coregister_swaths` selection):

| link | overlap | tile measured in |
|---|---|---|
| 128–129 | 1.459 km² / 1.443 km² | `4342-29-61` / `4342-28-61` |
| 129–130 | 1.395 km² / 1.378 km² | `4342-29-61` / `4342-28-61` |
| 130–131 | 1.310 km² | `4342-29-62` |
| 131–132 | 1.477 km² | `4342-29-62` |
| 132–133 | 1.383 km² | `4342-29-62` |
| 133–134 … 137–138 | solved | already in `data/derived/elbaext/corrections.json` |
| 143–144 | 1.421 km² | `4358-29-03` |
| 144–145 | 1.349 km² | `4358-29-03` |

No non-adjacent overlaps exist (checked: all zero), and no E–W tie lines appear in
any probed tile. The network is a **pure chain** — misclosure is identically zero and
carries no information, so the only real check on accumulated error is an
independent tie from the other side (§ bottom line 4).

**Why the corridor should sit at Elba's latitude.** Each swath's offset varies along
track (`fit_along_track_drift`, `coreg.py:380`, a spline in `gps_time`). All six
western overlaps and the Elba swaths lie in the same northing band
(N ≈ 4 882 650–4 886 240), so the chain is evaluated at essentially constant
along-track position and the drift term stays common rather than accumulating. A
corridor at a different latitude would not have that property. *This is a reasoned
argument from the code, not a measured result.*

**Note for the full pipeline.** `fit_along_track_drift` regresses against `Zref`, the
gen2 grid (`pipeline.py:718–723`). So a *full DoD* over a corridor would need gen2
there. The *chain tie* does not — it is `align_swaths` only, gen2-free.

---

## 4. What a tie would look like — and why it is not yet a number

gen1 returns within 10 m of each checkpoint, by line (measured, probe tiles):

| checkpoint | line | returns | class-2 | median |scan| | surveyed NAVD88(GEOID18) |
|---|---|---|---|---|---|
| 2210_NVA | 128 | 229 | 201 | 3° | 349.288 m |
| 3056_VVA | 128 | 209 | 187 | 4° | 353.259 m |
| 2024_NVA | 129 | 347 | 179 | 0° | 344.735 m |
| 2036_NVA | 144 | 201 | 173 | 5° | 353.119 m |

Sampling is ample and near-nadir, and exactly one line covers each point. **But a
local plane fit to the class-2 ground is strongly radius-dependent**, so any naive
difference is an artifact of the radius chosen:

| checkpoint | R=5 m | R=10 m | R=20 m | plane residual RMS at R=10 |
|---|---|---|---|---|
| 2210_NVA | −200 mm | −589 mm | −1 169 mm | 543 mm |
| 3056_VVA | −91 mm | −72 mm | −100 mm | 186 mm |
| 2024_NVA | −307 mm | −467 mm | −436 mm | 209 mm |
| 2036_NVA | −185 mm | −430 mm | −808 mm | 384 mm |

(RAW gen1 − surveyed, **no** geoid conversion, **no** swath alignment, **no**
slope-normal estimator.) A 543 mm plane residual on ground with a 2.4° fitted slope
is not a flat surface. At 2210, the class-2 elevations within 5 m span
348.83–349.26 m and the surveyed 349.288 m sits at their **p95** — the mark is on a
local high (road crown or shoulder) with ground falling away inside 10 m. These
checkpoints are sited for *survey* convenience, not for lidar patch extraction.

**Therefore: do not quote these numbers.** A defensible tie needs (1) the project's
slope-normal ground estimator at a small, justified radius rather than a plane fit,
(2) the GEOID03→GEOID18 conversion computed per point via
`references.geoid_difference` (not a constant), and (3) the per-swath alignment
applied. All three exist in the package. That is the next task, not this one.

---

## 5. gen2 (2021 3DEP)

All six checkpoints resolve to the **same** EPT project as Elba —
`MN_SEDriftless_2_2021` — from the cached boundary index, so there is no
project-boundary crossing to worry about:

```
https://s3-us-west-2.amazonaws.com/usgs-lidar-public/MN_SEDriftless_2_2021/ept.json
```

(2210 also falls inside `MN_SEDriftless_4_2021`; `MN_FullState` is the statewide
mosaic and is ranked last.) Access is `scripts/fetch_3dep_curl.py --auto`, the
recipe used for elbaext.

**Density and size, derived from the elbaext pull** (`ELBAEXT_BUILD.md`:
415 080 034 pts over the buffered extent E 575 450–580 200 × N 4 882 050–4 886 400 =
20.66 km² → **20.1 pts/m²**; bytes/point from `elbaext_3dep_fd_class2.laz`,
664 737 867 B / 131 988 280 pts = **5.04 B/pt**):

| what | area | points | LAZ |
|---|---|---|---|
| one 200 m box | 0.04 km² | 0.80 M | 4.1 MB |
| **all six boxes** | 0.24 km² | **4.8 M** | **24 MB** |
| row-29 corridor, full | 61.4 km² | 1.23 G | 6.2 GB |
| extent A | 58.3 km² | 1.17 G | 5.9 GB |
| extent B | 118.6 km² | 2.38 G | 12.0 GB |
| extent C | 214.5 km² | 4.31 G | 21.7 GB |

**Only the six 200 m boxes are actually needed**, and even those are a cross-check
on 3DEP rather than an input to the gen1 tie.

---

## 6. Cost table

gen1 per-tile figures are the **measured** mean of the ten tiles now on disk:
**8.09 M points, 25.5 MB** per tile (range 7.12–10.55 M; western tiles run denser).
Tile counts come from the cached centroid index and are boundary-sensitive by ±1–2.
CSF time scales from the brief's datum: elba gen1 7.7 M pts ≈ 15 min tiled 2×2 →
**1.95 min per M points**.

| option | gen1 tiles | to fetch | gen1 pts | gen1 LAZ | gen2 needed | CSF time | peak RAM | fits 46 GB? |
|---|---|---|---|---|---|---|---|---|
| **West chain** (`4342-29-61…64`) | 4 | **0 — on disk** | 32.9 M *(measured)* | 106 MB *(measured)* | 3 boxes, 12 MB | **none needed** | ~2 GB | yes, easily |
| **+ East closure** (`+4358-29-01/02/03`) | 7 | **2** (~51 MB) | 56.5 M | 177 MB | 4 boxes, 16 MB | none needed | ~3 GB | yes |
| A — west box | 16 | 12 | 129 M | 408 MB | 5.9 GB | ~4.2 h | ~31 GB gen2 load | **no** |
| B — A + 2099 | 26 | 22 | 210 M | 663 MB | 12.0 GB | ~6.8 h | ~63 GB gen2 load | **no** |
| C — all six | 37 | 33 | 299 M | 944 MB | 21.7 GB | ~9.7 h | ~114 GB gen2 load | **no** |

**Tiles per extent** (cached index; the elbaext row reproduces the six tiles of
`ELBAEXT_BUILD.md`, confirming the method):

- **A** (E 570 000–580 050, N 4 882 200–4 888 000), 16: `4342-28-61/62/63/64`,
  `4342-29-60/61/62/63/64`, `4342-30-60/61/62/63/64`, `4358-28-01`, `4358-29-01`
- **B** (+ N to 4 894 000), 26: A plus `4342-26-61/62/63/64`, `4342-27-61/62/63/64`,
  `4358-26-01`, `4358-27-01`
- **C** (E to 586 500, N to 4 895 200), 37: B plus `4358-26-02/03`, `4358-27-02/03`,
  `4358-28-02/03`, `4358-29-02/03`, `4358-30-01/02/03`

**RAM flags for A/B/C.** `ELBAEXT_BUILD.md` records two hard limits actually hit on
this machine: monolithic CSF was **OOM-killed** at 17.35 M points with ~18 GB free
(tiling to 4.3–5.0 M points per sub-tile fixed it, holding ~14 GB free); and the
non-streaming gen2 class-2 load of 113.7 M points took **9.4 GB RSS**. Scaling that
second rate (82.7 B/pt) to extent A's class-2 ground (~32 % of 1.17 G pts ≈ 375 M)
gives ~31 GB RSS. Measured this session, `free -h` shows 46 GB total with 8.6 GB
already in use and 37 GB available. **Extent A would not fit; B and C are far
beyond.** They would need a fundamentally different (out-of-core / COPC) gen2 path.

---

## 7. Recommendation

**Build the west chain corridor — tiles `4342-29-61`, `-62`, `-63`, `-64` — and drop
extents A, B and C.**

Reasons, in order of weight:

1. **A/B/C do not buy what they were meant to buy.** They were sized to bring a
   checkpoint inside the grid so lines 135–138 would cross it. Those lines run N–S
   and the checkpoints are 3–16 km east or west of them; no box fixes that (§2).
2. **The chain does buy it, and cheaply.** Five links, every one with 1.3–1.5 km² of
   measured overlap, ending on swath 133 which is already aligned to 134–138 in
   `data/derived/elbaext/corrections.json` (§3).
3. **It costs nothing further to fetch.** All four tiles are on disk. Include the
   E 575 032–575 600 strip so swath 133's nadir band is complete (§2 note).
4. **It needs no CSF and no gen2 corridor** — `coregister_swaths` uses the vendor
   classification, so the multi-hour CSF and the multi-GB gen2 pull that dominate
   A/B/C both vanish (§ bottom line 5, 6).
5. **It delivers two NVA checkpoints, not one.** 2210 (line 128) and 2024 (line 129)
   sit on adjacent links of the same chain, so the 128–129 link becomes
   over-determined — a first, cheap consistency check.

**Then, if the west tie looks sound, spend two more tiles** (`4358-29-01`,
`4358-29-02`) on the eastern chain to 2036 on line 144. West and east ties are
independent paths to the same answer at Elba; their disagreement is the only honest
error bar available on a chain with no internal redundancy.

**Do not chase 2099 or 3089.** 2099 needs `4342-26-61` and is on the same western
chain 9 km north — extra work for a tie the west corridor already provides. 3089 is
a VVA (vegetated) point with a 25–27 cm expected spread and no identified covering
line. Neither justifies the 10-tile jump from extent A to extent B.

---

## 8. Stated plainly: what could not be determined without more data

- Which line covers **3089_VVA** — no probed line comes within 2.8 km. Needs
  `4358-26-02`.
- Whether line **130 actually extends** to N 4 893 677 to cover **2099_NVA**. The
  +191 m cross-track is a 9 km extrapolation of a fit. Needs `4342-26-61`.
- The **PSIDs of lines 139–142**, which sit between 138 and 143 in the eastern chain.
  Inferred to be five lines at ~960 m spacing from the measured 138 (E 579 989) and
  144 (E 585 756) tracks, but not observed. Needs `4358-29-01` and `4358-29-02`.
- The **true along-track ends** of any flight line. Bounded only weakly by inter-line
  `gps_time` gaps, and not needed for this question.
- The **absolute offsets themselves.** §4 explains why the raw numbers are not
  trustworthy and what machinery the real computation requires.
- Point counts and sizes for `4358-29-01` and `4358-29-02` are the ten-tile measured
  mean, not measured for those tiles.
