# Absolute-elevation test at NGS leveled benchmark DG8385 ("9 DRL")

**Purpose.** Test the *absolute* NAVD88 vertical accuracy of the two lidar epochs
(gen1 = 2008 MN DNR SE-Minnesota; gen2 = 2021 USGS 3DEP SE-Driftless) against a
geoid-model-independent surveyed height, using **individual ground returns**
(classification 2) from each cloud at the one lidar-usable flat/stable NGS mark
in the study area.

**Verification status.** Every benchmark value below is quoted from the
authoritative NGS datasheet (retrieved 2026-08-22) or the NGS geoid API; both are
cited in the trail. Every lidar number is produced by the runnable script
`flagpole_absolute_test.py` in this directory and can be reproduced. Nothing is
fabricated; the dominant limitation (horizontal position + no resolvable slab) is
stated plainly and is in fact the headline caveat.

Compiled 2026-08-22.

---

## Bottom line (read this first)

1. **The epoch difference is the trustworthy result, and it is clean.**
   gen2 reads **+73 mm higher than gen1** (median over a ±3 m position grid, std
   11 mm). The GEOID03→GEOID18 update alone predicts **+68 mm** at this point.
   **The geoid-model change explains essentially all of the offset** (residual
   +5 mm). This epoch difference is *robust to horizontal mis-location* because
   both epochs are evaluated at the same (x,y) on the same graded surface, so a
   ±3 m position error shifts both identically and cancels.

2. **The absolute offsets are consistent with LOW but are position-limited.**
   gen1 − benchmark = **−180 mm**, gen2 − benchmark = **−102 mm**, each with a
   **±~180 mm** uncertainty *dominated by the ±3 m horizontal position* of the
   mark on a ~6% slope. Both epochs read low, and gen1 lower than gen2 by the
   geoid amount, but neither absolute number can be pinned to better than
   ~one-slope-step because we do not know where the disk is to better than 3 m.

3. **There is no resolvable concrete slab.** Measured, not assumed: the surface
   at the mark is a **smoothly graded lawn** (slope ~6%, residual scatter ~65 mm
   about a local plane = grass roughness, not a smooth slab). The flagpole-base
   slab is sub-meter and is averaged into the surrounding turf; discrete-return
   lidar cannot isolate it. We therefore fit a robust local plane and evaluate it
   *at the mark* rather than "isolating the slab."

4. **gen1 is genuinely sparse here** — 14 ground returns within 3 m, 39 within
   5 m — versus gen2's 210 / 593. The gen1 plane-at-mark fit uncertainty is still
   small (12 mm) because the fit averages many points over a plane, but the
   sparseness is real and is why we do not push to a tighter radius.

---

## 1. The benchmark: DG8385 "9 DRL" (NGS datasheet, verified 2026-08-22)

| Field | Value | Source |
|---|---|---|
| PID / designation | DG8385 / "9 DRL" | datasheet |
| NAD83(1986) position | 44°07′12.89″N, 092°00′12.73″W (**HD_HELD1**, ±3 m horiz) | datasheet |
| UTM 15N (EPSG:26915) | E 579 729.20, N 4 885 711.24 | pyproj (this work) |
| **NAVD88 ortho height** | **223.352 m** | datasheet |
| Height method | **LEVELED** — differential leveling, adjusted Feb 2005 | datasheet |
| Vert order | SECOND CLASS I | datasheet |
| Marker / setting | **survey disk (DD)** set in **mat foundation / concrete slab** (code 35), **flagpole base** | datasheet |
| **Mark setting detail** | **disk RECESSED 1 INCH** (25.4 mm) below the slab surface | datasheet |
| Stability | C — may hold, subject to surface motion | datasheet |
| Location | 3.7 km NE of Elba, at the MN DNR Whitewater WMA headquarters | datasheet |
| GEOID18 height (geoid N) | −30.080 m | datasheet & NGS geoid API |

**Leveled = geoid-independent.** The 223.352 m orthometric height was obtained by
running levels, not from GPS + a geoid model, so it is **true NAVD88 regardless
of geoid model**. This is the preferred kind of absolute control and is exactly
what lets it referee gen1 (GEOID03) against gen2 (GEOID18).

**Recessed-disk correction.** The disk sits 1 inch (25.4 mm) *below* the slab
surface. A lidar ground return lands on the **slab/turf surface**, not the disk
face, so the expected lidar surface height is

> H_surface = 223.352 m + 0.025 m = **223.377 m**.

All "epoch − benchmark" offsets below are taken against **223.377 m**, not the
raw 223.352 m mark height. (This correction moves both epochs 25 mm more negative;
it does not affect the gen1−gen2 difference.)

## 2. Geoid framing (NGS geoid API, verified 2026-08-22)

At DG8385 the NGS geoid models give:

- **GEOID03** N = **−30.012 m** (NGS model 3) — the gen1 (2008 MN DNR) geoid.
- **GEOID18** N = **−30.080 m** (NGS model 14) — the gen2 (2021 3DEP) geoid;
  matches the datasheet's −30.080 exactly.

A lidar reports orthometric height H = h_ellipsoid − N. For the **same physical
ground point**, the epoch difference from the geoid update alone is

> H_gen2 − H_gen1 = N_GEOID03 − N_GEOID18 = −30.012 − (−30.080) = **+0.068 m**.

So **on the geoid update alone, gen2 should read 68 mm higher than gen1.** This is
the quantity the epoch difference below is tested against.

## 3. Method

1. Convert DG8385 to UTM 15N; confirm coverage. The mark is covered by the
   existing gen1 tile `data/before/4342-29-64.laz` and by the gen2 full-density
   cloud `data/after/3dep2021_fulldensity.laz` (no fetch needed — 728 gen1 /
   ~9200 gen2 points within a 30 m box). The DNR HQ building (class 6) is >10 m
   away and is excluded by the small analysis radius; class 2 only is used.
2. Extract **individual class-2 (ground) returns** within radii of 3/5/8 m.
3. Because no flat slab is resolvable (§Bottom-line 3), fit a **robust local
   plane** to the ground returns and **evaluate it at the mark (dx=dy=0)** to
   remove the slope×position error. Report roughness (MAD-scaled residual std),
   a bootstrap SE on the fit, and the elevation SE implied by ±3 m horizontal
   position on the local slope.
4. Compare gen1 and gen2 plane-at-mark heights to 223.377 m and to each other;
   scan the center over a ±3 m grid to show the epoch difference is
   position-robust.

Primary radius = **5 m** (enough gen1 points; tight enough to stay off the
building and away from curbs).

## 4. Results (script output, r = 5 m)

| epoch | ground returns ≤3 m / ≤5 m | z@mark (plane) | local slope | roughness | SE (fit) | SE (horiz ±3 m) | SE total |
|---|---|---|---|---|---|---|---|
| gen1 (2008, GEOID03) | 14 / 39 | 223.198 m | 6.2% | 64 mm | 12 mm | 185 mm | **185 mm** |
| gen2 (2021, GEOID18) | 210 / 593 | 223.275 m | 5.7% | 61 mm | 2 mm | 171 mm | **171 mm** |

**Absolute offsets** (vs expected surface 223.377 m; ± dominated by horizontal position):

- **gen1 − benchmark = −180 mm (LOW), ±185 mm**
- **gen2 − benchmark = −102 mm (LOW), ±171 mm**

**Epoch difference** (robust to horizontal position; ±3 m center grid, n = 49):

- **gen2 − gen1 = +73 mm** (median; range +43…+96 mm, std 11 mm)
- GEOID03→GEOID18 prediction = **+68 mm**
- **residual (measured − geoid) = +5 mm**

## 5. Verdict

**Epoch difference (high confidence).** gen2 sits **+73 mm above gen1** at this
benchmark, and the GEOID03→GEOID18 model update predicts +68 mm of that. The
geoid change accounts for essentially the entire epoch offset (5 mm residual, well
inside the 11 mm scatter of the position scan and the surface roughness). This is
a clean, independent confirmation that the ~68–73 mm gen1↔gen2 vertical offset
seen elsewhere in this project is the **geoid-model difference**, not real ground
change and not a gross registration error. It also means that *after* accounting
for the geoid, the two epochs agree to within a few mm at an absolute leveled
mark — a strong internal consistency check.

**Absolute accuracy (position-limited, but directionally clear).** Both epochs
read **low** relative to the leveled surface — gen1 by ~180 mm, gen2 by ~100 mm —
and the gen2 low bias (~100 mm) is within the range expected once the ±180 mm
position uncertainty is admitted. These absolute numbers are **not** trustworthy
to better than ~one slope-step (±~180 mm) because the mark's horizontal position
is only surveyed to ±3 m (hand-held GPS) and it sits on a 6% grade: a 3 m along-
slope error is ±180 mm of height by itself. We report the *sign* (both low) and
the *magnitudes* with that caveat; we do **not** claim a mm-level absolute bias
from this single mark.

**On which epoch is "more right" absolutely:** gen2 (GEOID18, the current model)
is closer to the leveled truth than gen1 (GEOID03), consistent with the geoid
update having moved 2008-era heights ~68 mm in the correct direction. But given
the position limit, treat the absolute offsets as *order-of-magnitude* (both
epochs low by ~0.1–0.2 m at this site) and the **epoch difference as the
quantitative result**.

## 6. Honest limitations

- **No resolvable slab.** The flagpole-base slab is sub-meter; the lidar sees a
  graded lawn (65 mm roughness). This is a *lawn* absolute test, not a *slab*
  one. The plane-at-mark estimator handles the slope but inherits grass
  roughness and any turf-vs-slab micro-offset (unknown, ≤ a few cm).
- **Horizontal position dominates the absolute error budget** (±3 m → ±180 mm on
  this slope). A better absolute test needs a mark on flat ground *or* a mark
  whose (x,y) is known to sub-decimeter. Neither exists in this AOI (see
  `ABSOLUTE_ELEVATION_REFS.md`).
- **Single mark.** One point cannot separate a whole-collection vertical bias
  from a local artifact. The value here is the geoid-difference confirmation,
  which one clean mark *can* deliver because it is a within-mark, position-robust
  comparison.
- **Stability "C"** (subject to surface motion): a decade of slab settlement
  could add a few mm of real gen1↔gen2 change; it is not separable from the 5 mm
  residual and is not claimed.

---

## Appendix: sources & reproduction

**Benchmark & geoid (authoritative):**
- NGS datasheet DG8385: `https://www.ngs.noaa.gov/cgi-bin/ds_mark.prl?PidBox=DG8385`
  — leveled NAVD88 223.352 m; disk recessed 1 inch; setting 35 (mat/slab),
  flagpole base; GEOID18 N = −30.080 m; ±3 m hand-held horizontal.
- NGS geoid API `geodesy.noaa.gov/api/geoid/ght`: GEOID03 (model 3) N = −30.012 m;
  GEOID18 (model 14) N = −30.080 m at 44.120247 N, −92.003536 W.

**Lidar (this project):**
- gen1: `data/before/4342-29-64.laz` (2008 MN DNR, NAVD88/GEOID03), class 2.
- gen2: `data/after/3dep2021_fulldensity.laz` (2021 3DEP, NAVD88/GEOID18), class 2.

**Reproduce:**
```
cd <repo>
PROJ_DATA=/usr/share/proj ./lidar-icp/bin/python analysis/ridgelines/flagpole_absolute_test.py
```
Prints per-epoch return counts, plane-at-mark heights at r = 3/5/8 m, the three
offsets, and the position-robust epoch-difference grid vs the geoid prediction.
