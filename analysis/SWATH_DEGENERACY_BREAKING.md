# What breaks the per-line across-track degeneracy, and what does not

**Date:** 2026-08-26
**Scripts:** `analysis/degeneracy_flightline_inventory.py` (headings, spacings, water),
`analysis/degeneracy_crossline_geometry.py` (the cross line against the N-S lines),
`analysis/degeneracy_identifiability.py` (what the null space *is*, ground control, gen2),
`analysis/degeneracy_water_surface.py` (a level water surface as a single-line reference)
**Run:** `env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/<script>.py`
**Read, not redone:** `analysis/SWATH_ACROSS_TRACK_TEST.md`, `analysis/ELBAEXT2_SCOPE.md`,
`analysis/ABSOLUTE_BASIS_ELBA.md`, `docs/groundtruth.md`, `analysis/ADDITIONAL_GROUND_CONTROL.md`.
Every table is emitted by the script that computed it through `trust/provenance.py`, so each
column carries its own definition and each parameter names its source. Nothing was fetched;
all thirteen gen1 tiles read were already on disk. `coreg.py` and `pipeline.py` are imported,
not modified.

---

## The headline

**A cross line exists, it is already on disk, and it breaks the degeneracy.** Flight line
**`point_source_id` 10010, heading 270.98°** – due west – runs across lines 136, 137 and 138
in tile `4342-28-64`. In the cells it shares with them the two lines' scan-angle tangents
correlate at **−0.05 to −0.15**, against **−0.92 to −0.96** on every N-S pair, and the
standard error on the combination the N-S chain cannot see falls from **5–7× worse than the
symmetric one to 1.05–1.14×**. That is the degeneracy broken, not mitigated.

**But the degeneracy is not what it was called.** Sec 1 below proves, exactly and
numerically, that the null direction of the N-S network is **a global cross-track tilt of the
mosaic**: the "alternating ±v in `c_s`" and a linear ramp in the per-line nadir offsets are
the *same object*, and that object is `e(x) = g·(x − x₀)`. A tilt is identical for both lines
at a shared ground point, so it cancels in **every** between-line difference – adjacent,
non-adjacent or crossing. **The cross line therefore fixes the body-fixed per-line
coefficients and is blind to the tilt; the tilt needs an absolute external reference.** Those
are two different instruments for two different terms, and both exist in our data.

---

## 1. The null space is a global cross-track tilt (verified, not asserted)

With the per-line model `e_s(θ) = a_s + c_s·tan θ`, each adjacent overlap gives exactly two
numbers: a slope `P_s = (c_s + c_{s+1})/2` and an intercept
`K_s = (a_s − a_{s+1}) + S_s·(c_s − c_{s+1})/2`, where `S_s = tan θ_A + tan θ_B` is the
measured near-constant. Because the lines fly there-and-back, `S_s` alternates in sign with
the body-fixed side, and so does the fitted `h_s`
(`SWATH_ACROSS_TRACK_TEST.md` §0: −2424, +2564, −2626, +2473 m).

Apply `e(x) = g·(x − x₀)` to a synthetic six-line network. It reproduces
`c_s → c_s + g·h_s` (alternating, because `h_s` alternates) **together with**
`a_s → a_s + g·(x_s − x₀)` (a linear ramp, because the tracks are evenly spaced), and the
two conspire so that every line's error at a given easting is the same number:

| pair | P | K | ΔP under the direction | ΔK |
|---|---|---|---|---|
| 0 | +0.000 | +14.116 | +0.000e+00 | +0.000e+00 |
| 1 | +0.000 | −52.687 | +0.000e+00 | +0.000e+00 |
| 2 | +0.000 | +26.107 | +0.000e+00 | +0.000e+00 |
| 3 | +0.000 | +30.697 | +0.000e+00 | +0.000e+00 |
| 4 | +0.000 | −38.799 | +0.000e+00 | +0.000e+00 |

Largest change in **any** observation under a 10 mm/km tilt: **0.000e+00 mm**. The direction
line by line is `Δc = +25.62, −25.62, +25.62, …` and `Δa = 0.00, +9.16, +18.32, …` – the
alternating vector and the linear ramp, at once. At h = 2562 m and 916 m spacing,
`v = g·h` and the ramp step is `g·spacing`.

**Consequence.** The right statement is not "the `c_s` are unidentifiable" but "**the mosaic
tilt is unidentifiable from overlaps, and it projects onto the `c_s` as an alternating
vector**". That reframing is what makes the problem tractable: a cross-track tilt is a
classical airborne-lidar block-adjustment nuisance with a classical fix (ground control),
and it is not a per-line instrument term at all.

Adding a due-west cross line to the same synthetic recovers every line's own `c` exactly, and
returns the **same** value whether or not the tilt is present:

| line | c_true | c_fit, no tilt | c_fit, 10 mm/km tilt added | rms resid, no tilt | rms, tilted |
|---|---|---|---|---|---|
| 0 | +130.40 | +130.40 | +130.40 | 1.82e-14 | 1.94e-14 |
| 1 | +94.71 | +94.71 | +94.71 | 3.77e-15 | 8.96e-15 |
| 2 | −70.37 | −70.37 | −70.37 | 1.08e-14 | 1.10e-14 |
| 3 | −126.54 | −126.54 | −126.54 | 8.99e-15 | 9.43e-15 |
| 4 | −62.33 | −62.33 | −62.33 | 7.89e-15 | 8.30e-15 |
| 5 | +4.13 | +4.13 | +4.13 | 4.53e-15 | 5.53e-15 |

The cross line's own across-track direction is *northing*, so an easting tilt cannot be
absorbed into its two parameters and cannot masquerade as a per-line roll. It is blind to
the tilt rather than confused by it.

---

## 2. Candidate 1 – cross / tie lines. **They exist. Highest value, zero cost.**

Every `point_source_id` in the thirteen local gen1 tiles was enumerated and its nadir ground
track fitted from returns with `|scan_angle_rank| ≤ 1`, the method of
`ELBAEXT2_SCOPE.md` §2. **102,635,205 returns, 19 lines fitted.**

**Eighteen of the nineteen fly the boustrophedon** – per-tile headings 177.68–180.09°
southbound and 358.23–0.40° northbound, consistent with the known 133–138 pattern. (The
script's own summary line says seventeen, because averaging line 135's per-tile headings
across the 180°/−180° wrap corrupts its consolidated value; every one of its per-tile
headings is 179.4–180.1° and it is not a candidate.) **One line does not:**

| psid | tile | returns | heading | speed | class 2 | class 12 | gps_time (s of week) | easting span | northing span | scan_angle_rank |
|---|---|---|---|---|---|---|---|---|---|---|
| **10010** | `4342-28-64` | **747,107** | **270.98°** | 84.6 m/s | **0** | 745,942 | 238,319.99–238,344.42 | 578,000–579,960 | 4,888,888–4,889,709 | −17 … +2 |

**It is the same sortie, not a foreign dataset.** Its 24.43 s of `gps_time` sits between line
138 (237,404 s) and line 137 (243,608 s) in the single monotonic seconds-of-week progression
`ELBAEXT2_SCOPE.md` §2 documents, its fitted speed of 84.6 m/s sits inside the 68–85 m/s of the
N-S lines, and its 99.5th-percentile half-width of 719 m sits inside their 668–748 m. The
5-digit id marks it as the vendor's separate numbering block for tie lines.

**What it overlaps, and how well it separates the lines** (5 m cells, ≥3 returns per
cell-line, the `coregister_swaths` class selection):

| pair | kind | cells | area | tan range A | tan range B | corr(tanA,tanB) | sd(sum) | **SE(q)/SE(p)** |
|---|---|---|---|---|---|---|---|---|
| 135-136 | N-S/N-S | 41,463 | 103.7 ha | −0.287…−0.141 | −0.231…−0.052 | **−0.955** | 0.0111 | **6.46** |
| 136-137 | N-S/N-S | 63,976 | 159.9 ha | +0.070…+0.306 | +0.070…+0.306 | **−0.922** | 0.0219 | **4.95** |
| 137-138 | N-S/N-S | 69,004 | 172.5 ha | −0.287…−0.052 | −0.306…−0.070 | **−0.961** | 0.0161 | **7.09** |
| **136-10010** | **CROSS** | 22,093 | **55.2 ha** | +0.000…+0.287 | −0.306…+0.000 | **−0.056** | 0.1173 | **1.06** |
| **137-10010** | **CROSS** | 44,128 | **110.3 ha** | −0.287…+0.287 | −0.306…+0.035 | **−0.154** | 0.1675 | **1.14** |
| **138-10010** | **CROSS** | 22,903 | **57.3 ha** | −0.306…+0.000 | −0.287…+0.035 | **−0.052** | 0.1206 | **1.05** |

`q = (c_A − c_B)/2` is precisely the combination the N-S chain cannot see. On the cross pairs
it is determined as well as the symmetric combination. Line 137 is seen over its **full**
across-track range, −0.287 to +0.287.

**Coverage, stated plainly.** The cross line reaches lines 136, 137 and 138. It does **not**
reach 135 (its westernmost return is at E 578,000; line 135's swath ends at ≈ E 577,860), and
it is absent from tile `4342-28-63` immediately west, so it does not touch 133–135 anywhere
on disk. That is enough: with `c_136`, `c_137` and `c_138` fixed individually, the measured
pair sums `(c_135+c_136)/2`, `(c_134+c_135)/2` and `(c_133+c_134)/2` deliver 135, 134 and 133
in turn, and the three cross-line values are checked against the three N-S pair sums, so the
network gains real redundancy for the first time.

**What it costs, and what it does not yet give.**

* **Cost to try: none.** `4342-28-64` is on disk: 27,153,945 bytes, mtime 2026-08-15.
* **A CSF pass is required.** Line 10010 carries **zero** vendor class-2 returns – all 745,942
  are class 12 (overlap) – so its ground must be classified. That is the whole cost, and it
  is one 8 M-point tile.
* **The quick fit is not the answer, and must not be quoted as one.** Running the pair fit on
  the per-cell median of raw `z` (no CSF, no gen2) returns coefficients of ±30 to ±1850 mm per
  unit tangent, and the three independent estimates of `c_10010` disagree by a factor of four
  (−295.7, −74.3, −120.3 on open cells). The same crude estimator returns +383.1 for the
  136-137 pair where the pipeline's own reduction gives +169.7 and +184.6, so it is the
  estimator that is failing, not the geometry. **The geometry result above does not depend on
  any of it.**
* **Latitude.** The cross line sits at N 4,888,888–4,889,709, which is 2.6–3.5 km north of the
  elbaext grid's top edge (N 4,886,250). Whether `c_s` is constant along track is untested;
  the cross-line values against the Elba pair sums are the test.
* **More of it exists.** Its northing is cut by the tile top and its scan angles run −17…+2,
  so its nadir track and entire `+` side lie in `4342-27-64`; its easternmost return is 33 m
  from the tile's east edge, so it continues into `4358-28-01` over lines 139–140 and
  probably `4358-28-02` over 141–142. At the measured mean of 25.5 MB per tile
  (`ELBAEXT2_SCOPE.md` §6) that is ~25 MB each. Whether further cross lines (10011, …) exist
  cannot be answered from headers – `point_source_id` is not in the LAS header – so it costs
  one tile download per probe.

**Verdict: breaks the per-line degeneracy fully, at zero acquisition cost.**

---

## 3. Candidate 2 – ground control. **Does not reach `c_s`; it is the instrument for the tilt.**

Where each mark sits in its line's swath, from `ELBAEXT2_SCOPE.md` §2's measured off-nadir
distances and h = 2562 m, with the ties re-read from
`data/derived/groundtruth/elba_gen1_ties.json`:

| mark | line | easting | off-nadir | tan θ at the mark | tie (mm) | σ (mm) | move from c = 130 |
|---|---|---|---|---|---|---|---|
| 2210 | 128 | 570,492 | 147 m | **0.0574** | +21.3 | 12.4 | **+7.5 mm** |
| 3056 | 128 | 570,474 | 131 m | **0.0511** | −103.2 | 52.3 | +6.6 mm |
| 2024 | 129 | 571,244 | 32 m | **0.0125** | +156.6 | 54.5 | +1.6 mm |
| 2036 | 144 | 585,982 | 227 m | **0.0886** | +28.9 | 27.0 | **+11.5 mm** |

**Every mark sits essentially at nadir.** A per-line coefficient of 130 mm per unit tangent –
the representative value of `SWATH_ACROSS_TRACK_TEST.md` §2 – moves the tie by 1.6 to 11.5 mm,
against the marks' own σ of 12.4 to 54.5 mm. **No mark measures a `c_s`, and no line has two.**
The geometry that would be needed is two marks on the **same** line at **opposite ends of its
swath** – |off-nadir| ≳ 500 m each, on opposite sides – giving a tan lever arm of ~0.4 and a
`c_s` to ±(σ_mark·√2)/0.4 ≈ ±150 mm per unit tangent for a 40.8 mm mark. Even that is worse
than one cross pair. `ADDITIONAL_GROUND_CONTROL.md` reports 16 NVA+LCP points under a line of
the Elba network spanning −23.7 to +20.3 km in *northing*; its own ranking is siting, then
spread, then count, and it says nothing about across-track position because that is not what
those marks were sited for. **Marks will not supply a swath-edge pair.**

**But they are the right instrument for the tilt, and they already bound it.** The two anchors
are 15.49 km apart in easting. Their agreement converts directly into the mosaic tilt:

* **g = −0.485 ± 1.916 mm/km**, hence **v = g·h = −1.24 ± 4.91 mm per unit tangent**;
  2σ bound |v| ≤ **11.1**.
* Propagating instead the 42.6 mm *unmodelled* bound that `ABSOLUTE_BASIS_ELBA.md` deliberately
  keeps out of σ (validity of the lateral shift 7–16 km out) gives **σ(v) = 10.0**.

Against per-pair coefficients of **+34 to +193** mm per unit tangent, that is a tight bound:
the tilt component of the `c_s` is at most ~20 mm per unit tangent, ten per cent of the signal.
**The existing two anchors have already very nearly closed the part of the problem that the
cross line cannot touch.**

One tension to record rather than smooth over: the stable-ground DoD tilt at Elba is
**dE −14.19 ± 5.15 mm/km** (`ADDITIONAL_GROUND_CONTROL.md`), which as a gen1 mosaic tilt would
be −36.4 ± 13.2 mm per unit tangent and is **not** consistent with −1.24 ± 4.91. The two
measurements have different baselines (4.5 km within the tile against 15.5 km between marks)
and different epochs in play, so this says the Elba DoD tilt is either gen2's or local, not
that either measurement is wrong. It is worth settling.

**Verdict: does not break the per-line degeneracy at all; already bounds the tilt to
±5–10 mm per unit tangent.**

---

## 4. Candidate 3 – non-adjacent overlap. **None exists, anywhere.**

Measured across all thirteen tiles, from fitted nadir tracks and 99.5th-percentile half-widths:

| quantity | value |
|---|---|
| adjacent nadir-track spacing, 18 pairs | 848–1,060 m, median 974 m |
| swath half-width, 19 lines | 655–748 m, median 719 m |
| adjacent sidelap | +374 to +607 m (every pair overlaps) |
| **second-neighbour separation** | **1,822–1,974 m** |
| **second-neighbour swath sum** | **1,382–1,476 m** |
| **second-neighbour sidelap** | **−380 to −559 m** |

A second neighbour would need the line spacing to fall below ~719 m; the smallest observed is
848 m, between lines 143 and 144. **No second-neighbour overlap exists on any of the eighteen
adjacent pairs from line 128 to line 145**, confirming `ELBAEXT2_SCOPE.md` §3's "all zero" on
a wider footprint and with the numbers attached.

And even if one existed it would not help: §1 shows the null direction is a tilt, which cancels
in every between-line difference regardless of which two lines are differenced. **A loop in the
N-S network cannot break it.** This candidate is closed twice over.

**Verdict: does not exist, and would not break it if it did.**

---

## 5. Candidate 4 – terrain relief within a swath. **Real, but a factor of 3–5 short, and confounded.**

The degeneracy is exact only if `S = tan θ_A + tan θ_B` is exactly constant. It is not:
`S = spacing/(H − z)`, so relief modulates it. That prediction holds, on the Elba data, read
through the pipeline's own reduction (`swath_across_track_test.cell_swath_ground`, so this is
the identical population §1 of that document used):

| tile | pair | cells | sd(sum) | sd(diff) | p | SE(p) | **q** | **SE(q)** | SE(q)/SE(p) |
|---|---|---|---|---|---|---|---|---|---|
| elba | 135-136 | 7,394 | 0.0206 | 0.0675 | +124.2 | 17.9 | −19.1 | 52.7 | 2.9 |
| elba | 136-137 | 9,817 | 0.0202 | 0.1062 | +169.7 | 8.0 | +25.2 | 30.1 | 3.8 |
| elba | 137-138 | 10,411 | 0.0182 | 0.1117 | −41.1 | 21.2 | **+138.5** | 92.5 | 4.4 |
| elbaext | 133-134 | 7,218 | 0.0244 | 0.0637 | +185.2 | 17.1 | +179.7 | 45.9 | 2.7 |
| elbaext | 134-135 | 14,751 | 0.0185 | 0.1002 | +62.7 | 6.6 | −18.3 | 33.1 | 5.0 |
| elbaext | 135-136 | 14,582 | 0.0194 | 0.1049 | +131.1 | 9.3 | −114.7 | 42.0 | 4.5 |
| elbaext | 136-137 | 11,534 | 0.0198 | 0.1040 | +184.6 | 9.1 | +97.9 | 36.0 | 3.9 |
| elbaext | 137-138 | 12,030 | 0.0174 | 0.1129 | −36.1 | 17.9 | **−215.4** | 83.2 | 4.6 |

The leverage is **real geometry, not quantisation.** Regressing `sum_tan` on the gen2 ground
elevation of the same cell gives the predicted sign in **8 cases of 8** and the predicted
magnitude to within a factor of 2.4:

| tile | pair | d(sum_tan)/dz observed | predicted `sum_tan/(h−z)` | r² |
|---|---|---|---|---|
| elba | 135-136 | −3.782e-04 | −1.603e-04 | 0.431 |
| elba | 136-137 | +6.932e-05 | +1.669e-04 | 0.025 |
| elba | 137-138 | −3.263e-04 | −1.555e-04 | 0.615 |
| elbaext | 133-134 | −2.418e-04 | −1.632e-04 | 0.164 |
| elbaext | 134-135 | +3.271e-04 | +1.721e-04 | 0.327 |
| elbaext | 135-136 | −2.526e-04 | −1.607e-04 | 0.278 |
| elbaext | 136-137 | +4.920e-05 | +1.671e-04 | 0.015 |
| elbaext | 137-138 | −3.209e-04 | −1.549e-04 | 0.585 |

So the degeneracy is lifted **partially** – `q` is estimable, at 2.7–5.0× the standard error of
`p`, i.e. **±30 to ±93 mm per unit tangent**. Two things say not to trust it:

1. **It does not reproduce.** The two tiles estimate three pairs independently. For 137-138 they
   give **+138.5 ± 92.5 and −215.4 ± 83.2** – opposite signs, 2.8σ apart. The symmetric
   coefficient reproduces to <1.2σ on the same pairs (`SWATH_ACROSS_TRACK_TEST.md` §2).
2. **The regressor is terrain elevation.** With r² up to 0.615, fitting on `sum_tan` is very
   nearly fitting the between-line difference on `z`. Any elevation-correlated between-line
   difference – cover, valley against upland, along-track drift sampled differently by the two
   lines – lands on `q` directly. That is the same failure mode the along-track control exposed
   on the two 137-138 pairs in the original test, and 137-138 is exactly where `q` breaks.

**Verdict: lifts it partially – `q` to ±30–93 mm per unit tangent – but the estimate does not
reproduce across extents and its regressor is a terrain proxy. Do not build a correction on it.**

---

## 6. Candidate 5 – gen2 as an external reference. **Confounded with its own tilt, by construction.**

Within a flight line, easting and scan angle are the same variable
(`x = x_track + h·tan θ`, r² = 0.95–0.997, `SWATH_ACROSS_TRACK_TEST.md` §10). So a cross-track
tilt in **gen2** maps one-to-one onto the per-line coefficient it would be used to measure,
with coefficient `g·h_s` – **alternating in sign with the line's body-fixed side**, which is
precisely the null direction of §1. gen2 does not break the degeneracy; it substitutes gen2's
tilt for gen1's.

The size is not hypothetical. The measured stable-ground DoD tilt of **−14.19 ± 5.15 mm/km**
is **36.4 ± 13.2 mm per unit tangent** in the null direction – three times the 11.1 mm 2σ
bound the ground-control anchors give, and comparable to the spread among the pair
coefficients themselves.

What gen2 *would* resolve: each line's across-track ramp **relative to gen2**, one line at a
time, with no pairing and no chain. `SWATH_ACROSS_TRACK_TEST.md` §5 already ran that regression
and found it systematically *under*-reads: reading `(c_A + c_B)/2` off the gen2-referenced
within-swath slopes disagrees with the overlap value by −4.09 to +1.69 mm/deg, negative on five
of eight pairs. What it would not resolve: anything absolute. It makes gen1's correction a
function of gen2, which is acceptable for a DoD and **not** acceptable for an
absolutely-referenced gen1 DEM – and `analysis/ABSOLUTE_BASIS_ELBA.md` §4 has now built exactly
such a product.

**Verdict: does not break it. It relabels gen1's unknown tilt as gen2's, and the relabelled
quantity is three times larger than the bound ground control already provides.**

---

## 7. Candidate 6 (not on the list; found in the data) – a level water surface

A gravity-level surface is the one reference that sees a tilt. The tiles hold **562,664 vendor
class-9 (water) returns**, and **458,239 of them are in one tile**, `4358-26-03`, at a modal
elevation of **201.02 m** over E 584,890–587,373, N 4,893,254–4,896,747 – a 2.5 by 3.5 km
sheet of water at 44.203°N, 91.922°W, roughly 12 km north-east of the Elba reference point and
150 m below it. Its size, its flatness and its position in the Mississippi floodplain make it
a backwater of a navigation pool rather than a channel, *which is an inference from the
coordinates and the elevation, not a verified identification, and the level-surface premise is
exactly the thing that has to be checked before any number is used.* Three lines cross it
(144 and 145 over most of their scan range, 146 one-sided at the pool's edge), and the returns
are real rather than hydro-flattened (z sd 53–62 mm, not zero):

| psid | n (|Δz| < 0.15 m) | 50 m blocks | easting span | tan range | z sd | c (mm/unit tan) | SE | implied tilt |
|---|---|---|---|---|---|---|---|---|
| 144 | 96,558 | 158 | 1,078 m | −0.194…+0.231 | 0.062 m | **−125.8** | 14.9 | −49.1 mm/km |
| 145 | 168,670 | 391 | 1,330 m | −0.249…+0.231 | 0.053 m | **−133.9** | 8.4 | +52.3 mm/km |
| 146 | 120,192 | 232 | 541 m | −0.268…−0.052 | 0.060 m | +12.7 | 20.3 | – (one-sided; no nadir returns in the tile, so its heading sign is unfitted) |

**Read this carefully, because the obvious reading is wrong.** Lines 144 and 145 fly opposite
directions, so a real tilt of the water would give implied tilts of the *same* sign. It gives
**opposite** signs of nearly equal magnitude, which means the two lines agree in the **body**
frame – and that is exactly what a common water-return response to view angle (specular and
wave-facet geometry, identical on every line) would produce. **The absolute −130 must not be
quoted as `c_144` or `c_145`.**

What survives is the difference, in which any common water response cancels:

| Δz window | c_144 | c_145 | **q = (c_A − c_B)/2** | SE |
|---|---|---|---|---|
| 0.15 m | −125.8 | −133.9 | **+4.1** | 8.6 |
| 0.30 m | −193.5 | −152.6 | −20.5 | 15.5 |
| 0.60 m | −124.4 | −34.0 | −45.2 | 41.9 |

At the tightest window this is `q` to **±8.6 mm per unit tangent**, three to ten times better
than the N-S overlaps deliver (§5). Two things stand between that and a usable number, and
both are stated rather than hoped away:

* **It is window-sensitive.** The `Δz` window is **MINE**; over 0.15/0.30/0.60 m the point
  estimate moves from +4.1 to −45.2. A defensible run needs a real water mask, not a window.
* **The pool's own slope re-enters with the same alternating structure.** Each 1 mm/km of true
  easting slope adds **2.56 mm per unit tangent** to `q`, so a pool level to 10 mm/km leaves a
  26 mm systematic beside the quoted 8.6. A backwater 10 km above the dam is close to still,
  but that is a hydraulic argument, not a measurement.

A crude land-overlap control on the same tile gives +214.3 ± 30.6 for the 144-145 pair, which
disagrees with the water's −129.9. It does **not** adjudicate: it is the same raw-median
estimator that returns +383.1 for 136-137 where the pipeline gives +169.7 and +184.6. **Neither
reading is trustworthy without CSF.**

**Verdict: promising and cheap – the only gravity-referenced instrument in our data, and it
sits on the east chain (lines 144–146) that carries the 2036 anchor. Needs a proper water mask
and CSF before any number leaves it. Second priority, behind the cross line.**

---

## 8. Ranked recommendation

**Do this one thing: classify tile `4342-28-64` with the pipeline's own CSF and fit the
cross-line pairs 136-10010, 137-10010 and 138-10010.**

It is ranked first on every axis that matters:

1. **It is the only candidate that breaks the per-line degeneracy fully.** SE(q)/SE(p) of
   1.05–1.14 against 4.95–7.09; 55–110 ha of shared ground per pair; line 137 seen over its
   full across-track range.
2. **It costs one CSF pass on a tile already on disk.** No download, no gen2, no chain.
3. **It arrives with its own check.** Three independent estimates of `c_10010`, and three
   cross-line values of `c_136`, `c_137`, `c_138` to test against the three N-S pair sums the
   Elba tiles already measure. The network has had zero redundancy until now.
4. **It answers the along-track question at the same time.** If the cross-line `c_s` reproduce
   the Elba pair sums 3 km south, `c_s` is a line constant; if not, that is the finding.

Then, in order: **(2)** fetch `4342-27-64` and `4358-28-01` (~25 MB each) to complete the cross
line's own swath and carry it east over lines 139–142; **(3)** do the water surface properly on
`4358-26-03`, with a real water mask and a check on the pool's own slope, for an independent
read on the east chain; **(4)** leave
the tilt to ground control, where the two existing anchors already give
**g = −0.485 ± 1.916 mm/km** and where `ADDITIONAL_GROUND_CONTROL.md`'s marks are the way to
sharpen it; and **(5)** settle the −14.19 mm/km DoD tilt against that bound, since the two do
not presently agree.

**What is not worth doing:** looking for a second-neighbour overlap (none exists, and it would
not help), looking for a swath-edge pair of survey marks (the marks are all within 0.09 in
tangent of nadir and none sits on the same line twice), or using gen2 to fix the per-line
ramps (it substitutes a 36 mm per unit tangent confound for an 11 mm one).

---

## 9. The honest negative, stated plainly

**The global cross-track tilt of the gen1 mosaic is not recoverable from lidar geometry at
all** – not from adjacent overlaps, not from non-adjacent overlaps, not from cross lines, and
not from gen2 without inheriting gen2's own. §1 proves it exactly: the tilt is common to both
lines at a shared point, so it cancels in every between-line difference that can be formed.
Only an absolute external reference reaches it, and we have exactly two – the surveyed anchors
2210 and 2036, 15.49 km apart – which bound it to ±5 mm per unit tangent formally and ±10 mm
with the unmodelled term carried.

That is the part of the original claim that stands. The part that does not is the framing: the
per-line coefficients themselves **are** recoverable, from a cross line that has been sitting in
`data/before/4342-28-64.laz` since 2026-08-15.
