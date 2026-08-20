# Roughness characterization at Battle Creek

**Battle Creek is a characterization site, not a study site.** It was chosen for
its clean, ground-truthable surface types — sealed pavement, skinned softball
infields, mowed turf, prairie — which natural study sites lack. The purpose is to
separate, per epoch, the **instrumental ranging noise** from **real surface
roughness**, and to find where the two lidar generations agree versus diverge.

- **gen1** = 2008 Minnesota (Metro) lidar
- **gen2** = 2021 USGS 3DEP lidar

**Metric:** NMAD of the residuals to a per-cell least-squares plane (detrended
within-cell roughness). NMAD, not RMS — RMS is inflated by blunder returns near
vegetation. Reported as the **p20 "floor"** (smoothest cells ≈ the surface's clean
value) and the **median** (typical cell). Surface roughness `S` is recovered by
removing each epoch's instrumental floor in quadrature: `S = √(r² − σ²)`.

Frame (EPSG:26915): `498135 4975136 499365 4976876` (4× the pilot tile). Code:
`surface_roughness.py`. Data provenance and reacquisition commands are in that
script's header.

## Instrumental precision (sensor ranging noise)

Measured as the roughness floor on near-zero-roughness built surfaces:

| epoch | σ (instrumental floor) |
|---|---|
| **gen2 (3DEP 2021)** | **≈ 0.011 m** |
| **gen1 (2008)** | **≈ 0.017 m** |

These are **softball-infield floors = instrument ⊕ a few mm of surface irregularity**
(rake/grain/grading texture), so they are an *upper bound* on pure ranging noise.
gen1's 2008 sensor is ~1.6× noisier than modern 3DEP.

**Material-independence check** (this is what proves it is the *sensor*): the floor
matches across sealed pavement and clay/sand infield to within ~1 mm (gen2) / a few
mm (gen1), and after removing the floor the two epochs *agree*. gen2 σ = 0.011 m
reproduced five ways (pilot frame, 4× frame, pavement, infield, NDVI-built subsets).
Note the infield is *smoother* than the parking lot — a lot carries a drainage
crown, expansion joints, and painted stripes — so the ballfield, not the lot, is
the cleanest reference.

**Literature cannot supply this number.** Published specs are *absolute* accuracy
vs. GPS checkpoints (3DEP QL1/QL2 = 10 cm RMSEz NVA; MN 2008 ≈ 5 cm open, MnGeo
CVA) — a coarser quantity that includes georeferencing and systematic error at
~10 cm scale. Our *within-cell relative* precision (1.1 / 1.7 cm) sits well inside
the 3DEP within-swath relative ceiling. The built-surface measurement was
necessary; literature only bounds it.

## Surface-roughness ladder (median cell roughness, m; MODERN processing only)

**gen1 = CSF ground, gen2 = class-2 ground.** (Last-return is retired — it reads
canopy, not ground, in tall vegetation; see the historical note below.)

| surface | region (UTM15N box) | gen1 CSF | gen2 | gen1/gen2 |
|---|---|---|---|---|
| infield (skinned) | (498980, 499230, 4976180, 4976400) | 0.019 | 0.012 | 1.67 |
| parking lot       | (498870, 499160, 4976090, 4976145) | 0.025 | 0.014 | 1.81 |
| mowed grass       | (498980, 499230, 4976180, 4976400) | 0.022 | 0.015 | 1.52 |
| prairie grass     | (498780, 499180, 4975680, 4976030) | 0.063 | 0.033 | 1.94 |
| oak forest        | (498150, 498850, 4975200, 4976050) | 0.067 | 0.043 | 1.58 |

- **The epoch ratio is ~1.6–1.9 across the WHOLE ladder** — pavement to oak forest —
  ≈ the sensor ratio (σ₁/σ₂ = 0.017/0.011 ≈ 1.6). With modern (CSF) processing there
  is **no vegetation divergence**: both surveys, once on true ground, are
  sensor-limited everywhere.
- **What changes up the ladder is the ABSOLUTE roughness**, both epochs tracking it:
  gen2 rises 0.012 (hard) → 0.033 (prairie) → 0.043 (forest) — real surface roughness
  (microtopography, low-vegetation-supported ground) climbing above the ~0.011 m
  sensor floor, while the ratio stays pinned at ~1.7×.
- **All-returns check** (infield/parking): roughness from ALL returns equals the
  classified-ground value to ~1 mm (infield 0.019/0.012, parking 0.025/0.015),
  despite ~3× the points — on hard surfaces every return *is* the surface, so the
  floor is genuinely the sensor, not a classification artifact.

**Historical note (why last-return is retired):** with obsolete gen1 last-return,
prairie showed a spurious 8.7× ratio (gen1 "roughness" 0.283 m) — last-return grabs
grass tops, not ground. CSF gen1 pulls that to 0.063 m (ratio 1.94), in line with the
rest. The apparent "vegetation divergence" was entirely a ground-definition artifact
of an obsolete method.

## A caveat on floor-removed surface roughness

Recovering surface roughness by removing the floor in quadrature, `S = √(r² − σ²)`,
is **not reliable at the smooth end**: on the median cell it gives an infield ratio
of ~2.3 — larger than grass — even though the infield is hard ground. Near the floor,
S is dominated by floor-subtraction sensitivity and the sparser gen1's per-cell
estimation noise. Use the raw ratio (~1.7 throughout), not floor-removed S, at these
scales; S only becomes meaningful where roughness clears the floor (prairie, forest).

## Oak forest — extra error-structure detail (`oak_forest.py`)

The forest is the top rung of the ladder above (ratio 1.58, sensor-limited like the
rest, once gen1 uses CSF). `oak_forest.py` adds the detail that only the forest needs:

- **Penetration is comparable:** gen1 (CSF) median **167** ground returns/cell,
  gen2 **148**. Deciduous leaf-off oak lets both surveys reach ground; gen1 is not
  starved of ground here.
- **Error is near-white:** gen2 ground-roughness autocorrelation drops 0.48 (5 m) →
  0.26 (15 m), 1/e length ~7 m, with a weak ~0.25 canopy-scale plateau. Largely
  cell-independent, mild dependence on density (gen2 roughness↔count Spearman −0.24:
  rougher where sparser).

**Regime caveat:** this is deciduous **leaf-off** oak, where both surveys penetrate.
The forest roughness *breakdown* seen at Cook (roughness meaningless, DoD error
un-localizable) is a **different regime** — conifer, leaf-on, poor penetration — not
a universal property of "forest".

## Epoch roughness structure: gen1 ≈ 1.7× gen2, multiplicatively

Across the whole modern ladder, gen1 (CSF) roughness is ~1.5–1.9× gen2's, with
comparable or higher gen1 ground-return density:

| surface | ratio | gen1 n | gen2 n | slope |
|---|---|---|---|---|
| infield | 1.67 | 234 | 387 | 0.4° |
| parking | 1.81 | 248 | 383 | 0.9° |
| grass | 1.52 | 260 | 389 | 0.7° |
| prairie | 1.94 | 234 | 222 | 1.9° |
| oak forest | 1.58 | 167 | 148 | 7.3° |

**It is the instrument, not density or sampling.** In the oak forest gen1 has *more*
ground returns than gen2 (167 vs 148) yet *higher* roughness — at equal-or-higher
density gen1 is still rougher, so undersampling is excluded. Prairie is not
pulse-starved either (gen1 234/cell = the hard surfaces), so its marginally-top ratio
is not a density effect.

**But it is not a fixed additive noise floor.** A constant ranging-noise σ added in
quadrature (`r² = S² + σ²`) would make the ratio fall toward 1.0 where real roughness
is large — yet it *stays* ~1.6–1.9 at 3–7 cm prairie/forest roughness. So the gen1
excess **scales with the roughness being measured** → points to the 2008 sensor's
larger **footprint / beam divergence** (a bigger laser spot smears a rough surface
proportionally), not to ranging precision. Still "the instrument" — the footprint,
not the noise floor.

**Decimation test — CONFIRMED (density formally excluded).** Thinning gen2 to gen1's
exact per-cell count does not move gen2's roughness at all: infield 0.0117→0.0118,
parking 0.0139→0.0139, prairie 0.0326→0.0325, oak forest 0.0426→0.0423. At *identical*
density gen2 stays ~1.7× smoother than gen1, so the gap is purely the 2008 sensor
(footprint), not sampling — airtight.
