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

## Surface-roughness ladder (cell roughness, m)

| surface | region (UTM15N box) | gen1 floor / med | gen2 floor / med | gen1/gen2 (med) |
|---|---|---|---|---|
| infield (skinned) | (498980, 499230, 4976180, 4976400) | 0.017 / 0.019 | 0.011 / 0.012 | 1.65 |
| parking lot       | (498870, 499160, 4976090, 4976145) | 0.023 / 0.025 | 0.012 / 0.014 | 1.80 |
| mowed grass       | (498980, 499230, 4976180, 4976400) | 0.020 / 0.023 | 0.012 / 0.015 | 1.54 |

- **gen2 floor is material-independent** (0.011–0.012 across all three) → it is the
  sensor. gen1 floor is 0.017–0.023 (a bit more surface-dependent; infield cleanest).
- **All three are sensor-dominated:** the gen1/gen2 median ratio is ~1.5–1.8 on
  every surface, ≈ the sensor ratio (0.017/0.011 ≈ 1.6). These surfaces are too
  smooth for real surface roughness to rise above the ranging noise.

## What does NOT separate yet, and why the ladder must go rougher

Removing each epoch's floor in quadrature to recover surface roughness
`S = √(r² − σ²)` is **not reliable at these smooth scales**: on the median cell it
gives an infield ratio of 2.3 — *larger* than grass's 1.5 — even though the infield
is hard ground. Near the floor, S is dominated by floor-subtraction sensitivity and
the sparser gen1's per-cell estimation noise, not by real surface differences. So
mowed grass shows **no** clean vegetation-specific divergence; it is just another
smooth, sensor-dominated surface.

The vegetation divergence is real — under forest canopy roughness broke down
entirely (roughness failed to predict DoD error; see the roughness-covariate work).
But it needs a surface where real roughness clears the ranging noise. Mowed turf
does not; **prairie grass and, above it, brush/forest are the rungs where the gen1
over-read (sparse 2008 mixing canopy returns into "ground") should emerge.** That
ladder — smooth (sensor-limited) → prairie → brush → forest (breakdown) — is the
open characterization thread.
