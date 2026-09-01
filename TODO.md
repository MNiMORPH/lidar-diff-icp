# Open items

Each entry states the problem, what has been measured, and what remains.
Numbers here are pasted from a run, with the script that produced them.

---

## gen1 may not resolve ground under dense grass

**The problem.** gen1's per-cell ground elevation is the MEDIAN of its CSF-classified
ground returns. Under a dense sward, returns that stop in the grass enter that class and
the median floats above true ground. gen2 has roughly 12x the ground-return density
(5.78 vs 0.49 pts/m^2), so far more of its shots reach the real surface through the same
grass. The DoD is `gen2 - gen1`, so the result reads as **false erosion**.

This is a SAMPLING argument, not a leaf-state one. It is independent of, and opposite in
sign to, the leaf-on canopy bias that lifts gen2 and reads as false deposition, which the
DoD already corrects via `q2(cover)` (`analysis/ridgelines/q2_cover_fit.py`). Two epochs,
two mechanisms, opposite signs.

**What is measured** (`analysis/gen1_grass_lift.py`, commit `6e41db1`; gen1-INTERNAL,
`p50 - p10` of each cell's own gen1 ground returns from `gen1_csf_angles.npz`):

```
  group                        ncell   p50-p10   slope   cover    DoD mm
  floodplain, cover<0.15       38894      96.0     5.4   0.010     +18.8
  upland,     cover<0.15       99555      60.2     3.8   0.000      -5.2

  floodplain EROSION cells  n=12,789  lift 220.0 mm  slope 10.0  cover 0.170  DoD -188.6
  floodplain, not eroding  n=114,724  lift 106.7 mm  slope 13.8  cover 0.230  DoD  +28.5
```

Open floodplain against open upland: `96.0` vs `60.2 mm`, at near-zero canopy in both and
similar low slope -- meadow sward against cropland. The cells detected as erosion carry
twice the lift while being FLATTER and LESS canopied, so neither intra-cell relief nor
leaf-on explains them; both would push the other way. The lift, `220 mm`, is about the
size of the apparent erosion, `-189 mm`.

Independent support: Google Earth shows the eroding floodplain as open meadow in
September 2008 and a dense sapling thicket in April 2021.

Per-cell lift grid: `data/derived/elba_fulldensity/gen1_lift_p50_p10.npy`.

**Two ways forward, and the choice is open.**

1. *Correct gen1.* Replace the per-cell median with a lower statistic. **`p10` is not
   proposed** -- it was used to measure the spread and is an unjustified estimator.
   A defensible replacement needs either calibration against surveyed 2008 control on
   grass, or a floor derived from return count and footprint. Whatever is chosen must be
   applied the same way to both epochs, or its effect on gen2 measured, or it manufactures
   a difference.

2. *Blank the affected cells.* Mask cells whose gen1 lift exceeds a threshold and report
   them as no-data rather than as change. Cheaper, honest, and loses coverage exactly
   where the grass is thickest. The threshold is still a choice that has to be justified.

**Caution before any tile-wide version.** The lift also grows with intra-cell relief. The
flat floodplain is a clean control precisely because relief is low there; a global
correction needs the slope dependence separated first.

**Still to do.** Confirm the high-lift cells coincide with the 2008 open meadow; choose
between correction and masking; re-run the floodplain DoD against whichever is chosen.

---

# Decisions

Closed questions, with the reasoning and the conditions that would reopen them.

## gen2's horizontal accuracy is not checked, and does not need to be

**Decided 2026-08-31 (Andy).** Closes next-action 3 of
`ground_control/HANDOFF_FROM_GROUND_CONTROL.md`, which listed gen2's unchecked horizontal
accuracy as an open item on the grounds that step 4 registers gen1 *to* gen2, so a lateral
error in gen2 propagates into gen1 invisibly and becomes a vertical error on a slope.

**Why it is closed, and it is not "3DEP is probably fine".** The DoD is invariant to it.
Registration applies a Nuth-Kaab lateral shift that moves gen1 onto gen2, so a UNIFORM
horizontal error in gen2 is absorbed: both epochs end up in the same frame and the
difference of two clouds that agree with each other does not depend on where that shared
frame sits. The quantity we compute cannot see the error.

**Where it does not cancel** is the tie to control, because the marks carry independently
surveyed horizontal coordinates, so a lateral error samples gen2 at the wrong place. The
geometry makes that cheap, because control is sited flat by design:

```
  slope at the 389 gen2 control marks (deg): p50=2.6  p90=6.5  p99=15.5

  vertical error produced by a lateral error, at those marks (mm)
     lateral      p50      p90      p99
       10 cm      4.5     11.4     27.8
       20 cm      9.0     22.9     55.6
       50 cm     22.6     57.2    139.0

  for scale: the gen1 Elba datum SE is +/- 23.38 mm
```

A 20 cm lateral error costs 9 mm at the median mark. It also averages down over 389 marks,
because a fixed lateral shift produces vertical errors whose sign follows local aspect, so
with mixed aspects it behaves as noise rather than as an offset.

**Two conditions that would reopen it.**

1. *Non-uniform lateral error.* The cancellation assumes one shift fits the whole tile. A
   per-swath or along-track lateral drift would not be absorbed, and this acquisition has
   already required a per-swath VERTICAL tie, so a lateral analogue is not far-fetched. It
   would appear as edge-of-swath artefacts, not as a whole-tile offset.
2. *Comparing gen2 to anything it has not been registered to.* The flatness that makes the
   error cheap holds at the marks, NOT on the tile: Elba slope is p50 9.2 deg, p90 28.6, so
   a 20 cm lateral error there is 32.5 mm at the median and 108.9 mm at p90. Those numbers
   are harmless only while the error cancels.

Numbers from `analysis/control_lowveg_offset.py` mark slopes and
`data/derived/elba_fulldensity/slope.npy`.
