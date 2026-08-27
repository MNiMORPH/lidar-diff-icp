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
sign to, the leaf-on canopy bias that lifts gen2 and reads as false deposition
(`analysis/ngv.py`, `analysis/ngv_correct_dod.py`). Two epochs, two mechanisms, opposite
signs -- which is the likeliest reason a single vegetation correction overshot at Elba.

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
