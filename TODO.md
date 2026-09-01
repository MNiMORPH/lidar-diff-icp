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

## Generalization: what stops a second tile running end to end

**State as of 2026-09-01, after the two producers were parameterized.** elbaext now has
`floodplain_mask.npy` and `crest_mask.npy`; it is still missing `penetration.npy`,
`canopy_struct.npz`, `core_forest.npy` and `core_open.npy`.

### 1. DONE: both producers are parameterized by tile

`convexity_dod_landcover.py` (commit `94a3833`) and `strata_core.py` (commit `3f051d3`).
Grid origin and resolution are read from each tile's own `corrections.json`; every output
path is built from `--tile`. Both reproduce elba **byte-identically** (`floodplain_mask`,
`crest_mask`, `kappa_L10/20/30`, `core_forest`, `core_open`; `ridgecrest_pixels.npz` matches
in column order and in every column but `dod_m`, which differs because
`elba_refdatum/dod_geoid.npy` was regenerated after that table was written).

**The blocking gate is now clear.** `reference_cells` runs on both tiles: 28,533 cells at
elba_fulldensity, 78,270 at elbaext.

Two things this turned up, both recorded and neither acted on:

- **elbaext's DoD was chosen, not given.** The elbaext run read
  `data/derived/elbaext/dod_geoid.npy` -- the only `dod_geoid` on that tile and the
  name-for-name analogue of elba's `elba_refdatum/dod_geoid.npy`. `--dod` is required with
  no default precisely so this choice is visible. **Confirm it is the intended grid.**
- **`curvature_diffusion.py` has the same defect**, out of scope here: it takes a tile
  argument but hardcodes `data/derived/elba_refdatum/dod_geoid.npy` in two places, and it
  `sys.exit(0)`s silently when `penetration.npy` is absent. It also AUGMENTS
  `ridgecrest_pixels.npz` in place with `curv_xx`/`curv_yy`/`curv_laplacian`, so it must be
  re-run after the convexity producer or those three columns are silently dropped.

### 2. STILL OPEN: retire `penetration.npy`, or give it a producer

`penetration.npy` has **no producer anywhere in the repo** — every write form was searched,
twice. It exists only for elba, dated 2026-08-21. `src/lidar_diff_icp/canopy.py:19`
`ground_penetration()` computes exactly this quantity from the gen2 cloud, so writing a
producer is a short script, not research; the objection to it is scientific, not effort.

**`canopy_struct.npz` is a second orphan of the same kind** — no producer either, read by six
files, present only for elba. `strata_core.py` needs BOTH, so it is blocked twice over.

What is still gated on the decision: the forest/open crest split (step 4 of the convexity
producer) and `core_forest`/`core_open`. Nothing else.

**Recommendation: retire it.** `analysis/ridgelines/AUDIT_findings.md` flags `penetration.npy`
as a gen2-derived variable contaminating gen1-internal conclusions; it is referenced by 34
tracked files. `canopy_cover_pfs.npy`, `forest_pfs.npy` and `open_pfs.npy` already exist for
BOTH tiles and are the cover measure meant to replace it. The cost is not free: swapping the
crest split to the PFS layers CHANGES elba's step-4 answer, so it needs the two cover
definitions measured against each other on elba first, not substituted. Andy's call.

### 3. Then produce elbaext's cover-corrected products

`dod_cover_corrected.py` and `lod_cover_q2.py` are ALREADY tile-parameterized; they have
simply never been run there. Two commands. The q2 slope must be elbaext's OWN fit, because
the relation is per-site (it depends on each pair's phenology and undergrowth).

### 4. Loose end, not blocking

`nearground_cells_sn.npz`, `gen1_canopy_frac.npz` and `ridge_mask.npy` showed no producer in
the scan, but the scan cannot see filenames built from flags (`nearground_cells.py --out
..._sn.npz` is invisible to it). All exist for BOTH tiles, so they do not block; the
producers just are not identified. `curv_laplacian.npy` IS identified —
`curvature_diffusion.py` writes it, along with `curv_xx` and `curv_yy`.

### The rule these now follow

A missing input REFUSES; it does not run differently. Running without an optional layer is
allowed and often right, but must be stated -- `use_floodplain_mask=False`,
`--without penetration,core_forest`. Silent zero-fill turned an ABSENT layer into a
MEASURED EMPTY one, which is how a table came to look like a finding of "no forest", and how
two tiles' reference populations came to differ by 39,038 cells unnoticed.

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
