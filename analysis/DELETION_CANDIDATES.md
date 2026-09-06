# Deletion candidates — a report. Nothing here has been deleted.

Phase 3 of `analysis/REORGANIZATION_PLAN.md`, and it is deliberately not autonomous:
deleting a tracked record is Andy's call, and a classifier I measured at 17% wrong on a
hand-checked sample does not get to make it.

**The rule**: no workflow Step names it, no `TOOLS` entry declares it, no other script
imports it, and no document cites it — by path or by bare filename. 26 of 157.

## The list, dated by its last SUBSTANTIVE commit

Not by last commit: today's mechanical passes (`anchor:` and `reorg:`) touched 38 of these
files and would have made half of them look like current work. Those commits are skipped.

    2026-08-21  analysis/ridgelines/curvature_cdf.py
    2026-08-21  analysis/slope_bias/coarsen_resolution.py
    2026-08-21  analysis/slope_bias/ground_class_structure.py
    2026-08-21  analysis/slope_bias/ground_return_stats.py
    2026-08-21  analysis/slope_bias/recoverability_test.py
    2026-08-21  analysis/slope_bias/understory_from_lidar.py
    2026-08-23  analysis/ridgelines/elbaext_tie_offsets.py
    2026-08-23  analysis/slope_bias/andy_spot.py
    2026-08-23  analysis/slope_bias/blufftop_margin_forest.py
    2026-08-23  analysis/slope_bias/cliff_cells.py
    2026-08-23  analysis/slope_bias/flat_forest_offset.py
    2026-08-24  analysis/ridgelines/boresight_lateral_coupling.py
    2026-08-25  analysis/ridgelines/forest_on_flat_map.py
    2026-08-26  analysis/ridgelines/gen1_csf_pdal_defaults.py
    2026-09-01  analysis/geoid_plane_vs_extent.py
    2026-09-01  analysis/provenance_of.py
    2026-09-03  analysis/control_mixture_validation.py
    2026-09-03  analysis/control_percentile_fit.py
    2026-09-04  analysis/aspect_audit.py
    2026-09-04  analysis/control_mat_fit_nowindow.py
    2026-09-04  analysis/control_q_in_ground.py
    2026-09-04  analysis/control_shift_predictors.py
    2026-09-04  analysis/gen2_coverage_audit.py
    2026-09-04  analysis/plot_correction_method.py
    2026-09-04  analysis/steady_state/lod_on_steady_cells.py
    2026-09-05  analysis/swath_network_weighting.py

## Where the rule is WRONG, and I can name the cases

**Do not delete the last eleven without reading this.** Everything dated 2026-09-01 or later
was written during the current work, several at Andy's explicit request, and the rule misses
them because **a citation in a COMMIT MESSAGE is not a citation in a document**:

* `swath_network_weighting.py` — written 2026-09-05 to reproduce the variance-weighting
  result. Commit 9e78f4a says "reproduces every number above". It is the evidence for a
  change to the shipped swath network.
* `control_q_in_ground.py`, `control_mat_fit_nowindow.py` — the per-mark
  percentile-vs-spread producers, RECOVERED on 2026-09-04 precisely because they had gone
  missing once already.
* `plot_correction_method.py`, `gen2_coverage_audit.py`, `aspect_audit.py` — figure and
  audit producers written when Andy asked to see the correction and the coverage seam.
* `provenance_of.py` — trust tooling.

The pattern is the one this whole exercise keeps rediscovering: **evidence of use lives in
places a static scan does not read.** First it was scripts that compute a quantity instead
of loading it; then documents that name a file without its path; now work recorded in git
history rather than in a document.

## What I would actually propose

1. **The 14 dated 2026-08-2x are the real candidates.** They are the slope-bias and
   beam-geometry investigations, and their premise was refuted or superseded — the same
   class as the 33 penetration scripts already removed.
2. **Leave everything from 2026-09-01 on.** It is current work whose only record is a
   commit message.
3. **Better than deleting: cite them.** A script named in the findings document it produced
   stops being a candidate, and the citation lint then keeps that link honest. That converts
   this list from "what can we delete" into "what was never written up" — which is the more
   useful question, and the one the 53%-unclaimed figure has been pointing at all along.
