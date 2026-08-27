# ★★ ground_control FRAME — 2026-08-27. READ FIRST.

Verify every structural claim here against git and the files before acting.

## Why this subsystem exists

The DoD is `gen2 − gen1`, and each surface sits at its own unknown height. **Coregistration
— relative alignment — works and was never touched here.** What did not work is the
**absolute datum**, and it never did: the shipped `data/derived/elba_fulldensity/
z_before_absolute.json` carries **+22.7 mm with its own ±39.7 mm**, an uncertainty larger
than the value, built by chaining two *2021* marks 5–6 flight lines onto a *2008* surface.

That was not a design decision. git shows the product was written 2026-08-26 15:53 and
gen1's own 2008 control entered the repo at **17:17** — 83 minutes later. It used the only
control that existed.

**The rule this subsystem enforces: 2008 control for gen1, 2021 control for gen2, no
cross-epoch control, no chaining where the site's own line carries marks.**

## ★★ CONCLUSION 2026-08-27 — the epoch difference is settled; the absolute is not

**Stop chasing the absolute datum at Elba.** It is not measurable with the available
control, and it is not needed for the science.

### What IS settled: the epoch difference, by two independent instruments

| instrument | value |
|---|---|
| DoD on stable AND open ground (116,507 cells) | **−2.12 mm** |
| leveled benchmark DG8385, gen2 − gen1 raw | **+73.00 ± 11.00 mm** |
| ⤷ minus the pipeline's +67.38 mm geoid shift, predicts | **+5.62 mm** |
| **the two agree to** | **7.74 mm** |

Unrelated failure modes, same answer. **gen1 and gen2 are already level to ~8 mm after the
geoid**, so *the DoD is correct for change detection and always was.*

### What is NOT settled: either epoch's absolute level

Three defensible control-based estimators disagree by more than any of their error bars,
and all three contradict the two instruments above:

    all covers pooled                +31.01 ± 18.17    averages the canopy mix
    open only                        +62.74 ± 23.38    2 marks inside 10 km; extrapolates
    line + cover model, read at open  +72.08 ± 36.99    right structure, unidentifiable here

The measurement is not at fault — our ties match the vendor's own published residuals at
the same marks to **−6.71 mm**. The 8 open marks on Elba's lines simply sit **+54.46 mm**
above the 230-mark open population, and the near ones disagree with the far ones (+17.20
inside 10 km against +76.33 beyond).

The benchmark cannot supply the absolute either: **±180 mm**, dominated by the mark's ±3 m
horizontal position on a ~6% slope.

`products/ANSWER_gen1_elba.json` is marked **RETRACTED** and carries this reasoning.

### What follows

* The **absolute datum is a supporting measurement, not a blocker.** It matters only for
  tying our DEMs to external absolute data — not for erosion, not for the spatial pattern,
  and **not for the forest−open contrast, where a constant cancels exactly**.
* Ground control's payoff is **statewide per-line structure**, not a single site's level.
  `run_weighted_datum_test.py` quantifies it: the weighted sd goes from 22.75 mm toward
  11 mm as per-line constants improve, which only many marks per line can deliver.
* gen2's bridge remains unmeasurable (engineered checkpoint siting), which no longer
  blocks anything now that the DoD is known to be correct as it stands.

## Guardrails

* **Every parameter is swept or declared.** Cover treatment moves gen1 by 68.03 mm; gen2
  surface × cover by 15.16 mm; the collinearity σ changes the line grouping. None is chosen
  in code.
* **A failed computation is REPORTED, never dropped.** A silently skipped radius printed
  `spread = 0.0` on 29 marks and read as certainty.
* **Undo the geoid explicitly** when comparing our GEOID18-framed surface to GEOID03 control.
* Local commits only; nothing pushed. Write only inside `ground_control/`.
* Shared laptop: crop before CSF (~4 s vs minutes); peak so far 1932040 kB, swap untouched.

## NEXT ACTIONS

1. **Decide the cover treatment** (Andy leans open-only, which also resolves the catchment's
   urban-mark confound) and the **gen2 surface** (ql1 vs ql0, 12.55 mm).
2. **`REPORT.md` and `INTEGRATION.md`** — `HANDOFF.md` §8 deliverables, still unwritten.
3. Ask whether the absolute datum is needed at all: the DoD's spatial pattern does not need
   it, and the **forest−open contrast is a difference of differences in which a constant
   cancels exactly**. It matters only for net volume.

## Reproducibility

Every number above comes from a `trust/provenance.py` run with a ledger in `.trust/runs/`.
Tracks are committed (`data/gen1_line_tracks.json`, 67 passes) and so is the swath-constant
cache, so nothing depends on a scratchpad — which is exactly how `+53.6 ± 13.0` became
unreproducible.
