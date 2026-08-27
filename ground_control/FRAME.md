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

## Where the number stands — ADOPTED 2026-08-27

**`ground_control/products/ANSWER_gen1_elba.json` is the canonical gen1 surface offset.**

| surface | offset (mm, ADD) | uncertainty | of what |
|---|---|---|---|
| **DELIVERED** (what the vendor shipped; what the 2008 control measured) | **+62.74** | ± 23.38 | SE over the 5 flight lines, 8 marks |
| bridge, delivered → ours | −4.04 | ± 11.12 | SE of the mean over 29 open marks |
| **OURS** (CSF, swath-aligned, geoid-shifted 5 m grid the DoD runs on) | **+58.70** | ± 25.89 | the two in quadrature |

Producer: `run_datum_by_returns.py` (catchment-free) + `run_bridge_wide.py`, assembled by
`run_gen1_elba_answer.py`. Ledger `.trust/runs/20260827T220404-3482736.json`.

**Superseded, and why.** The `+55` lineage (2026-08-26 ledger +55.0; my catchment rebuild
+54.77) was inflated ~9 mm by marks on line segments that are NOT collinear with Elba's own
passes. A catchment search conflates "found near a track" with "belongs to that line"; the
catchment-free estimator tests each mark's own pass and rejects three (138.0 at 19.3σ,
138.3 at 4.1σ, 136.2 at 3.5σ), which drops line 136 entirely.

**Cover is now OPEN-ONLY (`L1O`), 2026-08-27**, and it bought more than it cost. It moves
the datum **+17.17 mm** (from +45.57 on `L1O+L5U`) and shrinks the `collinear_sigma`
sensitivity from **8.69 mm to 1.90 mm** — σ = 3 and σ = 5 now give the identical answer, and
no mark is rejected by collinearity at all. The σ ambiguity was entirely an urban-mark
phenomenon: every mark whose pass assignment was in doubt (136.2 at 3.5σ, 138.3 at 4.1σ,
138.0 at 19.3σ) was `L5U`. Cost: 14 marks → 8, SE 18.81 → 23.38.

The DoD absolute correction still cannot be assembled cleanly: it differences two constants
and only the gen1 side has been carried onto our surface.

## What was corrected, and what it cost

* **`+53.6 ± 13.0` is reconstructed.** The value reproduces (**+54.77** vs the 2026-08-26
  ledger's **+55.0**, −0.23 mm). The **±13.0 was the SE over MARKS**; over LINES it is
  ±18.21, a ratio of **1.40×** — exactly the documented design effect. Marks under one line
  share that line's unknown constant. The trust ledger reconstructed a number whose script no
  longer exists, which is what it was built for.
* **`point_source_id` is not a flight line — but my reason was wrong.** Across-track
  separation is NOT evidence of two lines: a near-N–S line drifts ~1.1 km in easting over
  94 km. Valid test = extrapolate and scale by the extrapolation's own prediction sd.
* **The catchment radius is a SEARCH bound, not a criterion.** Set it wide; the returns
  discriminate. A narrow one drops marks whose returns are unanimous.

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
