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

## ★★ CONCLUSION 2026-08-28 — everything closes; the 2026-08-27 retraction is WITHDRAWN

**gen1 at Elba: +58.70 ± 25.89 mm on our surface** (delivered +62.74 ± 23.38, bridge
−4.04 ± 11.12). Verified by an EXACT reconciliation with the DoD:

    PREDICTED DoD = c1_ours - c2 - g = +58.70 - (-6.56) - 67.38 = -2.12 mm
    MEASURED DoD  (stable AND open, 116,507 cells)              = -2.12 mm
    MISS                                                         =  0.0000 mm

**The DoD absolute correction is +2.12 mm** — negligible. gen1 sits in the DoD *after* the
+67.38 mm geoid shift, so its constant there is `c1_ours − g = −8.68 mm`; gen2's is
−6.56 mm; the difference is +2.12, and applying it puts stable open ground at **−0.00 mm**.

### What the 2026-08-27 retraction got wrong

It concluded the control was contradicted and the absolute unmeasurable. That rested on
**my own algebra error**: I asserted the geoid cancels between `c1` and `c2`. It does not.

**NAVD88 is the DATUM; GEOID03 and GEOID18 are MODELS** for converting GPS ellipsoidal
heights to orthometric. Both control sets publish NAVD88 and are directly comparable — but
each epoch's **lidar z** was converted with a different geoid model, so `c1` and `c2`
reference surfaces in different frames. The geoid does not cancel; it is exactly the term
the pipeline adds. The correct relation is `DoD = c1 − c2 − g`.

So nothing was contradicted. The DoD, the leveled benchmark DG8385 and the 2008/2021
control **all agree**, and gen2's unmeasured bridge is ~0, which is why it never blocked
anything.

### The signature I misread three times

The discrepancy equalled **67.38 mm — the geoid — to the digit**, and it surfaced three
separate times (the raw bridge at L1O101; the H1/H2 hypothesis test; the implied gen2
bridge). Each time I called it suspicious and then looked OUTWARD for the fault. **A
discrepancy that matches a known constant to the digit is a bookkeeping signature, not a
measurement failure** — check your own arithmetic first.

## ★ WHY THE CORRECTION IS REQUIRED, even at +2.12 mm

`align_swaths` is gauged on the lowest-numbered flight line. That choice touches no
swath-to-swath difference, but it sets the absolute level the mosaic inherits — which
becomes **the reference line's own vertical error**. Measured on elbaext:

    133 +0.00   134 +22.00   135 +6.20   136 -9.80   137 -18.40   138 -22.60
    => re-gauging on a different line moves EVERY elevation by up to 44.60 mm

**The gauge choice is worth 21x the correction.** So an uncorrected elevation is an
arbitrary implementation detail (`ref=int(ps.min())`), not a measurement, and the
correction's smallness at Elba is a property of line 133 having been a lucky pin.

Applying a control datum removes the dependence *exactly*: `corrected = z + c` with `c`
measured against the same gauged product, so re-gauging by `d` shifts `z` by `+d` and `c`
by `-d` and they cancel. `tests/test_apply_datum.py` demonstrates this — uncorrected
spread 44.60 mm across the six gauges, corrected spread < 1e-9.

**This is now a required pipeline step.** `pipeline.difference_dem` records
`swath_gauge_ref` and leaves `absolute_datum_mm` None until a constant is supplied, and
its docstring states that the absolute level is gauge-dependent. Applier:
`ground_control/apply_datum.py`; a constant is tied to its gauge and must be re-expressed
via `regauged_to()` if `ref` changes.

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
