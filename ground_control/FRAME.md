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

### The relation that governs it, as a closed level circuit

**Definitions.** Every constant below is a *tie*, `c = surveyed − z_lidar`, so **positive
means the surface reads LOW and the constant is what you ADD**.

| symbol | is | from | at Elba |
|---|---|---|---|
| `c1` | surveyed NAVD88 minus gen1's **delivered** surface | the 2008 MnGeo/MnDNR validation control, open ground, per flight line | **+62.74 ± 23.38 mm** |
| `c2` | surveyed NAVD88 minus gen2's **delivered** surface | the 2021 USGS **held-out** NVA checkpoints (the LCPs calibrated gen2 and are excluded) | **−6.56 mm** |
| `bridge` | delivered surface minus **our** reconstruction | per-mark comparison; carries a constant from the vendor's product onto ours | gen1 **−4.04 ± 11.12**, gen2 **unmeasured** |
| `g` | the geoid carry ADDED to gen1, GEOID03 → GEOID18 | `references.geoid_difference`, from the PROJ grids | **+67.38 mm** |

`c1` and `c2` are *not* interchangeable: each is measured against **its own epoch's**
control, and each describes that epoch's **delivered** product, not ours.

**The relation is a level circuit.** Walk from surveyed NAVD88 onto gen1, across to gen2,
and back to surveyed NAVD88. If every leg is right, the walk closes on zero.

```
            surveyed NAVD88  ......... START and END, the same datum both epochs
                   |                                              ^
       +c1 = +62.74|  2008 control                                | -c2 = +6.56
                   v                                              |  2021 held-out control
          gen1 DELIVERED                                   gen2 DELIVERED
                   |                                              ^
   +bridge1 = -4.04|  our reconstruction                          | +bridge2 = 0 +/- 26
                   v                                              |    ** UNMEASURED LEG **
            our gen1  (gen1's own geoid frame)                our gen2
                   |                                              ^
        -g = -67.38|  undo the geoid carry                        |
                   v                                              |
            our gen1  (gen2's frame) --- -DoD = +2.12 ------------+
                                          measured on 116,507 stable open cells

    leg                                                   mm
    +c1        2008 control -> gen1 delivered          +62.74
    +bridge1   gen1 delivered -> our gen1               -4.04
    -g         undo the geoid carry                    -67.38
    -DoD       our gen1 -> our gen2 (measured)          +2.12
    +bridge2   our gen2 -> gen2 delivered  UNMEASURED   +0.00
    -c2        gen2 delivered -> 2021 control           +6.56
    -------------------------------------------------------
    MISCLOSURE                                          0.0050
    expected from the legs' own uncertainties            26.06
```

Three things follow, and they are ordinary surveying practice:

1. **A misclosure far larger than the legs' combined uncertainty means a blunder or a
   mis-modelled circuit, not bad measurements.** An earlier version of this analysis missed
   by 71.42 mm because `g` had been left out of the relation. The measurements were fine;
   the *model of the traverse* was wrong.
2. **A misclosure far smaller than expected is luck, not precision.** 0.0050 mm against an
   expected 26.06 mm is a coincidence. The circuit confirms there is no gross error; it does
   not establish agreement to microns.
3. **Closure does not validate an individual leg.** With a 26 mm tolerance, a 26 mm error in
   `bridge2` is invisible — which is exactly why gen2's bridge remains open even though the
   circuit closes.

Because the constants come from three unrelated sources — two survey networks, the PROJ
geoid grids, and the point clouds — one relation among them closing to 1.92 mm on the
`c1 − c2` versus `g` reading is a real cross-validation of the *relation*. It is one
closure, not several.

## Open items — status 2026-08-31

- **Near versus far marks — RESOLVED 2026-08-31, and it was not a distance effect.** The
  split is not significant (Welch t = −1.435, p = 0.234) and is **confounded with line**:
  both near marks sit on the low-tie lines 137/138 while every far mark sits on 133/134/135.
  Within a line, where distance is unconfounded, the comparisons are small and opposite in
  sign (+153.7 over 44.7 km, −36.3 over 10.1 km, +19.0 over 25.6 km). It is therefore the
  per-line structure the estimator already averages over, not an uncorrected bias. Producer:
  `ground_control/run_nearfar_and_holdout.py`.
- **The unpublished vendor bias adjustments are ABSORBED, not a limitation on the number.**
  `c1` is measured against the *delivered* surface, which already carries the vendor's
  adjustment, so our constant is what remains after it and its value is never needed. The
  condition is that our marks were held out from that calibration, and they were: the 963
  published residuals have mean −43.41 mm, t = −10.17 against zero, p = 3.860e-23. Had they
  been the calibration set they would sit on zero by construction. This limits *explaining*
  `c1` — it cannot be decomposed into geoid error, lidar error and residual bias — but not
  *using* it. The one caveat that survives: the documentation does not say whether the
  adjustment was global, per lift or per line, so its spatial uniformity is unverified.
- **gen2's bridge** remains the one genuinely open term: bounded to 0 ± 26 mm by the
  closure, never measured directly (engineered checkpoint siting, radius spreads
  131–715 mm). gen2 is barely reprocessed, so there is little for it to be.

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
