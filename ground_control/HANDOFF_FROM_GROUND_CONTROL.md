# Handoff — the gen1/gen2 absolute intercomparison is closed

**2026-08-31, from the `ground_control/` session to whoever picks the repo up next.**
Read `FRAME.md` first for state; this is the narrative and the next actions.
53 commits, all on `main`. 268 tests pass. 142 provenance ledgers tracked.

---

## 1. What was wrong, and what is now true

**Coregistration was never broken.** The absolute datum never worked: the shipped
`z_before_absolute` carried **±39.7 mm on a value of +22.7 mm**, tied gen1 at two *2021*
marks chained 5–6 flight lines, and git shows that was chronology rather than choice — it
was written 83 minutes before gen1's own 2008 control entered the repo.

**Adopted, `ground_control/products/ANSWER_gen1_elba.json`:**

| quantity | value |
|---|---|
| gen1 at Elba, delivered surface | **+62.74 ± 23.38 mm** (open ground, 8 marks on 5 lines) |
| bridge, delivered → our reconstruction | **−4.04 ± 11.12 mm** (29 open marks) |
| gen1 at Elba, our surface | **+58.70 ± 25.89 mm** |
| gen2, its own held-out control | **−2.37 ± 2.37 mm** project-wide; **−6.83 ± 2.96** in QL1 |
| geoid carry added to gen1 | **+67.38 mm** |
| **DoD absolute correction** | **+2.18 mm**, putting stable open ground at −0.003 mm |

The level circuit closes at **0.0050 mm** against an expected **26.06** — which confirms
there is no gross error and does *not* establish agreement to microns. See the ASCII
diagram in `README.md` and `FRAME.md`.

## 2. Why the correction is required even though it is 2 mm

`align_swaths` is gauged on the lowest-numbered line, so the mosaic inherits **that line's
own vertical error** as its absolute level. Re-gauging moves every elevation by up to
**44.60 mm** at elbaext — **21× the correction**. Applying a control datum cancels that
exactly (`tests/test_apply_datum.py`: uncorrected spread 44.60 mm across six gauges,
corrected below 1e-9). **Judge a correction by what it removes, not by its magnitude.**

## 3. Changes outside `ground_control/` — please review

1. **`src/lidar_diff_icp/pipeline.py`.** `difference_dem` gained an optional
   `absolute_datum=` dict and now records `swath_gauge_ref`. `None` (the default) leaves
   every existing product byte-identical; 268 tests pass. Supplied, it corrects **both**
   epochs and shifts the DoD by the difference. It **raises** if the constant's `gauge_ref`
   does not match the run's gauge, because a constant belongs to the product it was
   measured against.
2. **`.gitignore`.** `.trust/` is now tracked, on Andy's instruction, because `README.md`,
   `FRAME.md`, `REPORT.md` and the adopted products cite ledgers by filename. **The counter
   argument is stated in the README as yours and not dismissed** — a clone generates its own
   ledgers, and the durable half is the `argv`, which now lives in the product.
3. **`README.md`.** Substantially extended: the ordered recipe, the three tiers
   (within-epoch / between-epoch registration / external tie), reductions versus
   corrections, the datum section, and the level-circuit diagram with `c1`/`c2` defined.

## 4. What is settled, and what is not

**Settled.** The relation `DoD = c1 − c2 − g` (the geoid does *not* cancel: NAVD88 is the
datum, GEOID03/GEOID18 are conversion models and each epoch's lidar used a different one).
Epoch-matched control. Open ground only. The flight line as the unit of replication
(design effect 1.40×). The returns assign the line, never the geometry.

**Resolved late, after having been carried as open:**
- *Near/far marks* — not a distance effect. Not significant (Welch t = −1.435, p = 0.234)
  and confounded with line; within a line the comparisons are small and opposite in sign.
  The estimator already handles it by aggregating per line.
- *The unpublished vendor bias adjustments* — **absorbed**. `c1` measures the delivered
  surface, which already carries them, so their values are never needed. Hold-out verified:
  963 published residuals, mean −43.41 mm, t = −10.17, **p = 3.860e-23**.

**Still open.**
- **gen2's bridge** is the one genuinely unmeasured leg: bounded to 0 ± 26 mm by the
  circuit, never measured directly because its four usable checkpoints sit on engineered
  ground (order-2 fit RMS 182–479 mm, radius spreads 131–715 mm). A *siting* problem; more
  checkpoints of the same kind will not fix it.
- **Spatial uniformity of the vendor adjustments** is unstated in the documentation — global,
  per lift or per line is not recorded.
- **The statewide per-line correction** is where control actually pays off: the weighted
  uncertainty falls from 22.75 mm toward 11 mm only as per-line constants improve, which
  needs many marks per line rather than more marks at one site.

## 5. Next actions, in the order I would take them

1. **Nothing is required to use the DoD.** It was correct for change detection all along;
   the absolute constant matters only for external tie-in, and it cancels exactly in the
   forest−open contrast.
2. **Apply the datum wherever an absolute elevation is quoted**, via
   `run_site_datum.py` → `difference_dem(absolute_datum=...)`. Not for accuracy, for
   gauge-invariance.
3. **Check gen2's horizontal accuracy.** Step 4 registers gen1 *to* gen2, so a lateral
   displacement in gen2 propagates into gen1 invisibly. The 2021 checkpoints carry surveyed
   *horizontal* coordinates, so this is checkable and has not been done. A lateral error on
   a slope becomes a vertical one.
4. **Take the datum to a second site** before trusting the method's generality. Everything
   here rests on one tile and 8 open marks.

## 6. Rules that cost the most to learn

- **A discrepancy equal to a known constant to the digit is a bookkeeping error in your own
  arithmetic.** The geoid's 67.38 mm appeared three times before I stopped looking outward;
  I retracted a *correct* result on the strength of it.
- **Never report a computation that silently dropped an input.** A skipped radius printed
  `spread = 0.0` on 29 marks and read as certainty.
- **Shared inputs make agreement cheap.** Two "independent" checks turned out to be the same
  equation twice; a cross-session comparison turned out to use the same control table.
- **Commit the producer, not just the number.** The headline gen1 datum was reconstructed
  after its script had been deleted, from a ledger alone.
