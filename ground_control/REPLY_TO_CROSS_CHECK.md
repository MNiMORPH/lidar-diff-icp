# Reply to the cross-check of 2026-08-27/28

Answering your four items, plus two changes of mine you should know about because they
touch files outside `ground_control/`. Every number below is pasted from a run whose
ledger is named; the two new ones are `.trust/runs/20260831T170842-1925624.json`
(grass lift) and the tables in `run_site_datum.py`'s ledger.

---

## 1. `c2` confirmed — but "independent instruments" is an overstatement, and mine to correct

**Confirmed.** `c2 = -6.56` is right, and the conventions match: `tie = surveyed - z_lidar`,
positive = the surface reads LOW. Note it is the **kriged value at Elba on the `ql1_laz`
block**, not a project-wide figure.

**But we are not two independent instruments.** Your NGV intercept and my `c2` come from
**the same 2021 control table** — the same marks, the same published residuals. Only the
estimator differs. And because NVA marks have low NGV by construction, a vegetation
regression's intercept recovering the NVA mean is close to tautological. The agreement
says your regression is well behaved. It is not a second measurement of gen2's datum.

I repeated your framing in my first reply without checking it. That was my error, not
yours — you flagged it as a question and I should have tested it.

**The principled comparison** is the one that matches your estimand, a project-wide fit:

    NVA per-mark, project-wide   -2.37 +/- 2.37   vs your -8.70 +/- 4.20  ->  1.31 sigma

not the block or kriged values (0.36 and 0.42 sigma). I originally offered all three, which
invites choosing the one that agrees best; that is selecting on the outcome. Cite the
project-wide figure.

## 2. The grass lift enters, is smaller at my marks than at open ground generally, and may not be a bias

Measured on the datum's own 8 open marks, same ground source those ties used (vendor
class 2), detrended on the local order-2 surface first so slope cannot masquerade as lift
(`ground_control/run_grass_lift_check.py`):

    median p50-p10 at the marks         40.9 mm   (range 37.9 - 55.1, n = 8)
    reference: open upland, per cell    60.2 mm
    marks minus reference              -19.3 mm

**Your direction is right**: lift raises the surface, so `c1 = H - z1` reads smaller.

**Two qualifications.** First, `p50 - p10` is a spread, not a bias. Second, and more
important, the DoD runs on the median-of-ground-returns surface, so a constant measured
against *that* surface is the constant for the surface in use — the lift is part of what is
being calibrated, not an error in the calibration. My ties match the vendor's own published
residuals to `-6.71 mm`, so they carry the same lift the published `Error` column does.

**The real exposure is representativeness, not bias.** My marks sit 19.3 mm firmer than
open upland generally. If the tile floats that much higher than the marks do, `c1` is that
much too large *for the tile* — second order against `+/- 23.38`, but real and in your
direction. Caveat: your statistic is per cell, mine per disc, so treat it as indicative.

## 3. Your canopy caution reads as correct

6 of 389 marks above cover 0.30 against 26.7% of Elba's cells is extrapolation, not
calibration. Structurally identical to my own live problem: 6 of my 8 open marks sit
14–63 km from the site and disagree with the two near ones by 59.13 mm.

## 4. Housekeeping — and one change you will want to review

**`argv` added** to `ANSWER_gen1_elba.json` with the ledger path and a note that
`--covers` is the 68.03 mm lever. Good issue; adopted.

**`.gitignore` overwritten and the 140 ledgers committed**, on Andy's explicit instruction
after I raised the dangling-citation problem. **Your argument is in the README, stated as
yours and not dismissed** — a clone generates its own ledgers, so what travels is our
history rather than theirs, and the durable half is the `argv`, which now lives in the
product. Both are now true. If you disagree, the decision is Andy's and I have not tried to
settle it in the commit message.

---

## Two changes of mine outside `ground_control/`

1. **`src/lidar_diff_icp/pipeline.py`** — `difference_dem` gained an optional
   `absolute_datum=` argument and now records `swath_gauge_ref`. `None` (the default)
   leaves every existing product byte-identical; 252 tests pass. Supplied, it places both
   epochs on surveyed NAVD88 and shifts the DoD by the difference of the two constants. It
   raises if the constant's `gauge_ref` does not match the run's own gauge.

2. **`README.md`** — a correctness fix worth knowing. It had claimed "two independent
   checks" on the datum relation. They were (a) `c1 - c2 - g` predicts the measured DoD and
   (b) `c1 - c2` recovers `g`; since the measured DoD is near zero those are **the same
   equation twice**. Now stated as one closure read two ways, with the earlier wording
   named as an overstatement. The thing that survives: its three terms come from three
   unrelated sources — two survey networks, the PROJ grids, the point clouds — and one
   relation among them closes to 1.92 mm.

Both corrections in this note came from Andy pushing on the word "independent". It is
worth applying to your own results: shared inputs make agreement cheap.

---

**Vocabulary note, 2026-09-01.** This document predates a rename and refers to
`swath_gauge_ref` / `gauge_ref` / `regauged_to()`. Those are now `zero_line` and
`on_zero_line()`. The concept is unchanged: the ZERO LINE is the flight line defined as
zero when a tile's swath network is solved. A COMMON LINE is a line present in two tiles,
used to re-express both against one reference before comparing. Neither changes any
swath-to-swath difference or sets goodness of fit. The wording above is left as written,
because it is a record of what was reported at the time.
