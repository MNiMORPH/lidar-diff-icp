# Session summary — ground control for the MN lidar vertical datum

**2026-08-27/28.** 35 commits, all in `ground_control/` plus one hook fix. 138 provenance
ledgers in `.trust/runs/`. `FRAME.md` is the state anchor; this is the narrative.

## The brief, and what actually happened

The brief asked for a reusable method giving a datum constant and uncertainty at any site.
What the work produced is that method **and** a correction to the project's understanding
of what it needs one for.

**Coregistration was never broken.** The absolute datum never worked: the shipped
`z_before_absolute` carried **±39.7 mm on a value of +22.7 mm**, built by chaining two
*2021* marks onto a *2008* surface. git shows that was chronology, not choice — the product
was written 83 minutes before gen1's own 2008 control entered the repo.

## Results

| quantity | value |
|---|---|
| gen1 at Elba, delivered surface | **+62.74 ± 23.38 mm** (open ground, per-line, 8 marks / 5 lines) |
| bridge, delivered → ours | **−4.04 ± 11.12 mm** (29 open marks) |
| gen1 at Elba, our surface | **+58.70 ± 25.89 mm** |
| gen2, its own control, delivered | **−2.37 ± 2.37 mm** project-wide; QL1 block **−6.83 ± 2.96** |
| geoid term the pipeline adds to gen1 | **+67.38 mm** |
| **DoD absolute correction** | **+2.12 mm**, uncertainty **±26.06** |

**The closure.** `DoD = c1_ours − c2 − g` predicts −2.12 mm against a measured −2.12 mm on
116,507 stable open cells, **miss 0.0050 mm**. Independently, the control's own epoch
separation `c1 − c2 = +69.30` recovers the PROJ geoid difference **+67.38** to **1.92 mm**
— two survey networks reproducing a geoid model they know nothing about.

**But the correction is not robust.** It is a small difference of large numbers, so it
inherits `c1`'s full ±25.89, and it swings **54.88 mm** (−7.22 to +47.66) across choices
still open. It is consistent with zero only because its uncertainty is 12× its value.

## What was built

`control.py` (epoch-agnostic adapter), `lines.py` + committed tracks (67 passes, 41 psids),
`same_line.py` (catchment-free estimator), `our_surface.py` (local reconstruction, gated at
+1.7 mm against the shipped grid), `datum.py`, `gen2_swath_deviation.py`, 14 drivers, and
12 tests including a regression **proven to bite** (`assert 209 == 230` → the L1O/L10 trap,
restored byte-identically).

Method decisions that changed answers: open ground only; the flight **line** as the unit of
replication (design effect 1.40×); candidates bounded by data not by a radius; passes merged
into physical lines by collinearity scaled by its own prediction sd.

## Errors made and corrected, in order

1. **Heading drift read as a gap** — an 800 m across-track separation is what one line at
   179.3° looks like over 94 km, not two lines.
2. **A silent estimator failure** — `radius_spread = 0.0` on all 29 marks read as certainty;
   in fact one radius enclosed 4 cells against the 6 an order-2 fit needs and was skipped.
3. **The geoid trap in the bridge** — raw gap −59.7 mm of which +69.3 was the conversion.
4. **Conflated "found near a track" with "belongs to that line"** — cost ~9 mm.
5. **Over-trusted my own along-track null**, which licensed distant marks.
6. **Retracted a correct result** on my own algebra error: I asserted the geoid cancels
   between `c1` and `c2`. It does not — `DoD = c1 − c2 − g`. The signature (67.38 mm, to
   the digit) appeared three times and I looked outward each time.
7. **Mis-arranged signs** in a decomposition (`−5.96` where `+2.12` was right).
8. **Promoted a rounding to exactness** (0.0000 for 0.0050).

The trust gate caught four of these classes; it fired on typed-not-pasted numbers three
times and on a missing provenance banner once, and was right every time.

## Open

* **Near vs far marks.** Six of eight open marks are 14–63 km away and disagree with the two
  near ones by **59.13 mm**. This is the live question; it changes the answer's sign.
* **Mechanism, not relation.** The relation is verified; "most of the difference was the
  geoid" is consistent but unproven, and competes with the **unpublished vendor bias
  adjustments** both epochs carry — neither value published, neither recoverable from data.
* **gen2's bridge** — bounded to 0 ± ~26 mm by the closure, never measured directly
  (engineered checkpoint siting, radius spreads 131–715 mm).
* **The statewide per-line correction** is where control actually pays off: the weighted sd
  falls from 22.75 mm toward 11 mm only as per-line constants improve, which needs many
  marks per line, not more marks at one site.

## What a reader should take

The DoD is correct for change detection and always was. The absolute constant matters only
for external tie-in — **not for erosion, not for the spatial pattern, and not for the
forest−open contrast, where a constant cancels exactly.** The honest state of the absolute
is +58.70 ± 25.89 mm at Elba, with a mechanism that is not established.
