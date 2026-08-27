# The gen1 headline +53.6 ± 13.0 mm does not reproduce from committed code

**2026-08-27.** Recorded because a load-bearing number whose producer is gone is the
continuity failure this project is built to avoid, and because the gap is in the
*uncertainty*, which is the half that matters for the DoD.

## What the frame claims

`analysis/FRAME_2026-08-26-PM.md`: *"Elba's offset: +53.6 ± 13.0 mm — the constant to ADD
to gen1 after swath alignment. 18 open/urban marks on Elba's own flight lines (133-138),
Elba's swath constants applied to put them in one frame, SE over marks."*
It names its producer as scratchpad `screen_marks.py` with results `screen_results.csv`.

**Neither exists.** `find . -name 'screen_marks*'` returns nothing. The figure appears only
in prose, in four `.md` files. No committed script computes it.

## What committed code gives at the same site

`analysis/groundtruth/gen1_datum_at_site.py --easting 578762.8 --northing 4884487.6
--radius-km 20 --tiles data/before --mode per_line --cover ...`, ledgers
`.trust/runs/20260827T173649-3294897.json` and the `--cover L1O,L5U` run beside it.
Restricted to marks whose `assign_line_from_returns` line is one of 133-138:

| selection | n marks | n lines | estimate |
|---|---|---|---|
| open (`L1O`) | 4 | 4 | mean **+39.03**, median +17.20 |
| open+urban (`L1O,L5U`), marks pooled | 7 | 4 | **+40.56 ± 23.21** |
| open+urban, per line (lines as unit of replication) | 7 | 4 | **+41.20 ± 25.94** |

Line means: 133 −9.9, 134 +111.2, 137 +17.9, 138 +45.7 mm.

**Sign agrees** — gen1 reads LOW at Elba on its own lines, as the frame says. **Magnitude
is lower** (~+40 against +53.6). **The uncertainty is about twice the frame's**: 23-26 mm
against ±13.0.

## Why the two differ, identified

Not a sign or estimator problem: a **different discovery mode**.

* The committed driver discovers marks near the **site** (`discover_near_point`,
  `--radius-km`). At 20 km only 7 open+urban marks sit on lines 133-138.
* The frame's 18 marks must have come from discovery near the **lines**
  (`gen1_datum.discover_near_lines`), which finds marks anywhere along a track, tens of
  kilometres from the site, and is the correct search for a per-line quantity.

`discover_near_lines` needs a track geometry per line. Its producer,
`analysis/groundtruth/gen1_line_tracks.py`, writes `line_tracks.json` **to the
scratchpad**, which is gone. Regenerating it means re-reading the 47 gen1 tiles in
`data/before/`.

Note `discover_near_lines`'s own docstring: *"This is a search, not an assignment."* The
line that actually hit a mark is settled by `assign_line_from_returns` from the returns;
a fitted track walks off by hundreds of metres over tens of kilometres.

## What this does and does not establish

* It does **not** show +53.6 ± 13.0 is wrong. It shows it is **unreproduced**.
* It does show that the only estimate reachable from committed code today is
  **+40.56 ± 23.21 mm** (open+urban, 7 marks, 4 lines), and that its uncertainty is
  roughly double the frame's.
* Because the DoD absolute correction is a difference of two ~10 mm-scale constants, a
  ±13 vs ±25 difference in gen1's term is not cosmetic.

## What it takes to settle it

1. Re-derive the flight-line tracks and **commit them** rather than write them to a
   scratchpad (~47 tiles read; one heavy job, shared laptop).
2. `discover_near_lines(control, {133..138}, half_width_m=481)` -- 481 m is the measured
   half of the line spacing, the vendor's own class-12 seam, from
   `analysis/GEN1_DATUM_MORE_MARKS.md` §1; it is a measurement, not a chosen radius.
3. Confirm each candidate with `assign_line_from_returns`, keep those truly on 133-138.
4. Combine with lines as the unit of replication.

Until then the frame's ±13.0 should not be carried into a product.

## The lesson, in the project's own terms

A number that drives a decision must be produced by code that is committed. The producer
here was a scratchpad script and its inputs were scratchpad JSON; both evaporated, and
with them the ability to check the number or extend it to another site. The statewide
goal needs exactly this estimator at every site.
