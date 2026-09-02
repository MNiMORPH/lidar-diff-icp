# A swath tie that no estimator produced — four rules tested, none adopted

## The observation

Battle Creek's `beam_offset_table` refused to build: swath 1102 has no along-track drift
curve. Chasing that turned up a second thing, which is the subject here.

`corrections.json` ships `1102 dz = +0.0137 m` as a solved constant. The pair that supplies
it reports:

    1014-1102   dz -3.4640 m   n 0   n_iter 1   dx_sigma NaN   converged False

Swath 1102 crosses the tile as a strip **2.8 m wide by 67.3 m long** — 89 of 4,166,880
returns. Nuth & Kääb regresses `dh/tan(slope)` on ASPECT, so a strip that narrow samples
essentially one aspect, the design matrix is rank-deficient, and there is no horizontal
solution to find. Under `tie="intercept"` the vertical is then extrapolated to
`tan(scan) = 0`; the same cells give **+0.0600 m** under `tie="overlap_median"`.

The edge enters the network with weight `sqrt(n) = 0`, so it constrains nothing — and the
swath, now unconstrained, takes a **minimum-norm value from `lstsq`**, which is shipped as
if measured. The tile's whole misclosure, 3.484 m against millimetres elsewhere, is this.

## Four rules tested against the real sites. Every one misclassifies a real case

    n == 0                  WRONG. The rigid vertical fallback legitimately reports n=0
                            with a valid dz -- the code says so in a comment, correctly.

    converged is False      TOO BLUNT. Carlton 89-90 is converged=False and sound:
                            22,293 cells, finite dx_sigma 0.0070, dz -0.0143 m, and it
                            IMPROVED the scatter (NMAD 0.1705 -> 0.1679). It hit max_iter,
                            which is slowness, not failure. Dropping it leaves 1.99% of
                            Carlton -- 161,842 returns -- unaligned.

    n_dz <= 0               TOO PERMISSIVE. The intercept fit reports a population of 36
                            for 1014-1102, so its -3.4640 m survives. (It also broke a
                            passing pipeline test.)

    dtan must bracket 0     WRONG. 1012-1014 spans dtan -0.9069..-0.7810, never reaching
                            zero, on 37,436 cells -- and gives a perfectly sound -0.0072 m.

## What actually separates them, and why it is not adopted

The lever arm of the extrapolation:

    pair          cells   dtan min   dtan max     span        dz
    1012-1013    53,599    -0.7537     0.2215   0.9752    +0.0287
    1013-1014    91,536    -0.7899     0.7717   1.5616    -0.0363
    1012-1014    37,436    -0.9069    -0.7810   0.1259    -0.0072
    1014-1102        36    -0.8083    -0.7899   0.0184    -3.4640
    1101-1102        36    -0.7357    -0.7357   0.0000    +0.0000

`1101-1102` has a dtan span of exactly zero: every cell at one across-track position, so the
slope is unidentifiable and the intercept is undefined. `1014-1102` extrapolates ~0.8 in
dtan from a span of 0.018 — a lever ratio near 43:1.

Turning that into a cut means inventing a ratio threshold, which is the trap this project
exists to avoid. **The threshold-free answer is to weight the network by each observation's
own variance** rather than by `sqrt(n)`: the extrapolation variance of an intercept grows
with `(distance / span)^2`, so a tie like these self-downweights to nothing with no constant
chosen by anyone. That is a change to the weighting scheme — which `coregister_swaths`
currently holds fixed ON PURPOSE, so the two tie modes differ in the vertical estimator
alone — and so it is a design decision, not a guard.

## What was changed

`Coreg` gains `n_dz`, REQUIRED with no default: the number of cells that actually determined
the reported `dz`. It is not `n` (the cosine fit's count) — the rigid fallback has n=0 while
its dz rests on a real overlap, and the intercept tie's dz comes from a different population
again, whose count `across_track_tie` computed and threw away.

Nothing filters on it yet. The diagnostic exists so the question can be asked; the answer is
open.

## Also found, not fixed

* `data/before_battlecreek/gen1_4tile.laz` has `point_source_id = 0` for all 24,642,050
  returns — the merge flattened every flight-line id. Per-swath work against that file would
  see one giant swath. `scripts/merge_gen1_tiles.py` exists precisely because an earlier
  merge dropped fields silently; this file predates it.
* A rigid-fallback edge with n=0 is kept but enters at weight `sqrt(n) = 0`, so the case the
  comment protects has never actually constrained anything.
