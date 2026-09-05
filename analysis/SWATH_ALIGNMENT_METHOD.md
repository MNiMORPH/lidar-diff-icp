# The gen1 swath-alignment method

The chain that takes a raw 2008 flight-line mosaic to one internally consistent, absolutely
levelled surface. Code: `pipeline.register_gen1` → `coreg.align_swaths` →
`coreg.coregister_swaths` → `coreg.nuth_kaab` / `coreg.across_track_tie`, then
`ground_control/apply_datum.py`.

It is **gen1-internal through step 6**: no gen2 data is read and no absolute level is set.
That separation is what stops the swath tie absorbing a cross-epoch correction, and a test
asserts `register_gen1`'s signature carries no gen2 argument.

---

## 1. Ground classification

CSF on the before cloud (cached; deterministic and slow, ~minutes/tile), or the raw
last-return heuristic. Ground points are `rn == nr`, **singles included** — dropping them
empties flat open ground.

## 2. Boresight roll — optional, and first

`z -= b · scan_angle`, removed per point **before** the empirical alignment, so the network
cannot absorb an instrumental term. Off by default. Calibrated from gen1 self-overlap, so it
is gen2-free and needs no iteration.

## 3. Pairwise co-registration

For every overlapping pair: Nuth & Kääb on density-robust per-cell median-Z surfaces at 2 m.
The horizontal shift comes from the `dh/tan(slope)` vs aspect cosine fit.

> Nuth, C. & Kääb, A. (2011), *The Cryosphere* **5**, 271–290.

## 4. The vertical tie

The horizontal solution is untouched by this choice — the aspect fit removes the median of
`dh` before fitting, so it never sees the vertical constant.

| mode | what it is | where used |
|---|---|---|
| `overlap_median` | `median(z_ref − z_src)` over the covered part of the overlap | the **library default** |
| `intercept` | the LAD intercept at `dtan = 0`, the sidelap centre where both lines view the ground at the same incidence | what the **pipeline** passes |

`intercept` **nests** `overlap_median`: `_lad` is median regression, so an intercept-only
design converges on the median. The only difference is the slope term.

Why it matters: the between-line offset is not a constant. It sweeps 15–48 mm with
across-track position; `c = +80.0 ± 5.7` mm per unit tangent pooled, per-pair +34.5 to
+192.7, homogeneity across pairs rejected at `p = 6e-64`
(`analysis/SWATH_ACROSS_TRACK_TEST.md`). So `median(dh) = k + c·mean(dtan)` depends on which
part of the overlap the tile covers, and two tiles of different extent get different ties for
the same pair of flight lines.

**Provenance of this estimator: no source is cited for it.** `coreg.py` cites Nuth & Kääb
(2011) and DeLong et al. (2022); `across_track_tie` cites nothing. As this repository records
it, the construction is ours. Whether prior art exists is an open literature question.

## 5. The free network, weighted by information

Least squares for per-swath `(Δx, Δy, Δz)` making all overlaps mutually consistent,
**weighted by `1/dz_var`** — not by cell count (changed 2026-09-05).

Each edge reports the variance of *its own* estimator: the median tie's sampling variance, or
the intercept's standard error including the extrapolation lever arm
`scale·sqrt(1/n + mean(t)²/Σ(t−mean(t))²)`. A tie fitted far from `dtan = 0` over a short span
downweights itself, with **no threshold chosen by anyone**.

Measured on Battle Creek (5 swaths, 7 edges, so loops): the sliver pair 1014–1102 —
36 cells, `dtan` span 0.0184, extrapolating ~0.8, a 43:1 lever ratio — returns
`dz = −3.4640 m` and now carries **0.0000%** of the network weight. Swath 1102 moves
+13.6513 → +27.0895 mm, into line with its neighbours. Largest misclosure excluding that edge:
**19.314 → 6.468 mm**.

Nothing is filtered. Non-finite observations are dropped; everything else enters at its
earned weight, and a bad observation surfaces as a large residual on its own edge rather than
moving the solution.

## 6. The zero line — and it is deliberately arbitrary

The network is solved **free**; only then is the reference swath's value subtracted
(`coreg.py:577`). So the zero line touches no swath-to-swath *difference*. It does set the
absolute level the whole mosaic inherits, because that level becomes the reference line's own
error: at elbaext the six per-swath dz span **44.60 mm**.

*Naming.* The object is the **zero line** — `corrections.json` records `zero_line`,
`align_swaths` takes `ref`, and `swath_gauge_ref` / `gauge_ref` / `regauged_to()` were
renamed away on 2026-09-01. The ADJECTIVE survives, deliberately and in current code:
"gauge-invariant" is the property step 7 delivers, and `apply_datum.gauge_invariance_residual`
is the function that demonstrates it. A COMMON LINE — a line present in two tiles, used to
re-express both before comparing — is the related term.

`corrections.json` records `zero_line`, `swath_tie` and
`absolute_level_depends_on_zero_line: true`, so two products can be related.

## 7. Correct the reference line — the ground-control datum

**This is what removes the arbitrariness of step 6, and it removes it exactly.** It is what
makes a product's elevation *gauge-invariant* — the one place the adjective is the right word.

Re-expressing on a different zero line shifts every elevation by `+d` and the measured
control constant by `−d`,
so the corrected surface is unchanged (`pipeline.py:1093`). Applied to BOTH epochs, so the DoD
moves by the *difference* of the two constants and true change on stable ground goes to zero.

Demonstrated rather than asserted, on elbaext's real per-swath dz
(`ground_control/tests/test_apply_datum.py`, 4 tests passing):

    spread across zero-line choices, uncorrected   44.60 mm
    spread across zero-line choices, corrected     < 1e-9 mm
    ratio                                          > 1e9
    on_zero_line() re-expression is reversible to 1e-9

`apply_datum` REFUSES if the datum's `zero_line` does not match the product's: a constant is
tied to the zero line it was measured on, and applying it to a product solved on a different
one would silently mis-level the tile.

Elba's constants, open ground only: gen1 delivered **+62.74 ± 23.38 mm**, ours
**+58.70 ± 25.89**. The datum is measured on **open ground only** (gen1 `--covers L1O`,
gen2 `--point-types NVA`), never pooled — pooling would bake canopy response into the datum
and pre-decide the canopy-versus-erosion question.

## 8. Coverage — recorded, not acted on

`swath_coverage` reports which lines contribute no exclusive cells on the analysis grid:
their constants cannot be checked against any cell they alone determine. Recorded because
exclusivity is grid-dependent and a swath can still be a shared cell's majority.

---

## Known weakness

Both Elba tiles take their **zero line from an edge-cut swath** whose nadir track falls
outside the tile
(elba/135 at 577,143 against x₀ = 577,492.8), so it is sampled over an 8–10° one-sided scan
range and its constant and across-track slope correlate at 0.99
(`SWATH_ACROSS_TRACK_TEST.md` §7). Gauging on the line with the widest two-sided scan
coverage, or using the zero-mean origin `align_swaths` already supports (`ref=None`), costs
nothing.

It matters for an **uncorrected** product and for relating two tiles that carry no datum.
Where step 7 is applied, it cancels.
