# gen1 flight-line (swath) alignment: what the method is, and what a more principled one would be

**Date:** 2026-09-18. **Status:** STUDY ONLY — no pipeline code, product or tile was changed.
**Tile:** `data/derived/elbaext_delong` / `data/derived/elbaext_independent` (identical swath
solution; verified below). **gen1 cloud:** `data/csf_cache/elbaext.las`.
**Python:** `./lidar-icp/bin/python`.

Every number in this document is either (a) quoted from a file with its path/line, or
(b) pasted from a command run in this session. The scripts are in the session scratchpad
(`edgeset.py`, `pairdetail.py`, `tilt_vs_shift.py`, `imposed_bias.py`,
`per_line_invariance.py`, `fp_minus_upland.py`, `star_tie.py`, `dod_by_line.py`,
`dod_by_line_variants.py`) and are reproduced in §8. Claims are tagged
**[CODE]** (verified by reading the source), **[MEASURED]** (run here), **[LIT]**
(retrieved literature, URL given), **[INFERENCE]** (reasoning, not measurement).

---

## 1. VERDICT on topology: it is a TREE — in fact a bare CHAIN — and worse than the brief assumed

**[MEASURED]** At the resolution the code actually uses (`res=2.0`, `align_swaths`' own
default; `register_gen1` passes no `res`), the network on elbaext is
**6 nodes, 5 edges, zero loops**. The non-adjacent pairs are not "dropped as NaN" — their
**bounding boxes are disjoint**, so `coregister_swaths` raises `ValueError` and they are
never even candidate edges:

```
      pair           bbox  grid cells  ref finite  BOTH finite
133-134      bbox-overlap      543750      430758       415771
133-135     disjoint-bbox           0           0            0
133-136     disjoint-bbox           0           0            0
133-137     disjoint-bbox           0           0            0
133-138     disjoint-bbox           0           0            0
134-135      bbox-overlap      578550      481306       423733
134-136     disjoint-bbox           0           0            0
134-137     disjoint-bbox           0           0            0
134-138     disjoint-bbox           0           0            0
135-136      bbox-overlap      635100      543548       455136
135-137     disjoint-bbox           0           0            0
135-138     disjoint-bbox           0           0            0
136-137      bbox-overlap      641625      531829       441072
136-138     disjoint-bbox           0           0            0
137-138      bbox-overlap      669900      528606       447692
```

Going to a **finer** grid cannot create edges (a coarser grid is what merges two lines'
points into one cell), and the missing pairs fail on **extent**, not on cell coincidence.
So 5 m did not miss anything: the 5 m edge set in the brief and the 2 m edge set the code
uses are the same five edges. **The tree verdict is confirmed and is not resolution-dependent.**

Three consequences, all measured:

**(a) The weighted least-squares solve is EXACTLY a chain traversal on this tile.**
`align_swaths` really is a weighted LS network solve **[CODE]** — `coreg.py:568-580`:

```python
    n, E = len(swaths), len(edges)
    A = np.zeros((E, n)); w = np.zeros(E); O = np.zeros((E, 3))
    for e, (a, b, dx, dy, dz, nn, ww) in enumerate(edges):
        A[e, idx[a]] = -1.0; A[e, idx[b]] = 1.0
        w[e] = ww; O[e] = (dx, dy, dz)
    sw = np.sqrt(w)
    corr = np.zeros((n, 3)); mis = np.zeros((E, 3))
    for k in range(3):
        c, *_ = np.linalg.lstsq(A * sw[:, None], O[:, k] * sw, rcond=None)
        c -= c[idx[ref]] if ref is not None else c.mean()   # gauge choice
```

— but on a tree the design matrix is square-equivalent and the fit is exact
interpolation. **[MEASURED]** the solved constants are the running sum of the edge
observations to the digit (intercept tie, mm): `21.987`; `21.987−15.800 = 6.187`;
`6.187−15.983 = −9.797`; `−9.797−8.602 = −18.399`; `−18.399−4.238 = −22.637` — identical
to the reported `+21.987 / +6.187 / −9.797 / −18.399 / −22.637`.

**(b) The edge weights have ZERO effect here.** The carefully reasoned `1/dz_var` weighting
(changed 2026-09-05, `coreg.py:489-505`, `analysis/SWATH_TIE_DEGENERACY.md`) and the whole
degeneracy apparatus it replaced **cannot change any number on this tile**, because a tree
has no redundancy for weights to arbitrate. It matters only where loops exist (Battle Creek:
5 swaths, 7 edges).

**(c) `misclosure` is identically zero, so the method's own quality metric is uninformative.**
**[MEASURED]**, both tie modes: `max |misclosure| in dz: 0.0000 mm` on all five edges.
The docstring says as much (`coreg.py:508-509`: "`~0` for a tree/chain; nonzero only where
redundant overlaps -- loops -- exist to check consistency"). **There is no error control on
this tile of any kind.** Not a small residual — none.

For the record, the full network, both tie modes, at the shipped settings **[MEASURED]**:

```
=== align_swaths(res=2.0, ref=133, tie='overlap_median') ===
      edge      dz_mm      se_mm        n         weight  w_share%  mis_dz_mm
133-134        18.877      0.238   257395      1.762e+07   18.1483     0.0000
134-135       -16.173      0.183   233481      2.983e+07   30.7255    -0.0000
135-136       -15.882      0.223   263447      2.005e+07   20.6566     0.0000
136-137        -6.921      0.241   255988      1.717e+07   17.6859    -0.0000
137-138        -3.540      0.284   234606      1.241e+07   12.7837    -0.0000
dz span across lines: 42.52 mm  (-23.64 .. 18.88)

=== align_swaths(res=2.0, ref=133, tie='intercept') ===       <- what the pipeline ships
      edge      dz_mm      se_mm        n         weight  w_share%  mis_dz_mm
133-134        21.987      0.195   257395      2.626e+07   16.8651    -0.0000
134-135       -15.800      0.141   233481      5.006e+07   32.1541     0.0000
135-136       -15.983      0.177   263447      3.204e+07   20.5792     0.0000
136-137        -8.602      0.195   255988      2.634e+07   16.9198    -0.0000
137-138        -4.238      0.218   234606      2.099e+07   13.4818    -0.0000
dz span across lines: 44.62 mm  (-22.64 .. 21.99)
```

The intercept-tie constants reproduce `data/derived/elbaext_delong/corrections.json`
`per_swath_internal_alignment_dxdydz_m` exactly (`0.022 / 0.0062 / −0.0098 / −0.0184 /
−0.0226` m), so this reproduction *is* the shipped solution.

---

## 2. RANKED OPTIONS, with cost and blast radius

Ranked by expected gain per unit implementation cost. The gain figures are stated as what
they would *target*, not as achieved improvements — none of these has been built.

### Rank 1 — Measure the star tie as a CHECK, and publish the disagreement. (cost: ~1 day; no pipeline change)

Tie every gen1 line **directly to the gen2 surface** and compare with the chain. This costs
nothing in the product and is the only thing on this list that gives the network an
independent second determination **today**. **[MEASURED]**, from
`beam_offset_table.parquet`, per-line median of `d_noswath = d_mm + dz_geoid + dz_lateral +
dz_drift` (the swath term removed), re-expressed on line 133 so it is comparable with the
chain, on the product's own stable mask:

```
=== star observation per line, STABLE cells ===
 line  n_returns   n_cells  median_mm   NMAD_mm  SE_med_mm |  dz_star_mm  dz_chain_mm   diff_mm
  133     430723     34197     -85.92    117.64      0.797 |       -0.00         0.00     -0.00
  134    2141073    155662    -120.46     95.05      0.302 |       34.53        22.00     12.53
  135    2341647    172671    -102.56     96.84      0.292 |       16.64         6.20     10.44
  136    2146197    154707     -72.81    121.51      0.387 |      -13.12        -9.80     -3.32
  137    1515640    114787     -72.70    127.88      0.473 |      -13.23        -18.40      5.17
  138     457333     36251     -62.51    224.86      1.480 |      -23.41        -22.60     -0.81
  spread of dz_star 57.95 mm; spread of dz_chain 44.60 mm; RMS(star-chain) 7.12 mm
```

Upland-only (non-floodplain) cells give `RMS(star−chain) 8.93 mm`, max `+14.64` on line 134.
**The two determinations disagree by 7–9 mm RMS and up to 15 mm on one line.** That is an
order of magnitude larger than either method's formal SE (0.14–0.80 mm) and is the same size
as the residuals the project is chasing.

**What it does NOT tell you, and this is the point.** The disagreement is an *upper bound* on
the chain's error, not a measurement of it: the star tie is contaminated by (i) real change,
(ii) gen2's own vegetation offset (memory: gen2 ground floats above true ground ∝ near-ground
vegetation), (iii) gen1's slope/scan-geometry bias — all of which differ per line because each
line's terrain and cover mix differs. So a 12.5 mm disagreement on line 134 could be the
chain, or could be that line 134's footprint is more forested. **The discriminating test,
which I did not run: does the star−chain disagreement correlate with the line's cover/slope
mix?** If it does, the star tie is measuring cover, not level.

**Cost:** one analysis script; no product rebuild; no parameter. **Blast radius:** none — it
is a diagnostic. **Ask first:** whether to restrict the star tie to open ground only, as the
standing rule for the ground-control datum already does (`datum-use-open-ground-only`).

### Rank 2 — Audit and constrain the per-line HORIZONTAL shift. (cost: ~2-3 days; touches register_gen1)

This is the largest single finding of the study and it is not about the vertical constant at
all. See §4.3. The per-line horizontal shifts `align_swaths` solves are large, they accumulate
**monotonically** along the chain to **1.39 m**, and because a lateral shift acts on sloping
ground through the surface gradient they impose a floodplain-minus-own-upland DoD difference of
**+34.7 mm (line 137)** and **+48.6 mm (line 138)** — while the vertical constants impose
**+1.9 and +0.4 mm**. The horizontal solution is ~20× more consequential for that comparison
than the vertical one, and it has had none of the scrutiny.

Concrete, cheap sub-steps, in order:
1. Report `dx_sigma`/`dy_sigma` (already computed, `coreg.py:167-168`) and the horizontal
   res-sensitivity alongside the constants in `corrections.json`. **[MEASURED]** the formal
   `dx_sigma` is 0.0023–0.0031 m while the res=2 → res=5 change on the same edge reaches
   0.0787 m (edge 134-135) — a factor of ~30, the same pathology `LOCAL_TIE_CHAINING.md` §4
   found for the vertical tie.
2. Weight the horizontal components by the **horizontal** variance. **[CODE]** `coreg.py:575-576`
   solves all three components with the *same* `sw`, derived from `dz_var` — the vertical tie's
   variance. On a tree this is inert; on Battle Creek's loops it is wrong.
3. Decide whether an unverifiable 1.39 m per-line lateral spread should be applied at all.
   **Ask Andy:** this is a parameter-level choice, not an implementation detail.

**Blast radius:** any change to `register_gen1` invalidates every tile's `corrections.json`,
every `beam_offset_table.parquet` (`dz_swath_mm` is derived from it), the `q2_cover_fit`,
and every number in the slope/cover analyses that used `d_mm_corr`. This is the expensive one.

### Rank 3 — Estimate a per-LINE across-track term instead of a per-PAIR one. (cost: ~1 week)

**[CODE]** `across_track_tie` already estimates the across-track slope `c` per edge and
**throws it away** — `coregister_swaths` uses only `at.k`, `at.n`, `at.se` (`coreg.py:788-801`).
`ACROSS_TRACK_TIE_LITERATURE.md` §C already proposed the fix (regress on `tan θ_ref` and
`tan θ_src` as two covariates, giving a per-strip rather than per-pair term) and did not run it.
**[LIT]** this is precisely the field's own parameterisation: Crombaghs et al. (2000) estimate
per strip {offset, along-track tilt, across-track tilt} in one integral adjustment
(https://www.isprs.org/proceedings/Xxxiii/congress/part3/230_XXXIII-part3.pdf).

**[MEASURED]** the size of the term being ignored. Refitting each edge with the library's
existing joint translation+tilt model (`coreg.tie_translation_tilt`, unmodified, its own
defaults) at res=2.0. Flight lines here run N-S (y extent 4350 m, x extent ~1500 m per line),
so `c1 = dz/dEast` **is** the across-track tilt:

```
      pair |    NK dx    NK dy  NK dz_mm  NK nmad | T+tilt dx T+tilt dy    c0_mm  c1 mm/km E  c2 mm/km N     nmad   conv
133-134    |   0.2543   0.1355     18.88    120.9 |    0.2246    0.1421    20.32      -162.0        -9.1    118.5   True
134-135    |   0.0787   0.0458    -16.17     88.5 |    0.0748    0.0850   -17.57        44.2        -3.0     86.8   True
135-136    |   0.3220  -0.1068    -15.88    114.6 |    0.3189   -0.1112   -12.34       -89.8         0.2    113.7   True
136-137    |   0.2663  -0.0848     -6.92    122.1 |    0.2889   -0.1124   -15.94       142.8         8.8    119.3   True
137-138    |   0.4676   0.2939     -3.54    137.5 |    0.4479    0.2690    -3.85       -49.5         1.8    137.9   True
```

Two readings. (i) Allowing a tilt moves the **vertical constant** by up to **9.0 mm**
(edge 136-137: `−6.92 → −15.94`), against that edge's reported formal `se = 0.195 mm`. The
constant is model-dependent at the 9 mm level and its error bar says 0.2 mm. (ii) The
across-track tilts are ±44 to ±162 mm/km and **alternate in sign** — over a ~300 m sidelap
that is a ±13 to ±49 mm vertical ramp inside the very overlap the constant is read from.
**[INFERENCE]** sign-alternation across adjacent pairs is what a *per-line* term looks like
when you fit it *per pair* (each pair sees the difference of two lines' terms).

**Blast radius:** as Rank 2, plus it changes `across_track_tie`'s contract.
**Expected gain:** it is the natural target for the 5.07 mm residual mean |per-line step|
that `align_swaths` leaves in shared overlap cells (`pipeline.py:45`).

### Rank 4 — Tie every line to gen2 in ONE simultaneous adjustment (a star, gauged on gen2). (cost: 2+ weeks; CONTRAINDICATED as specified)

This is the option the brief asks about directly. Structurally it does everything one wants:
6 star edges + 5 line-to-line edges = 11 observations on 5 free unknowns → **6 degrees of
freedom, real misclosure, no chain, no zero line, no gauge problem.** **[LIT]** it is
established practice, not an invention: Glira et al. (2015) list the datum options as fixing
one or more strips' trajectories **or** "ground-truth data … as a DEM to which the ALS block
should fit" (https://isprs-annals.copernicus.org/articles/II-3-W5/73/2015/isprsannals-II-3-W5-73-2015.pdf);
Wimmer et al. (2024) name and systematically evaluate exactly this design as **bi-temporal
adjustment against a reference block (REF)** and find it beats both control-patches-only and
all-epochs-simultaneous (https://repositum.tuwien.at/bitstream/20.500.12708/200265/1/Wimmer-2024-Remote%20Sensing-vor.pdf).
Crombaghs et al. (2000) reached for it for our exact reason: "From the height differences of
tie 'points' in the overlapping areas, it is not possible to say whether the additional
systematic errors originate from the left or right strip, or even from both strips."

**Why it is nonetheless ranked last here, and why I am not recommending it.** It inverts this
project's decided reference. Memory `reference-epoch-is-gen1`: Andy, 2026-08-27 — "trust gen1,
it was leaf-off"; gen2 is the epoch that gets corrected, because **gen2 reads high in
vegetation by 73.4 mm (measured)** and **gen2's own absolute level is unresolved (−5 ± 24 mm)**.
Gauging gen1's flight lines on gen2 would:
- **Absorb real change into the levels.** Each line's constant would soak up the mean gen2−gen1
  difference over *that line's footprint*, and the footprints are 1.5 km-wide N-S strips that
  sample different amounts of floodplain. Line 138 is 69.5% of the floodplain; line 134 is 0.02%.
  A real floodplain signal becomes a per-line constant.
- **Alias gen2's vegetation offset into gen1's geometry.** Lines differ in cover mix, so they
  would differ in inherited gen2 veg bias.
- **Re-introduce a circularity the pipeline deliberately removed.** **[CODE]**
  `SWATH_ALIGNMENT_METHOD.md`: "It is gen1-internal through step 6: no gen2 data is read and no
  absolute level is set. That separation is what stops the swath tie absorbing a cross-epoch
  correction, and a test asserts `register_gen1`'s signature carries no gen2 argument."
  Wimmer et al. (2024) hit the same wall from the other side and answered it with an explicit
  stability trim — their reference block is "trimmed to stable areas". That is the only version
  of this option I would consider: star edges computed **on stable, open ground only**, with the
  gen1-internal line-to-line edges kept so the two can be cross-checked rather than replaced.
- **Removes the gauge problem that is already removed.** `ground_control/apply_datum.py` makes
  the product gauge-invariant *exactly* (spread across zero-line choices: 42.40 mm uncorrected
  → `< 1e-9` mm corrected, `ground_control/tests/test_ground_control_apply_datum.py`). The
  zero-line problem is solved. Redundancy is the thing that is not.

**Net:** build the star as a **check** (Rank 1), not as the **gauge** (Rank 4). The redundancy
is the prize; the re-gauging is the risk.

### Not recommended — a per-line time-dependent term

`coreg.fit_along_track_drift` exists and is a per-swath `f(gps_time)` spline, the standard
strip-adjustment trajectory form. It is OFF by default and the reason is measured, not
stylistic (`pipeline.py:39-58`): out-of-sample skill 0.160/0.142/0.087 against the correction
surface's 0.444/0.384/0.283, and it **makes per-line agreement worse** — align_swaths brings
the mean |per-line step| in shared overlap cells from 15.15 to 5.07 mm and the drift returns
it to 19.91/16.97/21.46 mm. **[LIT]** Glira et al. (2015) reached the same restraint from
theory: they limit trajectory correction to the constant part per strip because
"the determinability of such parameters strongly depends on the terrain geometry and can not
be guaranteed in any case (e.g. over flat terrain)." Adding time-dependent terms to a network
with **zero** redundancy would be worse still.

---

## 3. What `align_swaths` actually estimates (question 1)

**[CODE], all of it.**

| | |
|---|---|
| Parameters per flight line | **exactly three**: `(Δx, Δy, Δz)`, a rigid 3-D translation. `corrections = {swath: (Δx, Δy, Δz)}`, `coreg.py:580` |
| Tilt per line | **none**. Not estimated, not available in `align_swaths`. |
| Time-dependent term per line | **none** in `align_swaths`. `coreg.fit_along_track_drift` is a separate, later, OFF-by-default stage in `apply_datum`. |
| Scan-angle term | **one COMMON term for all lines**, not per line: `coreg.estimate_boresight_roll` → `z -= b·scan_angle`, applied in `register_gen1` *before* the network, opt-in, and `boresight_roll_mm_per_deg: null` on this product. |
| Per-pair across-track slope | **estimated and discarded**. `across_track_tie` returns `ATTie(k, c, n, se)`; `coregister_swaths` keeps `k`, `n`, `se` and drops `c` (`coreg.py:788-801`). |
| Observation per edge | `c_b − c_a = s_ab` where `s_ab` aligns b onto a (`coreg.py:571`: `A[e, idx[a]] = -1.0; A[e, idx[b]] = 1.0`) |
| Solve | weighted least squares, `np.linalg.lstsq(A*sqrt(w), O[:,k]*sqrt(w))`, **the same `w` for all three components** (`coreg.py:575-576`) |
| Weight | `w = 1/dz_var` if finite and positive else `0.0` (`coreg.py:557`) — the **vertical** tie's variance, used for the horizontal components too |
| Gauge | `c -= c[idx[ref]]` if `ref` given, else `c -= c.mean()` (`coreg.py:577`) |
| Filtering | non-finite observations dropped (`coreg.py:528`); **nothing else** — the long comment at `coreg.py:530-548` records four candidate rules and why each was rejected |
| Grid the edges are built on | per-cell **median-Z** at `res` over the pair's bbox intersection, `exclude=(5,6,9)`. **[MEASURED]** on this tile `exclude` is a no-op: `data/csf_cache/elbaext.las` is `{2: 15581591}` — 100% class 2 — so the grids are ground medians. |

**`tie="intercept"` vs `"overlap_median"`.** **[CODE]** The horizontal solution is identical
either way (the Nuth & Kääb aspect fit removes `median(dh)` before fitting, so it never sees
the vertical constant) — confirmed **[MEASURED]**: `dx`/`dy` are bit-identical in the two
network runs in §1. Only `dz` and its variance change:
- `overlap_median` → `dz = median(z_ref − z_src)`, `dz_var = (nmad_after/√n)²`.
- `intercept` → the LAD (median-regression) intercept at across-track position `dtan = 0`,
  `dz_var = se²` with `se = scale·√(1/n + mean(t)²/Σ(t−mean(t))²)` — i.e. it **carries the
  extrapolation lever arm**, which is the whole point (`coreg.py:716-728`).

**[MEASURED]** the difference costs up to 3.4 mm on a constant (line 136: `−13.178` median vs
`−9.797` intercept) and 2.1 mm on the span (42.52 vs 44.62 mm).

**Is it least squares or a chain traversal?** It is genuinely least squares **[CODE]**, and on
this tile it is **numerically identical to a chain traversal** **[MEASURED]** (§1a). Both
statements are true and neither is the misleading one to omit.

---

## 4. Does the tree explain the offsets? (question 3) — a NO, then a different YES

### 4.1 The formal edge SEs say NO, by two orders of magnitude

**[MEASURED]** The five edge SEs (intercept tie) are `0.195, 0.141, 0.177, 0.195, 0.218` mm.
On a chain the level errors add in quadrature if independent:
`√(0.195² + 0.141² + 0.177² + 0.195² + 0.218²) = 0.42 mm` at the far end (line 138, five links
from the gauge). Against 25–28 mm, that is **~60× too small**. Taken at face value, the tree
cannot carry the observed offsets.

### 4.2 But the formal SE is measured to be wrong, and then the tree CAN carry it

The formal SE is the standard error of a median over ~250,000 cells. It is not wrong as
arithmetic; it is the wrong estimand. Three independent measurements of the honest per-link
error, two of them pre-existing and one from this session:

| source | what was varied | per-link scatter |
|---|---|---|
| `LOCAL_TIE_CHAINING.md` §5 (read, not redone) | window POSITION, 9 windows inside one tile, 400 m half-width | **6.0 mm** (median tie), **10.7 mm** (intercept) against formal 0.4 mm |
| `LOCAL_TIE_CHAINING.md` §3 | tile EXTENT (elba vs elbaext, same lines) | **8.42 mm RMS per link** |
| **this session [MEASURED]** | grid RESOLUTION, res 2.0 → 5.0, same code | `dz` moves `0.30, 1.80, 0.85, 0.33, 0.93` mm; `dx` moves up to `0.079 m` |

Take the in-tile window figure, 6.0–10.7 mm per link: five links in quadrature gives
**13.4–23.9 mm** at line 138. **That is the size of the observed offsets.** So the tree *does*
have the capacity to deliver 10–25 mm of level error to the far end — and it is verified
**by measurement of the estimator's repeatability**, not by its formal error bar.

**[INFERENCE], and it makes the chain worse than quadrature.** Quadrature assumes the five
edge errors are independent. They are not: **[MEASURED]** every edge on this tile is the same
geometry — the src line is always the eastern one (per-line x ranges increase monotonically
with line number), so each sidelap is sampled at high positive scan angle in the ref line and
high negative in the src line. Identically-structured observations share their systematic
error, and a chain **adds** correlated errors linearly instead of averaging them down. The
smoking gun is in the horizontal solution (§4.3): **all five `dx` observations have the same
sign.**

### 4.3 The offsets in the brief, however, are ALGEBRAICALLY IMMUNE to the vertical constants — and the lateral shift explains them instead

This is the finding I would put in front of everything else.

**The argument [INFERENCE, then verified].** A per-line vertical constant `c_L` shifts *every*
cell line L touches by the same amount. So it shifts line L's floodplain mean and line L's own
upland mean **equally**, and their **difference is invariant to `c_L`**. No matter how wrong the
chain's constants are, and no matter the topology, they cannot produce a per-line
floodplain-minus-own-upland difference. The only residual channel is that cells are *shared*,
so a line's floodplain cells and its upland cells can mix the *other* lines' constants in
different proportions.

**[MEASURED]** that residual channel, and it is negligible. Per line, the imposed DoD shift
(`= −v`, `v = 1000·(dz − gx·dx − gy·dy)`, vertical, return-weighted per cell, from
`corrections.json` + `z_after.npy` + the beam table's per-cell line membership), split into the
**constant** part `1000·dz` and the **lateral** part `−1000·(gx·dx + gy·dy)`:

```
 line |  const fp  const up     diff |   lat fp   lat up     diff | TOTAL diff
  133 |         -   -11.696        - |        -   -0.878        - |          -
  134 |   -22.000   -16.848   -5.152 |    2.195    0.557    1.638 |     -3.514
  135 |     1.359    -6.039    7.397 |    3.825    1.739    2.085 |      9.483
  136 |    11.941     8.465    3.475 |  -13.945  -14.270    0.324 |      3.800
  137 |    19.142    17.197    1.945 |   -9.821  -44.474   34.653 |     36.597
  138 |    21.464    21.084    0.381 |  -19.270  -67.893   48.624 |     49.004
```

And confirmed independently on the per-return table, by adding each component back to the
swath-free offset and re-measuring the floodplain-minus-upland difference
(`elbaext_independent`, the surface-free route the brief's numbers come from; the delong route
gives the same picture):

```
 line  const contrib  lateral contrib  total contrib
  134           0.02            -2.64          -2.38
  135           0.09            -1.98          -1.88
  136          -0.19             2.02           1.78
  137          -0.45            -9.32          -9.84
  138          -0.77           -21.92         -22.89
```

**The vertical constants contribute ≤ 0.8 mm. The per-line LATERAL shifts contribute 10–49 mm,
on the two lines that carry the floodplain.** The mechanism is not subtle: the floodplain is
flat and the upland is not — **[MEASURED]** median gen2 |gradient| is `0.0233` over floodplain
cells and `0.1535` over the rest, a factor of 6.6. A lateral shift changes a line's reported
elevation by `−(gx·dx + gy·dy)`, so it does almost nothing on the floodplain and a great deal
on the hillslope. Line 137's mean east-gradient over its upland cells works out to ≈ −0.048
(terrain falling east), and 0.9213 m east × 0.048 = **44 mm** of elevation change on the
upland alone.

**And the lateral shifts are the part of the solution with no redundancy, no sanity bound, and
a signature of systematic error.** **[MEASURED]**, edge by edge at res=2.0:

```
    133->134: dx +0.2543 m  dy +0.1355 m
    134->135: dx +0.0787 m  dy +0.0458 m
    135->136: dx +0.3220 m  dy -0.1068 m
    136->137: dx +0.2663 m  dy -0.0848 m
    137->138: dx +0.4676 m  dy +0.2939 m
```

**All five `dx` observations are positive** (`dy` alternates: `+,+,−,−,+`). Flight lines run
N-S, so east *is* the across-track direction. The chain integrates the five same-signed
increments into a **1.3889 m** spread across six lines (res=5.0 gives 1.2393 m, so it is not a
resolution artifact). Under a null of independent random signs, 5/5 is p = 0.0625 two-sided by
sign test — suggestive, not decisive on 5 observations, but the *accumulation* is what does the
damage and it does not care about the p-value.

**Two hypotheses for the eastward bias, one of which I tested and REFUTED.**
- *Refuted* **[MEASURED]**: it is not an across-track vertical tilt aliased into a horizontal
  shift. Refitting with `tie_translation_tilt` (translation **plus** planar vertical tilt) barely
  moves `dx`: `0.2543→0.2246`, `0.0787→0.0748`, `0.3220→0.3189`, `0.2663→0.2889`,
  `0.4676→0.4479` — max change 0.030 m against a 1.39 m accumulation. The horizontal signal
  survives the tilt. That is a clean negative and it kills the tidiest explanation.
- *Untested* **[INFERENCE]**: a scan-angle-dependent horizontal error (mirror-angle scale /
  boresight heading or pitch). Because every sidelap is viewed at + scan by the ref line and
  − scan by the src line, such an error produces opposite across-track displacements in the two
  lines and therefore the **same-signed** apparent relative shift on every edge. **[LIT]** this
  is exactly the term Glennie et al. (2014) had to add to repair a legacy survey — "a
  second-order polynomial in scan angle … for scanning-mirror/laser-beam misalignment, estimated
  simultaneously with boresight" — and Filin (2003) lists "a scan angle error λω_s that measures
  inaccurate scan angle determination that vary linearly" as a standard system parameter.
  I have not tested it; doing so is Rank 2 above.

**Caution on how much of the observed difference is thereby "explained".** The lateral
contribution (+34.7 / +48.6 mm imposed) is the same sign and larger than the brief's residuals
(+28 / +25 mm), but I could not reproduce the brief's per-line table (§6), so I cannot close
the arithmetic. What is solid is the **decomposition**: whatever the residual turns out to be,
the constants are not carrying it and the lateral shifts are carrying tens of millimetres.

### 4.4 The other thing the tree does: it biases the POOLED comparison

**[MEASURED]** Line composition is wildly different between the two populations — the
floodplain is the far end of the chain:

```
 line    dz_mm   fp cells    fp %   up cells    up %
  133     0.00          0    0.00      44889    7.32
  134    22.00         15    0.02     197898   32.27
  135     6.20       1877    2.05     224058   36.54
  136    -9.80       9574   10.45     218879   35.70
  137   -18.40      61479   67.09     166405   27.14
  138   -22.60      63643   69.46      55027    8.97
```

So although the constants are invariant *per line*, the **pooled** floodplain mean and the
**pooled** upland mean sit on different mixtures of them. **[MEASURED]** imposed pooled DoD
shift, floodplain minus upland: **constant part +19.385 mm**, lateral part −4.608 mm,
**total +14.777 mm**. (Statistic: each cell's value is the return-weighted mean of the
contributing lines' imposed shift; the population statistic is then the unweighted mean over
cells. Floodplain 91,631 cells, upland 613,170 cells.) **If any floodplain-vs-upland
comparison is made on pooled cells rather than within lines, ~19 mm of it is the swath network's
constants, not the ground.** That is a real, quantified, avoidable artifact.

**[MEASURED]** and confirming the brief: overlap availability is not the problem. Floodplain
cells seen by ≥2 lines **49.06%** (mean 1.491 lines), upland **47.95%** (mean 1.479) —
reproducing the brief's 49.06% exactly.

---

## 5. What the established methods do that this does not (question 4)

**[LIT]** — every source below was retrieved (PDF or page fetched) by a subagent in this
session; URLs given so any claim can be checked. I have read the quoted extracts, not the full
papers. Verification status is stated per item. The repo already holds a partial version of
this in `analysis/ACROSS_TRACK_TIE_LITERATURE.md` (2026-09-05), which explicitly flags its
section C as unverified; this extends it and verifies more of it.

### 5.1 The field's own taxonomy

Glira, Pfeifer, Briese & Ressl (2015), *ISPRS Annals* II-3/W5, 73–80
(https://isprs-annals.copernicus.org/articles/II-3-W5/73/2015/isprsannals-II-3-W5-73-2015.pdf,
**full text fetched**) states it cleanly: "Basically, two types of strip adjustment methods
exist: (a) rubber-sheeting coregistration solutions, which use the measured terrain points only
as input … and (b) rigorous solutions that start from the original scanner and trajectory
measurements." **`align_swaths` is type (a)** — and it must be, because the 2008 vendor supplied
no raw trajectory. Habib et al. (2008) is the explicit statement of that constraint: rigorous
calibration "demands the original observations (GPS, IMU and the laser measurements), which are
not usually available to the end-user" (https://www.isprs.org/proceedings/xxxvii/congress/1_pdf/64.pdf,
**full text fetched**).

### 5.2 Parameters per strip: everyone estimates more than three, and the tilts are the standard ones

| method | per-strip parameters |
|---|---|
| Crombaghs et al. (2000), the canonical data-driven vertical adjustment | **3, all vertical**: `ΔH(U,V) = a + bU + cV` — offset, **along-track tilt**, **across-track tilt**. Plus optional cross-strip parabolic effects, strip torsions, along-track periodic terms. |
| Habib et al. (2008), data-driven | **6**: three shifts + three rotations (rigid body) |
| Glira et al. (2015), rigorous | **6 per strip** (3 attitude + 3 position trajectory corrections) plus **12 block-wide** (range/angle offsets and scales, boresight misalignment, lever arm) = "12+6n parameters"; 627 parameters for 103 strips |
| Filin (2003) | mounting bias, range bias, **linearly varying scan-angle error**, GPS offset **and drift**, INS offset **and drift** (https://www.isprs.org/proceedings/XXXIV/3-W13/papers/Filin_ALSDD2003.pdf, **full text extracted; URL confirmed by search result**) |
| Wimmer et al. (2024), multi-temporal | scanner (constant + linear range/direction), mounting, **trajectory as "bias, linear strip-wise correction models or time-dependent corrections within a flight line using splines of arbitrary length"**, plus an explicit per-session datum shift (X,Y,Z) |

**`align_swaths` estimates 3 parameters per line, of which only 1 is vertical.** Crombaghs'
three are *all* vertical. The across-track tilt in particular is not an optional refinement in
the literature — it is one of the three baseline parameters, and Crombaghs reports its
"standard deviation of 33 cm/km" across strips, which brackets the ±44 to ±162 mm/km I measured
per pair in §2 Rank 3.

**[LIT]** Two restraints the field states that this project should keep: Glira deliberately
drops time-dependent trajectory terms because "the determinability of such parameters strongly
depends on the terrain geometry and can not be guaranteed in any case (e.g. over flat terrain)";
and Skaloud & Lichti (2006) found their estimated range bias `Δρ = −0.158 m` was not real,
being correlated ~0.4 with every boresight angle and limited by 1-dm GPS height bias
(https://infoscience.epfl.ch/server/api/core/bitstreams/cc1b0b2a-7fcd-410a-9d8c-08e97ee08d38/content,
**full text fetched**). Adding parameters to an unredundant network is how you get a confident
wrong number.

### 5.3 Observations: nobody uses a bare overlap median

Ranked by prevalence in the retrieved sources: (i) **point-to-plane distances on ICP
correspondences** between the original points of overlapping strips, re-established iteratively
(Glira 2015: 345,205 strip-to-strip correspondences; roughness gate σ_max = 10 cm; outliers
"typically does not exceed 5%"); (ii) **plane/patch conditioning**, where a group of points is
forced onto a common plane whose 4 parameters are co-estimated (Skaloud & Lichti 2006, Filin
2003); (iii) **conjugate linear features** from intersected planar patches (Habib 2008);
(iv) **averaged height differences over flat 50×50 m tie areas of ≥100 points** — Crombaghs
2000 explicitly rejects individual points as too noisy ("point noise of about 10–15 cm";
averaging 100 points gives 1.7 cm); (v) **ground control** as points, planes, TLS clouds, or a
DEM; (vi) **GNSS/IMU trajectory as observations with their own stochastic model**, in the
rigorous formulations.

The repo's estimator is closest to (iv). **[LIT]** Skaloud & Lichti's geometric requirement is
the one most directly relevant to our narrow N-S sidelaps: "to achieve sufficient de-correlation
between the plane and system parameters it is important that the planes vary in slope and
orientation," and horizontal planes "contribute very little to the estimates of bore-sight
angles in pitch and heading direction."

### 5.4 Redundancy and gauge: simultaneous, over-determined, and with cross-strips

This is the sharpest contrast with what we have.

- **Simultaneous, not chained.** Crombaghs: "Execute an integral least squares adjustment."
  Glira: one adjustment over all 748 overlapping strip pairs plus ground-truth correspondences.
  `align_swaths` *is* a simultaneous solve by construction — but with 5 observations on 5
  unknowns it has nothing to adjudicate.
- **An explicit redundancy rule.** Crombaghs: "The number of height differences from step 1 and
  step 2 together has to exceed **three times the number of strips** (per strip three parameters
  are estimated), to ensure redundancy," and tie points "must not lie on a straight line in
  order to avoid singularities." **Our tile: 5 observations, 6 strips. The rule asks for > 18.**
  We are at 28% of the *minimum* for a 3-parameter-per-strip model, and our five tie areas
  *are* collinear — five parallel N-S sidelaps.
- **Cross-strips are a determinability requirement, not a nicety.** Filin: "not all errors can
  be recovered under any flight configuration … to recover this error flight patterns that
  include cross strip(s) are required." Crombaghs uses them as deformation diagnostics and
  reports results "improve dramatically when cross strips are added." Glira's block: 95
  longitudinal + **8 cross strips, flown in both directions**. USGS institutionalizes it: the
  interswath check must cover "cross-tie swaths and a sample of intersecting project swaths **in
  both flight directions**."
  **We have no cross-strips, and 2008 MN DNR did not fly them for us to find.** That is a
  permanent, unfixable limitation of this dataset, and it is the reason an external reference
  surface (§2) is the only available substitute.
- **Datum/gauge.** Three verified options: fixing one or more strips' trajectories (our
  `ref=zero_line`, a minimum-constraint free network — legitimate and standard); external ground
  truth (points, planes, TLS, **or a DEM**); or an explicit per-session X/Y/Z datum-shift
  parameter (Wimmer 2024, specifically so reference-station error cannot leak into the other
  parameters). **Our gauge choice is defensible and the project has already made it exactly
  cancellable** via `apply_datum`. The gauge is not the problem.
- **Regularization when the normal matrix goes near-singular:** pseudo-observations pulling
  corrections toward zero with a stated σ (Wimmer 2024). Worth knowing before adding parameters.

### 5.5 Accuracies and quality metrics

**[LIT]** Post-adjustment residuals cluster at a few cm: 2–3 cm residual std.dev. and 3–5 cm
tie-area RMS, from 6–10 cm before (Crombaghs 2000; block histogram 7.1 → 4.6 cm); 5.3–5.4 cm
std.dev. of point-to-plane residuals in mountainous terrain (Glira 2015, 103 strips);
1–2 cm relative between epoch DSMs and 3–5 cm against independent checkpoints, "below 2–3 cm
for all involved epochs with a mean of 1.7 cm" (Wimmer 2024). Crombaghs' simulation says 1 cm
mean height over 0.01 km² is achievable **only with control in every strip**, degrading to
9.5 cm with control every 49th strip.

**Our 5.07 mm mean |per-line step| in shared overlap cells (`pipeline.py:45`) is better than
anything in that list**, and 5–10× better than USGS's acceptance limit. That deserves saying
plainly: the shipped method is not performing badly against published practice on its own
metric. What it lacks is not precision but **checkability**.

Standard metrics **[LIT]**: strip-difference (misclosure) rasters restricted to smooth,
low-slope, single-return, non-extrapolated cells, summarized by RMSDz or robust σ_MAD; percent
of cells over a client threshold; sliding-window least-squares matching to split the misclosure
into X/Y/Z (Ressl, Kager & Mandlburger 2008,
https://www.isprs.org/proceedings/XXXVII/congress/3b_pdf/54.pdf, **full text extracted**).
Note the last one: the field's standard QC **already** splits misclosure into horizontal and
vertical, which is the diagnostic §4.3 says we need.

**USGS Lidar Base Specification** (v2.1 PDF, **fetched**:
https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/atoms/files/Lidar-Base-Specification-version-2-1.pdf;
2022 rev. A copy **fetched**, tables identical):
- Table 2, relative vertical accuracy for swath data: swath-overlap-difference **RMSDz ≤ 0.04 m
  (QL0), ≤ 0.08 m (QL1/QL2), ≤ 0.16 m (QL3)**; smooth-surface repeatability ≤0.03/0.06/0.06/0.12 m.
- Table 4, absolute: RMSEz(non-veg) ≤0.100 m and NVA ≤0.196 m at 95% for QL1/QL2; VVA ≤0.30 m
  at the 95th percentile.
- Interswath consistency is assessed "in nonvegetated areas of only single returns and with
  slopes of less than 10 degrees"; current USGS page (**fetched**,
  https://www.usgs.gov/ngp-standards-and-specifications/interswath-consistency-check-ssis):
  "Δz ≤ 8 cm (green) in areas of swath overlap shall be considered acceptable", unacceptable if
  "Δz is greater than 16 cm, is correlated with the collection flightlines, contiguous and
  prevalent."
- Sampling instruction: "test areas should be located such that the full width of the overlap is
  represented" — which `ACROSS_TRACK_TIE_LITERATURE.md` §B3 already reads, correctly, as the
  spec conceding the extent-dependence and handling it by sampling design rather than modelling.

### 5.6 The direct precedent for adjusting an old survey against a new one

**[LIT]** This exists and is not marginal. Three precedents, all retrieved:
1. **Wimmer, Mandlburger, Ressl & Pfeifer (2024)**, *Remote Sensing* 16(15):2838
   (https://repositum.tuwien.at/bitstream/20.500.12708/200265/1/Wimmer-2024-Remote%20Sensing-vor.pdf,
   **full text fetched**): 21 epochs over 8 years; four strategies compared; **REF**
   (bi-temporal adjustment of each epoch against a whole reference block) **won**, beating both
   control-patches-only and the all-epochs-simultaneous FULL, which suffered "a vast number of
   highly correlated observations … near-singular equation system." Their four criteria for the
   reference block — **Coverage, Stability, Accuracy, Comparability** — are directly applicable,
   and their reference block is "trimmed to stable areas." Two counterintuitive findings worth
   carrying: temporal distance from the reference block did **not** degrade accuracy, and the
   reference need **not** be the densest epoch, only dense enough for reliable local normals.
2. **Glira et al. (2015)** already admits "a DEM to which the ALS block should fit" as a
   first-class datum option, so tying strips to an external reference surface is inside the
   standard model, not an extension of it.
3. **Glennie et al. (2014)**, *GRL* 41 (https://web.uvic.ca/~enissen/papers/glennie-etal-2014-GRL.pdf,
   **full text fetched**) is the methodological sibling of our problem: a sparse 2006 legacy
   survey (0.013 pts/m², ±32°, 6000 m AGL) repaired because differencing against a good 2010
   survey exposed flight-line-correlated artifacts. Their sequence: reprocess DGPS → fix the
   ITRF epoch → re-estimate boresight from flight-line overlap → add a second-order polynomial
   in scan angle. They also "cropped the pre-event point cloud data set to remove points at
   higher scan angles in areas of overlap," and note residual errors "are mainly at the edge of
   the flight line swaths (i.e., at larger scan angles)."

**No settled name for the technique.** The verified vocabulary is "reference block",
"ground-truth data", "bi-temporal / multi-temporal strip adjustment", and at the DEM level
"co-registration" (Nuth & Kääb 2011). "DEM-based strip adjustment" appeared only in an
unretrieved search snippet — **NOT VERIFIED** as a term of art.

**NOT VERIFIED, do not cite from this document:** Kilian/Ackermann/Englich (1996) — conflicting
author lists in snippets, paper never retrieved; Habib/Kersting/Bang "quasi-rigorous" strip
adjustment — publisher returned 403, snippet only; Friess, Burman, Kager (2004), Toth, Hebel &
Stilla, Kersting et al. (2012), Latypov, Glennie (2012) — present only as citations *inside*
sources that were read.

---

## 6. Things I could NOT reproduce, and one documentation contradiction

**Report these as negatives; do not build on my versions of them either.**

### 6.1 I could not reproduce the per-line floodplain/upland table in the brief

The brief gives, on `elbaext_independent`: 135 fp −43 / upland +10 / diff −52 (n=711);
136 +5/+3/+2 (1,409); 137 +34/+6/+28 (10,828); 138 +35/+11/+25 (9,354).

**[MEASURED]** from `elbaext_independent/dod.npy` (the product DoD, ×1000 mm) with per-cell
line membership from its own beam table, over cells each line contributes any in-grid return to:

```
 line |          ANY-return cells, MEAN         |        ANY-return cells, MEDIAN
      |        fp    upland     diff   n_fp     |       fp       up     diff     n_fp
  134 |     -66.9       7.3    -74.2     11     |    -86.1      6.3    -92.4       11
  135 |     -46.7       5.7    -52.4   1877     |    -45.0      5.7    -50.7     1877
  136 |      20.0      -2.6     22.6   9553     |      2.7     -1.2      3.8     9553
  137 |     -60.0      -0.7    -59.3  60909     |      2.3      0.6      1.7    60909
  138 |     -78.8       4.4    -83.3  63087     |    -10.0      1.1    -11.1    63087
```

Line 135 reproduces well (median: −45.0/+5.7/−50.7 against −43/+10/−52) and line 136 is close
(+2.7/−1.2/+3.8 against +5/+3/+2). **Lines 137 and 138 do not reproduce under any of the seven
cell-set/statistic combinations I tried** (any-return, majority, exclusive, each × with/without
`stable.npy`, mean and median): my floodplain values for 137/138 are near zero or negative where
the brief has +34/+35, and my floodplain cell counts are 5–6× larger (60,909 vs 10,828;
63,087 vs 9,354). Note also **[MEASURED]** `stable & floodplain == 0 cells` — the product's
stable mask excludes the floodplain entirely, so the brief's floodplain numbers cannot be on
stable ground.

I could not find the producing definition, and I am not going to guess a mask. **What I need:
the script or the exact cell-set definition behind that table.** Until then, treat those four
numbers as unverified by me — but note that **none of §1–§5 depends on them**; §4.3's
decomposition (constants ≤0.8 mm, lateral 10–49 mm) is measured independently of that table.

One relevant caution from the above: for line 137's floodplain the **mean is −60.0 mm and the
median +2.3 mm**. That distribution is heavily left-skewed, so mean-vs-median is worth tens of
millimetres there and any floodplain number must say which it is.

### 6.2 The documented 42.40 mm span is the span of the tie mode the pipeline does NOT use

`pipeline.py:595` (and :538, :1014, :1315) and `corrections.json`'s own `absolute_datum_note`
all say the six per-swath dz span **42.40 mm (−23.60 .. +18.80)**. **[MEASURED]** those are the
`tie="overlap_median"` constants (`−23.6389 .. +18.8767`, span 42.52 mm, truncated to 4 dp in
metres → 42.40). The pipeline ships `tie="intercept"`, and the constants actually written into
the same `corrections.json` are `−0.0226 .. +0.0220` m — **span 44.60 mm**. So the file's note
and the file's numbers disagree by 2.2 mm, and `pipeline.py`'s comment describes a
non-default configuration. Small, but it is a docs-vs-behaviour drift in exactly the place the
gauge argument is made. Not fixed here (study only).

### 6.3 What I did not measure, and would need to before any of §2 is built

- Whether the 1.39 m eastward `dx` accumulation is a **scan-angle-dependent horizontal error**.
  The tilt hypothesis is refuted; the scan-angle hypothesis is untested.
- Whether the star−chain disagreement (7–9 mm RMS) is the chain or is gen2's cover-dependent
  bias. The test is whether it correlates with each line's cover/slope mix.
- Whether a **per-line** (rather than per-pair) across-track term is even identifiable on five
  collinear parallel sidelaps. Glira's warning applies verbatim: "depending on the assembly of
  the sensors, the flight configuration, or the terrain geometry, some of these parameters may be
  completely correlated and therefore not estimable."
- Anything about Battle Creek, Carlton, Cook or mnrv. **Everything here is elbaext.** Battle
  Creek has 7 edges on 5 nodes and therefore loops, so its topology verdict is different and its
  weighting scheme is live. This study says nothing about it.

---

## 7. Summary judgement

**The vertical part of the method is defensible and is performing well against published
practice.** 5.07 mm mean |per-line step| in shared overlap cells beats every post-adjustment
figure in §5.5 and is 5–10× inside USGS's 8 cm acceptance limit. The gauge choice is standard
(minimum-constraint free network) and the project has made it *exactly* cancellable. The
weighted-LS network, the `1/dz_var` weighting, and the intercept tie are all well-reasoned.
I am not recommending they be replaced.

**What is genuinely missing is redundancy, and the tile cannot supply it.** Five collinear
parallel sidelaps, 6 nodes, 5 edges, zero loops, misclosure identically zero, no cross-strips
flown in 2008, and Crombaghs' minimum-redundancy rule asking for >18 observations where we have
5. There is no internal check available, at any effort. The only substitute this project has is
the gen2 surface, and using it as a **check** (Rank 1) costs one script and no product change.

**And the thing nobody was looking at is the horizontal solution.** It carries a 1.39 m per-line
spread built from five same-signed increments, its formal σ (2.5 mm) understates its
res-sensitivity (79 mm) by ~30×, it is weighted in the network by the *vertical* tie's variance,
and it imposes a floodplain-minus-own-upland DoD difference of +34.7 and +48.6 mm on the two
lines that carry the floodplain — against +1.9 and +0.4 mm from the vertical constants everyone
has been worrying about. Whether those shifts are right or wrong, they are ~20× more
consequential for the comparison in the brief than the levels are, and they have had none of the
scrutiny.

---

## 8. Commands run (all of it, in order)

```
./lidar-icp/bin/python <scratchpad>/edgeset.py            # §1: real edge set at res=2.0, both tie modes
./lidar-icp/bin/python <scratchpad>/pairdetail.py         # §2,§4: per-edge Coreg fields; res 2.0 vs 5.0
./lidar-icp/bin/python <scratchpad>/tilt_vs_shift.py      # §2,§4.3: translation+tilt refit per edge
./lidar-icp/bin/python <scratchpad>/fp_composition.py     # line composition + scan geometry of fp vs upland
./lidar-icp/bin/python <scratchpad>/imposed_bias.py       # §4.4: pooled imposed bias, constant vs lateral
./lidar-icp/bin/python <scratchpad>/per_line_invariance.py# §4.3: per-line imposed bias table
TILE=data/derived/elbaext_delong      ./lidar-icp/bin/python <scratchpad>/fp_minus_upland.py
TILE=data/derived/elbaext_independent ./lidar-icp/bin/python <scratchpad>/fp_minus_upland.py
./lidar-icp/bin/python <scratchpad>/star_tie.py           # §2 Rank 1: star tie vs chain
./lidar-icp/bin/python <scratchpad>/dod_by_line.py        # §6.1: product-grid reproduction attempt
./lidar-icp/bin/python <scratchpad>/dod_by_line_variants.py # §6.1: seven cell-set variants
```

Nothing was written to `data/`, `src/` or any product. No CSF was run, no tile rebuilt, no gen2
cloud streamed. The heaviest step was reading `data/csf_cache/elbaext.las` (561 MB,
15,581,591 points) three times, ~50 s each.

**Documents read, not redone:** `analysis/SWATH_ALIGNMENT_METHOD.md`,
`analysis/ACROSS_TRACK_TIE_LITERATURE.md`, `analysis/LOCAL_TIE_CHAINING.md` §§1-6,
`README.md` §"What step 5 replaced", `src/lidar_diff_icp/coreg.py`,
`src/lidar_diff_icp/localtie.py`, `src/lidar_diff_icp/pipeline.py` (`register_gen1`,
`apply_datum`), `src/lidar_diff_icp/registration.py`,
`analysis/modules/vegetation_correction/gen1_save_angles_slope.py` (for `d_mm`'s sign:
`d = (z_gen1 − gen2_plane)/nnorm × 1000`, line 120 — positive means gen1 ABOVE gen2).
