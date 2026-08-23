# Steady-state planar hillslope cells as a datum / DoD check

## Rationale

For linear hillslope diffusion `dz/dt = K * grad^2(z)`, a cell with **zero
planform Laplacian curvature** (`grad^2 z = 0`, locally planar) has zero flux
divergence: it delivers as much sediment downslope as it receives from upslope.
In a landscape near **topographic steady state** such a cell's surface elevation
should not change over time, *provided it is not mass-wasting*. Restricting to
**slope < 15 deg** conservatively excludes mass-wasting cells.

Planar, low-slope cells are therefore an **independent, geomorphically-motivated
set of expected-zero-change cells**. The gen2 - gen1 elevation difference (the
DoD) over them should center on ~0 mm if (a) the landscape is near steady state
and (b) the vertical datum tie between epochs is correct. Crucially, this check
uses **no hard (built) surfaces** — it is a purely geomorphic cross-check on the
datum, orthogonal to the roof/road tie-point machinery.

**Steady-state caveat (assumption, not fact):** topographic steady state is
*assumed*. A landscape that is net-aggrading or net-incising would move planar
cells coherently, and this test cannot separate that real signal from a datum
offset. Over the ~13-yr interval on soil-mantled hillslopes the expected
steady-state departure is tiny (diffusive-signal budget below), so the
assumption is defensible here — but it remains an assumption.

## Thresholds and why

| Threshold | Value | Justification |
|---|---|---|
| `max_slope` | 15 deg | Conservative mass-wasting exclusion. (Selected cells have median slope ~10 deg, max 15 deg — genuine hillslopes, not flats.) |
| `eps_curv` (primary) | 0.00302 1/m | Symmetric `|grad^2 z|` band keeping the **central 30%** of the curvature distribution over the base population — isolates the *genuinely planar* cells. |
| cover | core forest | Soil-mantled diffusive hillslopes; avoids built/agricultural surfaces. |
| finite | curv, slope, DoD all finite | usable data. |

**Curvature units** are 1/m (`grad^2 z = d2z/dx2 + d2z/dy2`; verified
`curv_laplacian == curv_xx + curv_yy`). Typical `|kappa|` over forest is
~0.005 1/m.

**Why the central-quantile band rather than the physical diffusion budget:** the
diffusion budget is *far more permissive* than "genuinely planar," so it does not
constrain the selection. With `dz/dt = K*kappa*dt`, `dt = 13 yr`:

| K (m2/yr) | eps keeping signal < 5 mm | diffusion signal at our band edge (0.00302 1/m) |
|---|---|---|
| 0.002 | 0.192 1/m | 0.08 mm |
| 0.010 | 0.039 1/m | 0.39 mm |
| 0.050 (stiff) | 0.0077 1/m | 1.96 mm |

Even for a stiff K, the residual diffusion signal at our chosen planar-band edge
is **< 2 mm** — far below the ~50 mm forest DoD NMAD. The steady-state read is
therefore **not limited by residual diffusion** but by datum + scan-geometry
noise. The tighter geomorphic constraint ("truly planar") is what the
central-quantile band enforces; the physics only confirms that even a loose band
is safe.

**Cell counts vs eps choice** (core forest, slope<15, finite DoD): central 20% ->
713 cells; **central 30% -> ~1070 cells** (used); central 40% -> 1426 cells. The
base population is 3566 cells.

## Results — gen2 - gen1 elevation difference PDF

`eps_curv = 0.00302 1/m`, `slope < 15 deg`, core forest:

| Set | n | median (mm) | NMAD (mm) | mean (mm) | IQR (mm) |
|---|---|---|---|---|---|
| **Steady-state, geoid datum** (primary) | 1070 | **+19.8** | **51.3** | +19.3 | 69.4 |
| **Steady-state, parabola datum** | 1069 | **+18.3** | **48.4** | +19.0 | 65.2 |
| all core forest, geoid datum | 3566 | +18.1 | 50.6 | +19.4 | 68.2 |
| all core forest, parabola datum | 3565 | +16.5 | 48.1 | +17.9 | 64.6 |

Figure: `figures/refdatum/steady_state_diff_pdf.png`

## Honest read

**The steady-state cells do NOT center on zero.** They carry a **+19.8 mm**
(geoid) / **+18.3 mm** (parabola) median offset — statistically
indistinguishable from the all-core-forest median (+18.1 / +16.5 mm). Restricting
to the flattest, most planar forest cells **did not remove the forest offset**.

The two datum choices agree to ~1.5 mm (geoid +19.8 vs parabola +18.3), so the
**datum choice is not the source** of the offset — the offset is common to both
ties.

Sizing it: the offset is ~**19-20 mm**, which is the same order as the "~20 mm
forest signal" and roughly a quarter of the "~67 mm datum" scale. The NMAD (~50
mm) is unchanged from the full forest population — planarity does not tighten the
spread either.

**This is the informative negative result:** the ~19 mm forest positive offset is
**not** a curvature / steady-state (diffusive-evolution) artifact. It survives on
the most planar, low-curvature cells and it is common to both datum ties. A
diffusion signal would (i) scale with curvature and (ii) make concave cells gain
and convex cells lose; instead, concave (k<0) and convex (k>0) selected cells
differ by only ~5 mm and in the *opposite* sense to diffusion (+18.8 concave vs
+23.5 convex). The offset therefore points at the **gen2 scan-geometry /
leaf-on-canopy incidence artifact** documented earlier — a systematic gen2-high
bias in forest — rather than at real landscape evolution or a datum-tie error.

## Confounds (flagged)

1. **Curvature is computed from the gen2 DEM.** The planar selection therefore
   shares the gen2 reference frame; a systematic gen2-side artifact is not
   removed by selecting on gen2 curvature. This is consistent with the offset
   surviving.
2. **The forest DoD carries the scan-geometry / incidence artifact** (gen2 flown
   2021-05-01 at green-up, leaf-on; documented +offset in forest). This is the
   most likely source of the residual +19 mm and is not addressed by planarity.
3. **Steady state is assumed**, not established. A uniform ~19 mm of net
   aggradation across all planar forest hillslopes over 13 yr (~1.5 mm/yr) is
   geomorphically implausible for these bluffs, which is itself evidence the
   offset is instrumental rather than geomorphic — but the assumption stands as
   a caveat.
4. **Mild concave skew:** the symmetric `|kappa|` band selects 69% concave / 31%
   convex cells (the forest population is concave-biased). Mean curvature over
   selected cells is -0.00075 1/m, i.e. ~0.5 mm of diffusion signal at stiff K —
   negligible, and it cannot explain a 19 mm offset.

## Reusable code

- `analysis/steady_state/steady_state_cells.py` — data-agnostic module:
  `steady_state_mask()`, `extract_diff()`, `extract_elevations()`,
  `diff_stats()`, `eps_curv_from_quantile()`, `eps_curv_from_diffusion()`.
- `analysis/steady_state/run_steady_state.py` — Elba driver (this analysis).
