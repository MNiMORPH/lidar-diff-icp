# Ridgeline curvature, hillslope diffusion, and the lidar canopy offset

Elba, MN (Whitewater valley). gen1 = 2008-11 lidar, gen2 = 2021-05 3DEP.
DoD = gen2 − gen1. All results below use the **clean OSM hard-surface vertical datum**
(const −82.5 mm, no ill-posed tilt, `resid_nmad` ≈ 30 mm; parabola tie **deactivated**),
full-density gen2 class-2 ground, cached CSF gen1, registered by
`align_swaths → Nuth & Kääb lateral (−65.9, −1.7 cm) → const datum → drift`.

## 1. Why ridgelines

On a **divergent ridge crest** there is no incoming material — nothing flows in from
upslope. So the measured elevation-change rate there is *local*: real lowering (erosion),
plus any *measurement* offset. That lets us separate a geomorphic signal from a lidar
artifact on two **independent axes** (on the crests, canopy cover ⟂ slope, corr ≈ −0.08):

    dz/dt(crest) = −E(geomorphic) + B(lidar) + noise

- **E (geomorphic):** hillslope soil transport. Linear diffusion (Culling 1960; hilltop
  curvature as erosion-rate proxy, Roering/Hurst) gives ∂z/∂t = K∇²z, so on a crest
  (where along-ridge curvature ≈ 0) **erosion rate ∝ curvature**. The nonlinear slope
  term (Roering et al. 1999) grows as gradients steepen.
- **B (lidar):** a cover-dependent measurement offset — gen2 (2021 leaf-on) reads high
  under canopy; ag reads low (2008 fall stubble high vs 2021 spring bare). Not geomorphic.

## 2. Curvature

At every ridgecrest pixel we compute grid-aligned curvature by a windowed 2nd-order
(Savitzky–Golay) fit over ±15 m: **d²z/dx² and d²z/dy²** (convex-up ⇒ negative);
Laplacian = d²z/dx² + d²z/dy². Saved per pixel (`ridgecrest_pixels.npz`) and as full-map
grids (`curv_xx/yy/laplacian.npy`). Ridgecrest curvature CDF: forest crests are
**sharper** (median κ 0.0078 /m) than open (0.0052) — important, because κ and canopy
then covary on forest crests.

## 3. Farmland — two-parameter best fit (solved)

Two unknowns: **K_ag** (diffusivity) and **c_ag** (constant error term).

    dz/dt = K_ag · d²z/dx² + c_ag

Method: ordinary least squares of per-cell dz/dt on [d²z/dx², 1] over all open crest
cells (n = 3202; canopy penetration ≥ 0.45), dt = 12.44 yr, on the OSM-datum DoD.

| curvature used | **K_ag** (m²/yr) | **c_ag** (mm/yr) | R² |
|---|---:|---:|---:|
| d²z/dx² | **0.0636** | **−0.96** | 0.001 |
| Laplacian (d²z/dx²+d²z/dy²) | 0.122 | −0.28 | — |

- **K_ag ≈ 0.06–0.12 m²/yr** — far above natural soil creep (10⁻³–10⁻²), squarely in the
  **tillage-diffusion** range: cultivated crests. Robust to the datum fix.
- **c_ag ≈ −0.96 mm/yr** — a curvature-*independent* lowering: the **seasonal residue
  offset** (2008 fall stubble sits high → gen2−gen1 negative), ag mirror of the forest
  canopy offset (opposite sign), possibly plus uniform tillage/wind lowering.
- **Caveat on the fit:** R² ≈ 0 — per-cell dz/dt is noise-dominated (~±4 mm/yr over
  12.4 yr), so K_ag and c_ag are *population-mean* estimates with wide per-cell scatter;
  c_ag (the mean) is well determined, K_ag (the curvature slope) is weak and the two
  curvature choices disagree by ~2×.

## 4. Forest — two-parameter best fit (offset solved; K not identifiable)

Two unknowns: **K_for** (diffusivity) and **β** (lidar-offset strength, with f = β·cover
so no canopy → no offset).

    dz/dt = K_for · d²z/dx² + β · cover

Method: OLS of per-cell dz/dt on [d²z/dx², cover] over forest crest cells
(n = 5691; penetration < 0.25).

| fit | **K_for** (m²/yr) | **β** (mm/yr per unit cover) | R² |
|---|---:|---:|---:|
| joint (both free) | **−0.030** (unphysical) | **+1.11** | 0.003 |
| K fixed to K_ag = 0.064 | 0.064 (imposed) | +1.82 | — |

Nonparametric check — shared K, one offset per cover bin (f rises monotonically with
canopy density, confirming the offset is real and cover-driven):

| canopy cover | f (mm/yr) |
|---|---:|
| 0.57 | +0.3 |
| 0.79 | +0.4 |
| 0.87 | +0.9 |
| 0.94 | +1.7 |

- **β (the offset) is clean and positive**, scaling with canopy density — the lidar
  canopy artifact, as found independently (density, not slope, drives it).
- **K_for is NOT identifiable.** The free fit returns a **negative (unphysical)** K, and
  it flips sign by elevation band. Two real causes: (1) on forest crests **κ and cover
  are correlated** (forest crests are both sharper-κ and denser-canopy — §2), so K·κ and
  β·cover trade off; (2) **R² ≈ 0** — the ~1 mm/yr curvature term is buried under the
  ~±4 mm/yr per-cell noise.

So the **density axis (β) is well-constrained; the curvature axis (K_for) is not** with
the data as-is.

## 5. Open items
- **K_for extraction:** either (a) fix K to a physical prior / K_ag and read f as the
  forest residual (assumes K_for ≈ K_ag, which we expect to differ), or (b) 2-D bin by
  (κ, cover) to break the confound, or (c) reduce noise (coarsen / more stable ground).
- **K may vary with elevation** — upper dolostone caprock vs lower slopes; the elevation
  bands hint at variation but are not yet reliable (same noise/confound problem).
- **Datum tilt:** the const is tight (3 clustered OSM hard surfaces near the town); it
  can't constrain a tilt. **Fix identified:** pull **MnDOT "Roadway Surface Material in
  Minnesota"** (paved roads, GeoJSON via `gis.data.mn.gov`) — roads thread the whole
  valley, so paved segments (buffered inward, resurfaced stretches excluded, road crown
  handled) give spatially-distributed flat-hard references that constrain a real tilt.
  Backups: FEMA / Microsoft building slabs and MnGeo's 2008-era lidar footprints as
  spatially-spread pad locators.
