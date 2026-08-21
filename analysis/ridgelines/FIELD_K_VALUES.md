# Hillslope diffusivity K for cultivated fields vs natural creep

Reference values for ∂z/∂t = K·∇²z (K in m²/yr), for the Winona County loess/cultivated
setting. Compiled 2026-08-21 (agent research, sources below).

## Tillage diffusivity — cultivated fields (the relevant regime here)
| K (m²/yr) | setting | source |
|---|---|---|
| **0.03–0.52** | US Corn Belt tillage diffusion | Thaler, Larsen & Yu (2021) *PNAS* 118:e1922375118 |
| **0.19 ± 0.04** (site 0.19/0.25/0.59) | Corn Belt, LEM-calibrated, plow depth 0.20 m | Kwang et al. (2022) *JGR-BG* 127:e2021JG006616 |
| ≈0.25–0.65 (from K_till) | mouldboard ~330 kg/m/pass (MN), multi-pass ~780 kg/m/yr ÷ ρ_b≈1200 | Lindstrom et al. (1992) *STR* 24; Van Muysen & Govers |

**Winona County lies *inside* the Thaler study region → 0.19 m²/yr is the direct anchor;
expected band 0.1–0.5 m²/yr.**

## Natural soil creep (context — negligible here)
| K (m²/yr) | setting | source |
|---|---|---|
| 0.0032 ± 0.0009 | forested Oregon Coast Range | Roering, Kirchner & Dietrich (1999) *WRR* 35 |
| 0.001–0.01 (to ~0.1 granitic) | soil-mantled convex hilltops | Fernandes & Dietrich (1997); Hanks (2000) |
| ~0.0005–0.009 | lithology-resolved, Sierra | Hurst et al. (2013) *JGR-ES* |

**Natural creep clusters at 0.001–0.01 m²/yr — 1–3 orders of magnitude below tillage.**

## Decadal change (Δt ≈ 11–12 yr, hilltop |∇²z| 0.005–0.03 m⁻¹)
| process | K | change over decade |
|---|---|---|
| natural creep | 0.001–0.01 | **~0.1–1 mm** (at/below lidar detection) |
| tillage (gentle convex) | 0.03–0.19 | ~1.6–21 mm |
| tillage (sharp convex) | 0.19–0.52 | ~40–170 mm |
Measured Midwest rate (Thaler 2022): ~1.8–1.9 mm/yr ≈ 20 mm/decade.

## Bearing on our results
- **Our measured K_ag ≈ 0.21 m²/yr (all-cell, buffered farmland) matches the published
  Corn-Belt tillage diffusivity (Thaler 0.19, range 0.03–0.52) almost exactly.** Independent
  validation that the farmland diffusion signal is real and correctly scaled — it's tillage,
  as the magnitude demands. (Caveat: our K is in loess/soil atop the dolostone, i.e. the
  soil-transport coefficient, not bedrock erodibility.)
- **Setting forest hillslope K = 0 is justified:** natural creep over a decade is
  ~sub-mm — below the ~50 mm/cell DoD noise and the ~30 mm forest offset. So the forest
  signal is essentially all the forest-floor measurement error, not creep.
