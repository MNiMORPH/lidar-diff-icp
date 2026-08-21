# Slope-dependent DoD bias — investigation record

**Status (2026-08-21): open. Best-motivated mechanism (slope-induced pulse
broadening) is not proven; a competing vegetation-penetration mechanism is now in
play. Bedrock-outcrop control is the decisive discriminator (pending).**

This directory preserves the *process*, not just the conclusion — including the
laser-beam / footprint-geometry work, which is reusable methodology regardless of
whether broadening turns out to be the cause.

---

## 1. The observation

Elba DoD (gen2 2021 − gen1 2008), bare-earth, 5 m grid, has a systematic
**aspect-independent bias that scales with tan(slope)**, gen1 reading low on slopes:

| slope | median DoD |
|---|---|
| 2–10° | +3.6 mm |
| 10–20° | +8.1 |
| 20–30° | +13.7 |
| 30–90° | +28.8 |

Robust harmonic decomposition of `dod/tan(slope)` vs aspect (slope 5–40°):
- **dipole (= residual horizontal shift): 0.5 cm** → NOT horizontal registration.
- quadrupole (E–W symmetric): 1.4 cm (secondary, unexplained).
- **aspect-independent term: 3.4 cm-equiv** (dominant); `DoD ≈ 35 mm × tan(S)`, R²≈0.77.
- Flat ground is pinned to ~0 by the flat-hard datum references → NOT a vertical shift.

## 2. Mechanisms tested and RULED OUT (the process)

| candidate | test | result |
|---|---|---|
| ground estimator | median vs plane vs poly2 (windowed quadratic) | all ~+29 mm at steep — **invariant** |
| curvature-aggregation | poly2 (exact on synthetic paraboloid) | no change → cancels between epochs |
| point density | decimate gen1→gen2 density, re-grid | steep bias mostly survives (+23 of +29) |
| extraction method | CSF on both epochs; last-return | same (+25.9) |
| scan-angle / incidence range-walk | residual vs scan angle at fixed slope | flat, wrong sign — **but this test controlled OUT the slope term** (see §4, methodological lesson) |
| boresight residual | opposing-swath overlap within gen1 | flat on flat ground; slope disagreement inconsistent-sign (sampling, not roll) |
| horizontal registration | Nuth–Kääb residual dipole | 0.5 cm, negligible |

## 3. Laser-beam / footprint physics (reusable work — the point of preserving this)

### Incidence angle per return (`compute_incidence.py`, `save_incidence_linked.py`)
- `cos i = cos θ cos S + sin θ sin S cos Δψ`, θ = scan angle, S = slope, Δψ = scan
  azimuth − upslope azimuth. Footprint ellipse semi-axes (ρ, ρ/cos i), ρ = R·γ/2
  (Baltsavias 1999; Sheng 2008; Zhang & Shen 2014).
- Computed for every gen1 return; **sanity verified: median(i − S) = 0.00° at nadir**.
- Saved LINKED TO THE PULSES: `data/las_local/gen1_geom.las` (incidence_deg / slope_deg
  / dpsi_deg as extra per-point dims) and `data/derived/elba/gen1_pulse_geometry.npz`
  (keyed by gps_time + point_source_id + return_number).

### Forward broadening model (`broadening_forward.py`, `broadening_waveform.py`)
- Return elevation spread `σ_geom = (D/2)/√2 · tan(S)` — aspect-independent, ∝ tan(S).
  **The measured bias has this exact form** (R²≈0.77).
- `DoD_bias = (k/2)(D₂₀₀₈ − D₂₀₂₁)·tan(S)`, k = detector offset below footprint centroid.
- Measured c = 35 mm ⇒ **required k ≈ 0.40 σ** — plausible for a below-centroid /
  walk-prone detector, **zero for a peak/centroid detector.** Spread ≠ bias.

### Footprint sizes (verified from survey reports + datasheets)
| | H AGL | divergence | **D (1/e)** |
|---|---|---|---|
| 2021 TerrainMapper (dominant) | 2000 m | 0.177 mrad | **0.35 m** |
| 2008 Gemini — narrow (standard) | 2400 m | 0.25 mrad | **0.60 m** |
| 2008 ALS50-II | 2400 m | 0.156 mrad | 0.37 m |
- **2008 footprint larger under all realistic assumptions** (higher flight + larger
  min divergence) → sign is right (gen1 low). The Gemini beam SETTING (narrow/wide)
  is NOT documented for the 2008 survey — pending (§6).

### Detector / range-walk physics (agent-verified)
- Leading-edge timing `t = b − c·ln(a/V_th)`: a broadened/slower pulse (larger c, from
  slope) crosses the fixed threshold LATER → range long → **ground low.** Right sign.
- CFD / peak / centroid are ~walk-free (caveat: CFD has residual rise-time walk;
  centroid biases if broadening is asymmetric).
- **Discriminator type is UNDISCLOSED for all four sensors.** 2008 = hardware-timed
  discrete-return (waveform optional); 2021 TerrainMapper = real-time full-waveform.
  Gemini likely CFD+AGC (least walk-prone). **Per-sensor susceptibility: undetermined.**

## 4. Methodological lesson (kept visible)

The first incidence "test" regressed the residual on incidence **controlling for
slope**. Because incidence ≈ slope (92% collinear) and the mechanism lives in the
slope term, that held fixed the very thing it should have varied — turning a
CONFIRMATION into an apparent refutation. The fix was the FORWARD model (derive the
form from geometry, test that form, reduce the remainder to an independent spec
check). See global CLAUDE.md "Price the wasted-rework time" and
memory `default-principled-not-fast`.

## 5. The vegetation cover-split (2026-08-21) — reframes toward vegetation

Splitting DoD-vs-slope by gen2 canopy fraction:

| slope | BARE (canopy<15%) | VEGETATED (canopy>50%) |
|---|---|---|
| 2–10° | **−2.4** (n=13535) | **+12.0** (n=67060) |
| 10–20° | +0.8 (n=306) | +8.2 |
| 20–30° | — (n=20) | +13.9 |
| 30–90° | — (n=14) | +29.1 |

- **Bare ground is ~unbiased even where sloped-ish; the bias tracks vegetation.** At
  matched slope (10–20°): bare +0.8 vs veg +8.2.
- **Even flat forest shows +12 mm** (gen2 higher) — a vegetation offset independent of
  slope. Points at a ground-under-vegetation penetration difference between epochs,
  possibly phenology (gen2 late-spring vs gen1 November-dormant — see §6).
- **This CHALLENGES broadening** (which is cover-independent — should bias bare rock too).
- **Limitation:** steep terrain here is almost entirely forested (only 19 bare cells
  >25°), so the data-derived proxy cannot test bare-STEEP.

## 6. Decisive next test + pending needs

- **BEDROCK OUTCROP as stable slope control** (agent finding the dataset). Bedrock is
  bare AND steep AND stable — the one surface that discriminates: bedrock DoD ≈ 0 on
  slopes ⇒ vegetation mechanism (broadening out); bedrock shows +tan(S) ⇒ sensor/broadening.
- **Raw sensor info (Andy flagged):** 2008 AeroMetric acquisition logs (beam setting)
  + discriminator specs / raw waveforms, from MnGeo / MN DNR / NV5-Quantum-Spatial.
- **2021 leaf-out / phenology:** vendor report states "leaf-off conditions", but the
  block flight window (May 31 2021 – May 17 2022) is borderline leaf-on in Minnesota;
  the flat-forest +12 mm suggests gen2 carried more low vegetation than gen1. Need the
  TILE-level acquisition date (gen2 gps_time is zeroed in the LAS) and a leaf-out check.
- **Correction fork (unresolved):** correct c·tan(S) fit on bedrock/stable control vs
  a slope/vegetation-aware LoD. Hold until bedrock resolves vegetation-vs-sensor.

## Citations
Baltsavias 1999 (ISPRS J 54:199); Sheng 2008 (IEEE GRSL 5:419); Zhang & Shen 2014
(IEEE GRSL 11:701) — footprint/incidence geometry. Gardner 1992 (IEEE TGRS 30:1061);
Kirchhof/Jutzi/Stilla 2008; Wagner 2006 — slope pulse broadening. Range walk: Sensors
17:2369 (2017), PMC5948876, PMC5677387 (t=b−c·ln(a/V_th)). Su & Bork 2006; Aguilar
2010; Hodgson & Bresnahan 2004 — slope×density DEM error. Kraus & Pfeifer 1998;
Viedma 2022; Bickel 2002/2006 + Bartels 2006 — robust/modal/correct-to-dense ground.
Nuth & Kääb 2011 — co-registration.

## Code & data
- Scripts (this dir): `save_incidence_linked.py`, `compute_incidence.py`,
  `broadening_forward.py`, `broadening_waveform.py`, `poly2_prototype.py`.
- Pipeline: additive `ground="plane"` and `ground="poly2"` options (uncommitted).
- Data: `data/derived/elba/{dod_refdatum,dod_poly2,z_after}.npy`,
  `gen1_pulse_geometry.npz`, `data/las_local/gen1_geom.las`.
- Figures: `figures/{broadening_forward,incidence_test,boresight_test}.png`.
