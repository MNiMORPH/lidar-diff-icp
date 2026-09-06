# Audit: methodological contamination in the gen1 forest-floor analysis

Scope: bare-earth DoD (gen2 − gen1) at Elba, MN. gen1 = 2008 lidar (leaf-off,
sparse ~0.87 pts/m²); gen2 = 2021 3DEP (leaf-on, ~24× denser). All numbers below
were re-derived by running code against the data, not taken from prior summaries.

Reproduce the corrected numbers with:
`env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/AUDIT_corrected_floor_signal.py`

---

## Verified baseline facts

- gen1 all-return density = **0.869 pts/m²** (7,728,747 returns / 8.89e6 m²).
  At this density gen1 **cannot resolve a 0.5–2 m understory as vertical
  structure** — returns are ~1 m apart. Any gen1 "understory"/"forward-scatter"
  mechanism is unsupported by the data.
- `d_mm` in `gen1_csf_angles.npz` is **intrinsically gen1-vs-gen2**: it is the
  gen1 return elevation minus the **gen2** reference plane `z_after.npy`,
  slope-normal (`gen1_save_angles_slope.py:45`). That is correct for a DoD; the
  problem is what is done *with* it.
- Incidence is beam-to-**gen2**-surface-normal; on flat ground it reduces to
  `|scan angle|`, so the reconstruction is geometrically valid regardless of epoch.
- The raw median `d` ≈ **−56 mm** (all gen1 ground) is dominated by a **constant
  ~+67 mm gen1↔gen2 geoid/co-registration datum** (GEOID03→GEOID18), which the
  prior analyst correctly identified in `ground_fine_datumshift.py` (`DATUM =
  0.067`). The **signal is only the ~20 mm of structure on top of that constant.**
- **Slope and achievable incidence are nearly collinear (corr ≈ 0.85 in forest).**
  The scanner is only ±17° off-nadir, so on steep slopes low-incidence
  (slope-perpendicular) beams are geometrically impossible. Verified:

  | slope band | min achievable incidence | median incidence |
  |---|---|---|
  | 0–3°   | 0.0° | 7.1° |
  | 15–20° | 0.2° | 18.5° |
  | 25–30° | **9.3°** | 27.5° |
  | 50–70° | **40.1°** | 55.4° |

---

## Problem 1 — hidden gen2→gen1 data mixing

Gen2-derived variables (`penetration.npy`, `canopy_struct.npz` fields
`understory_frac`/`canopy_height_p95`/`canopy_cover`/`veg_frac`, and the
`stratum`/`core_forest`/`core_open` labels, all cut on gen2 leaf-on penetration)
are used to explain **gen1-internal** behavior in the files below. `z_after.npy`
as the slope-normal datum is the accepted common reference and is *not* counted
as contamination; the contamination is using a gen2 **canopy magnitude** or a
**gen2-cut stratum** as the explanatory axis or filter for a gen1 question.

### Tier 1 — gen2 canopy/penetration is the explanatory axis for a gen1 result
| # | file:line | gen2 variable | gen1-internal conclusion it drives |
|---|---|---|---|
| 1 | `gen1_sink_vs_density.py:19` | `penetration.npy` | "gen1 ground median sinks as [gen2] canopy thickens" |
| 2 | `gen1_ground_mechanism.py:21` (+`core_forest.npy:22`) | `penetration.npy` | gen1 single-return floor "still sinks with [gen2] density" |
| 3 | `gen1_elev_vs_3.py:15,34,46` | `penetration.npy` | `corr(gen1 d, gen2 pen)` framed as gen1 covariate |
| 4 | `test_incidence_veg_hypothesis.py:18–19,45` (+`stratum` L17) | `canopy_height_p95`, `understory_frac` | gen1 floor incidence/d regressed on gen2 veg |

### Tier 2 — gen1 conclusion confined to / contrasted on the gen2-cut "core" stratum
`core_forest`/`core_open` are cut on gen2 penetration + gen2 `veg_frac`; an
all→core "sink" or a core-only form is therefore a gen2-density statement.

| # | file:line | gen2 variable |
|---|---|---|
| 5 | `gen1_nadir_elev_intensity.py:18` | `core_forest.npy` |
| 6 | `gen1_core_nadir_vs_oblique.py:18` | `core_forest.npy` |
| 7 | `gen1_csf_nadir_all_vs_core.py:18,20` | `penetration.npy`,`core_forest.npy` |
| 8 | `gen1_csf_nadir_oblique_hist.py:20` | `core_forest.npy` |
| 9 | `gen1_csf_vs_internal.py:20–21` | `penetration.npy`,`core_forest.npy` |
| 10 | `incidence_angle.py:27,29` | `penetration.npy`,`core_forest.npy` |
| 11 | `ground_fine_core.py:21–22` | `core_forest.npy`,`core_open.npy` |
| 12 | `core_vs_orig_forms.py:12` | (consumes gen2-cut core from #11) |

### Tier 3 — gen2 selection baked into a gen1 artifact used downstream
| # | file:line | gen2 variable |
|---|---|---|
| 13 | `gen1_save_angles_slope.py:21,23–24,49–54` | writes gen2-cut `stratum`,`core_forest`,`core_open` into `gen1_csf_angles.npz`; every consumer (e.g. #4) inherits gen2 selection |

### Tier 4 — gen1 form/distribution conclusions on gen2-cut forest/open strata
(No direct gen2 load in-file; the forest/open labels in the pooled/csf npz were
cut on gen2 penetration in the producer.)
`gen1_ground_dist.py:17`, `gen1_form_logderiv.py:19`, `gen1_tail_powerlaw.py:23`,
`ground_fine_semilog.py:11`, `ground_mixture_fit2.py:18`,
`ground_csf_all_analyze.py:14` — all via `ground_fine_pooled.npz` /
`ground_fine_csf_all.npz`.

### Selection-only (milder; measured axes are gen1, but "forest" label is gen2)
`gen1_pen_vs_intensity.py:15`, `gen1_penfrac_by_angle_stratum.py:14–18`,
`gen1_csf_oblique_strata_hist.py:18,22`, `gen1_intensity_fit.py:17`,
`gen1_combined_nadir_vs_oblique.py:17`, `ground_mixture_fit.py:38`,
`ground_fine_csf_all.py:27`, `all_returns_full.py:18`.

### Clean (gen1-only, or gen2 loaded only to *expose* the contamination)
- **`gen1_own_penetration.py`** — the corrective file: builds gen1's OWN leaf-off
  penetration and bins gen1 d against it, loading gen2 `penetration.npy` only to
  plot it as the misleading axis. Draws no gen1 conclusion from gen2.
- **`diagnose_06_split.py`** — gen1 penetration/scan-angle/structure from the gen1
  cloud only.

**Recommended gen1-only alternative to the strata:** define forest/open from a
gen1-only descriptor (e.g. gen1's own leaf-off penetration as in
`gen1_own_penetration.py`, or land cover from an independent source used only for
*labeling*, never for masking). Never bin a gen1 quantity against a gen2 canopy
magnitude.

---

## Problem 2 — constant offset vs signal, and r-vs-systematic-shift

- The raw `d` ≈ −56 mm is ~67 mm **constant datum** + ~20 mm **signal**. Any
  analysis that reports absolute `d` (~−56 to −70 mm) and correlates it against a
  covariate is conflating the 67 mm constant with the 20 mm signal.
- **Per-return scatter is ±272 mm** (std of `d`), *not* ±150 mm as previously
  stated — even wider, so correlation coefficients are tinier still. Verified in
  forest: `corr(d, slope) = −0.09`, `corr(d, incidence) = −0.003`. **These small
  r-values do not mean "no effect."** The *systematic median shift* across the
  covariate range is tens of mm (the whole budget and more). **Judging by r was a
  repeated error** — e.g. `gen1_save_angles_slope.py:69` reports `corr(d,slope)`;
  `test_incidence_veg_hypothesis.py:45–46` and `gen1_elev_vs_3.py:34,46` judge by
  `corrcoef`. Reframe every effect as **systematic median shift in mm relative to
  the 20 mm budget**, after removing the constant.

---

## Problem 3 — mixed-angle / correlation contamination (the central finding)

Per-cell median floor pools returns across incidence, and incidence is ~collinear
with slope (0.85). The floor reading depends **strongly and intrinsically** on
incidence — verified gen1-only at fixed slope (`AUDIT_corrected_floor_signal.py`
STEP 3):

| slope band | d/d(incidence) | floor span over achievable incidence |
|---|---|---|
| 6–9°   | +1.7 mm/deg | +29 mm (1.5× budget) |
| 12–15° | +2.1 mm/deg | +62 mm (3.1× budget) |
| 18–21° | +2.7 mm/deg | +80 mm (4.0× budget) |
| 25–30° | +4.0 mm/deg | +123 mm (6.2× budget) |

Grazing (high-incidence) beams read the floor **higher** (less negative d);
perpendicular beams read **lower**. This effect alone is several × the signal
budget, and its magnitude rises with slope. **Therefore any mixed-angle
median-vs-slope or median-vs-canopy comparison is contaminated by the shifting
angle composition.**

**The contamination reverses the apparent slope effect.** Mixed-angle forest floor
vs slope looks flat/slightly rising (−65→−63 mm, +2 mm over 4–23°), which reads as
"steep ground reads higher." But holding incidence fixed and removing the datum
(forest floor − open floor at the same incidence band):

| incidence held | differential slope trend | over 4–23° |
|---|---|---|
| 6–10°  | −2.0 mm/deg-slope | −37 mm (1.8× budget) |
| 10–14° | −2.8 mm/deg-slope | −52 mm (2.6× budget) |
| 14–18° | −2.3 mm/deg-slope | −43 mm (2.1× budget) |

The sign **flips to strongly negative** and is consistent across all three
incidence bands: the forest floor genuinely **deepens** with slope once angle is
controlled. The mixed-angle "flat" trend was an artifact of the composition
sweeping to higher incidence (which raises d) as slope steepens.

Note: linear residualization of d on incidence does **not** fix this
(`corr` stays −0.09→−0.087) because the incidence effect is nonlinear and
slope-dependent. **You must bin at matched incidence**, not partial out a line.

---

## Problem 4 — unsupported mechanisms / logic

- **gen1 understory / forward-scatter:** unsupported. At 0.869 pts/m² gen1 cannot
  measure a 0.5–2 m understory layer as structure. `understory_frac` is a **gen2**
  quantity (`slope_normal_returns.py:200–205`); using it to explain gen1's floor
  (`test_incidence_veg_hypothesis.py`) attributes a gen1 effect to a layer gen1
  never resolved. A leaf-**off** November canopy also does not provide the dense
  scattering medium such mechanisms assume.
- **r-driven "no effect" and "effect" calls:** repeatedly downweighted the 20 mm
  budget (Problem 2). Effects were declared present/absent by |r| ~ 0.03–0.2 when
  the systematic shift was the whole budget.
- **gen2-cut "core sharpens gen1" claims** (Tier 2): the sharpening is a statement
  about gen2 leaf-on density, not an intrinsic gen1 property.

---

## Correct analysis frame (going forward)

1. **Remove the ~67 mm constant** before any covariate work. Best practice here:
   anchor to the **open reference at matched incidence** (differential d), so both
   the datum and the common incidence effect cancel and the residual is genuinely
   the differential signal.
2. **The signal is the ~20 mm residual.** Report effects as **systematic median
   shift in mm relative to 20 mm**, never as Pearson r on ±272 mm per-return data.
3. **Hold incidence FIXED** (bin at matched incidence); do not use mixed-angle
   per-cell medians. Slope↔incidence collinearity (0.85) makes mixed-angle
   slope/canopy comparisons invalid.
4. **Use each epoch's own data, labeled.** For gen1 questions, use gen1-only
   penetration/labels; never a gen2 canopy magnitude or a gen2-cut stratum as the
   explanatory axis.

---

## Verdict

After removing the ~67 mm constant datum and controlling incidence with gen1
returns, **there is a real residual forest-floor signal, and it is large relative
to the 20 mm budget.**

- The dominant per-return nuisance is the **incidence effect** (+1.7 to +4.0
  mm/deg, up to +123 mm across the achievable range) — this is geometry, not
  landscape change, and must be controlled.
- Once incidence is fixed and the open datum removed, the forest floor **deepens
  with slope** at **≈ −2.0 to −2.8 mm/deg-slope**, a systematic differential of
  **≈ −40 to −52 mm over the 4–23° slope range = about 2–3× the 20 mm budget** —
  and this is the **opposite sign** from the mixed-angle appearance.
- The covariate the residual depends on is **slope** (with incidence held fixed);
  the previously-claimed dependence on **gen2 canopy/understory is not
  gen1-measurable** and is an artifact of gen2 stratification plus the
  slope↔incidence↔canopy collinearity.

Blunt summary: the prior "flat / floor-reads-higher-on-steep-slopes" reading was
an artifact of (a) leaving the 67 mm datum in, (b) reading r on ±272 mm scatter,
and (c) mixing incidence angles that are collinear with slope. Corrected, the real
gen1 forest-floor signal is a slope-dependent deepening of order 2–3× the signal
budget, on top of an even larger pure-geometry incidence effect that dwarfs both.
Whether that slope-dependent deepening is real ground change vs. a residual gen1
slope/footprint bias is the next question — but it is **not** explained by gen2
canopy, and it is real and much larger than r-values suggested.
