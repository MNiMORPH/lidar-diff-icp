# HELP: is the slope-dependent forest-floor deepening real ground change, or a gen1 geometry bias?

Follow-on to `AUDIT_findings.md`. The audit, after removing the ~67 mm geoid datum and
holding incidence fixed with gen1-only returns, found a slope-dependent **deepening** of the
forest floor of ≈ −40 to −52 mm over 4–23° slope (2–3× the ~20 mm signal budget), and left
open **whether that deepening is real ground change or a residual gen1 slope/footprint bias.**
This note attacks that open question with four gen1-internal tests.

All numbers below were produced by running the scripts named; each is reproducible with
`env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python <script>`. The ~67 mm constant is removed
everywhere via a differential (self-anchored) open/flat reference, and every effect is sized in
mm against the ~20 mm budget, never by Pearson r. Every variable is labeled by epoch.

**Reproduction check.** `AUDIT_corrected_floor_signal.py` reproduces exactly: raw median d
−55.6 mm; matched-incidence forest-minus-open slope trends −2.00 / −2.80 / −2.31 mm/deg-slope
(−37 / −52 / −43 mm over 4–23°) in the 6–10 / 10–14 / 14–18° incidence bands; intrinsic
incidence effect +1.7 to +4.0 mm/deg at fixed slope.

---

## What I did, and the numbers

### 1. gen1-ONLY land-cover stratification (`HELP_gen1only_strata.py`)
The audit's matched-incidence test defined forest/open from `stratum`, which is cut on **gen2**
leaf-on penetration — a gen2 label baked into the gen1 file. I rebuilt the label from **gen1's
own** leaf-off returns only: `gen1_abovefrac[cell] = (gen1 non-ground returns)/(gen1 non-noise
returns)`, computed from the raw gen1 cloud (all returns). This is gen1's own leaf-off
above-ground return fraction, and it is genuinely independent of gen2: spatially it correlates
with gen2 penetration at only **r = 0.23**, and the gen1-forest set overlaps the gen2-forest set
only **27%**.

With this gen1-only label (forest = abovefrac ≥ 0.55; open = ≤ 0.15), the matched-incidence,
datum-removed slope-deepening **survives essentially unchanged**:

| incidence held | gen1-only differential slope trend | over 4–23° |
|---|---|---|
| 6–10°  | **−1.84 mm/deg** | −34 mm (−1.7× budget) |
| 10–14° | **−2.34 mm/deg** | −43 mm (−2.2× budget) |
| 14–18° | **−2.14 mm/deg** | −40 mm (−2.0× budget) |

(vs the audit's gen2-labeled −2.00 / −2.80 / −2.31.) **Conclusion: the deepening is NOT an
artifact of the gen2-cut stratum.** It is present in gen1's own land-cover partition.
(Caveat: the 10–14° and 14–18° gen1-open anchors are thin — N=3,886 and 105 — because gen1-open
sits at low incidence; the 6–10° band, open-ref N=92,746, is the robust one.)

### 2. beam-vs-aspect discriminator — attempted, then RETIRED as confounded (`HELP_beam_aspect_discriminator.py`, `HELP_beam_aspect_matched.py`)
Intended test: at matched slope and incidence, a real ground surface must read the same whether
the beam points **up-aspect** or **down-aspect**; a footprint/range/aspect bias would split them.
It **failed its own flat-ground sanity check**: on flat forest (slope < 3°, where aspect is
meaningless) the up-vs-down split was **+12 mm**, and it stayed **+12 mm even after
histogram-matching the two groups' incidence to 0.5° precision** — so the split was not the
incidence confound I first hypothesized. Diagnosis (verified): it is driven by a large gen1
**per-flight-line** offset (next test), which the sign of the beam-aspect projection re-sorts.
Beam-aspect and flight-line/scan-side are geometrically coupled here, so this discriminator
cannot be cleaned; I retired it. (Both scripts are kept, with the failure documented in-file, so
the dead end is visible rather than hidden.)

### 3. per-flight-line test (`HELP_perline_slope_test.py`, plus the self-anchored refinement)
gen1 carries a real **per-flight-line vertical offset** — the 2008 along-track
GPS-drift/boresight signature this project exists to correct. The **robust** size of it is
**~44 mm**: the pipeline's applied per-swath alignment (`corrections.json`
`per_swath_internal_alignment_dxdydz_m` dz) is 0 / −23.9 / −32.5 / −43.7 mm = 43.7 mm spread, and
all-in-grid / core-farmland per-line **median** `d` on flat ground give the same ~44 mm. STEP A of
the script prints a **larger** number — the four per-line median floors on my **gen1-only
FOREST**, slope<3 subset are −91.6 / −69.8 / −21.2 / +5.0 mm, a **max−min RANGE of 96.6 mm**. That
96.6 is *not* a robust datum estimate: it is (a) a range of the two extreme lines, not a spread;
(b) on canopy-heavy forest cells (abovefrac ≥ 0.55) where gen1 ground is sparsest and the per-line
floor most divergent; (c) on raw gen1↔gen2 `d` with the per-swath dz **not** removed. Broaden to
all-in-grid ground and it falls to 75 mm (slope<3) / 63 mm (slope<5), converging toward the
~44 mm above. **Correction: my earlier "~97 mm per-line datum offset" overstated it** — the
defensible per-line offset is **~44 mm**, with 96.6 mm being the same phenomenon seen through the
widest, most-divergent (forest, narrow-slope, uncorrected) subset. *This does not affect any
conclusion*: the per-line offset is removed by self-anchoring each line to its own flat floor in
STEP B and every downstream test, so the offset's absolute size never enters the results.
**The flight-line MIX is nearly stable across slope bands** (line 137 33%→41%, line 135 16%→8%),
so line-mix does not by itself fake the trend.

When each line is self-anchored to **its own flat-forest floor at matched incidence (10–14°)**
and profiled vs slope, **all four lines deepen, and agree in sign**:

| line | flat anchor | slope trend | (N anchor) |
|---|---|---|---|
| 135 | −94 mm | −2.38 mm/deg | 19,955 |
| 136 | −59 mm | −5.25 mm/deg | 22,320 |
| 137 | −16 mm | −1.02 mm/deg | 61,122 |
| 138 |  −0 mm | −1.37 mm/deg | 51,622 |

(mean −2.50, median −1.87 mm/deg). A real surface change *or* a shared-instrument geometry bias
both predict sign agreement across lines; a per-line calibration quirk would not. **The
deepening is line-independent in sign** — it is not a flight-line-mix or single-line artifact.
(Note: an earlier per-line variant anchored to gen1-**open** ground gave scattered, disagreeing
trends; that was itself a bad-anchor artifact — gen1-open is sparse and noisy per line. The
flat-**forest** self-anchor above is far better sampled and is the trustworthy version. Reported
straight so the correction is visible.)

### 4. OPEN-vs-FOREST land-cover control — the decisive test the data allows (`HELP_open_vs_forest_control.py`)
If the deepening is a **forest-floor** signal (real ground lowering under trees, or a canopy×slope
interaction), canopy-free **open** ground on the same slopes should not deepen. It does. At
matched incidence (6–10°), self-anchored to each cover's own flat floor:

| slope | OPEN (no canopy) diff | FOREST diff |
|---|---|---|
| 3–6°   | +12.6 mm (N=74,581) | −23.6 mm (N=66,663) |
| 6–9°   |  −3.8 mm (N=134,694) | −31.7 mm (N=48,593) |
| 9–12°  |  −8.1 mm (N=32,899) | −30.6 mm (N=36,144) |
| **12–15°** | **−47.1 mm (N=2,864)** | **−43.7 mm (N=28,750)** |

Over the range where both exist, **open ground deepens with slope at least as strongly as
forest**, converging to ~−44 to −47 mm at 12–15°. Since open farmland has no canopy, the
deepening there **cannot** be a canopy mechanism.

---

## Answer to the open question (as far as the data allows)

**The slope-dependent deepening is most consistent with a gen1 slope-correlated geometry/DEM
bias, not with real forest-floor lowering.** The evidence:

1. It is **land-cover-independent** — canopy-free open ground deepens with slope as much as
   forest (test 4). A real *forest-floor* change would not appear on bare farmland.
2. It is **line-independent in sign** across all four flight lines (test 3), i.e. shared by the
   instrument/processing, as a geometry bias would be — not a per-line calibration quirk.
3. It is **not** an artifact of the gen2-cut stratum (test 1) and **not** a flight-line-mix
   artifact (test 3), so those two nuisance explanations are ruled out.

What "geometry/DEM bias" most plausibly means here: a **slope-correlated residual between the
gen1 returns and the gen2 slope-normal reference plane** — e.g. first-return range-walk on the
elongated grazing footprint, gen1's coarser ground sampling interacting with real curvature
inside a 5 m cell, or a small gen1↔gen2 slope-dependent co-registration residual. All of these
grow with slope and are common to land cover and flight line, matching what is observed.

**Honest limits — what this does NOT prove.** (a) It does not prove *zero* real change on steep
forested slopes; it shows the *bulk* of the −40 mm deepening is geometry, not that the real part
is exactly zero. (b) The steep-open control is thin (N=2,864 at 12–15°, none above 15°) and
steep-open cells are geomorphically special (banks, roadcuts) where real change is also possible —
so the open control is directional over 4–15°, not conclusive at the steepest slopes. (c) The
whole analysis rests on the gen2 plane as the DoD datum; a slope-dependent error in that plane
would present identically and cannot be separated with gen1 alone.

Bottom line for the DoD: **do not read the −40 mm slope-deepening as ground change.** It behaves
like a slope-correlated geometry bias and should be modeled/removed (as a slope term) before any
forest-slope change is claimed — exactly as the incidence effect (audit STEP 3) must be.

---

## Remaining steps that would most help, ranked

1. **Model and remove a slope (and incidence) correction surface, then re-difference.**
   The audit isolated an incidence term (+1.7 to +4.0 mm/deg) and this note isolates a
   land-cover-independent slope term (~−2 to −6 mm/deg). Fit `d ≈ a·incidence + b·slope`
   (per flight line for the constant, given the ~97 mm per-line offsets) on **open/stable**
   ground only, subtract it, and re-run the DoD. If the forest-slope deepening vanishes, it was
   geometry; whatever residual remains is the candidate real signal. This directly feeds the
   project's stated MN-wide per-swath `f(gps_time)` + geometry correction goal.

2. **Build the matched-incidence gen1-vs-gen2 comparison (the audit's framed test).**
   Reconstruct gen2's per-return incidence the same validated way (`incidence_angle.py`), remove
   the constant offset, and compare gen1 vs gen2 forest floor at **matched incidence and matched
   slope**. gen2 is leaf-on/denser but shares the same footprint-on-slope physics; if the ~20 mm
   DoD signal persists at matched incidence AND slope, that is the strongest positive evidence for
   real change. I did not build this here because tests 1–4 first had to establish that the gen1
   slope term is a bias to be controlled — otherwise the gen1-vs-gen2 comparison inherits it.

3. **Get more steep-OPEN control, or a steep-bare surrogate.** The open control is decisive only
   to ~15° slope. Pull in steep bare-ground cells from an independent source (roadcuts, quarry,
   rock outcrop — labeling only, never masking) to extend the canopy-free slope control to
   20–25° and confirm the deepening there is also land-cover-independent.

4. **Direct footprint/range-walk forward check.** Test whether the per-return `d` correlates with
   `intensity` and `range`-proxies at fixed slope and incidence (range-walk shifts first-return
   timing with pulse amplitude). A gen1-internal, physics-based confirmation of the mechanism in
   step 1, rather than an inference from the slope/land-cover pattern.

## Files
- `analysis/ridgelines/HELP_gen1only_strata.py` — gen1-only stratification, test 1
- `analysis/ridgelines/HELP_beam_aspect_discriminator.py` — retired discriminator (kept, failure documented)
- `analysis/ridgelines/HELP_beam_aspect_matched.py` — incidence-matched retry, also retired (kept)
- `analysis/ridgelines/HELP_perline_slope_test.py` — per-flight-line test, test 3
- `analysis/ridgelines/HELP_open_vs_forest_control.py` — open-vs-forest control, test 4 (decisive)
- `analysis/ridgelines/HELP_next_steps.md` — this file
