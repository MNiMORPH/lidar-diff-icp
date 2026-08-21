# Vertical return structure — the source of the false forest elevation rise

**Question (Andy):** From the per-cell slope-normal return histograms, can we identify
ground / undergrowth / trees, and see whether gen1 or gen2 has an issue that produces the
false forest elevation rise (gen2 reads higher than gen1)?

**Data:** `data/derived/elba_fulldensity/slope_normal_returns.npz` — per-cell histograms of
slope-normal height `d` (0.25 m bins, −1…40 m) for every return, split into *all* returns
and *ground-classified* returns, for both epochs. `d` is measured perpendicular to a
**single common reference plane: the gen2 bare-earth surface** (so gen1 and gen2 `d` are
directly comparable). gen2 carries ~24× more returns than gen1; all comparisons are per-epoch
**percentiles/shape**, never raw counts.

Scripts: `return_structure.py` (layer profiles + median driver), `return_penetration_test.py`
(classifier-independent deep-echo test). Figures in `figures/refdatum/`.

## Layers are cleanly separable
`return_structure_all.png`: both epochs show the textbook three-layer column — a sharp
**ground** spike at `d≈0`, an **understory** minimum near `d≈1.5–2 m`, and a **canopy** body
above. gen2 (leaf-on, May green-up) carries far more mid-canopy mass than gen1 (leaf-off,
Nov dormant) — expected, and the reason gen2's *median return* sits up in the canopy while
gen1's sits at the ground.

## The rise is a leaf-off vs leaf-on ground-EXPOSURE difference, driven by gen1's deeper penetration
All heights slope-normal, common frame:

| quantity | forest | open |
|---|---|---|
| ground-class median rise, per-cell (gen2−gen1) — *what the DoD differences* | **+80 mm** | **+50 mm** |
| → forest-specific differential (survives the open-tied datum) | **≈ +30 mm** | — |
| deepest echoes, physical floor p0.5 (gen1 / gen2) | −0.81 / −0.54 m | −0.55 / −0.37 m |
| gen1 penetrates deeper than gen2 at p0.5 | **+266 mm** | +179 mm |
| gen1−gen2 floor gap at p2 / p5 / p10 | +240 / +195 / +132 mm | +8 / +15 / +26 mm |

Reading `return_penetration.png`:
- **In forest, gen1's deep-echo CDF sits well below gen2's** across the entire low tail
  (0.5–12 %). gen1 (leaf-off, Nov 2008) reaches a *lower* ground surface than gen2
  (leaf-on, May 2021). **In open, the two CDFs coincide** — the penetration gap is a
  forest phenomenon, not a global sensor difference.
- The count asymmetry works *against* this finding: gen2 has ~24× more shots, so its low
  percentile samples far deeper into the tail by count alone — yet it *still* bottoms out
  ~2–27 cm **higher** than gen1 in forest. That makes the obscuration conclusion conservative.
- gen2 is not fully blocked: ~10 % of gen2 forest returns still reach below gen1's ground
  median (open: ~32 %). But gen2's ground *class* median sits +81 mm above gen1's — in
  leaf-on conditions the bulk of low returns land on the spring flush / fresh litter / low
  growth, so the median-of-ground surface is pulled up even where deeper points exist.

## Mechanism (what is actually happening)
This is **not a gen2 error and not a gen1 error** — it is a real difference in the surface
each survey *sees*, set by phenology:
- **Nov 2008 (gen1): leaf-off, dormant.** Pulses pass through bare deciduous canopy to the
  true forest floor → lower ground surface.
- **May 2021 (gen2): green-up, leaf-on** (verified: NDVI 0.49, gps_time flight date).
  Emerging canopy + spring herbaceous groundcover + fresh litter intercept pulses higher;
  even at 24× density the deepest returns stop above gen1's, and the ground-class median
  sits ~30 mm higher than open (beyond the common datum).

DoD = gen2 − gen1 therefore reads **positive = false aggradation** on forested hillslopes.
It is the same signal as the `f(veg_frac)` forest-floor offset (R²=0.72) already characterised,
now traced to its physical cause in the return column.

## What this rules in / out
- **Ruled out — gen2 understory leaking into the ground class:** gen2's upper tail
  (`frac(d>0.5 m)`, `frac(d>1 m)`) is *smaller* than gen1's, and the rise does not track
  `understory_frac` (r=+0.05). The whole ground distribution shifts up coherently; it is not
  a fat-tail contamination.
- **Ruled out — a global gen1-vs-gen2 sensor/noise offset:** in open ground the deep tails
  coincide (p2–p10 within 8–26 mm) and there is no forest-sized penetration gap. gen1's
  slightly greater raw noise (σ≈0.017 vs 0.011 m) cannot produce the forest-specific
  266 mm floor gap.
- **Ruled in — leaf-off/leaf-on ground exposure:** the penetration advantage is forest-only,
  in the direction and magnitude expected from Nov-dormant vs May-green-up phenology.

## Decomposition — groundcover shoulder vs leaf-off penetration (`groundcover_decomp.py`)
Andy's reading of `return_structure_ground.png`: gen2's broader near-ground peak = gen1's
true ground (still detected) *plus* an added upward shoulder from leaf-on herbaceous
groundcover. Tested directly on the forest ground-class column (unit area, no scaling):

- **Common ground MODE preserved: gen1 −0.125 m, gen2 −0.125 m (shift 0 mm).** The
  true-ground return sits in the identical place both epochs — gen2 does not shift off
  the ground. Confirms the model's core claim.
- **gen2 adds a +19.2-point upper shoulder** (mass above the mode: 30.1 % → 49.3 %),
  concentrated ~0–0.15 m — the leaf-on groundcover/litter layer. This is the dominant
  half of the asymmetry.
- **gen1 has +8.3 points more deep tail** (mass below the mode: 10.9 % vs 2.6 %) — the
  leaf-off penetration effect. Smaller, opposite side, *same* sign of DoD error.

So the forest rise is TWO mechanisms of the same leaf-state mismatch, both making
gen2−gen1 positive: (1) gen2 groundcover shoulder up (~⅔), (2) gen1 deeper penetration
down (~⅓).

Two consequences for correction design:
- **The median rise does not scale with groundcover amount** (corr = +0.007 vs the correct
  0.1–0.6 m near-ground band; *not* `understory_frac`, which is the 0.5–2 m shrub band and
  missed this). The median saturates — any shoulder past the 50 % point bumps it ~one step
  then stops. So a linear per-cell `f(groundcover)` off the median is not recoverable.
- **A lower ground percentile recovers the true ground under the cover — partially:**

  | estimator | forest rise | open rise | forest-specific anomaly |
  |---|---|---|---|
  | p50 (current) | +80 mm | +50 mm | **+30 mm** |
  | p25 | +54 mm | +39 mm | **+14.5 mm** |
  | p10 | +100 mm | +80 mm | +20 mm |

  p25 halves the forest anomaly by dropping beneath the groundcover shoulder — no external
  veg model. p10 overshoots: it falls into gen1's deeper-penetration tail and re-inflates
  the difference. A single percentile cannot cancel both mechanisms; ~p25 is the sweet spot.

## CORRECTION at 1 cm resolution — penetration dominates, not a groundcover shoulder
The decomposition above (0.25 m bins, `groundcover_decomp.py`) over-credited a "groundcover
shoulder." Re-streamed at **1 cm** (`ground_mixture_fit.py` → `ground_fine_pooled.npz`,
`ground_mixture_fit2.py`) the ground CLASS tells a different, cleaner story. Three
frame-invariant diagnostics (independent of the reference plane):
- **Upper tails coincide:** gen1 & gen2 forest ground-class share p90 (+0.123 m); the
  per-percentile rise is **+129 mm at p10 but 0 mm at p90** — it *falls* low→high.
- **Monotone-falling rise** (p10 +129, p25 +104, p50 +70, p75 +34, p90 0) is the signature
  of gen1's **deeper low tail (leaf-off penetration)**, not a gen2 upper shoulder (which
  would make the rise grow toward high percentiles).
- **gen2 forest is symmetric** (Bowley skew −0.004); if groundcover leaked upward into the
  ground class, gen2 forest would be the most right-skewed — instead gen1 *open* is.

This agrees with the earlier classifier-independent deep-echo test (gen1 forest floor −0.81
vs gen2 −0.54 m). The intermediate coarse "shoulder-dominant" claim was a 0.25 m binning +
cell-subset artifact and is **retracted**.

**Two-Gaussian fit on gen2 class-2 ground (1 cm):** the ground class is essentially ONE
symmetric Gaussian — μ_g = −0.001 m, σ_g = 0.080 m (forest) / 0.023 m (open); the 2nd
"plant" component collapses onto the same mean (symmetric halo, offset +0.000 m), i.e. no
resolvable upward plant mode. Extraction level = μ_g ≈ the median; nothing to strip from the
ground class. The groundcover returns are real (visible in the all-returns column) but 3DEP's
classifier files them as VEGETATION, so they don't enter the ground class — the "plant leak"
into ground is negligible.

**Corrected mechanism:** gen1 (leaf-off) reaches a lower, truer forest floor; gen2 (leaf-on)
ground sits ~70 mm higher (pulses can't penetrate as deep). gen1 is closer to true ground;
gen2 is biased high. Both ground classes are internally clean. The mixture decomposition is
the right tool but belongs on the **all-returns** near-ground column (where plant returns
actually live and gen2 still has deep ground returns under the canopy), applied identically
to both epochs to recover ground independent of vendor classifiers — the promising next step.

## Consequence for the workflow
1. The forest +dz is an artifact; do **not** count it as deposition. Apply the `f(veg_frac)`
   forest-floor correction (`dod_corrections.py`), or — the durable fix — pair **leaf-off
   with leaf-off** epochs so the exposed surface matches. This is the MN-wide lesson: gen1
   is dormant-season; gen2 3DEP tiles flown at green-up will carry this bias (cf. the
   mn-3dep-audit NDVI screen).
2. Open farmland's +50 mm is mostly the common datum plus a smaller seasonal
   stubble-vs-tillage surface difference; it is tied out by the hard-surface/geoid datum.
