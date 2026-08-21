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

## Consequence for the workflow
1. The forest +dz is an artifact; do **not** count it as deposition. Apply the `f(veg_frac)`
   forest-floor correction (`dod_corrections.py`), or — the durable fix — pair **leaf-off
   with leaf-off** epochs so the exposed surface matches. This is the MN-wide lesson: gen1
   is dormant-season; gen2 3DEP tiles flown at green-up will carry this bias (cf. the
   mn-3dep-audit NDVI screen).
2. Open farmland's +50 mm is mostly the common datum plus a smaller seasonal
   stubble-vs-tillage surface difference; it is tied out by the hard-surface/geoid datum.
