# Steady-state DoD check on two additional strata: all-forest & bluffland farmland

Zero-curvature ("steady-state") DoD check, reusing `steady_state_cells.py`
(`steady_state_mask`, `eps_curv_from_quantile`, `diff_stats`) exactly as
`run_steady_state.py` does. Driver: `run_steady_state_strata.py`.

**Datum: GEOID tie only** — `data/derived/elba_refdatum/dod_geoid.npy` (gen2 − gen1).
The parabola tie is a flexible warp that absorbs real hillslope change and is
deliberately not used.

Selection (identical to prior runs): locally planar `|∇²z| < eps_curv` AND
`slope < 15°` AND cover mask AND finite DoD. `eps_curv` = central-30% `|κ|` band
computed on each stratum's own low-slope population.

## Results (mm)

| stratum | cover_n | ss_n | median | NMAD | mean | eps_curv (1/m) |
|---|---:|---:|---:|---:|---:|---:|
| **ALL forest** (pen<0.25 & ¬floodplain) | 81474 | 10100 | **−10.3** | 76.8 | −21.5 | 0.00164 |
| **Bluffland farmland** (see def.) | 8951 | 1986 | **+31.1** | 37.4 | +26.9 | 0.00123 |
| core forest [context] | 16826 | 1070 | +19.8 | 51.3 | +19.3 | 0.00302 |
| upland farmland [context] | 42765 | 12611 | −2.2 | 37.0 | −1.9 | 0.00087 |

Context curves reproduce the known references: core forest **+19.8 mm** (exact),
upland farmland **−2.2 mm** (the prior −4.9 mm used a 66th-pct elevation split of
`core_open`; here `core_open` restricted to the plateau, z ≥ 330 m).

## Bluffland definition

Open/farmland = `penetration ≥ 0.45 & ¬floodplain_mask & isfinite`. The
open-farmland elevation distribution is strongly **bimodal**: a dense flat
**plateau** at ~330–348 m and a low **valley-floor / near-floodplain terrace** at
~215–225 m, with the **dissected bluff** terrain spanning the sparse middle
(~230–325 m) where the steep bluff faces are.

**Bluffland farmland = open cells at mid elevation, 230 ≤ z < 330 m** — below the
plateau base (excludes the flat upland) and above the valley-floor terrace
(excludes the high-scatter alluvial / floodplain-margin band, whose signal is
hydrologic, not hillslope: valley-floor z<230 gave median −40.6 mm, NMAD 98).

This is not an arbitrary elevation slice: **46.9%** of bluffland cells lie within
30 m of a steep (>20°) bluff face, versus **3.0%** of plateau farmland cells — the
band is genuinely dissected bluff terrain.

## Honest read vs the known offsets

**Neither prior expectation held.**

1. **ALL forest does NOT reproduce the core-forest +19.8 mm offset.** It comes in
   at **−10.3 mm** — slightly *negative*, and with a much broader spread (NMAD
   76.8 vs core-forest 51.3; mean −21.5 mm, a strong negative tail). The full
   forest population (adding edges and marginal/thin-canopy cells to the eroded
   interior) does not carry the core-forest's positive offset; it washes it out
   and tips slightly negative. So the +19.8 mm is **specific to the core-forest
   interior**, not a property of forest cover at large. This is the informative
   negative result: the offset does not generalize across the forest population.

2. **Bluffland farmland does NOT stay near zero like the upland plateau.** It sits
   at **+31.1 mm** (vs upland −2.2 mm, and even above core forest's +19.8 mm),
   with a tight NMAD (37.4). A clean, well-sampled positive offset in *open*
   farmland — with no canopy — means the offset seen on the bluffs is **not
   canopy-specific**. It tracks the **dissected-bluff terrain / slope setting**,
   not the vegetation. The upland-vs-bluffland farmland contrast (−2.2 vs +31.1
   mm) is a ~33 mm terrain-driven split within a single cover class.

**Net:** the data argue *against* "the offset is canopy-specific." A canopy-free
bluffland stratum shows a *larger* positive offset than core forest, while the
full forest population shows none. The common factor is dissected-bluff terrain,
not vegetation — consistent with a slope-geometry / scan-incidence origin rather
than a canopy-penetration origin.

## Confounds (flagged, not corrected)

- **Shared gen2 frame.** DoD, curvature, and slope are all derived from the gen2
  DEM. The "planar" cells are planar *in gen2*, not independently in gen1, so the
  selection and the difference share a frame.
- **Scan-incidence / point-density artifact.** Any per-cell difference in scan
  geometry or ground point density between the 2008 and 2021 epochs is not
  removed, and can bias a terrain-dependent (bluff-vs-plateau) contrast. The
  bluffland offset's alignment with steep terrain is exactly where such an
  incidence artifact would show up, so it cannot be cleanly separated here.
- **Steady-state assumption.** As in the core-forest run, topographic steady state
  is assumed; net incision/aggradation would move planar cells coherently. On a
  ~13 yr interval the diffusive budget is small, but bluff terrain is the setting
  most likely to violate it.

## Files

- Figure: `figures/refdatum/steady_state_allforest_bluffland_pdf.png` (836×505 px)
- Driver: `analysis/steady_state/run_steady_state_strata.py`
- Reused module: `analysis/steady_state/steady_state_cells.py`
