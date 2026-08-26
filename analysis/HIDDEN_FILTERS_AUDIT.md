# Audit: undisclosed count filters and covariate truncations in the analysis scripts

**Date:** 2026-08-26
**Scope:** every `.py` under `analysis/`, `scripts/` and `src/lidar_diff_icp/`.
**Trigger:** `SESSION_AUDIT_AND_ENFORCEMENT.md` §1.1 — a self-invented `--minn 20` in
`q2_cover_fit.py` deleted 35 of the 69 highest-cover cells and turned a filter artefact into
a reported power law. That one is fixed (`--minn` default is now 1). This audit finds the rest.

**Method.** Every number below was produced by re-running the script (a copy with its figure
output redirected to a scratch directory, so nothing in `figures/` was touched) or by
replicating its selection exactly against the same inputs. Nothing here is from recall.
**No script's default behaviour was changed.** Nothing was committed.

**Scan result:** 149 count-filter or truncation sites across 63 files. Most are one-line
`if b.sum() < N: continue` guards inside superseded, single-use diagnostic scripts. The ones
that matter — the ones behind a figure or a stated number Andy has used — are the 8 named in
the session audit plus **three the session audit missed**: a covariate truncation in
`incidence_correction_fit.py`, a slope truncation in `nearnadir_slope_dependence.py`, and a
hardcoded per-cell return filter still live in `q2_cover_fit.py` itself.

---

## Summary table — filters behind live results

| # | file | line | expression | exposed? | what it actually removes (measured) | damage |
|---|---|---|---|---|---|---|
| 1 | `analysis/ridgelines/nearnadir_slope_dependence.py` | 117, 128 | `FIT_MAX_SLOPE = 35.0`; `good = … & (ns >= 200) & (cen < FIT_MAX_SLOPE)` | **buried** | 2 of 13 slope bins (35–40°, n=14,647, r=−35.6 mm; >40°, n=3,623, r=+57.9 mm) from the tan-law/step fit | **HIGH — changes a reported statistic from significant to not** |
| 2 | `analysis/ridgelines/incidence_correction_fit.py` | 37 | `--inc-max` default `35.0` ("drop the sparse high-incidence tail") | CLI, but framed as a nuisance cut | 18 populated incidence bins, 36–70°, **52,797 returns (1.43%)** — the bins holding the *largest* deltas (+28 → +1138 mm) | **HIGH — changes the fitted correction and its shape** |
| 3 | `analysis/ridgelines/offset_by_beam_selection.py` | 53 | `MIN_BIN = 100` | **buried** | forest panel of `_sp5` figure: **16 of 22 non-empty slope bins, 45% of forest points, across the whole 0–44° range** | **HIGH for that figure** |
| 4 | `analysis/ridgelines/offset_model_slope_cover.py` | 28 | `--min-cells` default `30` | CLI | 5 of 72 grid boxes (n = 7–28 cells each), **and prints a false reason for their absence** | **MED-HIGH — a printed claim is untrue** |
| 5 | `analysis/ridgelines/offset_vs_angle.py` | 38 | `MIN_N = 300` | **buried** | 5–6 of ~105 bins per run, always the extreme x-tail (incidence/slope ≥ 36–39°, `|scan|` ≥ 16°); 0.01–1.5% of returns per panel | **MED-LOW — see §4** |
| 6 | `analysis/ridgelines/q2_cover_fit.py` | 52 | `ok = stable & (n1[cells] >= 5) & (ng >= 10) & …` | **buried, hardcoded** | 1,302 of 92,404 stable cells (1.41%) — but **cover-enriched**: median cover 0.346 vs 0.185 kept; top cover bin loses 33 of 102 cells (32%) | **MED — principle; the committed slope survives** |
| 7 | `analysis/ridgelines/offset_model_slope_cover.py` | 27 | `--min-n` default `3` (returns per cell) | CLI | 752 of 96,653 cells (0.78%); dropped cells median cover **0.356 vs 0.193**, median slope 30.5° vs 7.3°, median offset −114.5 vs −60.2 mm | **MED — small volume, 14.5× cover enrichment** |
| 8 | `analysis/ridgelines/offset_by_beam_selection.py` | 34 | `--min-n` default `3` (returns per cell) | CLI | 4,595 of 6,771,479 returns (0.07%) in 3,009 cells; dropped median cover **0.338 vs 0.197**, 16.6% vs 0.69% at cover > 0.50 (24× enrichment) | **LOW magnitude — but methodologically required here** |
| 9 | `analysis/ridgelines/cover_offset_isotonic.py` | 106 | `bs.binned_stats(…, min_n=5)` | **buried** | exactly **one** display bin: cover 0.85–1.00, n=3 (elbaext, −443 mm) / n=2 (elba, −245 mm) — the largest offset on the plot | **LOW-MED — display only, fit unaffected** |
| 10 | `analysis/ridgelines/nearnadir_slope_dependence.py` | 450 | `OPEN_RELIABLE = 2000` | **buried** | 1 bin (>40°, open n=291/1,028) — **but the prose it governs is false regardless** (see §1b) | **MED — false narrative, not the filter** |
| 11 | `analysis/ridgelines/cover_offset_reference.py` | 48 | `--min-n` default `200` (returns per cover bin) | CLI | **nothing** — sparsest bin is n=273 (elba divides+inc5) / n=400 (elbaext) | harmless today, **latent** |
| 12 | `analysis/ridgelines/cover_offset_regression.py` | 126 | `bs.binned_stats(…, min_n=200)` | **buried** | **nothing** — quantile bins are equipopulated by construction (13 bins with and without); check-display only, the OLS/LAD fits are unbinned | **harmless** |
| 13 | `analysis/ridgelines/nearnadir_slope_dependence.py` | 57 | `MIN_N = 30` | **buried** | **nothing** on either tile — sparsest bin over all four curves is n=219 | **dead filter** |
| 14 | `analysis/ridgelines/offset_vs_angle.py` | 136, 156 | `if bm.sum() < 500` (cover band); `if b.sum() < 50` (scatter check) | **buried** | **nothing** — no cover band skipped on either tile; scatter check is a table row only | harmless |
| 15 | `src/lidar_diff_icp/offset_model.py` | 141 | `median_surface(…, min_cells=30)` | library default | same boxes as #4; **its docstring asserts the removed boxes are "combinations the terrain does not supply"**, which is false | **doc error — fix the docstring** |
| 16 | `src/lidar_diff_icp/offset_model.py` | 172 | `matched_band_effects(…, min_n=200)` | library default, not exposed by the caller | **nothing** at current defaults — sparsest band n=572 | harmless, latent |
| 17 | `src/lidar_diff_icp/coreg.py` | 381, 419 | `min_pts=2000` per swath; `if ok.sum() < 10` | library default | **nothing** — 0 of 4 (elba) and 0 of 6 (elbaext) in-grid swaths got a zero drift curve | harmless |
| 18 | `src/lidar_diff_icp/pipeline.py` | 289 | `if m.sum() < 500: return None` | buried | guards the heteroscedastic-LoD model fit; returns `None` and falls back rather than deleting data | **harmless guard** |
| 19 | `analysis/ridgelines/percentile_float_fit.py` | 25–26 | `--min-gen1` 3, `--min-gen2` 5 | CLI | not separately quantified — but note these are **different values (3/5) from `q2_cover_fit.py`'s (5/10) on the same cells** | **comparability risk** |

### The rest (superseded / single-use diagnostics)

The remaining ~120 sites are `if <mask>.sum() < N: continue` guards, N ∈ {3, 20, 30, 50, 100,
200, 300, 400, 500}, inside scripts that produced a one-off table in a line of work the
`FRAME_2026-08-26.md` has since closed (the whole `gen1_*`, `HELP_*`, `*_pdf_by_slope`,
`glennie_*`, `test_incidence_veg_hypothesis` family — the slope/incidence hypothesis, now
attributed to per-swath misalignment). They are all buried and all invented, but none is
behind a live figure or a number in the FRAME. They are listed by file:line in the appendix.
**They should not be individually repaired; they should be left as historical record** — the
correct action is to stop the pattern in new work, which is what the enforcement hook in
`SESSION_AUDIT_AND_ENFORCEMENT.md` §M1 is for.

---

## §1 — `nearnadir_slope_dependence.py`: the truncation that made the 27° knee significant

### 1a. `FIT_MAX_SLOPE = 35` (line 117), applied at line 128

The fit uses `good = np.isfinite(meds) & (ns >= 200) & (cen < FIT_MAX_SLOPE)` — a *second*
count filter (`ns >= 200`) plus a covariate truncation, neither exposed.

Refit of exactly the script's own model on exactly its own binned medians, changing only
`FIT_MAX_SLOPE` (`gen1_csf_angles.npz`, near-nadir `|scan| < 5°`, datum-removed residual):

| tile | FIT_MAX_SLOPE | bins in fit | tan-law R² | **step@27° R²** | step amplitude | **F** |
|---|---|---|---|---|---|---|
| elba_fulldensity | **35 (default)** | 11 | 0.468 | **0.800** | −22.2 mm | **10.85** |
| elba_fulldensity | 90 (lifted) | 13 | −0.065 | **0.267** | −22.1 mm | **2.49** |
| elbaext | 35 (default) | 11 | −0.936 | 0.417 | −19.6 mm | 4.56 |
| elbaext | 90 (lifted) | 13 | −0.902 | 0.416 | −19.8 mm | 4.72 |

The two excluded bins are 35–40° (**n = 14,647**, median r = **−35.6 mm**) and >40°
(**n = 3,623**, median r = **+57.9 mm**). They are not empty and they are not negligible —
18,270 near-nadir returns.

**Verdict: the qualitative conclusion does not survive.** On elba the headline
"27° switch-on is real (R² = 0.80, F = 10.85)" becomes R² = 0.267, F = 2.49 — below the 5%
critical value F(1,10) ≈ 4.96. The step *amplitude* is stable (−22.2 → −22.1 mm); what
evaporates is the *evidence* for it. That statistic is recorded as established in
`NEARNADIR_SLOPE_DEPENDENCE.md` and in the `nearnadir-slope-27deg-knee.md` memory
("the ~27° 'switch-on' is REAL … R²=0.80 F=10.85 … NOT a tan-curve artifact"). It is an
artefact of where the fit was truncated. (Note `FRAME_2026-08-26.md` §2 already retired the
knee for an independent reason — per-swath misalignment — so this does not change the
current science; it does mean the *statistical* claim was never sound.)

### 1b. `OPEN_RELIABLE = 2000` (line 450) — the filter is minor, the prose is false

The filter drops one bin (>40°, open n = 291 elba / 1,028 elbaext). But the paragraph it
generates asserts:

> "above ~12 deg the 'open' median rests on n~30-3000 and swings wildly … the matched-slope
> canopy test can only be made where both are reliably populated"

**Measured, same run:** open (cc<0.2) near-nadir bin counts are 48,363 (12–15°), 39,335
(15–18°), 36,410 (18–21°), 37,345 (21–24°), 31,843 (24–27°), 24,199 (27–30°), 16,308
(30–35°), 2,302 (35–40°). Every bin to 40° passes the filter, and the script's own
"reliable" list runs 0–3° through 35–40°. The narrative restricting the comparison to
"~0–12 deg" is stale text left over from an earlier version of the data and contradicts the
table printed three lines below it. **This is a bigger problem than the filter.**

---

## §2 — `incidence_correction_fit.py`: `--inc-max 35` deletes the whole high-incidence regime

Run on `data/derived/elbaext`, `d_mm_corr`, divides, |Laplacian| ≤ 0.015, all covers
(3,688,609 returns). `--min-n 500` on its own **removes nothing** inside 0–35° — the sparsest
retained bin is n = 37,045. The entire bite is `--inc-max`.

Bins deleted by the 35° truncation (delta d relative to the incidence<5° anchor):

| bin (°) | n | delta d (mm) | robust SE |
|---|---|---|---|
| 36–38 | 22,317 | +27.7 | 1.85 |
| 38–40 | 13,931 | +16.1 | 2.42 |
| 40–42 | 8,126 | +37.7 | 3.51 |
| 42–44 | 4,101 | +44.5 | 5.53 |
| 44–46 | 2,135 | +65.6 | 7.46 |
| 46–48 | 1,145 | +109.3 | 11.40 |
| 48–50 | 448 | +162.7 | 27.64 |
| 50–70 | 594 (10 bins) | +15 to +1138 | 100–380 |
| **total** | **52,797 (1.43%)** | | |

The first six of these are individually significant at 5–15σ. They are not noise; they are
the rare regime, and they carry the largest signal in the dataset — exactly the case
`geoscience-data-scarcity-defaults.md` was written about.

**Effect on the fitted correction** (same AIC contest, same code, only the truncation moved):

| theta | correction with `--inc-max 35` | with the tail restored |
|---|---|---|
| 20° | −11.6 mm | −12.1 mm |
| 25° | −12.3 mm | −13.5 mm |
| **30°** | **−12.1 mm** | **−14.3 mm** |
| **35°** | **−10.9 mm** | **−14.4 mm** |

**Verdict: the qualitative shape changes.** With the truncation the quadratic turns over at
~28° and *decreases*, suggesting the geometric term saturates and reverses. With the tail
restored it is monotone rising to the edge of the data. The turnover is a boundary effect of
where the fit was cut. (Caveat, stated because it matters: the tail bins are also the
highest-cover ones — the script's own confound check reports corr(incidence, cover) = +0.93 —
so the tail is not clean geometry either. That is an argument for *showing* it with its
error bars, not for deleting it.)

---

## §3 — `offset_by_beam_selection.py`: `MIN_BIN = 100` guts the forest panel

`MIN_BIN = 100` (line 53, hardcoded, absent from the docstring and from the figure title)
filters slope bins in `curve()`, which draws three of the four panels.

Measured on `data/derived/elbaext`, `d_mm_corr`, `--min-n 3`:

| figure variant | panel | non-empty slope bins dropped | points dropped | slope span of dropped bins |
|---|---|---|---|---|
| `offset_by_beam_selection_elbaext.png` | all returns | 0 | 0 | — |
| " | open cc<0.10 | 3 | 81 (0.03%) | 38–44° |
| " | forest cc>0.50 | 3 | 156 (**2.19%**) | 38–44° |
| `offset_by_beam_selection_elbaext_sp5.png` | all returns | 0 | 0 | — |
| " | open cc<0.10 | 7 | 292 (0.39%) | 30–44° |
| " | **forest cc>0.50** | **16** | **832 (45.19% of 1,841)** | **0–44° — the whole range** |

**Verdict: the `_sp5` forest panel does not survive.** Requiring 5° of within-cell scan
spread leaves only 1,841 forest points; `MIN_BIN = 100` then deletes 45% of them, spread over
the entire slope axis, so what is plotted is a handful of surviving bins and not a curve.
The dropped bins' medians run −73 to +908 mm — i.e. the deleted points are not a tail, they
are half the panel. The default (no spread) version is fine except at 38–44°.

---

## §4 — `offset_vs_angle.py` and `MIN_N = 300` — full treatment

`MIN_N = 300` (line 38, comment: "drop bins below this (matches prior slope analysis)") is
hardcoded, is not exposed as an argument, is not in the docstring, and does not appear in the
figure title or filename. It is applied in `binned()`, which produces **every** curve in
**both** panels and every printed table.

### What it removes, per figure variant actually on disk

All runs replicate the script's own selection exactly (same columns, `in_grid`, same
`curv_laplacian` cut, same `ridge_mask`, same edges, same cover bands).

**(a) Default run — `offset_vs_incidence.png` / `offset_vs_scan_angle.png`, elba, raw `d_mm`,
all curvatures, no ridge:**

| panel | non-empty bins | dropped by MIN_N | returns dropped |
|---|---|---|---|
| all returns (incidence) | 15 | **0** | 0 |
| open cc<0.10 | 15 | **0** | 0 |
| forest cc>0.50 | 15 | **0** | 0 |
| all returns (|scan angle|) | 9 | **0** | 0 |
| forest cc>0.50 (|scan angle|) | 9 | 1 (16–18°, n=275, median **+28.8** mm vs kept range −119.6…−24.7) | 275 (0.59%) |

**(b) The headline figure — `offset_vs_incidence_curv0.015_reg_ridge_cbands.png`
(`--x incidence --curv-max 0.015 --offset corr --ridge --cover-bands`):**

elba_fulldensity, 1,928,149 returns:

| cover band | bins dropped | returns dropped | dropped-bin x span | dropped-bin medians | kept-bin median range |
|---|---|---|---|---|---|
| all returns | 0 | 0 | — | — | −6.1 … +59.2 mm |
| 0–5% | 2 | 123 (0.018%) | 39–45° | +139.2, +216.7 | −34.7 … +61.9 |
| 5–10% | 2 | 73 (0.110%) | 39–45° | +104.4, +328.3 | −7.2 … +43.3 |
| 10–20% | 1 | 259 (0.098%) | 42–45° | −4.2 | −14.5 … +22.4 |
| 20–35% | 0 | 0 | — | — | −16.7 … +16.0 |
| 35–50% | 0 | 0 | — | — | −24.4 … +21.6 |
| 50–100% | 1 | 201 (**1.499%**) | 42–45° | +131.3 | −84.0 … +150.6 |

elbaext, 3,688,609 returns: 5 bins / 772 returns dropped, same 39–45° span, same pattern
(0–5%: 2 bins/435; 5–10%: 2/165; 50–100%: 1/172 = 0.854%).

**(c) `offset_vs_slope_curv0.015_reg_ridge*.png`** (elba): 1 bin in the all-returns panel
(42–45°, n=269, median +106.7 mm); 2 bins in the forest panel (39–45°, n=122 = 0.909%,
medians **+154.8 and +880.4 mm** against a kept range of −45.3…+146.3).

### The pattern, stated plainly

`MIN_N = 300` **never** removes a mid-range bin, **never** removes a whole cover band, and
**never** touches the region where the finding lives (incidence 0–36°, cover 0–50%). It
removes only the last one or two bins of the covariate axis: incidence/slope ≥ 36–39°, or
`|scan angle|` ≥ 16°. In returns it is 0.01–1.5% per panel.

But those bins are not empty and they are not uninformative: they are the sparse extreme,
and their medians (+104 to +880 mm) sit far outside the range of everything retained. The
consequence is that **every curve in the right-hand panel stops before the data do, with no
mark saying so** — the left panel's hexbin still shows the underlying points, but the
binned-median line and the whole right panel end silently.

### Verdict on `offset_vs_angle.py`

**The figure's message survives.** The message Andy has been using — offset is organised by
canopy cover, ~0 in open ground, −50 to −85 mm at cover > 0.5 at low incidence, and flat
against incidence once cover is held fixed — is carried entirely by bins of n = 600 to
390,000, none of which the filter touches. Removing `MIN_N` changes no cover ordering, no
low-incidence level, and no cross-tile agreement. I re-ran both tiles and every variant on
disk to check this rather than assert it.

**But the filter should go, and it should not be replaced with another number.** It is
undisclosed, it is unmarked on the figure, and it does exactly the thing this project has
already been burned by once: it deletes the sparse extreme of a covariate. The honest form
is the one `binstats.binned_stats` already implements and documents — plot every non-empty
bin at its true span with its cluster-robust error bar, and let a 30-point bin show a large
error bar instead of vanishing. That also removes the current asymmetry where the left panel
shows the tail and the right panel silently drops it. I have not proposed a replacement
threshold; there should not be one.

---

## §5 — `offset_model_slope_cover.py` `--min-cells 30`: a printed claim that is untrue

Run: `--tile data/derived/elba_fulldensity --curv-max 0.015 --ridge` (matches
`offset_model_slope_cover_curv0.015_ridge.png`), raw `d_mm`, 95,901 cell medians.

`--min-cells 30` blanks 5 of 72 boxes in the (slope × cover) prediction table:

| slope band | cover band | cells | median offset (mm) — visible only with the filter lifted |
|---|---|---|---|
| 12–15° | 0.50–1.01 | 28 | −78 |
| 18–21° | 0.50–1.01 | 23 | +47 |
| 27–30° | 0.00–0.05 | 28 | −51 |
| 35–45° | 0.00–0.05 | 16 | −141 |
| 35–45° | 0.05–0.10 | 7 | −178 |

The script prints, at line 101 and in the summary line:

> `'.' = fewer than 30 cells: the terrain does not supply that combination, so it is UNSUPPORTED, not extrapolated.`
> `supported boxes: 67/72 (93%) -- the empty ones ARE the slope-cover covariance.`

**This is false.** The terrain does supply all five combinations, with 7–28 cells each. With
`--min-cells 1` the table is 72/72 and the covariance story disappears. The same false claim
is baked into the library docstring at `src/lidar_diff_icp/offset_model.py:141`
("those are combinations the terrain does not supply, and NaN is the honest answer there").

Downstream: the per-cell interaction model's fit to this table degrades from RMS 33.0 mm /
max |resid| 131.9 mm to RMS 42.6 mm / max 155.5 mm when the five boxes are restored — i.e.
the reported goodness of the interaction model is partly a product of hiding the boxes it
fits worst.

**Also in this script:** `--min-n 3` returns per cell removes 752 of 96,653 cells (0.78%),
whose median cover is **0.356 vs 0.193** for the kept cells, median slope 30.5° vs 7.3°, and
median offset **−114.5 vs −60.2 mm**. 17.2% of the dropped cells are cover > 0.50 against
1.2% of the kept — a 14.5× enrichment. Small in volume, but it is the steep-forested corner
of a table whose whole purpose is the steep-forested corner.

---

## §6 — `q2_cover_fit.py`: the `--minn` fix did not remove the cell filter

Line 52, hardcoded, no CLI, not in the docstring:

```python
ok = stable & (n1[cells] >= 5) & (ng >= 10) & np.isfinite(cover)
```

Measured (elba_fulldensity, `nearground_cells_sn.npz`, `d_mm_corr`, `reference_cells`):

- drops **1,302 of 92,404 stable cells (1.41%)**
- dropped cells: median cover **0.346**, mean 0.343 — against kept median **0.185**, mean 0.171
- top cover bin (0.65–0.93) grows **69 → 102 cells (+48%)** when lifted; 0.50–0.65 grows 688 → 833 (+21%)
- the two sparsest bins move: q2\* 0.199 ± 0.065 → **0.250 ± 0.106** (top bin), 0.388 ± 0.044 → **0.366 ± 0.037** (0.50–0.65)
- fitted linear coefficient: **b = −0.1949 → −0.1887**; χ²/dof improves for all four forms
  (linear 1.13 → 0.73, quadratic 0.77 → 0.61, power 0.79 → 0.58)

**Verdict: the committed result survives.** `q2 = 0.5 − 0.19·cover` (commit `5c88952`) and the
DoD rebuilt from it (commit `4cfd586`) are unchanged to the reported precision — both −0.195
and −0.189 round to −0.19, and the linear form is still the one Andy selected. But the filter
is undisclosed, hardcoded, and re-selects toward open canopy by exactly the mechanism
`FRAME_2026-08-26.md` §5 warns about ("A minimum-returns-per-cell filter does NOT filter
quality — it re-selects toward open canopy"). It should be exposed or removed, not retuned.

**Comparability note:** `percentile_float_fit.py` — whose code `q2_cover_fit.py` `exec`s —
applies `--min-gen1 3 --min-gen2 5` to the same cells. Two scripts in the same result chain
select two different cell populations from the same cube. Under "same-way or don't compare"
these should be one number defined in one place.

---

## §7 — `cover_offset_isotonic.py` `min_n=5`: one bin, the most extreme one

`bs.binned_stats(x, y, edges, block=blk, min_n=5)` at line 106 governs the plotted binned
medians only; the isotonic fit and its jackknife use every cell.

It removes exactly one bin on each tile — the top one:

| tile | bin | cells | median offset |
|---|---|---|---|
| elbaext | cover 0.85–1.00 | 3 | **−443 mm** |
| elba_fulldensity | cover 0.85–1.00 | 2 | **−245 mm** |

That is the single largest offset anywhere in the analysis, and it is the only bin removed.
The script's own docstring says "no bins, no bandwidth, no cover cut, no seed" and its figure
comment says the bins exist so that "each high-cover interval is drawn where it actually
sits" — both are contradicted by the line. The fit is unaffected, so no conclusion changes,
but the plot understates its own high-cover end. A 3-cell bin plotted with its (enormous)
jackknife error bar is the honest version.

---

## §8 — Filters verified to remove nothing (no action needed beyond disclosure)

Each was checked against the data, not assumed:

- `cover_offset_reference.py --min-n 200`: nothing removed in any of four runs (upland and
  divides+inc5, both tiles). **Latent risk:** the top cover bin is n = 273 (elba) / 400
  (elbaext) and carries the −140.3 / −130.0 mm point that is the whole finding. A smaller
  tile, or a min-n of 500, deletes it.
- `cover_offset_regression.py min_n=200`: nothing removed (13 bins either way, both tiles);
  quantile edges make bins equipopulated by construction, and the fits are unbinned anyway.
- `nearnadir_slope_dependence.py MIN_N = 30`: nothing removed on either tile across all four
  curves (sparsest bin n = 219). A dead filter.
- `offset_vs_angle.py` cover-band rule `bm.sum() < 500`: no band skipped on either tile.
- `offset_model.matched_band_effects(min_n=200)`: nothing removed (sparsest band n = 572).
- `coreg.py` drift `min_pts=2000`: 0 of 4 (elba) and 0 of 6 (elbaext) in-grid swaths skipped;
  every swath has a non-zero `dz_drift_mm`.
- `pipeline.py:289 if m.sum() < 500: return None`: a fallback guard on the xdem
  heteroscedastic model, not a data deletion.

---

## §9 — What needs regenerating

Ranked. Nothing here has been changed or regenerated; this is the list for Andy to pull from.

1. **`NEARNADIR_SLOPE_DEPENDENCE.md` and `NEARNADIR_SLOPE_DEPENDENCE_ELBAEXT.md`, and the
   `nearnadir-slope-27deg-knee.md` memory.** The "R² = 0.80, F = 10.85, the 27° step is real"
   statistic is conditional on `FIT_MAX_SLOPE = 35`; without it R² = 0.267, F = 2.49. Also
   the false "only reliable below ~12 deg" paragraph. *(The FRAME already retires the knee on
   other grounds — this is a correction to the record, not to current science.)*
2. **`figures/refdatum/incidence_correction_elbaext.png` and the correction table.** `--inc-max
   35` hides 52,797 returns whose deltas run to +1138 mm and creates the apparent turnover at
   ~28°. Re-run with the tail shown (and the cover confound plotted alongside).
3. **`figures/refdatum/offset_by_beam_selection_elbaext_sp5.png`.** The forest panel is drawn
   from 55% of its own points; re-run without `MIN_BIN`, or state the panel is not usable.
4. **`figures/refdatum/offset_model_slope_cover_*_ridge.png` and the printed grid table**
   (both tiles, all four curvature variants) — plus the false sentence at
   `offset_model_slope_cover.py:101` and the docstring at `offset_model.py:141`.
5. **`figures/refdatum/offset_vs_*.png`** (all 30 variants). Cosmetic-to-moderate: the curves
   are correct where they are drawn, but they stop before the data do without saying so.
   Regenerate with every non-empty bin plotted at its true span with an error bar.
6. **`figures/refdatum/cover_offset_isotonic*.png`** — restore the cover 0.85–1.00 point with
   its error bar; the fit line needs no change.
7. **No action:** `cover_offset_reference*.png`, `cover_offset_regression*.png`,
   `q2_vs_cover_fits*.png` and the committed `q2 = 0.5 − 0.19·cover` / `dod_cover_q2.npy` —
   all verified unchanged by their filters. The `q2_cover_fit.py:52` cell filter should still
   be exposed or removed on principle, but it does not move the number.

---

## Appendix — remaining filter sites (superseded or single-use)

`AUDIT_corrected_floor_signal.py` 91,105,127 · `HELP_beam_aspect_discriminator.py` 154 ·
`HELP_gen1only_strata.py` 129 · `HELP_open_vs_forest_control.py` 133 ·
`HELP_perline_slope_test.py` 134,140 · `canopy_struct_regress.py` 59 ·
`conifer_vs_deciduous.py` 50,263 · `core_vs_orig_forms.py` 20,27 · `cover_curve_select.py` 77 ·
`dod_cover_attribution.py` 71,81,89 · `forest_density_driver.py` 75,88,99 ·
`gen1_combined_nadir_vs_oblique.py` 44 · `gen1_core_nadir_vs_oblique.py` 43 ·
`gen1_elev_vs_3.py` 40 · `gen1_form_logderiv.py` 29,50 ·
`gen1_ground_mechanism.py` 53,60,79,91,101,117 · `gen1_nadir_elev_intensity.py` 44 ·
`gen1_own_penetration.py` 57 · `gen1_pen_vs_intensity.py` 42,52 ·
`gen1_save_angles_slope.py` 98 · `gen1_sink_vs_density.py` 28 · `gen1_vendor_vs_csf.py` 24,70 ·
`gen2_csf_compare.py` 25,70 · `gen2_form_logderiv.py` 21,45,69 ·
`gen2_incidence_test.py` 106,115,133,136,160,188,198 · `geology_forest_split.py` 66,84 ·
`glennie_scanangle_swath_test.py` 163,236,237 · `ground_csf_all_analyze.py` 23,37 ·
`groundcover_decomp.py` 86 · `grow_flat_surface.py` 54 · `hillslope_fits.py` 62 ·
`incidence_angle.py` 64 · `lowangle_slope_dependence.py` 35,56 ·
`nadir_vs_all_pdf_by_slope.py` 29,33,47,48 · `nearground_class_split.py` 108 ·
`nearground_profile.py` 96 · `nearground_q_for_gen2_median.py` 105 · `nearground_rank.py` 95 ·
`nearnadir_vs_perp_slope.py` 35 · `per_flightline_offset.py` 26,107,122,134 ·
`percentile_float_fit.py` 128 · `perp10_vs_all_pdf_by_slope.py` 30,34,50,51 ·
`perp_vs_all_pdf_by_slope.py` 30,34,50,51 · `return_structure.py` 109 · `roof_usability.py` 175 ·
`test_incidence_veg_hypothesis.py` 31,43,54,63 · `dep_internal_check.py` 31,38 ·
`m3c2_lesson.py` 93 · `mass_balance/elba.py` 47 ·
`roughness_characterization/oak_forest.py` 78 · `slope_bias/ground_class_structure.py` 25 ·
`slope_bias/ground_return_stats.py` 56,118 · `slope_bias/understory_from_lidar.py` 36,37,44 ·
`scripts/naip_cover_error.py` 105 · `src/detect.py` 81,84,129 · `src/coreg.py` 96,274,416 ·
`src/pipeline.py` 163,204

---

## FIXED — 2026-08-26

Every filter in the damaging list above has been removed or reduced to its definitional
floor and the affected results regenerated. Each fix is one commit on `datum-fix-hillslope`;
nothing was pushed. The rule applied throughout: **remove, do not retune**. Where a cut is
required for the arithmetic to be defined (a median needs a point, a two-parameter line fit
needs two), it is kept at that floor, exposed as an argument, and said so in the help text.
Every number below was measured by re-running, not recalled.

### Conclusions that changed

| what | before | after |
|---|---|---|
| **near-nadir 27° step, elba** | R² **0.800**, F **10.85** — reported as "the switch-on is REAL" | R² **0.267**, F **2.49**, p 0.146, below F(1,10)crit 4.96 — **UNRESOLVED** |
| near-nadir 27° step, elbaext | R² 0.417, F 4.56 | R² 0.416, F **4.72**, p 0.055 — still short of 4.96 |
| near-nadir step amplitude | −22.2 mm (elba) / −19.6 (elbaext) | −22.1 / −19.8 — **stable; the amplitude stands, the evidence for a break does not** |
| near-nadir forest-vs-open prose | "only reliably populated ~0–12°"; ">27° is ENTIRELY forest, no steep open control" | **both false** — open n = 48,363 (12–15°), 24,199 (27–30°), 16,308 (30–35°). Matched-bin median forest−open −39.3 → −38.3 mm |
| **incidence correction @30°, elbaext** | −12.1 mm, quadratic vertex at **26.2°**, so the curve *turned over* near 28° | −**14.3** mm, vertex at **33.6°**, rising to the edge of the reporting range — **the turnover was a boundary effect of the 35° truncation** |
| incidence correction @20/25/35° | −11.6 / −12.3 / −10.9 | −12.1 / −13.5 / −14.4 |
| **beam-selection `_sp5` forest panel** | 6 slope bins, all but one below zero (−68…−30 mm): reads "gen1 low in forest" | **22 bins over 0–44°**: crosses zero near 12°, +120 mm by 39°, +300…+900 mm at 41–44°. **The sign of the slope dependence was hidden, not absent** |
| **slope × cover grid, elba curv0.015** | 67/72 boxes; model-vs-table RMS **33.0** mm, max 131.9 | **72/72**; RMS **42.6** mm, max 155.5 — the blanked boxes were the ones the model fits worst |
| slope × cover grid, elba curv0.002389 | (20 of 72 blanked) RMS **22.5** mm | 72/72, RMS **54.9** mm |
| slope × cover grid, elbaext curv0.001 | (25 of 72 blanked) RMS **30.9** mm | 72/72, RMS **75.4** mm |
| median-surface slope term (elba curv0.015) | +0.037 mm/deg; ∂d/∂slope at cover 0.02 = **+0.02** | −0.017; **−0.03** — sign flip |
| matched-band cover effect, elbaext curv0.001 35–45° | "(sparse)", suppressed | **−181 mm/unit** (n=79), in line with the −143 and −250 below it |
| **isotonic top bin** | cover 0.85–1.00 not drawn | drawn: **−443 ± 386 mm** (n=3, elbaext), **−245 ± 186** (n=2, elba) — the largest offset in the analysis |
| **q2(cover) fit** | b = −0.1949; χ²/dof 1.13 (lin) / 0.77 (quad) / 0.79 (pow) | b = **−0.1887**; 0.73 / 0.61 / 0.58. **`q2 = 0.5 − 0.19·cover` and the DoD built from it SURVIVE** |
| q2 top cover bins | 0.65–0.93: 69 cells, q2* 0.199 ± 0.065; 0.50–0.65: 688, 0.388 ± 0.044 | **102** cells, 0.250 ± 0.106; **833**, 0.366 ± 0.037 |
| percentile_float_fit cells | 91,741 (gen1≥3, gen2≥5) | **92,359** — now identical to `q2_cover_fit.py`'s population on the same cube. q1* 0.5227 unchanged, q2* 0.4605 → 0.4606 |
| offset_vs_angle curves | stopped at incidence/slope ≥ 36–39°, \|scan\| ≥ 16°, unmarked | every non-empty bin drawn, each with a **cluster-robust (50 m block) error bar**. The headline — cover organises the offset, ~0 open, −84 mm at cover > 0.5 at low incidence, flat vs incidence at fixed cover — is **unchanged** |

### Filters removed, by site (numbering as in the summary table)

1, 10, 13 — `nearnadir_slope_dependence.py`: `FIT_MAX_SLOPE = 35`, the unexposed `ns >= 200`
inside the fit, `OPEN_RELIABLE = 2000`, and the dead `MIN_N = 30` are all gone; the
significance verdict now turns on the F-test against its own critical value rather than on a
raw R² increment, and the write-ups carry a dated correction-to-the-record.
2 — `incidence_correction_fit.py`: `--inc-max` deleted outright (bins now run to the largest
observed incidence); `--min-n` 500 → 1 with a degenerate-NMAD guard in its place; the three
count cuts in the cover-band shape check replaced by printed counts; figure rebuilt as a
full-range panel plus a 0–36° detail panel with median cover on a twin axis.
3 — `offset_by_beam_selection.py`: `MIN_BIN = 100` and the `m.sum() >= 50` hexbin cut gone;
bins now carry robust SEs in table and figure.
4, 15 — `offset_model_slope_cover.py` `--min-cells` 30 → 1 and `offset_model.median_surface`
`min_cells` 30 → 1; **the sentences claiming the blanked boxes are "combinations the terrain
does not supply" are deleted from both** — they were false. The script now prints, on every
run, what a `--min-cells 30` rule would cost, so the flattery is never invisible again.
5, 14 — `offset_vs_angle.py`: `MIN_N = 300`, the 500-return cover-band skip and the 50-cell
scatter-check cut gone; `--block-m` added for cluster-robust error bars.
6, 19 — `q2_cover_fit.py`'s hardcoded `n1 >= 5 & ng >= 10` exposed as `--min-gen1/--min-gen2`
at the definitional 1, and `percentile_float_fit.py`'s 3/5 matched to it, closing the
two-populations-from-one-cube comparability gap.
7 — `offset_model_slope_cover.py --min-n` 3 → 1 (the dropped cells were 14.5× enriched in
canopy cover and sat at 30.5° median slope).
9 — `cover_offset_isotonic.py` `min_n=5` → 1, and the y-axis now stretches to the restored
point instead of cropping it off.
11, 12, 16 — the three latent-but-inert minimums (`cover_offset_reference --min-n 200`,
`cover_offset_regression min_n=200`, `matched_band_effects min_n=200`) lowered to their
definitional floors; output verified byte-identical, except that `matched_band_effects` was
in fact biting on elbaext curv0.001 (see the table).

### Kept, with the reason stated

- 8 — `offset_by_beam_selection.py --min-n 3`. This one is **methodological, not tidying**:
  the script picks four *different* end-member beams within each cell, and a cell with one or
  two returns cannot supply four distinct picks. It is exposed, documented in the docstring,
  and named in the figure label.
- 17 — `coreg.py` `min_pts=2000` per swath and `if ok.sum() < 10`. Verified inert (0 of 4 and
  0 of 6 in-grid swaths skipped) and they *guard* an estimator rather than deleting data —
  the fallback is a zero drift curve, not a dropped observation. Left untouched deliberately:
  this is the registration path, and every `d_mm_corr` in the project descends from it, so a
  change here would invalidate the products all of the above are measured against. Flagged
  rather than altered.
- 18 — `pipeline.py:289 if m.sum() < 500: return None`. A fallback guard on the xdem
  heteroscedastic LoD model; returns None and falls back, deletes nothing.
- The ~120 sites in the appendix: left as historical record, per this audit's own §"The rest".

### Not committed

`analysis/ridgelines/percentile_float_fit.py` is **untracked** in this repo. Its fix (item 19)
is applied in the working tree and verified, but it has not been added to git — adding a new
file to the history is Andy's call, not a filter fix.

### Figures regenerated

`nearnadir_slope_dependence{,_elbaext}.png` and both `NEARNADIR_SLOPE_DEPENDENCE*.md` ·
`incidence_correction_elbaext.png` (newly built; it was not on disk) ·
`offset_by_beam_selection_elbaext{,_sp5}.png` · all 6 `offset_model_slope_cover_*_ridge.png`
and their 6 `offset_model_grid_*_ridge.png` siblings · all 38 `offset_vs_*.png` variants
(set verified unchanged: none added, none missed) · `cover_offset_isotonic{,_elbaext}.png` ·
`cover_offset_reference*.png` (5, byte-identical) · `cover_offset_regression{,_elbaext}.png`
(byte-identical) · `q2_vs_cover_fits_se.png`.

### Still to do, outside this repo

The memory note `nearnadir-slope-27deg-knee.md` still records "the ~27° 'switch-on' is REAL
… R²=0.80 F=10.85 … NOT a tan-curve artifact". That statistic is retracted above and the
memory needs the same correction; it was not edited here.
