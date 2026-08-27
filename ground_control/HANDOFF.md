# `ground_control/` — a generalizable ground-control method for MN lidar vertical datums

**This directory is the permanent home for ground-control code. This file is the brief for the
session that builds it.**

**Written 2026-08-27 by the session that produced the evidence below. Read this whole file
before doing anything.**

---

## 0. YOUR BOUNDARY — read first

**You may READ anything in this repository.** The existing modules and the `analysis/*.md`
reports are your evidence base and you should use them heavily.

**You may WRITE only inside `ground_control/`.** Do not create, edit, move or
delete any file outside this directory. That includes `src/`, `tests/`, `analysis/`,
`data/`, `docs/` and the repository root. Another session owns those and is working in them.

**`ground_control/` is the permanent home for ground-control code in this project** — not a
staging area. Build here as if the code lives here for good: a real module, real tests, a real
driver. It may later be promoted into `src/lidar_diff_icp/` or kept as its own importable
subsystem, but that is a separate decision and nothing you write should depend on it. Keep the
public API small and name things as they would be named in the package either way.

Standing project rules that apply to you:

- **Do not push, do not tag, do not create releases, do not close issues.** Local commits
  only, granular, one logical change each.
- **Invent no thresholds, filters, minimum counts, radii, bin widths or defaults.** Anything
  of that kind is a caller-supplied argument with no default, or is swept and reported as a
  sensitivity. This rule has cost real time on this project when broken: an undisclosed
  `FIT_MAX_SLOPE = 35` once turned a non-significant result into a reported finding.
- **Every number in your write-up is pasted from a command you actually ran.** Never retype a
  number from this document or from another report. If you cannot re-derive it, mark it
  `UNVERIFIED` and name its source.
- **Measurement and interpretation are separate acts.** Give numbers bare; mark any mechanism
  story as hypothesis.
- **A negative result is a result.** Report it straight.
- This is a **shared working laptop**: one heavy job at a time, watch `free -h`, stop if swap
  fills. Do not download 3DEP tiles (0.7–2.3 GB each). gen1 tiles are ~25 MB and 47 are
  already in `data/before/`.

---

## 1. WHAT YOU ARE BUILDING

A reusable method that answers, for **any** Minnesota lidar tile and epoch:

> What constant must be added to this surface to place it on surveyed NAVD88, and with what
> uncertainty *at this location*?

It must work for both generations, return a **prediction at the site** rather than an average
over whatever control happened to be nearby, and state what its uncertainty is the uncertainty
*of*.

Sign convention throughout, and do not deviate: `tie = surveyed − z_lidar`. **Positive means
the surface reads LOW** (add the constant to the surface). Verified for gen1 as
`Control Z − Surface Z`, exact on 1,004 of 1,004 rows; for gen2 as `zdiff = Z − DEMz`, exact
on 395 of 395 rows at 1e-9 m. In both, the reverse ordering misses by up to ~1.1 m.

---

## 2. THE ARCHITECTURE — residual-first, cloud-second

This ordering is the main thing this handoff exists to transmit. It is the opposite of how
the current code is built, and it follows from the measurements in §4.

1. **Look for PUBLISHED per-point residuals** for the site's acquisition. Both MN generations
   publish them, in different places, and this path needs no point cloud, no tiles, no CSF,
   no chaining and no geoid conversion. It is the strongest instrument found. Make it the
   primary route and make "no residuals published" an explicit, reported outcome.
2. **Fit the residual field and predict AT THE SITE**, with a prediction variance from the
   field's correlation structure — not an SE of a sample mean.
3. **Cover enters as a covariate, not a filter.** It is the largest single term in the whole
   problem (§4.3). A method that takes `--cover L1O` as a flag is hiding its dominant
   parameter.
4. **Use the point cloud for exactly two jobs**: (a) the bridge from the *published* surface
   to *our reconstructed* surface, and (b) sites where no residuals were published.
5. **When chaining between flight lines is unavoidable, use LOCAL ties**, not constants
   imported from another tile (§4.5).

---

## 3. COMPONENTS THAT ALREADY EXIST — reuse, do not reimplement

Read all of these. Import them; do not fork them. If one does not fit, say why in your report
rather than silently writing a second copy.

| path | lines | what it is |
|---|---|---|
| `src/lidar_diff_icp/groundtruth/tie.py` | 526 | per-mark tie estimator, radius ladder, ground crops. Solid, tested, epoch-agnostic. **The core primitive.** |
| `src/lidar_diff_icp/groundtruth/gen1_datum.py` | 1092 | discover → resolve tiles → assign line from returns → measure → combine. Good spine, wrong estimator (§4.2). |
| `src/lidar_diff_icp/groundtruth/residual_field.py` | 491 | variogram fit + kriged prediction + LOO/blocked CV on published residuals. **Closest to the right architecture.** |
| `src/lidar_diff_icp/localtie.py` | 839 | ties measured at a stated place; local chaining between lines. |
| `src/lidar_diff_icp/variogram.py` | — | Dowd estimator, `fit_spherical`, `detection_limit`. |
| `src/lidar_diff_icp/groundtruth/chain.py` | — | path planning across the swath graph. |
| `src/lidar_diff_icp/tiles.py` | — | `find_tile`, `download_tile`, `centroid_index`. |
| `analysis/groundtruth/gen1_datum_at_site.py` | — | the existing driver; a good CLI shape to imitate. |
| `analysis/groundtruth/parse_gen2_control.py` | — | how the gen2 tables were parsed. |

**These four modules are essentially unconnected** — `gen1_datum.py` carries `dnr_error_m` as
a dataclass field and does nothing with it, and imports neither `residual_field` nor
`localtie`. Assembling them is most of your job.

### Data already bundled (read-only for you)

- `src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv` — gen1 control. 1,004
  rows, **963 unique** marks after de-duplicating on exact `(easting, northing, elevation)`:
  41 rows in 39 groups, one a triple. Carries the vendor's own `dnr_error_m`.
- `src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv` — gen2 control. 534
  rows = 143 LCP + 227 NVA + 164 VVA, with USGS DEM and LAZ residuals per mark.
- `data/before/` — 47 gen1 tiles, ~25 MB each. 88 of the 963 gen1 marks fall inside them.
- `data/after/` — 11 GB of gen2, covering **only** Elba/elbaext. 6 of 390 gen2 checkpoints.

### Reports to read before designing anything

`analysis/CONTROL_RESIDUAL_FIELD.md` · `analysis/GEN2_ABSOLUTE_DATUM.md` ·
`analysis/GEN1_DATUM_MODULE.md` · `analysis/LOCAL_TIE_CHAINING.md` ·
`analysis/GEN1_OWN_CONTROL_TIE.md` · `analysis/GEN1_DATUM_MORE_MARKS.md` ·
`analysis/DATUM_FROM_MASS_BALANCE.md` · `analysis/ADDITIONAL_GROUND_CONTROL.md` ·
`analysis/LIDAR_DOCUMENTATION_MINE.md`

---

## 4. WHAT WAS MEASURED — the evidence your design must respect

Every number below is in a committed artifact. **Re-derive anything you build on.**

### 4.1 Both epochs publish per-point residuals

- **gen1**: the eight MnGeo county validation reports tabulate `Control Z`, `Surface Z`,
  `Error` for every checkpoint → 963 marks.
- **gen2**: `MN_SE_Driftless_2021_B21/Vertical_Accuracy/USGS/*.dbf`, fields
  `srcChkptId, X, Y, Z, Lndcover, DEMz, zdiff, LAZz, LAZzdiff` → 390 checkpoints, and it
  gives **two** surfaces (delivered DEM and delivered cloud) where gen1 gives one.
- **Both epochs also carry an unpublished vendor bias adjustment** — stated in the metadata,
  value published nowhere, in gen1 *and* gen2. The DoD therefore differences two unknown
  constants. This is a known, permanent floor; do not pretend otherwise.

### 4.2 Mark-averaging is the wrong estimator

The sample mean's SE falls with n while the mean itself moves further, because the quantity
varies spatially. Measured on gen1's published residuals near Elba:

```
  radius  n_all     mean     SE |  n_open+urb     mean     SE
     5km     11    -15.4   27.5 |           3     55.0   46.3
    10km     34    -18.7   12.9 |           9     -3.1   26.0
    15km     66    -31.2   11.5 |          18      1.8   20.8
    20km     99    -54.2   10.4 |          32    -36.6   17.1
    30km    207    -69.5    7.6 |          60    -37.2   12.6
```

Open+urban moves 33 mm between 10 and 20 km against SEs of 17–26 mm. **Your method must
return a prediction at a location, with a prediction variance.**

Related: marks are clustered by flight line, so treating them as independent understates the
SE by a design effect of **1.38–1.42×**. If you average marks anywhere, respect the clustering
and say what the unit of replication is.

### 4.3 Cover is the dominant term, and it is UNRESOLVED

- gen1, kriged to Elba and carried to our surface, by cover treatment: open **−20.4**,
  open+urban **+12.1**, cover-covariate **+50.7** mm — a **71 mm** spread against
  per-treatment prediction sds near 33 mm.
- gen2, delivered cloud vs its own held-out control: NVA (open) **−2.22 ± 2.35** mm,
  VVA (vegetated) **−74.84 ± 6.68** mm. **NVA − VVA = +72.62 ± 7.08, Welch t = +10.252,
  p = 3.94e-20.** Single epoch, so erosion cannot explain it.
- **Three incompatible taxonomies are in play**: MnDNR's five classes
  (`L1O`/`L2T`/`L3B`/`L4F`/`L5U`), USGS's binary NVA/VVA plus an unextracted `Lndcover` code,
  and our own continuous `canopy_cover_pfs`. Any grouping that equates them is a *choice* and
  must be flagged as one.

**Do not try to settle the open scientific question underneath this** (see §6). Your job is to
make the method carry cover explicitly and report its effect, not to pick a value.

### 4.4 The field is real but local

LOO skill against a global constant, scored on a common 230-mark set: 0.110–0.276 (open),
0.167–0.320 (open+urban), 0.369–0.400 (cover-covariate). **Under spatially blocked CV the
skill decays with block size and is effectively gone by 80 km**, one cell negative. The
variogram reaches no sill inside the data — the fitted range pins to the largest lag centre in
9 of 12 rows. **Range is not a quantity this control set determines**; report the sweep and
the empirical points, never a single fitted range.

Note the contrast: **gen2's open-ground residual has almost no spatial structure** (kriging
barely moves it from −2.3 mm), unlike gen1's. Do not assume both epochs need the same model.

### 4.5 Imported swath constants do not transport

- Per-link chaining error **8.42 mm**, verified by re-derivation.
- Pair 136-137 in 400 m windows over 90 km of track: **−57.1 to +28.8 mm, sd 28.3**. The same
  estimator at nine windows a few hundred metres apart inside one tile: **sd 6.0**.
  **F = 22.18 on (5,8), p = 0.0002.** `coreg`'s formal σ for that number is 0.4 mm.
- Local minus imported constants over 12 chained marks: mean +20.8, sd 36.6, **RMS 40.7 mm**.
- **Window size dominates the tie error**: ladder spreads 32.2 / 29.5 / 39.8 mm over
  half-widths 100–1200 m against formal σ of 0.4–4.2 mm. Quote `dz_sigma_window_m`, not the
  formal σ.
- On a window that does not sample `dtan = 0`, the intercept tie becomes a lever arm:
  **+71.7 mm at 100 m half-width vs −17.8 at 1200 m**. `LocalTie.extrapolated` flags it.

### 4.6 The bridge from published surface to ours is not zero

Published residuals describe the **delivered** surface; our DoD uses a CSF-reprocessed,
swath-aligned, geoid-shifted cloud. Measured: gen1 **−7.2 ± 10.8 mm** (18 marks, r = +0.807);
gen2 **+5.3 ± 22.6 mm** (6 marks, r = +0.453). **Estimate and carry this term; never assume
it is zero.** Note gen2's bridge is currently larger than the level it corrects.

### 4.7 The current headline numbers

Delivered-surface against its own contemporaneous control, project-wide:

```
stratum     |  gen1 n     mean     SE |  gen2 n     mean     SE |  gen1-gen2     SE
open        |     230    13.79   7.13 |     227    -2.22   2.35 |      16.01   7.51
vegetated   |     534   -86.51   5.84 |     163   -74.84   6.68 |     -11.67   8.87
```

`gen1 − gen2` on open ground, **+16.0 ± 7.5 mm**, is the current best DoD absolute correction.
The vegetated row depends on a grouping choice (`L2T`+`L3B`+`L4F` = vegetated) that was
flagged as a proposal, not adopted.

---

## 5. NEGATIVE RESULTS — do not rebuild these

- **The radius-spread siting screen does not work.** Three independent measurements: σ_site
  83–97 mm at every cut from 15 mm to none; sd 105.6 / 95.4 / 92.6 / 92.9 / 96.5 with the
  *tightest* cut having the *largest* scatter; 47.7–59.8 flat. Keep the columns as reported
  diagnostics. **Do not filter on them.** An earlier `screen_marks.py` implemented this
  filter; it is deliberately not in the repo and should not be resurrected.
- **The routed mass balance cannot estimate a datum.** Continuity is one-sided
  (`V_acc ≤ +zσ`); raising the surface only makes it more satisfied. It yields a **lower
  bound** and nothing more. Scene-integrated at 1σ: **−24.5 ± 15.2 / −5.8 ± 9.4 /
  −2.3 ± 23.8 mm** for three named aggregates. Consistent with any positive datum; not
  evidence for one.
- **LCPs cannot check gen2.** The 143 LCPs *calibrated* gen2; NVA/VVA were held out. Using
  LCPs is circular. They carry `role=calibration` in the bundled CSV and no residual.
- **Do not parameterise the field by county.** County structure is significant
  (open ground: F = 17.13, p = 4.704e-18 on 230 marks) but within-county sd (91.5 mm) exceeds
  the between-county spread (58.4 mm). County is a diagnostic, not a parameter.

---

## 6. THE OPEN SCIENTIFIC QUESTION — do not settle it, do not foreclose it

The forest−open contrast in the DoD is **degenerate** between (a) gen2 leaf-on canopy bias and
(b) forests genuinely eroding less than farmland. gen2 was flown **2021-05-01 at green-up,
NDVI 0.49** — leaf-on, contradicting its own vendor "leaf-off" spec.

The current state: the two control sets predict the DoD should read **−27.7 ± 11.6 mm** lower
in vegetation than in open; we measure **+23.4 mm** higher. Opposite sign, ~51 mm apart. That
gap runs in the direction of differential erosion — but it compares three taxonomies (§4.3) and
sets project-wide statistics against one tile, so it is **not established**.

Your method must let this question be *asked* — cover as a covariate, both epochs on
comparable strata, the bridge term carried — without answering it by a default.

---

## 7. TRAPS THAT COST TIME TODAY

1. **`L1O` vs `L10`.** 21 gen1 marks have a point_id starting `L10` with a **digit zero**.
   Parsing cover from the `point_id` prefix silently drops them (209 instead of 230 open
   marks). **The CSV's `point_type` column is correct — use it.**
2. **`corrections.json` vs `corrections_geoid.json`.** A tile rebuilt on the geoid datum
   writes the latter and leaves the obsolete `reference_plane` former beside it. Use
   `registration.read_corrections`, which applies the right precedence.
3. **Re-gauging swath constants: `dz` only.** Re-gauging `dx, dy` too shifts the cloud
   horizontally relative to marks at absolute surveyed coordinates. It cost 2.8 mm.
4. **Inverse-variance weighting is invalid when χ² rejects.** gen2's six marks gave a naive
   −25.2 ± 6.5 mm; χ² = 31.9 on 5 dof, p = 6.30e-06. Inflated for σ_site the answer is
   −5 ± 24 mm. **Always test consistency before weighting.**
5. **De-duplicate the gen1 control** on exact `(easting, northing, elevation)`: 1,004 rows,
   963 marks.
6. **`Lndcover` was dropped** when the gen2 CSV was built. It is in the source `.dbf` and is
   free to recover; it gives gen2 a finer cover taxonomy than NVA/VVA.

---

## 8. WHAT "DONE" LOOKS LIKE

Inside this directory:

- A module with a small public API, importable, with the four existing components wired
  together rather than reimplemented.
- Tests, including **at least one regression test you have PROVEN bites**: break the code,
  observe the failure, restore it byte-identically, and report that proof in your write-up.
- A driver with a CLI, in the shape of `analysis/groundtruth/gen1_datum_at_site.py`.
- `REPORT.md`: the API, every design decision with its reason, what you verified and how,
  what you could NOT verify, and **every parameter you chose flagged as YOURS with its
  measured effect**.
- A short `INTEGRATION.md`: how this subsystem is imported by the rest of the project, which
  of its pieces (if any) belong in `src/lidar_diff_icp/` instead, and what would need renaming
  if promoted — so that call is a decision later, not an archaeology exercise.

**Acceptance test:** given a site coordinate and an epoch, the method returns a constant and an
uncertainty *at that site*, states which route produced it (published residual / field
prediction / cloud measurement), names what the uncertainty is the uncertainty of, and reports
the cover treatment's effect on the answer. Run it on Elba (easting 578762.8, northing
4884487.6) for both epochs and show the output.
