# `ground_control` — how it is imported, and what belongs in `src/` instead

Deliverable for `HANDOFF.md` §8. The point is to make promotion a **decision** later
rather than an archaeology exercise.

> **AMENDED 2026-09-09 (Andy agreed 1-3). THE PURPOSE HAS CHANGED.** This was written
> 2026-08-31 as a tidying document. Two things since have made it something else: every
> pilot site now has gen1 control, and `absolute_datum_mm` is `None` at 6 of 6 while
> README step 6 and `FRAME.md` both call it required. So its job is now **to make the
> datum a routine pipeline step at any Minnesota site**, and the list below should be
> judged against that. The four promotions and three non-promotions still stand; what
> follows amends them.
>
> **1. Everything promoted must be ACQUISITION-AWARE. This document assumes gen1 is one
> survey. It is four.** `lines.py` opens "Flight-line ground tracks for the 2008 gen1
> acquisition", singular, and `data/gen1_line_tracks.json` is keyed by **psid alone**,
> holding only SE-MN lines (115-156, 1512, 1542, 10008-10012; 67 passes over 41 psids).
> elba 135-138 and whitewater 143-146 are in it; mnrv 6191-6251, cook 8-11, carlton 86-90
> and battlecreek 1012-1014/1101-1102 are not. **`psid` is a PER-PROJECT line number**, so
> extending that file statewide keyed by psid alone will eventually merge two different
> flight lines that share a number -- silently. `TrackSet` must be keyed by
> `(project_id, psid)`, mirroring `lidar_diff_icp.acquisitions`. Same for the geoid: gen1
> is four surveys on two geoid models, and applying one survey's frame to another's data
> cost +54.87 mm at Ramsey.
>
> **2. The `sys.path` item in section 3 is now half-done, and the promotions finish it.**
> `trust` was declared a package on 2026-09-09 (Andy: it is SEPARATE from the science
> library, so it keeps its own namespace and no `import trust.provenance` changed); it had
> been importable only because everything ran with the repo root as CWD. That removes the
> 34 hacks that existed only to reach `trust` or the installed package. Of the 35 that
> remain, all sibling imports, these four promotions account for 23: `lines` (9),
> `same_line` (7), `control` (5), `our_surface` (2). The cleanup is a consequence of the
> promotion, not a separate task.
>
> **3. "Do NOT promote `datum.py`" is now MORE right, not less.** Its stated reason was an
> unidentified `sd_field` (2.97-37.36 mm across the sweep), and it was kept as the
> documented fallback for sites whose own lines carry no control. As of 2026-09-09 no site
> is in that position, so the fallback has no remaining customer.
>
> **PROPOSED, not yet agreed -- strike this if unwanted.** Section 4's checklist is five
> procedural steps and contains no GATE. Add one: *the promoted code must reproduce Elba's
> adopted constant -- 58.70 +- 25.89 mm from 8 marks on 5 lines -- before it is trusted at
> a new site.* Reproduce the known case first, then generalise; that is what caught the
> parser's behaviour when it was extended to new acquisitions on 2026-09-09.

---

## 1. How the rest of the project uses it today

`ground_control/` is a flat, importable directory. Nothing in `src/lidar_diff_icp/`
imports it; the dependency runs one way only.

```python
import sys; sys.path.insert(0, "ground_control")
import control, lines, same_line, our_surface, datum, gen2_swath_deviation
```

Each module inserts `src/` on `sys.path` itself, so it works from the repo root without
installation. **That is a development convenience and the first thing to remove on
promotion.**

Consumed from `src/` (imported, never modified):
`groundtruth.tie`, `groundtruth.gen1_datum`, `groundtruth.residual_field`,
`groundtruth.chain`, `coreg`, `io`, `references`, `registration`, `variogram`,
plus `trust.provenance`.

Artifacts other work can use directly, no import needed:

| file | what it is |
|---|---|
| `data/gen1_line_tracks.json` | 67 passes over 41 psids; the object whose scratchpad predecessor evaporated |
| `data/swath_constants_cache.json` | per-tile `align_swaths` constants, so nothing re-reads 25 tiles |
| `products/*.json` | per-site datum, bridge and diagnostic products, each naming its ledger |

---

## 2. What belongs in `src/lidar_diff_icp/`, and what does not

### Promote — these are general capabilities

| module | destination | why | on promotion |
|---|---|---|---|
| `control.py` | **fold into `groundtruth/residual_field.py`** | it exists only because that module's *edges* are gen1-schema-bound while its estimators are already epoch-agnostic. Promotion means making `load_residuals` take an epoch, then **deleting this module** | `load_control` → `residual_field.load_residuals(epoch, surface=)`; keep `verify_sign_convention` |
| `lines.py` | `groundtruth/lines.py` | flight-line tracks are a property of an acquisition, needed at every site | rename `Pass.key` → `pass_id`; keep `INHERITED_PARAMS` and its provenance note |
| `same_line.py` | `groundtruth/same_line.py` | the site-local estimator; the statewide goal needs it everywhere | rename `estimate_by_returns` → `estimate` and drop the old catchment `estimate` (superseded, see §3) |
| `our_surface.py` | `groundtruth/reconstruct.py` | rebuilding our surface at a point is generally useful — it is how any bridge is measured | rename to say what it does; `SurfacePoint.z_geoid18_m/z_geoid03_m` → `z_after_frame_m` / `z_native_frame_m`, which do not hard-code geoid names |

### Do NOT promote

* **`datum.py`** — the kriged residual-field route. Superseded for gen1 by the same-line
  estimator, and its `sd_field` is not identified (2.97–37.36 mm across the sweep). Keep it
  here as the documented fallback for sites whose own lines carry no control; do not put an
  unidentified uncertainty into the package.
* **`gen2_swath_deviation.py`** — a one-off measurement that answered its question (gen2
  ties ≤ 4.8 mm, ~20× tighter than gen1). Its *finding* belongs in the docs; the code does
  not need a home in `src/`.
* **every `run_*.py`** — drivers belong beside `analysis/groundtruth/`, not in the package.

---

## 3. What would need renaming or deleting

* `same_line.estimate` (catchment-based) is **superseded** by `estimate_by_returns` and
  should be deleted, not carried forward — it conflates "found near a track" with "belongs
  to that line" and cost ~9 mm. `discover`, `Scope`, `SEAM_HALF_SPACING_M` and
  `site_scope`'s `pass`/`psid` scopes go with it; only `scope="track"`'s collinearity logic
  survives, as `collinear_groups` / `collinearity_sigma`.
* `SurfacePoint.csf_half_width_m` is gen1-specific and meaningless on the gen2 path, where
  it is reused as a window half-width. Split, or rename to `window_half_width_m`.
* The `sys.path.insert` preambles in every module.
* `run_gen1_elba_answer.py` takes the σ sweep as CLI numbers because it assembles committed
  products rather than recomputing. On promotion it should read the sweep from the product.

---

## 4. Promotion checklist

1. Make `residual_field.load_residuals` epoch-aware; delete `control.py`; keep its tests.
2. Move `lines.py`, `same_line.py` (returns-based path only), `our_surface.py` into
   `groundtruth/`, dropping the `sys.path` preambles.
3. Move `ground_control/data/*.json` to a package-data or `data/derived/` location and
   update `run_*.py`.
4. Run `tests/test_control.py` against the new import paths; the L1O/L10 regression must
   still bite (`assert 209 == 230`).
5. Re-run `run_datum_by_returns.py` and confirm it reproduces its ledgered value.

---

## 5. Standing constraints any consumer must respect

* **Open ground only** for a datum (`L1O` / `NVA`), never pooled — pooling bakes canopy
  response into the level and pre-decides the canopy-vs-erosion question.
* **The line is the unit of replication.** SE over lines, never over marks; the design
  effect is 1.40×.
* **`surface=` is required for gen2** — four delivered surfaces, worth 12.55 mm.
* **The delivered surface and ours are different objects.** A constant measured against
  control applies to the delivered product until the bridge is added.
* **The absolute level is not measurable at one site**; the epoch difference is. See
  `FRAME.md`.
