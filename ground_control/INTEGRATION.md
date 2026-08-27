# `ground_control` — how it is imported, and what belongs in `src/` instead

Deliverable for `HANDOFF.md` §8. The point is to make promotion a **decision** later
rather than an archaeology exercise.

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
