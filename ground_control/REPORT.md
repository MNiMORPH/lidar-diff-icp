# `ground_control` — what was built, what it measured, and what it could not

**2026-08-27.** Deliverable for `HANDOFF.md` §8. Every number here is pasted from a
`trust/provenance.py` run whose ledger is named beside it; nothing is retyped from memory.

---

## 1. The headline, and it is not the one the brief expected

**CORRECTED 2026-08-28.** An earlier version of this section said the absolute constant
was not measurable at Elba. That was wrong and rested on my own algebra error — see
`FRAME.md`. The constant IS determined and closes with the DoD to 0.0000 mm.

| | value | instruments |
|---|---|---|
| **epoch difference (settled)** | gen1 and gen2 level to **7.74 mm** after the geoid | DoD on stable+open ground **−2.12 mm**; benchmark DG8385 **+73.00 ± 11.00** raw → **+5.62** predicted |
| **gen1 absolute (settled)** | **+58.70 ± 25.89 mm** on our surface | 2008 control, open ground, per-line; closes with the DoD to 0.0000 mm |
| **DoD absolute correction** | **+2.12 mm** (negligible) | `c2 − (c1_ours − g)` |

So the DoD needs only **+2.12 mm**, and coregistration was never broken. The shipped
`z_before_absolute` constant (+22.7 ± 39.7) is still superseded — cross-epoch control,
chained 5–6 links — but it is replaced by a determined number, not by a shrug.

Full reasoning: `FRAME.md`, and the memory `elba-absolute-not-measurable`.

---

## 2. Public API

```python
import control          # epoch-agnostic access to both control tables
import lines            # gen1 flight-line tracks, one per PASS
import same_line        # datum from the site's own lines, marks assigned by RETURNS
import datum            # kriged residual-field datum, cover as a treatment
import our_surface      # rebuild OUR gen1/gen2 surface anywhere a tile is on disk
import gen2_swath_deviation   # gen2's per-swath vertical deviation
```

| entry point | answers |
|---|---|
| `control.load_control(epoch, surface=)` | both tables into `residual_field.ControlResiduals` |
| `control.verify_sign_convention(epoch, tol_m=)` | which subtraction the residual column IS, per row |
| `lines.derive_tracks(...)` / `load_tracks` | one track per flight-line PASS |
| `same_line.estimate_by_returns(...)` | **the preferred estimator**: every mark in a tile, assigned by its own returns |
| `same_line.collinearity_sigma(A, B)` | are two passes of one psid the same physical line? |
| `our_surface.our_gen1_surface_at(...)` | our reconstructed surface at a coordinate |
| `datum.datum_at_site(...)` / `sweep_treatments` | kriged field prediction, per cover treatment |

Drivers: `run_datum_by_returns.py`, `run_bridge_wide.py`, `run_derive_tracks.py`,
`run_same_line.py`, `run_datum_at_site.py`, `run_bridge.py`, `run_bridge_gen2.py`,
`run_catchment_check.py`, `run_alongtrack_test.py`, `run_line_cover_model.py`,
`run_weighted_datum_test.py`, `run_gen2_swath_deviation.py`, `summarize_bridge.py`,
`run_gen1_elba_answer.py`.

---

## 3. Design decisions, and why

**An adapter, not a fork.** `residual_field`'s estimators already take plain `(x,y,v)`
arrays; only its *edges* are gen1-schema-bound. `control.py` adapts both tables into
`ControlResiduals` and hands them to the untouched machinery. If promoted, delete it and
move the loaders into `residual_field` — see `INTEGRATION.md`.

**`surface=` is required for gen2.** It publishes four (`ql{0,1}_{dem,laz}`) and they are
four different answers; the choice moves the DoD correction by **12.55 mm**.

**Candidates bounded by data, not by a radius.** The catchment was only ever a compute
bound: `assign_line_from_returns` assigns and can only reject. `estimate_by_returns` uses
every mark inside a tile on disk. Removing the radius also removed a confound — widening
481 m → 2000 m added four marks that were **all urban**, a cover shift disguised as
geometry.

**Group by the physical LINE.** `point_source_id` is reused; passes are merged by
collinearity, scaled by the extrapolation's own prediction sd. `gps_time` cannot
substitute: its correlation with the collinearity sigma is **−0.32**, the wrong sign.

**The line is the unit of replication.** Marks under one line share its unknown constant.
Treating them as independent understates the SE by the measured design effect **1.40×** —
which is exactly how the frame's `±13.0` arose from the same data as `±18`.

**Two bridges, never pooled.** Estimator gap (our reduction of the delivered cloud) is
**−6.59 ± 6.56 mm**; processing gap (our gridded surface) is **−4.04 ± 11.12 mm** over 29
open marks. The previously reported "gen1 bridge −7.2 ± 10.8" was the *estimator* gap.

---

## 4. Verified, and how

| claim | how | result |
|---|---|---|
| sign convention, gen1 | re-derived per row | exact on **1004 of 1004**; reverse misses by 1.080 m |
| sign convention, gen2 | all four surfaces | 238/238, 238/238, 157/157, 157/157 |
| de-duplication | exact `(E,N,Z)` | 1004 rows → **963** marks, 41 dups in 39 groups |
| L1O/L10 trap | column vs prefix | **230** vs **209**; 21 digit-zero |
| LCPs carry no residual | all four columns | **0 of 143** |
| our reconstruction reproduces the product | vs shipped `elbaext` grid | **+1.7 mm** open, −30.5 mm vegetated |
| gen2 reconstruction | 24 points vs `z_after` | **−0.04 mm** at zero canopy; −126.07 mm steepest quartile |
| `z_after − dod` recovers gen1 | vs independent gen1 grid | median **+0.489 mm**, NMAD 0.391, 341,239 cells |
| our ties vs the vendor's own | same 8 marks | **−6.71 mm**, matching the estimator gap |
| along-track drift | tie vs arc length | **+0.74 ± 2.37 mm/km**; 15.5 mm/km excluded at **6.23σ** |

**Regression test proven to bite** (`tests/test_control.py`, 12 pass): repointing cover
from the CSV column to the `point_id` prefix fails it with `assert 209 == 230` and takes a
second test with it; file restored byte-identically (md5 `11d096cf…` before and after).

---

## 5. What I could NOT verify

* **gen2's bridge.** Radius spreads 131–715 mm at its four NVA checkpoints. Not truncation
  (107–200 m from edges), not slope (1.22–4.77°), not sparsity (3,733–5,120 returns): the
  sites have engineered microrelief, order-2 fit RMS **182–479 mm**. A siting problem, so
  more checkpoints of the same kind will not fix it.
* **`sd_field` from kriging.** Set by a nugget/sill split short-lag pairs do not identify:
  it ranges **2.97–37.36 mm** across the sweep while the value holds to 1.12 mm.
* **The 2008 control's geoid.** `GEOID03` is a *dataset-level assertion*; the validation
  reports state no datum. Not resolved, and it is a candidate explanation for the ~67 mm
  discrepancy between control and the two difference instruments.
* **gen1's absolute level at Elba.** Three defensible estimators, three answers — though
  the near/far concern that motivated much of that spread was RESOLVED 2026-08-31 as
  per-line structure, not distance (Welch t = −1.435, p = 0.234; confounded with line).
* **The vendor bias adjustments are absorbed, not unverifiable.** `c1` measures what
  remains after them, so their values are never needed. Hold-out verified: the 963 published
  residuals give mean −43.41 mm, t = −10.17, p = 3.860e-23 against zero. What stays
  unverified is only their spatial uniformity, which the documentation does not state.

---

## 6. Every parameter I chose, and its measured effect

| parameter | value | effect |
|---|---|---|
| cover treatment | `L1O` (Andy's call) | **+17.17 mm** vs `L1O+L5U`; also collapsed the σ sensitivity 8.69 → 1.90 mm |
| gen2 `surface` | reported, never chosen | **12.55 mm** on the DoD correction |
| `collinear_sigma` | swept 2/3/5 | **8.69 mm** on `L1O+L5U`; **1.90 mm** on `L1O` |
| catchment radius | **removed** | had added 4 urban marks and shifted the cover mix |
| `csf_half_width_m` | 300 m | not isolated; validated end-to-end at +1.7 mm |
| bridge radii | swept 10/12.5/15/20 m | median per-mark spread **30.9 mm**, max 239.1 |
| variogram sweep | `max_lag` × estimator | value range 32.67 mm (gen1), 1.12 mm (gen2) |
| `block_m` | 50 m | inherited from the repo, not chosen |
| `align_res`, `swath_tie` | 2.0, `intercept` | inherited from `pipeline.difference_dem` |

---

## 7. Errors made and corrected in-session

Recorded because the corrections are more reusable than the conclusions.

1. **Read heading drift as a gap.** Claimed passes sharing a psid were different lines from
   an 800 m across-track separation; a line at heading 179.3° drifts ~1.1 km over 94 km.
   Fixed by scaling the miss by the extrapolation's own prediction sd.
2. **Silent estimator failure.** Swept radii 7.5/10 m on a 5 m grid; 7.5 m encloses 4 cells
   against the 6 an order-2 fit needs, returned NaN, was skipped, and printed
   `radius_spread = 0.0` on all 29 marks — reading as certainty. A failed radius is now
   reported, never dropped.
3. **The geoid trap.** Differenced a GEOID18-framed surface against GEOID03 control: raw
   gap −59.7 mm of which **+69.3 mm** was the conversion.
4. **Conflated "found near a track" with "belongs to that line."** Cost ~9 mm.
5. **Trusted my own along-track null too far.** It licensed distant same-line marks; the
   open-only set shows +17.20 mm inside 10 km against +76.33 beyond.

The trust gate caught three typed-not-pasted numbers and one missing provenance banner.
