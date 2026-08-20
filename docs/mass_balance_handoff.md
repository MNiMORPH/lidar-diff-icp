# Routed sediment mass-balance — tool handoff

Internal, physics-based validation of a lidar **DEM of Difference (DoD)** by routing
its sediment volume down a flow network and enforcing sediment continuity. Written for
whoever builds the **floodplain-masking component** in parallel, and to seed the tool's
own repository (extracted from `MNiMORPH/lidar-diff-icp`, where it was developed).

## What it does (physics)

Treat the DoD as a morphological sediment budget. Each cell's net volume change is
`V(c) = DoD(c) · area` (+ deposition, − erosion). Accumulating `V` down the flow
network gives, at every cell,

    V_acc(c) = Σ over the upstream area (incl. c) of DoD · area .

By Exner continuity the sediment flux leaving a cell is `Q_out(c) = −V_acc(c)`, and a
channel cannot carry negative sediment, so the **physical constraint** is

    V_acc(c) ≤ 0   ⟺   cumulative erosion ≥ cumulative deposition, everywhere downstream.

A cell where `V_acc` climbs positive **beyond its error** has deposited more than any
upstream erosion can supply → external input (bank collapse, aeolian, anthropogenic
fill, or an **off-map trunk river** — see below) **or** DoD error. Used forward this is
the "morphological method" (Ashmore & Church; Lane et al. 2003; Wheaton 2010 GCD;
Heckmann & Vericat 2018). Here it is **inverted** as a QC / error detector on the DoD
itself — which appears novel (Heckmann & Vericat is the closest, forward, precedent).

**Validity rests on a closed catchment.** In small/headwater catchments sediment is
sourced locally, so the budget can close (Dietrich & Dunne 1978; Walling 1983). The
main threats are **storage** (colluvial/valley-floor sinks legitimately hold deposition
below eroding cells — Trimble 1983, Coon Creek, *same Driftless region as our pilots*:
~94% of eroded sediment goes to storage, ~6% exported) and **off-map input**.

## The code (`massbalance.py`)

Pure-NumPy core, routing injected (backend-agnostic), mass-conserving by construction:

- `dinf_proportions(dem, *, breach=True)` → `(props, valid)`. Thin RichDEM helper:
  D-infinity flow proportions (Tarboton 1997) after **breaching** (Lindsay 2016,
  preferred over filling on noisy lidar). `props[:,:,0]` is a flag (0 = resolved
  interior, −1 = edge/NoData); `props[:,:,1:9]` are downslope fractions (−1 = none).
- `weighted_accumulation(weight, props, valid, *, exponent=1)` → `(acc, exited)`.
  Accumulates `weight` downslope in the flow fractions via a Kahn topological order;
  mass to an invalid/off-grid neighbour EXITS (returned). `exponent=2` for variance.
- `mass_balance(dod, perror, props, valid, res, *, z=1.96)` → dict:
  - `V_acc` — signed accumulated volume (m³), the budget.
  - `sigma_Vacc` — 1σ error envelope. **INDEPENDENT-error LOWER BOUND** — omits flow-
    reconvergence cross-terms and, crucially, **spatial correlation** of the DoD error
    (the systematic bias that dominates a long accumulation). Treat surplus flags as
    candidates until the correlated/`N_eff` envelope (variogram tooling exists in-repo)
    is wired in.
  - `contaminated` — cells whose upstream area touches the domain edge or a data hole
    (excluded; budget cannot close).
  - `surplus` — `V_acc > z·sigma_Vacc & ~contaminated` (unphysical deposition flag).

Tests in `tests/test_massbalance.py` (pure-NumPy mass-conservation + physics; RichDEM-
gated end-to-end). Scope choices in this build: **no bulking** (soil-dominated erosion,
equal density); external inputs ignored; off-map-contaminated cells excluded.

## Why it needs a floodplain mask (the reason for the parallel build)

The check **false-positives on floodplains of large / off-map trunk rivers**: overbank
deposition there is fed by a river entering from off-tile, and **overbank flow is not
steered by local topographic drainage** (Lewin & Ashworth 2014), so the local routing
calls that deposition "unsupported." On the Elba pilot, after switching to the median
ground estimator, the residual surplus concentrates in the **largest drainage (main
valley floor)** — exactly this signature.

**Two hard-won constraints for the mask (from the crude placeholder):**

1. **Handle it at the ROUTING stage, in the right ORDER OF OPERATIONS — not by
   masking the reported output.** The floodplain is the *downstream sink*, not a
   source: it does not feed the accumulation, it *receives* it, so `V_acc` there is a
   real accumulated budget whose sediment came from an **off-map trunk** the local
   routing never saw. Masking floodplain cells from the *output* therefore does
   nothing — the budget already arrived. The fix is to **account for the off-map
   inflow where the routing happens** (a boundary-inflow / off-map-source term at the
   trunk's tile entry), so the sink's budget can close. Order of operations first,
   then the mask.
2. **It must discriminate by SCALE / sourcing, not flatness.** A flat-valley or single-
   HAND cutoff over-masks the headwater valleys we need to keep. The discriminator is
   **sediment sourcing** (locally-sourced vs off-map), i.e. channel-network topology:
   an **on-map channel head** (channel initiates inside the tile) → keep; a channel
   that **enters the boundary already large** (off-map head, trunk) → mask its
   floodplain. Detect off-map trunks by large accumulation at boundary crossings.

**Method options (researched):**
- **Clubb et al. 2017** (*Earth Surf. Dynam.* 5:369–385, doi:10.5194/esurf-5-369-2017)
  — objective floodplain/terrace extraction: dual threshold on local gradient + channel
  relief above a ≥N-Strahler-order channel (HAND-type), thresholds set by Q–Q plots;
  built on DrEICH channel-head extraction (Clubb et al. 2014). LSDTopoTools (C++). The
  Strahler-order reference is the natural trunk-vs-headwater lever. **Recommended
  concept** — pairs with the on/off-map channel-head classification.
- **MRVBF** (Gallant & Dowling 2003, doi:10.1029/2002WR001426) — valley-bottom flatness
  index whose value *is a scale coordinate* (trunk valleys score high, headwaters low);
  DEM-only (no channel network — useful since the trunk is off-tile). Good coarse
  first-pass / cross-check, but flatness ≠ sourcing.
- **HAND** (Rennó 2008; Nobre 2011) — the underlying relief metric; needs a channel
  network and an order/sourcing rule to avoid over-masking.
- **GFPLAIN / V-BET** — drainage-area-scaled; need the trunk's (off-tile) area, so seed
  the network from NHD or a tile-extended DEM.

## Integration points

- The mask feeds two places: (1) the **accumulation** in `weighted_accumulation` /
  `mass_balance` (exclude floodplain/off-map DoD from the sum), and (2) the existing
  **off-map `contaminated`** logic (refine it: exclude only cells fed by an off-map
  *channel*, not every boundary-touching hillslope — recovering evaluable area).
- Routing (`dinf_proportions`) already gives the flow graph the channel-head/HAND work
  needs; reuse it rather than re-routing.

## Open items

- Correlated-error / `N_eff` envelope for `sigma_Vacc` (variogram/xdem tooling exists
  in the origin repo) — to make the surplus % trustworthy.
- The floodplain component (this handoff's subject).
- Validate on a known-change / known-gullied and a known-trunk-floodplain site.
- RichDEM as a clean dependency (dev-checkout version-lookup bug; monkeypatch in use).

## Provenance

Developed in `MNiMORPH/lidar-diff-icp` (Elba 2008-vs-2021 lidar differencing). Key
commits: `massbalance.py` + tests + the Elba driver; extracted here with history.

## References

Ashmore & Church; Lane, Westaway & Hicks (2003) *ESPL* 28:249, doi:10.1002/esp.483 ·
Wheaton, Brasington, Darby & Sear (2010) *ESPL* 35:136, doi:10.1002/esp.1886 (GCD) ·
Paola & Voller (2005) *JGR-ES* 110:F04014, doi:10.1029/2004JF000274 (generalized Exner)
· Heckmann & Vericat (2018) *ESPL* 43:1547, doi:10.1002/esp.4334 (forward DoD routing) ·
Erwin, Schmidt, Wheaton & Wilcock (2012) *WRR* 48:W10512, doi:10.1029/2011WR011035
(budget closure vs measured flux) · Dietrich & Dunne (1978); Trimble (1983) *Am. J.
Sci.* 283:454 (Driftless storage) · Walling (1983) *J. Hydrol.* 65:209 · Borselli et
al. (2008) *CATENA* 75:268; Cavalli et al. (2013) *Geomorph.* 188:31 (connectivity) ·
Tarboton (1997) D-infinity; Lindsay (2016) breaching · Clubb et al. (2017)
doi:10.5194/esurf-5-369-2017; Gallant & Dowling (2003) MRVBF doi:10.1029/2002WR001426;
Rennó (2008)/Nobre (2011) HAND; Lewin & Ashworth (2014) *ESR* 129:1 (trunk-built
floodplain relief).
