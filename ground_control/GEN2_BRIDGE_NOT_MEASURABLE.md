# gen2's bridge is not measurable at the checkpoints we hold

**2026-08-27.** A negative result, recorded because the DoD absolute correction needs
BOTH epochs' bridges and this one is missing for a reason that will not go away by
re-running.

## What was attempted

`ground_control/run_bridge_gen2.py`, ledger `.trust/runs/20260827T195723-3396806.json`.
Our gen2 surface rebuilt from the six per-checkpoint crops already in
`data/after/checkpoints/` (nothing downloaded) and read at the mark, against the
delivered `usgs_ql1_laz_z_m`.

## The producer was validated first, and it passed

24 points inside elbaext where the shipped `z_after.npy` gives an independent answer:

| stratum | n | median local − grid |
|---|---|---|
| zero canopy cover | 8 | **−0.04 mm** |
| slope 1.1–4.0 deg | 6 | −3.70 mm |
| slope 21.1–32.5 deg | 6 | −126.07 mm |

So the reconstruction is faithful on open, low-slope ground and unusable on steep or
vegetated ground.  NVA marks — non-vegetated vertical accuracy — are sited on exactly the
former, so the producer is the right instrument for them.

## It still fails, and not for any of the reasons checked

| point_id | bridge_mm | radius_spread_mm | n_radii |
|---|---|---|---|
| 2024_2021_MN | +133.5 | 131.1 | 4/4 |
| 2036_2021_MN | −22.7 | 202.4 | 4/4 |
| 2099_2021_MN | +45.1 | 293.0 | 4/4 |
| 2210_2021_MN | +153.4 | **714.8** | 4/4 |

A spread of 714.8 mm is not a measurement with an error bar; it is a fit that is not
determined.  Ruled out by direct check:

* **not crop truncation** — every mark is 107–200 m from its nearest crop edge;
* **not steep terrain** — fitted slopes are 1.22–4.77 deg;
* **not sparse ground** — 3,733–5,120 class-2 returns within 15 m.

The cause is **local relief at the sites**: the order-2 surface fit RMS is **182.4, 238.4,
395.7 and 479.2 mm**.  On open low-slope ground gen2's class-2 returns should fit a local
surface to tens of millimetres.  These do not, because the marks are sited on engineered
ground — 2210 is a road shoulder (`ground_control` memory and
`analysis/ABSOLUTE_BASIS_ELBA.md`) — where a curb, ditch or embankment sits inside the
fitting window.  An order-2 surface cannot represent that, so the answer moves with the
radius.

## What this means

**gen2's bridge stays unmeasured, and `bridge_mm` stays `None` on the gen2 products.**
The pooled figure (+77.29 ± 40.81 mm over 4 marks) is NOT reported as gen2's bridge; its
per-mark radius spreads exceed the quantity being measured.

This is a *siting* problem, not a sample-size problem.  More checkpoints of the same kind
would not fix it.  What would:

1. checkpoints on open ground without engineered microrelief, or
2. an estimator that does not assume an order-2 surface over the fitting window — the
   same limitation `tie.py`'s `radius_spread_mm` column exists to expose, and which the
   repo already declines to filter on.

## Consequence for the DoD

The DoD absolute correction differences two constants and needs both bridges. gen1's is
measured (**−4.04 ± 11.12 mm** over 29 open marks,
`.trust/runs/20260827T194848-3391056.json`). gen2's is not. So the correction can be
carried onto our surface on the gen1 side only, and the gen2 side keeps a declared,
unquantified offset from the delivered surface.
