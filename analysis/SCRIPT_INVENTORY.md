# The 155 scripts under `analysis/`: classified, and why I did not delete by class

Andy asked for "analyse and delete by class". The analysis is done and is reproducible:

    ./lidar-icp/bin/python analysis/classify_scripts.py /tmp/classes.json

**I did not do the deletion, and the reason is evidence, not caution.** The classes that can
be computed do not separate live from dead reliably enough to delete from. I found two
demonstrable false positives by hand in a twelve-item sample — a 17% error rate on a
destructive operation.

## The classification, 2026-09-05

| class | n | rule |
|---|---|---|
| GRAPH | 15 | named in a workflow `Step`'s command or `code` |
| IMPORTED | 5 | another `.py` imports it as a module |
| PRODUCER | 9 | writes a file another script reads |
| CITED | 99 | named in a `.md` |
| ORPHAN | 27 | none of the above |

## Why ORPHAN is not a delete list

**`analysis/ridgelines/epoch_canopy_frac.py` — classified ORPHAN, is a live producer.**
It writes `gen1_canopy_frac.npz` / `gen2_canopy_frac.npz`, which
`refcells.reference_cells` reads for the building and clear-cut cuts. Both files exist on
elbaext and elba_fulldensity. It was missed because it takes its output path as a
command-line argument instead of hard-coding it, so no static scan can attribute the
artifact to it.

**`analysis/ridgelines/plot_dod.py` — classified ORPHAN, is a general utility.** *"Render a
saved DoD in the project's standard figure form."* Nothing imports a plotting tool; that is
what a plotting tool is.

Fifteen of the 27 were last touched on 2026-09-01 or later — they are the figure producers
and control-mark tools from the current ground-q work, several written because Andy asked
for them. They are ORPHAN only because **the graph declares 15 scripts out of 155**.

## What would make this answerable

The reason 140 scripts cannot be classified is that nothing declares them. Fix the
declaration, not the classifier:

1. **Declare the hand-run tools.** Figure producers, audits and one-off validators get a
   `Step` (optional) or a `TOOLS` registry entry naming what they produce and what they are
   for. Then *undeclared* means dead **by construction**, and the sweep is mechanical.
2. **Then sweep.** Anything still undeclared after that is genuinely unreferenced, and
   deleting it is a bookkeeping act rather than a judgement.

Until (1) exists, every deletion here is a per-file judgement — which is exactly what the
penetration sweep was, and it worked because the criterion was a single named artifact, not
a class.

## The lesson the penetration sweep already taught, restated

That sweep inventoried by the ARTIFACT path `penetration.npy` and missed three scripts that
computed the quantity in memory by calling `ground_penetration()`. Same failure shape as
`epoch_canopy_frac` here: **a static scan finds the scripts that name a thing, not the
scripts that do a thing.** Any future sweep must key on the function and the import as well
as the artifact, and must still be read by a human before it deletes.
