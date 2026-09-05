# Data completeness of each site's gen2 cloud

Measured 2026-09-05, `analysis/ept_coverage_check.py --write`, which asks the SOURCE rather
than the file: the 3DEP EPT hierarchy JSON carries a point count per node, so the points
available over a bounding box are summed without downloading any of them.

`ratio` = our file / the AREA-WEIGHTED share of the node set inside the bbox. Nodes clip
the box, so a raw node-count sum overstates what a complete fetch would hold; that
overstatement once read two complete sites as short.

| site | project | in-bbox estimate | our file | ratio |
|---|---|---|---|---|
| battlecreek | MN_CentralMissRiver_5_B22 | 8,455,995 | 8,437,357 | 0.998 |
| carlton | MN_LakeSuperior_2_2021 | 201,587,049 | 199,332,626 | 0.989 |
| elba | MN_SEDriftless_2_2021 | 187,106,542 | 182,923,322 | 0.978 |
| whitewater | MN_SEDriftless_2_2021 | 152,399,301 | 148,050,625 | 0.971 |
| cook | MN_RainyLake_1_2020 | 257,373,163 | 240,897,623 | 0.936 |
| **mnrv** | — | — | — | **UNMEASURABLE, see below** |

**No threshold is applied to these numbers, here or in the code.** What counts as complete
enough for a particular question is a decision to be made out loud. The numbers are recorded
so that decision can be made at all.

## mnrv cannot be measured by this route

`ept_coverage_check.py` fails on mnrv with

    LookupError: the newest project(s) here do not fully cover the bbox

from `threedep.resolve_reference`, which raises rather than silently downgrading to an older
epoch. So mnrv's gen2 bbox spans more than one 3DEP project and no single project covers it.
Its completeness record is therefore ABSENT, and `completeness.check` refuses on it — which
is the correct outcome: unknown is not a pass.

**Correction to the record.** A figure of 0.67 for mnrv's completeness was carried in this
session's notes and repeated in commit 695d580's message. It was relayed from an earlier
summary and NOT recomputed. It does not reproduce: the measurement does not complete for
mnrv at all. Treat mnrv's completeness as unmeasured, not as 0.67.

## What is open

- mnrv needs either a multi-project version of the availability sum, or a stated decision
  that it is measured some other way. It is the only site with no record.
- cook is the lowest of the five that measured, at 0.936. Whether the missing ~16.5 M points
  are a truncated fetch or node clipping at the bbox edge is not established here.
