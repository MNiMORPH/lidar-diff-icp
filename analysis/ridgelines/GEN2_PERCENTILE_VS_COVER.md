# gen2 percentile that matches gen1's ground, as a function of canopy cover

Run 2026-08-26, elba (`data/derived/elba_fulldensity`).

## Inputs
- **gen1**: `beam_offset_table.parquet` (ground from `data/csf_cache/elba.las`, CSF
  rigidness 1 / threshold 1.5), per-cell **median** of `d_mm_corr` — i.e. the four
  registration terms (geoid, lateral, swath, drift) applied per return. q1 fixed at 0.50.
- **gen2**: per-cell near-ground column of the **vendor class-2** returns,
  `nearground_gen2_class_split.npz` (`Hg`), 2 cm bins, quantile interpolated within the bin.
- **cells**: stable reference cells (`lidar_diff_icp.refcells`, all slopes) with >= 5 gen1
  returns and >= 10 gen2 class-2 returns in the -1..+2 m slope-normal window.
- `q2*` solved per bin by Brent root-find on `median(gen1_q50 - gen2_q(q2)) = 0`.

## Table

| cover bin | cells | mean cover | q2* | residual at q2*=fit | residual at q2=0.50 | mm per 0.01 rank |
|---|---|---|---|---|---|---|
| 0.00-0.02 | 31,139 | 0.001 | 0.489 | -0.0 |  -0.7 | 0.7 |
| 0.02-0.05 |  1,604 | 0.035 | 0.496 |  0.0 |  -0.5 | 1.6 |
| 0.05-0.10 |  3,114 | 0.076 | 0.508 |  0.0 |  +1.6 | 1.8 |
| 0.10-0.15 |  4,757 | 0.127 | 0.477 | -0.0 |  -4.5 | 1.8 |
| 0.15-0.20 |  7,537 | 0.177 | 0.484 | -0.0 |  -3.2 | 1.8 |
| 0.20-0.25 | 10,039 | 0.226 | 0.461 | -0.0 |  -6.7 | 1.8 |
| 0.25-0.30 | 11,077 | 0.275 | 0.442 |  0.0 |  -9.8 | 1.7 |
| 0.30-0.35 | 10,256 | 0.324 | 0.430 |  0.0 | -11.7 | 1.7 |
| 0.35-0.40 |  6,920 | 0.372 | 0.428 | -0.0 | -11.9 | 1.6 |
| 0.40-0.50 |  3,921 | 0.432 | 0.424 |  0.0 | -12.4 | 1.7 |
| 0.50-0.65 |    688 | 0.550 | 0.388 | -0.0 | -18.1 | 1.7 |
| 0.65-1.01 |     69 | 0.718 | 0.199 | -0.0 | -59.7 | 2.6 |

Residuals in mm, `gen1 - gen2`; negative means gen1 reads low.

Weighted linear fit: **`q2 = 0.4946 - 0.1749 * cover`**.

## Notes
- Monotone from cover 0.10 up: **0.49 at zero cover -> 0.39 at 0.55 -> 0.20 in the top bin.**
- The linear fit tracks the shape to ~cover 0.5 and then under-corrects badly: it predicts
  0.37 for the top bin, which needs 0.20.
- The first four bins wobble in the wrong direction (0.508 at cover 0.076). Below cover 0.1
  the residual at q=0.50 is only a couple of mm, i.e. at the level of the registration noise.
- `mm per 0.01 rank` is flat at 1.6-1.8 across the whole canopy range, so the curve is not
  an artefact of changing sensitivity -- the required correction genuinely grows with cover.
- Top two bins are 688 and 69 cells. Kept and shown rather than pooled away; treat their
  q2* as indicative.

## Context
The classification side is exhausted: rigidness 1/2/3, threshold 1.5/0.5, cloth resolution
1.0/2.5, ELM, `filters.outlier`, CSF vs vendor TerraScan on gen1 (0.0 mm in every stratum),
and CSF vs class-2 on gen2 all leave this cover-dependent difference intact. It is in the
returns, not in how ground is selected.
