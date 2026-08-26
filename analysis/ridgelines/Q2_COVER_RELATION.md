# q2 = 0.5 − 0.19 · cover

The gen2 percentile whose elevation matches gen1's median ground, as a function of canopy
cover. Elba, 2026-08-26.

## The relationship

    q2(c) = 0.5 - 0.1922 * c          c = canopy cover fraction (0-1)

One parameter. `q2(0) = 0.50` is imposed, not fitted: at zero cover both epochs see the
true ground, so their medians must agree.

Fitted on uniform 0.05 cover bins, weighted by cell count, all 91,121 stable reference
cells, no count threshold. chi2/dof = 0.72.

## Why linear, and why one parameter

Four forms were fitted, all pinned to 0.5 at c = 0:

| form | params | chi2/dof | max residual (mm) |
|---|---|---|---|
| linear `0.5 + b*c` | 1 | 0.72 | 118.5 |
| quadratic `0.5 + b*c + d*c^2` | 2 | 0.72 | 106.1 |
| power `0.5 - b*c^k` | 2 | 0.70 | 108.3 |
| saturating `0.5 - b(1-e^-kc)` | 2 | 0.76 | 118.5 (degenerate -> linear) |

The extra parameters buy nothing. And the linear coefficient is stable under every binning
and weighting we tried, while the power exponent is not:

| variant | b (linear) | k (power) |
|---|---|---|
| original bin edges, SE weights | -0.1949 | 1.508 |
| uniform 0.05 bins, SE weights, min 20 cells | -0.1881 | 1.222 |
| uniform 0.05 bins, cell weights, min 20 cells | -0.1911 | 1.191 |
| uniform 0.05 bins, cell weights, ALL cells | **-0.1922** | 1.246 |
| uniform 0.05 bins, SE weights, ALL cells | -0.5325 | 2.159 |

b sits at -0.19 in every case except the last, where 1- and 2-cell bins with SE of 0.000
and NaN take over the fit. k moves from 1.19 to 2.16 on the same data. **The linear term
is the thing that is actually measured; the curvature is an artefact of binning and
weighting choices.**

The apparent strong curvature came from one point: the old pooled 0.65-1.01 bin, 69 cells
at q2* = 0.199, with enough leverage to set an exponent. Binned uniformly it splits into
0.321 / 0.012 / 0.389 / 0.182 / 0.030 / 0.877 across six bins -- scatter, not a trend.

## Method that produced it

1. **Anchor at theory.** Bare ground -> q2 = 0.50, imposed on every form.
2. **Bin the covariate uniformly**, so no bin gains leverage from how the edges were drawn.
3. **Weight by the data behind each point** (cell count), and keep every cell -- no minimum
   count, which would cut exactly the sparse high-cover regime.

## Inputs

- **gen1**: `beam_offset_table.parquet` (ground from `data/csf_cache/elba.las`), per-cell
  MEDIAN of `d_mm_corr` -- the four registration terms (geoid, lateral, swath, drift)
  applied per return. q1 fixed at 0.50.
- **gen2**: per-cell near-ground column of the vendor class-2 returns,
  `nearground_gen2_class_split.npz`, 2 cm bins, quantile interpolated within the bin.
- **cells**: stable reference cells (`lidar_diff_icp.refcells`, all slopes), >= 5 gen1
  returns and >= 10 gen2 class-2 returns in the -1..+2 m slope-normal window.
- **cover**: `canopy_cover_pfs.npy` (PyForestScan, >2 m, gen2).
- q2* per bin solved by Brent root-find on `median(gen1_q50 - gen2_q(q2)) = 0`.

Reproduce: `./lidar-icp/bin/python analysis/ridgelines/q2_cover_fit.py --binw 0.05 --weight cells`

## Scope and caveats

- One tile (elba). Not yet tested on elbaext.
- Cover range sampled: 0 to 0.93, but 99.9% of cells are below 0.65.
- Sensitivity is ~1.6-1.8 mm of elevation per 0.01 of rank across the canopy range, so the
  full correction from c = 0 to c = 0.5 is about 17 mm.
- This is a matching relation, not a claim about which epoch is correct. It brings gen1 and
  gen2 onto a common surface for differencing; whether gen2 reads high under leaf-on canopy
  remains open.
