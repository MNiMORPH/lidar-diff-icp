# Density-decimation test (M3C2 median)

2021 decimated to 2008 point count (same area => 2008 density). Positive = 2021 higher (deposition). Units metres.

| position | n | dense median | decim median | dense-decim | LoD dense | LoD decim |
|---|--:|--:|--:|--:|--:|--:|
| floodplain (reed canary) | 1,006,979 | +0.007 | +0.008 | -0.000 | 0.074 | 0.101 |
| upland flat | 188,676 | -0.016 | -0.016 | +0.000 | 0.045 | 0.056 |
| steep >15deg | 798,778 | +0.006 | +0.008 | -0.000 | 0.078 | 0.109 |
| ALL | 2,121,908 | -0.003 | -0.003 | -0.000 | 0.063 | 0.085 |

The M3C2 (median) difference is unchanged by decimation everywhere, including the
floodplain (dense - decim = -0.000). Decimation only raises the LoD (fewer points
per cylinder => noisier), so it is not worth doing.

## Why: density effect is real for a LOW ground pick, but not for the median

Per floodplain cell, the 2021 elevation at a given percentile, dense minus
decimated (to 2008 density):

| percentile | dense - decimated (m) |
|---|--:|
| min  | -0.080 |
| 5th  | -0.021 |
| 10th | -0.013 |
| 25th | -0.005 |
| median | +0.000 |

The denser 2021 cloud finds lower ground through the grass -- 8 cm lower at the
minimum, ~2 cm at the 5th percentile -- because more shots have more chances to
reach the surface. The effect decays up the percentiles and is zero at the median.
So a min / low-percentile ground surface WOULD carry a density bias between the
two epochs; the median aggregation used in the product does not.

Interpretation is the reader's; the numbers are what changed when the densities were matched.
