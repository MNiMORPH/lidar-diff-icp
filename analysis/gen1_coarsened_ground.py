#!/usr/bin/env python3
"""Does COARSENING gen1's ground estimate fix the floodplain, or is it non-penetration?

Andy, 2026-09-23: "Expanding the cellsize for gen1 for ground-return estimate: yes. Let's
do it. And we can use this to help make comparisons with gen2."

THE PREDICTION, STATED BEFORE RUNNING. gen1 carries ~16 near-ground returns per 5 m cell.
A mode estimator needs more: measured in heavy vegetation it is biased +19.7 mm and
scatters 56.4 mm at n=16, against +7.0 and 30.0 at n=64. If gen1's floodplain problem is
SPARSITY, pooling a 3x3 neighbourhood (~150 returns) should let a mode find the ground and
should CUT the DoD's dependence on gen1's own return spread, currently -0.466 mm/mm. If the
slope barely moves, the problem is NOT sparsity but genuine non-penetration, and masking is
the honest answer rather than a better estimator.

WHY POOLING IS LEGITIMATE HERE: ``d_mm_corr`` is gen1's offset from the gen2 surface, taken
SLOPE-NORMAL. It is already detrended, so pooling neighbouring cells mixes their offsets,
not their elevations, and does not smear real relief. That is exactly why the pooling is
done on d and not on z.

WHAT IT COSTS: the ground estimate's effective resolution drops to 15 m while the grid stays
at 5 m. Acceptable on a flat floodplain, NOT on a bank -- so this is reported for the flat
floodplain only and must not be applied where slope is steep.

    ./lidar-icp/bin/python analysis/gen1_coarsened_ground.py
"""
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import uniform_filter

S = "data/derived/elbaext_delong_star_nosurf"
LO, HI, DZ = -1500.0, 1500.0, 20.0            # mm; the d range to histogram, and bin width
edges = np.arange(LO, HI + DZ, DZ)
ctr = 0.5 * (edges[:-1] + edges[1:])
NB = ctr.size

dod = np.load(f"{S}/dod.npy") * 1000.0
sl = np.load(f"{S}/slope.npy")
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sp = np.load(f"{S}/gen1_spread_mm.npy")
st = np.load(f"{S}/stable.npy").astype(bool)
ny, nx = dod.shape

t = pq.read_table(f"{S}/beam_offset_table.parquet",
                  columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]
d = t["d_mm_corr"].to_numpy()[g].astype(float)
ok = np.isfinite(d) & (d > LO) & (d < HI)
cell, d = cell[ok], d[ok]
print(f"gen1 CSF ground returns in range: {cell.size:,}")

H = np.zeros((ny * nx, NB), np.int16)
np.add.at(H, (cell, np.clip(((d - LO) / DZ).astype(int), 0, NB - 1)), 1)
H = H.reshape(ny, nx, NB)
n_own = H.sum(2)
# 3x3 SUM of histograms = pooling the neighbourhood's returns. uniform_filter gives the
# mean, so multiply by 9; only the shape matters for quantiles, but keep counts honest.
Hp = (uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9.0)
n_pool = Hp.sum(2)


def qtl(Hx, q):
    c = np.cumsum(Hx, -1); tot = c[..., -1]
    out = np.full(tot.shape, np.nan)
    m = tot > 0
    out[m] = ctr[(c[m] >= q * tot[m][:, None]).argmax(-1)]
    return out


def hsm_from_hist(Hx):
    """Half-sample mode read off a histogram: recurse on the densest half."""
    c = np.cumsum(Hx, -1); tot = c[..., -1]
    out = np.full(tot.shape, np.nan)
    m = tot > 0
    idx = np.flatnonzero(m.ravel())
    flat = Hx.reshape(-1, NB)
    res = np.empty(idx.size)
    for j, i in enumerate(idx):
        v = np.repeat(ctr, flat[i].astype(int))
        while v.size > 3:
            k = v.size // 2
            w = v[k:] - v[:v.size - k]
            a = int(np.argmin(w)); v = v[a:a + k + 1]
        res[j] = np.median(v)
    out.ravel()[idx] = res
    return out


flat = fp & (sl <= 2.0) & np.isfinite(dod) & np.isfinite(sp) & (n_own >= 3)
print(f"flat floodplain cells: {flat.sum():,}")
print(f"  gen1 returns per cell  own {np.median(n_own[flat]):.0f}   "
      f"pooled 3x3 {np.median(n_pool[flat]):.0f}")

est_own_med = qtl(H, 0.50)
est_pool_med = qtl(Hp, 0.50)
est_pool_hsm = hsm_from_hist(np.where(flat[..., None], Hp, 0))

print("\nDoD-vs-gen1-SPREAD SLOPE -- the prediction is that pooling CUTS it")
print(f"  {'gen1 ground estimator':<34}{'n':>8}{'slope':>9}{'r':>8}{'medianDoD':>11}")
base = dod - np.median(dod[st & np.isfinite(dod)])
for lab, est in (("current: median of own cell", est_own_med),
                 ("pooled 3x3: median", est_pool_med),
                 ("pooled 3x3: half-sample mode", est_pool_hsm)):
    k = flat & np.isfinite(est)
    # DoD = gen2 - gen1. Lowering gen1's ground by delta RAISES the DoD by delta.
    newdod = base[k] + (est_own_med[k] - est[k])
    s_, _ = np.polyfit(sp[k], newdod, 1)
    print(f"  {lab:<34}{k.sum():>8,}{s_:>9.3f}"
          f"{np.corrcoef(sp[k], newdod)[0, 1]:>8.3f}{np.median(newdod):>11.1f}")

print("\nBY gen1-SPREAD QUARTILE (mm) -- where the erosion lived")
k = flat & np.isfinite(est_pool_hsm)
qs = np.percentile(sp[k], [25, 50, 75])
print(f"  {'estimator':<34}" + "".join(f"{f'Q{i+1}':>10}" for i in range(4)))
for lab, est in (("current: median of own cell", est_own_med),
                 ("pooled 3x3: median", est_pool_med),
                 ("pooled 3x3: half-sample mode", est_pool_hsm)):
    row = []
    for lo, hi in zip([0] + list(qs), list(qs) + [1e9]):
        m = k & (sp >= lo) & (sp < hi)
        row.append(np.median(base[m] + (est_own_med[m] - est[m])))
    print(f"  {lab:<34}" + "".join(f"{v:>10.1f}" for v in row))
