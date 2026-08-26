#!/usr/bin/env python3
"""WHY the best-fitting percentiles sit at or below the median -- measured, not fitted.

The pairs (q1, q2) that match gen1's near-ground column to gen2's are nearly all below
0.5.  That is the shape of the column, not a fitting quirk.  Nothing scatters from BELOW
the ground: the terrain is a hard floor and every non-ground return in the -1..+2 m window
sits ABOVE it.  So contamination is one-sided-up, and it can only push the ground's rank
DOWN.  Rank 0.5 is the ceiling, reached only when the window holds nothing but ground.

This script measures the rank of each epoch's own ground inside its own column, with no
fitting anywhere:

  * gen2: h = 0 in the cube IS gen2's gridded ground, so F2(0) is the rank directly.
  * gen1: gen1's ground sits at h = -reg, where reg is the per-cell registration
    correction (geoid + lateral + swath + drift) aggregated from beam_offset_table.parquet
    -- established independently of anything in this analysis.

F1 - F2 should then reproduce the fitted q1 - q2 without fitting it.

    ./lidar-icp/bin/python analysis/ridgelines/nearground_rank.py
"""
import numpy as np
import pyarrow.parquet as pq

STRATA = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
          ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01)]


def per_cell_reg_mm(tile, ncell):
    """Mean per-cell registration correction (d_mm_corr - d_mm), mm, gen1 -> gen2 frame."""
    t = pq.read_table(f"data/derived/{tile}/beam_offset_table.parquet",
                      columns=["cell", "d_mm", "d_mm_corr", "in_grid"])
    ing = t["in_grid"].to_numpy().astype(bool)
    cell = t["cell"].to_numpy()[ing]
    reg = (t["d_mm_corr"].to_numpy() - t["d_mm"].to_numpy())[ing]
    s = np.bincount(cell, weights=reg, minlength=ncell)
    c = np.bincount(cell, minlength=ncell)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


for tile in ("elba_fulldensity", "elbaext"):
    D = f"data/derived/{tile}"
    A = np.load(f"{D}/nearground_cells_sn.npz")
    cells, H1, H2 = A["cells"], A["H1"], A["H2"]
    dz = float(A["dz"]); zlo = float(A["zlo"]); NZ = H1.shape[1]
    cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()[cells]
    zf = np.load(f"{D}/z_after.npy")
    reg = per_cell_reg_mm(tile, zf.size)[cells] / 1000.0        # m; gen1 ground at h = -reg

    C1 = np.cumsum(H1, 1).astype(float); C2 = np.cumsum(H2, 1).astype(float)
    n1 = C1[:, -1]; n2 = C2[:, -1]
    k2 = np.full(len(cells), int(round((0.0 - zlo) / dz)))      # gen2 ground bin
    k1 = np.clip(np.round((-reg - zlo) / dz), 1, NZ - 1)        # gen1 ground bin
    k1 = np.where(np.isfinite(k1), k1, 0).astype(int)
    idx = np.arange(len(cells))

    def F(C, n, k):                       # rank of the ground = fraction strictly below it
        return C[idx, k - 1] / np.maximum(n, 1)

    def band(C, n, k, a, b):              # fraction in bins [k+a, k+b-1] i.e. h in [a*dz, b*dz)
        lo = np.clip(k + a, 1, NZ); hi = np.clip(k + b, 1, NZ)
        return (C[idx, hi - 1] - C[idx, lo - 1]) / np.maximum(n, 1)

    ok = (n1 > 0) & (n2 > 0) & np.isfinite(reg) & np.isfinite(cover)
    F1 = F(C1, n1, k1); F2 = F(C2, n2, k2)
    up1 = band(C1, n1, k1, 5, NZ); up2 = band(C2, n2, k2, 5, NZ)      # above ground + 0.10 m
    hi1 = band(C1, n1, k1, 25, NZ); hi2 = band(C2, n2, k2, 25, NZ)    # above ground + 0.50 m
    lo1 = band(C1, n1, k1, -5, 0); lo2 = band(C2, n2, k2, -5, 0)      # 0.10 m just below ground

    print(f"\n=== {tile}  ({cells.size:,} divide cells, {n1.sum()/1e6:.1f}M gen1 + "
          f"{n2.sum()/1e6:.1f}M gen2 returns; median reg {np.nanmedian(reg)*1000:+.1f} mm) ===")
    print("rank of EACH EPOCH'S OWN ground inside its OWN near-ground column (medians over cells)")
    print(f"{'stratum':14s} {'cells':>7s} | {'F1':>5s} {'F2':>5s} {'F1-F2':>6s} |"
          f" {'g1 >+.1':>7s} {'g2 >+.1':>7s} | {'g1 >+.5':>7s} {'g2 >+.5':>7s} |"
          f" {'g1 -.1..0':>9s} {'g2 -.1..0':>9s} | {'ret/cell g1':>11s} {'g2':>5s}")
    rows = [(nm, ok & (cover > lo) & (cover <= hi)) for nm, lo, hi in STRATA]
    rows.append(("ALL", ok))
    for nm, m in rows:
        if m.sum() < 50:
            continue
        print(f"{nm:14s} {m.sum():7,d} | {np.median(F1[m]):5.2f} {np.median(F2[m]):5.2f} "
              f"{np.median(F1[m])-np.median(F2[m]):+6.2f} |"
              f" {np.median(up1[m]):7.2f} {np.median(up2[m]):7.2f} |"
              f" {np.median(hi1[m]):7.2f} {np.median(hi2[m]):7.2f} |"
              f" {np.median(lo1[m]):9.2f} {np.median(lo2[m]):9.2f} |"
              f" {np.median(n1[m]):11.0f} {np.median(n2[m]):5.0f}")
