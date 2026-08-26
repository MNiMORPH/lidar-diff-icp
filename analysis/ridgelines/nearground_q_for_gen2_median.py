#!/usr/bin/env python3
"""Which quantile of gen1's near-ground return distribution reproduces gen2's MEDIAN?

Same object and same frame as everything else here: the per-cell near-ground cube
(`nearground_cells_sn.npz`, slope-normal residual from the local gen2 plane), on the
low-concavity ridgeline reference cells, with gen1's four registration terms applied
per cell (`reg`, from `beam_offset_table.parquet`).

For a level q we compare, per cell,

    r_i(q) = Q1_i(q) + reg_i  -  T_i

against two targets T, reported side by side because they are different surfaces:
  * `p50(gen2)` -- the median of gen2's own near-ground column.  THE LITERAL TARGET.
  * `z_after`   -- gen2's gridded ground (h = 0 by construction), i.e. the DoD reference.
    gen2's ground sits at rank ~0.38 of its column, so these differ by design.

q* = argmin median|r(q)|.  The objective is flat, so a single argmin is not a measurement:
the flat range (within 2 mm of the minimum) is reported with it, and q* is recomputed on
two DISJOINT spatial halves (a 50 m block checkerboard) as an independent spread.

    ./lidar-icp/bin/python analysis/ridgelines/nearground_q_for_gen2_median.py [--slope-max 12]
"""
import argparse
import numpy as np
import pyarrow.parquet as pq

from lidar_diff_icp.binstats import block_ids
from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--slope-max", type=float, default=12.0, help="cell slope cut, deg")
ap.add_argument("--gross-change-mm", type=float, default=500.0,
                help="gross-change guard; lowering it toward the LoD is CIRCULAR (see refcells)")
ap.add_argument("--block-m", type=float, default=50.0)
A = ap.parse_args()

STRATA = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
          ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01)]
QS = np.arange(0.02, 0.995, 0.02)
def per_cell_reg_mm(tile, ncell):
    t = pq.read_table(f"data/derived/{tile}/beam_offset_table.parquet",
                      columns=["cell", "d_mm", "d_mm_corr", "in_grid"])
    ing = t["in_grid"].to_numpy().astype(bool)
    cell = t["cell"].to_numpy()[ing]
    reg = (t["d_mm_corr"].to_numpy() - t["d_mm"].to_numpy())[ing]
    s = np.bincount(cell, weights=reg, minlength=ncell)
    c = np.bincount(cell, minlength=ncell)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def quantiles(C, n, zlo, dz):
    """Per-cell quantile of a histogram row, mm, at the bin centre. (len(QS), ncell)."""
    out = np.empty((QS.size, C.shape[0]), np.float32)
    for i, q in enumerate(QS):
        k = (C >= q * n).argmax(1)
        out[i] = (zlo + (k + 0.5) * dz) * 1000.0
    return out


def summarise(res, blk):
    """q*, its flat range, q* on two disjoint spatial halves, and the fit quality at q*."""
    med = np.median(np.abs(res), axis=1)
    i = int(med.argmin())
    flat = QS[med <= med[i] + 2.0]
    ub, inv = np.unique(blk, return_inverse=True)
    half = (inv % 2) == 0                       # 50 m block checkerboard -> disjoint halves
    qa = QS[np.median(np.abs(res[:, half]), axis=1).argmin()]
    qb = QS[np.median(np.abs(res[:, ~half]), axis=1).argmin()]
    r = res[i]
    return (QS[i], flat.min(), flat.max(), min(qa, qb), max(qa, qb),
            med[i], np.sqrt(np.mean(r**2)), np.median(r), res.shape[1], ub.size)


for tile in ("elba_fulldensity", "elbaext"):
    D = f"data/derived/{tile}"
    A_ = np.load(f"{D}/nearground_cells_sn.npz")
    cells, H1, H2 = A_["cells"], A_["H1"], A_["H2"]
    dz = float(A_["dz"]); zlo = float(A_["zlo"]); curv_max = float(A_["curv_max"])
    zf = np.load(f"{D}/z_after.npy"); NX = zf.shape[1]
    cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()[cells]
    reg = per_cell_reg_mm(tile, zf.size)[cells]

    C1 = np.cumsum(H1, 1).astype(np.float32); C2 = np.cumsum(H2, 1).astype(np.float32)
    n1 = C1[:, -1:]; n2 = C2[:, -1:]
    Q1 = quantiles(C1, n1, zlo, dz)
    P50 = quantiles(C2, n2, zlo, dz)[np.searchsorted(QS, 0.50)]        # gen2 column median
    stable, rep = reference_cells(D, cells=cells, curv_max=curv_max, slope_max=A.slope_max,
                                  gross_change_mm=A.gross_change_mm)
    keep = stable & (n1[:, 0] > 0) & (n2[:, 0] > 0) & np.isfinite(reg) & np.isfinite(cover)
    blk = block_ids(cells, NX, 5.0, A.block_m)

    print(f"\n=== {tile}: gen1 quantile that reproduces gen2 ===")
    print(f"stable reference cells: {keep.sum():,} of {cells.size:,} in the cube")
    print("  " + "  ".join(f"{k}: -{v:,}" if k not in ("start", "kept") else f"{k}: {v:,}"
                           for k, v in rep.items()))
    for tgt_name, T in (("p50(gen2) column median", P50), ("z_after (gen2 gridded ground)",
                                                           np.zeros(cells.size, np.float32))):
        print(f"\n  target = {tgt_name}")
        print(f"  {'stratum':14s} {'cells':>7s} {'blk':>5s} | {'q*':>5s} {'flat range':>11s} "
              f"{'halves':>11s} | {'med|r|':>7s} {'RMS':>6s} {'med r':>7s}")
        rows = [(nm, keep & (cover > lo) & (cover <= hi)) for nm, lo, hi in STRATA]
        rows.append(("ALL", keep))
        for nm, m in rows:
            if m.sum() < 100:
                continue
            res = Q1[:, m] + reg[m][None, :].astype(np.float32) - T[m][None, :]
            q, f0, f1, b0, b1, mr, rms, mdr, nc, nb = summarise(res, blk[m])
            print(f"  {nm:14s} {nc:7,d} {nb:5,d} | {q:5.2f} {f0:5.2f}-{f1:4.2f} "
                  f"{b0:5.2f}-{b1:4.2f} | {mr:7.1f} {rms:6.1f} {mdr:+7.1f}")
