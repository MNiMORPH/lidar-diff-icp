"""Scanner boresight (roll) self-calibration from flight-line overlap.

A boresight roll -- the fixed angular misalignment between the laser scanner and the IMU --
leaves an elevation error proportional to **scan angle** that is *the same on every flight
line* (a mounting constant), unlike per-line registration (a translation that differs
line-to-line) or along-track drift (time-varying). Because the vendor supplied no raw
GPS/IMU/trajectory, we cannot re-do the direct georeferencing; we recover the **residual**
boresight empirically from swath overlap -- a strip self-calibration -- and apply it as a
per-point term.

Method. In a grid cell covered by two flight lines the between-line offset difference to a
common reference (e.g. the other epoch) cancels terrain, real change, and datum -- all shared
by the cell -- leaving
        d_A - d_B  =  (reg_A - reg_B)  +  b * (scan_A - scan_B),
so a regression of the per-cell between-line offset difference on the between-line scan-angle
difference gives the common roll ``b`` (slope) and the pairwise vertical registration
(intercept). A least-squares solve over all pairwise intercepts yields per-line registration
offsets (mean-zero gauge). ``b`` is a sensor constant -- calibrate once, reuse across a lift.

Modular by design. :func:`estimate_boresight` reports ``b`` and its uncertainty WITHOUT
touching the data; :func:`apply_boresight` builds the per-point correction from a ``b``.
:func:`lateral_sensitivity` is the convergence self-check -- how much ``b`` moves under a
per-swath horizontal shift -- so we *measure* whether iterating the boresight/lateral
coupling is warranted rather than assuming it.

Units. ``d`` (and ``b``, the correction) are in the same units throughout -- millimetres in
our per-beam table (``d_mm``); ``b`` is then mm per degree of scan angle. Grid-agnostic:
cells are only grouping keys.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class BoresightSolution:
    """Result of :func:`estimate_boresight`."""

    b: float                       # common boresight roll (d-units per degree of scan angle)
    b_se: float                    # bootstrap standard error of b (over overlap cells)
    registration: dict             # flight-line id -> vertical offset (d-units), mean-zero gauge
    pairs: list                    # per-pair dicts: a, b, n_cells, slope, intercept, median_dd
    n_overlap_cells: int
    units: str = "mm"


def _cellline_means(cell, psid, scan, d, min_cell_line):
    """Per-(cell, flight-line) mean offset and mean scan angle, keeping only groups with at
    least ``min_cell_line`` returns (so a line's per-cell mean is stable)."""
    t = pd.DataFrame({"cell": np.asarray(cell), "psid": np.asarray(psid),
                      "sc": np.asarray(scan, float), "d": np.asarray(d, float)})
    t = t[np.isfinite(t.d) & np.isfinite(t.sc)]
    g = t.groupby(["cell", "psid"])
    clm = g.agg(d=("d", "mean"), sc=("sc", "mean"), n=("d", "size")).reset_index()
    return clm[clm.n >= min_cell_line]


def _pair_rows(clm):
    """One row per (cell, unordered line-pair): between-line offset diff dd and scan diff dsc."""
    m = clm.merge(clm, on="cell", suffixes=("_a", "_b"))
    m = m[m.psid_a < m.psid_b].copy()
    m["dd"] = m.d_a - m.d_b
    m["dsc"] = m.sc_a - m.sc_b
    return m


def estimate_boresight(cell, point_source_id, scan_angle, d, *,
                       min_cell_line=3, min_pair_cells=50, n_boot=200, seed=0,
                       units="mm"):
    """Estimate the common boresight roll and per-line registration from flight-line overlap.

    Parameters mirror the per-beam table: ``cell`` (grid cell id per return),
    ``point_source_id`` (flight line), ``scan_angle`` (deg, signed), ``d`` (offset to the
    common reference, e.g. gen2, in ``units``). Returns a :class:`BoresightSolution`.

    ``b`` is the pooled slope of the between-line offset difference on the scan-angle
    difference; ``b_se`` is a bootstrap SE resampling overlap cells (block bootstrap, since
    a cell may contribute several pairs). Pairs with fewer than ``min_pair_cells`` shared
    cells are still pooled but reported for inspection.
    """
    clm = _cellline_means(cell, point_source_id, scan_angle, d, min_cell_line)
    m = _pair_rows(clm)
    if len(m) < 2:
        raise ValueError("insufficient flight-line overlap to estimate boresight")
    lines = sorted(clm.psid.unique())

    # pooled roll + per-pair breakdown
    b_pool, ic_pool = np.polyfit(m.dsc, m.dd, 1)
    pairs, pair_int = [], []
    for (a, b), gg in m.groupby(["psid_a", "psid_b"]):
        if len(gg) < 2:
            continue
        sl, ic = np.polyfit(gg.dsc, gg.dd, 1)
        pairs.append(dict(a=a, b=b, n_cells=int(len(gg)),
                          slope=float(sl), intercept=float(ic), median_dd=float(gg.dd.median())))
        if len(gg) >= min_pair_cells:
            pair_int.append((a, b, float(ic)))

    # per-line registration: least squares on pairwise intercepts, mean-zero gauge
    idx = {p: i for i, p in enumerate(lines)}
    rows, rhs = [], []
    for a, b, ic in pair_int:
        r = np.zeros(len(lines)); r[idx[a]] = 1.0; r[idx[b]] = -1.0
        rows.append(r); rhs.append(ic)
    rows.append(np.ones(len(lines))); rhs.append(0.0)      # sum(reg) = 0 gauge
    reg_vec, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    registration = {p: float(reg_vec[idx[p]]) for p in lines}

    # bootstrap SE of the pooled roll, resampling overlap cells (block bootstrap)
    rng = np.random.default_rng(seed)
    cells = m.cell.to_numpy()
    uniq = np.unique(cells)
    by_cell = {c: np.where(cells == c)[0] for c in uniq}
    dsc = m.dsc.to_numpy(); dd = m.dd.to_numpy()
    boot = np.empty(n_boot)
    for k in range(n_boot):
        pick = rng.choice(uniq, size=uniq.size, replace=True)
        ridx = np.concatenate([by_cell[c] for c in pick])
        boot[k] = np.polyfit(dsc[ridx], dd[ridx], 1)[0]
    b_se = float(np.std(boot))

    return BoresightSolution(b=float(b_pool), b_se=b_se, registration=registration,
                             pairs=pairs, n_overlap_cells=int(uniq.size), units=units)


def apply_boresight(scan_angle, b):
    """Per-return correction to SUBTRACT from the measured elevation to remove the boresight
    roll: ``correction = b * scan_angle`` (same units as ``b``). For z in metres and ``b`` in
    mm/deg, subtract ``apply_boresight(scan_angle, b) / 1000.0``."""
    return b * np.asarray(scan_angle, float)


def lateral_sensitivity(cell_of, x, y, scan_angle, point_source_id, d, shifts, *,
                        res, x0, y0, nx, ny, min_cell_line=3, seed=0):
    """Convergence self-check: how much the estimated roll moves under a per-swath HORIZONTAL
    shift. Re-cells the returns after shifting each flight line by ``shifts[psid] = (dx, dy)``
    (metres) and re-estimates ``b``; returns ``(b0, b_shift, delta)``.

    This probes the boresight<->lateral coupling that the estimation order ("lateral first")
    is meant to handle: if ``|delta|`` is below the estimate's ``b_se``, one pass suffices and
    iterating to convergence is unnecessary. Approximation: it redistributes returns across
    cells but keeps each return's ``d`` fixed (it does not re-evaluate the reference surface at
    the shifted position -- a second-order effect), so it bounds the dominant coupling term.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    psid = np.asarray(point_source_id)
    b0 = estimate_boresight(cell_of, psid, scan_angle, d, min_cell_line=min_cell_line,
                            seed=seed).b
    xs, ys = x.copy(), y.copy()
    for s, (dx, dy) in shifts.items():
        m = psid == s
        xs[m] += dx; ys[m] += dy
    ix = ((xs - x0) / res).astype(np.int64); iy = ((ys - y0) / res).astype(np.int64)
    ing = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    cell2 = np.where(ing, iy * nx + ix, -1)
    ok = cell2 >= 0
    b_shift = estimate_boresight(cell2[ok], psid[ok], np.asarray(scan_angle)[ok],
                                 np.asarray(d)[ok], min_cell_line=min_cell_line, seed=seed).b
    return float(b0), float(b_shift), float(b_shift - b0)
