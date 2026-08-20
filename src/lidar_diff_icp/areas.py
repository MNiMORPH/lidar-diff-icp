"""Area building for spatial-coherence change detection.

**Area building is its own thing, decoupled from the count/significance math.**
Wheaton et al. (2010) assess coherence over a fixed 5x5 *square* window; but their
statistic is COUNT-based -- the coherence weight depends only on how many of the
``n`` neighbourhood cells share the change sign, not on the window's geometry (see
:func:`lidar_diff_icp.coherence.spatial_coherence_probability`). So the *shape* of
the neighbourhood is a free design choice: this module builds the shape, and the
Wheaton count math scores it, unchanged.

The isotropic square is one shape (:func:`lidar_diff_icp.coherence.isotropic_counts`).
The other, here, is the **flow corridor**: a footprint grown from a seed cell along
the DEM's flow line -- ``k`` cells downstream and ``k`` upstream -- optionally
widened perpendicular to flow so a channel wider than one cell is captured. This
gathers coherence *along flow*, recovering gullies/rills/channel change that the
isotropic square suppresses (a 1-cell-wide gully has too few same-sign neighbours in
a square window, so its Wheaton posterior collapses regardless of amplitude).

The area is a *footprint over which same-sign cells are counted* -- it is NOT sign-
gated, so the count is a genuine coherence test (a corridor with mixed signs fails).
Which cells to build footprints for is the caller's choice (``cells``); the natural
choice is the amplitude-significant seeds (``|DoD/perror| > z``), leaving the flow
connectivity to supply the coherence -- the flow analogue of Wheaton's logic.

Routing is injected: ``flowdown``/``flowup`` are per-cell flat indices of the
downstream / dominant-upstream neighbour (-1 where none), e.g. from a RichDEM D8 or
D-infinity routing. This module never imports a router, so it stays backend-agnostic.

References: Wheaton, Brasington, Darby & Sear (2010), Earth Surf. Process. Landforms
35(2):136-156, doi:10.1002/esp.1886 (spatial coherence, count-based). Growing the
neighbourhood along flow, and the sediment-connectivity intuition that change is
organised down drainage lines, follow Borselli et al. (2008) / Cavalli et al. (2013)
connectivity work; the flow-corridor footprint itself is our construction, faithful
to Wheaton's geometry-free count.
"""
from __future__ import annotations

import numpy as np


def _corridor(seed, flowdown, flowup, rows, cols, ny, nx, k, width):
    """Flat-index set of the flow-corridor footprint around cell ``seed``:
    ``k`` steps down (``flowdown``) + ``k`` up (``flowup``), plus ``width`` lateral
    cells each side perpendicular to the local flow step (D8-quantised)."""
    cells = {int(seed)}
    for nb in (flowdown, flowup):
        p = int(seed)
        for _ in range(k):
            q = int(nb[p])
            if q < 0:
                break
            cells.add(q)
            p = q
    if width > 0:
        for p in list(cells):
            q = int(flowdown[p])
            if q < 0:
                continue
            dr = int(rows[q] - rows[p]); dc = int(cols[q] - cols[p])   # local flow step
            if dr == 0 and dc == 0:
                continue
            pr, pc = -dc, dr                                            # perpendicular unit
            for sgn in (1, -1):
                for j in range(1, width + 1):
                    rr = int(rows[p]) + sgn * j * pr
                    cc = int(cols[p]) + sgn * j * pc
                    if 0 <= rr < ny and 0 <= cc < nx:
                        cells.add(rr * nx + cc)
    return cells


def flow_corridor_counts(dod, valid, flowdown, flowup, *, k=12, width=0, cells=None):
    """Same-sign COUNT over a FLOW-CORRIDOR footprint, for Wheaton's coherence math.

    For each cell in ``cells`` (default: every valid cell) build the flow corridor
    (:func:`_corridor`) and count how many of its cells are same-sign change. Returns
    ``(ndepos, neros, n)`` grid arrays -- deposition/erosion same-sign counts and the
    footprint size ``n`` -- ready to pass as ``counts=`` to
    :func:`lidar_diff_icp.coherence.spatial_coherence_probability`. ``n`` is per-cell
    (footprints vary in size at grid edges and headwaters), which that function
    handles: its ``low``/``up`` weight bounds scale elementwise with ``n``.

    ``dod``       : DEM of Difference (m), 2-D.
    ``valid``     : boolean mask of usable cells (finite DoD and error).
    ``flowdown``  : per-cell flat index of the downstream neighbour (-1 = none).
    ``flowup``    : per-cell flat index of the dominant-upstream neighbour (-1 = none).
    ``k``         : cells followed each way along flow (corridor length 2k+1).
    ``width``     : lateral cells each side (0 = single-cell-wide thalweg).
    ``cells``     : which cells to build footprints for -- a boolean mask (grid or
                    flat) or an array of flat indices; ``None`` = all valid cells.
                    Cells not built keep count 0 (Wheaton weight 0 -> not promoted).
    """
    dod = np.asarray(dod, float)
    valid = np.asarray(valid, bool)
    ny, nx = dod.shape
    N = nx * ny
    dodf = dod.ravel(); vf = valid.ravel()
    rows = np.arange(N) // nx
    cols = np.arange(N) % nx

    if cells is None:
        idx = np.where(vf)[0]
    else:
        cells = np.asarray(cells)
        idx = np.where(cells.ravel())[0] if cells.dtype == bool else cells.ravel().astype(np.int64)

    ndep = np.zeros(N); nero = np.zeros(N); nn = np.full(N, 2 * k + 1, float)
    for s in idx:
        fp = np.fromiter(_corridor(s, flowdown, flowup, rows, cols, ny, nx, k, width),
                         np.int64)
        vv = vf[fp]; dd = dodf[fp]
        ndep[s] = np.count_nonzero(vv & (dd > 0))
        nero[s] = np.count_nonzero(vv & (dd < 0))
        nn[s] = fp.size
    return ndep.reshape(ny, nx), nero.reshape(ny, nx), nn.reshape(ny, nx)
