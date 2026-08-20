"""Mass-conserving internal validation of a DEM of Difference.

Treat the DoD as a morphological sediment budget and route it downhill. Each cell's
net volume change is ``V(c) = DoD(c) * area``; positive = deposition, negative =
erosion. Accumulating ``V`` down the flow network gives, at every cell,

    V_acc(c) = sum over the upstream area (incl. c) of DoD * area .

Sediment continuity (Exner) makes the flux leaving a cell ``Q_out(c) = -V_acc(c)``.
Because a channel cannot carry negative sediment, the PHYSICAL constraint is

    Q_out(c) >= 0   <=>   V_acc(c) <= 0   <=>   cumulative erosion >= cumulative
                                                deposition, everywhere downstream.

A cell where ``V_acc`` climbs positive (beyond its error) has deposited more sediment
than any upstream erosion can supply -- unphysical, i.e. either external input (bank
collapse, aeolian, anthropogenic fill; ignored here by choice) or ERROR in the DoD.
So the routed budget is an internal consistency check that flags artefact-like
deposition and, via the down-network ``V_acc`` profile, diagnoses systematic bias
(a monotonic climb down-network is the signature of a residual tilt/offset).

This is the continuity core of the **morphological method** for sediment budgeting
from repeat topography (Ashmore & Church; Lane et al. 2003; Brasington; Wheaton et
al. 2010, GCD) -- there used forward, to estimate transport; here inverted, as a
quality check on the DoD itself.

**Two decoupled parts.** The core (:func:`weighted_accumulation`, :func:`mass_balance`)
is pure NumPy and takes the flow *proportions* injected, so it is testable without a
router. :func:`dinf_proportions` is a thin RichDEM helper (D-infinity, Tarboton 1997;
breaching preferred over filling for noisy lidar, Lindsay 2016).

**Assumptions / scope (this build):** no bulking (density of eroded ~ deposited --
soil-dominated erosion); external inputs ignored; any cell whose upstream area touches
the domain edge or a data hole is off-map-CONTAMINATED and excluded from the check
(its budget cannot close). The error envelope is the INDEPENDENT-error propagation --
a LOWER BOUND: it omits (a) flow-reconvergence cross-terms and, more importantly, (b)
SPATIAL CORRELATION of the DoD error (the systematic bias that dominates a long
accumulation). Treat flags as candidates, not proof, until the correlated/N_eff
envelope is wired in.
"""
from __future__ import annotations

import numpy as np

# RichDEM D-infinity FlowProportions neighbour convention (band k = 1..8), verified
# empirically (East-ramp test): DX = column offset, DY = row offset.
_DX = np.array([0, -1, -1, 0, 1, 1, 1, 0, -1])
_DY = np.array([0, 0, -1, -1, -1, 0, 1, 1, 1])


def dinf_proportions(dem, *, breach=True, nodata=-9999.0):
    """D-infinity flow proportions for a DEM, via RichDEM. Returns ``(props, valid)``:
    ``props`` is (ny, nx, 9) -- band 0 a flag (0 = resolved interior, -1 = edge/NoData),
    bands 1..8 the downslope fractions (>0 = fraction that way, -1 = none); ``valid``
    is the boolean mask of resolved interior cells (band 0 == 0).

    Depressions are removed by BREACHING (carving through noise dams; Lindsay 2016)
    when ``breach`` -- more faithful than filling on fine noisy lidar -- else by
    epsilon fill. Non-finite DEM cells are set to ``nodata`` before routing.
    """
    import richdem as rd
    rd._RichDEMVersion = lambda: "dev"                      # dev checkout lacks pkg metadata
    z = np.asarray(dem, float).copy()
    z[~np.isfinite(z)] = nodata
    rdem = rd.rdarray(z, no_data=nodata)
    if breach:
        rd.BreachDepressions(rdem, in_place=True)
        rd.FillDepressions(rdem, epsilon=True, in_place=True)   # mop up any residual pits
    else:
        rd.FillDepressions(rdem, epsilon=True, in_place=True)
    props = np.array(rd.FlowProportions(rdem, method="Dinf"))
    valid = props[:, :, 0] == 0
    return props, valid


def weighted_accumulation(weight, props, valid, *, exponent=1):
    """Accumulate ``weight`` down the D-infinity flow network (:func:`dinf_proportions`).

    Each cell passes its accumulated value to downslope neighbours in the flow
    fractions; a cell's accumulation is its own weight plus everything routed in from
    upstream. Processed high-to-low by resolved order so a cell is finalised before it
    distributes. Mass flowing to an invalid cell or off the grid EXITS (returned).

    ``exponent`` scales the routing fraction: 1 for a linear sum (V_acc), 2 for
    variance propagation (fractions enter squared). Returns ``(acc, exited)``.
    """
    props = np.asarray(props, float)
    valid = np.asarray(valid, bool)
    ny, nx = valid.shape
    N = nx * ny
    P = props.reshape(N, 9)
    vf = valid.ravel()
    I = np.arange(N) // nx
    J = np.arange(N) % nx
    acc = np.asarray(weight, float).ravel().copy()
    acc[~vf] = 0.0
    # resolved cells, high elevation first: use each cell's total outgoing fraction as
    # no proxy needed -- topological order = descending of the routing surface. We do
    # not have z here, so order by the fact that flow is a DAG: process by a stable
    # descending sort of a potential. Instead, sort by number of downstream steps is
    # unavailable; use Kahn on the DAG built from props.
    # Build downstream edges and in-degrees for a topological order (robust, no z).
    indeg = np.zeros(N, np.int64)
    down_idx = [[] for _ in range(N)]
    down_frac = [[] for _ in range(N)]
    for c in np.where(vf)[0]:
        for k in range(1, 9):
            f = P[c, k]
            if f <= 0:
                continue
            ni = I[c] + _DY[k]; nj = J[c] + _DX[k]
            if 0 <= ni < ny and 0 <= nj < nx and vf[ni * nx + nj]:
                d = ni * nx + nj
                down_idx[c].append(d); down_frac[c].append(f)
                indeg[d] += 1
    from collections import deque
    q = deque(c for c in np.where(vf)[0] if indeg[c] == 0)
    exited = 0.0
    while q:
        c = q.popleft()
        a = acc[c]
        for d, f in zip(down_idx[c], down_frac[c]):
            acc[d] += (f ** exponent) * a
            indeg[d] -= 1
            if indeg[d] == 0:
                q.append(d)
        # fraction leaving the domain (to invalid / off-grid) -- accounted as exit
        stay = sum(f for f in down_frac[c])
        exited += (1.0 - stay) * a if exponent == 1 else 0.0
    return acc.reshape(ny, nx), exited


def _adjacent_to_invalid(valid):
    """Valid cells that touch (8-connectivity) an invalid cell or the grid edge --
    the places where off-map / data-hole flow would enter unaccounted."""
    v = np.asarray(valid, bool)
    padded = np.zeros((v.shape[0] + 2, v.shape[1] + 2), bool)
    padded[1:-1, 1:-1] = v
    inv = ~padded
    touch = np.zeros_like(v)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            touch |= inv[1 + di:1 + di + v.shape[0], 1 + dj:1 + dj + v.shape[1]]
    return v & touch


def mass_balance(dod, perror, props, valid, res, *, z=1.96):
    """Routed sediment-continuity budget of a DoD (see module docstring).

    ``dod``    : DEM of Difference (m), 2-D (gen2 - gen1).
    ``perror`` : per-cell 1-sigma DoD error (m) = lod / 1.96.
    ``props``, ``valid`` : D-infinity routing from :func:`dinf_proportions`.
    ``res``    : cell size (m); cell area = res**2.

    Returns a dict:
      ``V_acc``       accumulated net volume upstream incl. self (m^3), SIGNED --
                      negative = net erosion (physical), positive = deposition surplus.
      ``sigma_Vacc``  1-sigma independent-error envelope on V_acc (m^3) -- LOWER BOUND.
      ``contaminated`` cells whose upstream area touches the edge / a data hole
                      (excluded: budget cannot close).
      ``surplus``     unphysical-deposition flag: V_acc > z * sigma_Vacc AND not
                      contaminated (deposition beyond what upstream erosion + error allow).
    """
    dod = np.asarray(dod, float); perror = np.asarray(perror, float)
    valid = np.asarray(valid, bool) & np.isfinite(dod) & np.isfinite(perror)
    area = res * res

    vol = np.where(valid, dod * area, 0.0)
    V_acc, _ = weighted_accumulation(vol, props, valid, exponent=1)

    var_cell = np.where(valid, (perror * area) ** 2, 0.0)
    var_acc, _ = weighted_accumulation(var_cell, props, valid, exponent=2)
    sigma_Vacc = np.sqrt(np.maximum(var_acc, 0.0))

    seed = _adjacent_to_invalid(valid).astype(float)
    contam_acc, _ = weighted_accumulation(seed, props, valid, exponent=1)
    contaminated = valid & (contam_acc > 0)

    surplus = valid & ~contaminated & (V_acc > z * sigma_Vacc)
    return {"V_acc": np.where(valid, V_acc, np.nan),
            "sigma_Vacc": np.where(valid, sigma_Vacc, np.nan),
            "contaminated": contaminated,
            "surplus": surplus}
