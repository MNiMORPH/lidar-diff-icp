"""Per-swath 3-D registration as a STAR against gen2, instead of a swath-to-swath CHAIN.

``coreg.align_swaths`` runs Nuth & Kaeaeb on each overlapping swath PAIR and chains the
results through a free network. On elbaext that network is 6 nodes and 5 edges -- a bare
chain, no loops, because non-adjacent swaths have disjoint bounding boxes -- so every
swath's constants are a running sum with NO redundancy and ``max |misclosure|`` is
0.0000 mm by construction. Three independent measurements say its HORIZONTAL half is
wrong:

  * the global Nuth & Kaeaeb removes -0.618 m of the cross-track shift it applies;
  * its horizontal alone makes per-line VERTICAL agreement 70% worse (63.9 -> 108.1 mm);
  * regressing the residual gen1-gen2 difference on the gen2 gradient says 74% of its
    1389 mm cross-track ramp should be undone (corr -0.969, slope -0.738).

WHY THE CHAIN'S HORIZONTAL IS BADLY POSED. It estimates the shift from swath OVERLAP
STRIPS -- a few hundred metres of cross-track extent, the narrowest footprint available --
which is exactly where Nuth & Kaeaeb was measured to lose conditioning. It then weights
all three components by the VERTICAL tie's variance (``coreg.py`` builds one weight per
edge and reuses it for k in range(3)), so the horizontal is weighted by how well the
vertical was determined.

WHAT THIS DOES INSTEAD. For each swath independently, on STABLE cells, fit

    d(cell) = dx * dz/deast + dy * dz/dnorth + c

where ``d`` is that swath's own gridded ground minus the gen2 reference. This is the same
physics as Nuth & Kaeaeb -- a lateral shift shows up as an elevation difference
proportional to the surface gradient -- but as a three-parameter LINEAR fit over the whole
swath footprint rather than an aspect-binned fit over an overlap strip. Measured SEs are
1-4 mm against the chain's spurious ramp.

Being a STAR it has redundancy the chain cannot have: every swath is determined
independently, so swaths CAN disagree, and each carries a standard error that means
something. The absolute level stays gauge -- as it must, since a vertical offset is not an
error -- but the RELATIVE levels no longer propagate through five links from a zero line.

THE COST, stated because it is real: ``align_swaths`` is gen2-FREE, the only consistency
check in the pipeline that does not reference the other epoch. This does reference it. Keep
the chain solve running as a DIAGNOSTIC so that check survives; only stop APPLYING it.

SIGN CONVENTION, established by synthetic recovery (scripts/test_starframe.py), not by
reasoning -- and the test caught that the two halves DIFFER, which the first draft of this
docstring got wrong. Displacing a known cloud and reading the estimator back:

    horizontal   the fitted coefficient comes back as MINUS the applied shift, because a
                 cloud displaced by s reads z_true(X - s) ~ z_true(X) - s*grad
    vertical     the fitted constant comes back as PLUS the applied shift, because the
                 elevations were simply raised

So the raw fit is not uniformly "the correction". This function NEGATES the vertical
coefficient before returning, so that all three components are consistently READY TO ADD
to the coordinates. Recovery on synthetic data is within 2.4% on every component.
"""
from __future__ import annotations

import numpy as np

__all__ = ["fit_star_shifts"]


def fit_star_shifts(x, y, z, source_id, ground, Zref, ground_of, grid, stable,
                    *, smooth_cells=1.2, min_cells=500, iterations=2, verbose=True):
    """Per-swath ``(dx, dy, dz)`` fitted independently against ``Zref`` on stable cells.

    Returns ``(corrections, report)``. ``corrections`` maps swath id -> (dx, dy, dz) in
    metres, ready to ADD to the coordinates. ``report`` carries n, the standard errors and
    R^2 per swath, so a badly determined swath is visible rather than silently chained.

    A swath with fewer than ``min_cells`` usable stable cells gets (0, 0, 0) and is named
    in the report: refusing to move it is honest, inventing a shift from 50 cells is not.

    ``iterations`` -- the relation ``d = s . grad z`` is a FIRST-ORDER linearisation, so one
    pass leaves a second-order residual that scales with the shift and the surface
    curvature. Measured on synthetic recovery of a 600 mm shift, one pass leaves 33.3 mm
    (5.5%) and two leave ~1 mm. Real per-swath shifts here are the same +/-500 mm size, so
    a single pass is not adequate. Each extra iteration costs one re-grid per swath.
    """
    from scipy.ndimage import distance_transform_edt, gaussian_filter

    X0, Y0, res, nx, ny = grid
    Zf = np.asarray(Zref, float).copy()
    bad = ~np.isfinite(Zf)
    if bad.any():                      # fill only to take a gradient; never to compare
        Zf = Zf[tuple(distance_transform_edt(bad, return_distances=False,
                                             return_indices=True))]
    Zs = gaussian_filter(Zf, smooth_cells)
    gnorth, geast = np.gradient(Zs, res)          # dz/dnorth (rows), dz/deast (cols)
    st = np.asarray(stable, bool) & np.isfinite(Zref)

    g = np.asarray(ground, bool)
    sid = np.asarray(source_id)
    Zr = np.asarray(Zref, float)
    corrections, rows = {}, []
    for s in sorted(np.unique(sid[g])):
        m = g & (sid == s)
        if m.sum() < 10:
            corrections[int(s)] = (0.0, 0.0, 0.0)
            rows.append(dict(swath=int(s), n=0, reason="no ground points"))
            continue
        cx = cy = cz = 0.0
        info = None
        for _ in range(max(int(iterations), 1)):
            G = ground_of(x[m] + cx, y[m] + cy, z[m] + cz)   # SAME estimator that made Zref
            use = st & np.isfinite(G)
            n = int(use.sum())
            if n < min_cells:
                info = dict(swath=int(s), n=n, reason=f"< min_cells={min_cells}")
                break
            d = G[use] - Zr[use]
            A = np.column_stack([geast[use], gnorth[use], np.ones(n)])
            coef, *_ = np.linalg.lstsq(A, d, rcond=None)
            r = d - A @ coef
            s2 = float(r @ r) / max(n - 3, 1)
            se = np.sqrt(np.diag(s2 * np.linalg.pinv(A.T @ A)))
            # x, y: fitted = -applied, so the fitted value IS the correction to add.
            # z:    fitted = +applied, so it must be NEGATED to become a correction.
            # Verified by scripts/test_starframe.py; do not "simplify" this asymmetry away.
            cx += float(coef[0]); cy += float(coef[1]); cz += -float(coef[2])
            info = dict(swath=int(s), n=n, dx=cx, dy=cy, dz=cz,
                        se_dx=float(se[0]), se_dy=float(se[1]), se_dz=float(se[2]),
                        r2=float(1.0 - r.var() / d.var()) if d.var() > 0 else float("nan"),
                        last_step_mm=(1000*float(coef[0]), 1000*float(coef[1]),
                                      -1000*float(coef[2])))
        if info is not None and "dx" not in info:
            corrections[int(s)] = (0.0, 0.0, 0.0)
        else:
            corrections[int(s)] = (cx, cy, cz)
        rows.append(info)
    if verbose:
        print("  star frame: per-swath 3-D fit against the gen2 reference, stable cells")
        print(f"    {'swath':>7}{'cells':>9}{'dx (mm)':>10}{'SE':>7}{'dy (mm)':>10}"
              f"{'SE':>7}{'dz (mm)':>10}{'SE':>7}{'R2':>8}")
        for r_ in rows:
            if "dx" not in r_:
                print(f"    {r_['swath']:>7}{r_['n']:>9,}   NOT MOVED: {r_['reason']}")
                continue
            print(f"    {r_['swath']:>7}{r_['n']:>9,}{1000*r_['dx']:>10.0f}"
                  f"{1000*r_['se_dx']:>7.0f}{1000*r_['dy']:>10.0f}{1000*r_['se_dy']:>7.0f}"
                  f"{1000*r_['dz']:>10.1f}{1000*r_['se_dz']:>7.1f}{r_['r2']:>8.3f}", flush=True)
    return corrections, rows
