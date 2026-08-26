#!/usr/bin/env python3
"""Does the 2008 cross line (point_source_id 10010) break the alternating degeneracy?

``analysis/SWATH_ACROSS_TRACK_TEST.md`` Sec 2 shows the per-line across-track coefficients
``c_s`` are recoverable from N-S swath overlaps only up to an ALTERNATING vector, because
consecutive lines fly there-and-back and every overlap therefore fixes ``tan th_A + tan th_B``
to a constant, so the pair fit returns only the SUM ``(c_A + c_B)/2``.

``analysis/degeneracy_flightline_inventory.py`` found one line whose heading is not N-S:
**psid 10010, heading 271.0 deg**, in the already-local tile ``4342-28-64``.  A due-west
line crosses each N-S line at an across-track position set by the N-S coordinate, which is
INDEPENDENT of the N-S line's own across-track coordinate.  This run measures whether the
two regressors are actually decorrelated in the shared cells, and by how much that shrinks
the standard error on the ANTISYMMETRIC combination ``q = (c_A - c_B)/2`` -- the exact
quantity the N-S chain cannot see.

For every pair the design is ``D = k + p*(tanA - tanB) + q*(tanA + tanB)``, which is the
same fit as ``c_A*tanA - c_B*tanB``; ``p`` is what the N-S overlaps measure and ``q`` is
what they cannot.  ``se_ratio = SE(q)/SE(p)`` from the design matrix alone is the whole
story: it is ~1 when the pair separates the two lines and blows up when it does not.

Ground selection is ``coreg.coregister_swaths``'s own: classification not in (5, 6, 9).
Cell size is the pipeline's 5 m (``data/derived/elbaext/corrections.json`` res_m).
Nothing is fetched and nothing in ``coreg.py`` or ``pipeline.py`` is modified.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/degeneracy_crossline_geometry.py
"""
from __future__ import annotations

import argparse
import os
import sys
import itertools

import numpy as np
import pandas as pd
import laspy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from trust.provenance import Run

EXCLUDE = (5, 6, 9)


def ols_cluster(X, y, groups):
    """OLS with a cluster-robust sandwich covariance -- the same estimator, and the same
    finite-sample factor, as ``analysis/swath_across_track_test.ols_cluster``."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ (X.T @ y)
    e = y - X @ beta
    g = np.asarray(groups)
    o = np.argsort(g, kind="stable")
    gs, Xs, es = g[o], X[o], e[o]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    meat = np.zeros((k, k))
    for i, j in zip(st, en):
        sc = Xs[i:j].T @ es[i:j]
        meat += np.outer(sc, sc)
    G = len(st)
    corr = (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
    return beta, corr * (XtXi @ meat @ XtXi)


def load(path, res, chunk):
    xs, ys, zs, ss, ps, cl = [], [], [], [], [], []
    veg = []          # (x, y) of class-5 returns, for the per-cell vegetation fraction
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(chunk):
            c = np.asarray(pts.classification)
            veg.append((np.asarray(pts.x), np.asarray(pts.y), c == 5))
            m = ~np.isin(c, EXCLUDE)
            xs.append(np.asarray(pts.x)[m].astype(np.float64))
            ys.append(np.asarray(pts.y)[m].astype(np.float64))
            zs.append(np.asarray(pts.z)[m].astype(np.float64))
            ss.append(np.asarray(pts.scan_angle_rank)[m].astype(np.float32))
            ps.append(np.asarray(pts.point_source_id)[m].astype(np.int32))
            cl.append(c[m].astype(np.int8))
    vx = np.concatenate([v[0] for v in veg]); vy = np.concatenate([v[1] for v in veg])
    vm = np.concatenate([v[2] for v in veg])
    return (np.concatenate(xs), np.concatenate(ys), np.concatenate(zs),
            np.concatenate(ss), np.concatenate(ps), np.concatenate(cl), vx, vy, vm)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/before/4342-28-64.laz")
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--min-per-cellline", type=int, default=3)
    ap.add_argument("--block-m", type=float, default=50.0)
    ap.add_argument("--open-vegfrac", type=float, default=0.05)
    A = ap.parse_args()

    R = Run("does the due-west cross line 10010 share ground with the N-S lines at across-track "
            "positions that separate their individual across-track coefficients?")
    R.input(A.tile, role="raw 2008 gen1 LAZ tile containing the cross line 10010 and N-S lines "
                         "135-138, vendor classification")
    R.param("res_m", A.res, src="repo", why="the pipeline grid, corrections.json res_m = 5.0")
    R.param("exclude_classes", EXCLUDE, src="repo",
            why="coreg.coregister_swaths selects ~isin(classification,(5,6,9))")
    R.param("min_per_cellline", A.min_per_cellline, src="repo",
            why="boresight.estimate_boresight's own min_cell_line = 3")
    R.param("block_m", A.block_m, src="repo",
            why="the 50 m cluster-robust block size used by swath_across_track_test.py")
    R.param("open_vegfrac", A.open_vegfrac, src="repo",
            why="cover < 0.05 is the 'open' stratum boundary of swath_across_track_test.py; "
                "here it is the per-cell fraction of returns the VENDOR classified 5, the only "
                "cover proxy available in this tile (no gen2, no PyForestScan cover)")
    R.param("chunk", A.chunk, src="MINE",
            why="points per streamed read, a memory ceiling only; it changes no result because "
                "the per-cell reduction is done after the whole tile is assembled")

    for c, d in [
        ("pair", "unordered flight-line pair (point_source_id)"),
        ("kind", "N-S/N-S = both lines fly the boustrophedon; CROSS = one is the due-west line"),
        ("cells", f"5 m cells covered by BOTH lines with >= {A.min_per_cellline} returns each"),
        ("area_ha", "area of those shared cells, hectares"),
        ("tanA_lo", "minimum per-cell mean tan(scan angle) of the lower-numbered line"),
        ("tanA_hi", "maximum of the same"),
        ("tanB_lo", "minimum per-cell mean tan(scan angle) of the higher-numbered line"),
        ("tanB_hi", "maximum of the same"),
        ("corr_AB", "Pearson correlation of the two lines' per-cell tan(scan angle)"),
        ("sd_sum", "standard deviation of tanA+tanB over the shared cells"),
        ("sd_dif", "standard deviation of tanA-tanB over the shared cells"),
        ("se_ratio", "SE(q)/SE(p) from the design matrix [1, tanA-tanB, tanA+tanB]: how many "
                     "times worse the ANTISYMMETRIC coefficient is determined than the "
                     "symmetric one. 1 = the pair separates the two lines"),
        ("cond", "condition number of that design matrix, column-standardised"),
        ("cls", "vendor LAS classification code"),
        ("count", "returns of the line carrying that classification in the tile"),
        ("psid", "point_source_id, i.e. the flight line, dimensionless integer"),
        ("n_ret", "returns of that line in the tile after the (5,6,9) exclusion"),
        ("frac2", "fraction of that line's returns classified 2 (vendor bare earth)"),
        ("frac12", "fraction classified 12 (vendor overlap/withheld)"),
        ("stratum", "all = every shared cell; open = cells whose vendor class-5 return "
                    "fraction is below open_vegfrac"),
        ("k_mm", "fitted intercept at tanA=tanB=0, mm; positive = line A reads above line B"),
        ("p", "fitted coefficient on tanA-tanB = (c_A+c_B)/2, mm per unit tangent"),
        ("p_se", "cluster-robust standard error of p, mm per unit tangent"),
        ("q", "fitted coefficient on tanA+tanB = (c_A-c_B)/2, mm per unit tangent: the "
              "combination a N-S overlap CANNOT see"),
        ("q_se", "cluster-robust standard error of q, mm per unit tangent"),
        ("c_A", "p+q, the lower-numbered line's own across-track coefficient, mm per unit tangent"),
        ("c_A_se", "cluster-robust standard error of c_A, mm per unit tangent"),
        ("c_B", "p-q, the higher-numbered line's own across-track coefficient, mm per unit tangent"),
        ("c_B_se", "cluster-robust standard error of c_B, mm per unit tangent"),
    ]:
        R.column(c, d)
    R.banner()
    print()

    x, y, z, s, p, c, vx, vy, vm = load(A.tile, A.res, A.chunk)
    x0, y0 = x.min(), y.min()
    nx = int((x.max() - x0) / A.res) + 1
    cell = (((y - y0) / A.res).astype(np.int64)) * nx + ((x - x0) / A.res).astype(np.int64)
    tan = np.tan(np.radians(s.astype(np.float64)))

    vcell = (((vy - y0) / A.res).astype(np.int64)) * nx + ((vx - x0) / A.res).astype(np.int64)
    ncell = int(max(cell.max(), vcell.max())) + 1
    tot = np.bincount(vcell, minlength=ncell).astype(float)
    nveg = np.bincount(vcell[vm], minlength=ncell).astype(float)
    vegfrac = np.divide(nveg, tot, out=np.zeros_like(tot), where=tot > 0)
    blk = ((cell // nx) // int(round(A.block_m / A.res))) * 10_000 + \
          ((cell % nx) // int(round(A.block_m / A.res)))
    df = pd.DataFrame(dict(cell=cell, psid=p, z=z, tan=tan, blk=blk))
    g = df.groupby(["cell", "psid"], sort=True)
    clm = g.agg(z=("z", "median"), tan=("tan", "mean"), n=("z", "size"),
                blk=("blk", "first")).reset_index()
    clm = clm[clm.n >= A.min_per_cellline]

    print("## per-line vendor classification in this tile (after the 5/6/9 exclusion)")
    rows = []
    for ps in sorted(np.unique(p)):
        m = p == ps
        rows.append([int(ps), int(m.sum()),
                     f"{(c[m] == 2).mean():.3f}", f"{(c[m] == 12).mean():.3f}"])
    R.table(["psid", "n_ret", "frac2", "frac12"], rows)
    print()

    print("## pair geometry: can the pair separate c_A from c_B?")
    ids = sorted(clm.psid.unique())
    rows, keep = [], {}
    for a, b in itertools.combinations(ids, 2):
        ca = clm[clm.psid == a][["cell", "tan", "z", "blk"]]
        cb = clm[clm.psid == b][["cell", "tan", "z"]]
        m = ca.merge(cb, on="cell", suffixes=("_a", "_b"))
        if len(m) < 200:
            continue
        ta = m.tan_a.to_numpy()
        tb = m.tan_b.to_numpy()
        X = np.c_[np.ones_like(ta), ta - tb, ta + tb]
        XtXi = np.linalg.pinv(X.T @ X)
        se_ratio = float(np.sqrt(XtXi[2, 2] / XtXi[1, 1]))
        Xs = X[:, 1:] - X[:, 1:].mean(0)
        Xs = Xs / np.linalg.norm(Xs, axis=0)
        cond = float(np.linalg.cond(np.c_[np.ones(len(Xs)) / np.sqrt(len(Xs)), Xs]))
        kind = "CROSS" if 10010 in (a, b) else "N-S/N-S"
        keep[(a, b)] = m
        rows.append([f"{a}-{b}", kind, len(m), f"{len(m)*A.res**2/1e4:.1f}",
                     f"{ta.min():+.3f}", f"{ta.max():+.3f}", f"{tb.min():+.3f}", f"{tb.max():+.3f}",
                     f"{np.corrcoef(ta, tb)[0,1]:+.3f}", f"{(ta+tb).std():.4f}",
                     f"{(ta-tb).std():.4f}", f"{se_ratio:.2f}", f"{cond:.1f}"])
    R.table(["pair", "kind", "cells", "area_ha", "tanA_lo", "tanA_hi", "tanB_lo", "tanB_hi",
             "corr_AB", "sd_sum", "sd_dif", "se_ratio", "cond"], rows)

    cr = [r for r in rows if r[1] == "CROSS"]
    ns = [r for r in rows if r[1] == "N-S/N-S"]
    if cr and ns:
        print(f"\nSE(q)/SE(p): cross pairs {min(float(r[11]) for r in cr):.2f}"
              f"-{max(float(r[11]) for r in cr):.2f}   "
              f"N-S pairs {min(float(r[11]) for r in ns):.2f}"
              f"-{max(float(r[11]) for r in ns):.2f}")
    # ------------------------------------------------------------------ the fit
    print()
    print("## the fit each pair supports: D = k + p*(tanA-tanB) + q*(tanA+tanB),")
    print("## so c_A = p + q and c_B = p - q, in mm per unit tangent.")
    print("## Vertical estimator: per-cell, per-line MEDIAN raw z (the pipeline's ground_q=0.50),")
    print("## on the coregister_swaths class selection. NO CSF and NO gen2 here -- see caveats.")
    fitrows = []
    for (a, b), m in keep.items():
        for lab, sel in (("all", np.ones(len(m), bool)),
                         ("open", vegfrac[m.cell.to_numpy()] < A.open_vegfrac)):
            if sel.sum() < 200:
                continue
            mm = m[sel]
            ta, tb = mm.tan_a.to_numpy(), mm.tan_b.to_numpy()
            D = (mm.z_a.to_numpy() - mm.z_b.to_numpy()) * 1000.0
            X = np.c_[np.ones_like(ta), ta - tb, ta + tb]
            beta, V = ols_cluster(X, D, mm.blk.to_numpy())
            se = np.sqrt(np.diag(V))
            cq = np.array([0.0, 1.0, 1.0]); cb_ = np.array([0.0, 1.0, -1.0])
            fitrows.append([f"{a}-{b}", "CROSS" if 10010 in (a, b) else "N-S/N-S", lab,
                            len(mm), f"{beta[0]:+.1f}", f"{beta[1]:+.1f}", f"{se[1]:.1f}",
                            f"{beta[2]:+.1f}", f"{se[2]:.1f}",
                            f"{cq @ beta:+.1f}", f"{np.sqrt(cq @ V @ cq):.1f}",
                            f"{cb_ @ beta:+.1f}", f"{np.sqrt(cb_ @ V @ cb_):.1f}"])
    R.table(["pair", "kind", "stratum", "cells", "k_mm", "p", "p_se", "q", "q_se",
             "c_A", "c_A_se", "c_B", "c_B_se"], fitrows)
    print()
    print("CAVEATS this fit carries, stated before any number is used:")
    print(" * line 10010 has ZERO vendor class-2 returns (all class 12), so its per-cell median")
    print("   is a mixed ground/vegetation population; a CSF pass is needed before these")
    print("   coefficients are trustworthy. The GEOMETRY result above does not depend on this.")
    print(" * this tile sits N 4,886,209-4,889,709, i.e. 0-3.5 km NORTH of the elbaext grid")
    print("   (N 4,882,200-4,886,250). Whether c_s is constant along track is untested.")
    print(" * no gen2 surface exists here, so no lateral (Nuth & Kaeaeb) or drift term is applied.")
    R.done(headline="cross-line 10010 overlap geometry against the N-S lines")


if __name__ == "__main__":
    main()
