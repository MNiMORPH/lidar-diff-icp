#!/usr/bin/env python3
"""Do the cross line's overlaps give each gen1 flight line its OWN across-track coefficient?

``analysis/SWATH_ACROSS_TRACK_TEST.md`` Sec 2 measures, on eight N-S swath overlaps, that the
between-line difference is a function of across-track position: ``D = k + c*(tan th_A -
tan th_B)``. Under the per-line model ``err_s(th) = a_s + c_s*tan th``, the there-and-back
flight pattern fixes ``tan th_A + tan th_B`` inside every N-S overlap, so that fit returns
only the SUM ``(c_A + c_B)/2``. Individual ``c_s`` are unidentifiable from the N-S chain.

``analysis/SWATH_DEGENERACY_BREAKING.md`` found the instrument that breaks it: **psid 10010**,
a due-west cross line in ``data/before/4342-28-64.laz``, whose across-track direction is
NORTHING and therefore independent of the N-S lines' own. In the cells it shares with lines
136, 137 and 138 the two tangents correlate at -0.05 to -0.15 rather than -0.92 to -0.96, so
the three-parameter design

    D = k + p*(tan th_A - tan th_B) + q*(tan th_A + tan th_B),   c_A = p + q,  c_B = p - q

is estimable in BOTH coefficients. This run fits it, and then closes the network: with
``c_136``, ``c_137``, ``c_138`` measured individually, the N-S pair sums give ``c_135``,
``c_134``, ``c_133`` in turn, and the pairs 136-137 and 137-138 -- measured independently on
two tiles -- become REDUNDANCY CHECKS the network has never had.

**The one thing that had to change, stated before any number.** The N-S fits reduce with the
pipeline's own per-return quantity ``d_mm_corr`` (the slope-normal residual of a gen1 CSF
ground return to the gen2 reference surface, after the pipeline's registration terms). No
gen2 exists over the cross-line tile, which sits 2.6-3.5 km NORTH of the elbaext grid, so
``d_mm_corr`` cannot be formed there. This run therefore uses the identical estimator with
the reference surface built from the POOLED gen1 ground of the same tile instead of gen2:
per 5 m cell take the median of all lines' ground returns, smooth it by the pipeline's
``sn_smooth_cells``, and read each line's ground as the median of its returns' residuals to
that common tilted plane. In a between-line difference a COMMON reference cancels, which is
the whole reason the overlap test works at all -- but that is an argument, so Sec 2 MEASURES
it: both estimators are run on elbaext's N-S pairs and the coefficients compared.

Nothing in ``coreg.py`` or ``pipeline.py`` is modified; both are imported.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/crossline_fit.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import laspy
from scipy.ndimage import distance_transform_edt as edt, gaussian_filter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lidar_diff_icp import binstats as bs                              # noqa: E402
from trust.provenance import Run                                       # noqa: E402
from swath_across_track_test import ols_cluster, cell_swath_ground, pair_rows  # noqa: E402

CROSS_PSID = 10010


# --------------------------------------------------------------- the estimator
def scan_angle_deg(f):
    """Scan angle in DEGREES from either LAS spelling.

    PDAL's writers.las promotes LAS 1.1/PF1 to 1.4/PF7, which renames ``scan_angle_rank``
    (integer degrees) to ``scan_angle`` (int16 in 0.006 deg units). Reading the wrong one
    silently returns zeros, which has happened twice in this project, so this refuses a
    field that is all zero rather than returning it.
    """
    names = set(f.point_format.dimension_names)
    if "scan_angle" in names:
        sa = np.asarray(f.scan_angle).astype(np.float64) * 0.006
    elif "scan_angle_rank" in names:
        sa = np.asarray(f.scan_angle_rank).astype(np.float64)
    else:
        raise RuntimeError("no scan angle field in this LAS")
    if not (sa != 0).any():
        raise RuntimeError("scan angle field is identically zero -- the LAS lost it")
    return sa


def gen1_cellline_ground(las_path, res, sn_smooth_cells, ground_q, *, grid=None):
    """Per-(cell, flight-line) gen1 ground, in the pipeline's slope-normal form, with the
    POOLED gen1 ground as the common reference surface.

    Mirrors ``pipeline.difference_dem``'s ``groundg`` for ``ground="slope_normal"``
    (pipeline.py:595-606) term for term:

        resid = z - (Zreg[cell] + dE*dzde[cell] + dN*dzdn[cell])
        ground(cell) = the ``ground_q`` quantile of resid in that cell

    with ``Zreg = gaussian_filter(gap-filled per-cell median, sn_smooth_cells)``. The only
    substitution is the surface the residual is taken to: the pipeline uses the gen2 grid
    ``Z21``; here it is the same construction on the gen1 ground itself, pooling every
    flight line. Returns ``(clm, nx, X0, Y0, n_returns)`` where ``clm`` has one row per
    (cell, point_source_id) with the ground estimate in mm, the mean tangent, and n.

    ``grid=(X0, Y0, nx, ny)`` pins the cell indexing to an existing product's grid so the
    cells are the same objects as that product's; otherwise the cloud's own extent is used.
    """
    f = laspy.read(str(las_path))
    cl = np.asarray(f.classification)
    rn = np.asarray(f.return_number)
    nr = np.asarray(f.number_of_returns)
    keep = (cl == 2) & (rn == nr)          # pipeline's `be`: last returns of the CSF ground
    x = np.asarray(f.x)[keep]
    y = np.asarray(f.y)[keep]
    z = np.asarray(f.z)[keep]
    ps = np.asarray(f.point_source_id)[keep]
    sa = scan_angle_deg(f)[keep]
    del f

    if grid is None:
        X0 = np.floor(x.min() / res) * res
        Y0 = np.floor(y.min() / res) * res
        nx = int(np.ceil((x.max() - X0) / res)) + 1
        ny = int(np.ceil((y.max() - Y0) / res)) + 1
    else:
        X0, Y0, nx, ny = grid
    ix = ((x - X0) / res).astype(np.int64)
    iy = ((y - Y0) / res).astype(np.int64)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    x, y, z, ps, sa, ix, iy = (a[ok] for a in (x, y, z, ps, sa, ix, iy))
    cell = iy * nx + ix

    # pooled gen1 ground grid -> gap-filled -> smoothed -> gradients (pipeline.py:505-524)
    med = pd.Series(z).groupby(cell).quantile(ground_q)
    Z = np.full(nx * ny, np.nan)
    Z[med.index.values] = med.values
    Z = Z.reshape(ny, nx)
    nan = np.isnan(Z)
    Zf = Z if not nan.any() else Z[tuple(edt(nan, return_distances=False, return_indices=True))]
    Zreg = gaussian_filter(Zf, sn_smooth_cells)
    dzde = np.gradient(Zreg, res, axis=1).ravel()
    dzdn = np.gradient(Zreg, res, axis=0).ravel()
    Zreg_f = Zreg.ravel()

    dxe = x - (X0 + (ix + 0.5) * res)
    dyn = y - (Y0 + (iy + 0.5) * res)
    resid_mm = (z - (Zreg_f[cell] + dxe * dzde[cell] + dyn * dzdn[cell])) * 1000.0

    df = pd.DataFrame(dict(cell=cell, point_source_id=ps.astype(np.int64),
                           r=resid_mm, tan=np.tan(np.radians(sa)), sc=sa))
    g = df.groupby(["cell", "point_source_id"], sort=True)
    # groupby.quantile is vectorised; an .agg(lambda) over ~1e6 groups is minutes slower
    clm = g[["r"]].quantile(ground_q).rename(columns={"r": "med_corr"})
    clm["tan"] = g.tan.mean()
    clm["sc"] = g.sc.mean()
    clm["n"] = g.r.size()
    return clm.reset_index(), nx, X0, Y0, int(len(df))


def joint_network(clm, nx, res, block_m, ref_line=None):
    """Solve every line's own ``a_s`` and ``c_s`` AT ONCE, with per-cell fixed effects.

    The pairwise fits answer one pair at a time on that pair's own overlap, so a
    disagreement between two of them could always be the two different regions rather than
    the model. This estimator removes that escape: it fits

        g(cell, line) = mu_cell + a_line + c_line * tan(scan)

    to every (cell, line) ground estimate on the tile at once. ``mu_cell`` absorbs terrain,
    cover, and any spatially varying field exactly as differencing within a cell does -- it
    IS differencing within a cell (Frisch-Waugh), generalised from two lines to all of them
    -- so the ``c_s`` are estimated on ONE shared population. Cells covered by a single line
    are dropped: they demean to zero and carry no information about any difference.

    The ``a_s`` are identified only up to a constant (they trade against ``mu_cell``), so one
    line's constant is dropped as the gauge; the ``c_s`` are not affected by that choice.
    """
    keep_cells = clm.groupby("cell").point_source_id.transform("size") >= 2
    d = clm[keep_cells]
    lines = sorted(d.point_source_id.unique())
    L = len(lines)
    ix = {s: i for i, s in enumerate(lines)}
    r = d.point_source_id.map(ix).to_numpy()
    n = len(d)
    ref = ix[ref_line] if ref_line in ix else 0
    keep_col = [j for j in range(L) if j != ref]
    X = np.zeros((n, L + len(keep_col)))
    X[np.arange(n), r] = d.tan.to_numpy()                    # per-line c_s
    for k, j in enumerate(keep_col):                         # per-line a_s, gauge dropped
        X[r == j, L + k] = 1.0
    y = d.med_corr.to_numpy()
    cell = d.cell.to_numpy()
    # within-cell demean (Frisch-Waugh): identical to fitting mu_cell explicitly
    u, inv = np.unique(cell, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    y = y - np.bincount(inv, y)[inv] / cnt[inv]
    for j in range(X.shape[1]):
        X[:, j] -= np.bincount(inv, X[:, j])[inv] / cnt[inv]
    blk = bs.block_ids(cell, nx=nx, res=res, block_m=block_m)
    beta, V, r2, nn, G = ols_cluster(X, y, blk)
    se = np.sqrt(np.diag(V))
    Xs = X / np.maximum(np.linalg.norm(X, axis=0), 1e-30)
    cond = float(np.linalg.cond(Xs))
    return dict(lines=lines, c=beta[:L], c_se=se[:L], V=V, r2=r2, n=nn, G=G,
                cond=cond, cells=int(len(u)), L=L)


def cross_pair_rows(clm, nx, res, block_m):
    """One row per (cell, unordered line pair), the same construction as
    ``swath_across_track_test.pair_rows`` but on the tangent this script already carries."""
    m = clm.merge(clm, on="cell", suffixes=("_a", "_b"))
    m = m[m.point_source_id_a < m.point_source_id_b].copy()
    m["D"] = m.med_corr_a - m.med_corr_b
    m["ta"] = m.tan_a
    m["tb"] = m.tan_b
    m["dtan"] = m.ta - m.tb
    m["stan"] = m.ta + m.tb
    m["blk"] = bs.block_ids(m.cell.to_numpy(), nx=nx, res=res, block_m=block_m)
    m["cx"] = (m.cell.to_numpy() % nx) * res
    m["cy"] = (m.cell.to_numpy() // nx) * res
    return m


def fit_pair(g, *, controls=False):
    """``D = k + p*dtan + q*stan`` (+ optional position controls), cluster-robust on blocks.

    Returns a dict with p, q, their SEs, the derived ``c_A = p+q`` and ``c_B = p-q`` with
    SEs propagated through the full covariance, the SE ratio from the design alone, and
    the two-parameter (N-S-comparable) coefficient ``c_pair``.
    """
    y = g.D.to_numpy()
    blk = g.blk.to_numpy()
    X2 = np.c_[np.ones(len(g)), g.dtan.to_numpy()]
    b2, V2, r2_2, _, _ = ols_cluster(X2, y, blk)
    X = np.c_[X2, g.stan.to_numpy()]
    if controls:
        # QUADRATIC position terms only. Within a flight line the ground coordinate and
        # the scan angle are the same variable (x = x_track + h*tan th, r2 = 0.95-0.997,
        # SWATH_ACROSS_TRACK_TEST Sec 10), so on a cross pair tan_A is linear in EASTING
        # and tan_B is linear in NORTHING: adding LINEAR position controls would be exactly
        # collinear with the two coefficients being estimated and would delete the signal
        # rather than control anything. The quadratic terms are not collinear with them,
        # and they are what a curved residual field would land on.
        u = (g.cx.to_numpy() - g.cx.mean()) / 1000.0
        v = (g.cy.to_numpy() - g.cy.mean()) / 1000.0
        X = np.c_[X, u ** 2, v ** 2, u * v]
    beta, V, r2, n, G = ols_cluster(X, y, blk)
    se = np.sqrt(np.diag(V))
    XtXi = np.linalg.pinv(X.T @ X)
    ca = np.zeros(X.shape[1]); ca[1] = 1.0; ca[2] = 1.0
    cb = np.zeros(X.shape[1]); cb[1] = 1.0; cb[2] = -1.0
    return dict(k=beta[0], k_se=se[0], p=beta[1], p_se=se[1], q=beta[2], q_se=se[2],
                c_A=float(ca @ beta), c_A_se=float(np.sqrt(ca @ V @ ca)),
                c_B=float(cb @ beta), c_B_se=float(np.sqrt(cb @ V @ cb)),
                se_ratio=float(np.sqrt(XtXi[2, 2] / XtXi[1, 1])),
                r2=r2, n=n, G=G, c_pair=b2[1], c_pair_se=float(np.sqrt(V2[1, 1])),
                r2_pair=r2_2)


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cross-las", default="data/csf_cache/4342-28-64.las",
                    help="CSF ground of the cross-line tile (analysis/slope_bias/csf_tiled.py)")
    ap.add_argument("--cross-raw", default="data/before/4342-28-64.laz",
                    help="the raw tile, read only for the per-line class census of Sec 1")
    ap.add_argument("--valid-las", default="data/csf_cache/elbaext.las",
                    help="CSF ground of a tile that ALSO has a beam_offset_table, so the two "
                         "estimators can be compared on the same N-S pairs")
    ap.add_argument("--valid-tile", default="data/derived/elbaext")
    ap.add_argument("--ns-tiles", nargs="+",
                    default=["data/derived/elba_fulldensity", "data/derived/elbaext"])
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--block-m", type=float, default=50.0)
    ap.add_argument("--min-cell-line", type=int, default=1)
    ap.add_argument("--sn-smooth-cells", type=float, default=1.2)
    ap.add_argument("--ground-q", type=float, default=0.50)
    ap.add_argument("--skip-validation", action="store_true")
    A = ap.parse_args()

    R = Run("does the due-west cross line 10010 give each gen1 flight line its OWN "
            "across-track coefficient, and does the resulting network close?")
    R.param("res_m", A.res, src="repo", why="corrections.json res_m = 5.0")
    R.param("ground_q", A.ground_q, src="repo", why="pipeline.difference_dem default")
    R.param("sn_smooth_cells", A.sn_smooth_cells, src="repo",
            why="pipeline.difference_dem default")
    R.param("block_m", A.block_m, src="repo",
            why="the cluster block size of swath_across_track_test.py")
    R.param("min_cell_line", A.min_cell_line,
            src="MINE" if A.min_cell_line == 1 else "repo",
            why="1 is the definitional floor -- a median needs one value, so it imposes NO "
                "cut. It is swath_across_track_test.py's own default; Sec 6 re-runs the "
                "cross-line fits at 3 and 5 and reports how far every coefficient moves.")
    R.param("reference surface for the cross tile", "pooled gen1 ground of the same tile",
            src="MINE",
            why="gen2 does not cover this tile, so the pipeline's gen2 Zreg cannot be formed. "
                "A common per-cell reference cancels in a between-line difference; Sec 2 "
                "measures the substitution on elbaext's N-S pairs instead of asserting it. "
                "It excludes nothing: every cell with ground in it keeps its ground.")
    R.param("last-return selection rn == nr", True, src="repo",
            why="pipeline.difference_dem's `be` for ground_source='csf'")

    for c, d in [
        ("psid", "point_source_id, i.e. the flight line, dimensionless integer"),
        ("n_raw", "returns of that line in the raw vendor tile, all classes"),
        ("frac2", "fraction of them the VENDOR classified 2 (bare earth)"),
        ("frac12", "fraction the vendor classified 12 (overlap/withheld)"),
        ("n_csf", "returns of that line our CSF classified as ground (class 2)"),
        ("n_csf_last", "of those, last returns (rn == nr) -- the pipeline's ground selection"),
        ("csf_rate", "n_csf / n_raw: what fraction of the line's returns CSF calls ground"),
        ("dens", "n_csf_last per square metre over the cells that line covers"),
        ("cells", "5 m cells the pair shares, with >= min_cell_line ground returns each"),
        ("area_ha", "area of those shared cells, hectares"),
        ("pair", "unordered flight-line pair A-B (point_source_id), A < B"),
        ("kind", "N-S/N-S = both lines fly the boustrophedon; CROSS = one is line 10010"),
        ("tile", "derived-product directory or LAS the estimate comes from"),
        ("estimator", "which per-(cell,line) ground quantity was differenced"),
        ("corr_AB", "Pearson correlation of the two lines' per-cell tan(scan angle)"),
        ("sd_sum", "standard deviation of tanA+tanB over the shared cells"),
        ("sd_dif", "standard deviation of tanA-tanB over the shared cells"),
        ("se_ratio", "SE(q)/SE(p) from the design matrix alone: how many times worse the "
                     "ANTISYMMETRIC coefficient is determined than the symmetric one. "
                     "1 = the pair separates the two lines"),
        ("blocks", "independent spatial blocks spanned by the shared cells"),
        ("k_mm", "fitted intercept at tanA = tanB = 0, mm; + = line A reads above line B"),
        ("k_se", "cluster-robust standard error of k_mm, mm"),
        ("c_pair", "two-parameter coefficient on tanA-tanB = (c_A+c_B)/2, mm per unit tangent "
                   "-- the SAME model and estimator as SWATH_ACROSS_TRACK_TEST Sec 2"),
        ("c_pair_se", "cluster-robust standard error of c_pair, mm per unit tangent"),
        ("p", "three-parameter coefficient on tanA-tanB, mm per unit tangent"),
        ("p_se", "cluster-robust standard error of p, mm per unit tangent"),
        ("q", "coefficient on tanA+tanB = (c_A-c_B)/2, the combination a N-S overlap "
              "CANNOT see, mm per unit tangent"),
        ("q_se", "cluster-robust standard error of q, mm per unit tangent"),
        ("c_A", "p+q: the lower-numbered line's OWN across-track coefficient, mm/unit tangent"),
        ("c_A_se", "cluster-robust standard error of c_A, mm per unit tangent"),
        ("c_B", "p-q: the higher-numbered line's OWN across-track coefficient, mm/unit tangent"),
        ("c_B_se", "cluster-robust standard error of c_B, mm per unit tangent"),
        ("r2", "coefficient of determination of that fit"),
        ("t", "coefficient divided by its cluster-robust standard error"),
        ("line", "flight line (point_source_id) whose own coefficient is being reported"),
        ("c_own", "that line's own across-track coefficient, mm per unit tangent"),
        ("c_own_se", "its standard error, mm per unit tangent"),
        ("source", "how that line's coefficient was obtained: measured on a cross pair, or "
                   "propagated along the N-S chain from a measured one"),
        ("predicted", "the pair sum (c_A+c_B)/2 implied by the individually measured "
                      "coefficients, mm per unit tangent"),
        ("observed", "the pair sum measured directly on that N-S overlap, mm per unit tangent"),
        ("resid", "predicted - observed, mm per unit tangent: the redundancy residual"),
        ("resid_sig", "resid divided by the quadrature sum of the two standard errors"),
        ("r_tA_E", "correlation of line A's per-cell tan(scan) with cell EASTING"),
        ("r_tA_N", "correlation of line A's per-cell tan(scan) with cell NORTHING"),
        ("r_tB_E", "correlation of line B's per-cell tan(scan) with cell EASTING"),
        ("r_tB_N", "correlation of line B's per-cell tan(scan) with cell NORTHING"),
        ("kind_check", "whether the row is a genuine redundancy check, i.e. whether both "
                       "lines' coefficients were measured on the cross line rather than "
                       "propagated from this very pair sum"),
        ("variant", "which sensitivity run produced the row"),
        ("d_c136", "change in c_136 from the headline run, mm per unit tangent"),
        ("d_c137", "change in c_137 from the headline run, mm per unit tangent"),
        ("d_c138", "change in c_138 from the headline run, mm per unit tangent"),
        ("d_c10010", "change in the mean c_10010 from the headline run, mm per unit tangent"),
    ]:
        R.column(c, d)

    # ------------------------------------------------------------- inputs exist
    R.input(A.cross_raw, role="raw 2008 gen1 LAZ tile carrying the cross line 10010 and the "
                              "N-S lines 135-138, vendor classification")
    R.input(A.cross_las, role="our CSF ground (class 2) for that tile, from "
                              "analysis/slope_bias/csf_tiled.py at classify_ground_csf defaults")
    if not A.skip_validation:
        R.input(A.valid_las, role="our CSF ground for elbaext -- the same cache the elbaext "
                                  "pipeline product was built from")
        R.input(f"{A.valid_tile}/corrections.json",
                role="elbaext grid geometry, so the two estimators index the SAME cells")

    # ---------------------------------------------------------------- LOAD FIRST
    # Every input is declared before R.banner(), so the banner covers the whole run:
    # cell_swath_ground() declares the parquet, z_after and canopy inputs itself.
    raw = laspy.read(A.cross_raw)
    rps = np.asarray(raw.point_source_id)
    rcl = np.asarray(raw.classification)
    del raw
    csf = laspy.read(A.cross_las)
    cps = np.asarray(csf.point_source_id)
    ccl = np.asarray(csf.classification)
    crn = np.asarray(csf.return_number)
    cnr = np.asarray(csf.number_of_returns)
    cx_ = np.asarray(csf.x)
    cy_ = np.asarray(csf.y)
    csa = scan_angle_deg(csf)
    del csf
    sec1_rows = []
    for ps in sorted(np.unique(rps)):
        mr = rps == ps
        mc = cps == ps
        ml = mc & (crn == cnr)
        if ml.sum():
            cellu = np.unique((np.floor(cy_[ml] / A.res).astype(np.int64) << 20)
                              + np.floor(cx_[ml] / A.res).astype(np.int64))
            dens = ml.sum() / (len(cellu) * A.res ** 2)
        else:
            dens = float("nan")
        sec1_rows.append([int(ps), int(mr.sum()), f"{(rcl[mr] == 2).mean():.3f}",
                          f"{(rcl[mr] == 12).mean():.3f}", int(mc.sum()), int(ml.sum()),
                          f"{mc.sum() / mr.sum():.3f}", f"{dens:.2f}"])
    _c = csa[cps == CROSS_PSID]
    sec1_note = (f"cross line {CROSS_PSID} scan angle after the CSF round trip: "
                 f"{_c.min():+.2f} to {_c.max():+.2f} deg, "
                 f"{100 * (_c != 0).mean():.1f}% nonzero")
    del rps, rcl, cps, ccl, crn, cnr, cx_, cy_, csa, _c

    # --------------------------------- 2. does the substituted reference change c?
    ns_meas = {}          # (tile, a, b) -> (c_pair, se) from the PIPELINE estimator
    sec2_rows = []
    if not A.skip_validation:
        rows = sec2_rows
        for tile in A.ns_tiles:
            meta = json.load(open(R.input(f"{tile}/corrections.json",
                                          role="pipeline corrections: grid geometry, per-swath "
                                               "alignment, geoid datum, drift curves")))
            res = float(meta["res_m"])
            b = meta["bounds"]
            nx = int(round((b[2] - b[0]) / res))
            clm, _, _, _, _ = cell_swath_ground(tile, R, A.min_cell_line, res)
            m = pair_rows(clm)
            m["blk"] = bs.block_ids(m.cell.to_numpy(), nx=nx, res=res, block_m=A.block_m)
            for (a, b_), g in m.groupby(["point_source_id_a", "point_source_id_b"]):
                X = np.c_[np.ones(len(g)), g.dtan.to_numpy()]
                beta, V, r2, n, G = ols_cluster(X, g.D.to_numpy(), g.blk.to_numpy())
                se = float(np.sqrt(V[1, 1]))
                ns_meas[(os.path.basename(tile), int(a), int(b_))] = (float(beta[1]), se)
                rows.append([os.path.basename(tile), f"{a}-{b_}", "pipeline d_mm_corr",
                             len(g), int(np.unique(g.blk).size),
                             f"{beta[1]:+.1f}", f"{se:.1f}", f"{r2:.4f}"])
            del clm, m
        vmeta = json.load(open(f"{A.valid_tile}/corrections.json"))
        vb = vmeta["bounds"]
        vres = float(vmeta["res_m"])
        vnx = int(round((vb[2] - vb[0]) / vres))
        vny = int(round((vb[3] - vb[1]) / vres))
        vclm, _, _, _, vn = gen1_cellline_ground(A.valid_las, vres, A.sn_smooth_cells,
                                                 A.ground_q, grid=(vb[0], vb[1], vnx, vny))
        vclm = vclm[vclm.n >= A.min_cell_line]
        vm = cross_pair_rows(vclm, vnx, vres, A.block_m)
        for (a, b_), g in vm.groupby(["point_source_id_a", "point_source_id_b"]):
            if len(g) < 200:
                continue
            f = fit_pair(g)
            rows.append([os.path.basename(A.valid_las), f"{a}-{b_}", "pooled-gen1 reference",
                         len(g), f["G"], f"{f['c_pair']:+.1f}", f"{f['c_pair_se']:.1f}",
                         f"{f['r2_pair']:.4f}"])
        del vclm, vm

    clm, nx, X0, Y0, nret = gen1_cellline_ground(A.cross_las, A.res, A.sn_smooth_cells,
                                                 A.ground_q)

    def run_variant(min_cl, controls):
        c = clm[clm.n >= min_cl]
        m = cross_pair_rows(c, nx, A.res, A.block_m)
        out = {}
        for (a, b_), g in m.groupby(["point_source_id_a", "point_source_id_b"]):
            if len(g) < 200:
                continue
            out[(int(a), int(b_))] = (fit_pair(g, controls=controls), g)
        return out

    fits = run_variant(A.min_cell_line, False)

    # ------------------------------------------------------------------- REPORT
    R.banner()
    print()
    print("## 1. Did CSF produce usable ground on the cross line?\n")
    R.table(["psid", "n_raw", "frac2", "frac12", "n_csf", "n_csf_last", "csf_rate", "dens"],
            sec1_rows)
    print("\n  " + sec1_note)

    if sec2_rows:
        print("\n## 2. The substituted reference surface, measured rather than argued\n")
        print("  Upper block: the pipeline quantity d_mm_corr (gen2 reference, registration\n"
              "  terms applied, reference cells only) through swath_across_track_test's own\n"
              "  code. Lower block: this script's estimator (pooled-gen1 reference, all cells)\n"
              "  on the same tile's CSF cache. Same model, same cluster blocks, same pairs.\n")
        R.table(["tile", "pair", "estimator", "cells", "blocks", "c_pair", "c_pair_se", "r2"],
                sec2_rows)

    print("\n## 3. The cross-line tile: geometry, then the three-parameter fit\n")
    print(f"  {nret:,} CSF last-return ground points on the cross tile, "
          f"{len(clm):,} (cell, line) ground estimates\n")
    rows = []
    for (a, b_), (fi, g) in sorted(fits.items()):
        rows.append([f"{a}-{b_}", "CROSS" if CROSS_PSID in (a, b_) else "N-S/N-S",
                     len(g), f"{len(g) * A.res ** 2 / 1e4:.1f}", fi["G"],
                     f"{np.corrcoef(g.ta, g.tb)[0, 1]:+.3f}", f"{g.stan.std():.4f}",
                     f"{g.dtan.std():.4f}", f"{fi['se_ratio']:.2f}",
                     f"{np.corrcoef(g.ta, g.cx)[0, 1]:+.2f}",
                     f"{np.corrcoef(g.ta, g.cy)[0, 1]:+.2f}",
                     f"{np.corrcoef(g.tb, g.cx)[0, 1]:+.2f}",
                     f"{np.corrcoef(g.tb, g.cy)[0, 1]:+.2f}"])
    R.table(["pair", "kind", "cells", "area_ha", "blocks", "corr_AB", "sd_sum", "sd_dif",
             "se_ratio", "r_tA_E", "r_tA_N", "r_tB_E", "r_tB_N"], rows)

    print("\n### 3a. two-parameter fit -- the quantity the N-S overlaps measure\n")
    rows = [[f"{a}-{b_}", "CROSS" if CROSS_PSID in (a, b_) else "N-S/N-S",
             f"{fi['c_pair']:+.1f}", f"{fi['c_pair_se']:.1f}",
             f"{fi['c_pair'] / fi['c_pair_se']:+.1f}", f"{fi['r2_pair']:.4f}"]
            for (a, b_), (fi, g) in sorted(fits.items())]
    R.table(["pair", "kind", "c_pair", "c_pair_se", "t", "r2"], rows)

    print("\n### 3b. three-parameter fit -- each line's OWN coefficient\n")
    rows = [[f"{a}-{b_}", "CROSS" if CROSS_PSID in (a, b_) else "N-S/N-S",
             f"{fi['k']:+.1f}", f"{fi['k_se']:.1f}",
             f"{fi['p']:+.1f}", f"{fi['p_se']:.1f}", f"{fi['q']:+.1f}", f"{fi['q_se']:.1f}",
             f"{fi['c_A']:+.1f}", f"{fi['c_A_se']:.1f}",
             f"{fi['c_B']:+.1f}", f"{fi['c_B_se']:.1f}", f"{fi['r2']:.4f}"]
            for (a, b_), (fi, g) in sorted(fits.items())]
    R.table(["pair", "kind", "k_mm", "k_se", "p", "p_se", "q", "q_se",
             "c_A", "c_A_se", "c_B", "c_B_se", "r2"], rows)

    # ----------------------------------------------- 4. per-line coefficients
    print("\n## 4. Each line's own coefficient, and the cross line's three estimates\n")
    own = {}
    rows = []
    cross_est = []
    for (a, b_), (fi, g) in sorted(fits.items()):
        if CROSS_PSID not in (a, b_):
            continue
        ns = a if b_ == CROSS_PSID else b_
        c_ns, se_ns = ((fi["c_A"], fi["c_A_se"]) if a == ns else (fi["c_B"], fi["c_B_se"]))
        c_cr, se_cr = ((fi["c_A"], fi["c_A_se"]) if a == CROSS_PSID else (fi["c_B"], fi["c_B_se"]))
        own[ns] = (c_ns, se_ns, f"cross pair {a}-{b_}")
        cross_est.append((c_cr, se_cr, f"{a}-{b_}"))
        rows.append([ns, f"{c_ns:+.1f}", f"{se_ns:.1f}", f"cross pair {a}-{b_}"])
    for c_cr, se_cr, lab in cross_est:
        rows.append([CROSS_PSID, f"{c_cr:+.1f}", f"{se_cr:.1f}", f"cross pair {lab}"])
    if cross_est:
        w = np.array([1.0 / s ** 2 for _, s, _ in cross_est])
        cs = np.array([c for c, _, _ in cross_est])
        cbar = float((w * cs).sum() / w.sum())
        cbar_se = float(np.sqrt(1.0 / w.sum()))
        chi2 = float((w * (cs - cbar) ** 2).sum())
        rows.append([CROSS_PSID, f"{cbar:+.1f}", f"{cbar_se:.1f}",
                     f"inverse-variance mean of the 3 (chi2={chi2:.2f}, dof=2)"])
    R.table(["line", "c_own", "c_own_se", "source"], rows)

    # --------------------------------- 5. the joint within-cell solve (all lines at once)
    print("\n## 5. The joint within-cell solve: every line's c_s on ONE shared population\n")
    print("  The pairwise fits of Sec 3-4 each use their own pair's own overlap, so any\n"
          "  disagreement between two of them could be the two REGIONS rather than the model.\n"
          "  This fit removes that escape: g(cell,line) = mu_cell + a_line + c_line*tan, every\n"
          "  (cell, line) estimate on the tile at once, cell fixed effects absorbing terrain,\n"
          "  cover and any spatial field exactly as differencing within a cell does.\n")
    J = joint_network(clm[clm.n >= A.min_cell_line], nx, A.res, A.block_m, ref_line=137)
    print(f"  {J['n']:,} (cell, line) rows on {J['cells']:,} multiply-covered cells, "
          f"{J['G']:,} blocks, within-cell r2 = {J['r2']:.4f}, "
          f"design condition number {J['cond']:.1f}\n")
    jc = {int(ln): (float(J["c"][i]), float(J["c_se"][i]))
          for i, ln in enumerate(J["lines"])}
    R.table(["line", "c_own", "c_own_se", "source"],
            [[ln, f"{c:+.1f}", f"{se:.1f}", "joint within-cell solve, all 5 lines"]
             for ln, (c, se) in sorted(jc.items())])

    # ---------------------------- 6. redundancy: does the per-line model reproduce the sums?
    print("\n## 6. Redundancy -- the first real test the model has ever had\n")
    print("  A per-line coefficient set is only meaningful if it reproduces the pair sums\n"
          "  measured directly on the N-S overlaps. The joint solve USES this tile's own N-S\n"
          "  overlaps, so its same-tile rows are in-sample residuals of an overdetermined fit\n"
          "  (5 coefficients against 9 informative combinations); the elba and elbaext rows\n"
          "  are a different region entirely and are out of sample.\n")
    rows = []
    Vc = J["V"][:J["L"], :J["L"]]
    li = {int(ln): i for i, ln in enumerate(J["lines"])}

    def _sum(c, V, idx, a, b_):
        v = np.zeros(len(c)); v[idx[a]] += 0.5; v[idx[b_]] += 0.5
        return float(v @ c), float(np.sqrt(v @ V @ v))

    for (a, b_), (fi, g) in sorted(fits.items()):
        if CROSS_PSID in (a, b_) or a not in li or b_ not in li:
            continue
        pr, pse = _sum(J["c"], Vc, li, a, b_)
        obs, ose = fi["c_pair"], fi["c_pair_se"]
        rows.append(["4342-28-64, joint solve", f"{a}-{b_}", f"{pr:+.1f}", f"{obs:+.1f}",
                     f"{pr - obs:+.1f}", f"{(pr - obs) / np.sqrt(pse ** 2 + ose ** 2):+.2f}",
                     "in-sample residual of the overdetermined joint fit, not a prediction"])
    for (tile, a, b_), (obs, ose) in sorted(ns_meas.items()):
        if a not in li or b_ not in li:
            continue
        pr, pse = _sum(J["c"], Vc, li, a, b_)
        rows.append([tile + ", joint solve", f"{a}-{b_}", f"{pr:+.1f}", f"{obs:+.1f}",
                     f"{pr - obs:+.1f}", f"{(pr - obs) / np.sqrt(pse ** 2 + ose ** 2):+.2f}",
                     "OUT OF SAMPLE, and 2.6-3.5 km SOUTH: also tests c_s along track"])
    R.table(["tile", "pair", "predicted", "observed", "resid", "resid_sig", "kind_check"], rows)

    # ------------------------------------ 7. propagate down the chain to 135, 134, 133
    print("\n## 7. Propagating the chain to lines 135, 134 and 133\n")
    if ns_meas:
        chain = {ln: c for ln, (c, se) in jc.items() if ln != CROSS_PSID}
        chain_se = {ln: se for ln, (c, se) in jc.items() if ln != CROSS_PSID}
        prop_rows = [[ln, f"{chain[ln]:+.1f}", f"{chain_se[ln]:.1f}",
                      "joint within-cell solve on the cross tile"] for ln in sorted(chain)]
        sums = {}
        for (tile, a, b_), (c, se) in ns_meas.items():
            sums.setdefault((a, b_), []).append((c, se))
        summ = {k: (float(sum(c / s ** 2 for c, s in v) / sum(1 / s ** 2 for c, s in v)),
                    float(np.sqrt(1.0 / sum(1 / s ** 2 for c, s in v))))
                for k, v in sums.items()}
        cur = min(chain)
        while (cur - 1, cur) in summ:
            P, Pse = summ[(cur - 1, cur)]
            chain[cur - 1] = 2 * P - chain[cur]
            chain_se[cur - 1] = float(np.sqrt((2 * Pse) ** 2 + chain_se[cur] ** 2))
            prop_rows.append([cur - 1, f"{chain[cur-1]:+.1f}", f"{chain_se[cur-1]:.1f}",
                              f"2*(elba+elbaext pair {cur-1}-{cur} sum {P:+.1f}) - c_{cur}"])
            cur -= 1
        R.table(["line", "c_own", "c_own_se", "source"], prop_rows)

    # ----------------------------------------------- 8. sensitivities
    print("\n## 8. Sensitivity of every reported coefficient\n")
    print("  Each variant re-runs BOTH the pairwise cross fits and the joint solve, and\n"
          "  reports how far each line's coefficient moves from the headline run.\n")
    base = {ln: own[ln][0] for ln in own}
    base_cross = float(np.mean([c for c, _, _ in cross_est])) if cross_est else float("nan")
    rows = []
    for lab, mcl, ctl in [("min_cell_line=3", 3, False), ("min_cell_line=5", 5, False),
                          ("quadratic position controls (E2,N2,EN)", A.min_cell_line, True)]:
        v = run_variant(mcl, ctl)
        vo, vc = {}, []
        for (a, b_), (fi, g) in v.items():
            if CROSS_PSID not in (a, b_):
                continue
            ns = a if b_ == CROSS_PSID else b_
            vo[ns] = fi["c_A"] if a == ns else fi["c_B"]
            vc.append(fi["c_A"] if a == CROSS_PSID else fi["c_B"])
        rows.append(["pairwise: " + lab]
                    + [f"{vo.get(ln, float('nan')) - base.get(ln, float('nan')):+.1f}"
                       for ln in (136, 137, 138)]
                    + [f"{np.mean(vc) - base_cross:+.1f}"])
    for lab, mcl in [("min_cell_line=3", 3), ("min_cell_line=5", 5)]:
        Jv = joint_network(clm[clm.n >= mcl], nx, A.res, A.block_m, ref_line=137)
        jv = {int(ln): float(Jv["c"][i]) for i, ln in enumerate(Jv["lines"])}
        rows.append(["joint: " + lab]
                    + [f"{jv.get(ln, float('nan')) - jc[ln][0]:+.1f}" for ln in (136, 137, 138)]
                    + [f"{jv.get(CROSS_PSID, float('nan')) - jc[CROSS_PSID][0]:+.1f}"])
    R.table(["variant", "d_c136", "d_c137", "d_c138", "d_c10010"], rows)

    print("\n## 9. What this run does NOT determine\n")
    print("  The GLOBAL CROSS-TRACK TILT of the gen1 mosaic is untouched by any of this.\n"
          "  SWATH_DEGENERACY_BREAKING Sec 1 shows the null direction of the overlap network\n"
          "  is e(x) = g*(x - x0), a tilt that is IDENTICAL for both lines at a shared ground\n"
          "  point and therefore cancels in EVERY between-line difference -- N-S, non-adjacent\n"
          "  or crossing -- to 0.000e+00 mm. Every coefficient above is a between-line\n"
          "  difference, and the joint solve's cell fixed effects remove exactly the same\n"
          "  thing. The tilt needs an absolute external reference (ground control); it is NOT\n"
          "  resolved here and the problem is NOT closed.\n")
    R.done(headline="per-line across-track coefficients from the cross line, and the network")


if __name__ == "__main__":
    main()
