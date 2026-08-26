#!/usr/bin/env python3
"""Is the gen1 per-swath vertical offset a CONSTANT, or a function of across-track position?

The pipeline models each gen1 flight line's disagreement with its neighbours as a constant
``(dx, dy, dz)`` (``coreg.align_swaths``). If the true error instead varies with **scan
angle** across the swath, that constant is an average over whatever across-track range the
tile happened to sample -- which would explain, in one mechanism, why the elba and elbaext
tiles disagree by 8.0/9.8/17.4 mm about the same swaths while sampling them at very
different mean scan angles.

Regressing the registered per-return offset ``d_mm_corr`` on ``scan_angle`` within a swath
CANNOT decide this: scan angle is collinear with across-track position, and there is a real
spatially varying residual field on these divides (38.0 mm real spatial sd, block medians
-167..+122 mm at 250 m -- ``analysis/MISSION_TIME_DRIFT.md``). A spatial field projects onto
scan angle and masquerades as beam geometry.

**The decisive test is the swath OVERLAP.** Where two flight lines cover the same cell they
observe the same ground at the same epoch, so terrain, land cover, vegetation, the gen2
reference and any spatial error field cancel exactly in the between-line difference, while
an across-track instrument error does not:

    D_cell = ground_A - ground_B = (k_A - k_B) + [ across-track error_A - across-track error_B ]

Three hypotheses, distinguished by the FORM of that bracket:

* **constant offset**       -> no scan-angle dependence at all; D is flat.
* **roll boresight** (odd)  -> ``c * (tan th_A - tan th_B)``, ONE ``c`` shared by every pair
  and both tiles, because a mounting roll is a sensor constant. Geometry: a pointing error
  ``delta`` at range ``R = h/cos th`` displaces the point by ``R*delta`` perpendicular to the
  beam, whose vertical component is ``R*delta*sin th = h*delta*tan th``. So ``c = h*delta``.
* **per-pair idiosyncratic** -> different ``c`` per pair; not an instrument term, and would
  point back at something spatial.

**A geometric warning this test must carry.** These lines are flown there-and-back, so in
every overlap the two swaths see the shared ground on the SAME body-fixed side: both scan
angles have the same sign, and their tangents sum to the constant line-spacing/height ratio
(``tan th_A + tan th_B = ~ +-S/h``). Two consequences, both measured and printed below:
(1) individual coefficients on ``tan th_A`` and ``tan th_B`` are NOT separately identifiable
-- only their antisymmetric combination is; (2) an EVEN (symmetric) across-track error
``c2 * tan^2 th`` reduces inside one pair to ``c2 * (tan th_A + tan th_B) * (tan th_A -
tan th_B)``, i.e. it is proportional to the odd predictor with a coefficient that FLIPS SIGN
between pairs of opposite side. Sign-alternation of the per-pair coefficient is therefore
the discriminator between an odd (roll) and an even (range/pointing) error, and it is
exactly the same quantity as "does one shared coefficient fit all pairs".

Reduction: the per-cell, per-swath ground estimate is the **median of ``d_mm_corr``**, which
is the pipeline's own ground estimator (``pipeline.difference_dem`` ``ground_estimator =
"slope_normal"``, ``ground_q = 0.50`` -- ``pipeline.py:603`` takes the ``ground_q`` quantile
of the slope-normal residual per cell) evaluated on one flight line at a time. No new
reduction is invented, and no minimum count is imposed by default.

Nothing in ``coreg.py`` or ``pipeline.py`` is modified; both are imported.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/swath_across_track_test.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lidar_diff_icp import binstats as bs
from lidar_diff_icp.refcells import reference_cells
from lidar_diff_icp.boresight import estimate_boresight
from lidar_diff_icp.registration import surface_gradients
from trust.provenance import Run

STRATA = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
          ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01)]

COLS = ["cell", "point_source_id", "scan_angle", "d_mm_corr", "dz_swath_mm", "in_grid",
        "slope", "aspect_deg"]


# ----------------------------------------------------------------- statistics
def ols_cluster(X, y, groups):
    """OLS with a cluster-robust (sandwich) covariance, clustered on ``groups``.

    Returns ``(beta, V, r2, n, n_clusters)``. The finite-sample factor is the usual
    ``G/(G-1) * (n-1)/(n-k)``. Clustering matters here because neighbouring overlap cells
    share woodlots, swath stripes and the residual spatial field; a return- or cell-counting
    SE is several times too small on these tiles (``binstats`` docstring).
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    g = np.asarray(groups)
    order = np.argsort(g, kind="stable")
    gs = g[order]
    Xs = X[order]
    es = e[order]
    starts = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    ends = np.r_[starts[1:], len(gs)]
    meat = np.zeros((k, k))
    for a, b in zip(starts, ends):
        s = Xs[a:b].T @ es[a:b]
        meat += np.outer(s, s)
    G = len(starts)
    corr = (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
    V = corr * (XtX_inv @ meat @ XtX_inv)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((e ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return beta, V, r2, n, G


def wald(beta, V, R):
    """Wald statistic and dof for ``R @ beta = 0`` under cluster-robust ``V``."""
    R = np.asarray(R, float)
    Rb = R @ beta
    M = R @ V @ R.T
    W = float(Rb @ np.linalg.pinv(M) @ Rb)
    return W, int(np.linalg.matrix_rank(M))


def chi2_sf(x, df):
    """Upper tail of chi^2 without scipy.stats (survival function)."""
    from math import exp, lgamma, log

    if x <= 0:
        return 1.0
    # regularised upper incomplete gamma Q(df/2, x/2) by continued fraction / series
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1.0:                                    # series for P, then Q = 1 - P
        term = 1.0 / a
        s = term
        ap = a
        for _ in range(2000):
            ap += 1.0
            term *= xx / ap
            s += term
            if abs(term) < abs(s) * 1e-14:
                break
        return max(0.0, 1.0 - s * exp(-xx + a * log(xx) - lgamma(a)))
    b, c, d, h = xx + 1.0 - a, 1e300, 1.0 / (xx + 1.0 - a), 1.0 / (xx + 1.0 - a)
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-14:
            break
    return h * exp(-xx + a * log(xx) - lgamma(a))


# ------------------------------------------------------------------- the data
def _aggregate(df, grad):
    """The pipeline's own ground estimator, one flight line at a time.

    ``pipeline.difference_dem`` runs ``ground_estimator="slope_normal"`` with
    ``ground_q=0.50``: per cell it takes the ``ground_q`` quantile of the slope-normal
    residual to the reference plane (``pipeline.py:603``). ``d_mm_corr`` IS that
    slope-normal residual per return, so the per-(cell, line) MEDIAN of it is that
    estimator evaluated on one line. Nothing new is invented here.

    ``gu``/``gv`` are the gen2 surface gradient from ``registration.surface_gradients`` --
    the same array the pipeline's own ``lateral_term`` projects a horizontal shift onto.
    They control for the one error that does NOT cancel in a between-line difference by
    being shared ground: a residual LATERAL misregistration BETWEEN the two lines, which
    lands on the difference as ``-(gx*ddx + gy*ddy)`` and would be mistaken for a scan-angle
    term if terrain gradient happened to trend across the overlap strip.
    """
    c = df.cell.to_numpy()
    df = df.assign(d_nosw=df.d_mm_corr.to_numpy() - df.dz_swath_mm.to_numpy(),
                   gu=grad[0].ravel()[c], gv=grad[1].ravel()[c])
    g = df.groupby(["cell", "point_source_id"], sort=True)
    return g.agg(med_corr=("d_mm_corr", "median"), med_nosw=("d_nosw", "median"),
                 sc=("scan_angle", "mean"), gu=("gu", "mean"), gv=("gv", "mean"),
                 n=("d_mm_corr", "size")).reset_index()


def cell_swath_ground(tile, R, min_cell_line, res):
    """``(clm_reference, clm_all_in_grid, per_return_reference, n_cellline, n_returns)``."""
    z = np.load(R.input(f"{tile}/z_after.npy",
                        role="gen2 gridded ground (the DoD reference surface), m"))
    gx, gy, _ = surface_gradients(z, res)
    grad = (gx, gy)
    df = pd.read_parquet(R.input(f"{tile}/beam_offset_table.parquet",
                                 role="per-return gen1 CSF ground offsets to the gen2 surface, "
                                      "slope-normal mm, with the four registration terms"),
                         columns=COLS)
    df = df[df.in_grid.to_numpy()]
    clm_all = _aggregate(df, grad)
    m, rep = reference_cells(tile)
    R.cuts(f"{os.path.basename(tile)} reference_cells", rep)
    R.mask(f"{os.path.basename(tile)} reference cells", m, of=m.size,
           defn="refcells.reference_cells() at its repo defaults: divide cells, |curv|<=0.015, "
                "slope<12 deg, no building returns, not clear-cut, not blufftop margin, |DoD|<=500 mm")
    sel = m[df.cell.to_numpy()]
    df = df[sel]
    clm = _aggregate(df, grad)
    ret = df[["cell", "point_source_id", "scan_angle", "d_mm_corr"]].copy()
    del df
    return clm[clm.n >= min_cell_line], clm_all, ret, len(clm), len(ret)


def pair_rows(clm):
    """One row per (cell, unordered line pair) -- the same construction as
    ``boresight._pair_rows``, on the median rather than the mean."""
    m = clm.merge(clm, on="cell", suffixes=("_a", "_b"))
    m = m[m.point_source_id_a < m.point_source_id_b].copy()
    m["D"] = m.med_corr_a - m.med_corr_b
    m["D_nosw"] = m.med_nosw_a - m.med_nosw_b
    m["ta"] = np.tan(np.radians(m.sc_a.to_numpy()))
    m["tb"] = np.tan(np.radians(m.sc_b.to_numpy()))
    m["dtan"] = m.ta - m.tb
    m["stan"] = m.ta + m.tb
    m["dlin"] = m.sc_a - m.sc_b
    m["dtan2"] = m.ta ** 2 - m.tb ** 2
    m["dgu"] = m.gu_a - m.gu_b
    m["dgv"] = m.gv_a - m.gv_b
    m["gu"] = 0.5 * (m.gu_a + m.gu_b)
    m["gv"] = 0.5 * (m.gv_a + m.gv_b)
    return m


# ------------------------------------------------------------------------ run
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", nargs="+",
                    default=["data/derived/elba_fulldensity", "data/derived/elbaext"])
    ap.add_argument("--min-cell-line", type=int, default=1,
                    help="minimum returns for a (cell, flight-line) ground estimate. 1 is the "
                         "definitional floor (a median needs one value). The repo's "
                         "boresight.estimate_boresight uses 3; --min-cell-line 3 reproduces it.")
    ap.add_argument("--block-m", type=float, default=50.0,
                    help="spatial block size (m) for the cluster-robust standard errors; the "
                         "repo default used by offset_vs_angle.py and cover_offset_regression.py")
    A = ap.parse_args()

    R = Run("is the gen1 per-swath vertical offset a constant, or a function of across-track "
            "(scan-angle) position that the constant averages over?")
    R.param("min_cell_line", A.min_cell_line, src="MINE" if A.min_cell_line == 1 else "repo",
            why="1 is the definitional floor -- a median needs one value, so this imposes NO "
                "cut at all. The repo's boresight.estimate_boresight uses 3; that value is "
                "run as a sensitivity in section 9 and both answers are reported.")
    R.param("block_m", A.block_m, src="repo")
    R.param("ground_q", 0.50, src="repo")
    R.param("reference_cells defaults", "curv_max=0.015, slope_max=12, gross_change_mm=500, "
                                        "clearcut_drop=0.30", src="repo")

    for c, d in [
        ("tile", "derived-product directory the returns and cells come from"),
        ("pair", "unordered flight-line pair A-B (point_source_id), A < B"),
        ("cells", "overlap cells: reference cells covered by BOTH lines of the pair"),
        ("ret_A", "gen1 returns of line A in those cells"),
        ("ret_B", "gen1 returns of line B in those cells"),
        ("blocks", f"independent {A.block_m:g} m spatial blocks spanned by the overlap cells"),
        ("scan_A", "mean scan angle of line A over the overlap cells, deg (signed, body frame)"),
        ("scan_B", "mean scan angle of line B over the overlap cells, deg (signed, body frame)"),
        ("sum_tan", "mean of tan(scan_A)+tan(scan_B) over the overlap; ~ +-S/h, the line "
                    "spacing over flying height, and the sign of the side the pair shares"),
        ("dtan_lo", "minimum of tan(scan_A)-tan(scan_B) over the overlap cells"),
        ("dtan_hi", "maximum of tan(scan_A)-tan(scan_B) over the overlap cells"),
        ("dtan_sd", "standard deviation of tan(scan_A)-tan(scan_B): the leverage available"),
        ("k_mm", "fitted intercept: between-line offset at tan(scan_A)=tan(scan_B), mm, "
                 "positive = line A reads above line B"),
        ("k_se", "cluster-robust standard error of k_mm, mm"),
        ("c_mm", "fitted coefficient on tan(scan_A)-tan(scan_B), mm per unit tangent"),
        ("c_se", "cluster-robust standard error of c_mm, mm per unit tangent"),
        ("c_mmdeg", "c_mm converted to mm per degree at small angle (c_mm * pi/180)"),
        ("t", "c_mm / c_se, cluster-robust"),
        ("r2", "coefficient of determination of the pair's own fit"),
        ("medD", "median between-line difference over the overlap cells, mm"),
        ("nmadD", "normalised MAD of the between-line difference, mm"),
        ("model", "which functional form was fitted"),
        ("c_shared", "single coefficient shared by all pairs, mm per unit tangent"),
        ("W", "cluster-robust Wald statistic for the stated restriction"),
        ("df", "degrees of freedom of the Wald statistic"),
        ("p", "upper-tail probability of W under chi-squared with df"),
        ("rss", "residual sum of squares of the pooled fit, mm^2"),
        ("swath", "flight line (point_source_id)"),
        ("returns", "gen1 returns of that swath on the reference cells"),
        ("slope_mmdeg", "OLS slope of the response on scan_angle within the swath, mm/deg"),
        ("slope_se", "cluster-robust standard error of that slope, mm/deg"),
        ("scan_sd", "standard deviation of scan_angle within the swath, deg"),
        ("scan_mean", "mean scan_angle within the swath, deg"),
        ("corr_ks", "correlation between the fitted constant and slope of that swath's "
                    "regression; |corr| -> 1 means the two are not separable"),
        ("resid_mmdeg", "the same slope after subtracting the shared roll term, mm/deg"),
        ("dz_elba", "elba per-swath dz re-referenced to swath 135, mm"),
        ("dz_ext", "elbaext per-swath dz re-referenced to swath 135, mm"),
        ("disagree", "dz_elba - dz_ext, mm: the tile-to-tile disagreement to be explained"),
        ("mean_dtan_elba", "mean tan(scan_A)-tan(scan_B) sampled by elba in that pair's overlap"),
        ("mean_dtan_ext", "mean tan(scan_A)-tan(scan_B) sampled by elbaext in the same overlap"),
        ("pred_mm", "disagreement predicted by the shared roll: c_shared * (mean_dtan_elba - "
                    "mean_dtan_ext), mm"),
        ("bin_lo", "lower edge of the tan(scan_A)-tan(scan_B) bin"),
        ("bin_hi", "upper edge of the bin"),
        ("bin_n", "overlap cells in the bin"),
        ("bin_med", "median between-line difference in the bin, mm"),
        ("bin_se", "cluster-robust standard error of that median, mm"),
        ("b_mmdeg", "repo boresight.estimate_boresight pooled roll, mm/deg"),
        ("b_se_mmdeg", "its bootstrap SE (cell resampling only -- optimistic), mm/deg"),
        ("b_pairsd", "std of its per-pair roll estimates -- the honest uncertainty, mm/deg"),
    ]:
        R.column(c, d)

    tiles = {}
    for tile in A.tiles:
        meta = json.load(open(R.input(f"{tile}/corrections.json",
                                      role="pipeline corrections: grid geometry, per-swath "
                                           "internal alignment, geoid datum, drift curves")))
        res = float(meta["res_m"])
        clm, clm_all, ret, n_all, n_ret = cell_swath_ground(tile, R, A.min_cell_line, res)
        b = meta["bounds"]
        nx = int(round((b[2] - b[0]) / res))
        cover = np.load(R.input(f"{tile}/canopy_cover_pfs.npy",
                                role="PyForestScan canopy cover from the gen2 cloud, "
                                     "plant-area based, min_height 2 m")).ravel()
        m = pair_rows(clm)
        m["blk"] = bs.block_ids(m.cell.to_numpy(), nx=nx, res=res, block_m=A.block_m)
        m["cover"] = cover[m.cell.to_numpy()]
        m["cx"] = (m.cell.to_numpy() % nx) * res      # across-track here: these lines fly N-S
        m["cy"] = (m.cell.to_numpy() // nx) * res     # along-track
        m_all = pair_rows(clm_all[clm_all.n >= A.min_cell_line])
        m_all["blk"] = bs.block_ids(m_all.cell.to_numpy(), nx=nx, res=res, block_m=A.block_m)
        tiles[tile] = dict(meta=meta, clm=clm, m=m, m_all=m_all, ret=ret, res=res, nx=nx,
                           n_cellline=n_all, n_ret=n_ret)

    R.banner()
    print()
    np.save  # keep numpy referenced for linters

    # ---------------------------------------------------------------- 1. geometry
    print("## 1. The overlap geometry, and what leverage it gives\n")
    rows = []
    for tile, T in tiles.items():
        for (a, b_), g in T["m"].groupby(["point_source_id_a", "point_source_id_b"]):
            rows.append([os.path.basename(tile), f"{a}-{b_}", len(g),
                         int(g.n_a.sum()), int(g.n_b.sum()),
                         int(np.unique(g.blk).size),
                         f"{g.sc_a.mean():+.2f}", f"{g.sc_b.mean():+.2f}",
                         f"{g.stan.mean():+.3f}", f"{g.stan.std():.4f}",
                         f"{g.dtan.min():+.3f}", f"{g.dtan.max():+.3f}", f"{g.dtan.std():.4f}"])
    R.column("sum_tan_sd", "standard deviation of tan(scan_A)+tan(scan_B) within the pair. "
                           "Near zero means the two lines see the shared ground on the same "
                           "body-fixed side with their tangents summing to a fixed S/h -- the "
                           "there-and-back flight pattern, and the reason only the "
                           "antisymmetric combination of the two angles is identifiable.")
    R.table(["tile", "pair", "cells", "ret_A", "ret_B", "blocks", "scan_A", "scan_B",
             "sum_tan", "sum_tan_sd", "dtan_lo", "dtan_hi", "dtan_sd"], rows)

    # ---------------------------------------------------------- 2. per-pair fits
    print("\n## 2. Per-pair fit  D = k + c * (tan scan_A - tan scan_B)\n")
    rows = []
    percc = {}
    for tile, T in tiles.items():
        for (a, b_), g in T["m"].groupby(["point_source_id_a", "point_source_id_b"]):
            X = np.c_[np.ones(len(g)), g.dtan.to_numpy()]
            beta, V, r2, n, G = ols_cluster(X, g.D.to_numpy(), g.blk.to_numpy())
            se = np.sqrt(np.diag(V))
            Xg = np.c_[X, g.dgu.to_numpy(), g.dgv.to_numpy(), g.gu.to_numpy(), g.gv.to_numpy()]
            bg, Vg, _, _, _ = ols_cluster(Xg, g.D.to_numpy(), g.blk.to_numpy())
            yy = g.cy.to_numpy()
            yy = (yy - yy.mean()) / 1000.0
            Xa = np.c_[X, yy, yy ** 2]
            ba, _, _, _, _ = ols_cluster(Xa, g.D.to_numpy(), g.blk.to_numpy())
            ry = np.corrcoef(g.dtan.to_numpy(), g.cy.to_numpy())[0, 1]
            rx = np.corrcoef(g.dtan.to_numpy(), g.cx.to_numpy())[0, 1]
            percc[(tile, a, b_)] = (beta[1], se[1], float(g.stan.mean()), ba[1])
            rows.append([os.path.basename(tile), f"{a}-{b_}",
                         f"{beta[0]:+.2f}", f"{se[0]:.2f}",
                         f"{beta[1]:+.1f}", f"{se[1]:.1f}", f"{beta[1] * np.pi / 180:+.3f}",
                         f"{beta[1] / se[1]:+.1f}", f"{r2:.4f}",
                         f"{np.median(g.D):+.2f}", f"{bs.nmad(g.D.to_numpy()):.1f}",
                         f"{bg[1]:+.1f}", f"{ba[1]:+.1f}", f"{rx:+.2f}", f"{ry:+.2f}"])
    R.column("c_ctrl", "the same coefficient with the terrain-gradient controls added "
                       "(between-line and mean gu, gv): what survives a residual LATERAL "
                       "misregistration between the two lines, mm per unit tangent")
    R.column("c_along", "the same coefficient with an along-track control (cell northing, "
                        "linear + quadratic) added: these lines fly N-S, so across-track is "
                        "easting and any residual along-track error is absorbed here")
    R.column("r_dtan_x", "correlation of tan(scan_A)-tan(scan_B) with cell EASTING (across-track)")
    R.column("r_dtan_y", "correlation of tan(scan_A)-tan(scan_B) with cell NORTHING (along-track)")
    R.table(["tile", "pair", "k_mm", "k_se", "c_mm", "c_se", "c_mmdeg", "t", "r2",
             "medD", "nmadD", "c_ctrl", "c_along", "r_dtan_x", "r_dtan_y"], rows)

    # ------------------- 2b. is the across-track term a CANOPY effect or a per-line one?
    print("\n## 2b. The same coefficient inside canopy-cover strata\n")
    print("  A per-line navigation/pointing error is a property of the LINE: its coefficient\n"
          "  should be the same in open and in vegetated cells of the same overlap. A canopy\n"
          "  effect is a property of the GROUND: path length through the canopy grows as\n"
          "  1/cos(scan), so the ground-return selection differs with scan angle -- and that\n"
          "  coefficient should be ~0 in open cells and grow with cover, with the SAME shape\n"
          "  in every pair. The two hypotheses are told apart here, not by argument.\n")
    R.column("stratum", "canopy-cover band (canopy_cover_pfs), the repo's four strata")
    rows = []
    for tile, T in tiles.items():
        for (a, b_), g in T["m"].groupby(["point_source_id_a", "point_source_id_b"]):
            for nm, lo, hi in STRATA:
                gg = g[(g.cover > lo) & (g.cover <= hi)]
                if not len(gg):
                    continue
                Xs = np.c_[np.ones(len(gg)), gg.dtan.to_numpy()]
                bs_, Vs, _, _, _ = ols_cluster(Xs, gg.D.to_numpy(), gg.blk.to_numpy())
                ses = float(np.sqrt(Vs[1, 1]))
                rows.append([os.path.basename(tile), f"{a}-{b_}", nm, len(gg),
                             f"{bs_[1]:+.1f}", f"{ses:.1f}", f"{bs_[1] / ses:+.1f}"])
    R.table(["tile", "pair", "stratum", "cells", "c_mm", "c_se", "t"], rows)

    # -------------------------------------------- 3. does ONE coefficient fit all?
    print("\n## 3. Does ONE shared coefficient fit every pair, on both tiles?\n")
    keys = sorted(percc)
    allm = []
    for tile, T in tiles.items():
        gg = T["m"].copy()
        gg["pairkey"] = [f"{os.path.basename(tile)}:{a}-{b_}" for a, b_ in
                         zip(gg.point_source_id_a, gg.point_source_id_b)]
        gg["blkkey"] = [f"{os.path.basename(tile)}:{v}" for v in gg.blk]
        allm.append(gg)
    P = pd.concat(allm, ignore_index=True)
    pk = pd.Categorical(P.pairkey)
    Dm = pd.get_dummies(pk).to_numpy(float)          # pair fixed effects (one per pair)
    npair = Dm.shape[1]
    y = P.D.to_numpy()
    blk = P.blkkey.to_numpy()

    rows = []
    fits = {}
    ctrl = np.c_[P.dgu.to_numpy(), P.dgv.to_numpy(), P.gu.to_numpy(), P.gv.to_numpy()]
    for name, pred, extra in [("constant only (no scan term)", None, None),
                              ("shared c * (tanA - tanB)", P.dtan.to_numpy(), None),
                              ("shared c * (scanA - scanB) deg", P.dlin.to_numpy(), None),
                              ("shared c2 * (tan^2A - tan^2B)", P.dtan2.to_numpy(), None),
                              ("shared c, + terrain-gradient controls", P.dtan.to_numpy(), ctrl)]:
        X = Dm if pred is None else np.c_[Dm, pred]
        if extra is not None:
            X = np.c_[Dm, extra, pred]
        beta, V, r2, n, G = ols_cluster(X, y, blk)
        e = y - X @ beta
        rss = float((e ** 2).sum())
        fits[name] = (beta, V, r2, rss, X)
        if pred is None:
            rows.append([name, "-", "-", "-", f"{r2:.4f}", f"{rss:.4g}"])
        else:
            se = float(np.sqrt(V[-1, -1]))
            rows.append([name, f"{beta[-1]:+.1f}", f"{se:.1f}", f"{beta[-1] / se:+.1f}",
                         f"{r2:.4f}", f"{rss:.4g}"])
    # pair-specific along-track controls: the two 137-138 pairs are the only ones whose
    # dtan correlates appreciably with northing, so the shared fit is re-run with a linear
    # and quadratic northing term PER PAIR, absorbing any residual along-track per-line error.
    ycen = P.groupby("pairkey")["cy"].transform(lambda v: (v - v.mean()) / 1000.0).to_numpy()
    Y1 = Dm * ycen[:, None]
    Y2 = Dm * (ycen ** 2)[:, None]
    Xy = np.c_[Dm, Y1, Y2, P.dtan.to_numpy()]
    by, Vy, r2y, _, _ = ols_cluster(Xy, y, blk)
    sey = float(np.sqrt(Vy[-1, -1]))
    fits["shared c, + per-pair along-track controls"] = (by, Vy, r2y, float(
        ((y - Xy @ by) ** 2).sum()), Xy)
    rows.append(["shared c, + per-pair along-track controls", f"{by[-1]:+.1f}", f"{sey:.1f}",
                 f"{by[-1] / sey:+.1f}", f"{r2y:.4f}", f"{float(((y - Xy @ by) ** 2).sum()):.4g}"])

    # the same shared-c fit on EVERY in-grid overlap cell, not only the reference cells:
    # the between-line difference is terrain-free by construction, so the divide/slope cuts
    # are not needed for it, and dropping them multiplies the sample ~6x.
    allp = []
    for tile, T in tiles.items():
        gg = T["m_all"].copy()
        gg["pairkey"] = [f"{os.path.basename(tile)}:{a}-{b_}" for a, b_ in
                         zip(gg.point_source_id_a, gg.point_source_id_b)]
        gg["blkkey"] = [f"{os.path.basename(tile)}:{v}" for v in gg.blk]
        allp.append(gg)
    PA = pd.concat(allp, ignore_index=True)
    DmA = pd.get_dummies(pd.Categorical(PA.pairkey)).to_numpy(float)
    bA, VA, r2A, nA, GA = ols_cluster(np.c_[DmA, PA.dtan.to_numpy()], PA.D.to_numpy(),
                                      PA.blkkey.to_numpy())
    seA = float(np.sqrt(VA[-1, -1]))
    rows.append([f"shared c, ALL in-grid overlap cells ({len(PA):,})", f"{bA[-1]:+.1f}",
                 f"{seA:.1f}", f"{bA[-1] / seA:+.1f}", f"{r2A:.4f}",
                 f"{float(((PA.D.to_numpy() - np.c_[DmA, PA.dtan.to_numpy()] @ bA) ** 2).sum()):.4g}"])
    R.column("c_shared_se", "cluster-robust standard error of c_shared")
    R.table(["model", "c_shared", "c_shared_se", "t", "r2", "rss"], rows)

    # per-pair coefficients in one model, then test they are all equal
    Inter = Dm * P.dtan.to_numpy()[:, None]
    X = np.c_[Dm, Inter]
    beta, V, r2, n, G = ols_cluster(X, y, blk)
    Rm = np.zeros((npair - 1, X.shape[1]))
    for i in range(npair - 1):
        Rm[i, npair + i] = 1.0
        Rm[i, npair + i + 1] = -1.0
    W, dfree = wald(beta, V, Rm)
    print(f"\n  Homogeneity of the {npair} per-pair coefficients (cluster-robust Wald, "
          f"H0: all c_pair equal): W = {W:.1f}, df = {dfree}, p = {chi2_sf(W, dfree):.3g}")
    Xi = np.c_[Dm, Y1, Y2, Inter]
    bi, Vi, r2i, _, _ = ols_cluster(Xi, y, blk)
    Ri = np.zeros((npair - 1, Xi.shape[1]))
    for i in range(npair - 1):
        Ri[i, 3 * npair + i] = 1.0
        Ri[i, 3 * npair + i + 1] = -1.0
    Wi, dfi = wald(bi, Vi, Ri)
    print(f"  Same test WITH the per-pair along-track controls: W = {Wi:.1f}, df = {dfi}, "
          f"p = {chi2_sf(Wi, dfi):.3g}")
    print(f"  Per-pair model r2 = {r2:.4f} against shared-c r2 = "
          f"{fits['shared c * (tanA - tanB)'][2]:.4f} "
          f"({npair - 1} extra parameters).")

    # cross-tile repeatability: three pairs are measured by BOTH tiles, on different extents
    print("\n  Cross-tile repeatability of the per-pair coefficient. Three pairs are measured\n"
          "  independently by both tiles, over different extents and different ground:\n")
    R.column("c_elba", "per-pair coefficient from the elba tile, mm per unit tangent")
    R.column("c_ext", "per-pair coefficient from the elbaext tile, mm per unit tangent")
    R.column("diff_se", "difference c_elba - c_ext divided by the standard error of that "
                        "difference (the two tiles share ground, so this is conservative "
                        "only if they were independent -- see the text)")
    rows = []
    te0, tx0 = A.tiles[0], A.tiles[1]
    for (t, a, b_) in keys:
        if t != te0:
            continue
        k2 = (tx0, a, b_)
        if k2 not in percc:
            continue
        c1, s1, _, _ = percc[(t, a, b_)]
        c2, s2, _, _ = percc[k2]
        sd = np.hypot(s1, s2)
        rows.append([f"{a}-{b_}", f"{c1:+.1f}", f"{s1:.1f}", f"{c2:+.1f}", f"{s2:.1f}",
                     f"{(c1 - c2) / sd:+.2f}"])
    R.column("c_ext_se", "cluster-robust standard error of c_ext, mm per unit tangent")
    R.table(["pair", "c_elba", "c_se", "c_ext", "c_ext_se", "diff_se"], rows)

    # sign alternation: odd (roll) vs even (symmetric) error
    print("\n  Per-pair coefficient against the side the pair shares (sum_tan). An ODD (roll)\n"
          "  error gives the SAME c on every pair; an EVEN error gives c proportional to\n"
          "  sum_tan, i.e. sign-flipping between pairs of opposite side:\n")
    R.column("c_over_sumtan", "c_mm divided by the pair's mean sum_tan; constant across pairs "
                              "if the error is EVEN in scan angle")
    rows = [[os.path.basename(t), f"{a}-{b_}", f"{v[2]:+.3f}", f"{v[0]:+.1f}", f"{v[1]:.1f}",
             f"{v[0] / v[2]:+.1f}"] for (t, a, b_), v in
            ((k, percc[k]) for k in keys)]
    R.table(["tile", "pair", "sum_tan", "c_mm", "c_se", "c_over_sumtan"], rows)

    # -------------------------------------------------- 4. tan vs linear in angle
    print("\n## 4. tan(scan) against linear-in-scan\n")
    rt = np.corrcoef(P.dtan.to_numpy(), np.radians(P.dlin.to_numpy()))[0, 1]
    tmax = float(np.abs(P.dtan).max())
    print(f"  corr[ tan(scan_A)-tan(scan_B) , radians(scan_A-scan_B) ] = {rt:.6f}")
    print(f"  Over the observed |scan| <= 17 deg the two predictors differ by at most "
          f"{100 * abs(np.tan(np.radians(17.0)) - np.radians(17.0)) / np.tan(np.radians(17.0)):.1f}% "
          f"in shape; max |dtan| observed = {tmax:.3f}.")
    print(f"  rss: tan {fits['shared c * (tanA - tanB)'][3]:.6g}  vs  "
          f"linear {fits['shared c * (scanA - scanB) deg'][3]:.6g}  "
          f"(difference {100 * (fits['shared c * (scanA - scanB) deg'][3] / fits['shared c * (tanA - tanB)'][3] - 1):+.4f}%)")

    # ---------------------------------------------- 5. the data, binned (no fit)
    print("\n## 5. The data behind the fit: between-line difference binned by dtan\n")
    rows = []
    for tile, T in tiles.items():
        for (a, b_), g in T["m"].groupby(["point_source_id_a", "point_source_id_b"]):
            e = bs.quantile_edges(g.dtan.to_numpy(), 6)
            st = bs.binned_stats(g.dtan.to_numpy(), g.D.to_numpy(), e, block=g.blk.to_numpy())
            for i in range(len(st.x)):
                rows.append([os.path.basename(tile), f"{a}-{b_}",
                             f"{st.lo[i]:+.3f}", f"{st.hi[i]:+.3f}", int(st.n[i]),
                             f"{st.y[i]:+.1f}", f"{st.se[i]:.1f}"])
    R.table(["tile", "pair", "bin_lo", "bin_hi", "bin_n", "bin_med", "bin_se"], rows)

    # --------------------------------------- 6. the elba / elbaext disagreement
    print("\n## 6. What an across-track term would do to the elba / elbaext disagreement\n")
    c_shared = fits["shared c * (tanA - tanB)"][0][-1]
    c_shared_se = float(np.sqrt(fits["shared c * (tanA - tanB)"][1][-1, -1]))
    al = {t: {int(k): v for k, v in tiles[t]["meta"]["per_swath_internal_alignment_dxdydz_m"].items()}
          for t in tiles}
    te, tx = A.tiles[0], A.tiles[1]
    ze = {s_: 1000 * v[2] for s_, v in al[te].items()}
    zx = {s_: 1000 * v[2] for s_, v in al[tx].items()}
    ze = {s_: v - ze[135] for s_, v in ze.items()}
    zx = {s_: v - zx[135] for s_, v in zx.items()}
    print("  align_swaths ties each line to its neighbour over their OVERLAP. If the\n"
          "  between-line difference varies across that overlap, the tie is that difference\n"
          "  AVERAGED over whatever part of the overlap the tile covers -- so two tiles of\n"
          "  different extent get ties differing by c * (mean dtan_1 - mean dtan_2), and the\n"
          "  disagreement in a swath's dz ACCUMULATES along the chain from swath 135.\n"
          "  Means are over ALL in-grid overlap cells (coregister_swaths uses the whole\n"
          "  overlap, not the reference cells; it grids at 2 m on the vendor terrain classes\n"
          "  rather than 5 m on CSF ground, so this is a close proxy, not identical).\n")
    R.column("c_pair", "the pair's own coefficient, averaged over the tiles that measure it, "
                       "mm per unit tangent")
    R.column("step_mm", "c_pair * (mean_dtan_elba - mean_dtan_ext) for this link of the chain, mm")
    R.column("cum_mm", "running sum of step_mm from swath 135: the predicted disagreement, mm")
    R.column("cum_shared", "the same running sum using the single shared coefficient, mm")
    rows = []
    cum = cum_shared = cum_alt = 0.0
    for s_ in (136, 137, 138):
        a_, b_ = s_ - 1, s_
        me = tiles[te]["m_all"]
        me = me[(me.point_source_id_a == a_) & (me.point_source_id_b == b_)]
        mx = tiles[tx]["m_all"]
        mx = mx[(mx.point_source_id_a == a_) & (mx.point_source_id_b == b_)]
        cs = [percc[(t, a_, b_)][0] for t in (te, tx) if (t, a_, b_) in percc]
        ca = [percc[(t, a_, b_)][3] for t in (te, tx) if (t, a_, b_) in percc]
        cp, cpa = float(np.mean(cs)), float(np.mean(ca))
        dme, dmx = float(me.dtan.mean()), float(mx.dtan.mean())
        step = cp * (dme - dmx)
        cum += step
        cum_alt += cpa * (dme - dmx)
        cum_shared += c_shared * (dme - dmx)
        rows.append([str(s_), f"{ze[s_]:+.1f}", f"{zx[s_]:+.1f}", f"{ze[s_] - zx[s_]:+.1f}",
                     f"{dme:+.4f}", f"{dmx:+.4f}", f"{cp:+.1f}", f"{step:+.2f}",
                     f"{cum:+.2f}", f"{cum_alt:+.2f}", f"{cum_shared:+.2f}"])
    R.column("cum_along", "the same running sum using the along-track-controlled per-pair "
                          "coefficients, mm")
    R.table(["swath", "dz_elba", "dz_ext", "disagree", "mean_dtan_elba", "mean_dtan_ext",
             "c_pair", "step_mm", "cum_mm", "cum_along", "cum_shared"], rows)
    print(f"\n  c_shared = {c_shared:+.1f} +- {c_shared_se:.1f} mm per unit tangent "
          f"({c_shared * np.pi / 180:+.3f} mm/deg at small angle).")

    # ------------------------------- 7. the within-swath residual slopes, before/after
    print("\n## 7. Within-swath regression of d_mm_corr on scan_angle, per RETURN\n")
    print("  Per return, on the same reference cells, so the `before` column reproduces the\n"
          "  tabulated residual slopes this test was set to explain. `after` subtracts the\n"
          "  single shared roll c_shared * tan(scan) from every return first.\n")
    R.column("fit0_se", "standard error of the fitted offset extrapolated to scan_angle = 0, "
                        "mm -- the quantity a per-swath CONSTANT claims to be")
    R.column("fitm_se", "standard error of the fitted offset at the swath's own mean scan "
                        "angle, mm -- the quantity the data actually pin down")
    rows = []
    for tile, T in tiles.items():
        ret = T["ret"]
        for s_, g in ret.groupby("point_source_id"):
            blkg = bs.block_ids(g.cell.to_numpy(), nx=T["nx"], res=T["res"], block_m=A.block_m)
            sc = g.scan_angle.to_numpy(float)
            X = np.c_[np.ones(len(g)), sc]
            yv = g.d_mm_corr.to_numpy(float)
            beta, V, r2, n, G = ols_cluster(X, yv, blkg)
            se = np.sqrt(np.diag(V))
            rho = V[0, 1] / np.sqrt(V[0, 0] * V[1, 1])
            se0 = float(np.sqrt(V[0, 0]))
            xm = np.array([1.0, sc.mean()])
            sem = float(np.sqrt(xm @ V @ xm))
            b2, V2, _, _, _ = ols_cluster(X, yv - c_shared * np.tan(np.radians(sc)), blkg)
            rows.append([os.path.basename(tile), str(s_), int(len(g)),
                         f"{sc.mean():+.2f}", f"{sc.std():.2f}",
                         f"{beta[1]:+.3f}", f"{se[1]:.3f}", f"{beta[1] / se[1]:+.1f}",
                         f"{rho:+.3f}", f"{se0:.1f}", f"{sem:.1f}", f"{b2[1]:+.3f}"])
    R.table(["tile", "swath", "returns", "scan_mean", "scan_sd", "slope_mmdeg", "slope_se",
             "t", "corr_ks", "fit0_se", "fitm_se", "resid_mmdeg"], rows)

    # --------------------------- 7b. per-swath coefficients: what the overlaps can and cannot say
    print("\n## 7b. Per-SWATH across-track coefficients: identifiable only up to an alternating sign\n")
    print("  Within one pair tan th_A + tan th_B is fixed at the line-spacing/height ratio S/h,\n"
          "  so a per-swath model  err_s = c_s * tan th_s  collapses to a pair coefficient\n"
          "  c_pair = (c_A + c_B)/2. A chain of n lines gives n-1 such sums, so the c_s are\n"
          "  determined only up to adding an ALTERNATING vector (+v, -v, +v, ...): overlaps\n"
          "  alone cannot separate them, and only a loop (a non-adjacent overlap) or an\n"
          "  external reference can. Reported: the pair sums, against the same quantity read\n"
          "  from the gen2-referenced within-swath slopes of section 7, which ARE contaminated\n"
          "  by the spatial residual field this test was designed to remove.\n")
    R.column("c_pair_ovl", "(c_A + c_B)/2 from the OVERLAP, mm/deg at small angle")
    R.column("c_pair_gen2", "the same quantity from section 7's gen2-referenced per-swath "
                            "slopes, (slope_A + slope_B)/2, mm/deg")
    R.column("gap", "c_pair_gen2 - c_pair_ovl, mm/deg: what the gen2-referenced regression "
                    "attributes to scan angle that the overlap does not")
    sl = {}
    for tile, T in tiles.items():
        for s_, g in T["ret"].groupby("point_source_id"):
            blkg = bs.block_ids(g.cell.to_numpy(), nx=T["nx"], res=T["res"], block_m=A.block_m)
            sc = g.scan_angle.to_numpy(float)
            b3, _, _, _, _ = ols_cluster(np.c_[np.ones(len(g)), sc],
                                         g.d_mm_corr.to_numpy(float), blkg)
            sl[(tile, s_)] = float(b3[1])
    rows = []
    for (t, a_, b_), v in ((k, percc[k]) for k in keys):
        o = v[0] * np.pi / 180
        gsum = 0.5 * (sl[(t, a_)] + sl[(t, b_)])
        rows.append([os.path.basename(t), f"{a_}-{b_}", f"{o:+.3f}", f"{gsum:+.3f}",
                     f"{gsum - o:+.3f}"])
    R.table(["tile", "pair", "c_pair_ovl", "c_pair_gen2", "gap"], rows)

    # ------------------------------------------ 8. the repo's own boresight estimator
    print("\n## 8. The repo's own estimator (boresight.estimate_boresight) on the same offsets\n")
    rows = []
    for tile, T in tiles.items():
        df = pd.read_parquet(f"{tile}/beam_offset_table.parquet",
                             columns=["cell", "point_source_id", "scan_angle", "d_mm", "d_mm_corr",
                                      "in_grid"])
        df = df[df.in_grid.to_numpy()]
        for tag, col in (("d_mm (raw)", "d_mm"), ("d_mm_corr (registered)", "d_mm_corr")):
            sol = estimate_boresight(df.cell.values, df.point_source_id.values,
                                     df.scan_angle.values, df[col].values,
                                     min_cell_line=3, min_pair_cells=50)
            rows.append([os.path.basename(tile), tag, f"{sol.b:+.3f}", f"{sol.b_se:.3f}",
                         f"{sol.b_pair_std:.3f}", sol.n_overlap_cells])
        del df
    R.column("offset", "which offset column the estimator was run on")
    R.column("n_cells", "overlap cells the estimator used")
    R.table(["tile", "offset", "b_mmdeg", "b_se_mmdeg", "b_pairsd", "n_cells"], rows)

    # -------------------------------------------- 9. sensitivity to the one exposed cut
    print("\n## 9. Sensitivity of the shared coefficient to min_cell_line\n")
    R.column("min_cell_line", "minimum gen1 returns required for a (cell, flight-line) median")
    R.column("kept_pairs", "cell-pair observations surviving the cut, both tiles pooled")
    rows = []
    for mcl in (1, 3, 5):
        parts = []
        for tile, T in tiles.items():
            c2 = T["clm"][T["clm"].n >= mcl]
            mm = pair_rows(c2)
            if not len(mm):
                continue
            mm["blkkey"] = [f"{os.path.basename(tile)}:{v}" for v in
                            bs.block_ids(mm.cell.to_numpy(), nx=T["nx"], res=T["res"],
                                         block_m=A.block_m)]
            mm["pairkey"] = [f"{os.path.basename(tile)}:{a}-{b_}" for a, b_ in
                             zip(mm.point_source_id_a, mm.point_source_id_b)]
            parts.append(mm)
        Q = pd.concat(parts, ignore_index=True)
        Dq = pd.get_dummies(pd.Categorical(Q.pairkey)).to_numpy(float)
        Xq = np.c_[Dq, Q.dtan.to_numpy()]
        bq, Vq, r2q, nq, Gq = ols_cluster(Xq, Q.D.to_numpy(), Q.blkkey.to_numpy())
        seq = float(np.sqrt(Vq[-1, -1]))
        rows.append([mcl, len(Q), f"{bq[-1]:+.1f}", f"{seq:.1f}",
                     f"{bq[-1] / seq:+.1f}", f"{bq[-1] * np.pi / 180:+.3f}"])
    R.table(["min_cell_line", "kept_pairs", "c_shared", "c_shared_se", "t", "c_mmdeg"], rows)

    # ---------------------------------- 10. flying height, and c as an angle
    print("\n## 10. Flying height from the scan geometry, and the coefficient as an angle\n")
    print("  Within a flight line the ground easting satisfies x = x_track + h * tan(scan),\n"
          "  so regressing cell easting on tan(scan_angle) recovers the nadir track and the\n"
          "  flying height h directly from the delivered points -- no metadata needed (the\n"
          "  2008 vendor supplied no trajectory). The across-track coefficient c has units of\n"
          "  mm per unit tangent, and a pointing error delta produces a vertical error\n"
          "  h * delta * tan(scan), so delta = c / h is the equivalent roll.\n")
    R.column("x_track", "fitted nadir-track easting of the flight line, m (EPSG:26915)")
    R.column("h_m", "fitted flying height above ground, m, from x = x_track + h*tan(scan). "
                    "The SIGN alternates line to line: scan_angle is body-fixed and the lines "
                    "are flown there-and-back, so the same ground side is +scan on one line "
                    "and -scan on the next. |h| is the height.")
    R.column("h_r2", "r2 of that regression")
    rows = []
    hs = []
    for tile, T in tiles.items():
        for s_, g in T["clm"].groupby("point_source_id"):
            xc = (g.cell.to_numpy() % T["nx"]) * T["res"] + T["meta"]["bounds"][0]
            tt = np.tan(np.radians(g.sc.to_numpy(float)))
            b4 = np.polyfit(tt, xc, 1)
            pr = np.corrcoef(tt, xc)[0, 1] ** 2
            hs.append(abs(b4[0]))
            rows.append([os.path.basename(tile), str(s_), int(len(g)),
                         f"{b4[1]:.1f}", f"{b4[0]:.1f}", f"{pr:.4f}"])
    R.column("cells_used", "cell-line units the height fit used")
    R.table(["tile", "swath", "cells_used", "x_track", "h_m", "h_r2"], rows)
    hm = float(np.median(hs))
    xt = sorted(float(r[3]) for r in rows)
    print(f"\n  median fitted |flying height| = {hm:.0f} m. Nadir-track spacings: "
          + ", ".join(f"{b - a:.0f}" for a, b in zip(xt, xt[1:]) if b - a > 100) + " m.")
    print(f"  Cross-check: with sum_tan = S/h ~ 0.37 and S ~ 940 m this gives h ~ "
          f"{940 / 0.37:.0f} m, and a full swath 2*h*tan(17 deg) = "
          f"{2 * hm * np.tan(np.radians(17.0)):.0f} m against a {940:.0f} m spacing -- "
          f"a {100 * (2 * hm * np.tan(np.radians(17.0)) - 940) / (2 * hm * np.tan(np.radians(17.0))):.0f}% sidelap, "
          f"which is the overlap this test uses.")
    R.column("c_arcsec", "the pair coefficient expressed as an equivalent pointing (roll) "
                         "error, arcseconds: 206265 * c_mm / (1000 * h_m)")
    rows = [[os.path.basename(t), f"{a}-{b_}", f"{percc[(t, a, b_)][0]:+.1f}",
             f"{206265 * percc[(t, a, b_)][0] / (1000 * hm):+.1f}"] for (t, a, b_) in keys]
    R.table(["tile", "pair", "c_mm", "c_arcsec"], rows)

    R.done(headline=f"shared across-track coefficient c = {c_shared:+.1f} +- {c_shared_se:.1f} "
                    f"mm per unit tan(scan)")


if __name__ == "__main__":
    main()
