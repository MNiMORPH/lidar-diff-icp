#!/usr/bin/env python3
"""What EXACTLY is unidentifiable in the gen1 swath network, and what would fix it?

``analysis/SWATH_ACROSS_TRACK_TEST.md`` Sec 2 established that the per-line across-track
coefficients ``c_s`` are recoverable from N-S overlaps only up to an ALTERNATING vector
``(+v, -v, +v, ...)``, and left it there.  This run does four things:

**1. Names the null space.**  With the per-line model ``e_s(th) = a_s + c_s*tan th``, each
adjacent overlap yields exactly two numbers: a slope ``P_s = (c_s + c_{s+1})/2`` and an
intercept ``K_s = (a_s - a_{s+1}) + S_s*(c_s - c_{s+1})/2``, where ``S_s = tan th_A +
tan th_B`` is the measured constant.  Because the lines fly there-and-back, ``S_s`` ALTERNATES
IN SIGN, and the claim tested here is that the null direction is not merely "alternating
``c``" but is **exactly a global cross-track TILT of the mosaic**: ``e(x) = g*(x - x0)``,
which reproduces ``c_s = g*h_s`` (alternating, because the fitted ``h_s`` alternates) together
with a LINEAR RAMP in ``a_s``.  That is verified numerically, not asserted -- a synthetic
network is built, the direction applied, and every observation checked for invariance.

**2. Measures how much the terrain lifts it in practice** (the "relief" candidate).  The
degeneracy is exact only if ``S_s`` is exactly constant.  It is not: the flying height above
ground varies with relief, so ``S = spacing/(H - z)`` varies with terrain elevation.  This
run tests that prediction against the gen2 surface, then reports the standard error the real
overlaps actually give on the antisymmetric coefficient ``q = (c_A - c_B)/2``.

**3. Prices the ground control** (the "marks" candidate).  Under the tilt reading, two
surveyed ties at different EASTINGS -- not different line parities -- constrain ``g``.  The
existing anchors 2210 and 2036 are 15.3 km apart in easting, and this run turns their
agreement into a bound on ``g`` and hence on ``v``.

**4. Prices the gen2 reference** (the "external surface" candidate).  Within a line, easting
and scan angle are the same variable, so a cross-track TILT IN GEN2 maps one-to-one onto the
per-line coefficients it would be used to measure -- with the alternating sign that is
precisely the null direction.  The size of that confound is computed from the DoD tilt the
repo has already measured.

Nothing in ``coreg.py`` or ``pipeline.py`` is modified; ``swath_across_track_test`` is
imported for its data reduction so the two runs read the identical population.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/degeneracy_identifiability.py
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lidar_diff_icp import binstats as bs
from lidar_diff_icp.registration import surface_gradients
from swath_across_track_test import cell_swath_ground, pair_rows, ols_cluster
from trust.provenance import Run

CP_CSV = "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_ql1_near_elba.csv"


# --------------------------------------------------------------- 1. null space
def null_space_check(nlines, h_mag, spacing, rng):
    """Build a synthetic N-S network, apply the candidate null direction, and report the
    largest change in ANY observation.  Zero (to machine precision) proves the direction is
    exactly unobservable; a nonzero value disproves it."""
    s = np.arange(nlines)
    sgn = (-1.0) ** s                                  # body-fixed side alternates
    h = sgn * h_mag                                    # the fitted h of Sec 10, alternating
    x_track = spacing * s
    a = rng.normal(0, 30, nlines)                      # true per-line nadir offsets, mm
    c = rng.normal(0, 100, nlines)                     # true per-line across-track coeffs

    def observations(a, c):
        """(P_s, K_s) for every adjacent pair, from the exact geometry."""
        out = []
        for i in range(nlines - 1):
            # a shared ground point at easting x sees tan th_A = (x-x_A)/h_A etc.
            x = np.linspace(max(x_track[i], x_track[i + 1]) - 0.5 * spacing,
                            min(x_track[i], x_track[i + 1]) + 0.5 * spacing, 401)
            tA = (x - x_track[i]) / h[i]
            tB = (x - x_track[i + 1]) / h[i + 1]
            D = (a[i] + c[i] * tA) - (a[i + 1] + c[i + 1] * tB)
            X = np.c_[np.ones_like(tA), tA - tB]
            beta = np.linalg.lstsq(X, D, rcond=None)[0]
            out.append((beta[0], beta[1], float((tA + tB).mean())))
        return np.array(out)

    O0 = observations(a, c)
    g = 0.01                                           # a global tilt of 10 mm/km, mm per m
    a2 = a + g * (x_track - x_track[0])
    c2 = c + g * h
    O1 = observations(a2, c2)
    dc = c2 - c

    # --- the same direction, now seen by a PERPENDICULAR (due-west) line.
    # The cross line's across-track coordinate is NORTHING, so a tilt in EASTING is not
    # expressible as (a_c + c_c*tan th_c): it lands on the cross-pair fit instead.
    def cross_obs(a, c, i, a_c, c_c):
        """Fit D = k + p*(tA-tC) + q*(tA+tC) on the patch line i shares with the cross line."""
        xx = np.linspace(x_track[i] - 0.28 * abs(h[i]), x_track[i] + 0.28 * abs(h[i]), 121)
        yy = np.linspace(-0.28 * h_mag, 0.0, 121)          # the cross line's own across-track
        X2, Y2 = np.meshgrid(xx, yy)
        tA = (X2.ravel() - x_track[i]) / h[i]
        tC = Y2.ravel() / h_mag
        D = (a[i] + c[i] * tA) - (a_c + c_c * tC)
        M = np.c_[np.ones_like(tA), tA - tC, tA + tC]
        beta = np.linalg.lstsq(M, D, rcond=None)[0]
        e = D - M @ beta
        return beta, float(np.sqrt((e ** 2).mean()))

    a_c, c_c = 12.0, 55.0
    cross0 = [cross_obs(a, c, i, a_c, c_c) for i in range(nlines)]
    # under the tilt, the cross line's TRUE error at a point is g*(x-x0): it varies along the
    # cross line's own ALONG-track direction, which its two parameters cannot absorb.
    def cross_obs_tilted(i):
        xx = np.linspace(x_track[i] - 0.28 * abs(h[i]), x_track[i] + 0.28 * abs(h[i]), 121)
        yy = np.linspace(-0.28 * h_mag, 0.0, 121)
        X2, Y2 = np.meshgrid(xx, yy)
        tA = (X2.ravel() - x_track[i]) / h[i]
        tC = Y2.ravel() / h_mag
        eA = a[i] + c[i] * tA + g * (X2.ravel() - x_track[0])
        eC = a_c + c_c * tC + g * (X2.ravel() - x_track[0])
        D = eA - eC
        M = np.c_[np.ones_like(tA), tA - tC, tA + tC]
        beta = np.linalg.lstsq(M, D, rcond=None)[0]
        e = D - M @ beta
        return beta, float(np.sqrt((e ** 2).mean()))
    cross1 = [cross_obs_tilted(i) for i in range(nlines)]

    return dict(O0=O0, O1=O1, dc=dc, da=a2 - a, g=g, h=h, x_track=x_track,
                cross0=cross0, cross1=cross1, c=c, a=a, a_c=a_c, c_c=c_c,
                dmax=float(np.abs(O1[:, :2] - O0[:, :2]).max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", nargs="+",
                    default=["data/derived/elba_fulldensity", "data/derived/elbaext"])
    ap.add_argument("--min-cell-line", type=int, default=1)
    ap.add_argument("--block-m", type=float, default=50.0)
    ap.add_argument("--h-m", type=float, default=2562.0)
    A = ap.parse_args()

    R = Run("what exactly is unidentifiable in the gen1 swath network, and which available "
            "observation would determine it?")
    R.param("min_cell_line", A.min_cell_line, src="repo",
            why="the same definitional floor swath_across_track_test.py runs at, so the two "
                "runs read the identical population")
    R.param("block_m", A.block_m, src="repo")
    R.param("h_m", A.h_m, src="repo",
            why="flying height, the median |h| of the per-line fits in "
                "analysis/SWATH_ACROSS_TRACK_TEST.md Sec 0")
    ties_p = R.input("data/derived/groundtruth/elba_gen1_ties.json",
                     role="gen1 absolute ties at the surveyed 3DEP checkpoints, mm to ADD to "
                          "gen1, with each mark's sigma and chain path")
    inv_p = R.input("analysis/ELBAEXT2_SCOPE.md",
                    role="the fitted nadir-track eastings and per-mark cross-track distances "
                         "quoted below (Sec 2), read not recomputed")

    for c, d in [
        ("pair", "adjacent flight-line pair index in the synthetic network, dimensionless"),
        ("P", "the pair's fitted slope (c_A+c_B)/2, mm per unit tangent"),
        ("K", "the pair's fitted intercept, mm"),
        ("dP", "change in P under the candidate null direction, mm per unit tangent"),
        ("dK", "change in K under the candidate null direction, mm"),
        ("line", "flight line index in the synthetic network, dimensionless"),
        ("dc", "change in that line's across-track coefficient under the direction, "
               "mm per unit tangent"),
        ("da", "change in that line's nadir offset under the direction, mm"),
        ("c_true", "the synthetic line's true across-track coefficient, mm per unit tangent"),
        ("c_fit_notilt", "what the cross-line pair fit returns for it with no tilt present, "
                         "mm per unit tangent"),
        ("c_fit_tilted", "the same with a 10 mm/km global cross-track tilt added to both "
                         "lines, mm per unit tangent"),
        ("rms_notilt", "root-mean-square residual of that cross-line fit, mm"),
        ("rms_tilted", "the same with the tilt present, mm"),
        ("tile", "derived-product directory the overlap cells come from"),
        ("pairname", "unordered flight-line pair (point_source_id)"),
        ("cells", "overlap reference cells used"),
        ("sd_sum", "standard deviation of tan(scan_A)+tan(scan_B) over those cells"),
        ("sd_dif", "standard deviation of tan(scan_A)-tan(scan_B) over those cells"),
        ("p_mm", "fitted (c_A+c_B)/2, mm per unit tangent"),
        ("p_se", "its cluster-robust standard error, mm per unit tangent"),
        ("q_mm", "fitted (c_A-c_B)/2, mm per unit tangent -- the degenerate combination"),
        ("q_se", "its cluster-robust standard error, mm per unit tangent"),
        ("se_ratio", "q_se / p_se: how many times worse the degenerate combination is known"),
        ("dsum_dz", "OLS slope of tan(scan_A)+tan(scan_B) on the gen2 ground elevation of the "
                    "same cell, per metre"),
        ("dsum_dz_pred", "the value predicted by geometry, mean(sum_tan)/(h - mean z), per metre"),
        ("r2_z", "r-squared of that regression, dimensionless"),
        ("mark", "surveyed 3DEP checkpoint id"),
        ("psid_mark", "the gen1 flight line that overflies the mark"),
        ("east", "the mark's easting, m (UTM 15N)"),
        ("offnadir", "the mark's distance from that line's fitted nadir track, m"),
        ("tan_mark", "tan of the scan angle at which that line sees the mark = offnadir / h"),
        ("tie", "surveyed - gen1, mm (the constant to ADD to gen1)"),
        ("tie_sigma", "that tie's own standard error, mm"),
        ("c_effect", "how much a per-line coefficient of 130 mm per unit tangent moves the "
                     "tie at that mark: 130 * tan_mark, mm"),
    ]:
        R.column(c, d)
    R.banner()
    print()

    # ------------------------------------------------------------- 1. null space
    print("## 1. Is the null direction exactly a global cross-track tilt? (synthetic, exact)")
    rng = np.random.default_rng(0)
    N = null_space_check(6, A.h_m, 916.0, rng)
    R.table(["pair", "P", "K", "dP", "dK"],
            [[i, f"{N['O0'][i,1]:+.3f}", f"{N['O0'][i,0]:+.3f}",
              f"{N['O1'][i,1]-N['O0'][i,1]:+.3e}", f"{N['O1'][i,0]-N['O0'][i,0]:+.3e}"]
             for i in range(len(N["O0"]))])
    print(f"\nlargest change in ANY pair observation under a {N['g']*1000:.0f} mm/km global "
          f"tilt: {N['dmax']:.3e} mm")
    print("the direction, line by line -- note dc ALTERNATES and da is a LINEAR RAMP:")
    R.table(["line", "dc", "da"],
            [[i, f"{N['dc'][i]:+.2f}", f"{N['da'][i]:+.2f}"] for i in range(len(N["dc"]))])
    print(f"\nso v (the alternating amplitude) = g*|h| = {N['g']:.4f} mm/m * {A.h_m:.0f} m = "
          f"{N['g']*A.h_m:.1f} mm per unit tangent,")
    print(f"and the ramp step per line = g*spacing = {N['g']*916.0:.2f} mm.")
    print("VERDICT: the alternating-c mode and a global cross-track tilt are THE SAME OBJECT.")
    print("A tilt is identical for every line at a shared point, so it cancels in EVERY")
    print("between-line difference -- which is why no set of overlaps, adjacent or not, and")
    print("no cross line, can ever see it. Only an external absolute reference can.")
    print()
    print("## 1b. The SAME synthetic, now with a due-west cross line over each N-S line.")
    print("Fitting D = k + p*(tA-tC) + q*(tA+tC) recovers c_A = p+q for the N-S line.")
    R.table(["line", "c_true", "c_fit_notilt", "c_fit_tilted", "rms_notilt", "rms_tilted"],
            [[i, f"{N['c'][i]:+.2f}",
              f"{N['cross0'][i][0][1]+N['cross0'][i][0][2]:+.2f}",
              f"{N['cross1'][i][0][1]+N['cross1'][i][0][2]:+.2f}",
              f"{N['cross0'][i][1]:.2e}", f"{N['cross1'][i][1]:.2e}"]
             for i in range(len(N['c']))])
    print("The cross line returns EACH line's own c exactly, with no alternating freedom --")
    print("that is the degeneracy broken. It returns the SAME value with the tilt added,")
    print("i.e. it is blind to the tilt rather than confused by it: a cross line fixes the")
    print("BODY-FIXED per-line term and leaves the mosaic tilt to ground control.")
    print()

    # -------------------------------------- 2. what the real overlaps give on q
    print("## 2. What the real Elba overlaps give on the degenerate combination q")
    rows, zrows = [], []
    for tile in A.tiles:
        meta = json.load(open(R.input(f"{tile}/corrections.json",
                                      role="pipeline corrections: grid geometry, per-swath "
                                           "alignment, datum, drift")))
        res = float(meta["res_m"])
        clm, _, _, _, _ = cell_swath_ground(tile, R, A.min_cell_line, res)
        b = meta["bounds"]
        nx = int(round((b[2] - b[0]) / res))
        z = np.load(f"{tile}/z_after.npy").ravel()
        m = pair_rows(clm)
        m["blk"] = bs.block_ids(m.cell.to_numpy(), nx=nx, res=res, block_m=A.block_m)
        m["z"] = z[m.cell.to_numpy()]
        for (a, bb), gg in m.groupby([m.point_source_id_a, m.point_source_id_b]):
            u = gg.dtan.to_numpy()
            S = gg.stan.to_numpy()
            X = np.c_[np.ones_like(u), u, S]
            beta, V, _, _, _ = ols_cluster(X, gg.D.to_numpy(), gg.blk.to_numpy())
            se = np.sqrt(np.diag(V))
            rows.append([os.path.basename(tile), f"{a}-{bb}", len(gg), f"{S.std():.4f}",
                         f"{u.std():.4f}", f"{beta[1]:+.1f}", f"{se[1]:.1f}",
                         f"{beta[2]:+.1f}", f"{se[2]:.1f}", f"{se[2]/se[1]:.1f}"])
            zz = gg.z.to_numpy()
            fin = np.isfinite(zz)
            sl = np.polyfit(zz[fin], S[fin], 1)[0]
            pred = float(np.sign(S.mean()) * abs(S.mean()) / (A.h_m - np.nanmean(zz)))
            pred = float(S.mean() / (A.h_m - np.nanmean(zz)))
            r2 = float(np.corrcoef(zz[fin], S[fin])[0, 1] ** 2)
            zrows.append([os.path.basename(tile), f"{a}-{bb}", int(fin.sum()),
                          f"{sl:+.3e}", f"{pred:+.3e}", f"{r2:.3f}"])
    R.table(["tile", "pairname", "cells", "sd_sum", "sd_dif", "p_mm", "p_se", "q_mm", "q_se",
             "se_ratio"], rows)
    print()
    print("## 2b. Is that residual sd_sum REAL relief, or scan-angle quantisation noise?")
    print("Geometry predicts sum_tan = spacing/(h - z), so d(sum_tan)/dz = sum_tan/(h - z).")
    R.table(["tile", "pairname", "cells", "dsum_dz", "dsum_dz_pred", "r2_z"], zrows)
    print()

    # ------------------------------------------------- 3. ground control on the tilt
    print("## 3. What the surveyed ties say about the tilt g (hence about v)")
    T = json.load(open(ties_p))
    # eastings and off-nadir distances: read from ELBAEXT2_SCOPE.md Sec 2, not recomputed
    cps = pd.read_csv(R.input(CP_CSV, role="the surveyed 3DEP checkpoints near Elba: id, type, "
                                            "UTM 15N easting/northing, NAVD88(GEOID18) elevation"))
    EAST = dict(zip(cps.point_id, cps.easting.astype(float)))
    OFFNADIR = {"2210_2021_MN": 147.0, "3056_2021_MN": 131.0,
                "2024_2021_MN": 32.0, "2036_2021_MN": 227.0}    # ELBAEXT2_SCOPE.md Sec 2
    trows, anchors = [], []
    for cp in T["checkpoints"]:
        pid = cp["point_id"]
        if not cp.get("attempted") or pid not in EAST:
            continue
        tn = OFFNADIR[pid] / A.h_m
        trows.append([pid, cp["line"], f"{EAST[pid]:.0f}", f"{OFFNADIR[pid]:.0f}",
                      f"{tn:.4f}", f"{cp['tie_mm']:+.1f}", f"{cp['sigma_mm']:.1f}",
                      f"{130.0*tn:+.1f}"])
        if pid in ("2210_2021_MN", "2036_2021_MN"):
            anchors.append((EAST[pid], cp["tie_mm"], cp["sigma_mm"]))
    R.table(["mark", "psid_mark", "east", "offnadir", "tan_mark", "tie", "tie_sigma",
             "c_effect"], trows)
    (x1, t1, s1), (x2, t2, s2) = anchors
    dx = x2 - x1
    # tie = surveyed - lidar, so the gen1 error field is e = -tie; g = de/dx
    g = -(t2 - t1) / dx
    sg = np.hypot(s1, s2) / abs(dx)
    print(f"\nbaseline between the two anchors: {abs(dx)/1000:.2f} km in easting")
    print(f"gen1 cross-track tilt implied by the two ties: g = {g*1000:+.3f} +- {sg*1000:.3f} "
          f"mm/km")
    print(f"equivalent alternating amplitude v = g*h = {g*A.h_m:+.2f} +- {sg*A.h_m:.2f} "
          f"mm per unit tangent")
    print(f"2-sigma bound on |v|: {abs(g*A.h_m) + 2*sg*A.h_m:.1f} mm per unit tangent, against "
          f"per-pair coefficients of +34 to +193 (SWATH_ACROSS_TRACK_TEST Sec 1)")
    print("CAVEAT, from ABSOLUTE_BASIS_ELBA.md's own budget: the ties carry a 42.6 mm")
    print("UNMODELLED bound (validity of the lateral shift 7-16 km out) that is deliberately")
    print("kept out of sigma. Propagated the same way that bound alone gives")
    print(f"sigma(v) = {42.6*np.sqrt(2)/abs(dx)*A.h_m:.1f} mm per unit tangent.")
    print()

    # ------------------------------------------------- 4. gen2 as external reference
    print("## 4. gen2 as the external reference: the size of its own tilt confound")
    tilt = json.load(open(R.input("data/derived/elba_fulldensity/z_before_absolute.json",
                                  role="the gen1 absolute datum product and its budget")))
    print(json.dumps({k: v for k, v in tilt.items()
                      if k in ("constant_mm", "sign_convention")}, indent=1))
    for name, val, sd in (("DoD tilt dE, stable ground (ADDITIONAL_GROUND_CONTROL.md)",
                           -14.19, 5.15),):
        print(f"{name}: {val:+.2f} +- {sd:.2f} mm/km")
        print(f"   -> as a per-line coefficient it is g*h = {val/1000*A.h_m:+.1f} +- "
              f"{sd/1000*A.h_m:.1f} mm per unit tangent, ALTERNATING in sign with line parity")
        print(f"   -> i.e. exactly the null direction, at {abs(val/1000*A.h_m):.0f} mm per unit "
              f"tangent")
    R.done(headline="the swath-network null space is a global cross-track tilt; a cross line "
                    "fixes the body-fixed per-line coefficients, ground control fixes the tilt")


if __name__ == "__main__":
    main()
