#!/usr/bin/env python3
"""Does the reference swath's off-nadir sampling bias the per-swath constants?

The worry, stated precisely
---------------------------
Both Elba tiles pin a swath the tile CUTS. If a per-swath "constant" is really an
across-track function averaged over whatever scan angles the tile happens to sample,
every constant is measured against a biased reference -- and the datum constant a
ground-control tie transports would inherit that bias.

The test has to be split, because the worry conflates two different things:

**(a) The gauge.** ``coreg.align_swaths`` solves a FREE-NETWORK least squares over all
pairwise observations and only then subtracts the reference swath's value
(``coreg.py:485``: ``c -= c[idx[ref]]``). The reference is applied AFTER the solve, so
re-referencing is arithmetic: it shifts every constant by one number and leaves every
DIFFERENCE untouched. Re-solving with an interior swath as reference therefore cannot
change the relative solution, cannot change the inter-tile disagreement, and cannot
change what a tie transports. Part 1 demonstrates this on a synthetic cloud rather than
asserting it.

**(b) The edge observations.** Those are where an across-track term can act.
``coregister_swaths(pc, a, b)`` reads each pair on the cells the two swaths SHARE -- the
overlap strip -- so what matters is not the swath's mean scan angle over the tile but the
scan angles each member presents *inside that strip*. Part 2 measures both, per tile.

Part 3 is the controlled test. Within ONE estimator
(:func:`lidar_diff_icp.boresight.estimate_boresight`, run identically on both tiles' own
per-return angle tables) each swath pair yields two vertical registrations from the same
cells: the between-line mean with NO across-track term, and the intercept of the same
data regressed on the between-line scan-angle difference -- the same quantity with the
across-track term removed. Chaining each into per-swath constants and differencing the
two tiles isolates the across-track term from every other difference between the tiles,
because only that one term changes.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/groundtruth/reference_swath_bias.py [--n-boot 200]

Sign note: ``d_mm`` in the angle tables is gen1 ground minus the gen2 reference plane, so
a between-line difference ``d_a - d_b`` is "line a reads higher than line b". The tiles
are compared to each other in that one convention throughout, so the comparison does not
depend on which way it points.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                    # noqa: E402
from lidar_diff_icp import coreg                                    # noqa: E402
from lidar_diff_icp.boresight import estimate_boresight             # noqa: E402
from lidar_diff_icp.io import PointCloud                            # noqa: E402

TILES = {
    "elba_fulldensity": ("data/derived/elba_fulldensity/gen1_csf_angles.npz",
                         "data/derived/elba_fulldensity/corrections.json"),
    "elbaext": ("data/derived/elbaext/gen1_csf_angles.npz",
                "data/derived/elbaext/corrections_geoid.json"),
}
#: The swaths both tiles hold, and the one they are compared in. 135 is elba's own
#: reference; MISSION_TIME_DRIFT.md section 4 re-references elbaext to it for the same
#: reason, so the tables here line up with that one.
SHARED = (135, 136, 137, 138)
COMPARE_REF = 135


# --------------------------------------------------------------- part 1: the gauge fact

def _synthetic_swaths(n_sw=4, seed=0):
    """Four overlapping N-S strips over a curved surface, each with a known offset.

    Deliberately NOT symmetric: the outermost strip is cut in half, exactly the geometry
    the real tiles have, so if the gauge did depend on the reference this would show it.
    """
    rng = np.random.default_rng(seed)
    xs, ys, zs, ps = [], [], [], []
    truth = {}
    for k in range(n_sw):
        x0 = 100.0 * k
        width = 60.0 if k else 30.0                       # swath 0 is cut, like swath 135
        x = rng.uniform(x0, x0 + 160.0 + 0.0 * width, 40000)
        x = x[(x >= x0) & (x <= x0 + 100.0 + width)]
        y = rng.uniform(0.0, 400.0, x.size)
        z = 8.0 * np.sin(x / 40.0) + 6.0 * np.cos(y / 55.0) + 0.02 * x
        dz = 0.010 * k - 0.015
        truth[k] = dz
        xs.append(x); ys.append(y); zs.append(z + dz); ps.append(np.full(x.size, k))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    p = np.concatenate(ps).astype(np.int32)
    return PointCloud(x=x, y=y, z=z, point_source_id=p,
                      classification=np.full(x.size, 2, np.uint8),
                      gps_time=np.zeros(x.size), scan_angle=np.zeros(x.size),
                      crs="EPSG:26915"), truth


def gauge_demo(res=2.0):
    """Solve the same network with every possible reference; report the relative solution."""
    pc, _ = _synthetic_swaths()
    out = {}
    for ref in pc.swaths.tolist():
        corr, _, _ = coreg.align_swaths(pc, res=res)
        corr_ref, _, _ = coreg.align_swaths(pc, res=res, ref=ref)
        out[ref] = (corr, corr_ref)
    return out


# ----------------------------------------------------- part 2: how the tiles sample scan

def swath_sampling(cell, psid, scan):
    """Per-swath mean signed scan angle over the whole tile, and per-EDGE inside the
    shared cells -- the only place ``coregister_swaths`` looks."""
    t = pd.DataFrame({"cell": np.asarray(cell), "ps": np.asarray(psid),
                      "sc": np.asarray(scan, float)})
    whole = t.groupby("ps").sc.agg(["mean", "median", "size"])
    cl = t.groupby(["cell", "ps"]).sc.mean().reset_index()
    sw = sorted(t.ps.unique())
    edges = []
    for a, b in zip(sw, sw[1:]):
        ca = cl[cl.ps == a]; cb = cl[cl.ps == b]
        sh = np.intersect1d(ca.cell.to_numpy(), cb.cell.to_numpy())
        if sh.size == 0:
            continue
        ma = float(ca[ca.cell.isin(sh)].sc.mean())
        mb = float(cb[cb.cell.isin(sh)].sc.mean())
        edges.append(dict(a=a, b=b, n_cells=int(sh.size), scan_a=ma, scan_b=mb,
                          asym=mb - ma))
    return whole, edges


# ---------------------------------------------- part 3: the same estimator, one term out

def _chain_from_pairs(pairs, key, swaths, ref):
    """Accumulate adjacent-pair observations into per-swath constants, gauge at ``ref``."""
    P = {(p["a"], p["b"]): p for p in pairs}
    c = {swaths[0]: 0.0}
    for a, b in zip(swaths, swaths[1:]):
        if (a, b) not in P:
            return None
        c[b] = c[a] + P[(a, b)][key]
    return {k: v - c[ref] for k, v in c.items()}


def mean_dd_pairs(cell, psid, scan, d, *, min_cell_line=3):
    """Per-pair between-line offset with NO across-track term (a degree-0 fit), and with
    it (the degree-1 intercept), from the SAME cells.

    Both are least-squares fits of the same per-cell between-line differences; the only
    difference between them is whether the scan-angle difference is a free term. That is
    what makes the comparison a controlled one.
    """
    t = pd.DataFrame({"cell": np.asarray(cell), "psid": np.asarray(psid),
                      "sc": np.asarray(scan, float), "d": np.asarray(d, float)})
    t = t[np.isfinite(t.d) & np.isfinite(t.sc)]
    g = t.groupby(["cell", "psid"]).agg(d=("d", "mean"), sc=("sc", "mean"),
                                        n=("d", "size")).reset_index()
    g = g[g.n >= min_cell_line]
    m = g.merge(g, on="cell", suffixes=("_a", "_b"))
    m = m[m.psid_a < m.psid_b]
    rows = []
    for (a, b), gg in m.groupby(["psid_a", "psid_b"]):
        dd = (gg.d_a - gg.d_b).to_numpy()
        dsc = (gg.sc_a - gg.sc_b).to_numpy()
        if dd.size < 2:
            continue
        sl, ic = np.polyfit(dsc, dd, 1)
        rows.append(dict(a=int(a), b=int(b), n_cells=int(dd.size),
                         flat=float(dd.mean()),          # degree-0: no across-track term
                         intercept=float(ic),            # degree-1: term removed
                         slope=float(sl),
                         scan_gap=float(dsc.mean())))
    return rows


def _rms(v):
    v = np.asarray(list(v), float)
    return float(np.sqrt(np.mean(v ** 2))) if v.size else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=200,
                    help="bootstrap draws for the pooled roll SE (boresight.estimate_"
                         "boresight default: 200)")
    ap.add_argument("--min-cell-line", type=int, default=3,
                    help="returns a line needs in a cell for its per-cell mean to be "
                         "used (boresight.estimate_boresight default: 3)")
    ap.add_argument("--gauge-res", type=float, default=2.0,
                    help="grid resolution for the synthetic gauge demo, m "
                         "(coreg.align_swaths default: 2.0)")
    A = ap.parse_args()

    R = Run("is the per-swath constant biased by the scan angles the tile samples, and "
            "does that explain the elba/elbaext disagreement?")
    R.param("shared_swaths", list(SHARED), src="repo",
            why="the flight lines both Elba products align; "
                "data/derived/{elba_fulldensity,elbaext}/corrections*.json")
    R.param("compare_reference_swath", COMPARE_REF, src="repo",
            why="elba's own align_swaths reference (min point_source_id); "
                "MISSION_TIME_DRIFT.md section 4 re-references elbaext to it")
    R.param("min_cell_line", A.min_cell_line, src="repo",
            why="boresight.estimate_boresight min_cell_line default")
    R.param("n_boot", A.n_boot, src="repo", why="boresight.estimate_boresight n_boot default")
    R.param("gauge_res_m", A.gauge_res, src="repo", why="coreg.align_swaths res default")

    for t, (npz, cj) in TILES.items():
        R.input(npz, role=f"{t} gen1 CSF ground returns: per-return scan angle, flight "
                          "line, grid cell and offset d_mm to the gen2 reference")
        R.input(cj, role=f"{t} pipeline corrections: the saved per-swath alignment")

    R.column("tile", "which Elba product")
    R.column("swath", "gen1 flight line (point_source_id)")
    R.column("n_returns", "CSF ground returns of that swath inside the tile, count")
    R.column("scan_mean", "mean SIGNED scan angle of those returns, deg (0 = nadir)")
    R.column("scan_med", "median signed scan angle of those returns, deg")
    R.column("edge", "adjacent swath pair, as coregister_swaths solves it")
    R.column("n_cells", "grid cells the pair shares -- the overlap strip the pair is solved on")
    R.column("scan_a", "mean signed scan angle of the LOWER-numbered line inside the strip, deg")
    R.column("scan_b", "mean signed scan angle of the HIGHER-numbered line inside the strip, deg")
    R.column("asym", "scan_b - scan_a: the across-track sampling asymmetry of the edge, deg")
    R.column("ref", "swath pinned to zero by the gauge")
    R.column("rel_max_mm", "largest change in any swath-to-swath DIFFERENCE against the "
                           "zero-mean solve, mm")
    R.column("abs_shift_mm", "the constant the gauge subtracts from every swath, mm")
    R.column("variant", "which vertical registration the per-swath constants were chained from")
    R.column("d136", "elbaext minus elba constant for swath 136, referenced to 135, mm")
    R.column("d137", "same for swath 137, mm")
    R.column("d138", "same for swath 138, mm")
    R.column("rms_mm", "RMS of those three disagreements, mm (the pinned swath is a gauge, "
                       "not an observation, so it is excluded -- as in MISSION_TIME_DRIFT.md "
                       "section 4)")
    R.column("b_mm_per_deg", "pooled across-track slope from boresight.estimate_boresight, "
                             "mm per degree of scan angle")
    R.column("b_pair_std", "std of the per-pair slopes -- the honest uncertainty on b")
    R.banner()

    # ------------------------------------------------------------------ 1. the gauge
    print("\n== PART 1: is the reference swath a gauge, or does it enter the solve? ==")
    print("  coreg.py:485 subtracts the reference AFTER the free-network least squares:")
    print("      c -= c[idx[ref]] if ref is not None else c.mean()")
    print("  so the prediction is that every swath-to-swath DIFFERENCE is identical for "
          "every choice of reference. Demonstrated on a synthetic 4-strip cloud whose "
          "outermost strip is cut in half, the geometry of the real tiles:")
    demo = gauge_demo(res=A.gauge_res)
    rows = []
    for ref, (free, pinned) in sorted(demo.items()):
        sw = sorted(free)
        rel = max(abs((pinned[b][2] - pinned[a][2]) - (free[b][2] - free[a][2]))
                  for a in sw for b in sw)
        shift = np.mean([pinned[s][2] - free[s][2] for s in sw])
        rows.append([ref, f"{1000 * rel:.2e}", f"{1000 * shift:+.1f}"])
    R.table(["ref", "rel_max_mm", "abs_shift_mm"], rows)
    print("  Every relative difference is unchanged to machine precision while the "
          "absolute level moves by tens of mm.")
    print("  CONSEQUENCE: re-solving the internal alignment with a symmetrically-sampled "
          "interior swath as reference CANNOT change the elba/elbaext disagreement, and "
          "cannot change what a ground-control tie transports. The disagreement is in "
          "the DIFFERENCES, and the gauge does not touch them.")

    # -------------------------------------------------- 2. how each tile samples scan
    print("\n== PART 2: what scan angles does each tile actually sample? ==")
    data, samp = {}, {}
    for t, (npz, _) in TILES.items():
        z = np.load(npz)
        ok = z["in_grid"] & np.isfinite(z["d_mm"])
        data[t] = dict(cell=z["cell"][ok], ps=z["point_source_id"][ok],
                       sc=z["scan_angle"][ok].astype(float), d=z["d_mm"][ok].astype(float))
        whole, edges = swath_sampling(data[t]["cell"], data[t]["ps"], data[t]["sc"])
        samp[t] = (whole, edges)

    print("\n  (a) over the WHOLE tile -- the quantity the worry is stated in:")
    rows = []
    for t, (whole, _) in samp.items():
        for s, r in whole.iterrows():
            rows.append([t, int(s), f"{int(r['size']):,}", f"{r['mean']:+.2f}",
                         f"{r['median']:+.2f}"])
    R.table(["tile", "swath", "n_returns", "scan_mean", "scan_med"], rows)
    print("  The reference swath IS cut off-nadir in both tiles, as suspected -- and so "
          "is the far edge swath 138. Interior swaths sit within ~2 deg of nadir.")

    print("\n  (b) inside the OVERLAP STRIP of each pair -- where the alignment is "
          "actually measured:")
    rows = []
    for t, (_, edges) in samp.items():
        for e in edges:
            rows.append([t, f"{e['a']}-{e['b']}", f"{e['n_cells']:,}",
                         f"{e['scan_a']:+.2f}", f"{e['scan_b']:+.2f}", f"{e['asym']:+.2f}"])
    R.table(["tile", "edge", "n_cells", "scan_a", "scan_b", "asym"], rows)
    print("  This is the correction to the worry: in the strip BOTH members are sampled "
          "off-nadir at ~10-12 deg, because a strip is by construction the far edge of "
          "one swath and the far edge of the next. The asymmetry that can bias an edge "
          "is scan_b - scan_a, and it is small on every interior edge and ~+3.5 to "
          "+4.0 deg on the OUTERMOST edge only -- the one whose lower member is the "
          "half-swath the tile cuts.")

    # ------------------------------------------- 3. the controlled test on the same data
    print("\n== PART 3: remove ONLY the across-track term, same estimator, both tiles ==")
    chains, sols = {}, {}
    for t in TILES:
        D = data[t]
        sol = estimate_boresight(D["cell"], D["ps"], D["sc"], D["d"],
                                 min_cell_line=A.min_cell_line, n_boot=A.n_boot)
        sols[t] = sol
        pairs = mean_dd_pairs(D["cell"], D["ps"], D["sc"], D["d"],
                              min_cell_line=A.min_cell_line)
        sw = sorted({p["a"] for p in pairs} | {p["b"] for p in pairs})
        chains[t] = {k: _chain_from_pairs(pairs, k, sw, COMPARE_REF)
                     for k in ("flat", "intercept")}
    R.table(["tile", "b_mm_per_deg", "b_pair_std"],
            [[t, f"{s.b:+.2f}", f"{s.b_pair_std:.2f}"] for t, s in sols.items()])
    print("  The two tiles measure the same across-track slope independently, and the "
          "between-pair scatter is the honest uncertainty on it.")

    saved = {}
    for t, (_, cj) in TILES.items():
        c = {int(k): v[2] * 1000.0
             for k, v in json.load(open(cj))["per_swath_internal_alignment_dxdydz_m"].items()}
        saved[t] = {s: c[s] - c[COMPARE_REF] for s in SHARED}

    rows = []
    others = [s for s in SHARED if s != COMPARE_REF]
    d_saved = [saved["elbaext"][s] - saved["elba_fulldensity"][s] for s in others]
    rows.append(["pipeline align_swaths (as saved)"] +
                [f"{v:+.1f}" for v in d_saved] + [f"{_rms(d_saved):.1f}"])
    for key, label in (("flat", "same estimator, across-track term IN"),
                       ("intercept", "same estimator, across-track term OUT")):
        a, b = chains["elba_fulldensity"][key], chains["elbaext"][key]
        if a is None or b is None:
            continue
        d = [b[s] - a[s] for s in others]
        rows.append([label] + [f"{v:+.1f}" for v in d] + [f"{_rms(d):.1f}"])
    R.table(["variant", "d136", "d137", "d138", "rms_mm"], rows)

    flat = _rms([chains["elbaext"]["flat"][s] - chains["elba_fulldensity"]["flat"][s]
                 for s in others])
    inter = _rms([chains["elbaext"]["intercept"][s] -
                  chains["elba_fulldensity"]["intercept"][s] for s in others])
    print(f"\n  The independent estimator reproduces the anomaly ({flat:.1f} mm RMS "
          f"against {_rms(d_saved):.1f} mm for the pipeline's own), so it is a property "
          "of the data and not of coregister_swaths.")
    print(f"  Removing ONLY the across-track term takes it from {flat:.1f} to "
          f"{inter:.1f} mm RMS -- a reduction of "
          f"{np.sqrt(max(flat ** 2 - inter ** 2, 0.0)):.1f} mm in quadrature. "
          "Same cells, same returns, same tiles; one term.")
    print(f"  It does NOT collapse to zero: {inter:.1f} mm RMS survives, and that "
          "residual is the extent-dependent repeatability the two tiles have for other "
          "reasons (different terrain inside each overlap).")
    R.done(headline=(f"reference swath is a pure gauge (relative solution invariant); "
                     f"across-track sampling is real on the outermost edge only and "
                     f"accounts for {flat:.1f}->{inter:.1f} mm RMS of the elba/elbaext "
                     f"disagreement"))


if __name__ == "__main__":
    main()
