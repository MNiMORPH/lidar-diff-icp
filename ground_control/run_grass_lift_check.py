"""Does the ground-return median float on grass AT THE CONTROL MARKS the datum uses?

Raised by the other session (issue #1): gen1's per-cell ground is the MEDIAN of its CSF
ground returns, and that median appears to float on dense grass, which would push the
apparent surface UP and make the tie look SMALLER.

Measured here on the marks the datum actually rests on, with the SAME ground source those
ties used (vendor class 2), detrended on the local order-2 surface first so terrain slope
cannot masquerade as lift.

Note on interpretation: ``p50 - p10`` is a SPREAD, not a bias. And the DoD runs on the
median-of-ground-returns surface, so a constant measured against that surface is the
constant for the surface in use -- the lift is part of what is being calibrated, not an
error in the calibration. The risk it does create is REPRESENTATIVENESS: if the marks sit
on firmer ground than the tile at large, the constant measured at them does not apply to
the tile.
"""
from __future__ import annotations
import argparse, glob, json, sys
from pathlib import Path
import numpy as np
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent)); sys.path.insert(0, str(_HERE.parent / "src"))
from lidar_diff_icp.groundtruth import gen1_datum as G, tie as T  # noqa: E402
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", required=True)
    p.add_argument("--tiles", required=True)
    p.add_argument("--radius-m", type=float, required=True)
    p.add_argument("--read-half-width-m", type=float, required=True)
    p.add_argument("--reference-open-upland-mm", type=float, required=True)
    a = p.parse_args(argv)
    import laspy

    R = Run("Does the ground-return median float on grass at the control marks the gen1 "
            "datum rests on, and by how much relative to open ground generally?")
    R.input(a.product, role="the marks the adopted gen1 constant is measured from")
    R.param("radius_m", a.radius_m, src="MINE",
            why="disc within which returns are taken; matches the tie estimator's own "
                "pipeline-scale radius so the statistic describes the ties' own support")
    R.param("ground_source", "vendor class 2", src="repo",
            why="the source those ties used; measuring the lift on a DIFFERENT ground "
                "source would not describe them")
    R.param("detrend", "order-2 surface removed first", src="MINE",
            why="without it terrain slope inside the disc masquerades as vertical spread")
    R.param("reference_open_upland_mm", a.reference_open_upland_mm, src="andy",
            why="the other session's gen1-internal open-upland figure, for comparison; "
                "their statistic is PER CELL and this one is per disc, so the two are "
                "indicative rather than strictly like-for-like")
    R.column("point_id", "control mark id")
    R.column("n_ground", "vendor class-2 returns within radius_m, count")
    R.column("p50_p10_mm", "median minus 10th percentile of the DETRENDED return "
                           "residuals, mm -- a spread, NOT a bias")
    R.column("p95_p05_mm", "95th minus 5th percentile of the same, mm")
    R.column("stat", "name of the statistic")
    R.column("value", "its value; mm unless stated")
    R.banner()

    ctl = {m.aliases[0]: m for m in G.load_control()}
    marks = json.loads(Path(a.product).read_text())["marks"]
    boxes = {}
    for q in glob.glob(f"{a.tiles}/*.laz"):
        with laspy.open(q) as fh:
            h = fh.header
            boxes[q] = (h.mins[0], h.mins[1], h.maxs[0], h.maxs[1])
    rows, lifts = [], []
    for m in marks:
        mk = ctl[m["point_id"]]
        E, N = mk.checkpoint.easting, mk.checkpoint.northing
        tp = next((q for q, b in boxes.items()
                   if b[0] <= E <= b[2] and b[1] <= N <= b[3]), None)
        if tp is None:
            continue
        g = T.vendor_ground_near(tp, E, N, a.read_half_width_m)
        x, y, z = np.asarray(g.x), np.asarray(g.y), np.asarray(g.z)
        s = np.hypot(x - E, y - N) <= a.radius_m
        if s.sum() < 10:
            continue
        u, v, zz = x[s] - E, y[s] - N, z[s]
        A = np.column_stack([np.ones(u.size), u, v, u * u, v * v, u * v])
        c, *_ = np.linalg.lstsq(A, zz, rcond=None)
        r = (zz - A @ c) * 1000.0
        lift = float(np.percentile(r, 50) - np.percentile(r, 10))
        lifts.append(lift)
        rows.append([m["point_id"], int(s.sum()), f"{lift:.1f}",
                     f"{np.percentile(r,95)-np.percentile(r,5):.1f}"])
    R.table(["point_id", "n_ground", "p50_p10_mm", "p95_p05_mm"], rows)
    L = np.array(lifts)
    R.table(["stat", "value"], [
        ["marks measured (count)", f"{L.size}"],
        ["median p50-p10 at the marks", f"{np.median(L):.1f}"],
        ["minimum / maximum", f"{L.min():.1f} / {L.max():.1f}"],
        ["reference: open upland, per cell", f"{a.reference_open_upland_mm:.1f}"],
        ["marks minus reference", f"{np.median(L) - a.reference_open_upland_mm:+.1f}"]])
    R.done(headline=f"median p50-p10 {np.median(L):.1f} mm at the datum's marks against "
                    f"{a.reference_open_upland_mm:.1f} mm for open upland generally")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
