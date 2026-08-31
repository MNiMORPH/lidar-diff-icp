"""Is the near/far mark split a DISTANCE effect, and were our marks held out from the
vendor's bias adjustment?

Two questions that between them govern how much the open items in FRAME.md matter.

1. **Near/far.** Six of the datum's eight open marks sit far from the site and disagree
   with the two near ones. If that is a distance effect the estimator is wrong; if it is
   the per-line structure the estimator already averages over, it is handled. The two are
   confounded, so the only clean test is WITHIN a line.
2. **Hold-out.** Both epochs carry an unpublished vendor bias adjustment. Our constant
   measures what remains AFTER it, so its value is never needed -- PROVIDED our marks were
   held out from the calibration. Had they been the calibration set, their residuals would
   sit on zero by construction.
"""
from __future__ import annotations
import argparse, json, sys, collections
from pathlib import Path
import numpy as np
from scipy import stats
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE)); sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))
import control  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G  # noqa: E402
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", required=True)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--near-km", type=float, required=True)
    a = p.parse_args(argv)

    R = Run("Is the near/far mark split a distance effect or the per-line structure, and "
            "were our marks held out from the vendor's bias adjustment?")
    R.input(a.product, role="the datum's own open marks, with the line their returns assign")
    R.input(str(control.GEN1_CSV), role="gen1's own 2008 control; dnr_error_m is the "
                                        "vendor's published residual on the DELIVERED "
                                        "surface, i.e. AFTER its bias adjustment")
    R.param("near_km", a.near_km, src="MINE",
            why="the near/far cut; it is REPORTED as confounded rather than used to "
                "select, and the within-line comparison below does not depend on it")
    R.column("line", "flight line point_source_id")
    R.column("km", "distance from the site, km")
    R.column("tie_mm", "surveyed - lidar at the mark, mm")
    R.column("stat", "name of the statistic")
    R.column("value", "its value; mm unless the name says otherwise")
    R.banner()

    ctl = {m.aliases[0]: m for m in G.load_control()}
    rows = []
    for m in json.loads(Path(a.product).read_text())["marks"]:
        mk = ctl[m["point_id"]]
        km = float(np.hypot(mk.checkpoint.easting - a.easting,
                            mk.checkpoint.northing - a.northing) / 1000)
        rows.append((int(m["line"]), km, float(m["tie_mm"])))
    R.table(["line", "km", "tie_mm"],
            [[l, f"{k:.1f}", f"{t:+.1f}"] for l, k, t in sorted(rows)])

    near = np.array([t for l, k, t in rows if k <= a.near_km])
    far = np.array([t for l, k, t in rows if k > a.near_km])
    tt, pp = stats.ttest_ind(near, far, equal_var=False)
    by = collections.defaultdict(list)
    for l, k, t in rows:
        by[l].append((k, t))
    within = [(l, sorted(v)) for l, v in by.items()
              if len(v) > 1 and sorted(v)[-1][0] - sorted(v)[0][0] > 5]
    out = [["near n / far n", f"{near.size} / {far.size}"],
           ["near mean / far mean", f"{near.mean():+.2f} / {far.mean():+.2f}"],
           ["Welch t", f"{tt:+.3f}"],
           ["p", f"{pp:.3f}"],
           ["lines carrying a NEAR mark", str(sorted({l for l, k, t in rows if k <= a.near_km}))],
           ["lines carrying only FAR marks", str(sorted({l for l, k, t in rows if k > a.near_km}
                                                        - {l for l, k, t in rows if k <= a.near_km}))]]
    for l, v in sorted(within):
        out.append([f"WITHIN line {l}: far minus near over {v[-1][0]-v[0][0]:.1f} km",
                    f"{v[-1][1]-v[0][1]:+.1f}"])
    L = control.load_control("gen1"); r = L.residuals
    t0, p0 = stats.ttest_1samp(r.resid_mm, 0.0)
    out += [["hold-out check: all 963 published residuals, mean", f"{r.resid_mm.mean():+.2f}"],
            ["hold-out check: t against zero", f"{t0:+.2f}"],
            ["hold-out check: p", f"{p0:.3e}"]]
    R.table(["stat", "value"], out)
    R.done(headline=f"near/far p={pp:.3f} and confounded with line; hold-out p={p0:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
