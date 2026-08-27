"""Assemble the canonical gen1-at-Elba answer: both surfaces, with every choice declared.

Reads the committed products rather than recomputing, so this file IS the answer and its
inputs are named:

  * ``gen1_datum_by_returns_elba.json``  -- the catchment-free estimate
  * ``bridge_wide_L1O.json``             -- the delivered -> our-surface bridge

Two surfaces are reported and they are NOT interchangeable:

  DELIVERED  what the vendor shipped; what the 2008 control measured directly.
  OURS       the CSF-reclassified, swath-aligned, geoid-shifted 5 m grid the DoD runs on.
             = delivered + bridge.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datum-product", required=True)
    p.add_argument("--bridge-product", required=True)
    p.add_argument("--sigma-sweep", type=float, nargs="+", required=True)
    p.add_argument("--sigma-values-mm", type=float, nargs="+", required=True)
    p.add_argument("--sigma-se-mm", type=float, nargs="+", required=True)
    p.add_argument("--adopted-sigma", type=float, required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    R = Run("What is gen1's datum at Elba, on the delivered surface and on ours, with "
            "every choice that moves it declared?")
    R.input(a.datum_product, role="catchment-free per-line datum; marks assigned by their "
                                  "own ground returns, reused psids disambiguated by "
                                  "track collinearity")
    R.input(a.bridge_product, role="per-mark bridge on open ground; its population mean "
                                   "carries a delivered-surface constant onto ours")
    R.param("adopted_collinear_sigma", a.adopted_sigma, src="MINE",
            why="the estimate is flat across sigma 2-3 (0.47 mm) and jumps 8.22 mm at "
                "sigma 5 where the disambiguation readmits marks on non-collinear passes; "
                "3 sits at the top of the stable range")
    R.param("catchment_radius_m", "NONE", src="andy",
            why="removed: it was a compute bound, not a criterion, and it admitted marks "
                "the pass test rejects")
    R.param("cover", "L1O+L5U as run", src="andy",
            why="NOT yet narrowed to open-only; doing so is an open decision that will "
                "move this number")
    R.param("gen2_side", "NOT carried onto our surface", src="repo",
            why="gen2's bridge is not measurable at its checkpoints (engineered siting), "
                "so a DoD correction built from these two constants is asymmetric")
    R.column("quantity", "what is being reported")
    R.column("value_mm", "millimetres; positive = the surface reads LOW, so ADD it")
    R.column("uncert_mm", "its uncertainty, mm")
    R.column("uncert_of", "what that uncertainty is the uncertainty OF")
    R.banner()

    d = json.loads(Path(a.datum_product).read_text())
    b = np.array([m["bridge_mm"] for m in
                  json.loads(Path(a.bridge_product).read_text())["marks"]])
    bm, bse = float(b.mean()), float(b.std(ddof=1) / np.sqrt(b.size))
    v, se = float(d["value_mm"]), float(d["se_mm"])
    tot = float(np.hypot(se, bse))
    sw = np.array(a.sigma_values_mm)

    rows = [
        ["gen1 at Elba, DELIVERED surface", f"{v:+.2f}", f"{se:.2f}",
         f"SE over the {d['n_lines']} flight lines ({d['n_marks']} marks)"],
        ["bridge, delivered -> ours", f"{bm:+.2f}", f"{bse:.2f}",
         f"SE of the mean over {b.size} open marks"],
        ["gen1 at Elba, OUR surface", f"{v + bm:+.2f}", f"{tot:.2f}",
         "the two in quadrature"],
        ["sigma sensitivity", f"{np.ptp(sw):+.2f}", "--",
         f"range of the datum over collinear_sigma {a.sigma_sweep}"],
    ]
    R.table(["quantity", "value_mm", "uncert_mm", "uncert_of"], rows)

    out = dict(
        site="elba", easting=578762.8, northing=4884487.6,
        delivered_mm=v, delivered_se_mm=se,
        bridge_mm=bm, bridge_se_mm=bse, n_bridge_marks=int(b.size),
        our_surface_mm=v + bm, our_surface_se_mm=tot,
        n_marks=d["n_marks"], n_lines=d["n_lines"],
        adopted_collinear_sigma=a.adopted_sigma,
        sigma_sweep={str(s): {"value_mm": vv, "se_mm": ss} for s, vv, ss in
                     zip(a.sigma_sweep, a.sigma_values_mm, a.sigma_se_mm)},
        sign_convention="tie = surveyed - z_lidar; POSITIVE = surface reads LOW, so ADD",
        open_choices=["cover is L1O+L5U, not yet narrowed to open-only",
                      "collinear_sigma moves the answer 8.22 mm between 3 and 5",
                      "gen2's bridge is not measurable, so the DoD correction is "
                      "asymmetric between epochs"],
        producers=["ground_control/run_datum_by_returns.py",
                   "ground_control/run_bridge_wide.py"])
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(f"\n  wrote {a.out}")
    R.done(headline=f"gen1 at Elba: delivered {v:+.2f} +/- {se:.2f}; "
                    f"our surface {v+bm:+.2f} +/- {tot:.2f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
