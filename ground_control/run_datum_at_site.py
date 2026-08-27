"""Driver: the datum constant at a site, per epoch, from each epoch's OWN control.

    ./lidar-icp/bin/python ground_control/run_datum_at_site.py \
        --easting 578762.8 --northing 4884487.6 --site elba \
        --gen2-surfaces ql1_laz ql0_laz \
        --max-lags-m 20000 40000 80000 160000 --n-lags 25 --n-pairs 800000 \
        --estimators dowd matheron --seed 0 \
        --out-dir ground_control/products

No defaults for the variogram sweep: the control does not determine a range, so a single
max_lag would privilege one unfalsifiable fit.  ``--gen2-surfaces`` is required because
gen2 publishes four surfaces and they are four different answers.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))

import control  # noqa: E402
import datum  # noqa: E402
sys.path.insert(0, str(_HERE.parent))
from trust.provenance import Run  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--site", required=True, help="name, used in the product filenames")
    p.add_argument("--gen2-surfaces", nargs="+", required=True,
                   choices=control.GEN2_SURFACES)
    p.add_argument("--max-lags-m", type=float, nargs="+", required=True)
    p.add_argument("--n-lags", type=int, required=True)
    p.add_argument("--n-pairs", type=int, required=True)
    p.add_argument("--estimators", nargs="+", required=True,
                   choices=("dowd", "matheron"))
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--gen1-bridge-mm", type=float, default=None)
    p.add_argument("--gen1-bridge-sd-mm", type=float, default=None)
    p.add_argument("--gen2-bridge-mm", type=float, default=None)
    p.add_argument("--gen2-bridge-sd-mm", type=float, default=None)
    return p.parse_args(argv)


def banner(a):
    """Declare every input, parameter and column through trust/provenance.py.

    The repo's rule: the label on a number is written by the code that made the number,
    not typed afterwards.  A run that prints a results table must print its provenance.
    """
    R = Run("What constant places each epoch's surface on surveyed NAVD88 at this site, "
            "using ONLY that epoch's own contemporaneous control?")
    R.input(str(control.GEN1_CSV),
            role="gen1's OWN 2008 MnGeo/MnDNR validation control; the column "
                 "dnr_error_m is Control Z - Surface Z, verified exact on 1004 of 1004 rows")
    R.input(str(control.GEN2_CSV),
            role="gen2's OWN 2021 USGS NVA/VVA held-out checkpoints; the usgs_*_error_m "
                 "columns are Z - <surface>z, verified exact on all four columns")
    R.param("site_easting", a.easting, src="andy", why="the site asked for")
    R.param("site_northing", a.northing, src="andy")
    R.param("epoch_matching", "gen1<-2008 control, gen2<-2021 control", src="andy",
            why="Andy's directive: adjust each epoch internally, no cross-epoch control")
    R.param("gen2_surfaces", tuple(a.gen2_surfaces), src="MINE",
            why="gen2 publishes FOUR surfaces and they are four different answers; "
                "these are reported side by side, not reduced to one")
    R.param("max_lags_m", tuple(a.max_lags_m), src="MINE",
            why="swept, not fitted once: the control does not determine a variogram "
                "range and single fits pin to the largest lag centre")
    R.param("n_lags", a.n_lags, src="MINE",
            why="lag bins in the empirical variogram; excludes nothing -- it sets the "
                "resolution at which the sweep can see structure, not which marks enter")
    R.param("n_pairs", a.n_pairs, src="MINE",
            why="random pairs drawn per variogram; excludes no marks, but a small value "
                "would add sampling noise to the fit and so to sd_field")
    R.param("estimators", tuple(a.estimators), src="MINE",
            why="dowd is robust, matheron is classical; both are run so the "
                "nugget/sill partition's instability is visible")
    R.param("seed", a.seed, src="MINE",
            why="fixes the random pair draw so the run is reproducible; changing it "
                "re-draws pairs and moves the fit slightly, not the marks")
    R.param("cover_treatment", "ALL swept; none chosen", src="MINE",
            why="the largest single lever in the problem; choosing one here would "
                "silently set the datum")
    R.param("bridge_mm", (a.gen1_bridge_mm, a.gen2_bridge_mm), src="MINE",
            why="delivered-surface -> our reconstructed surface; None means NOT "
                "applied and NOT assumed zero")
    R.column("epoch", "which acquisition the constant is for")
    R.column("selection", "gen2 surface (if any) and cover treatment")
    R.column("n_marks", "control marks entering the kriging system")
    R.column("constant_mm", "constant to ADD to the DELIVERED surface at the site, mm; "
                            "median over the variogram sweep; positive = surface reads LOW")
    R.column("sweep_min_mm", "smallest constant over the variogram sweep, mm")
    R.column("sweep_max_mm", "largest constant over the variogram sweep, mm")
    R.column("sd_field_mm", "sd of the error in predicting the spatially correlated "
                            "component of the surface offset at this coordinate, mm; "
                            "median over the sweep; NOT an SE of a mean over marks")
    R.column("sd_field_min_mm", "smallest sd_field over the sweep, mm")
    R.column("sd_field_max_mm", "largest sd_field over the sweep, mm")
    R.notes.append("gen1 reads ONLY the 2008 table and gen2 ONLY the 2021 table. No "
           "cross-epoch control is read and no flight-line chaining is performed.")
    R.notes.append("sd_field is set by the fitted nugget/sill partition, which short-lag pairs "
           "do not identify on a control set spread over ~200 km; its sweep range is "
           "reported for that reason.")
    R.banner()
    return R


HEADER = ["epoch", "selection", "n_marks", "constant_mm", "sweep_min_mm",
          "sweep_max_mm", "sd_field_mm", "sd_field_min_mm", "sd_field_max_mm"]


def row(est):
    sel = est.cover_treatment if est.surface is None else f"{est.surface}/{est.cover_treatment}"
    return [est.epoch, sel, est.n_marks,
            f"{est.constant_mm:+.2f}", f"{est.constant_min_mm:+.2f}",
            f"{est.constant_max_mm:+.2f}", f"{est.sd_field_mm:.2f}",
            f"{est.sd_field_min_mm:.2f}", f"{est.sd_field_max_mm:.2f}"]


def main(argv=None):
    a = parse_args(argv)
    R = banner(a)
    kw = dict(easting=a.easting, northing=a.northing, max_lags_m=a.max_lags_m,
              n_lags=a.n_lags, n_pairs=a.n_pairs, estimators=a.estimators, seed=a.seed)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    products, rows = {}, []

    g1 = datum.sweep_treatments("gen1", bridge_mm=a.gen1_bridge_mm,
                                bridge_sd_mm=a.gen1_bridge_sd_mm, **kw)
    for t_, e in g1.items():
        products[f"gen1__{t_}"] = e
        rows.append(row(e))
    for s in a.gen2_surfaces:
        g2 = datum.sweep_treatments("gen2", surface=s, bridge_mm=a.gen2_bridge_mm,
                                    bridge_sd_mm=a.gen2_bridge_sd_mm, **kw)
        for t_, e in g2.items():
            products[f"gen2__{s}__{t_}"] = e
            rows.append(row(e))
    R.table(HEADER, rows)

    v1 = [e.constant_mm for e in g1.values()]
    v2 = [e.constant_mm for k, e in products.items() if k.startswith("gen2")]
    print()
    print(f"  cover treatment alone moves gen1 by {np.ptp(v1):.2f} mm "
          f"({min(v1):+.2f} .. {max(v1):+.2f}) -- not chosen here")
    print(f"  surface x cover moves gen2 by {np.ptp(v2):.2f} mm "
          f"({min(v2):+.2f} .. {max(v2):+.2f}) -- not chosen here")

    g1o = products["gen1__open"]
    print()
    print("  OPEN-GROUND, the treatment that leaves forest-open in the DoD rather than")
    print("  absorbing it into the datum:")
    for s in a.gen2_surfaces:
        e = products[f"gen2__{s}__open"]
        d = g1o.constant_mm - e.constant_mm
        sd = float(np.hypot(g1o.sd_field_mm, e.sd_field_mm))
        print(f"    gen1 {g1o.constant_mm:+.2f} - gen2({s}) {e.constant_mm:+.2f}"
              f"  =  DoD absolute correction {d:+.2f} +/- {sd:.2f} mm")

    for name, e in products.items():
        (out / f"{a.site}__{name}__datum.json").write_text(json.dumps(e.to_dict(), indent=2))
    print(f"\n  wrote {len(products)} products to {out}/")
    R.done(headline=(f"gen1 open {g1o.constant_mm:+.2f} mm; "
                     + "; ".join(f"gen2 {s} open "
                                 f"{products[f'gen2__{s}__open'].constant_mm:+.2f} mm"
                                 for s in a.gen2_surfaces)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
