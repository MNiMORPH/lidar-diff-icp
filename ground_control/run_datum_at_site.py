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
    print("=" * 78)
    print("DATUM AT SITE -- each epoch from its OWN contemporaneous control")
    print("=" * 78)
    print(f"  site              {a.site}  E {a.easting}  N {a.northing}")
    print(f"  gen1 control      2008 MnGeo/MnDNR validation, {control.GEN1_CSV.name}")
    print(f"  gen2 control      2021 USGS NVA/VVA checkpoints, {control.GEN2_CSV.name}")
    print(f"  NO cross-epoch control is read.  NO chaining between flight lines.")
    print(f"  variogram sweep   max_lag_m {a.max_lags_m}  x  estimator {a.estimators}")
    print(f"                    n_lags {a.n_lags}  n_pairs {a.n_pairs}  seed {a.seed}")
    print(f"                    [CALLER-supplied; there are no defaults, because the")
    print(f"                     control does not determine a variogram range]")
    print(f"  sign convention   tie = surveyed - z_lidar; POSITIVE = surface reads LOW")
    print()


def show(tag, est):
    print(f"  {tag:34s} n={est.n_marks:4d}  "
          f"constant {est.constant_mm:+8.2f} mm  "
          f"[sweep {est.constant_min_mm:+8.2f} .. {est.constant_max_mm:+8.2f}]  "
          f"sd_field {est.sd_field_mm:6.2f} [{est.sd_field_min_mm:.2f}..{est.sd_field_max_mm:.2f}]")


def main(argv=None):
    a = parse_args(argv)
    banner(a)
    kw = dict(easting=a.easting, northing=a.northing, max_lags_m=a.max_lags_m,
              n_lags=a.n_lags, n_pairs=a.n_pairs, estimators=a.estimators, seed=a.seed)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    products = {}

    print("GEN1  (2008 surface, 2008 control)")
    g1 = datum.sweep_treatments("gen1", bridge_mm=a.gen1_bridge_mm,
                                bridge_sd_mm=a.gen1_bridge_sd_mm, **kw)
    for t, e in g1.items():
        show(f"cover treatment: {t}", e)
        products[f"gen1__{t}"] = e
    v = [e.constant_mm for e in g1.values()]
    print(f"  --> the COVER TREATMENT alone moves gen1 by {np.ptp(v):.2f} mm "
          f"({min(v):+.2f} .. {max(v):+.2f}); this file does not choose one")

    print("\nGEN2  (2021 surface, 2021 control)")
    allv = []
    for s in a.gen2_surfaces:
        g2 = datum.sweep_treatments("gen2", surface=s, bridge_mm=a.gen2_bridge_mm,
                                    bridge_sd_mm=a.gen2_bridge_sd_mm, **kw)
        for t, e in g2.items():
            show(f"surface {s} / cover {t}", e)
            products[f"gen2__{s}__{t}"] = e
            allv.append(e.constant_mm)
    print(f"  --> surface x cover together move gen2 by {np.ptp(allv):.2f} mm "
          f"({min(allv):+.2f} .. {max(allv):+.2f})")

    print("\nOPEN-GROUND HEADLINE (the treatment that does not foreclose the")
    print("canopy-vs-erosion question, because it leaves forest-open in the DoD):")
    g1o = products["gen1__open"]
    print(f"  gen1  ADD {g1o.constant_mm:+8.2f} mm   sd_field {g1o.sd_field_mm:.2f} mm")
    for s in a.gen2_surfaces:
        e = products[f"gen2__{s}__open"]
        print(f"  gen2 {s:8s} ADD {e.constant_mm:+8.2f} mm   sd_field {e.sd_field_mm:.2f} mm")
        d = g1o.constant_mm - e.constant_mm
        sd = float(np.hypot(g1o.sd_field_mm, e.sd_field_mm))
        print(f"       implied DoD absolute correction gen1-gen2 = {d:+.2f} +/- {sd:.2f} mm")

    for name, e in products.items():
        p = out / f"{a.site}__{name}__datum.json"
        p.write_text(json.dumps(e.to_dict(), indent=2))
    print(f"\nwrote {len(products)} products to {out}/")
    if g1o.notes:
        print("\nNOTES")
        for n in g1o.notes:
            print(f"  * {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
