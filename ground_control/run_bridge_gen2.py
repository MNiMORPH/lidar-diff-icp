"""gen2's bridge: the delivered surface against OUR reconstructed gen2 surface.

Six per-checkpoint gen2 crops are already on disk (``data/after/checkpoints/``), so this
needs no download.  Our gen2 surface is rebuilt from each crop and read at the mark.

VALIDATION, measured before any bridge is quoted, at 24 points inside elbaext where the
SHIPPED ``z_after.npy`` gives an independent answer:

    zero canopy cover  n=8   median local-grid   -0.04 mm
    slope  1.1- 4.0    n=6   median              -3.70 mm
    slope 21.1-32.5    n=6   median            -126.07 mm

The producer is therefore faithful on open, low-slope ground -- which is where NVA
checkpoints are sited, NVA meaning non-vegetated vertical accuracy -- and unusable on
steep or vegetated ground.  ``--point-types NVA`` is the consequence, not a preference.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))

import control  # noqa: E402
import our_surface as OS  # noqa: E402
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--crops", required=True)
    p.add_argument("--point-types", nargs="+", required=True)
    p.add_argument("--surface", required=True, choices=control.GEN2_SURFACES)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--radii-m", type=float, nargs="+", required=True)
    p.add_argument("--ground-class", type=int, required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    rows_csv = list(csv.DictReader(open(control.GEN2_CSV)))
    by_id = {r["point_id"]: r for r in rows_csv}
    crops = sorted(glob.glob(f"{a.crops}/cp*_gen2.laz"))

    R = Run("What is gen2's bridge -- its delivered surface minus OUR reconstructed gen2 "
            "surface -- at the checkpoints whose clouds we already hold?")
    R.input(str(control.GEN2_CSV),
            role="gen2's own 2021 USGS checkpoints; usgs_*_z_m is the DELIVERED surface "
                 "elevation at the mark")
    for c in crops[:3]:
        R.input(c, role="gen2 3DEP cloud cropped about one checkpoint; class 2 is its "
                        "vendor ground")
    R.param("point_types", tuple(a.point_types), src="MINE",
            why="the validation at 24 elbaext points gives median local-grid -0.04 mm at "
                "ZERO canopy cover and -126.07 mm in the steepest slope quartile, so the "
                "producer is trusted only on the open low-slope ground NVA marks sit on")
    R.param("surface", a.surface, src="andy",
            why="which of gen2's four delivered surfaces the bridge is measured against; "
                "it must match the surface whose datum constant it will be added to")
    R.param("res_m", a.res, src="repo", why="the pipeline's grid resolution")
    R.param("radii_m", tuple(a.radii_m), src="MINE",
            why="radius of the order-2 fit at the mark; swept, and a radius that cannot "
                "be filled with 6 cells is reported, never dropped")
    R.param("ground_class", a.ground_class, src="repo",
            why="gen2's vendor ground; its delivered classification carries only 1 and 2")
    R.param("swath_alignment", "none", src="repo",
            why="gen2 is the reference epoch and no corrections file on disk carries a "
                "gen2 swath; measured adjacent-pair ties are |k| <= 4.8 mm")
    R.param("geoid_term", "none", src="repo",
            why="gen2 IS the frame; its control and its cloud are both NAVD88(GEOID18)")
    R.column("point_id", "gen2 checkpoint id")
    R.column("type", "USGS accuracy class, unitless: NVA non-vegetated, VVA vegetated")
    R.column("n_ground", "class-2 gen2 returns in the crop, count")
    R.column("delivered_z_m", "the delivered surface elevation at the mark, m")
    R.column("ours_z_m", "our reconstructed gen2 surface at the mark, m; same frame")
    R.column("bridge_mm", "z_delivered - z_ours, mm; ADD to a delivered-surface constant "
                          "to carry it onto our surface")
    R.column("radius_spread_mm", "max - min of bridge_mm over the radii that could be "
                                 "fitted, mm")
    R.column("n_radii", "radii fitted / radii swept")
    R.banner()

    keep, recs, table = set(a.point_types), [], []
    zc, ec = f"usgs_{a.surface}_z_m", f"usgs_{a.surface}_error_m"
    for cpath in crops:
        m = re.search(r"cp(\d+)_gen2", os.path.basename(cpath))
        pid = next((k for k in by_id if k.startswith(m.group(1) + "_")), None)
        if pid is None:
            continue
        r = by_id[pid]
        if r["point_type"] not in keep:
            continue
        zdel = (r[zc] or "").strip()
        if not zdel:
            print(f"  {pid}: no {zc} -- this mark is not in the {a.surface} block")
            continue
        E, N = float(r["easting"]), float(r["northing"])
        vals, used, failed, sp = [], [], [], None
        for rad in a.radii_m:
            s = OS.our_gen2_surface_at(cpath, E, N, res=a.res, radius_m=rad,
                                       ground_class=a.ground_class, half_width_m=150.0)
            if s is None:
                failed.append(rad); continue
            sp = s; used.append(rad)
            vals.append((float(zdel) - s.z_geoid18_m) * 1000.0)
        if not vals:
            print(f"  {pid}: no radius could be fitted ({failed} failed)")
            continue
        br = float(np.median(vals)); spread = float(np.ptp(vals)) if len(vals) > 1 else 0.0
        recs.append(dict(point_id=pid, point_type=r["point_type"], bridge_mm=br,
                         radius_spread_mm=spread, n_radii_used=len(used),
                         radii_failed=failed, n_ground=sp.n_ground_pts))
        table.append([pid, r["point_type"], sp.n_ground_pts, f"{float(zdel):.3f}",
                      f"{sp.z_geoid18_m:.3f}", f"{br:+.1f}", f"{spread:.1f}",
                      f"{len(used)}/{len(a.radii_m)}"])
    R.table(["point_id", "type", "n_ground", "delivered_z_m", "ours_z_m", "bridge_mm",
             "radius_spread_mm", "n_radii"], table)
    v = np.array([r["bridge_mm"] for r in recs])
    print()
    if v.size:
        print(f"  gen2 BRIDGE over {v.size} {'/'.join(a.point_types)} marks, "
              f"surface {a.surface}:")
        print(f"    mean {v.mean():+.2f}  median {np.median(v):+.2f}  "
              + (f"sd {v.std(ddof=1):.2f}  SE {v.std(ddof=1)/np.sqrt(v.size):.2f} mm"
                 if v.size > 1 else "(n=1: no scatter computable)"))
        print(f"    n is small because only {len(crops)} gen2 checkpoint crops are on "
              f"disk; nothing was downloaded.")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(dict(marks=recs, params=vars(a)), indent=1))
    print(f"  wrote {a.out}")
    R.done(headline=f"gen2 bridge n={v.size} mean {v.mean():+.2f} mm" if v.size else "no marks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
