#!/usr/bin/env python3
"""Compare the continuity-derived gen1 raise against the control-derived one, and state
the erosion rate the difference between them implies.

Sediment continuity constrains ``delta`` only through the scene budget, and the budget
and a uniform raise trade off exactly one for one:

    scene net volume(delta) = sum over the evaluated cells of (DoD - delta) * area .

So any claim that ``delta`` is larger than the continuity match point is the same claim
that the hillslopes lost that much sediment over the epoch. This prints that equivalence
in mm/yr, so the two instruments can be compared in the units geomorphology is measured
in rather than argued about in the abstract.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/mass_balance/compare_to_control.py --delta-mb <mm> --se-mb <mm>
"""
import argparse
import datetime as dt

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--delta-mb", type=float, required=True, help="continuity delta, mm")
ap.add_argument("--se-mb", type=float, nargs="+", required=True,
                help="the 1-sigma COMPONENTS of that delta, mm; combined in quadrature "
                     "here so the total is derived rather than written down")
ap.add_argument("--se-name", nargs="+", default=None, help="what each component is")
ap.add_argument("--label", default="")
ap.add_argument("--delta-ctrl", type=float, default=53.6,
                help="control-derived raise, mm (analysis/FRAME_2026-08-26-PM.md)")
ap.add_argument("--se-ctrl", type=float, default=13.0)
ap.add_argument("--gen1-date", default="2008-11-25")
ap.add_argument("--gen2-date", default="2021-05-01")
A = ap.parse_args()

yrs = ((dt.date.fromisoformat(A.gen2_date) - dt.date.fromisoformat(A.gen1_date)).days
       / 365.2425)
se_mb = float(np.sqrt(np.sum(np.square(A.se_mb))))
d = A.delta_ctrl - A.delta_mb
sd = float(np.hypot(A.se_ctrl, se_mb))
print(f"{A.label}")
for i, c in enumerate(A.se_mb):
    nm = A.se_name[i] if A.se_name and i < len(A.se_name) else f"component {i+1}"
    print(f"    sigma component: {nm} = {c:.2f} mm")
print(f"  continuity (mass balance)      {A.delta_mb:+.1f} +/- {se_mb:.1f} mm "
      f"(quadrature of the above)")
print(f"  control marks                  {A.delta_ctrl:+.1f} +/- {A.se_ctrl:.1f} mm")
print(f"  difference                     {d:+.1f} +/- {sd:.1f} mm  ->  {d/sd:.2f} sigma")
print(f"  epoch {A.gen1_date} -> {A.gen2_date} = {yrs:.2f} yr")
print(f"  net hillslope lowering the control value implies, on top of the observed "
      f"budget: {d:.1f} mm over the epoch = {d/yrs:.2f} mm/yr")
