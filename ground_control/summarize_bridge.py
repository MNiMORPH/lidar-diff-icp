"""Summarise a committed bridge product, through trust/provenance.py.

The bridge table is produced by ``run_bridge_wide.py``, which carries its own banner.
This reads that committed product back and prints the population statistics, so the
summary numbers quoted in prose come from a provenance-carrying run rather than from an
ad-hoc script.

    ./lidar-icp/bin/python ground_control/summarize_bridge.py \
        --product ground_control/products/bridge_wide_L1O.json
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
    p.add_argument("--product", required=True)
    a = p.parse_args(argv)
    d = json.loads(Path(a.product).read_text())

    R = Run("What is the bridge on open ground, over all marks, and is any single "
            "mark's value resolvable?")
    R.input(a.product, role="per-mark bridge table from run_bridge_wide.py; bridge_mm is "
                            "z_delivered - z_ours, both on NAVD88(GEOID03)")
    R.param("covers", tuple(d["covers"]), src="andy",
            why="open ground only; the validation gate reproduces the shipped grid to "
                "+1.7 mm on the open mark and -30.5 mm on the vegetated one")
    R.param("estimator", "mean and median, both reported, neither chosen", src="MINE",
            why="the per-mark values are heavy-tailed, so a mean alone would be a claim "
                "about the tails; nothing is trimmed or filtered")
    R.column("stat", "name of the statistic")
    R.column("value_mm", "its value in mm, or a count where the name says count")
    R.banner()

    v = np.array([m["bridge_mm"] for m in d["marks"]])
    s = np.array([m["radius_spread_mm"] for m in d["marks"]])
    rows = [
        ["n marks (count)", f"{v.size}"],
        ["mean", f"{v.mean():+.2f}"],
        ["SE of the mean", f"{v.std(ddof=1)/np.sqrt(v.size):.2f}"],
        ["median", f"{np.median(v):+.2f}"],
        ["NMAD about the median", f"{1.4826*np.median(np.abs(v-np.median(v))):.2f}"],
        ["sd", f"{v.std(ddof=1):.2f}"],
        ["minimum", f"{v.min():+.1f}"],
        ["maximum", f"{v.max():+.1f}"],
        ["per-mark radius spread, median", f"{np.median(s):.1f}"],
        ["per-mark radius spread, max", f"{s.max():.1f}"],
        ["marks with |bridge| < own radius spread (count)", f"{int((np.abs(v) < s).sum())}"],
    ]
    R.table(["stat", "value_mm"], rows)
    R.done(headline=f"bridge {v.mean():+.2f} +/- {v.std(ddof=1)/np.sqrt(v.size):.2f} mm "
                    f"over {v.size} open marks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
