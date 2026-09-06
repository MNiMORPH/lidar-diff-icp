"""Where does the vendor's 2008 class-2 GROUND end, across the flight line?

This is the measurement that explains why widening the control catchment cannot add
marks.  The delivered 2008 tiles carry an explicit **overlap class (12)**: in the region
where two adjacent swaths see the same ground, the vendor kept ONE line's returns in the
ground class and moved the other line's into class 12.  So the vendor's bare-earth ground
is cut at a SEAM at roughly half the line spacing, not at the edge of the swath.

Printed per line, from the returns themselves:
  * ``all``     -- every return of that line: the true swath half-width
  * ``class12`` -- the overlap class: the outer band the vendor discarded from ground
  * ``class2``  -- the bare-earth ground: the band that survives, i.e. the seam

Cross-track distance is measured from the whole-line straight ground track fitted by
``gen1_line_tracks.py``.  Percentiles are p0.05/p99.95 rather than min/max so a handful
of stray points do not set the number.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

SCRATCH = os.environ.get("SCRATCH", ".")
LINES = (133, 134, 135, 136, 137, 138)
# tiles chosen to span the six lines along their whole length; any tile holding the line
# would do, and the extremes below are taken over all of them
TILES = ["4342-29-64", "4342-30-64", "4342-29-63", "4342-30-63", "4342-28-63",
         "4342-28-64", "5142-15-63", "4342-20-63"]
OVERLAP_CLASS = 12
GROUND_CLASS = 2


def main():
    import laspy
    tracks = json.load(open(f"{SCRATCH}/line_tracks.json"))["whole"]

    def cross(ln, x, y):
        t = tracks[str(ln)]
        return (x - (np.polyval(t["w"], y - t["ym"]) + t["xm"])) / math.hypot(1.0, t["w"][0])

    acc = {}
    for name in TILES:
        p = f"data/before/{name}.laz"
        if not os.path.exists(p):
            print(f"  {name}: not on disk, skipped"); continue
        f = laspy.read(p)
        cl = np.asarray(f.classification)
        s = np.asarray(f.point_source_id).astype(int)
        x, y = np.asarray(f.x), np.asarray(f.y)
        for ln in LINES:
            m = s == ln
            if not m.any():
                continue
            for lbl, mm in (("all", m), (f"class{OVERLAP_CLASS}", m & (cl == OVERLAP_CLASS)),
                            (f"class{GROUND_CLASS}", m & (cl == GROUND_CLASS))):
                if not mm.any():
                    continue
                c = cross(ln, x[mm], y[mm])
                acc.setdefault((ln, lbl), []).append(
                    (np.percentile(c, 0.05), np.percentile(c, 99.95), c.size))
        del f
        print(f"  read {name}", flush=True)

    print(f"\n{'line':>5} {'subset':>8} {'p0.05_m':>9} {'p99.95_m':>9} {'half_m':>8} {'n':>13}")
    for k in sorted(acc):
        A = np.array(acc[k])
        lo, hi = A[:, 0].min(), A[:, 1].max()
        print(f"{k[0]:>5} {k[1]:>8} {lo:>9.1f} {hi:>9.1f} {max(abs(lo), abs(hi)):>8.1f} "
              f"{int(A[:, 2].sum()):>13,}")

    # line spacing, at the Elba reference northing
    NREF = 4883677.71
    e = [np.polyval(tracks[str(l)]["w"], NREF - tracks[str(l)]["ym"]) + tracks[str(l)]["xm"]
         for l in LINES]
    d = np.diff(e)
    print(f"\nline spacing at N = {NREF:.0f}: " + ", ".join(f"{v:.0f}" for v in d) +
          f" m; mean {d.mean():.0f} m, half-spacing {d.mean()/2:.0f} m")


if __name__ == "__main__":
    main()
