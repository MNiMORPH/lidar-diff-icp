#!/usr/bin/env python3
"""Filter a point cloud to LAST RETURNS, streaming and memory-light.

Bare earth (the standard mapmaker definition) = the last return of each pulse:
``return_number == number_of_returns``. This INCLUDES single-return pulses
(for a single return, return_number == number_of_returns == 1), which dominate
flat, open ground (bare soil, short crop) where a pulse reflects once.

Caution -- the trap this script exists to avoid: PDAL's ``filters.returns``
``groups=last`` keeps only the last of MULTI-return pulses and silently drops
single returns; use ``last,only``. An earlier build here used the multi-only
form and emptied the agricultural fields (single returns), leaving 72.7% cell
coverage instead of ~99.7%. Here the mask is explicit: rn == nr.

Streams in chunks so a 180 M-point file never fully lands in memory.

    python scripts/filter_last_return.py \
        data/after/3dep2021_fulldensity.laz data/after/3dep2021_last.laz
"""
import argparse
from pathlib import Path
import numpy as np
import laspy


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--chunk", type=int, default=5_000_000)
    a = ap.parse_args()

    Path(a.dst).parent.mkdir(parents=True, exist_ok=True)
    kept = tot = 0
    with laspy.open(a.src) as fh:
        writer = laspy.open(a.dst, mode="w", header=fh.header)
        try:
            for pts in fh.chunk_iterator(a.chunk):
                rn = np.asarray(pts.return_number)
                nr = np.asarray(pts.number_of_returns)
                m = rn == nr                       # last return, INCLUDING singles
                tot += rn.size; kept += int(m.sum())
                writer.write_points(pts[m])
        finally:
            writer.close()
    print(f"{a.src}: {tot:,} pts -> {kept:,} last returns ({100*kept/tot:.1f}%)")
    print(f"wrote {a.dst}")


if __name__ == "__main__":
    main()
