#!/usr/bin/env python3
"""DeLong's diagnostic, done densely: does FLIGHT-LINE structure survive our corrections?

Andy, 2026-09-11: "Give the DeLong method a try."

WHY DENSE, AND NOT AT THE MARKS. The 14-mark test (`bridge_by_line.py`) returned a null --
F(4,9)=0.308, ICC(1)=-0.328 -- but it could not have found the effect it was looking for.
Within-line scatter at the marks is 99.2 mm, so one line mean at n=3 carries SE 57.2 mm and
the 95% upper bound on a true between-line sd is 56.1 mm; meanwhile the per-line dz our own
`align_swaths` applies has median sd 28.3 mm across 60 cached tiles. The real effect sits at
HALF the detection ceiling. Going to all 88 marks does not fix it: they touch 38 lines, 2.3
marks per line.

DeLong et al. (2022) faced exactly this and did NOT use checkpoints. They characterised the
spurious field from the dense stable-terrain DoD -- over 1.337e9 pixels, mean 0.002 m,
sd 0.103 m -- and never tied to a mark at all. Marks are for the ABSOLUTE LEVEL, which no
amount of dense data can supply; dense stable ground is for SHAPE and per-line structure,
which no realistic number of marks can resolve. This script is the second half of that.

WHAT IT MEASURES. `beam_offset_table.parquet` carries, per return, the raw offset of gen1
CSF ground to the gen2 surface (`d_mm`), the fully corrected offset (`d_mm_corr`), and the
four registration terms separately -- so the SAME returns can be grouped by
`point_source_id` before and after our correction chain. That is DeLong's "are there
flightline-parallel streaks left in the DoD" question, asked quantitatively.

DELONG'S STABLE MASK, used as they defined it: slope <= 3 deg and |dz| <= 0.7 m. They also
buffer streams by 100 m because real low-slope floodplain deposition would otherwise be
absorbed into "stable"; we substitute the pipeline's own `floodplain_mask.npy`, which cuts
the floodplain by ELEVATION (flow routing being unreliable on flats). Each mask's kept count
is printed.

NO SE IS QUOTED PER LINE. With millions of returns a per-return standard error would be
meaningless -- returns within a flight line are strongly spatially correlated, so it would
understate by orders of magnitude. The between-line SPREAD of the line means is reported
instead, which is the quantity the mark-based test could not pin down.

    ./lidar-icp/bin/python ground_control/dense_line_structure.py --tile elba_fulldensity
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from trust.provenance import Run

COLS = ["point_source_id", "d_mm", "d_mm_corr", "slope", "in_grid", "cell",
        "dz_swath_mm", "dz_drift_mm", "dz_geoid_mm", "dz_lateral_mm"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", required=True, help="a data/derived/<tile> holding the table")
    ap.add_argument("--slope-max-deg", type=float, default=3.0)
    ap.add_argument("--dz-max-mm", type=float, default=700.0)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    d = f"data/derived/{a.tile}"
    tab = os.path.join(d, "beam_offset_table.parquet")
    fp = os.path.join(d, "floodplain_mask.npy")

    R = Run("After our full correction chain, does the gen1-to-gen2 offset still organise "
            "by FLIGHT LINE on stable ground -- DeLong's flightline-streak question, asked "
            "with millions of returns instead of 14 marks?")
    R.input(tab, role="per-return gen1 CSF ground offset to the gen2 surface, with the four "
                      "registration terms and point_source_id per return; supplies BOTH the "
                      "raw d_mm and the corrected d_mm_corr for the SAME returns")
    R.input(fp, role="the pipeline's floodplain mask, cut by ELEVATION; stands in for "
                     "DeLong's 100 m stream buffer")
    R.param("slope_max_deg", a.slope_max_deg, src="MINE",
            why="DeLong et al. (2022) s4.4: they mask likely-real change by 'areas with "
                "slope >3 deg' before interpolating the stable residual. VERIFIED verbatim "
                "in the extracted text, not relayed. Marked MINE, not 'repo', because a "
                "published constant has no category in this harness and adopting it was my "
                "choice -- but the VALUE is checkable in the paper rather than invented. "
                "Its effect on the kept count is printed as its own mask")
    R.param("dz_max_mm", a.dz_max_mm, src="MINE",
            why="DeLong's companion cut, same sentence: 'DoD values >0.7 m in flat areas'. "
                "Same verification and same reasoning as slope_max_deg. Note it is applied "
                "to the RAW d_mm, so the same returns are kept for the before and after "
                "columns and the comparison is like with like")
    R.param("floodplain", "pipeline floodplain_mask.npy (elevation-cut)", src="MINE",
            why="DeLong buffer streams by 100 m so real low-slope deposition is not absorbed "
                "into 'stable'. We substitute the pipeline's own elevation-cut floodplain "
                "mask, because flow routing is unreliable on flats. A SUBSTITUTION for their "
                "buffer, not a reproduction of it")
    R.param("tile", a.tile, src="MINE",
            why="the tile carrying a built beam_offset_table; chosen for availability, not "
                "by any property of its answer")
    R.column("point_source_id", "gen1 flight line, unitless")
    R.column("n", "stable returns on that line, count")
    R.column("raw_mm", "mean d_mm -- gen1 CSF ground minus the gen2 surface BEFORE our "
                       "corrections, mm")
    R.column("corr_mm", "mean d_mm_corr -- the SAME returns AFTER geoid, lateral, swath and "
                        "drift, mm")
    R.column("swath_mm", "mean VERTICAL EFFECT of the swath transform on that line's "
                         "returns, mm. NOT a per-line constant: align_swaths solves "
                         "(dx, dy, dz), and the horizontal part drags each return across "
                         "sloped ground, so this varies return to return (on line 136 it "
                         "takes 173009 distinct values, sd 62.1 mm). Its MEAN therefore "
                         "depends on which returns are selected")
    R.column("drift_mm", "mean of the per-swath along-track drift spline, mm")
    R.column("between_line_sd", "sd over the LINE MEANS, ddof=1, mm. No per-return SE is "
                                "quoted: returns within a line are strongly spatially "
                                "correlated, so one would understate by orders of magnitude")

    t = pd.read_parquet(tab, columns=COLS)
    n_all = len(t)
    fpm = np.load(fp)
    ny, nx = fpm.shape
    iy, ix = np.divmod(t["cell"].to_numpy(), nx)
    ok = (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx)
    in_fp = np.zeros(n_all, bool)
    in_fp[ok] = fpm[iy[ok], ix[ok]].astype(bool)

    m_grid = t["in_grid"].to_numpy().astype(bool)
    m_slope = t["slope"].to_numpy() <= a.slope_max_deg
    m_dz = np.abs(t["d_mm"].to_numpy()) <= a.dz_max_mm
    R.mask("in_grid", m_grid, defn="return falls inside the product grid", of=n_all)
    R.mask("slope_le_3deg", m_grid & m_slope,
           defn=f"DeLong: local slope <= {a.slope_max_deg} deg", of=n_all)
    R.mask("dz_le_0p7m", m_grid & m_slope & m_dz,
           defn=f"DeLong: |raw dz| <= {a.dz_max_mm} mm, on top of the slope cut", of=n_all)
    stable = m_grid & m_slope & m_dz & ~in_fp
    R.mask("stable", stable,
           defn="DeLong stable ground: in grid, slope <= 3 deg, |raw dz| <= 0.7 m, and NOT "
                "floodplain (elevation-cut, standing in for their 100 m stream buffer)",
           of=n_all)

    s = t[stable]
    g = s.groupby("point_source_id").agg(
        n=("d_mm", "size"), raw_mm=("d_mm", "mean"), corr_mm=("d_mm_corr", "mean"),
        swath_mm=("dz_swath_mm", "mean"), drift_mm=("dz_drift_mm", "mean"),
        geoid_mm=("dz_geoid_mm", "mean"), lateral_mm=("dz_lateral_mm", "mean"))
    print(f"\n{'psid':>6}{'n':>11}{'raw_mm':>10}{'corr_mm':>10}{'swath_mm':>10}"
          f"{'drift_mm':>10}{'geoid_mm':>10}{'lateral_mm':>11}")
    for ln, r in g.iterrows():
        print(f"{ln:>6}{int(r['n']):>11,}{r['raw_mm']:>10.1f}{r['corr_mm']:>10.1f}"
              f"{r['swath_mm']:>10.1f}{r['drift_mm']:>10.1f}{r['geoid_mm']:>10.1f}"
              f"{r['lateral_mm']:>11.1f}")

    # *** SPATIAL CONFOUND -- report it before any per-line claim. ***
    # If the flight lines do not overlap spatially, "between-line" and "a tilt across the
    # tile" are the SAME quantity and this script cannot tell them apart. Measured, never
    # assumed: where each line's stable returns actually sit.
    fy, fx = np.divmod(s["cell"].to_numpy(), nx)
    print(f"\n  SPATIAL SPREAD OF EACH LINE (grid {ny} rows x {nx} cols):")
    print(f"  {'psid':>6}{'col_mean':>10}{'col_p5':>8}{'col_p95':>8}   <- if these barely "
          f"overlap, line number IS a spatial coordinate")
    spans = {}
    for ln in g.index:
        m = s["point_source_id"].to_numpy() == ln
        lo, hi = np.percentile(fx[m], 5), np.percentile(fx[m], 95)
        spans[ln] = (lo, hi)
        print(f"  {ln:>6}{fx[m].mean():>10.1f}{lo:>8.0f}{hi:>8.0f}")
    order = sorted(spans)
    ov = [(a_, b_, max(0.0, min(spans[a_][1], spans[b_][1]) - max(spans[a_][0], spans[b_][0])))
          for a_, b_ in zip(order, order[1:])]
    print("  adjacent-pair 5-95%% column overlap: "
          + ", ".join(f"{a_}-{b_}: {w:.0f} cols" for a_, b_, w in ov))

    sd_raw = float(g["raw_mm"].std(ddof=1))
    sd_corr = float(g["corr_mm"].std(ddof=1))
    rng_raw = float(g["raw_mm"].max() - g["raw_mm"].min())
    rng_corr = float(g["corr_mm"].max() - g["corr_mm"].min())
    print(f"\n  stable returns: {int(stable.sum()):,} of {n_all:,} "
          f"({100*stable.sum()/n_all:.1f}%), on {len(g)} flight lines")
    print(f"  BETWEEN-LINE sd of the line means:   raw {sd_raw:7.2f}  ->  "
          f"corrected {sd_corr:7.2f} mm")
    print(f"  BETWEEN-LINE range of the line means: raw {rng_raw:7.2f}  ->  "
          f"corrected {rng_corr:7.2f} mm")
    print(f"  within-line sd of single returns (pooled): "
          f"{float(s['d_mm_corr'].std(ddof=1)):.1f} mm")
    print(f"\n  for scale, the mark-based test could not exclude a between-line sd below "
          f"56.1 mm;\n  here it is measured at {sd_corr:.2f} mm after correction.")
    print("\n  READ THIS BEFORE CALLING IT PER-LINE STRUCTURE: the number above is only a\n"
          "  FLIGHT-LINE effect if the lines overlap spatially. Where they do not, the same\n"
          "  number is equally a smooth tilt across the tile, and this grouping cannot\n"
          "  separate them. Resolve it in the OVERLAP cells, where two lines measure the\n"
          "  SAME ground -- a step there is per-line; agreement there means it is spatial.")

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(dict(tile=a.tile, n_all=n_all, n_stable=int(stable.sum()),
                           between_line_sd_raw_mm=sd_raw,
                           between_line_sd_corr_mm=sd_corr,
                           between_line_range_raw_mm=rng_raw,
                           between_line_range_corr_mm=rng_corr,
                           by_line=g.reset_index().to_dict("records")), fh, indent=1)
            fh.write("\n")
        print(f"    wrote {a.out}")

    R.done(headline=(f"{int(stable.sum()):,} stable returns on {len(g)} lines; between-line "
                     f"sd of the line means {sd_raw:.2f} -> {sd_corr:.2f} mm after our "
                     f"corrections (mark-based ceiling was 56.1 mm)"))


if __name__ == "__main__":
    main()
