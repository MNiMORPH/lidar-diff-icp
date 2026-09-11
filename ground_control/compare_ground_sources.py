#!/usr/bin/env python3
"""Do the marks prefer the VENDOR's ground surface, our CSF surface, or our class-2 surface?

Andy, 2026-09-10: "run that comparison. vendor and CSF, both on vendor's elevation."

THE THREE SURFACES, all evaluated at the same marks:

    vendor   the MnDNR validation report's own `Surface Z` -- the delivered ASPRS class-2
             ground, interpolated at the mark by the DNR, not by us
    csf      our reconstruction with `ground_source="csf"` (cloth simulation filter)
    class2   our reconstruction with `ground_source="class2"` (the vendor's own class 2,
             re-gridded by US)

`run_bridge_wide.py` was run TWICE, changing ONLY `--ground-source`. Everything else --
swath re-solve, the absent lateral Nuth-Kaaeb term, 5 m gridding, the order-2 fit, the
window -- is byte-identical between the arms, so

    z_class2 - z_csf  =  bridge_csf - bridge_class2

is the ground-source contrast with NOTHING else in it. That isolation is the point: an
earlier comparison of "ours" against the vendor bundled FIVE differences at once and could
not attribute any of them.

OFFSETS AGAINST THE MARKS. `dnr_error_m` is the report's own Control Z - Surface Z, so

    offset_vendor = dnr_error_m                       (surveyed - z_delivered)
    offset_ours   = dnr_error_m + bridge              (surveyed - z_ours)

because bridge = z_delivered - z_ours. Positive = the surface reads LOW at the mark.

WHAT THE COMPARISON TURNS ON -- the SPREAD, not the mean. A constant offset common to all
marks is a datum constant and is absorbed downstream (that is what the whole gen1 datum
thread is for). Scatter about it is not absorbable. So the surface that tracks the marks
better is the one with the smaller sd, even if its mean sits further from zero.

SPLIT BY COVER, and that is the whole reason this run exists. The earlier comparison was
L1O-only -- open ground, exactly where a vendor classifier has the easiest job -- so it was
structurally incapable of seeing a classifier difference. Cover is a GROUPING here, never a
filter: every mark in both arms is reported, and the "other" class is shown too.

NO OUTLIER RULE IS APPLIED. Skew and the extremes are printed so the choice of whether to
cut, and where, stays Andy's.

    ./lidar-icp/bin/python ground_control/compare_ground_sources.py \
        --csf <cmp_csf.json> --class2 <cmp_class2.json> --project lidar_semn2008
"""
import argparse
import json

import numpy as np

from lidar_diff_icp.groundtruth import gen1_datum as G
from trust.provenance import Run


def stats(v):
    v = np.asarray(v, float)
    if v.size == 0:
        return dict(n=0)
    m = v.mean()
    sd = v.std(ddof=1) if v.size > 1 else float("nan")
    return dict(n=v.size, mean=m, median=float(np.median(v)), sd=sd,
                se=sd / np.sqrt(v.size) if v.size > 1 else float("nan"),
                nmad=1.4826 * float(np.median(np.abs(v - np.median(v)))),
                skew=float(((v - m) ** 3).mean() / sd ** 3) if v.size > 2 and sd > 0
                else float("nan"),
                lo=float(v.min()), hi=float(v.max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csf", required=True, help="cmp_csf.json (ground_source=csf arm)")
    ap.add_argument("--class2", required=True, help="cmp_class2.json (class2 arm)")
    ap.add_argument("--project", required=True, help="the gen1 acquisition, e.g. lidar_semn2008")
    ap.add_argument("--out", default="", help="optional JSON to write")
    a = ap.parse_args()

    csf = json.load(open(a.csf))
    cl2 = json.load(open(a.class2))
    control = G.control_for_survey(a.project)

    R = Run("Do the 2008 control marks prefer the VENDOR's delivered ground surface, our "
            "CSF reconstruction, or our class-2 reconstruction -- and does the answer "
            "depend on land cover?")
    R.input(a.csf, role="bridge arm with ground_source=csf; supplies bridge_mm per mark")
    R.input(a.class2, role="bridge arm with ground_source=class2; IDENTICAL in every other "
                           "parameter, so the arms differ only in ground selection")
    R.input(control.origin, role="gen1's own 2008 control; dnr_error_m is the REPORT's own "
                                 "Control Z - Surface Z against the DELIVERED surface")

    for k, v in csf["params"].items():
        other = cl2["params"].get(k)
        if k == "ground_source":
            R.param("ground_source", f"{v} vs {other}", src="andy",
                    why="THE variable under test. Andy: 'CSF might be better than the "
                        "vendor's surface.' Every other parameter is held identical")
            continue
        if k == "out":
            continue          # the two arms must write to different files, by construction
        if v != other:
            raise SystemExit(f"REFUSED: the arms differ in {k!r} ({v!r} vs {other!r}). "
                             f"They must differ ONLY in ground_source or the contrast is "
                             f"not attributable to the ground source.")
    R.param("covers", tuple(csf["params"]["covers"]), src="MINE",
            why="widened from the standing open-ground-only datum rule SO THAT the canopy "
                "case is visible -- a vendor classifier's hard cases are vegetated, and an "
                "L1O-only comparison cannot see them. A GROUPING here, not a filter: every "
                "mark in both arms is reported")
    R.param("csf_half_width_m", csf["params"]["csf_half_width_m"], src="MINE",
            why="CSF window half-width; wide enough that the cloth is not dominated by a "
                "crop edge. Crops the class2 arm IDENTICALLY, so the windows match")
    R.param("radii_m", tuple(csf["params"]["radii_m"]), src="MINE",
            why="order-2 fit radius, swept not chosen; a 5 m grid needs 6+ cells to fill "
                "one, which sets the floor. Identical in both arms")
    R.param("lateral_nuth_kaab", "NOT applied", src="MINE",
            why="gen1->gen2 registration, needing gen2 which is not on disk away from "
                "Elba; measured effect on a tie 10.0 mm. It applies EQUALLY to both arms "
                "so it CANCELS in the csf-vs-class2 contrast, and remains in each arm's "
                "offset against the marks")
    R.param("outlier_rule", "NONE APPLIED", src="andy",
            why="no mark is cut. Skew and the extremes are printed so the decision of "
                "whether to cut, and where, stays Andy's")

    R.column("offset_vendor_mm", "surveyed - z_delivered = dnr_error_m*1000; the DNR's own "
                                 "number against the DELIVERED class-2 surface. Positive = "
                                 "the surface reads LOW at the mark")
    R.column("offset_csf_mm", "surveyed - z_ours(csf) = offset_vendor + bridge_csf, mm")
    R.column("offset_class2_mm", "surveyed - z_ours(class2) = offset_vendor + bridge_class2, mm")
    R.column("z_class2_minus_z_csf_mm", "bridge_csf - bridge_class2, mm. The ground-source "
                                        "contrast with nothing else in it. Positive = the "
                                        "class-2 surface sits ABOVE the CSF surface")
    R.column("sd", "sample standard deviation over marks, ddof=1, mm -- the NON-absorbable "
                   "part. A constant mean offset is a datum constant and is removed "
                   "downstream; scatter is not")
    R.column("nmad", "1.4826 * median absolute deviation, mm; shown beside sd so a "
                     "difference driven by one mark is visible without cutting it")

    by_id_c2 = {}
    for r in cl2["marks"]:
        by_id_c2[r["point_id"]] = r
    err = {}
    for m in control.marks:
        if m.dnr_error_m is not None:
            for al in m.aliases:
                err[al] = m.dnr_error_m * 1000.0

    rows, missing_arm, missing_err = [], [], []
    for r in csf["marks"]:
        pid = r["point_id"]
        o = by_id_c2.get(pid)
        if o is None:
            missing_arm.append(pid)
            continue
        if pid not in err:
            missing_err.append(pid)
            continue
        ov = err[pid]
        rows.append(dict(point_id=pid, cover=r["cover"], tile=r["tile"],
                         offset_vendor_mm=ov,
                         offset_csf_mm=ov + r["bridge_mm"],
                         offset_class2_mm=ov + o["bridge_mm"],
                         z_class2_minus_z_csf_mm=r["bridge_mm"] - o["bridge_mm"],
                         spread_csf_mm=r["radius_spread_mm"],
                         spread_class2_mm=o["radius_spread_mm"]))

    R.mask("paired", np.ones(len(rows), bool),
           defn=f"marks present in BOTH arms AND carrying the report's own dnr_error_m. "
                f"{len(missing_arm)} absent from the class2 arm, {len(missing_err)} "
                f"carry no dnr_error_m; both counts are reported, never silently dropped",
           of=len(csf["marks"]))
    if missing_arm:
        print(f"    absent from class2 arm: {', '.join(missing_arm)}")
    if missing_err:
        print(f"    no dnr_error_m in the report: {', '.join(missing_err)}")

    covers = sorted({r["cover"] for r in rows})
    print(f"\n{'cover':<8}{'n':>4}  {'surface':<9}{'mean':>9}{'median':>9}{'sd':>9}"
          f"{'nmad':>9}{'SE':>8}{'skew':>7}{'min':>9}{'max':>9}")
    out = {}
    for cov in covers + ["ALL"]:
        sel = rows if cov == "ALL" else [r for r in rows if r["cover"] == cov]
        if not sel:
            continue
        out[cov] = {}
        for lab, key in (("vendor", "offset_vendor_mm"), ("csf", "offset_csf_mm"),
                         ("class2", "offset_class2_mm")):
            s = stats([r[key] for r in sel])
            out[cov][lab] = s
            print(f"{cov if lab=='vendor' else '':<8}{s['n'] if lab=='vendor' else '':>4}  "
                  f"{lab:<9}{s['mean']:>9.1f}{s['median']:>9.1f}{s['sd']:>9.1f}"
                  f"{s['nmad']:>9.1f}{s['se']:>8.1f}{s['skew']:>7.2f}"
                  f"{s['lo']:>9.1f}{s['hi']:>9.1f}")
        d = stats([r["z_class2_minus_z_csf_mm"] for r in sel])
        out[cov]["z_class2_minus_z_csf"] = d
        print(f"{'':<8}{'':>4}  {'  ^contrast':<9}{d['mean']:>9.1f}{d['median']:>9.1f}"
              f"{d['sd']:>9.1f}{d['nmad']:>9.1f}{d['se']:>8.1f}{d['skew']:>7.2f}"
              f"{d['lo']:>9.1f}{d['hi']:>9.1f}")
        print()

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(dict(rows=rows, summary=out,
                           params=csf["params"],
                           missing_from_class2_arm=missing_arm,
                           missing_dnr_error=missing_err), fh, indent=1)
            fh.write("\n")
        print(f"    wrote {a.out}")

    al = out.get("ALL", {})
    if al:
        R.done(headline=(
            f"n={al['vendor']['n']} paired marks; sd of the offset about its own mean: "
            f"vendor {al['vendor']['sd']:.1f}, csf {al['csf']['sd']:.1f}, "
            f"class2 {al['class2']['sd']:.1f} mm; ground-source contrast "
            f"(z_class2-z_csf) {al['z_class2_minus_z_csf']['mean']:+.1f} "
            f"+/- {al['z_class2_minus_z_csf']['sd']:.1f} mm"))


if __name__ == "__main__":
    main()
