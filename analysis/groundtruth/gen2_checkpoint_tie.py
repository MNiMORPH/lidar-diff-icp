#!/usr/bin/env python3
"""Test gen2 (2021 3DEP) DIRECTLY against its own surveyed QA checkpoints.

Why this one needs no chain
---------------------------
The six marks are the 2021 3DEP project's own vertical-accuracy checkpoints, surveyed in
NAVD88(GEOID18) -- the same vertical datum as gen2 (``analysis/ridgelines/
ABSOLUTE_ELEVATION_REFS.md`` section 1a, verified there from the shapefile ``.prj`` and
the project report). So every term the gen1 tie needs disappears:

* no swath chain -- gen2 covers the marks directly;
* no geoid shift -- the mark and the survey are on the same model;
* no lateral extrapolation -- gen2 is the horizontal frame everything else was shifted
  into.

What is left is the estimator and the terrain, read through the SAME
:func:`lidar_diff_icp.groundtruth.tie.estimate_tie` and the SAME radius ladder as gen1,
so the two epochs' absolute offsets are comparable rather than merely both present.

This has never been run. ``docs/groundtruth.md`` section 7 lists it as a gap: "Whether
the *2021* surface hits its own checkpoints is untested here; six 200 m boxes (~24 MB)
would settle it."

Ground source. gen2's own bare earth is the vendor ASPRS class 2 that
:func:`lidar_diff_icp.pipeline.read_after_ground` uses in its default ``mode="class2"``,
so that is the headline; ``--cross-check-ground`` re-reads each mark with CSF, the gen1
side's source, so the choice is visible rather than assumed away.

Sign convention, the same as the gen1 tie: ``tie_mm = surveyed - z_lidar``, i.e. the
constant to ADD to the lidar. A POSITIVE tie means the survey reads HIGH of the lidar --
the lidar sits LOW.

    scripts/fetch_3dep_checkpoint_boxes.sh          # one 400 m box per mark
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/groundtruth/gen2_checkpoint_tie.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                    # noqa: E402
from lidar_diff_icp.groundtruth import checkpoints as C             # noqa: E402
from lidar_diff_icp.groundtruth import tie as T                     # noqa: E402

BOXES = "data/after/checkpoints"
CACHE = "data/derived/groundtruth"
#: gen1's ties at the two marks that anchor the datum, from
#: ``analysis/groundtruth/elba_absolute_tie.py`` (re-run 2026-08-26, csf ground). Used
#: only to print the epoch difference beside each gen2 result; nothing here depends on
#: them.
GEN1_TIE_MM = {"2210_2021_MN": (21.3, 12.4), "2036_2021_MN": (28.9, 27.0),
               "3056_2021_MN": (-103.2, 52.3), "2024_2021_MN": (156.6, 54.5)}
GEN1_SRC = ("analysis/groundtruth/elba_absolute_tie.py, csf ground, run 2026-08-26; "
            "docs/groundtruth.md section 6")
#: 3DEP's own published accuracy for this block, used ONLY as the verdict tolerance --
#: gen2 is being tested against the very marks that certify it, so its own specification
#: is the right bar and it is stated rather than chosen.
NVA_RMSE_MM = 35.0
VVA_95_MM = 270.0
ACC_SRC = ("USGS project report / VA text files for the QL1 block covering Winona "
           "County: NVA 3.5 cm RMSEz, VVA 27 cm at the 95th percentile; quoted in "
           "analysis/ridgelines/ABSOLUTE_ELEVATION_REFS.md section 1b")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ground", choices=("class2", "csf"), default="class2",
                    help="gen2 ground source (repo default: class2, "
                         "pipeline.read_after_ground mode)")
    ap.add_argument("--half-width", type=float, default=200.0,
                    help="half-width of the fetched gen2 box, m")
    ap.add_argument("--json", default=None,
                    help="also write every per-checkpoint result to this JSON, so a "
                         "downstream step consumes the artifact instead of a transcription")
    ap.add_argument("--cross-check-ground", action="store_true", default=True,
                    help="also read each tie with the other ground source")
    ap.add_argument("--no-cross-check-ground", dest="cross_check_ground",
                    action="store_false")
    A = ap.parse_args()

    R = Run("does the 2021 gen2 surface hit its own surveyed QA checkpoints, and by how "
            "much does it sit low or high at each one?")
    cps = C.load_bundled()
    R.input(cps.origin, role="surveyed 3DEP QA checkpoints near Elba, NAVD88(GEOID18) -- "
                             "the SAME datum as gen2, so no geoid term is applied")
    R.param("checkpoint_set", os.path.basename(cps.origin), src="repo",
            why="the six 3DEP vertical-accuracy marks within ~12 km of Elba")
    R.param("ground_source", A.ground, src="repo",
            why="pipeline.read_after_ground default mode='class2' -- gen2's own vendor "
                "bare earth, the ground the product is built from")
    R.param("geoid_shift_mm", 0.0, src="repo",
            why="checkpoints are NAVD88(GEOID18) and gen2 is NAVD88(GEOID18) "
                "(ABSOLUTE_ELEVATION_REFS.md section 1a, verified from the shapefile "
                ".prj and the project report), so NO geoid term applies on this side. "
                "The measurement below is itself the check on that claim: a wrong geoid "
                "model would show as a ~70 mm offset at every mark.")
    R.param("swath_shift_m", (0.0, 0.0, 0.0), src="repo",
            why="gen2 is the frame gen1 was registered INTO; it carries no chain, no "
                "internal swath alignment and no lateral shift of its own")
    R.param("box_half_width_m", A.half_width, src="MINE",
            why="the fetched gen2 box around each mark. The widest fitting radius on the "
                "ladder is 25 m, so this excludes nothing from any estimate; it is "
                "sized so the 3DEP node tiles fetch in seconds.")
    R.param("nva_tolerance_mm", NVA_RMSE_MM, src="repo", why=ACC_SRC)
    R.param("vva_tolerance_mm", VVA_95_MM, src="repo", why=ACC_SRC)
    R.param("gen1_ties_for_comparison_mm", GEN1_TIE_MM, src="repo", why=GEN1_SRC)

    for k, v in T.TieEstimate.table_columns().items():
        R.column(k, v)
    R.column("checkpoint", "surveyed 3DEP QA mark id")
    R.column("type", "NVA = open ground, VVA = under vegetation (3DEP convention)")
    R.column("n", "gen2 ground returns inside the report radius, count")
    R.column("gen2_mm", "surveyed minus gen2 ground at the mark -- the constant to ADD "
                        "to gen2; POSITIVE means gen2 sits LOW, mm")
    R.column("sigma_mm", "half the tie's spread across the pipeline-scale radii, mm")
    R.column("gen2_med_mm", "the same tie taken as the median over the pipeline-scale radii, mm")
    R.column("gen1_mm", "gen1's tie at the same mark, for comparison, mm")
    R.column("epoch_mm", "gen1 tie minus gen2 tie: how much MORE the 2008 surface sits "
                         "low than the 2021 one at this mark, mm")
    R.column("usable", "does the radius spread fall inside 3DEP's own accuracy for the "
                       "point's class")
    R.column("note", "why a mark has no gen2 tie, where it has none")
    R.column("density", "gen2 ground returns per square metre inside the box")

    boxes = {}
    for cp in cps.usable():
        p = os.path.join(BOXES, f"cp{cp.point_id.split('_')[0]}_gen2.laz")
        if os.path.exists(p):
            boxes[cp.point_id] = R.input(
                p, role=f"gen2 2021 3DEP full-density box around checkpoint "
                        f"{cp.point_id} (fetched at EPT depth 12)")
    R.banner()

    results = []
    for cp in cps.usable():
        print(f"\n{'=' * 78}\n== {cp.point_id}  ({cp.point_type})  surveyed "
              f"{cp.elevation_m:.3f} m {cp.vertical_datum}/{cp.geoid_model} ==")
        if cp.point_id not in boxes:
            print("  no gen2 box on disk for this mark. Reported, not dropped.")
            results.append(dict(cp=cp, est=None, reason="no gen2 box on disk"))
            continue
        path = boxes[cp.point_id]
        loader = (T.vendor_ground_near if A.ground == "class2" else T.csf_ground_near)
        kw = dict(cache_dir=os.path.join(CACHE, "csf_gen2")) if A.ground == "csf" else {}
        g = loader(path, cp.easting, cp.northing, A.half_width, **kw)
        dens = len(g) / (2 * A.half_width) ** 2
        print(f"  ground: {g.source}, {len(g):,} returns from {g.n_input:,} in the "
              f"{2 * A.half_width:.0f} m box  ({dens:.2f} ground pts/m^2)")

        # No chain, no geoid, no lateral shift -- gen2 IS the frame and IS the datum.
        est = T.estimate_tie(cp, g, line=None, geoid_shift_m=0.0,
                             swath_shift_m=(0.0, 0.0, 0.0))
        R.table(list(T.TieEstimate.table_columns()), est.table_rows())
        tol = NVA_RMSE_MM if cp.point_type.upper() == "NVA" else VVA_95_MM
        ok, why = est.verdict(tol, tolerance_source=ACC_SRC)
        print(f"  gen2 TIE = {est.tie_mm:+.1f} +/- {est.sigma_mm:.1f} mm "
              f"(n={est.n_report} at R={est.report_radius_m:g} m); "
              f"{'gen2 sits LOW' if est.tie_mm > 0 else 'gen2 sits HIGH'} of the survey")
        print(f"  median over the pipeline-scale radii: {est.tie_median_mm:+.1f} mm")
        print(f"  radius spread {est.radius_spread_mm:.1f} mm over the pipeline-scale "
              f"radii, {est.radius_spread_all_mm:.1f} mm over the whole ladder")
        print(f"  usable: {ok} -- {why}")
        other_delta = None
        if A.cross_check_ground:
            other = "csf" if A.ground == "class2" else "class2"
            ldr = T.csf_ground_near if other == "csf" else T.vendor_ground_near
            kw2 = dict(cache_dir=os.path.join(CACHE, "csf_gen2")) if other == "csf" else {}
            try:
                g2 = ldr(path, cp.easting, cp.northing, A.half_width, **kw2)
                e2 = T.estimate_tie(cp, g2, line=None, geoid_shift_m=0.0,
                                    swath_shift_m=(0.0, 0.0, 0.0))
                other_delta = e2.tie_mm - est.tie_mm
                print(f"  ground-source check: {other} gives {e2.tie_mm:+.1f} mm "
                      f"({other_delta:+.1f} mm from {A.ground})")
            except Exception as exc:                                  # reported, not hidden
                print(f"  ground-source check ({other}) FAILED: {exc}")
        results.append(dict(cp=cp, est=est, dens=dens, ok=ok, reason="",
                            other_ground_mm=other_delta))

    # -------------------------------------------------------------------- every mark
    print(f"\n{'=' * 78}\n== gen2 against every checkpoint, separately ==")
    rows = []
    for r in results:
        cp = r["cp"]
        g1 = GEN1_TIE_MM.get(cp.point_id)
        if r["est"] is None:
            rows.append([cp.point_id, cp.point_type, "-", "-", "-", "-",
                         f"{g1[0]:+.1f}" if g1 else "-", "-",
                         f"NOT ATTEMPTED: {r['reason']}"])
            continue
        e = r["est"]
        rows.append([cp.point_id, cp.point_type, e.n_report, f"{e.tie_mm:+.1f}",
                     f"{e.sigma_mm:.1f}", f"{e.tie_median_mm:+.1f}",
                     f"{g1[0]:+.1f}" if g1 else "-",
                     f"{g1[0] - e.tie_mm:+.1f}" if g1 else "-", str(r["ok"])])
    R.table(["checkpoint", "type", "n", "gen2_mm", "sigma_mm", "gen2_med_mm", "gen1_mm",
             "epoch_mm", "usable" if all(r["est"] is not None for r in results) else "note"],
            rows)

    got = [r for r in results if r["est"] is not None]
    nva = [r for r in got if r["cp"].point_type.upper() == "NVA"]
    print("\n== strata, declared by the survey's own accuracy class ==")
    for label, sel in (("all marks with a gen2 box", got),
                       ("NVA only (open ground)", nva),
                       ("VVA only (under vegetation)",
                        [r for r in got if r["cp"].point_type.upper() == "VVA"])):
        if not sel:
            continue
        v = [r["est"].tie_mm for r in sel]
        print(f"  {label:32s} n={len(sel)}  median {np.median(v):+7.1f} mm  "
              f"spread {max(v) - min(v):6.1f} mm  RMS {np.sqrt(np.mean(np.square(v))):6.1f} mm")

    if nva:
        v = [r["est"].tie_mm for r in nva]
        print(f"\n  gen2 on open ground is within {max(np.abs(v)):.0f} mm of the survey at "
              f"every NVA mark, against its own published {NVA_RMSE_MM:.0f} mm RMSEz.")
    pairs = [(r["cp"].point_id, GEN1_TIE_MM[r["cp"].point_id][0] - r["est"].tie_mm)
             for r in nva if r["cp"].point_id in GEN1_TIE_MM]
    if pairs:
        print("\n== the epoch difference, mark by mark (gen1 tie minus gen2 tie) ==")
        print("  gen1's tie already carries the GEOID03->GEOID18 shift and the chain into "
              "the elbaext swath frame, so this difference is the epoch residual AFTER "
              "the geoid term, for the gen1 surface AS THE PRODUCT USES IT -- not for the "
              "raw 2008 returns at the mark.")
        for pid, d in pairs:
            print(f"    {pid:16s} {d:+7.1f} mm")
        print(f"  median {np.median([d for _, d in pairs]):+.1f} mm over {len(pairs)} "
              "open-ground marks. Compare the leveled benchmark DG8385, which put the "
              "residual after the GEOID03->GEOID18 shift at +5 mm with a std of 11 mm "
              "(analysis/ridgelines/FLAGPOLE_ABSOLUTE_TEST.md).")
    if A.json:
        os.makedirs(os.path.dirname(A.json) or ".", exist_ok=True)
        v_nva = [r["est"].tie_mm for r in nva]
        rec = dict(
            produced_by="analysis/groundtruth/gen2_checkpoint_tie.py",
            ground_source=A.ground,
            sign_convention=("tie_mm = surveyed - z_gen2: the constant to ADD to gen2. "
                             "POSITIVE means gen2 sits LOW of the survey. No geoid, no "
                             "chain and no lateral term are applied -- gen2 and the marks "
                             "share NAVD88(GEOID18)."),
            nva_rms_mm=float(np.sqrt(np.mean(np.square(v_nva)))) if v_nva else None,
            nva_median_mm=float(np.median(v_nva)) if v_nva else None,
            nva_spread_mm=float(max(v_nva) - min(v_nva)) if len(v_nva) > 1 else None,
            n_nva=len(nva),
            checkpoints=[])
        for r in results:
            cp = r["cp"]
            if r["est"] is None:
                rec["checkpoints"].append(dict(point_id=cp.point_id,
                                               point_type=cp.point_type,
                                               attempted=False, reason=r["reason"]))
                continue
            e = r["est"]
            rec["checkpoints"].append(dict(
                point_id=cp.point_id, point_type=cp.point_type, attempted=True,
                n_report=e.n_report, report_radius_m=e.report_radius_m,
                ground_density_per_m2=r["dens"], tie_mm=e.tie_mm, sigma_mm=e.sigma_mm,
                tie_median_mm=e.tie_median_mm, radius_spread_mm=e.radius_spread_mm,
                radius_spread_all_mm=e.radius_spread_all_mm, fit_se_mm=e.fit_se_mm,
                other_ground_delta_mm=r["other_ground_mm"], usable=bool(r["ok"])))
        with open(A.json, "w") as f:
            json.dump(rec, f, indent=2)
        print(f"\n  wrote {A.json}")
    R.done(headline=("gen2 vs its own checkpoints: " +
                     ", ".join(f"{r['cp'].point_id.split('_')[0]}="
                               f"{r['est'].tie_mm:+.0f}mm" for r in got)))


if __name__ == "__main__":
    main()
