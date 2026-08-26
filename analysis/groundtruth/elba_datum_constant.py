#!/usr/bin/env python3
"""One datum constant for gen1 at Elba, its uncertainty budget, and the absolutely-
referenced gen1 surface.

What this does, in order
------------------------
1. Reads both epochs' checkpoint ties as ARTIFACTS (``elba_absolute_tie.py --json`` and
   ``gen2_checkpoint_tie.py --json``), never as constants transcribed into this file.
2. Declares which marks anchor the datum and, for every mark that does not, the reason --
   the exclusion is an argument with its source, not a hidden filter.
3. Builds the uncertainty budget term by term, MEASURING each term here rather than
   quoting it: the per-link chain repeatability is derived from the two Elba products'
   own per-swath alignments; the estimator-and-mark scatter is derived from the gen2 run,
   where no chain, no geoid and no lateral term exist to confound it.
4. Combines the anchors with :func:`lidar_diff_icp.groundtruth.datum.combine_ties`, which
   keeps common-mode terms at full size.
5. Writes the absolutely-referenced gen1 surface for ``elba_fulldensity`` as a NEW array
   with a JSON sidecar. Nothing existing is overwritten.
6. Measures the term the datum constant does NOT touch: the spatially varying residual on
   stable divide cells, block by block.

Sign convention throughout: the datum constant is what you **ADD to gen1**.
``tie = surveyed - z_lidar_corrected``. The gen1 surface is reconstructed as
``z_after - dod`` (``pipeline.py``: ``dod = Zref - Z08c``), so adding the constant to it
raises gen1 and therefore LOWERS the DoD by the same amount.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/groundtruth/elba_datum_constant.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from trust.provenance import Run                                     # noqa: E402
from lidar_diff_icp.groundtruth.datum import (                        # noqa: E402
    BudgetTerm, DatumConstant, combine_ties)
from lidar_diff_icp.refcells import reference_cells                  # noqa: E402

GEN1_JSON = "data/derived/groundtruth/elba_gen1_ties.json"
GEN2_JSON = "data/derived/groundtruth/elba_gen2_ties.json"
FD = "data/derived/elba_fulldensity"
ELBA_CORR = f"{FD}/corrections.json"
EXT_CORR = "data/derived/elbaext/corrections_geoid.json"
OUT_NPY = f"{FD}/z_before_absolute.npy"
OUT_JSON = f"{FD}/z_before_absolute.json"


def per_link_repeatability(elba_corr, ext_corr, ref=135):
    """Measure the alignment estimator's per-link repeatability from the two Elba tiles.

    The same four flight lines are aligned twice, from two different extents. Referenced
    to a common swath, the two solutions disagree; a swath *n* links away from the
    reference has accumulated *n* independent link errors, so ``|disagreement| / sqrt(n)``
    is one estimate of the per-link error and the RMS over the shared swaths is the
    headline. This is the same comparison ``analysis/MISSION_TIME_DRIFT.md`` section 4
    reports as 12.4 mm RMS on dz; here it is decomposed per link, because the corridor
    chains are 5 and 6 links long.
    """
    a = {int(k): v[2] * 1000.0 for k, v in
         json.load(open(elba_corr))["per_swath_internal_alignment_dxdydz_m"].items()}
    b = {int(k): v[2] * 1000.0 for k, v in
         json.load(open(ext_corr))["per_swath_internal_alignment_dxdydz_m"].items()}
    shared = sorted(set(a) & set(b))
    rows, per = [], []
    for s in shared:
        if s == ref:
            continue
        d = (b[s] - b[ref]) - (a[s] - a[ref])
        n = abs(s - ref)
        rows.append([s, n, f"{d:+.1f}", f"{abs(d) / np.sqrt(n):.1f}"])
        per.append(abs(d) / np.sqrt(n))
    tot = [((b[s] - b[ref]) - (a[s] - a[ref])) for s in shared if s != ref]
    return (float(np.sqrt(np.mean(np.square(per)))),
            float(np.sqrt(np.mean(np.square(tot)))), rows, shared, ref)


def block_residual(tile_dir, corr, block_m, *, slope_max=12.0):
    """Per-block median DoD on stable divide cells: the spatially varying term.

    Uses the repo's own :func:`lidar_diff_icp.refcells.reference_cells` with its own
    defaults, so the stable set is the one the pipeline already trusts. The between-block
    standard deviation is corrected for the within-block sampling noise it inevitably
    contains: ``sd_real^2 = var(block medians) - mean(se_of_a_block_median^2)``, with the
    median's standard error taken as ``1.2533 * sd / sqrt(n)``. Without that correction a
    handful of thin blocks would manufacture spatial structure.
    """
    dod = np.load(os.path.join(tile_dir, "dod.npy")) * 1000.0
    mask, rep = reference_cells(tile_dir, slope_max=slope_max)
    mask = mask.reshape(dod.shape)
    ny, nx = dod.shape
    res = corr["res_m"]
    step = int(round(block_m / res))
    gy, gx = np.mgrid[0:ny, 0:nx]
    blk = (gy // step) * (nx // step + 1) + (gx // step)
    ok = mask & np.isfinite(dod)
    meds, ses, ns = [], [], []
    for b in np.unique(blk[ok]):
        v = dod[ok & (blk == b)]
        if v.size < 2:                     # a median's SE is undefined below 2 samples
            continue
        meds.append(float(np.median(v)))
        ses.append(1.2533 * float(np.std(v, ddof=1)) / np.sqrt(v.size))
        ns.append(int(v.size))
    meds = np.array(meds); ses = np.array(ses)
    var_obs = float(np.var(meds, ddof=1)) if meds.size > 1 else np.nan
    var_noise = float(np.mean(ses ** 2)) if ses.size else np.nan
    sd_real = float(np.sqrt(max(var_obs - var_noise, 0.0)))
    return dict(n_blocks=int(meds.size), n_cells=int(np.sum(ns)),
                lo=float(meds.min()) if meds.size else np.nan,
                hi=float(meds.max()) if meds.size else np.nan,
                sd_obs=float(np.sqrt(var_obs)), sd_noise=float(np.sqrt(var_noise)),
                sd_real=sd_real, report=rep, block_m=block_m,
                median_n_per_block=int(np.median(ns)) if ns else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band-km", type=float, default=1.0,
                    help="a mark anchors the datum only if it lies within this distance "
                         "of the corridor band the chain links were solved in")
    ap.add_argument("--anchor-types", default="NVA",
                    help="comma-separated 3DEP accuracy classes allowed to anchor")
    ap.add_argument("--block-m", type=float, default=500.0,
                    help="block size for the spatially varying residual, m")
    ap.add_argument("--ref-slope-max", type=float, default=12.0,
                    help="reference_cells slope cut, deg (repo default 12.0)")
    ap.add_argument("--no-write", action="store_true",
                    help="compute and report, but do not write the product")
    A = ap.parse_args()
    types = tuple(t.strip().upper() for t in A.anchor_types.split(","))

    R = Run("what single constant places gen1 at Elba on the surveyed datum, how well is "
            "it known, and what does that fix in the DoD?")
    g1 = json.load(open(R.input(GEN1_JSON, role="gen1 checkpoint ties, chained into the "
                                                "elbaext swath frame and geoid-shifted")))
    g2 = json.load(open(R.input(GEN2_JSON, role="gen2 checkpoint ties, direct: no chain, "
                                                "no geoid, no lateral term")))
    R.input(ELBA_CORR, role="elba_fulldensity pipeline corrections: datum method, "
                            "per-swath alignment, grid definition")
    R.input(EXT_CORR, role="elbaext pipeline corrections: the SAME swaths aligned from a "
                           "different extent -- the repeatability measurement")
    R.input(f"{FD}/z_after.npy", role="gen2 gridded ground, 5 m -- the surface dod is "
                                      "measured from")
    R.input(f"{FD}/dod.npy", role="gen2 minus gen1 on the 5 m grid; gen1 = z_after - dod")

    corr = json.load(open(ELBA_CORR))
    R.param("band_km", A.band_km, src="repo",
            why="the on-band stratum elba_absolute_tie.py declares before any tie is "
                "computed; off-band marks carry unmodelled along-track drift at the "
                f"measured {g1['drift_median_mm_per_km']:.1f} mm/km")
    R.param("anchor_types", types, src="repo",
            why="3DEP's own accuracy classes: NVA is open ground at 3.5 cm RMSEz, VVA is "
                "under vegetation at 27 cm @95%. Both are reported; only NVA anchors.")
    R.param("block_m", A.block_m, src="MINE",
            why="block size for the spatially varying residual. Chosen to be far larger "
                "than the 5 m grid and far smaller than the 2.5 x 3.5 km tile so there "
                "are tens of blocks; the result is reported with its block count and its "
                "sampling-noise correction so the choice is visible.")
    R.param("ref_slope_max_deg", A.ref_slope_max, src="repo",
            why="refcells.reference_cells slope_max default")
    R.param("datum_method", corr["cross_epoch_datum"]["method"], src="repo",
            why="elba_fulldensity/corrections.json -- the geoid-only datum the product "
                "was built on; the tie constant is applied ON TOP of it, not instead")
    R.param("grid_bounds", corr["bounds"], src="repo", why="elba_fulldensity corrections.json")

    for k, v in DatumConstant.table_columns().items():
        R.column(k, v)
    R.column("mark", "surveyed 3DEP QA mark id")
    R.column("type", "NVA = open ground, VVA = under vegetation")
    R.column("tie_mm", "constant to ADD to gen1 from this mark, mm")
    R.column("sig_mm", "half its spread across the pipeline-scale radii, mm")
    R.column("links", "cross-swath links from the covering line to the elbaext frame")
    R.column("dN_km", "distance from the corridor band the links were solved in, km")
    R.column("role", "anchor, or the reason this mark does not anchor")
    R.column("swath", "gen1 flight line whose two independent alignments are compared")
    R.column("nlink", "links from the comparison reference swath")
    R.column("disagree_mm", "elbaext minus elba_fulldensity alignment for that swath, mm")
    R.column("per_link_mm", "that disagreement divided by sqrt(nlink), mm")
    R.column("quantity", "what is being reported about the spatially varying residual")
    R.column("value", "its value, with its own units stated in the cell (mm, or a count)")
    R.column("source", "the array, mask or formula the value was computed from")
    R.banner()

    # ------------------------------------------------------- 1. who anchors, and why not
    print("\n== every gen1 control point, and its role ==")
    rows, anchors = [], []
    for c in g1["checkpoints"]:
        if not c["attempted"]:
            rows.append([c["point_id"], c["point_type"], "-", "-", "-", "-",
                         f"NOT ATTEMPTED: {c['reason']}"])
            continue
        why = []
        if c["point_type"].upper() not in types:
            why.append(f"type {c['point_type']} not in {types}")
        if abs(c["dN_from_band_km"]) > A.band_km:
            why.append(f"{abs(c['dN_from_band_km']):.2f} km off band "
                       f"(> {A.band_km:g}); ~{c['unmodelled_drift_mm']:.0f} mm of "
                       "unmodelled drift")
        role = "ANCHOR" if not why else "; ".join(why)
        if not why:
            anchors.append(c)
        rows.append([c["point_id"], c["point_type"], f"{c['tie_mm']:+.1f}",
                     f"{c['sigma_mm']:.1f}", c["links"],
                     f"{c['dN_from_band_km']:+.2f}", role])
    R.table(["mark", "type", "tie_mm", "sig_mm", "links", "dN_km", "role"], rows)
    if len(anchors) < 2:
        raise SystemExit("fewer than two anchors; the two-chain check is what makes this "
                         "a measurement rather than a single reading")
    print(f"  {len(anchors)} anchors on {len({tuple(c['path'][:1]) for c in anchors})} "
          "independent chains: " +
          "; ".join(f"{c['point_id'].split('_')[0]} via "
                    f"{'-'.join(str(n) for n in c['path'])}" for c in anchors))
    spread = max(c["tie_mm"] for c in anchors) - min(c["tie_mm"] for c in anchors)
    sig_list = ", ".join(f"{c['sigma_mm']:.1f}" for c in anchors)
    print(f"  They agree to {spread:.1f} mm, which is SMALLER than either mark's own "
          f"sigma ({sig_list} mm). That is a consistency check, not the accuracy.")

    # ------------------------------------------------- 2. the terms, each measured here
    print("\n== the per-link chain term, measured from the two Elba products ==")
    per_link, tot_rms, lrows, shared, cref = per_link_repeatability(ELBA_CORR, EXT_CORR)
    R.table(["swath", "nlink", "disagree_mm", "per_link_mm"], lrows)
    print(f"  swaths {shared} aligned twice from two extents, referenced to {cref}: "
          f"{tot_rms:.1f} mm RMS overall (MISSION_TIME_DRIFT.md section 4 reports the "
          f"same comparison), {per_link:.1f} mm RMS PER LINK.")
    chain_terms = {c["point_id"]: per_link * np.sqrt(c["links"]) for c in anchors}
    for pid, v in chain_terms.items():
        print(f"    {pid}: {[c['links'] for c in anchors if c['point_id'] == pid][0]} "
              f"links -> {v:.1f} mm")
    chain_mm = float(np.sqrt(np.mean(np.square(list(chain_terms.values())))))

    g2_nva = g2["nva_rms_mm"]
    print(f"\n== the estimator-and-mark term, measured on gen2 where nothing else acts ==")
    print(f"  gen2 shares the marks' datum, covers them directly, and needs no chain, no "
          f"geoid and no lateral shift. Over {g2['n_nva']} open-ground marks it still "
          f"scatters {g2['nva_spread_mm']:.0f} mm, {g2_nva:.1f} mm RMS, median "
          f"{g2['nva_median_mm']:+.1f} mm.")
    print("  Whatever makes a surveyed monument disagree with the lidar ground at that "
          "size -- siting on a crown, the monument against the surface, the survey's own "
          "target definition -- is a property of the MARK and applies just as much when "
          "gen1 is read there. It is an UPPER bound: part of it may be real gen2 spatial "
          "error, which this cannot separate.")

    lat = float(np.sqrt(np.mean([c["lateral_effect_mm"] ** 2 for c in anchors])))
    lat_list = ", ".join(f"{c['lateral_effect_mm']:+.1f}" for c in anchors)
    grd_list = ", ".join(f"{c['other_ground_delta_mm']:+.1f}" for c in anchors)
    lat_all = max(abs(c["lateral_effect_mm"]) for c in g1["checkpoints"] if c["attempted"])
    grd = float(np.sqrt(np.mean([c["other_ground_delta_mm"] ** 2 for c in anchors])))
    ext_alignment = max(abs(float(r[2])) for r in lrows)
    drift = max(c["unmodelled_drift_mm"] for c in anchors)

    terms = [
        BudgetTerm("chain, per-link alignment", chain_mm, "random",
                   f"measured here: the same swaths aligned from two extents disagree by "
                   f"{tot_rms:.1f} mm RMS, {per_link:.1f} mm per link; propagated over "
                   f"each anchor's {[c['links'] for c in anchors]} links",
                   "each anchor's own chain (the two paths are disjoint)"),
        BudgetTerm("estimator + mark scatter", g2_nva, "random",
                   f"measured here on gen2, where no chain, geoid or lateral term exists: "
                   f"{g2_nva:.1f} mm RMS over {g2['n_nva']} open-ground marks "
                   f"(analysis/groundtruth/gen2_checkpoint_tie.py). UPPER bound -- it "
                   f"cannot be separated from real gen2 spatial error.",
                   "each mark independently"),
        BudgetTerm("alignment extent-dependence into the elbaext frame", ext_alignment,
                   "common",
                   f"measured here: elbaext's swath alignment repeats to {ext_alignment:.1f} "
                   "mm at worst against elba_fulldensity's. Treated as common, which is "
                   "CONSERVATIVE: the east chain ends on swath 138 and carries it, the "
                   "west chain ends on 133, the gauge swath, and carries none.",
                   "the east tie via elbaext's swath-138 alignment"),
        BudgetTerm("lateral (Nuth & Kaeaeb) shift, extrapolated", lat, "common",
                   f"measured here: applying/withholding the elbaext shift "
                   f"{tuple(g1['elbaext_lateral_shift_m'])} m moves the anchors by "
                   f"{lat_list} mm; one shift vector, one validity assumption",
                   "all ties (one shift, applied everywhere)"),
        BudgetTerm("ground source (CSF vs vendor class 2)", grd, "common",
                   f"measured here: re-reading the anchors with the other ground source "
                   f"moves them by {grd_list} mm",
                   "all ties (one ground source)"),
        BudgetTerm("horizontal datum EPSG:6344 vs EPSG:26915", 0.0, "common",
                   "measured with pyproj: the transform is null to 0.1 mm at all six "
                   "marks, so using the surveyed coordinates directly costs nothing. The "
                   "NAD83(1986)-to-NAD83(2011) realization difference is not modelled by "
                   "PROJ here and is therefore not in this number.",
                   "all ties"),
        BudgetTerm("along-track drift, not applied", drift, "unmodelled",
                   f"pipeline.fit_along_track_drift regresses against the gen2 grid, "
                   f"which does not exist at the marks. Measured scale "
                   f"{g1['drift_median_mm_per_km']:.1f} mm/km; the anchors sit "
                   f"<= {max(abs(c['dN_from_band_km']) for c in anchors):.2f} km off band, "
                   f"so <= {drift:.0f} mm. A bound, not a distribution.",
                   "both anchors"),
        BudgetTerm("validity of the lateral shift 7-16 km away", lat_all, "unmodelled",
                   f"the shift was measured against gen2 at Elba and the marks are "
                   f"3.1-15.6 km away. Its measured effect across all four attempted "
                   f"marks reaches {lat_all:.1f} mm; whether it holds out there is an "
                   "assumption, so this is a bound on being wrong about it.",
                   "all ties"),
    ]

    dc = combine_ties([(c["point_id"], c["tie_mm"], c["sigma_mm"]) for c in anchors],
                      terms,
                      notes=[f"anchors: {[c['point_id'] for c in anchors]}",
                             f"sign: {g1['sign_convention']}"])

    print("\n== the uncertainty budget ==")
    R.table(list(dc.table_columns()), dc.table_rows())
    print(f"\n  DATUM CONSTANT = {dc.value_mm:+.1f} mm, to be ADDED to gen1.")
    print(f"  sigma_total {dc.sigma_total_mm:.1f} mm (random {dc.random_mm:.1f}, "
          f"common {dc.common_mm:.1f}); largest unmodelled bound "
          f"{dc.unmodelled_mm:.1f} mm, NOT included.")
    print(f"  The two anchors agree to {dc.spread_mm:.1f} mm. The budget is "
          f"{dc.sigma_total_mm / dc.spread_mm:.0f}x that, and the honest reading is that "
          "the agreement is luck within a much wider distribution -- the gen2 run above "
          "shows four marks of the same class spanning "
          f"{g2['nva_spread_mm']:.0f} mm with every transport term removed.")
    print("  The per-mark radius sigmas are already carried by the inverse-variance "
          "weighting and are not repeated as budget rows; the estimator-and-mark term "
          "partly overlaps them, so sigma_total is conservative by that overlap.")

    # ------------------------------------------------ 3. the absolutely-referenced grid
    print("\n== the absolutely-referenced gen1 surface ==")
    z_after = np.load(f"{FD}/z_after.npy")
    dod = np.load(f"{FD}/dod.npy")
    z_gen1 = z_after - dod                       # pipeline.py: dod = Zref - Z08c
    z_abs = z_gen1 + dc.value_mm / 1000.0
    print(f"  gen1 grid reconstructed as z_after - dod: {np.isfinite(z_gen1).sum():,} "
          f"finite cells of {z_gen1.size:,}")
    print(f"  + {dc.value_mm:+.1f} mm -> {OUT_NPY}")
    print(f"  The DoD implied by this surface is the existing dod.npy MINUS "
          f"{dc.value_mm:.1f} mm everywhere: raising gen1 lowers gen2-minus-gen1. "
          "dod.npy itself is NOT modified.")
    side = dc.to_dict()
    side.update(
        product="elba_fulldensity gen1 (2008 MN DNR) bare-earth surface on the surveyed "
                "NAVD88(GEOID18) datum",
        array=os.path.basename(OUT_NPY), units="m", crs=corr["crs"],
        bounds=corr["bounds"], res_m=corr["res_m"], shape=list(z_abs.shape),
        built_from=dict(z_after=f"{FD}/z_after.npy", dod=f"{FD}/dod.npy",
                        formula="z_before_absolute = (z_after - dod) + datum_constant"),
        pipeline_datum=corr["cross_epoch_datum"],
        chains=[dict(point_id=c["point_id"], point_type=c["point_type"],
                     covering_line=c["line"], path=c["path"], links=c["links"],
                     tie_mm=c["tie_mm"], sigma_mm=c["sigma_mm"],
                     chain_mm=c["chain_mm"], geoid_mm=c["geoid_mm"])
                for c in anchors],
        not_anchors=[dict(point_id=c["point_id"], point_type=c["point_type"],
                          tie_mm=c.get("tie_mm"), reason=r[6])
                     for c, r in zip(g1["checkpoints"], rows) if r[6] != "ANCHOR"],
        implied_dod_shift_mm=-dc.value_mm,
        dod_note=("dod.npy is unchanged. The DoD on this absolute basis is "
                  "dod.npy - datum_constant/1000 in metres; a positive datum constant "
                  "raises gen1 and therefore lowers gen2-minus-gen1."),
        gen2_absolute_offsets=g2["checkpoints"],
        produced_by="analysis/groundtruth/elba_datum_constant.py")
    if not A.no_write:
        for p in (OUT_NPY, OUT_JSON):
            if os.path.exists(p):
                raise SystemExit(f"{p} exists; refusing to overwrite. Move it aside first.")
        np.save(OUT_NPY, z_abs)
        with open(OUT_JSON, "w") as f:
            json.dump(side, f, indent=2)
        print(f"  wrote {OUT_NPY} and {OUT_JSON}")
    else:
        print("  --no-write: nothing written")

    # ------------------------------------ 4. what the constant does NOT fix
    print("\n== what a datum constant does NOT fix: the spatially varying residual ==")
    br = block_residual(FD, corr, A.block_m, slope_max=A.ref_slope_max)
    print("  reference_cells cuts, in order: " +
          ", ".join(f"{k}={v}" for k, v in br["report"].items()))
    R.table(["quantity", "value", "source"], [
        ["blocks with >= 2 stable cells", f"{br['n_blocks']}",
         f"{A.block_m:g} m blocks over the elba_fulldensity grid"],
        ["stable cells in them", f"{br['n_cells']:,}",
         f"refcells.reference_cells(slope_max={A.ref_slope_max:g})"],
        ["median cells per block", f"{br['median_n_per_block']}", "same"],
        ["per-block median DoD, range", f"{br['lo']:+.0f} to {br['hi']:+.0f} mm",
         "median of dod.npy inside each block"],
        ["between-block sd, observed", f"{br['sd_obs']:.1f} mm", "sd of those medians"],
        ["within-block sampling noise", f"{br['sd_noise']:.1f} mm",
         "root-mean-square standard error of a block median, 1.2533*sd/sqrt(n)"],
        ["between-block sd, real", f"{br['sd_real']:.1f} mm",
         "observed variance minus the sampling variance"],
    ])
    print(f"  A constant moves every cell by the same {dc.value_mm:+.1f} mm and therefore "
          f"removes NONE of this. The two are separate error terms: the datum fixes the "
          f"LEVEL of the DoD, the block residual is its SHAPE, and only the first is "
          f"what a ground-control tie can reach.")
    R.done(headline=(f"gen1 datum constant at Elba = {dc.value_mm:+.1f} mm "
                     f"+/- {dc.sigma_total_mm:.1f} mm (random {dc.random_mm:.1f}, common "
                     f"{dc.common_mm:.1f}), unmodelled bound {dc.unmodelled_mm:.0f} mm; "
                     f"spatial residual sd {br['sd_real']:.1f} mm is untouched by it"))


if __name__ == "__main__":
    main()
