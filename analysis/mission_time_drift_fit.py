#!/usr/bin/env python3
"""Fit the single-curve mission-time drift model F(t) = (Fx, Fy, Fz) and compare it,
on evidence, with the per-swath constants the pipeline uses now.

The comparison is deliberately made against the SAVED per-swath constants themselves,
not against a re-run of ``align_swaths``. Five of the seven sites were aligned on a
CSF-classified cloud whose cache (``data/csf_cache/<site>.las``) no longer exists, and
``coregister_swaths`` excludes classes (5, 6, 9) -- which CSF overwrites. Re-running on
the vendor classification would produce different numbers, so the two models would no
longer be compared the same way. Fitting F(t) TO the saved constants keeps the method
identical, and the residual is then directly the per-swath vertical bias, in mm, that
substituting F(t) would introduce.

Sections
  1. per-site fit of F linear in gps_time, per axis, and its residual to the saved
     constants; alongside F linear in flight-line ordinal, which has the same 2 parameters.
  2. parameter count: free parameters after the gauge, current model vs F(t).
  3. what "small residual" has to beat -- the repeatability of the per-swath constant
     itself, measured from the two Elba tiles, and each site's stable-ground 1-sigma.
  4. does F(t) make elba and elbaext agree about the same swaths? Fit each tile
     separately, evaluate both at the shared swaths, re-reference to swath 135, and see
     whether the 8-17 mm disagreement shrinks.
  5. IDENTIFIABILITY and the two timescales: what the between-swath gaps cannot
     constrain, and the amplitude of the WITHIN-swath along-track drift already in
     corrections.json, which the same curve would also have to carry.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/mission_time_drift_fit.py

Reads only the cached swath statistics written by analysis/mission_time_drift.py and the
saved corrections.json files. No thresholds, no count filters.
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from trust.provenance import Run
from analysis.mission_time_drift import SITES, CACHE, ols, fmt


def main():
    R = Run("does one linear mission-time curve F(t) reproduce the gen1 per-swath constants, "
            "and does it make the two Elba tiles agree about the same flight lines?")
    R.input(CACHE, role="per-swath gps_time window, midpoint and return count, from the gen1 LAZ files")
    for name, (_laz, ddir) in SITES.items():
        R.input(f"{ddir}/corrections.json",
                role=f"{name}: saved per-swath (dx,dy,dz), along-track drift curves, stable 1-sigma")
    R.param("swath_reference", "min(point_source_id)", src="repo",
            why="align_swaths gauge in pipeline.difference_dem; that swath is pinned to zero "
                "and is excluded from every fit below")
    R.column("site", "study site = one gen1 tile/merge")
    R.column("axis", "dz (mm), dx (m) or dy (m) -- the alignment component being modelled")
    R.column("nsw", "swaths in the fit (the forced-zero reference swath excluded)")
    R.column("rate", "fitted F rate: mm/hour for dz, m/hour for dx and dy")
    R.column("rms_time", "RMS residual of F linear in gps_time to the saved constants (mm or m)")
    R.column("rms_ord", "RMS residual of F linear in flight-line ordinal, same units")
    R.column("max_time", "largest absolute residual of the gps_time model (mm or m)")
    R.column("p_now", "free parameters of the current model on this axis: n swaths minus the gauge")
    R.column("p_Ft", "free parameters of linear F on this axis after the gauge: 1 (the rate)")
    R.column("stable_sigma_mm", "stable_1sigma_m from the site's corrections.json, mm")
    R.column("swath", "point_source_id shared by the two Elba tiles")
    R.column("obs_mm", "observed elbaext-minus-elba dz for that swath, both re-referenced to 135, mm")
    R.column("Ft_mm", "the same difference predicted by each tile's own fitted F(t), mm")
    R.column("Fk_mm", "the same difference predicted by each tile's own fit linear in swath id, mm")
    R.column("model", "which model the row summarises")
    R.column("sw135", "elbaext minus elba for swath 135, mm (zero by the shared re-reference)")
    R.column("sw136", "elbaext minus elba for swath 136, re-referenced to swath 135, mm")
    R.column("sw137", "elbaext minus elba for swath 137, re-referenced to swath 135, mm")
    R.column("sw138", "elbaext minus elba for swath 138, re-referenced to swath 135, mm")
    R.column("rms_mm", "RMS over the shared swaths, mm")
    R.column("pair", "consecutive pair of swaths in time order")
    R.column("gap_h", "elapsed time between consecutive swath midpoints, hours")
    R.column("obs_h", "hours of the gap in which returns actually exist (the two 45-60 s windows)")
    R.column("frac_obs", "obs_h / gap_h -- the fraction of the interval F(t) is constrained on")
    R.column("drift_pp_mm", "peak-to-peak of the saved WITHIN-swath along-track drift curve, mm")
    R.column("drift_rate_mm_h", "drift_pp_mm divided by the swath's own duration, mm/hour")
    R.column("between_rate_mm_h", "the site's fitted BETWEEN-swath rate, mm/hour, for comparison")
    R.banner()

    db = json.load(open(CACHE))
    per = {}
    for name, (_laz, ddir) in SITES.items():
        C = json.load(open(f"{ddir}/corrections.json"))
        corr = C["per_swath_internal_alignment_dxdydz_m"]
        ref = min(int(s) for s in corr)
        rec = []
        for s, (dx, dy, dz) in corr.items():
            st = db[name]["swaths"].get(s)
            if st is None:
                continue
            rec.append(dict(swath=int(s), dx=dx, dy=dy, dz_mm=1e3 * dz,
                            t=st["t_mid"], span=st["t_max"] - st["t_min"], n=st["n"],
                            is_ref=(int(s) == ref)))
        rec.sort(key=lambda r: r["t"])
        per[name] = dict(rec=rec, ref=ref, C=C)

    # -------------------------------------------------------------- 1. the fits
    print("\n### 1. F linear in gps_time, fitted to the saved per-swath constants\n")
    rows = []
    for name in SITES:
        obs = [r for r in per[name]["rec"] if not r["is_ref"]]
        t = np.array([r["t"] for r in obs]); k = np.arange(len(obs), dtype=float)
        for ax, key in (("dz (mm)", "dz_mm"), ("dx (m)", "dx"), ("dy (m)", "dy")):
            v = np.array([r[key] for r in obs])
            b, a, _, n = ols(t, v)
            res_t = v - (a + b * t)
            bk, ak, _, _ = ols(k, v)
            res_k = v - (ak + bk * k)
            rows.append([name, ax, n, fmt(3600 * b, 3),
                         fmt(np.sqrt(np.mean(res_t ** 2)), 2),
                         fmt(np.sqrt(np.mean(res_k ** 2)), 2),
                         fmt(np.max(np.abs(res_t)), 2)])
    R.table(["site", "axis", "nsw", "rate", "rms_time", "rms_ord", "max_time"], rows)
    print("  With 3 observations a 2-parameter line has 1 residual degree of freedom, so a")
    print("  small rms there is nearly automatic; elbaext (5) and mnrv (6) carry the weight.")

    # -------------------------------------------------------------- 2. parameters
    print("\n### 2. free parameters per axis, after the gauge\n")
    rows = []
    for name in SITES:
        nsw = len(per[name]["rec"])
        rows.append([name, "dz/dx/dy each", nsw, str(nsw - 1), "1",
                     fmt(1e3 * per[name]["C"]["stable_1sigma_m"], 1)])
    R.table(["site", "axis", "nsw", "p_now", "p_Ft", "stable_sigma_mm"], rows)

    # -------------------------------------------------------------- 3. what to beat
    print("\n### 3. the bar: repeatability of a per-swath constant, measured")
    ce = json.load(open("data/derived/elba/corrections.json"))["per_swath_internal_alignment_dxdydz_m"]
    cx = json.load(open("data/derived/elbaext/corrections.json"))["per_swath_internal_alignment_dxdydz_m"]
    shared = sorted(set(ce) & set(cx), key=int)
    free = [s for s in shared if s != "135"]
    print("\n  The same four flight lines, aligned twice from two different gen1 extents")
    print("  (elba 2.5x3.5 km single tile; elbaext 4.45x4.05 km merge), re-referenced to")
    print("  swath 135. All three axes, in mm:\n")
    rows = []
    axd = {}
    for ax, j in (("dx (m)", 0), ("dy (m)", 1), ("dz (mm)", 2)):
        axd[ax] = {s: 1e3 * ((cx[s][j] - cx["135"][j]) - (ce[s][j] - ce["135"][j])) for s in shared}
        rows.append([ax] + [fmt(axd[ax][s], 1) for s in shared]
                    + [fmt(np.sqrt(np.mean([axd[ax][s] ** 2 for s in free])), 1)])
    R.table(["axis"] + [f"sw{s}" for s in shared] + ["rms_mm"], rows)
    obs_d = axd["dz (mm)"]
    print("  That is this estimator's own repeatability under a change of extent, and it is")
    print("  the number any residual in section 1 has to be judged against. dx -- the axis")
    print("  carrying the metre-scale offsets -- repeats best of the three.")

    # -------------------------------------------------------------- 4. does F(t) resolve it?
    print("\n### 4. does F(t) make the two Elba tiles agree about the same swaths?\n")
    pred = {}
    for mdl, xkey in (("F(t)", "t"), ("F(swath id)", "swath")):
        vals = {}
        for name in ("elba", "elbaext"):
            obs = [r for r in per[name]["rec"] if not r["is_ref"]]
            u = np.array([r[xkey] for r in obs], float)
            v = np.array([r["dz_mm"] for r in obs])
            b, a, _, _ = ols(u, v)
            lut = {str(r["swath"]): r[xkey] for r in per[name]["rec"]}
            base = a + b * float(lut["135"])
            vals[name] = {s: (a + b * float(lut[s])) - base for s in shared}
        pred[mdl] = {s: vals["elbaext"][s] - vals["elba"][s] for s in shared}
    rows = [[s, fmt(obs_d[s], 1), fmt(pred["F(t)"][s], 1), fmt(pred["F(swath id)"][s], 1)]
            for s in shared]
    R.table(["swath", "obs_mm", "Ft_mm", "Fk_mm"], rows)
    rows = [["as saved (per-swath constants)",
             fmt(np.sqrt(np.mean([obs_d[s] ** 2 for s in free])), 1)],
            ["F(t) linear in gps_time",
             fmt(np.sqrt(np.mean([pred["F(t)"][s] ** 2 for s in free])), 1)],
            ["F linear in swath id",
             fmt(np.sqrt(np.mean([pred["F(swath id)"][s] ** 2 for s in free])), 1)]]
    R.table(["model", "rms_mm"], rows)

    # -------------------------------------------------------------- 5. identifiability
    print("\n### 5. what the gaps cannot constrain\n")
    rows = []
    for name in SITES:
        rec = per[name]["rec"]
        for a, b in zip(rec[:-1], rec[1:]):
            gap = (b["t"] - a["t"]) / 3600.0
            obsh = (a["span"] + b["span"]) / 2.0 / 3600.0
            rows.append([name, f"{a['swath']}->{b['swath']}", fmt(gap, 3), fmt(obsh, 4),
                         fmt(obsh / gap, 4)])
    R.table(["site", "pair", "gap_h", "obs_h", "frac_obs"], rows)

    print("\n  The two timescales the SAME curve would have to carry: the within-swath")
    print("  along-track drift already fitted by coreg.fit_along_track_drift, against the")
    print("  between-swath rate fitted above.\n")
    rows = []
    for name in SITES:
        C = per[name]["C"]
        obs = [r for r in per[name]["rec"] if not r["is_ref"]]
        brate, _, _, _ = ols([r["t"] for r in obs], [r["dz_mm"] for r in obs])
        span = {str(r["swath"]): r["span"] for r in per[name]["rec"]}
        for s, cur in sorted(C["along_track_drift_gpsTime_to_m"].items(), key=lambda kv: int(kv[0])):
            d = 1e3 * np.asarray(cur["drift_m"], float)
            pp = float(np.nanmax(d) - np.nanmin(d))
            sp = span.get(s, np.nan)
            rows.append([name, s, fmt(pp, 1), fmt(3600 * pp / sp, 0), fmt(3600 * brate, 1)])
    R.table(["site", "swath", "drift_pp_mm", "drift_rate_mm_h", "between_rate_mm_h"], rows)

    R.done(headline="single-curve mission-time model fitted and compared with the per-swath constants")


if __name__ == "__main__":
    main()
