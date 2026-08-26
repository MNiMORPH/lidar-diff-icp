#!/usr/bin/env python3
"""Is the gen1 per-swath (dx, dy, dz) constant a SAMPLE OF ONE CONTINUOUS DRIFT IN
MISSION TIME, or is it a per-flight-line quantity that mission time merely correlates with?

Andy's hypothesis: the decomposition the pipeline currently uses -- a constant
(dx, dy, dz) per point_source_id from ``coreg.align_swaths`` PLUS a within-swath
along-track curve from ``coreg.fit_along_track_drift`` -- is artificial. There is one
position error F(t) over mission time; the per-swath constants are just F(t) sampled in
each line's 45-second window, and the windows are separated by hours.

The hypothesis makes a FALSIFIABLE prediction that a bare correlation with gps_time does
not test. Within a tile, flight-line ORDINAL, ACROSS-TRACK POSITION and MISSION TIME are
all collinear -- lines are flown in order and laid down side by side -- so "dz rises with
gps_time" is equally "dz rises with line number". The two separate only through the GAP
STRUCTURE: consecutive-line gaps here range from 688 s to 70,966 s. A continuous F(t)
predicts the step between two lines scales with the ELAPSED TIME between them; a per-line
quantity predicts it scales with the NUMBER OF LINES. This script measures both.

Sections
  1. per-swath table: gps_time window + midpoint + return count (from the gen1 LAZ the
     pipeline actually aligned), joined to the saved (dx, dy, dz).
  2. dz, dx, dy vs gps_time per site: Pearson r and slope, reference swath excluded
     (it is forced to zero by construction) and, for disclosure, included.
  3. the map-frame (dx, dy) rotated into the flight's own ALONG/ACROSS-track frame,
     since a north-south line and an east-west line put the same physical error on
     different map axes.
  4. THE TEST: consecutive-pair step vs elapsed gap (per hour) against step per line.
  5. two-parameter linear-in-time model vs two-parameter linear-in-ordinal model.
  6. the elba / elbaext anomaly: the two tiles disagree about the same swaths by
     8-17 mm. Under the hypothesis that is because each tile averages F(t) over its own
     along-track segment -- so it is a quantitative claim about the tiles' gps_time
     midpoints, and this section checks the mm-per-second it would require.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/mission_time_drift.py
    ... --rebuild        # re-read every gen1 LAZ instead of using the cached swath stats

No count filter is applied anywhere. Swaths with very few returns (mnrv 6251: 1,087;
battlecreek 1102: 90) are PRINTED with their counts and kept in the tables; --min-pts
additionally reports every fit with them removed, using coreg.py's own per-swath value.
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from trust.provenance import Run

# ---------------------------------------------------------------- site inventory
# before-epoch (gen1) file + derived dir. The gen1 paths are the ones the products were
# actually made from: scripts/run_all_sites.py SITES for the six standard sites, and
# data/derived/elbaext/meta.json["before"] for elbaext.
SITES = {
    "elba":        ("data/before/4342-29-64.laz",              "data/derived/elba"),
    "elbaext":     ("data/before/elbaext_gen1_merged.laz",     "data/derived/elbaext"),
    "whitewater":  ("data/before/4358-26-03.laz",              "data/derived/whitewater"),
    "mnrv":        ("data/before_mnrv/4342-23-01.laz",         "data/derived/mnrv"),
    "cook":        ("data/before_ne/1158-31-59.laz",           "data/derived/cook"),
    "carlton":     ("data/before_carlton/2742-12-53.laz",      "data/derived/carlton"),
    "battlecreek": ("data/before_battlecreek/4342-03-32_b_a.laz", "data/derived/battlecreek"),
}
# derived dirs that re-use another site's alignment verbatim (same gen1 file, same
# per_swath block) -- listed so the reader knows they are not independent evidence.
ALIASES = {"elba_fulldensity": "elba", "elba_refdatum": "elba", "final": "elba",
           "carlton_density": "carlton", "ne": "cook"}

CACHE = "data/derived/mission_time/swath_gps_stats.json"


# ---------------------------------------------------------------- LAZ pass
def swath_gps_stats(laz_path, bounds):
    """Per point_source_id, over the WHOLE gen1 file (which is the extent
    ``coreg.align_swaths`` sees -- pipeline.difference_dem reads the file uncropped and
    aligns before any bounds clip) and, separately, over the analysis bounds (which is
    what ``fit_along_track_drift``'s stable points are limited to).

    Also fits x(t), y(t) per swath by least squares to get the flight heading; t and x,y
    are referenced to the first point seen for that swath so the normal equations stay
    conditioned (adjusted-standard gps_time is ~3.6e7 s).
    """
    import laspy
    X0, Y0, X1, Y1 = bounds
    acc = {}
    f = laspy.open(laz_path)
    for ch in f.chunk_iterator(3_000_000):
        ps = np.asarray(ch.point_source_id)
        gt = np.asarray(ch.gps_time, float)
        x = np.asarray(ch.x); y = np.asarray(ch.y)
        inb = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)
        for s in np.unique(ps):
            m = ps == s
            a = acc.get(int(s))
            if a is None:
                a = acc[int(s)] = dict(n=0, tmin=np.inf, tmax=-np.inf, tsum=0.0,
                                       n_in=0, tmin_in=np.inf, tmax_in=-np.inf, tsum_in=0.0,
                                       t0=float(gt[m][0]), x0=float(x[m][0]), y0=float(y[m][0]),
                                       St=0.0, Stt=0.0, Sx=0.0, Sy=0.0, Stx=0.0, Sty=0.0)
            t = gt[m] - a["t0"]
            a["n"] += int(m.sum())
            a["tmin"] = min(a["tmin"], float(gt[m].min()))
            a["tmax"] = max(a["tmax"], float(gt[m].max()))
            a["tsum"] += float(gt[m].sum())
            a["St"] += float(t.sum()); a["Stt"] += float((t * t).sum())
            a["Sx"] += float((x[m] - a["x0"]).sum()); a["Sy"] += float((y[m] - a["y0"]).sum())
            a["Stx"] += float(t @ (x[m] - a["x0"])); a["Sty"] += float(t @ (y[m] - a["y0"]))
            mb = m & inb
            if mb.any():
                a["n_in"] += int(mb.sum())
                a["tmin_in"] = min(a["tmin_in"], float(gt[mb].min()))
                a["tmax_in"] = max(a["tmax_in"], float(gt[mb].max()))
                a["tsum_in"] += float(gt[mb].sum())
    out = {}
    for s, a in acc.items():
        n = a["n"]
        det = n * a["Stt"] - a["St"] ** 2
        vx = (n * a["Stx"] - a["St"] * a["Sx"]) / det if det > 0 else np.nan   # m/s East
        vy = (n * a["Sty"] - a["St"] * a["Sy"]) / det if det > 0 else np.nan   # m/s North
        out[str(s)] = dict(
            n=n, t_min=a["tmin"], t_max=a["tmax"], t_mid=a["tsum"] / n,
            n_in=a["n_in"],
            t_min_in=(a["tmin_in"] if a["n_in"] else np.nan),
            t_max_in=(a["tmax_in"] if a["n_in"] else np.nan),
            t_mid_in=(a["tsum_in"] / a["n_in"] if a["n_in"] else np.nan),
            vx=vx, vy=vy, heading_deg=float(np.degrees(np.arctan2(vx, vy))),
        )
    hdr = f.header
    return out, dict(gps_time_type=int(hdr.global_encoding.gps_time_type),
                     point_count=int(hdr.point_count),
                     generating_software=str(hdr.generating_software))


def build_cache(rebuild=False):
    if os.path.exists(CACHE) and not rebuild:
        return json.load(open(CACHE))
    db = {}
    for name, (laz, ddir) in SITES.items():
        bounds = json.load(open(f"{ddir}/corrections.json"))["bounds"]
        st, meta = swath_gps_stats(laz, bounds)
        db[name] = dict(laz=laz, derived=ddir, bounds=bounds, meta=meta, swaths=st)
        print(f"[cache] {name}: {len(st)} swaths from {laz}", flush=True)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    json.dump(db, open(CACHE, "w"), indent=1)
    return db


# ---------------------------------------------------------------- little stats
def ols(t, v):
    """(slope, intercept, r, n) of v on t. r is Pearson. NaN where undefined."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    n = t.size
    if n < 2 or np.ptp(t) == 0:
        return np.nan, np.nan, np.nan, n
    b, a = np.polyfit(t, v, 1)
    r = np.nan if (np.std(v) == 0) else float(np.corrcoef(t, v)[0, 1])
    return float(b), float(a), r, n


def fmt(v, p=2):
    return "nan" if (v is None or not np.isfinite(v)) else f"{v:.{p}f}"


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true",
                    help="re-read every gen1 LAZ rather than using the cached swath stats")
    ap.add_argument("--min-pts", type=int, default=2000,
                    help="secondary disclosure only: also report every fit with swaths below "
                         "this return count removed. Default is coreg.fit_along_track_drift's "
                         "own min_pts (src/lidar_diff_icp/coreg.py:381); NOTHING is removed "
                         "from the primary tables.")
    A = ap.parse_args()

    R = Run("are the gen1 per-swath dx/dy/dz constants samples of one continuous mission-time "
            "drift, or per-flight-line quantities that gps_time merely correlates with?")
    for name, (laz, ddir) in SITES.items():
        R.input(laz, role=f"{name}: gen1 (2008-2012 MN) cloud align_swaths was run on -- "
                          "point_source_id + gps_time per return")
        R.input(f"{ddir}/corrections.json",
                role=f"{name}: saved per_swath_internal_alignment_dxdydz_m (m) and grid bounds")
    R.param("min_pts", A.min_pts, src="repo",
            why="coreg.fit_along_track_drift default; used for a SECONDARY disclosure row only")
    R.param("swath_reference", "min(point_source_id)", src="repo",
            why="pipeline.difference_dem calls align_swaths(ref=int(ps8.min())) -- that swath's "
                "(dx,dy,dz) is pinned to zero and is NOT an observation")
    R.column("site", "study site = one gen1 tile/merge processed by scripts/run_all_sites.py")
    R.column("swath", "point_source_id (flight line) in the gen1 file")
    R.column("n_file", "returns of that swath in the whole gen1 file (the extent align_swaths saw)")
    R.column("n_grid", "returns of that swath inside the analysis bounds (what the drift fit saw)")
    R.column("t_mid", "mean gps_time of the swath over the whole file, s (raw file units)")
    R.column("span_s", "t_max - t_min of the swath over the whole file, s")
    R.column("gap_s", "elapsed gps_time from the previous swath's t_mid, s (time order)")
    R.column("dx_m", "saved per-swath alignment East shift, m (added to gen1 x)")
    R.column("dy_m", "saved per-swath alignment North shift, m (added to gen1 y)")
    R.column("dz_mm", "saved per-swath alignment vertical shift, mm (added to gen1 z)")
    R.column("hdg_deg", "flight heading from a least-squares fit of x(t), y(t), deg clockwise from North")
    R.column("axis", "which of dz / dx / dy / across-track / along-track is being regressed")
    R.column("nsw", "number of swaths entering the fit (reference swath excluded unless stated)")
    R.column("r", "Pearson correlation of the axis against gps_time across swaths")
    R.column("slope_per_h", "ordinary-least-squares slope against gps_time, mm/hour (dz) or m/hour (dx,dy)")
    R.column("r_incl_ref", "same Pearson r with the forced-zero reference swath INCLUDED")
    R.column("R2_time", "R-squared of a 2-parameter fit linear in elapsed gps_time")
    R.column("R2_ordinal", "R-squared of a 2-parameter fit linear in flight-line ordinal (time order)")
    R.column("cv_per_h", "coefficient of variation of the consecutive-pair step expressed PER HOUR")
    R.column("cv_per_line", "coefficient of variation of the same steps expressed PER LINE")
    R.column("pair", "consecutive pair of swaths in time order")
    R.column("d_dz_mm", "dz of the later swath minus dz of the earlier, mm")
    R.column("rate_mm_h", "d_dz_mm divided by the elapsed gap, mm/hour")
    R.column("pred_mm", "step the site's own fitted continuous rate predicts over that gap, mm")
    R.column("dow", "GPS day of week of the swath midpoint (0 = Sunday); only defined for week-seconds files")
    R.column("date_utc", "UTC date from adjusted-standard gps_time; blank where the file stores week-seconds")
    R.column("t_mid_elba", "elba mean gps_time for that swath, in-grid, s")
    R.column("t_mid_ext", "elbaext mean gps_time for the SAME swath, in-grid, s")
    R.column("dt_s", "t_mid_ext - t_mid_elba, s")
    R.column("disagree_mm", "elbaext dz minus elba dz for the same swath, both re-referenced to swath 135, mm")
    R.column("need_mm_s", "mm/s of drift that would be required for dt_s to explain disagree_mm")
    R.column("gps_time_type", "LAS global encoding bit 0: 0 = GPS week seconds, 1 = adjusted standard GPS time")
    R.banner()

    db = build_cache(A.rebuild)

    # ------------------------------------------------------------ 1. per-swath table
    print("\n### 1. per-swath table (time order within each site)\n")
    print("  aliases sharing another site's alignment verbatim (not independent evidence): "
          + ", ".join(f"{k}->{v}" for k, v in ALIASES.items()))
    print("  gps_time units per site: " + ", ".join(
        f"{k}={'adjusted-standard' if db[k]['meta']['gps_time_type'] else 'week-seconds'}"
        for k in SITES))
    rows = []
    per_site = {}
    for name in SITES:
        d = db[name]
        corr = json.load(open(f"{d['derived']}/corrections.json"))["per_swath_internal_alignment_dxdydz_m"]
        rec = []
        for s, (dx, dy, dz) in corr.items():
            st = d["swaths"].get(s)
            if st is None:
                print(f"  WARNING {name} swath {s} has a saved correction but no returns in the gen1 file")
                continue
            rec.append(dict(swath=int(s), dx=dx, dy=dy, dz_mm=1e3 * dz, **st))
        rec.sort(key=lambda r: r["t_mid"])
        ref = min(int(s) for s in corr)          # align_swaths(ref=min(point_source_id))
        prev = None
        for r in rec:
            r["gap_s"] = np.nan if prev is None else r["t_mid"] - prev
            prev = r["t_mid"]
            r["is_ref"] = (r["swath"] == ref)
            rows.append([name, f"{r['swath']}{'*' if r['is_ref'] else ''}", r["n"], r["n_in"],
                         f"{r['t_mid']:.1f}", fmt(r["t_max"] - r["t_min"], 1), fmt(r["gap_s"], 0),
                         fmt(r["dx"], 4), fmt(r["dy"], 4), fmt(r["dz_mm"], 1), fmt(r["heading_deg"], 1)])
        per_site[name] = (rec, ref)
    R.table(["site", "swath", "n_file", "n_grid", "t_mid", "span_s", "gap_s",
             "dx_m", "dy_m", "dz_mm", "hdg_deg"], rows)
    print("  * = reference swath, pinned to (0,0,0) by align_swaths; not an observation.")

    # ------------------------------------------------------------ 1b. mission calendar
    print("\n### 1b. when each site was flown, from the file's own gps_time encoding\n")
    import datetime as _dt
    GPS_EPOCH = _dt.datetime(1980, 1, 6, tzinfo=_dt.timezone.utc)
    rows = []
    for name in SITES:
        rec, _ = per_site[name]
        typ = db[name]["meta"]["gps_time_type"]
        t0 = min(r["t_min"] for r in rec); t1 = max(r["t_max"] for r in rec)
        if typ == 1:      # adjusted standard GPS time: seconds since GPS epoch, minus 1e9
            d0 = GPS_EPOCH + _dt.timedelta(seconds=t0 + 1e9)
            d1 = GPS_EPOCH + _dt.timedelta(seconds=t1 + 1e9)
            dow, date = f"{d0.weekday():d}(iso)", f"{d0:%Y-%m-%d %H:%M}Z -> {d1:%H:%M}Z"
        else:             # GPS week seconds: day of week is known, the WEEK is not recorded
            dow = "/".join(sorted({str(int(r["t_mid"] // 86400)) for r in rec}))
            h0 = _dt.timedelta(seconds=t0 % 86400); h1 = _dt.timedelta(seconds=t1 % 86400)
            date = f"week unknown; {str(h0)[:8]} -> {str(h1)[:8]} GPS-time-of-day"
        rows.append([name, len(rec), typ, dow, date])
    R.table(["site", "nsw", "gps_time_type", "dow", "date_utc"], rows)
    print("  Leap seconds (15-16 s in 2011-2012) are not removed, so the adjusted-standard")
    print("  timestamps are UTC to within ~16 s -- irrelevant for dating the flight.")
    print("  battlecreek straddles a day boundary: lines 1012-1014 sit on GPS day 0 and")
    print("  1101-1102 on day 1, a 70,974 s (19.7 h) gap. That is a genuine mission break and")
    print("  the single largest lever this dataset has on a continuous-time model.")
    print("  elba, elbaext and whitewater all fall on GPS day 2 with disjoint windows")
    print("  (whitewater 226,000-229,475 s; elba/elbaext 237,350-249,711 s). IF they are the")
    print("  same GPS week they are one flight day and one F(t) must carry both; the week is")
    print("  NOT recorded in these LAS 1.1 files, so that remains conditional, not established.")

    # ------------------------------------------------------------ 2/3. slopes
    print("\n### 2-3. drift rate per site, per axis, reference swath EXCLUDED\n")
    rows = []
    for name in SITES:
        rec, ref = per_site[name]
        obs = [r for r in rec if not r["is_ref"]]
        allr = rec
        t = np.array([r["t_mid"] for r in obs])
        ta = np.array([r["t_mid"] for r in allr])
        # across/along-track: the TRACK AXIS, from a circular mean of the doubled headings.
        # Alternating flight lines run +1 deg and -179 deg (there and back); a plain mean of
        # those is meaningless, but doubling folds heading modulo 180 so the two directions of
        # the same axis agree, and halving recovers the axis.
        h2 = np.radians(2.0 * np.array([r["heading_deg"] for r in rec]))
        hd = 0.5 * np.arctan2(np.nanmean(np.sin(h2)), np.nanmean(np.cos(h2)))
        along = np.array([np.sin(hd), np.cos(hd)])      # (East, North)
        across = np.array([np.cos(hd), -np.sin(hd)])
        series = {
            "dz (mm)":       ([r["dz_mm"] for r in obs],  [r["dz_mm"] for r in allr], 1.0),
            "dx (m)":        ([r["dx"] for r in obs],     [r["dx"] for r in allr],    1.0),
            "dy (m)":        ([r["dy"] for r in obs],     [r["dy"] for r in allr],    1.0),
            "across (m)":    ([r["dx"] * across[0] + r["dy"] * across[1] for r in obs],
                              [r["dx"] * across[0] + r["dy"] * across[1] for r in allr], 1.0),
            "along (m)":     ([r["dx"] * along[0] + r["dy"] * along[1] for r in obs],
                              [r["dx"] * along[0] + r["dy"] * along[1] for r in allr], 1.0),
        }
        for ax, (v, va, _) in series.items():
            b, _, r_, n = ols(t, v)
            _, _, ra, _ = ols(ta, va)
            rows.append([name, ax, n, fmt(r_, 3), fmt(3600 * b, 3), fmt(ra, 3)])
    R.table(["site", "axis", "nsw", "r", "slope_per_h", "r_incl_ref"], rows)
    print("  Sign convention: the shift is ADDED to gen1. A site whose reference swath is the")
    print("  LAST in time (elba, elbaext, whitewater, carlton) has its zero at the late end;")
    print("  re-referencing changes every value but not r or the slope, which are gauge-free.")

    # ------------------------------------------------------------ 4. the test
    print("\n### 4. THE FALSIFIABLE TEST -- does the step between two lines scale with the")
    print("###    ELAPSED TIME between them (continuous F(t)) or with the LINE COUNT?\n")
    print("    pred_mm is what the site's OWN best-fit continuous rate (section 2, reference")
    print("    swath excluded) predicts for that gap. Under one continuous F(t) every pair")
    print("    should sit near pred_mm; under a per-line quantity the steps are alike in SIZE")
    print("    regardless of gap.\n")
    rows = []
    for name in SITES:
        rec, ref = per_site[name]
        obs = [r for r in rec if not r["is_ref"]]
        rate, _, _, _ = ols([r["t_mid"] for r in obs], [r["dz_mm"] for r in obs])
        for a, b in zip(rec[:-1], rec[1:]):
            gap = b["t_mid"] - a["t_mid"]
            d = b["dz_mm"] - a["dz_mm"]
            rows.append([name, f"{a['swath']}->{b['swath']}", fmt(gap, 0),
                         fmt(d, 1), fmt(3600 * d / gap, 1), fmt(rate * gap, 1)])
    R.table(["site", "pair", "gap_s", "d_dz_mm", "rate_mm_h", "pred_mm"], rows)

    print("\n  Consistency of the step, expressed per hour vs per line (lower = the better")
    print("  description). CV = std/|mean| over the site's consecutive pairs.\n")
    rows = []
    for name in SITES:
        rec, _ = per_site[name]
        for ax, key in (("dz (mm)", "dz_mm"), ("dx (m)", "dx"), ("dy (m)", "dy")):
            steps, gaps = [], []
            for a, b in zip(rec[:-1], rec[1:]):
                steps.append(b[key] - a[key]); gaps.append(b["t_mid"] - a["t_mid"])
            steps = np.array(steps); gaps = np.array(gaps)
            if steps.size < 2:
                continue
            per_h = 3600 * steps / gaps
            cvh = np.std(per_h) / abs(np.mean(per_h)) if np.mean(per_h) else np.inf
            cvl = np.std(steps) / abs(np.mean(steps)) if np.mean(steps) else np.inf
            rows.append([name, ax, steps.size, fmt(cvh, 2), fmt(cvl, 2)])
    R.table(["site", "axis", "nsw", "cv_per_h", "cv_per_line"], rows)

    # ------------------------------------------------------------ 5. model comparison
    print("\n### 5. two-parameter linear-in-TIME vs two-parameter linear-in-ORDINAL\n")
    rows = []
    for name in SITES:
        rec, _ = per_site[name]
        obs = [r for r in rec if not r["is_ref"]]
        t = np.array([r["t_mid"] for r in obs])
        k = np.arange(len(obs), dtype=float)
        for ax, key in (("dz (mm)", "dz_mm"), ("dx (m)", "dx"), ("dy (m)", "dy")):
            v = np.array([r[key] for r in obs])
            if v.size < 3 or np.std(v) == 0:
                rows.append([name, ax, v.size, "n/a", "n/a"]); continue
            def r2(u):
                b, a, _, _ = ols(u, v)
                return 1 - np.sum((v - (a + b * u)) ** 2) / np.sum((v - v.mean()) ** 2)
            rows.append([name, ax, v.size, fmt(r2(t), 3), fmt(r2(k), 3)])
    R.table(["site", "axis", "nsw", "R2_time", "R2_ordinal"], rows)

    # secondary disclosure: the same with tiny swaths removed
    print(f"\n  Secondary disclosure -- same table with swaths of fewer than {A.min_pts} returns")
    print("  removed (coreg.fit_along_track_drift's own min_pts; nothing above uses it):\n")
    rows = []
    for name in SITES:
        rec, _ = per_site[name]
        obs = [r for r in rec if not r["is_ref"] and r["n"] >= A.min_pts]
        dropped = [r["swath"] for r in rec if r["n"] < A.min_pts]
        if dropped:
            print(f"    {name}: drops swath(s) {dropped}")
        t = np.array([r["t_mid"] for r in obs])
        for ax, key in (("dz (mm)", "dz_mm"), ("dx (m)", "dx"), ("dy (m)", "dy")):
            v = np.array([r[key] for r in obs])
            if v.size < 3 or np.std(v) == 0:
                rows.append([name, ax, v.size, "n/a", "n/a"]); continue
            k = np.arange(v.size, dtype=float)
            def r2(u):
                b, a, _, _ = ols(u, v)
                return 1 - np.sum((v - (a + b * u)) ** 2) / np.sum((v - v.mean()) ** 2)
            rows.append([name, ax, v.size, fmt(r2(t), 3), fmt(r2(k), 3)])
    R.table(["site", "axis", "nsw", "R2_time", "R2_ordinal"], rows)

    # ------------------------------------------------------------ 6. elba vs elbaext
    print("\n### 6. the elba / elbaext anomaly: is the disagreement a gps_time-midpoint effect?\n")
    ce = json.load(open("data/derived/elba/corrections.json"))["per_swath_internal_alignment_dxdydz_m"]
    cx = json.load(open("data/derived/elbaext/corrections.json"))["per_swath_internal_alignment_dxdydz_m"]
    rows = []
    shared = sorted(set(ce) & set(cx), key=int)
    base = "135"
    for s in shared:
        dz_e = 1e3 * (ce[s][2] - ce[base][2])
        dz_x = 1e3 * (cx[s][2] - cx[base][2])
        te = db["elba"]["swaths"][s]["t_mid_in"]; tx = db["elbaext"]["swaths"][s]["t_mid_in"]
        dt = tx - te
        dis = dz_x - dz_e
        rows.append([s, f"{te:.1f}", f"{tx:.1f}", fmt(dt, 2), fmt(dis, 1),
                     fmt(dis / dt, 2) if dt else "inf"])
    R.table(["swath", "t_mid_elba", "t_mid_ext", "dt_s", "disagree_mm", "need_mm_s"], rows)
    rec, _ = per_site["elba"]
    obs = [r for r in rec if not r["is_ref"]]
    b, _, _, _ = ols([r["t_mid"] for r in obs], [r["dz_mm"] for r in obs])
    print(f"\n  elba's own fitted drift rate is {3600*b:.1f} mm/hour = {b:.5f} mm/s.")
    print("  Compare that with the need_mm_s column: the ratio is the factor by which the")
    print("  midpoint-difference explanation falls short.")

    R.done(headline="mission-time drift: correlation vs gap-structure test across 7 gen1 sites")


if __name__ == "__main__":
    main()
