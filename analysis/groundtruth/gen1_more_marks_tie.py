"""gen1's datum on ITS OWN flight lines 133-138, from every 2008 control mark those
lines actually illuminated.

What is new here against ``gen1_own_control_tie.py`` / the first screening pass
--------------------------------------------------------------------------------
1. **The line assignment comes from the data, not from geometry.**  The first pass
   assigned each mark to the flight line whose fitted centreline was nearest.  Line
   spacing here is ~960 m and the measured swath half-width is ~775 m, so the swaths
   OVERLAP: a mark can be, and often is, inside two lines' swaths at once, and the
   nearest centreline is then simply the wrong label for at least one of the returns.
   Here the line is read off the ``point_source_id`` of the ground returns AT the mark,
   and a separate tie is estimated from each line's own returns
   (``tie.estimate_tie(..., line=N)``, which is what that argument is for).
   The centreline model (``gen1_line_tracks.py``) is used ONLY to decide which tiles to
   fetch.

2. **Marks seen by two lines give a within-mark line difference.**  That difference
   cancels the site term -- the local relief, the vegetation, the survey error, the mark
   itself -- and is therefore the sharp measurement of the RELATIVE line offsets, which
   is the quantity ``corrections_geoid.json``'s ``per_swath_internal_alignment`` holds.

3. **The catchment is widened and the widening is measured**, not assumed: every
   open/urban mark within 2 km of any of the six centrelines is read, and the answer to
   "how far out is it worth going" is the measured fraction of marks at each distance
   that actually carry 133-138 returns.

Sign convention, unchanged and the same on both sides:
``tie = surveyed - z_lidar``, positive = the lidar reads BELOW the mark.
``geoid_shift_m = 0``: the 2008 control is NAVD88(GEOID03) and so is the raw gen1
cloud, so this is not an approximation, it is the same frame -- the geoid cancels and
there is nothing to convert.  ``swath_shift_m = (0,0,0)``: no gen2-derived lateral or
vertical term belongs in a gen1-against-its-own-control comparison.

Ground source: VENDOR class 2.  CSF is ~460 s/tile and the 2026-08-26 report measured
the ground-source dependence of this comparison at 6.5 mm median absolute over 16 marks.

Usage
-----
    SCRATCH=... env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/groundtruth/gen1_more_marks_tie.py [--radius-m 2000]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lidar_diff_icp.groundtruth import tie as T                       # noqa: E402
from lidar_diff_icp.groundtruth.checkpoints import Checkpoint         # noqa: E402

SCRATCH = os.environ.get("SCRATCH", ".")
CONTROL_CSV = "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv"
LINES = (133, 134, 135, 136, 137, 138)
RES_M = 5.0                 # corrections_geoid.json
CROP_HALFWIDTH_M = 300.0    # carried unchanged from elba_absolute_tie.py so the runs are one method
SITING_RADIUS_M = 5.0       # ADDITIONAL_GROUND_CONTROL.md section 7.1


def cross_track(track, x, y):
    """Signed perpendicular distance from a whole-line straight fit, metres."""
    w = track["w"]
    xp = np.polyval(w, y - track["ym"]) + track["xm"]
    return (x - xp) / math.hypot(1.0, w[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius-m", type=float, default=2000.0,
                    help="candidate catchment: |cross-track| to any of lines 133-138")
    ap.add_argument("--covers", default="L1O,L5U",
                    help="MnDNR land-cover classes to read; 'all' reads every class")
    ap.add_argument("--tag", default="",
                    help="suffix for the output CSV names")
    a = ap.parse_args()
    covers = None if a.covers.strip().lower() == "all" else a.covers.split(",")

    tracks = json.load(open(f"{SCRATCH}/line_tracks.json"))["whole"]
    inv = pd.read_csv(f"{SCRATCH}/control_inventory.csv")
    ctrl = pd.read_csv(CONTROL_CSV).drop_duplicates("point_id").set_index("point_id")

    x = inv.easting.to_numpy(); y = inv.northing.to_numpy()
    C = {}
    for ln in LINES:
        t = tracks[str(ln)]
        c = np.abs(cross_track(t, x, y))
        C[ln] = np.where((y >= t["nmin"]) & (y <= t["nmax"]), c, np.inf)
    A = np.stack([C[l] for l in LINES], 1)
    inv["d_min"] = A.min(1)
    inv["line_nearest"] = [LINES[i] for i in A.argmin(1)]
    for ln in LINES:
        inv[f"d_{ln}"] = C[ln]

    sel_cov = np.ones(len(inv), bool) if covers is None else inv.cls.isin(covers).to_numpy()
    cand = inv[sel_cov & (inv.d_min <= a.radius_m).to_numpy()].copy()
    cand = cand.sort_values("d_min").reset_index(drop=True)
    print(f"candidates: {len(cand)} marks (covers={a.covers}) within {a.radius_m:.0f} m "
          f"of a 133-138 centreline")

    rows, per_line_rows = [], []
    for k, r in cand.iterrows():
        p = f"data/before/{r.tile}.laz"
        base = dict(point=r.point_id, cls=r.cls, county=r.county, km=float(r.km),
                    tile=r.tile, d_min=float(r.d_min), line_nearest=int(r.line_nearest),
                    **{f"d_{ln}": float(r[f"d_{ln}"]) for ln in LINES})
        if not os.path.exists(p):
            rows.append({**base, "status": "tile not on disk"}); continue
        c = ctrl.loc[r.point_id]
        cp = Checkpoint(point_id=r.point_id, point_type=str(c.point_type),
                        easting=float(c.easting), northing=float(c.northing),
                        elevation=float(c.elevation), elevation_units="m",
                        horizontal_crs="EPSG:26915", vertical_datum="NAVD88",
                        geoid_model="GEOID03", project_id="lidar_semn2008",
                        collected="2008", source=str(c.source), verified=str(c.verified))
        try:
            g = T.vendor_ground_near(p, cp.easting, cp.northing, CROP_HALFWIDTH_M)
        except Exception as ex:
            rows.append({**base, "status": f"FAIL {type(ex).__name__}: {ex}"[:70]}); continue

        report_R = 1.5 * RES_M
        d = np.hypot(g.x - cp.easting, g.y - cp.northing)
        sel = d <= report_R
        psid, cnt = np.unique(g.point_source_id[sel].astype(int), return_counts=True)
        comp = {int(a_): int(b_) for a_, b_ in zip(psid, cnt)}
        ours = {k_: v for k_, v in comp.items() if k_ in LINES}
        # siting screen, section 7.1
        zs = g.z[d <= SITING_RADIUS_M]
        spread = 1000.0 * float(np.percentile(zs, 95) - np.percentile(zs, 5)) if zs.size >= 5 else np.nan
        pct = 100.0 * float((zs < cp.elevation_m).mean()) if zs.size >= 5 else np.nan

        est = T.estimate_tie(cp, g, res=RES_M, swath_shift_m=(0., 0., 0.), geoid_shift_m=0.0)
        rr = {e.radius_m: e for e in est.curve}[est.report_radius_m]
        row = {**base, "status": "ok", "n_report": rr.n,
               "psid_at_mark": ",".join(f"{k_}:{v}" for k_, v in sorted(comp.items())),
               "n_psid": len(comp), "lines_ours": ",".join(str(k_) for k_ in sorted(ours)),
               "n_lines_ours": len(ours),
               "line_dominant": (max(ours, key=ours.get) if ours else -1),
               "spread_mm": spread, "ctrl_pct": pct,
               "slope_deg": rr.slope_deg, "relief_mm": rr.relief_mm,
               "ladder_mm": est.radius_spread_mm, "sigma_mm": est.sigma_mm,
               "tie_mixed_mm": est.tie_mm, "dnr_mm": 1000.0 * float(c.dnr_error_m)}
        for ln in sorted(ours):
            e1 = T.estimate_tie(cp, g, line=ln, res=RES_M, swath_shift_m=(0., 0., 0.),
                                geoid_shift_m=0.0)
            r1 = {e.radius_m: e for e in e1.curve}[e1.report_radius_m]
            row[f"tie_{ln}_mm"] = e1.tie_mm
            row[f"n_{ln}"] = r1.n
            row[f"sig_{ln}_mm"] = e1.sigma_mm
            per_line_rows.append(dict(point=r.point_id, cls=r.cls, county=r.county,
                                      km=float(r.km), tile=r.tile, line=ln,
                                      d_line=float(r[f"d_{ln}"]), d_min=float(r.d_min),
                                      n=r1.n, n_at_mark=ours[ln],
                                      tie_mm=e1.tie_mm, sigma_mm=e1.sigma_mm,
                                      ladder_mm=e1.radius_spread_mm,
                                      spread_mm=spread, ctrl_pct=pct,
                                      slope_deg=r1.slope_deg,
                                      tie_mixed_mm=est.tie_mm,
                                      line_nearest=int(r.line_nearest),
                                      dnr_mm=1000.0 * float(c.dnr_error_m)))
        rows.append(row)
        print(f"  [{k+1:3d}/{len(cand)}] {r.point_id:26s} d_min={r.d_min:7.1f} "
              f"nearest={r.line_nearest} ours={sorted(ours)} psids={sorted(comp)}", flush=True)

    marks = pd.DataFrame(rows)
    marks.to_csv(f"{SCRATCH}/more_marks{a.tag}.csv", index=False)
    pl = pd.DataFrame(per_line_rows)
    pl.to_csv(f"{SCRATCH}/more_marks_perline{a.tag}.csv", index=False)
    print(f"\nwrote more_marks{a.tag}.csv ({len(marks)} marks) and "
          f"more_marks_perline{a.tag}.csv ({len(pl)} mark-line ties)")


if __name__ == "__main__":
    main()
