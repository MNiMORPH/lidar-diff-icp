"""Parse the MnGeo per-county validation reports into a bundled checkpoint CSV.

The 2008 SE-Minnesota lidar (MN DNR / AeroMetric -- our "gen1") was validated by MnDNR
against its own ground control. ``analysis/ridgelines/ABSOLUTE_ELEVATION_REFS.md`` §1b
called those coordinates a dead end; they are not. The per-county **validation reports**
in the open MnGeo lidar tree carry the full checkpoint tables, one row per mark:

    Name, Control X, Control Y, Control Z, Surface Z, Error, Z-Diff Squared, Absolute Error

    https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/<county>/
        <County>_county_validation_report.pdf     (capitalisation varies by county)

**Sign convention, pinned by arithmetic, not by the column name.** On every parsed row of
every report ``Error == Control Z - Surface Z`` exactly (max |residual| 0.000000 m over
1147 rows); the other order misses by up to 1.08 m. So a **negative** ``Error`` means the
delivered 2008 surface reads **above** the surveyed mark. This is the *same* sign family
as ``groundtruth.tie``: ``tie = surveyed - z_lidar``. ``--check`` re-runs that test and
reproduces each report's own published RMSE.

``Surface Z`` is the **delivered 2008 DNR DEM**, not our reconstruction from the point
cloud, so the ``Error`` column is not our residual; it is carried through as a column for
comparison. ``Control Z`` is the useful number.

Land cover is encoded in the point name -- L1O open, L2T tall weeds/crops, L3B brush,
L4F forest, L5U urban -- and ``pdftotext`` renders "L1O" as "L10" in some reports, so the
class is taken from the digit.

Usage
-----
    python analysis/groundtruth/parse_mndnr_2008_control.py --pdf-dir <dir> [--check]
        [--out src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess

# The eight counties the 2008 SE-MN project covers, i.e. the ones for which the dataset
# metadata (lidar_semn2008.html) publishes a per-county RMSE and sample count.  Freeborn
# has a validation report in the same tree but is NOT in that list and is not gen1.
SEMN_COUNTIES = ("dodge", "fillmore", "houston", "mower", "olmsted", "steele",
                 "wabasha", "winona")

# InPort 68818 / lidar_semn2008.html, "MnDNR's Tests": county -> (RMSE m, n).
PUBLISHED = {"dodge": (0.129, 121), "fillmore": (0.155, 128), "houston": (0.110, 134),
             "mower": (0.161, 115), "olmsted": (0.117, 125), "steele": (0.125, 137),
             "wabasha": (0.106, 97), "winona": (0.161, 176)}

COVER = {"1": "L1O", "2": "L2T", "3": "L3B", "4": "L4F", "5": "L5U"}
COVER_NAME = {"L1O": "open terrain", "L2T": "tall weeds and crops",
              "L3B": "brush lands and low trees", "L4F": "forested",
              "L5U": "urban", "other": "unclassed (FSA photo target etc.)"}

ROW = re.compile(r"^\s*(?P<name>\S.*?)\s+(?P<x>\d{6}\.\d+)\s+(?P<y>\d{7}\.\d+)\s+"
                 r"(?P<cz>-?\d+\.\d+)\s+(?P<sz>-?\d+\.\d+)\s+(?P<err>-?\d+\.\d+)\s+"
                 r"(?P<sq>-?\d+\.?\d*)\s+(?P<abs>\d+\.\d+)")
SUMMARY = re.compile(r"^\s*(?:Z-)?RMSE\s+(?P<v>\d\.\d+)\s*$")


def cover_of(name: str) -> str:
    m = re.match(r"^L\s?([1-5])\s?[O0TBFUotbfu]", name)
    return COVER[m.group(1)] if m else "other"


def parse_report(pdf: str, *, all_tables: bool = False):
    """Rows of one county report, plus the RMSE it prints under its first table.

    Each report leads with an **overall** table (every mark) and then repeats the marks
    broken out per land-cover class.  Stopping at the first printed RMSE takes the
    overall table only, which is what reproduces the published per-county figure;
    ``all_tables=True`` keeps everything, for the sign test.
    """
    txt = subprocess.run(["pdftotext", "-layout", pdf, "-"],
                         capture_output=True, text=True, check=True).stdout
    rows, printed = [], None
    for line in txt.splitlines():
        s = SUMMARY.match(line.strip())
        if s and rows and not all_tables:
            printed = float(s.group("v"))
            break
        m = ROW.match(line)
        if m:
            d = m.groupdict()
            rows.append(dict(name=d["name"].strip(), x=float(d["x"]), y=float(d["y"]),
                             cz=float(d["cz"]), sz=float(d["sz"]), err=float(d["err"]),
                             sq=float(d["sq"]), aerr=float(d["abs"])))
    seen, out = set(), []
    for r in rows:                                   # a mark can repeat across tables
        k = (r["name"], r["x"], r["y"])
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out, printed


def find_pdfs(pdf_dir: str) -> dict:
    got = {}
    for f in sorted(os.listdir(pdf_dir)):
        if not f.lower().endswith(".pdf"):
            continue
        for c in SEMN_COUNTIES:
            if f.lower().startswith(c):
                got[c] = os.path.join(pdf_dir, f)
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf-dir", required=True,
                    help="directory holding <county>*_validation_report.pdf")
    ap.add_argument("--out", default="src/lidar_diff_icp/groundtruth/data/"
                                     "mn_dnr_2008_control_semn.csv")
    ap.add_argument("--check", action="store_true",
                    help="run the sign test and the published-RMSE reproduction")
    A = ap.parse_args()

    pdfs = find_pdfs(A.pdf_dir)
    missing = [c for c in SEMN_COUNTIES if c not in pdfs]
    if missing:
        raise SystemExit(f"no report found for {missing} in {A.pdf_dir}")

    if A.check:
        import numpy as np
        print("=== sign test: Error == ControlZ - SurfaceZ ?  (every table, every county)")
        n_tot = n_cs = n_sc = 0
        rcs = rsc = 0.0
        for c in SEMN_COUNTIES:
            rows, _ = parse_report(pdfs[c], all_tables=True)
            e = np.array([r["err"] for r in rows])
            cz = np.array([r["cz"] for r in rows])
            sz = np.array([r["sz"] for r in rows])
            n_tot += len(rows)
            n_cs += int((np.abs(e - (cz - sz)) < 5e-4).sum())
            n_sc += int((np.abs(e - (sz - cz)) < 5e-4).sum())
            rcs = max(rcs, float(np.abs(e - (cz - sz)).max()))
            rsc = max(rsc, float(np.abs(e - (sz - cz)).max()))
        print(f"  Control - Surface : {n_cs}/{n_tot} rows agree, max |resid| {rcs:.6f} m")
        print(f"  Surface - Control : {n_sc}/{n_tot} rows agree, max |resid| {rsc:.6f} m")
        print("\n=== published-RMSE reproduction (overall table of each report)")
        print(f"{'county':10s} {'n':>5s} {'n_pub':>6s} {'RMSE':>7s} {'in-report':>10s} "
              f"{'published':>10s}")
        for c in SEMN_COUNTIES:
            rows, printed = parse_report(pdfs[c])
            e = np.array([r["err"] for r in rows])
            pub, npub = PUBLISHED[c]
            print(f"{c:10s} {len(rows):5d} {npub:6d} {np.sqrt((e ** 2).mean()):7.4f} "
                  f"{(printed if printed else float('nan')):10.3f} {pub:10.3f}")
        print()

    os.makedirs(os.path.dirname(A.out), exist_ok=True)
    fields = ["point_id", "point_type", "easting", "northing", "elevation",
              "elevation_units", "horizontal_crs", "vertical_datum", "geoid_model",
              "project_id", "collected", "source", "verified",
              "county", "dnr_surface_z_m", "dnr_error_m"]
    n = 0
    with open(A.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fields)
        w.writeheader()
        for c in SEMN_COUNTIES:
            rows, _ = parse_report(pdfs[c])
            for r in rows:
                w.writerow(dict(
                    point_id=r["name"], point_type=cover_of(r["name"]),
                    easting=r["x"], northing=r["y"], elevation=r["cz"],
                    elevation_units="m", horizontal_crs="EPSG:26915",
                    vertical_datum="NAVD88", geoid_model="GEOID03",
                    project_id="lidar_semn2008", collected="2008",
                    source=f"{os.path.basename(pdfs[c])} table 'Control Z'",
                    verified="datum is a DATASET-level assertion (lidar_semn2008.html: "
                             "'Vertical datum: NAVD88 (Geoid03)'); the validation "
                             "reports themselves state no datum and no geoid",
                    county=c, dnr_surface_z_m=r["sz"], dnr_error_m=r["err"]))
                n += 1
    print(f"wrote {n} control points to {A.out}")


if __name__ == "__main__":
    main()
