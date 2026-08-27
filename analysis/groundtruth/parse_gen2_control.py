"""Parse the 2021 MN_SE_Driftless (gen2) ground control into a bundled checkpoint CSV.

**The 2021 project publishes per-point residuals.** They are not in the survey report
and not in the contractor's checkpoint shapefile; they are in the *USGS* vertical-accuracy
shapefiles written by the NGTOC "VATool":

    .../metadata/MN_SE_Driftless_2021_B21/Vertical_Accuracy/USGS/
        USGS_MN_SE_Driftless_2021_B21_QL{0,1}.{shp,shx,dbf,prj}

whose ``.dbf`` carries, one row per checkpoint:

    proj_name, srcChkptId, X, Y, Z, Lndcover, DEMz, zdiff, zdiffSq, LAZz, LAZzdiff,
    LAZzdiffSq

``Z`` is the surveyed orthometric height; ``DEMz`` is the delivered OPR DEM read at the
mark and ``LAZz`` the delivered classified point cloud read at the mark; ``zdiff`` and
``LAZzdiff`` are the residuals. This is the gen2 analogue of the MnGeo county validation
reports' ``Control Z / Surface Z / Error`` that gave gen1 its 963-point residual field --
with two delivered surfaces instead of one.

**Sign, pinned by arithmetic and not by the column name.** ``zdiff == Z - DEMz`` and
``LAZzdiff == Z - LAZz`` on every row of both blocks. So a **negative** residual means the
delivered 2021 surface reads **above** the surveyed mark, which is the same sign family as
``groundtruth.tie``'s ``tie = surveyed - z_lidar``: positive = the surface reads LOW.
``--check`` re-runs that test and reproduces all eight published aggregate statistics
(NVA RMSEz and VVA 95th percentile, for the DEM and the LAZ, for each of QL0 and QL1)
from the parsed rows.

**Which points calibrated the data, and which were held out.** The vendor FGDC metadata
(``MN_SEDriftless_2_2021_Classified_Point_Cloud_Metadata.xml``, Ground Conditions) states
that the LiDAR Control Points were used to "calibrate the lidar to known ground
locations", while the accuracy checkpoints "were not used to calibrate or post process
the data". The 143 LCPs are therefore **not independent** of gen2's own vertical
adjustment and must not enter an accuracy statement about gen2; the 227 NVA and 164 VVA
checkpoints are. The CSV records this in the ``role`` column, and only the NVA/VVA rows
carry residuals -- the LCPs were never tested by the VATool.

Sources, and what each one supplies
-----------------------------------
1. ``MN_SE_Driftless_2021_B21_Ground_Control_Survey_Report.pdf`` (Woolpert, January 2023),
   §2.2 / §2.3 / §2.4 -- coordinate tables for **533** points, the *only* public source of
   the 143 LCPs. Its §1.3 text says 143 + 227 + 164 = 534; the tables hold 143 + 227 +
   163. The 164th VVA (``3000_2021_MN``) is missing from the tables and is recovered from
   the USGS shapefile, which completes the set at 534.
2. ``USGS/USGS_..._QL{0,1}.dbf`` -- the per-point residuals, 390 unique checkpoints
   (5 of them tested against both blocks, so 395 rows).
3. ``contractor_provided/MN_Driftless_NVA_VVA_UTM15_QL{0,1}.dbf`` -- per-mark
   ``geoid`` = "Geoid 18" and the per-mark collection date. The per-mark geoid attribute
   is the reason a gen2 tie needs no geoid conversion, and it is asserted here per mark;
   for the LCPs, which appear in no shapefile, the assertion is the report's own
   per-table header ("Geoid Model: Geoid18").

Note the shapefile's ``source_v_1`` field reads "US Feet". That label is wrong and was
already disproved (see ``groundtruth/data/README.md``); the values are metres, and the
report's §1.8.4 says so in words.

Usage
-----
    python analysis/groundtruth/parse_gen2_control.py \
        --report-pdf <survey_report.pdf> --usgs-dir <dir> --contractor-dir <dir> \
        [--check --tol-m 1e-6] \
        [--out src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv]

``--tol-m`` has no default: the sign test needs a tolerance and the caller states it.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import struct
import subprocess

# Point-ID / northing / easting / orthometric height / code, in that column order.
# One 2022 NVA mark carries a letter suffix (2198A_2022_MN); a regex without [A-Z]? drops
# it silently and parses 532 rows instead of 533.
ROW = re.compile(r"^\s*(?P<pid>\d{4}[A-Z]?_20\d\d_MN)\s+(?P<n>\d{7}\.\d+)\s+"
                 r"(?P<e>\d{6}\.\d+)\s+(?P<z>\d+\.\d+)\s+(?P<code>LCP|NVA|VVA)\s*$",
                 re.M)

# Ground Control Survey Report §1.3, verbatim: "143 lidar control points, 227
# non-vegetated check points, and 164 vegetated check points".
REPORT_COUNTS = {"LCP": 143, "NVA": 227, "VVA": 164}

# Published aggregates, quoted from the USGS VATool report files in the same directory.
# (block, surface) -> (NVA n, NVA RMSEz cm, VVA n, VVA 95th percentile cm), from
#   USGS_MN_SE_Driftless_2021_B21_QL{0,1}_VA.txt              -- the OPR DEM ("DEMz")
#   USGS_MN_SE_Driftless_2021_B21_QL{0,1}_las_checkpoint_report.txt -- the LAZ ("LAZz")
PUBLISHED = {
    ("QL1", "DEM"): (139, 3.51, 99, 25.48),
    ("QL1", "LAZ"): (139, 3.54, 99, 27.14),
    ("QL0", "DEM"): (91, 3.51, 66, 13.33),
    ("QL0", "LAZ"): (91, 3.55, 66, 12.64),
}

ROLE = {"LCP": "calibration", "NVA": "check", "VVA": "check"}


def read_dbf(path):
    """Minimal .dbf reader -- returns a list of dicts, numerics as float or None."""
    with open(path, "rb") as fh:
        b = fh.read()
    nrec, hlen, rlen = struct.unpack("<IHH", b[4:12])
    fields, off = [], 32
    while b[off] != 0x0D:
        fields.append((b[off:off + 11].split(b"\0")[0].decode(),
                       chr(b[off + 11]), b[off + 16]))
        off += 32
    out = []
    for i in range(nrec):
        p = hlen + i * rlen
        if b[p:p + 1] == b"*":          # deleted record
            continue
        p += 1
        r = {}
        for name, typ, ln in fields:
            v = b[p:p + ln].decode("latin-1").strip()
            p += ln
            r[name] = (float(v) if v else None) if typ == "N" else v
        out.append(r)
    return out


def report_text(pdf_path):
    out = pdf_path + ".layout.txt"
    subprocess.run(["pdftotext", "-layout", pdf_path, out], check=True)
    with open(out, errors="replace") as fh:
        return fh.read()


def parse_report(pdf_path):
    """{point_id: (easting, northing, elevation, code)} from the report's UTM tables."""
    txt = report_text(pdf_path)
    rows = {}
    for m in ROW.finditer(txt):
        rows[m.group("pid")] = (float(m.group("e")), float(m.group("n")),
                                float(m.group("z")), m.group("code"))
    return rows


def load_usgs(usgs_dir):
    """{point_id: {block: row}} from the two VATool shapefiles."""
    out = {}
    for block in ("QL0", "QL1"):
        path = os.path.join(usgs_dir,
                            f"USGS_MN_SE_Driftless_2021_B21_{block}.dbf")
        for r in read_dbf(path):
            out.setdefault(r["srcChkptId"], {})[block] = r
    return out


def load_contractor(contractor_dir):
    """{point_id: row} from the two contractor NVA/VVA shapefiles."""
    out = {}
    for block in ("QL0", "QL1"):
        path = os.path.join(contractor_dir,
                            f"MN_Driftless_NVA_VVA_UTM15_{block}.dbf")
        for r in read_dbf(path):
            out.setdefault(r["unique_ind"], {})[block] = r
    return out


def _date(v):
    if v is None:
        return ""
    s = f"{int(v):08d}"
    return f"{s[:4]}-{s[4:6]}-{s[6:]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-pdf", required=True)
    ap.add_argument("--usgs-dir", required=True)
    ap.add_argument("--contractor-dir", required=True)
    ap.add_argument("--out", default="src/lidar_diff_icp/groundtruth/data/"
                                     "mn_se_driftless_2021_control.csv")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tol-m", type=float, default=None,
                    help="tolerance for the sign test and the report/shapefile "
                         "coordinate agreement test; required with --check, no default")
    A = ap.parse_args()

    rep = parse_report(A.report_pdf)
    usgs = load_usgs(A.usgs_dir)
    contr = load_contractor(A.contractor_dir)

    print("=== inputs")
    print(f"  survey report        {A.report_pdf}")
    print(f"  USGS VA shapefiles   {A.usgs_dir}")
    print(f"  contractor shapefiles{A.contractor_dir}")
    counts = {}
    for pid, (_, _, _, code) in rep.items():
        counts[code] = counts.get(code, 0) + 1
    print(f"\n=== report coordinate tables: {len(rep)} points {counts}")
    print(f"    report section 1.3 text says            {REPORT_COUNTS}")
    print(f"    USGS VATool shapefiles: {len(usgs)} unique checkpoints, "
          f"{sum(len(v) for v in usgs.values())} rows "
          f"({sum(1 for v in usgs.values() if len(v) > 1)} tested against both blocks)")
    only_usgs = sorted(set(usgs) - set(rep))
    only_rep = sorted(p for p in rep if rep[p][3] != "LCP" and p not in usgs)
    print(f"    in the VATool files but NOT in the report tables: {only_usgs}")
    print(f"    NVA/VVA in the report tables but NOT tested:      {only_rep}")

    if A.check:
        import numpy as np
        if A.tol_m is None:
            raise SystemExit("--check requires --tol-m; there is no default tolerance")
        tol = A.tol_m

        print(f"\n=== sign test at tolerance {tol:g} m: is zdiff == Z - surface?")
        for block in ("QL1", "QL0"):
            rows = [v[block] for v in usgs.values() if block in v]
            for surf, zc, dc in (("DEM", "DEMz", "zdiff"), ("LAZ", "LAZz", "LAZzdiff")):
                Z = np.array([r["Z"] for r in rows])
                S = np.array([r[zc] for r in rows])
                D = np.array([r[dc] for r in rows])
                a = np.abs(D - (Z - S))
                b = np.abs(D - (S - Z))
                print(f"  {block} {surf}:  surveyed - surface  "
                      f"{int((a < tol).sum())}/{len(rows)} rows, max |resid| {a.max():.3e} m"
                      f"   |   surface - surveyed  {int((b < tol).sum())}/{len(rows)}, "
                      f"max {b.max():.3e} m")

        print(f"\n=== report table vs VATool shapefile, same mark "
              f"(tolerance {tol:g} m / same code)")
        de = dn = dz = 0.0
        nmatch = mism = 0
        for pid, blocks in usgs.items():
            if pid not in rep:
                continue
            E, N, Z, code = rep[pid]
            r = next(iter(blocks.values()))
            de = max(de, abs(r["X"] - E)); dn = max(dn, abs(r["Y"] - N))
            dz = max(dz, abs(r["Z"] - Z))
            nmatch += 1
            mism += (code != r["Lndcover"])
        print(f"  {nmatch} marks in both: max |dE| {de:.4f} m, max |dN| {dn:.4f} m, "
              f"max |dZ| {dz:.4f} m, land-cover code mismatches {mism}")
        print(f"  agree within {tol:g} m in all three coordinates: "
              f"{'YES' if max(de, dn, dz) < tol else 'NO'}")

        print("\n=== published-aggregate reproduction (VATool's own report files)")
        print(f"{'block':6s} {'surf':4s} {'n_NVA':>6s} {'RMSEz cm':>9s} {'pub':>6s} "
              f"{'n_VVA':>6s} {'95pct cm':>9s} {'pub':>6s}")
        for block in ("QL1", "QL0"):
            rows = [v[block] for v in usgs.values() if block in v]
            lc = np.array([r["Lndcover"] for r in rows])
            for surf, dc in (("DEM", "zdiff"), ("LAZ", "LAZzdiff")):
                d = np.array([r[dc] for r in rows])
                nvn, nvr, vvn, vvp = PUBLISHED[(block, surf)]
                rmse = 100 * np.sqrt((d[lc == "NVA"] ** 2).mean())
                p95 = 100 * np.percentile(np.abs(d[lc == "VVA"]), 95, method="linear")
                print(f"{block:6s} {surf:4s} {int((lc=='NVA').sum()):6d} {rmse:9.4f} "
                      f"{nvr:6.2f} {int((lc=='VVA').sum()):6d} {p95:9.4f} {vvp:6.2f}")
        print("  (the VVA figure is the 95th percentile of |residual|; numpy's default")
        print("   'linear' interpolation is the one that reproduces the published value.")
        print("   Alternatives, QL1 LAZ: ", end="")
        rows = [v["QL1"] for v in usgs.values() if "QL1" in v]
        lc = np.array([r["Lndcover"] for r in rows])
        d = np.abs(np.array([r["LAZzdiff"] for r in rows]))[lc == "VVA"]
        print(", ".join(f"{m}={100*np.percentile(d,95,method=m):.2f}"
                        for m in ("lower", "higher", "nearest", "midpoint")))

        print("\n=== geoid, asserted per mark")
        bad = [pid for pid, bl in contr.items()
               for r in bl.values() if r["geoid"] != "Geoid 18"]
        if bad:
            raise SystemExit(f"geoid attribute is not 'Geoid 18' on {len(bad)} marks: "
                             f"{bad[:5]}")
        print(f"  contractor shapefile 'geoid' == 'Geoid 18' on all "
              f"{sum(len(v) for v in contr.values())} NVA/VVA rows "
              f"({len(contr)} unique marks)")
        no_shp = sorted(p for p in rep if p not in contr)
        lcp = [p for p in no_shp if rep[p][3] == "LCP"]
        other = [p for p in no_shp if rep[p][3] != "LCP"]
        print(f"  {len(no_shp)} marks in the report tables appear in no shapefile: "
              f"the {len(lcp)} LCPs, plus {other}; for those the geoid comes from the "
              f"report's own per-table header, 'Geoid Model: Geoid18' (section 2.2), "
              f"and from section 1.8.4")

    fields = ["point_id", "point_type", "role", "easting", "northing", "elevation",
              "elevation_units", "horizontal_crs", "vertical_datum", "geoid_model",
              "project_id", "collected", "source", "verified",
              "va_blocks",
              "usgs_ql1_dem_z_m", "usgs_ql1_dem_error_m",
              "usgs_ql1_laz_z_m", "usgs_ql1_laz_error_m",
              "usgs_ql0_dem_z_m", "usgs_ql0_dem_error_m",
              "usgs_ql0_laz_z_m", "usgs_ql0_laz_error_m"]

    # Union of the two coordinate sources, report first.
    merged = {}
    for pid, (E, N, Z, code) in rep.items():
        merged[pid] = dict(easting=E, northing=N, elevation=Z, code=code,
                           src="Ground_Control_Survey_Report.pdf sections 2.2/2.3/2.4")
    for pid, blocks in usgs.items():
        if pid in merged:
            continue
        r = next(iter(blocks.values()))
        merged[pid] = dict(easting=r["X"], northing=r["Y"], elevation=r["Z"],
                           code=r["Lndcover"],
                           src="USGS_MN_SE_Driftless_2021_B21_QL*.dbf fields X/Y/Z "
                               "(this mark is MISSING from the survey report's tables)")

    os.makedirs(os.path.dirname(A.out) or ".", exist_ok=True)
    n_by_code = {}
    with open(A.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fields)
        w.writeheader()
        for pid in sorted(merged):
            m = merged[pid]
            blocks = usgs.get(pid, {})
            cb = contr.get(pid, {})
            collected = ""
            for r in cb.values():
                collected = _date(r["collection"])
                break
            row = dict(
                point_id=pid, point_type=m["code"], role=ROLE[m["code"]],
                easting=m["easting"], northing=m["northing"], elevation=m["elevation"],
                elevation_units="m", horizontal_crs="EPSG:6344",
                vertical_datum="NAVD88", geoid_model="GEOID18",
                project_id="MN_SE_Driftless_2021_B21", collected=collected,
                source=m["src"],
                verified=("geoid asserted PER MARK from the contractor shapefile "
                          "'geoid' attribute = 'Geoid 18'" if cb else
                          "geoid from the survey report's per-table header 'Geoid "
                          "Model: Geoid18' (sections 1.8.4, 2.2); this mark is in no "
                          "shapefile"),
                va_blocks="+".join(sorted(blocks)),
            )
            for b in ("QL1", "QL0"):
                r = blocks.get(b)
                p = f"usgs_{b.lower()}_"
                row[p + "dem_z_m"] = "" if r is None else r["DEMz"]
                row[p + "dem_error_m"] = "" if r is None else r["zdiff"]
                row[p + "laz_z_m"] = "" if r is None else r["LAZz"]
                row[p + "laz_error_m"] = "" if r is None else r["LAZzdiff"]
            w.writerow(row)
            n_by_code[m["code"]] = n_by_code.get(m["code"], 0) + 1

    print(f"\nwrote {sum(n_by_code.values())} control points to {A.out}  {n_by_code}")
    print(f"  of which {sum(1 for p in merged if p in usgs)} carry per-point residuals; "
          f"{n_by_code.get('LCP', 0)} LCPs carry none and are labelled "
          f"role=calibration -- they calibrated gen2 and must be excluded from any "
          f"accuracy statement about it")


if __name__ == "__main__":
    main()
