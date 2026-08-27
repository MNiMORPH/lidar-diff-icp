# Bundled surveyed control

## `mn_se_driftless_2021_ql1_near_elba.csv`

The six USGS 3DEP vertical-accuracy checkpoints within ~12 km of Elba, MN, transcribed
from the **contractor-provided** QA shapefile of the 2021 SE-Minnesota Driftless QL1
block:

```
https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/
    MN_SE_Driftless_2021_B21/Vertical_Accuracy/contractor_provided/
    MN_Driftless_NVA_VVA_UTM15_QL1.shp
```

238 checkpoints in the full file; these six are the ones near the pilot site. The
transcription, and the checks behind it, are recorded in
`analysis/ridgelines/ABSOLUTE_ELEVATION_REFS.md` §1a.

**Units.** The shapefile `.dbf` labels `source_ele` as *US Feet*. That label is wrong --
the values are **metres**, verified two ways (a valley-floor "feet" reading would be
~106 m against a ~225 m NAVD88 valley floor; and the USGS EPQS DEM matches `source_ele`
to the centimetre at these coordinates). The CSV therefore records `elevation_units = m`
and `read_3dep_va_shapefile` has **no default** for the units: a caller must state them,
so the mislabel cannot pass through silently.

**Datum.** NAVD88 orthometric height on **GEOID18**, horizontal NAD83(2011) UTM 15N
(EPSG:6344, which agrees with EPSG:26915 to a few cm). gen1 is on **GEOID03**
(InPort 68818), so a tie must convert -- see `lidar_diff_icp.references.geoid_difference`.

This CSV is the offline fallback. Where the shapefile itself is available, read it with
`read_3dep_va_shapefile` and prefer it: it is the authoritative source and carries all
238 points.

## `mn_dnr_2008_control_semn.csv`

The **1 004** MnDNR ground control points used to validate the 2008 SE-Minnesota lidar --
**gen1's own control**, on gen1's own geoid -- parsed from the nine per-county validation
reports in the open MnGeo lidar tree:

```
https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/<county>/
    <County>_county_validation_report.pdf        (capitalisation varies by county)
```

Eight counties: Dodge, Fillmore, Houston, Mower, Olmsted, Steele, Wabasha, Winona -- the
ones for which the dataset metadata publishes a per-county RMSE. A ninth report exists
for Freeborn; Freeborn is **not** in that list and is not part of this acquisition, so it
is excluded. Regenerate with
`analysis/groundtruth/parse_mndnr_2008_control.py --pdf-dir <dir> --check`; `--check`
reproduces each report's own printed RMSE from the parsed rows.

**Sign.** The reports' `Error` column is **`Control Z - Surface Z`**, established by
arithmetic on 1 022 of 1 022 rows (the other order misses by up to 1.08 m), not from the
column name. Negative therefore means the delivered 2008 surface reads **above** the
mark. That is the same sign family as `groundtruth.tie`'s `tie = surveyed - z_lidar`.
`Surface Z` is the **delivered 2008 DNR DEM**, a different surface from our
reconstruction; it is carried as `dnr_surface_z_m` / `dnr_error_m` for comparison and is
never the answer.

**Datum.** NAVD88 on **GEOID03**, UTM 15N metres. The validation reports state *no* datum
and *no* geoid; the linkage is a **dataset-level** assertion in
`lidar_semn2008.html` ("Vertical datum: NAVD88 (Geoid03)"), recorded in each row's
`verified` field. The names carry `RTK`/`VRS`, so these are GPS-derived, not levelled, and
inherit GEOID03 in full. The horizontal realization is not stated anywhere; `EPSG:26915`
is recorded as the best available reading and the NAD83 realization difference against
the 2011-epoch 2021 marks is **not** modelled.

**Land cover** is encoded in the point name and carried in `point_type`: `L1O` open
terrain, `L2T` tall weeds and crops, `L3B` brush and low trees, `L4F` forested, `L5U`
urban. It matters: pooled over all classes the residual is dominated by vegetation in the
ground class. See `analysis/GEN1_OWN_CONTROL_TIE.md`.

## `mn_se_driftless_2021_control.csv`

The **534** surveyed control points of the 2021 3DEP MN_SE_Driftless project -- **gen2's
own control** -- with, for 390 of them, the **per-point residual USGS publishes**. This
supersedes the six-mark file above for everything except the offline-fallback role that
file's docstring describes.

Regenerate with `analysis/groundtruth/parse_gen2_control.py --check --tol-m <tol>`; the
script's docstring names all three sources and `--check` reproduces every published
aggregate from the parsed rows.

**Where the residuals come from.** Not the survey report and not the contractor's
checkpoint shapefile -- neither carries a lidar elevation. They are in the *USGS* NGTOC
"VATool" output shapefiles,

```
https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/
    MN_SE_Driftless_2021_B21/Vertical_Accuracy/USGS/
    USGS_MN_SE_Driftless_2021_B21_QL{0,1}.dbf
```

whose fields are `srcChkptId, X, Y, Z, Lndcover, DEMz, zdiff, zdiffSq, LAZz, LAZzdiff,
LAZzdiffSq`. `Z` is the surveyed height; **`DEMz` is the delivered OPR DEM** read at the
mark and **`LAZz` the delivered classified point cloud** read at the mark. Two delivered
surfaces, where gen1's validation reports give one.

**Sign.** `zdiff == Z - DEMz` and `LAZzdiff == Z - LAZz`, exact on all 395 rows (max
residual 5e-16 m; the other order misses by up to 1.098 m). Negative therefore means the
delivered 2021 surface reads **above** the mark -- the same sign family as
`groundtruth.tie`'s `tie = surveyed - z_lidar`, and the same as gen1's `dnr_error_m`.

**`role` is the column that matters.** The vendor FGDC metadata states that the 143
**LCP**s were used to *calibrate* the lidar and that the NVA/VVA checkpoints "were not
used to calibrate or post process the data". LCP rows therefore carry `role=calibration`
and **no residual** (the VATool never tested them); using them to check gen2 would be
circular. The 227 NVA and 164 VVA rows carry `role=check`.

**Counts, and two discrepancies carried rather than hidden.** The survey report's §1.3
text says 143 + 227 + 164 = 534. Its coordinate tables hold **533**: the 164th VVA,
`3000_2021_MN`, is missing from them and is recovered from the USGS shapefile, whose X/Y/Z
agree with the report to 0.0000 m on all 389 marks the two share. One NVA id carries a
letter suffix (`2198A_2022_MN`) and a regex without `[A-Z]?` silently drops it. One
report-table VVA, `3021_2021_MN`, was never tested by the VATool and so has no residual.

**Datum.** NAVD88 on **GEOID18**, NAD83(2011) epoch 2010.00, UTM 15N, metres. gen2 is
delivered on the same geoid, so **a gen2 tie needs no geoid conversion** -- the reason
these marks are a direct absolute check where gen1's need `references.geoid_difference`.
The geoid is asserted **per mark** for the 390 shapefile marks (`geoid` attribute =
"Geoid 18" on all 395 rows) and from the report's own per-table header
("Geoid Model: Geoid18", §1.8.4 and §2.2) for the 144 that are in no shapefile.

**`va_blocks`** records which QL block(s) tested the mark. Five marks were tested against
both, and the two blocks' DEMs differ at them by up to 3 cm, so both are carried in
separate columns rather than averaged.

**What the residuals do not remove.** gen2 carries a vendor vertical **bias adjustment**
of unpublished magnitude (`MN_SE_Driftless_2_2021_Lidar_Mapping_Report.pdf` p. 15), tuned
against the LCPs. These checkpoints measure gen2 *after* that adjustment.
