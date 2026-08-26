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
