# Absolute-elevation reference points near Elba, MN

**Purpose.** Find surveyed points with known NAVD88 heights near Elba, MN
(Whitewater River valley, Winona County) that could serve as an *absolute*
vertical test of the two airborne-lidar epochs used in this project
(gen1 = 2008 MN DNR SE-Minnesota lidar; gen2 = 2021 USGS 3DEP SE-Driftless).

**Area of interest (AOI):** lat 44.09–44.13 N, lon 92.00–92.05 W
(UTM 15N / EPSG:26915: E 575 600–580 050, N 4 882 200–4 886 250).

**Verification status of this document.** Every coordinate, elevation, and datum
below was pulled from an authoritative source (USGS 3DEP metadata tree, NGS
datasheets, USGS EPQS) and is tagged with how it was checked. Aggregate lidar
accuracy numbers are quoted from the official project reports. Dead ends are
stated plainly. Nothing here is fabricated; where the raw data is not public,
that is said outright.

Compiled 2026-08-22.

---

## Bottom line (read this first)

1. **gen2 (2021 3DEP) surveyed QA checkpoints ARE downloadable** as a point
   shapefile with UTM 15N coordinates and surveyed NAVD88(GEOID18) elevations.
   **Zero checkpoints fall inside the AOI box, but 6 lie within ~12 km of Elba.**
   The published product gives checkpoint *locations + surveyed elevations* and
   *project-level* RMSE, but **not** the per-point lidar−checkpoint residual.
   These are still the cleanest absolute tie available and can be tested directly
   against the gen2 point cloud/DEM. **This is the Priority-1 win.**

2. **gen1 (2008 SE-MN) surveyed QA checkpoints are NOT publicly downloadable as
   coordinates.** Only aggregate accuracy statistics are published (per-county
   RMSE, CVA/FVA). The full "Lidar Accuracy Assessment Report" exists as a
   deliverable but is not posted online in the USGS/MnGeo metadata I could reach;
   the USGS `MN_SEMN_2008` metadata directory returns 404. To get raw 2008
   checkpoints you would have to request the report from MnGeo / MN DNR.

3. **NGS benchmarks:** 6 marks fall inside the AOI box, clustered along TH 74 /
   the Whitewater valley around Elba. **Only one (DG8385 "9 DRL") is set in a
   flat concrete slab** that discrete-return lidar could plausibly measure, and
   even that is a small flagpole-base slab. It is **leveled** (geoid-independent),
   which is the preferred kind. The other five are on bridge abutments, rod caps,
   or a concrete monument top — narrow/vertical settings lidar cannot resolve —
   and two of those are additionally GPS-derived. **The box is poor in
   lidar-usable NGS marks.**

4. **Airport/runway:** **None in or near the box.** Nearest is Winona Municipal
   (KONA), ~30 km E, well outside the AOI. No usable runway flat test here.

---

## 1. Lidar QA/QC checkpoints (Priority 1)

### 1a. gen2 — 2021 USGS 3DEP "MN_SE_Driftless_2021_B21" (DOWNLOADABLE)

**Datum / geoid (verified from the shapefile `.prj` and project report):**
NAVD88 height, **GEOID18**, meters. Horizontal NAD83(2011) UTM 15N (EPSG:6344;
= EPSG:26915 to within a few cm). Vertical EPSG:5703.

**Project-level accuracy (verified — USGS project report + VA text files):**
For the QL1 block covering Winona County (which contains Elba):
- NVA (non-vegetated), lidar point cloud: **RMSEz = 3.54 cm**, 95%-conf = 6.94 cm, n = 139 checkpoints.
- NVA, bare-earth DEM: RMSEz = 3.51 cm, 95%-conf = 6.88 cm.
- VVA (vegetated), point cloud: **27.14 cm at 95th percentile**, n = 99.
- Sensor: Optech Galaxy T2000. Collection 2021-04 to 2022-05. Contractor: Woolpert.

**Where the raw checkpoints live (verified — files downloaded):**
`https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/MN_SE_Driftless_2021_B21/Vertical_Accuracy/`
- `contractor_provided/MN_Driftless_NVA_VVA_UTM15_QL1.shp` (+ `.dbf/.shx/.prj`, also `.gpkg`) — **238 checkpoint points**, UTM 15N, with surveyed elevation in field `source_ele`.
- `USGS/USGS_..._QL1_las_checkpoint_report.txt` and `..._QL1_VA.txt` — aggregate stats only (RMSE, 95th pct); **no per-point residuals**.

**Important unit note (verified):** the shapefile `.dbf` labels `source_ele`
units as "US Feet", but this is **wrong / an inherited field label**. The values
are **meters (NAVD88)**. Verified two ways: (1) physical plausibility — a valley-
floor "feet" reading would be ~106 m, far below the Whitewater floor (~225 m
NAVD88), impossible; (2) cross-check against the USGS 3DEP DEM via EPQS at the
checkpoint coordinates — e.g. checkpoint 2210 `source_ele` = 349.288, EPQS DEM =
349.288 m; checkpoint 3089 `source_ele` = 205.737, EPQS DEM = 205.69 m. The raw
values are meters and tie the 3DEP surface to the cm.

**Checkpoints in / near the AOI (verified — computed from the shapefile):**
In-box count = **0**. Within ~12 km of Elba = **6** (all `prj_id` 222538,
collected 2021-09-11, geoid GEOID18). Coordinates below are the *surveyed*
checkpoint position (UTM 15N) and surveyed NAVD88 elevation; lat/lon are my
conversion (EPSG:6344→WGS84):

| unique_ind | type | UTM15 E | UTM15 N | lat | lon | elev NAVD88 (m) | dist to Elba |
|---|---|---|---|---|---|---|---|
| 2210_2021_MN | NVA | 570492.1 | 4884126.1 | 44.10693 | −92.11918 | 349.288 | 7.2 km |
| 3056_2021_MN | VVA | 570473.8 | 4884246.4 | 44.10801 | −92.11939 | 353.259 | 7.2 km |
| 2024_2021_MN | NVA | 571243.8 | 4887693.8 | 44.13897 | −92.10931 | 344.735 | 7.4 km |
| 2036_2021_MN | NVA | 585982.1 | 4884249.8 | 44.10638 | −91.92564 | 353.119 | 8.3 km |
| 2099_2021_MN | NVA | 571954.9 | 4893677.0 | 44.19276 | −92.09960 | 355.053 | 11.3 km |
| 3089_2021_MN | VVA | 582881.0 | 4894754.2 | 44.20130 | −91.96273 | 205.737 | 11.9 km |

- **NVA points** (2210, 2024, 2036, 2099) sit on open/non-vegetated ground — the
  best absolute flat targets. All are ridge-top/upland (~345–355 m), i.e. on the
  Driftless uplands surrounding the valley, not on the valley floor. They test
  the *gen2* surface directly at a surveyed height.
- **VVA points** (3056, 3089) are under vegetation — less clean for an absolute
  ground test; expect the 25–27 cm VVA-scale spread.

**Usable as an absolute test?** **Yes, for gen2.** These are surveyed
NAVD88(GEOID18) heights on the same datum as the gen2 lidar, and their locations
are known to sub-meter. Extract gen2 ground points within a small radius of each
NVA checkpoint and compare to `source_ele`. Caveat: none are *inside* the AOI box,
so this is a regional (7–12 km) absolute check of the gen2 collection, not a
point test at the study site. It does **not** by itself test gen1.

### 1b. gen1 — 2008 MN DNR "Lidar Elevation, Southeast Minnesota, 2008" (checkpoints NOT public)

**Datum / geoid (verified — InPort record 68818):** NAVD88, **GEOID03**,
NAD83(NSRS2007). Collected Nov 2008 by AeroMetric (now Quantum Spatial); DNR QA
from Apr 2009.

**Note:** gen1 uses **GEOID03**, gen2 uses **GEOID18**. Any absolute comparison
of the two epochs' checkpoints must account for the geoid-model difference (a
small, spatially smooth offset), or work in ellipsoidal heights.

**Aggregate accuracy (verified — InPort record 68818):**
- AeroMetric: FVA = 0.161 m @95% (open, 26 pts); CVA = 0.36 m @95% (127 pts).
- MN DNR independent: CVA = 0.287 m @95% (1,009 ground-control points, 5 land classes).
- **Per-county vertical RMSE — Winona = 0.161 m** (the county containing Elba).
  (Others: Dodge 0.129, Fillmore 0.155, Houston 0.110, Mower 0.161, Olmsted
  0.117, Steele 0.125, Wabasha 0.106 m.)

**Raw checkpoint coordinates:** **Not found online.** The InPort metadata lists
"One paper copy of the Lidar Accuracy Assessment Report" as a deliverable but
gives no download link and no checkpoint file. The USGS metadata directory
`.../metadata/MN_SEMN_2008/` returns **HTTP 404**, so there is no USGS-side
`Vertical_Accuracy` shapefile analogous to gen2's. **Dead end for public raw
coordinates.**

**How to get them:** request the "Lidar Accuracy Assessment Report" (and any
checkpoint shapefile/CSV) from **MnGeo / MN DNR** — the 1,009-point DNR control
set and the AeroMetric 127-point set both had surveyed coordinates that are not
posted publicly. This is the missing piece for a *gen1* absolute test.

**Usable as an absolute test?** Not yet — no coordinates in hand. If the report's
checkpoint list is obtained, several of the 1,009 DNR points likely fall near
Elba and would give a genuine gen1 absolute tie.

---

## 2. NGS control marks in the AOI box (Priority 2 — filtered to flat, stable surfaces)

All six marks below fall **inside** the target box (verified — WGS84→UTM15N).
Heights, settings, and adjustment method are **verified from individual NGS
datasheets** (ngs.noaa.gov `ds_mark.prl`). All ortho heights are on GEOID18.

| PID | Desig. | lat | lon | NAVD88 (m) | height source | setting | lidar-usable? |
|---|---|---|---|---|---|---|---|
| DG8385 | 9 DRL | 44.1202 | −92.0035 | 223.352 | **LEVELED** (adj. 2005) | **concrete mat/slab** (flagpole base) | **maybe — see below** |
| DO2105 | 8511 A 4 | 44.0916 | −92.0131 | 226.279 | LEVELED (adj. 2013) | bridge abutment/pier | no (bridge; classified out) |
| DG8388 | TT 20 E | 44.0940 | −92.0165 | 226.871 | LEVELED (adj. 2005) | top of concrete monument | no (narrow monument top) |
| DG8387 | 8508 G | 44.1002 | −92.0092 | 228.973 | LEVELED | disk on 10-ft+ aluminum rod | no (rod cap, not a surface) |
| DR9927 | 8508 X | 44.0951 | −92.0158 | 229.99 | **GPS-derived** | flange-encased rod | no (rod + GPS) |
| DR9925 | 8508 Y | 44.1217 | −92.0047 | 218.75 | **GPS-derived** | aluminum rod (10 ft+) | no (rod + GPS) |

**Leveled vs GPS-derived:** four marks (DG8385, DO2105, DG8388, DG8387) have
**leveled** orthometric heights — geoid-model-independent, the preferred kind.
Two (DR9927, DR9925) are **GPS-derived** and carry geoid ambiguity — note this if
ever used.

**The one candidate flat surface — DG8385 "9 DRL":**
- Set in a **concrete mat foundation / slab** (NGS setting code, not pavement),
  0.4 ft south of a **flagpole**, ~150 ft W of TH 74, near a **DNR headquarters
  building** — i.e. the Whitewater Wildlife Management Area / DNR area HQ,
  **2.3 mi (3.7 km) NE of Elba**. Established USGS 1972; recovered "GOOD" 2019-10-08
  and 2023-04-05.
- **Leveled** NAVD88 = 223.352 m, stability "may hold, subject to surface motion."
- **Usable as an absolute flat test? Marginally.** It is the only in-box mark on
  a flat horizontal slab, and it is leveled — good. But a flagpole-base slab is
  small (a few m²), so a discrete-return pass may put few clean ground points
  squarely on it, and it sits beside a flagpole and building (edge/multipath
  risk). Worth checking whether both epochs actually have ground returns on the
  slab before trusting it. Its "subject to surface motion" stability also weakens
  it as a decade-spanning datum for the *DoD*.

**Honest summary of §2:** the box is **poor in lidar-usable NGS marks**. Five of
six are on rods, a bridge, or a monument top — settings discrete-return lidar
cannot measure. Only DG8385 is on a flat slab, and it is a small one. The NGS
network here is a weak absolute-flat resource; the gen2 3DEP NVA checkpoints
(§1a) are far better.

---

## 3. Airport / runway note

**No airport, airstrip, or runway in or near the AOI box.** Searches for a
public or private strip near Elba / Altura / the Whitewater valley returned
nothing in the box. The nearest is **Winona Municipal Airport (KONA / Max Conrad
Field)**, ~3 mi NW of Winona and roughly **30 km east** of Elba — well outside the
AOI and not useful as a local flat test. (FAA airport-elevation / runway-end data
exist for KONA if a distant large-flat check is ever wanted, but it is out of
area.) **No runway flat test is available here.**

---

## Appendix: sources and verification trail

**gen2 3DEP (Priority 1 win):**
- InPort 70275 — 2021-2022 USGS Lidar (QL1) SE-MN Driftless (aggregate accuracy, datum).
  https://www.fisheries.noaa.gov/inport/item/70275
- Project report PDF (parsed locally with pdftotext; WebFetch could not read it):
  https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/MN_SE_Driftless_2021_B21/USGS_MN_SE_Driftless_2021_B21_Project_Report.pdf
- Checkpoint shapefile (downloaded, read with pyshp; 238 pts):
  `.../MN_SE_Driftless_2021_B21/Vertical_Accuracy/contractor_provided/MN_Driftless_NVA_VVA_UTM15_QL1.shp`
- Elevation-unit cross-check: USGS EPQS (`epqs.nationalmap.gov`) DEM at checkpoint
  coords matched `source_ele` to the cm → values are meters, not feet.

**gen1 2008 (checkpoints not public):**
- InPort 68818 — 2008 MN DNR Lidar SE Minnesota (aggregate accuracy, per-county
  RMSE, GEOID03). https://www.fisheries.noaa.gov/inport/item/68818
- MnGeo first-gen lidar page: https://mn.gov/mngeo/gis-data-and-maps/info-by-topic/elevation/lidar/lidar-2008-2012.jsp
- `.../metadata/MN_SEMN_2008/` → HTTP 404 (no public checkpoint tree).

**NGS marks (verified per-datasheet):**
- ngs.noaa.gov `ds_mark.prl?PidBox=<PID>` for DG8385, DO2105, DG8388, DG8387, DR9927, DR9925.

**Airport:** web search (no strip in box); Winona Municipal (KONA) is the nearest, out of area.

**Local files retained (scratchpad, session-scoped — not committed):**
`vaqc/MN_Driftless_NVA_VVA_UTM15_QL1.*`, `vaqc/USGS_QL1_VA.txt`,
`vaqc/USGS_QL1_checkpoint_report.txt`, `report_3dep.txt`.
