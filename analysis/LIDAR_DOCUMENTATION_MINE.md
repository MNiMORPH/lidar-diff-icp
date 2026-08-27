# Mining the lidar documentation sets: 2008 gen1 and 2021 gen2

**Compiled 2026-08-26.** Companion to `analysis/ADDITIONAL_GROUND_CONTROL.md`, which
mined the *control* documents. This one mines the *acquisition and processing*
documents. Nothing here is committed to git; downloads went to the session
scratchpad.

**Why this exists.** Twice the answer sat in a report rather than in the
machine-readable deliverable. So the method here is to walk the document trees –
directory by directory – rather than the download endpoints. That method paid
again: the 2021 project turns out to publish a per-work-unit **Lidar Mapping
Report**, **sensor calibration certificates**, **flight logs**, and a
**swath-polygon geodatabase whose `PS_ID` is exactly our `point_source_id`**; and
the 2008 project turns out to publish a **statewide flight-line layer** that
carries every SE Minnesota line, dated, on a plan laid out at a constant
43.31 arcsecond longitude interval (962.8 m at Elba's latitude).

Every number below is quoted verbatim from a named document with its page, or is
marked as measured by us, or is marked **unverified**. Nothing is inferred from
plausibility.

---

## Bottom line

1. **The 2008 vertical bias adjustment is not published, and the negative is now
   well covered.** The sentence asserting it exists appears in AeroMetric's own
   FGDC record (found under USGS `legacy/`, a path I had first written off as a
   404), in MnGeo's derivative of that record, and in NOAA's copy of it. **The
   value appears in none of them, nor in the 2012 FEMA certification letter, the
   Digital Elevation Committee minutes, the trade-press write-up, or the archived
   2009 project page.** It survives only in the paper accuracy assessment report
   and must be requested. §1.
2. **Both epochs had a vendor vertical bias adjustment applied, and neither value
   is published.** This is new: the 2021 Lidar Mapping Report says so in as many
   words (p. 15). Our DoD therefore differences two independently
   bias-adjusted surfaces, and *neither* offset is recoverable from the data.
   §1.2.
3. **Cross lines: the 2008 documentation records none.** The published
   `flight_lines_and_dates` layer holds 204 SE Minnesota lines and **every one is
   north–south** – zero lines within 45° of east–west. Yet `point_source_id`
   10010, heading 270.98°, is in the delivered point cloud, on the same sortie
   (`SWATH_DEGENERACY_BREAKING.md` §2). **At least one cross line was flown, and
   the documentation does not record it.** §2.
4. **The 2021 lines over Elba run east–west** – bearing 90.0–90.8°, i.e.
   **perpendicular to gen1's north–south lines**. Verified from the vendor swath
   polygons, whose `PS_ID` 3040/3041/3042 are the `point_source_id` values in our
   own gen2 file. §2.4.
5. **Swath adjustment.** gen1: one set of TerraMatch orientation corrections
   (heading, pitch, roll, scale) "applied to the entire dataset" – *no* per-line
   or per-swath adjustment is described, and no relative-accuracy statistic is
   published. gen2: automated line-to-line calibration plus published interswath
   RMSDz, **one of whose 15 test polygons sits 1.8 km from the Elba reference
   point at RMSDz 0.009 m**. §3.
6. **The published flight plan pins gen1's line spacing at 962.8 m at Elba** (sd
   0.16 m over 145 gaps), and shows the plan was laid out at a constant longitude
   interval, so the metric spacing varies with latitude – 967.8 m at N 4 850 000,
   960.5 m at N 4 900 000. Fed back through our own `S/h` it brings our two
   independent flying-height estimates into agreement. §2.2.
7. **A documentation error worth knowing:** the 2008 metadata says "Data was
   acquired between November 18 - 24, 2008", and the FEMA certification letter
   repeats it. **Elba's own lines were flown 2008-11-25** – outside that window –
   on the flight-line layer's own date field, and our gen1 GPS week-seconds
   confirm a Tuesday. §2.2.
8. **The published 2008 control counts do not agree with each other.** The FEMA
   certification letter is the only source that gives **Freeborn County (0.144 m,
   126 points)**, and it transposes Houston relative to MnGeo (0.134 m / 110 vs
   0.110 m / 134). Neither list sums to the "1009 control points" both assert.
   §1.4b.

---

## 1. The 2008 vertical bias adjustment (Priority 1)

### 1.1 What the documentation says, verbatim

`https://resources.gisdata.mn.gov/pub/data/elevation/lidar/documentation/lidar_semn2008.xml`
(17,360 bytes; identical wording in the `.html`), Lineage, Lidar Processing,
process step 8:

> "AeroMetric provided Quality Assurance and Quality Control (QA/QC) data for this
> project. AeroMetric captured 127 QA/QC points in multiple land cover categories
> that were used to test the accuracy of the lidar ground surface. TerraScan's
> Output Control Report (OCR) was used to compare the QA/QC data to the lidar
> data. This routine searches the lidar dataset by X and Y coordinate, finds the
> closest lidar point and compares the vertical (Z) values to the known data
> collected in the field. **Based on the QA/QC data, a bias adjustment was
> determined, and the results were applied to the lidar data.** A final OCR was
> performed with a resulting RMSE of 0.109 meters."

**That is the whole of it. The magnitude of the bias is not stated, and no other
public document I reached states it.** Note also what the sentence does *not*
say: it does not say the adjustment was a single constant, it does not say whether
it was applied per lift, per line or globally, and it does not give its sign
convention.

Two details in the same record bear on how to interpret it:

- The comparison was **point-to-point**, not surface-to-point: "searches the lidar
  dataset by X and Y coordinate, **finds the closest lidar point** and compares the
  vertical (Z) values". So the bias was determined against the *point cloud*, not
  against a DEM or TIN.
- The final residual after adjustment was **RMSE 0.109 m** on those 127 points.

### 1.2 The same thing happened to gen2, and it is stated more plainly

`MN_SE_Driftless_2_2021_Lidar_Mapping_Report.pdf` (Woolpert, November 2022;
1,123,464 bytes), **p. 15**, §2.3 Lidar Data Classification:

> "Statistical absolute accuracy was assessed by direct comparisons of ground
> classified points to ground RTK survey data. **Based on the statistical analysis,
> the lidar data was then adjusted to reduce the vertical bias when compared to the
> survey ground control of higher accuracy.**"

And the vendor FGDC metadata
(`MN_SEDriftless_2_2021_Classified_Point_Cloud_Metadata.xml`, 12,522 bytes,
Ground Conditions):

> "Woolpert established ground control points that were used to **calibrate the
> lidar to known ground locations** established throughout the entire project area.
> An additional independent accuracy checkpoints were collected throughout the
> entire project area and used to assess the vertical accuracy of the data. These
> checkpoints were not used to calibrate or post process the data."

**This is the load-bearing find of the whole search.** It identifies which of the
534 surveyed 2021 points did what: the **143 LiDAR Control Points calibrated the
data**; the **227 NVA and 164 VVA checkpoints were held out**. So a tie made at an
LCP (`ADDITIONAL_GROUND_CONTROL.md` §1.5, §1c – LCPs 1078, 1079, 1080 near Elba)
is **not independent of gen2's own vertical calibration**, whereas a tie made at
an NVA/VVA mark is. That distinction was not previously drawn and it changes how
those marks should be weighted.

**Neither bias value is published.** For gen1 that is the single number that would
resolve the sign tension in `ADDITIONAL_GROUND_CONTROL.md` §5. For gen2 it is the
matching unknown on the other side of the difference. Neither is recoverable by
any processing of ours: both are constants already folded into the delivered
elevations.

### 1.3 Where I looked and did not find it

- MnGeo `documentation/` tree – 39 files listed; the only SE MN items are
  `lidar_semn2008.html` / `.xml`. No vendor report.
- MnGeo `county/<county>/` – validation report PDF, tile list, tile index map,
  licence, READMEs, LAZ and geodatabase. No vendor report.
- MnGeo `q250k/q4342/` (the quad containing Elba) – tile index PDF, LAZ,
  geodatabase, buildings. No reports.
- MnGeo `projects/` – ten project folders, **none of them SE MN**.
- MnGeo `lidar/` root – `LAS_File_Processing_Using_LASTOOLS.pdf`,
  `about_this_data.rtf`, `raw_LiDAR_Data_README.rtf`, `county_mosaic_readme_first.rtf`,
  `readme_first.rtf`, `readme_first.txt`. I read all five text documents: they are
  Tim Loesch's 2011 *statewide* delivery-disk documentation ("the Minnesota State
  LiDAR collect circa 2011"), not SE MN 2008, and none mentions bias, calibration
  or flight lines.
- USGS – **and here I was wrong at first, so the correction is recorded rather
  than quietly fixed.** `.../Elevation/LPC/Projects/MN_SEMN_2008/` and
  `.../Elevation/metadata/MN_SEMN_2008/` do both return 404, and the project is
  in neither current directory index (14 MN projects, earliest
  `MN_BlueEarth_2011`), so I concluded USGS had dropped it. **It has not been
  dropped; it moved under `legacy/`**, and the correct path is
  `https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/legacy/MN_SEMN_2008/`
  (HTTP 200, verified by me). See §1.4. NOAA InPort 68818 still advertises the
  old non-legacy path; **that link is dead**.
- NOAA: InPort 68818's lineage is a verbatim copy of the MnGeo record – same
  sentence, no number. The NOAA S3 bucket `noaa-nos-coastal-lidar-pds` under
  `laz/geoid18/9684/` holds **only COPC LAZ tiles**; the bucket's top level has
  four prefixes (`dem/`, `entwine/`, `laz/`, `testing/`) and no supplemental or
  report prefix. `chs.coast.noaa.gov/htdata/lidar1_z/{geoid12b,geoid18}/data/9684/`
  is 404.

The metadata's own Final Deliverables list explains why: deliverable 1 is

> "One paper copy of the Lidar Accuracy Assessment Report."

**One paper copy.** That is very likely why nothing is online. See §6 for the
request.

### 1.4 The documents that do exist, and what each does not say

Five 2008-era documents turned up that were not previously in play. **None states
the bias.** They are recorded here so the search does not get repeated.

**(a) The USGS legacy metadata – the vendor's own FGDC record.**
`.../metadata/legacy/MN_SEMN_2008/` holds `semn_metadata.xml` (14,356 bytes) and
`semn_metadata.docx` (19,319 bytes), plus per-county tile-index and hydro-breakline
`.lyr` files and a `DobbinsCreek_tile_index` shapefile. The LAZ sits at
`.../LPC/Projects/legacy/MN_SEMN_2008/LAZ/`. I downloaded and read both metadata
files. This is **AeroMetric's own wording**, of which the MnGeo record is a lightly
edited derivative: same process steps 1–10, same "a bias adjustment was determined,
and the results were applied to the LIDAR data. A final OCR was performed with a
resulting RMSE of 0.109m." **Still no value.** It does add three things the MnGeo
copy drops:

> "The Consolidated Vertical Accuracy (CVA) of the TIN achieved a 0.205 meters at a
> 95% confidence level in all land cover categories. 127 control points were used
> in this evaluation."

– a second, tighter CVA on the same 127 points than the 0.36 m ASPRS figure quoted
elsewhere in the same record; the processing tiling scheme, "a ILDOT specified 2000
meter by 2000 meter tiling scheme"; and the horizontal datum written as
"NAD83(2007)" rather than MnGeo's "NAD83 (NSRS2007)".

**A pointed contrast in the same tree.** `.../metadata/legacy/MN_PINECO_2006/`
publishes **`Pine_LIDAR_Final_Report.doc`** – the vendor final report for the 2006
Pine County project – alongside its control shapefiles. Verified by listing the
directory. **There is no analogous file under `MN_SEMN_2008/`.** So USGS does host
vendor final reports for Minnesota legacy projects when it has them; SE MN's simply
was not deposited. That is the strongest precedent to cite when asking for it, and
it sits beside the `lidar_checkpts_pine2007` precedent already noted in
`ADDITIONAL_GROUND_CONTROL.md` §5.

**(b) The FEMA certification letter – and it contains a number nobody else
publishes.** `https://files.dnr.state.mn.us/waters/watermgmt_section/floodplain/mn_lidar_certification-se_mn_10-22-2012.pdf`,
81,131 bytes, 2 pages, 22 October 2012, from **Peter W. Jenkins, PLS, CFedS,
Photogrammetric Unit Supervisor** (MnDOT) to Eric Ratcliffe, STARR MT-1 Project
Manager, Atkins Global. Downloaded and read. Nothing on bias, cross lines, lifts or
boresight – but its per-county table is **not** the MnGeo table:

> "The vertical RMSE and sample count per county as tested by the State of
> Minnesota is as follows: Dodge a.129m, 121; Fillmore 0.155m, 128; **Freeborn
> a.144m, 126**; **Houston 0.134m, 110**; Mower a.161m, 115; Olmsted a.117m, 125;
> Steele 0.125m, 137; Wabasha a.106m, 97; Winona a.161m, 176."

(The `a.` are a font-encoding artefact of the text extraction; they are `0.`.)

Two things follow. **(i) This is the only public source for Freeborn County –
0.144 m over 126 points.** MnGeo's record lists eight counties and omits Freeborn
entirely. **(ii) Houston is transposed between the two sources**: MnGeo says
"Houston 0.110m, 134", the certification letter says "Houston 0.134m, 110". One of
them has swapped the RMSE and the count. **Which is right is not resolvable from
these two documents** – parsing `Houston_county_validation_report.pdf` would settle
it in one step, and that has not been done. Note that every other county agrees
exactly between the two sources.

**(iii) And neither count list reaches the stated total.** Both documents say
"1009 control points". MnGeo's eight counties sum to **1 033**; the certification
letter's nine sum to **1 135** (or 1 159 if Houston is 134). **The published
control counts are internally inconsistent in all three copies.** That is worth
raising in the same email as everything else.

The letter also repeats "Date of acquisition: November 18-24, 2008" – the same
window the flight-line dates contradict (§2.2).

**(c) The 2009 MnDNR project page, from the Wayback Machine.** A single capture,
2009-02-11, of `dnr.state.mn.us/mis/gis/semn_lidar/index.html`, states that an RFP
was advertised in July 2008, nine proposals were received, and AeroMetric was
awarded the contract in August 2008. **The RFP is the document that would say
whether cross lines were required**, and no copy of it was found anywhere.
*(Sub-agent finding; I have not re-fetched the capture.)*

**(d) MnGeo Digital Elevation Committee minutes, 26 February 2009**
(`https://www.mngeo.state.mn.us/committee/elevation/elev_09feb.pdf`, 13,514 bytes).
I downloaded and searched it. The whole of what it says is:

> "DNR LiDAR project in southeast Minnesota – Tim Loesch gave an update about
> ongoing activities including quality control / accuracy testing and data delivery
> issues."

No numbers. Adjacent 2008–2009 minutes were checked by the sub-agent with the same
result.

**(e) MN GIS/LIS News, Fall 2009, issue 58**, "Southeast Minnesota LiDAR Project
Completed!" by Tim Loesch – the per-county RMSE table and a narrative, no bias and
no flight plan. *(Sub-agent finding; the live site is behind Cloudflare and I have
not re-fetched the Wayback copy.)*

**Conclusion on Priority 1: a clean, well-covered negative.** The bias value is
absent from the vendor's own FGDC record, MnGeo's derivative of it, NOAA's copy of
that, the FEMA certification letter, the committee minutes, the trade-press
article, and the archived project page. It exists only in the paper accuracy
assessment report. **It must be requested.**

---

## 2. Flight-line geometry, and cross lines (Priority 2)

### 2.1 The document that exists

`https://resources.gisdata.mn.gov/pub/data/elevation/lidar/tile_index/flight_lines_and_dates.zip`
– **759,353 bytes**, last modified 2012-06-26. Shapefile `flight_lines.shp`
(2,312,424 bytes), EPSG:26915, **2,893 features** across seven Minnesota
first-generation projects, fields `OBJECTID, Date_Flown, Project, Shape_Leng`.
A `.kmz` twin (1,050,400 bytes) was not downloaded.

`Project = 'SE Minnesota'` selects **204 lines**. Median line length **77.3 km**
(range 1.6–108.4 km).

### 2.2 What it says about Elba

**Headings.** Folded to an axis, all 204 lines fall in two bins: 54 at 0–5° and
150 at 175–180°. **Zero lines lie within 45° of east–west.** The layer is a pure
north–south boustrophedon plan.

**Spacing, and it is not a metric constant.** At N 4 883 678 (Elba's northing),
146 lines cross; the 145 consecutive gaps between their crossing eastings are
**mean 962.84 m, median 962.76 m, range 962.20–963.35 m, sd 0.16 m**. That is a
plan, not a fit. But repeat it at other northings and the number moves
systematically:

| northing | gaps | mean spacing | sd |
|---|---|---|---|
| 4 900 000 | 56 | **960.49 m** | 0.15 m |
| **4 883 678 (Elba)** | 145 | **962.84 m** | 0.16 m |
| 4 850 000 | 199 | **967.81 m** | 0.22 m |

**Because the plan was laid out in longitude, not in metres.** Converting the
crossing eastings to geographic coordinates, the interval between adjacent lines is
**43.3063 arcseconds of longitude (0.0120295°), with a standard deviation of
0.003 arcsec, and it is the same number at both northings tested.** The metric
spacing therefore shrinks northward as a degree of longitude shortens – 7.3 m over
50 km of northing.

**Use 962.8 m at Elba**, not a single project-wide figure. (This is the source of
a discrepancy worth naming: a spacing of "967 m" is also correct, at N 4 850 000.
Neither is wrong; the question is where.)

**Line identification.** The layer has no line-number field, but the geometry
identifies itself. Taking our own fitted nadir tracks (`ELBAEXT2_SCOPE.md` §2,
easting where each track crosses N 4 884 126) against the shapefile's crossing
easting at the same northing:

| `point_source_id` | our fitted track E | shapefile track E | difference |
|---|---|---|---|
| 128 | 570 345 | 570 357.8 | +12.8 m |
| 129 | 571 311 | 571 320.7 | +9.7 m |
| 130 | 572 275 | 572 283.6 | +8.6 m |
| 131 | 573 286 | 573 246.5 | −39.5 m |
| 132 | 574 222 | 574 209.3 | −12.7 m |
| 134 | 576 129 | 576 135.0 | +6.0 m |
| 135 | 577 119 | 577 097.9 | −21.1 m |
| 136 | 578 054 | 578 060.8 | +6.8 m |
| 137 | 579 046 | 579 023.7 | −22.3 m |
| 138 | 579 989 | 579 986.6 | −2.4 m |
| 144 | 585 756 | 585 764.2 | +8.2 m |
| 145 | 586 764 | 586 727.1 | −36.9 m |

n = 12, mean **−6.9 m**, sd **18.9 m**, max **39.5 m** – against a 962.8 m line
spacing. The identification is unambiguous:

> **`point_source_id` = `OBJECTID` − 2644**, for the SE Minnesota block, verified
> over lines 128–138 and 144–145.

That mapping makes the whole layer usable: **every gen1 flight line now has a
date**, per `point_source_id`, for the entire nine-county project.

**Elba's dates.** OBJECTIDs 2777–2782 = `point_source_id` 133–138, **all
`Date_Flown` 2008-11-25**. Lines 129–147 are all that date; line 128 and westward
are 2008-11-26.

**Independently confirmed from the data.** `data/las_local/4342-29-64.las` carries
GPS week-seconds (global encoding `gps_time_type = 0`, so seconds from Sunday
00:00 GPS time). Lines 135–138 give day-of-week **2.75–2.86** – **Tuesday**,
18:00–20:38 GPS time, i.e. 17:59:45–20:37:45 UTC (GPS ran 15 s ahead of UTC in
November 2008) or roughly noon to 14:40 local. **2008-11-25 was a Tuesday.** Two
independent sources agree.

**A second, sharper confirmation, at the date boundary.** The layer says line 128
was flown 2008-11-26 while lines 129–147 were flown 2008-11-25 – so it predicts
the date changes *between* 128 and 129. Reading `data/before/4342-29-61.laz`
(10,549,493 returns, all three lines present) gives GPS day-of-week **2.964** for
line 129, **2.951** for line 130 and **3.068** for line 128: Tuesday 23:08 GPS for
129 and **Wednesday 01:38 GPS for 128**. The date flips exactly where the layer
says it does. The layer's `Date_Flown` is therefore a **UTC** date – line 128 was
flown in the small hours of 26 November UTC, still the evening of the 25th local –
and lines 128–130 are a later sortie of the same night than Elba's 135–138. Their
nadir tracks in this tile (E 570 342, 571 308, 572 272) also reproduce the
`ELBAEXT2_SCOPE.md` values to 3 m.

**And the metadata is wrong about the window.** `lidar_semn2008.xml` states:

> "Data was acquired between November 18 - 24, 2008."

The flight-line layer's own dates for SE Minnesota run **2008-11-18 through
2008-11-27** (44 lines on 11-26, 37 on 11-25, 23 on 11-27, and 9 with a null
date; those nine are `OBJECTID` 2702–2710, at E 503–511 km, some 70 km west of
Elba, so the gap does not touch us). Elba's lines are on 11-25, **outside the
stated window**. Treat the
metadata's date range as unreliable and the per-line `Date_Flown` as the source of
truth.

**Flying height and sidelap, reconciled.** The metadata gives:

> "All these data products were acquired at 2400 meters above mean terrain (AMT)
> and have a horizontal accuracy of 0.40 meters, with a nominal point spacing of
> 1.0 meters."

and, under Acquisition parameters: "2. Flight Height - 2400 meters above mean
terrain; 3. Swath Width - 32 degrees; 4. Sidelap - 60%; 5. Nominal Post Spacing -
1.0 meter."

- **The height checks out.** `SWATH_ACROSS_TRACK_TEST.md` derives |h| = 2562 m
  from the fitted per-line coefficients. Its second route, `S / mean|sum_tan|`,
  used a spacing of 916 m and gave 2482 m; **with the published 962.8 m the same
  arithmetic gives 2609 m**, so the two routes now agree with each other (2562 and
  2610 m) far better than before. Both sit 7–9% above the documented 2400 m, which
  is a height above *mean* terrain: Elba sits in a dissected valley, so the height
  above the ground actually ranged is expected to exceed the project mean. That
  last sentence is reasoning, not measurement. The rest is arithmetic on our own
  fits with a documented constant substituted.
- **The 32° swath checks out.** At 2400–2560 m and ±16°, a swath is 1376–1470 m;
  our measured 99.5th-percentile half-widths are 696–742 m, i.e. swaths of
  1392–1484 m.
- **The "Sidelap - 60%" does not check out.** With a 962.8 m spacing and a
  1376–1470 m swath, the sidelap is **30–35%** (`SWATH_ACROSS_TRACK_TEST.md`
  independently says "41% sidelap 650 m wide" using its 916 m spacing). A 60%
  sidelap would require a 2408 m swath, which needs a 53° field of view – flatly
  inconsistent with the 32° in the line above it. **Do not use the 60% figure.**

**Lifts, missions, and line numbering: not documented.** The layer gives one date
per line and nothing else. There is no lift identifier, no sortie grouping, no
altitude per line, no aircraft, no sensor assignment. The metadata names two
sensors – "Optech ALTM Gemini (AeroMetric) and Leica ALS50-2 (Surdex)" – and says
Surdex's data was "post-processed to a raw point cloud by Surdex and then
delivered to AeroMetric to be merged into one raw point cloud dataset", but
**nothing says which lines belong to which sensor.** That is a real gap for us:
two different scanners in one merged cloud is a candidate systematic, and the
delivered data gives no flag for it.

### 2.3 Cross lines

**Documented: none.** Zero east–west lines in the 204. No document I found
mentions a cross line, a tie line, or a calibration line for the 2008 project.
The metadata's calibration step describes only "sampling the data collected across
all flight lines".

**Flown: yes.** `point_source_id` **10010**, heading **270.98°**, 747,107 returns
in tile `4342-28-64`, gps_time 238,320–238,344 s of week – between line 138's pass
(237,404 s) and line 137's (243,608 s), so **the same sortie on 2008-11-25**. It
carries **zero class-2 returns**; all 745,942 classified points are class 12,
overlap. Full analysis in `analysis/SWATH_DEGENERACY_BREAKING.md` §2, which is
where this was found; I re-read its gps_time here and confirm the day-of-week is
2.758–2.759, the same Tuesday.

**So the state of knowledge is: at least one cross line exists in the data and it
is absent from the documentation.** Its five-digit numbering block, its
all-class-12 classification, and its absence from the published plan are all
consistent with the vendor treating cross lines as QA product rather than
deliverable coverage. **How many there are, and where, is answerable only by
enumerating `point_source_id` across tiles** – there is no index, and one line in
one tile is all we have seen. That is worth doing systematically: at Elba,
one cross line reaches lines 136, 137 and 138 and breaks the across-track
degeneracy outright.

### 2.4 gen2's geometry at Elba, which is perpendicular to gen1's

From `SwathPolygon.gdb` (see §4): the swaths covering the Elba reference point are
`PS_ID` **3042** and **3043**, both from lift **`Day12121_TM515`**. Fitting the
long axis of each swath polygon by SVD:

| `PS_ID` | lift | axis bearing | swath width | northing centre |
|---|---|---|---|---|
| 3040 | Day12121_TM515 | 89.8° | 1 453 m | 4 885 888 |
| 3041 | Day12121_TM515 | 90.0° | 1 444 m | 4 884 932 |
| 3042 | Day12121_TM515 | 90.8° | 1 432 m | 4 883 964 |
| 3043 | Day12121_TM515 | 90.1° | 1 441 m | 4 882 984 |

Spacing between swath centres **~965 m**, swath width **~1 440 m**, sidelap
**~33%** – almost exactly gen1's geometry, rotated 90°.

**And the two swaths over Elba were flown in opposite directions**, a
boustrophedon like gen1's. The swath polygons' `START_TIME` / `END_TIME` are
adjusted-standard GPS seconds; converting (GPS epoch 1980-01-06, GPS − UTC = 18 s
in 2021) gives `PS_ID` 3043 at **14:55:29–15:07:14 UTC** and `PS_ID` 3042 at
**15:07:14–15:20:17 UTC** on 2021-05-01. The flight log for that lift lists
**line 43, direction E, 14:57–15:06** and **line 42, direction W, 15:09–15:19**.
So `PS_ID` = 3000 + the log's block-3 line number, and 3042/3043 are an
east-then-west pair. (The swath times are contiguous – each swath's `END_TIME` is
the next one's `START_TIME` – so they are segmentation boundaries rather than
on-line windows; the identification rests on the direction and the ten-minute
duration, both of which match.)

**Verified against our own data:** `data/after/3dep2021_subpatch.laz` contains
`point_source_id` **3041, 3042, 3040** and nothing else. The vendor's `PS_ID` is
our `point_source_id` with no remapping.

`Day12121_TM515` decodes as day 121 of 2021 = **2021-05-01**, sensor Leica
TerrainMapper serial 90515 – and the flight logs (§4) contain a
2021-05-01 / TM 90515 log whose lines are all E and W. **Elba's gen2 was flown
east–west on 2021-05-01 by TM 90515**, which is consistent with the
green-up/leaf-on finding already in the project memory.

**This matters for the analysis.** gen1's across-track axis at Elba is east–west;
gen2's is north–south. Any error that lives in the across-track coordinate of one
epoch appears as an *along-track-invariant* pattern in the other. The two epochs
are, geometrically, each other's cross lines.

---

## 3. Swath adjustment methodology, both epochs (Priority 3)

### 3.1 gen1, verbatim

`lidar_semn2008.xml`, Lineage, Lidar Processing, steps 5–6:

> "5. Inspected for calibration errors in the dataset using the TerraMatch
> software. This was accomplished by sampling the data collected across all flight
> lines and classifying the individual lines to ground. The software used the
> ground-classified lines to compute corrections (Heading, Pitch, Roll, and Scale).
> 6. Orientation corrections (i.e., calibration corrections) were then applied to
> the entire dataset."

Software versions, same record: "TerraSolid TerraScan (version 009.010),
TerraModeler (version 009.002) and TerraMatch (version 009.003) and Intergraph
MicroStation (version.08.01.02.15)."

**Read this precisely.**

- The corrections are **orientation** corrections – heading, pitch, roll, scale.
  There is **no vertical shift** among them, and no mention of a per-line dz.
- They were "applied to the **entire dataset**". The wording describes a single
  solution, not a per-line or per-lift block adjustment. It does not explicitly
  exclude a per-lift solution, so I will not claim it did not happen – but nothing
  in the record says it did.
- **No boresight or misalignment values are published**, for either sensor.
- **No relative-accuracy statistic is published at all** – no interswath RMSDz, no
  swath separation figure, no smooth-surface repeatability. The only accuracy
  numbers in the whole 2008 record are absolute, against control.

**Consequence for our work, stated as a consequence and not as a finding:** the
per-swath constants we solve are not re-deriving something the vendor documented
solving. The vendor documents solving *angles*, once, for everything. A residual
per-line vertical offset is exactly what that procedure leaves behind. Whether an
across-track term should be expected is not settled by the documentation either
way: a global scale and roll correction would remove the *mean* across-track ramp
while leaving per-line departures from it.

### 3.2 gen2, verbatim

`MN_SE_Driftless_2_2021_Lidar_Mapping_Report.pdf`, **p. 12**, §2.1.6 Geometric
Calibration:

> "Laser point position was calculated by associating the SBET position to each
> laser point return time, scan angle, intensity, etc. Raw laser point cloud data
> was created for the whole project area in LAS format. **Automated line-to-line
> calibrations were then performed for system attitude parameters (pitch, roll,
> heading), mirror flex (scale) and GPS/IMU drift.** Statistical reports were
> generated for comparison and used to make the necessary adjustments to remove any
> residual systematic error. […] Software used included proprietary software,
> TerraMatch v20, and Leica CloudPro 1.2.4., Optech's LiDAR Mapping Suite (LMS)
> v4.5"

So gen2 is **line-to-line** where gen1 is dataset-wide, and gen2 additionally
solves **GPS/IMU drift**. Same software family, thirteen years apart.

**And gen2 publishes the resulting relative accuracy.** Same report, p. 12, §2.1.7:

> "Interswath (overlap) consistency was assessed at multiple locations within
> overlap in non-vegetated areas containing only single returns and located in
> areas with slopes of less than 10-degrees. […] These overlap areas include
> adjacent, overlapping parallel swaths within a project, **cross-tie swaths**, and
> a sample of intersecting project swaths in both flight directions, and adjacent,
> overlapping lifts. […] This project required the interswath accuracy to meet ≤
> 8-cm RMSDz."

Table 3.4.1 gives 15 test polygons, RMSDz **0.006–0.018 m**. The shapefile
`interswath.shp` (1,420 bytes) carries the same 15 records with geometry, so each
number has a location. **Test polygon FID 14 sits 1.8 km from the Elba reference
point** (centroid E 579 929, N 4 885 413): MIN −0.022 m, MAX +0.021 m,
**RMSDz 0.009 m**. That is the vendor's own between-flight-line agreement
essentially at our site.

Intraswath (p. 14, §2.2): 3 polygons, RMSDz 0.017–0.033 m, requirement ≤ 6 cm.

**Boresight values: published for gen2, at the factory.**
`Attachment1_Sensor_Calibration_Reports.pdf` (1,833,652 bytes) is three Leica
TerrainMapper calibration certificates. For serial **90515** – the sensor that flew
Elba – p. 3, §5.1 LiDAR Geometric Calibration Results, dated 12 December 2018:

> IMU Misalignment ω −0.022555°, Φ 0.056357°, κ 0.000504°
> Boresight Θ 0.015419°, Φ −0.001923°
> Wedge 0: Δ Alpha −0.043014°, Wedge Position Δ Offset 0.442789°, Position
> Correction X −0.012826° Y 0.000012°, Mount Roll 0.045379° Pitch 0.210132°,
> Rotation Axis Roll 0.031087° Pitch 0.076675°
> Wedge 1: Δ Alpha −0.005517°, Wedge Position Δ Offset 0.559649°, Position
> Correction X 0.030760° Y −0.001169°, Mount Roll 0.012366° Pitch 0.054254°,
> Speed Pitch 1.50E-06 °/rps², Rotation Axis Roll 0.032485° Pitch −0.029191°

Certificates are also present for serials **91511** (3 July 2019) and **91557**
(1 July 2020). **These are factory values, not the project's in-flight TerraMatch
solution**, and the report's §1.8 lists the serials as "90515 … 90557 … 90511",
which does not match the certificates' 90515 / 91511 / 91557 – a transcription
slip in the report, worth not tripping over. The **project's** line-to-line
correction values are not published for gen2 either.

---

## 4. The 2021 project's documents (Priority 4)

Root: `https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/MN_SE_Driftless_2021_B21/`

The project splits into five work units. **Elba is in `MN_SEDriftless_2_2021`,
Work Unit ID 222535, Quality Level 1** – verified two ways: our fetch script pins
`--base .../usgs-lidar-public/MN_SEDriftless_2_2021`
(`scripts/fetch_3dep_curl.py:19`), and the cached EPT boundary index has Elba
inside `MN_SEDriftless_2_2021` and outside the other four.

Each work unit `MN_SEDriftless_<n>_2021/` holds `breaklines/`, `reports/` and
`spatial_metadata/`. The `reports/` tree, which was not previously opened:

| file | what it holds |
|---|---|
| `MN_SE_Driftless_<n>_2021_Lidar_Mapping_Report[_WU…].pdf` | **the main document** – flight planning, sensors, planned flight specifications, timeline, GNSS/IMU, trajectory quality, geometric calibration, interswath and intraswath results, classification, hydro-flattening, DEM, accuracy |
| `Attachment1_Sensor_Calibration_Reports.pdf` | factory calibration certificates per sensor, with misalignment and boresight angles (WU4 has none) |
| `Attachment2_Flight_Logs_….pdf` | **21 Woolpert Lidar Acquisition Logs** for WU2 – date, aircraft, tail number, pilot, operator, sensor serial, airport, weather, air speed, altitude AGL/MSL, point spacing, density, FOV, scan frequency, pulse rate, laser power, then per-line: line number, direction, start/end UTC, satellites, PDOP, comments |
| `Attachmnet3_GPS IMU Images_WU…/` | trajectory plots (sic, the directory name is misspelt on the server) |
| `USGS_MN_SEDriftless_<n>_2021_Summary_Report.pdf` | NGTOC data validation: accept/reject per deliverable class, no new numbers |
| `USGS_…_FINAL_{DEM,LPC}_Report.txt` | per-file validation logs |
| `vendor_provided_xml/*.xml` | FGDC metadata for point cloud, DEM, intensity, breaklines |

And `spatial_metadata/contractor_provided/`: `SwathPolygon.gdb`,
`interswath.shp`, `intraswath.shp`, tile index or DPA, `intensity_imagery/`,
`maximum_surface_height_rasters/`, `swath_separation_images/`.

Project-level: `USGS_MN_SE_Driftless_2021_B21_Project_Report.pdf` (501,262 bytes,
3 pages) – accuracy table per quality level, classification list, "Sensor(s) Used:
**Optech Galaxy T2000 - Aerial Oscillating Mirror**", and the five work units with
their IDs, geoid model (GEOID18), GSD (0.5 m) and collection dates. **That sensor
line is incomplete**: the work-unit report and the flight logs show three Leica
TerrainMappers doing most of the flying, and Elba specifically was flown by
TerrainMapper 90515, not by the Galaxy. Do not take the project report's sensor
field as the sensor over any given tile.

### 4.1 The swath polygons are the key artefact

`SwathPolygon.gdb` is an Esri file geodatabase served as loose files. I fetched 25
of its component files (**737,088 bytes total**) and GDAL's OpenFileGDB reads the
local copy directly. One layer, `SwathPolygon`, **181 features** – matching the
report's "A total of 181 individual flight lines were collected" (p. 9) – with
fields:

`WP_ID, WU_ID, LIFT_ID, PS_ID, TYPE, START_TIME, END_TIME`

**`PS_ID` is `point_source_id`.** Verified: the two swaths containing the Elba
reference point are `PS_ID` 3042 and 3043, and our own
`data/after/3dep2021_subpatch.laz` contains `point_source_id` 3040, 3041, 3042.

`LIFT_ID` groups swaths into sorties (`Day12121_TM515` = day-of-year 121 of 2021,
sensor TM 90515), which joins straight to the flight logs. `START_TIME` /
`END_TIME` are adjusted-standard GPS seconds.

**`TYPE` is `PROJECT` on all 181 features, and on every swath of all five work
units.** The field carries a coded-value domain, and reading that domain out of
the geodatabase's `GDB_Items` table gives its allowed values:

> `<DomainName>TYPE</DomainName> … <CodedValue><Name>PROJECT</Name></CodedValue>
> <CodedValue><Name>CROSS-TIE</Name></CodedValue>
> <CodedValue><Name>FILL-IN</Name></CodedValue>
> <CodedValue><Name>CALIBRATION</Name></CodedValue>
> <CodedValue><Name>OTHER</Name></CodedValue>`

So the vendor's own schema has a **`CROSS-TIE`** code, and **it is never used**. I
fetched the feature table (`a00000009.gdbtable`) from each of the other four work
units – 557,422 / 711,266 / 711,266 / 403,217 bytes – and scanned each for the
domain strings: `PROJECT` only, in all four. **The 2021 project delivered no
swath labelled as a cross-tie, in any work unit**, despite §2.1.7's boilerplate
listing cross-tie swaths among the kinds of overlap assessed. That boilerplate
should not be read as evidence that cross ties were flown here.

### 4.2 Flight geometry for gen2, from the report and the logs

Planned specifications, report p. 8–9, §1.8, for the Leica TerrainMapper sensors
(90515, 90557, 90511 as printed):

> "Maximum Number of Returns: 15; Nominal Point Spacing: 0.35-m; Nominal Point
> Density: 8 ppsm; **Flying Height Above Ground Level: 2,000-m**; Flight Speed:
> 160-knots; **Scan Angle: 40°**; Scan Rate Used: 150-Hz; Pulse Rate Used:
> 1600-kHz; Multi-Pulse in Air: Enabled; **Overlap: Minimum 25%**"

and for the Optech Galaxy: 8 returns, 0.35 m, 8 ppsm, **1,400 m AGL**, 120 knots,
**45°**, 90 Hz, 800 kHz, minimum 20% overlap.

Report p. 9, §1.9: "Lidar data was collected from May 31, 2021, through May 17,
2022. A total of 181 individual flight lines were collected." **That start date
contradicts the same report's own acquisition window** ("April 25, 2021, through
May 15, 2022", p. 6), the vendor XML ("April 25, 2021, through May 15, 2022"), the
USGS project report (collection start 2021-04-25) and **the flight logs, which
begin 04/25/2021**. Elba's lift is 2021-05-01. **Treat the p. 9 sentence as an
error**; the logs and the swath `START_TIME` are the reliable sources.

The 21 flight logs for WU2 cover 04/25, 04/26, 04/29, 04/30, 05/01, 05/02, 05/11,
05/12, 05/26 of 2021 and 05/06, 05/07, 05/15 of 2022, on sensors TM 90515, TM
90511 and TM 557, aircraft Cessna 404 Titan N475RC and Reims 406 N406SD, out of
KRST (Rochester). Line directions across the block: **N 39, S 42, E 64, W 67** –
the work unit is flown north–south in some sub-blocks and east–west in others.
The 2021-05-01 TM 90515 log is 19 lines, all E or W – Elba's lift.

**Elba's own lift, from the log** (log 10 of 21): Project # 81926, "Minnesota SE
Driftless_Block 3", Unique ID **`Day121_90515_A`**, flight date **05/01/2021**
(day of year 121), Cessna 404 Titan **N404CP**, pilot Dar Perl, operator
Galambos, out of and back to **KRST**, 14:35–19:17 UTC. Settings: FOV **40°**,
scan frequency **150 Hz**, pulse rate **1600 kHz**, laser power 100%, air speed
160 kt, **altitude 6,562 ft AGL / 7,083 ft MSL** (2,000 m AGL). Conditions clear,
18 °C, wind 230° at 18 kt gusting 24. Twenty-one lines numbered 25–45, alternating
E and W, PDOP **1.1–1.3**, 18–23 satellites; the operator noted "93-95% Return
rate". Elba sits under lines 43 and 42 of that log.

GNSS base stations, p. 10, Table 2.5.2: **MNLS_CORS** (44°26'28.13634",
−93°54'24.61955", ellipsoid height L1 phase centre 239.951 m) and **MNSV_CORS**
(43°54'09.04784", −92°28'55.92990", 362.281 m), plus Applanix PP-RTX. Trajectory
software POSPac 5.3, IPAS Pro 1.35, Novatel Inertial Explorer 8.60.6129; goals
combined separation < 10 cm and average PDOP < 3.0. Elba's log shows PDOP 1.1–1.6.

---

## 5. Everything else load-bearing (Priority 5)

### 5.1 Sensors, rates, heights

| | gen1 (2008) | gen2 (2021, Elba) |
|---|---|---|
| sensor | Optech ALTM Gemini (AeroMetric) **and** Leica ALS50-2 (Surdex) | Leica TerrainMapper, **serial 90515** |
| serial | **not published** | 90515, factory-calibrated 12 Dec 2018 |
| flying height | 2 400 m above mean terrain | 2 000 m AGL planned; Elba's log gives 6 562 ft AGL / 7 083 ft MSL |
| field of view | 32° full | 40° full |
| pulse rate | **not published** | 1 600 kHz |
| scan frequency | **not published** | 150 Hz |
| nominal spacing | 1.0 m | 0.35 m, 8 ppsm |
| overlap | "Sidelap - 60%" **(inconsistent, see §2.2)** | minimum 25% planned; ~33% measured at Elba |
| line spacing | 962.8 m at Elba (planned, published; 43.31" of longitude) | ~965 m (measured from swath polygons) |
| lines at Elba | 133–138, plus cross line 10010 | 3040, 3041, 3042 |
| flown | 2008-11-25 | 2021-05-01 |
| trajectory software | Applanix POSGPS/POSProc (AeroMetric); GravNav GNSS + Leica IPAS (Surdex) | POSPac 5.3, IPAS Pro 1.35, Inertial Explorer 8.60 |
| base stations | ≥2 MnDOT CORS per lift, from BLUE, CLDN, DDGC, ELKT, EYTA, LCHI, LCRS, NALB, PRSP, REDW, RSHF, STWV, TWNL, WBSH, WINO, WSCA; max baseline "Not greater than 30km" | MNLS_CORS, MNSV_CORS, plus Applanix PP-RTX |

### 5.2 Classification routine and parameters

**gen1**, `lidar_semn2008.xml` step 7, in full:

> "Automatic ground classification was performed using algorithms with customized
> parameters to best fit the project area. Several areas of varying relief and
> planimetric features were inspected to verify the final ground surface."

**No algorithm is named and no parameter is given.** Delivered classes are 0, 2,
5, 6, 8, 9, 10, 12. This is the whole of what is documented about gen1's ground
classification – which is worth stating plainly, because our CSF reconstruction of
gen1 ground has no vendor recipe to be compared against.

**gen2**, report p. 15, §2.3: initial classification on "first and only" as well
as "last of many" returns, then TerraScan v20 / TerraModeler v20 / GeoCue LP360 /
Global Mapper v20, manual review, a **0.7 m buffer** around hydro-flattened
features reclassified to Ignored Ground (class 20). Classes 1, 2, 7, 9, 17, 18, 20.

### 5.3 The definition of "Surface Z" in the 2008 validation tables

**This one is not settled by the documents, and it matters for §1.**

The per-county validation report PDFs contain **no narrative at all** – only the
point tables (`Name, Control X, Control Y, Control Z, Surface Z, Error, Z-Diff
Squared, Absolute Error`), per-class RMSE and NSSDA(95%) summaries, and two charts.
I re-read `winona_val.pdf` end to end looking for a method statement; there is
none. "Surface Z" is never defined in the document that uses it.

The only evidence is the dataset metadata's wording, which describes both
accuracy tests as being **of the TIN**:

> "The Fundamental Vertical Accuracy (FVA) **of the TIN** achieved 0.161 meters …"
> "The Consolidated Vertical Accuracy (CVA) **of the TIN** as tested by MnDNR
> achieved 0.287 meters at a 95% confidence level of all land cover categories.
> 1009 control points covering the 5 land classes were used in this evaluation."

So the best available reading is that the DNR's `Surface Z` is **interpolated from
a TIN of the delivered ground points**, not sampled from the 1 m DEM and not the
nearest point. **I mark this as inferred from the metadata's wording, not stated
in the validation report.** It is worth pinning down before the −22 mm intercept
of `ADDITIONAL_GROUND_CONTROL.md` §5 is reconciled against our own +22.7 mm
anchor, because a TIN interpolation and a nearest-point comparison have different
biases on sloping ground – and note that AeroMetric's *own* bias determination
used the **nearest point** (§1.1), so the two 2008 numbers are not even
like-for-like with each other.

The 1 009 / 1 033 discrepancy already flagged in `ADDITIONAL_GROUND_CONTROL.md` §5
is confirmed here: the metadata's narrative says "1009 control points" while its
per-county list sums to 1 033 (121 + 128 + 134 + 115 + 125 + 137 + 97 + 176 = 1 033
for eight counties; Freeborn is in the coverage list but has no per-county RMSE).
**§1.4b makes it worse, not better:** the FEMA certification letter publishes a
ninth county (Freeborn, 0.144 m / 126) and a transposed Houston, and its nine
counties sum to 1 135. Three published copies, three different count lists, one
shared claim of 1 009.

### 5.4 A note on what the 2008 project was for

> "The Southeast Minnesota lidar project's goal was to provide high accuracy,
> bare-earth processed lidar data suitable for the FEMA National Flood Insurance
> Program. In August 2007, the nine counties had experienced a major rainfall
> event with over 10 inches of rain…"

Relevant because FEMA NFIP work drives the accuracy specification and the paper
deliverable list, and because it suggests FEMA Region V may hold a copy of the
accuracy assessment report.

---

## 6. What is NOT in the documentation, and how to get it

Ranked by what it would buy us.

1. **The value of the 2008 bias adjustment.** Not published anywhere reachable.
   Ask MN DNR / MnGeo. This is the number that gates the absolute datum.
2. **The value of the 2021 bias adjustment** – the same quantity for gen2, and I
   had not previously realised it existed. Ask USGS NGTOC or Woolpert. Without it,
   the DoD's absolute level is the difference of two unknown vendor constants.
3. **The 2008 "Lidar Accuracy Assessment Report"** – "One paper copy" per the
   metadata's deliverable list, so plausibly never digitised. Ask MN DNR.
4. **The 2008 cross-line plan** – how many, where, which lifts. The data proves
   at least one exists; no index lists any. Ask MN DNR, or enumerate
   `point_source_id` across tiles ourselves.
5. **Which 2008 lines were flown by AeroMetric's Optech ALTM Gemini and which by
   Surdex's Leica ALS50-2.** Two sensors merged into one cloud with no flag. Ask
   MN DNR; the AeroMetric report would say.
6. **The 2008 TerraMatch correction values** (heading, pitch, roll, scale) and
   whether they were solved once or per lift. Ask MN DNR.
7. **The 2021 project's in-flight line-to-line calibration values.** The factory
   certificates are published; the project solution is not. Ask USGS NGTOC.
8. **The 127 AeroMetric QA/QC checkpoint coordinates** – already on the request
   list in `ADDITIONAL_GROUND_CONTROL.md` §5.2.
9. **The definition of "Surface Z"** in the DNR validation tables (§5.3). Ask MN
   DNR; it may be answerable in one sentence by whoever ran the test.
10. **The July 2008 RFP and the AeroMetric contract.** The RFP is the document that
    would say whether cross lines were *required*, and what the flight plan had to
    deliver. No copy was found. Ask MN DNR.
11. **Which per-county control count is right** – the MnGeo record's "Houston
    0.110m, 134" or the certification letter's "Houston 0.134m, 110" – and why
    neither list sums to the stated 1 009 (§1.4b). Parsing
    `Houston_county_validation_report.pdf` settles the first half ourselves; the
    second half needs asking.

**The precedent to lead with has got stronger.** USGS's own legacy tree publishes
`.../metadata/legacy/MN_PINECO_2006/Pine_LIDAR_Final_Report.doc` – the vendor final
report for a 2006 Minnesota project – in exactly the directory where SE MN's would
sit if it had been deposited. Together with the published
`lidar_checkpts_pine2007` shapefile, that is two Pine County precedents for
publishing precisely what is being asked for here.

`ADDITIONAL_GROUND_CONTROL.md` §5.2 already drafts an email to Sean Vaughn
(MN DNR lidar data steward, sean.vaughn@state.mn.us, 763-284-7223) and MnGeo
(gisinfo.mngeo@state.mn.us, 651-201-2499) covering items 1, 3 and 8. **Add items
4, 5, 6 and 9 to it** – the paragraph below is drafted to drop in after its
numbered list:

> Four further things about the 2008 acquisition would help, all of which the
> AeroMetric report or the delivery records would answer:
>
> 4. Whether cross or tie lines were flown, and if so how many and where. The
>    delivered point cloud contains at least one east–west line over Winona County
>    (point source id 10010, flown 25 November 2008), but it does not appear in the
>    published flight-line layer, so I have no way to find the others.
> 5. Which flight lines were flown by AeroMetric's Optech ALTM Gemini and which by
>    Surdex's Leica ALS50-2. The metadata says the two point clouds were merged,
>    but the delivered files carry no sensor flag.
> 6. The TerraMatch orientation corrections (heading, pitch, roll and scale) that
>    process step 6 says were applied to the entire dataset, and whether they were
>    solved once for the project or separately per lift.
> 9. How "Surface Z" in the county validation tables was obtained – interpolated
>    from a TIN of the ground points, sampled from the one-metre DEM, or taken as
>    the nearest lidar point. The reports give the column but not the method, and
>    the answer changes how I read the residuals on sloping ground.
>
> Two small things you may want to know. First, the metadata says the data was
> acquired between 18 and 24 November 2008, but the flight-line layer you publish
> dates lines to 25, 26 and 27 November as well, including all of the lines over
> Winona County that I am working with. Second, the per-county control counts do
> not agree between sources: the Geospatial Commons metadata lists eight counties
> summing to 1,033 points and gives Houston as 0.110 m over 134 points, while the
> 2012 FEMA certification letter lists nine counties summing to 1,135 and gives
> Houston as 0.134 m over 110 points. Both say the total was 1,009. The
> certification letter is also the only place I can find a figure for Freeborn
> County.

**For the 2021 side**, a separate request to USGS (tnm_help@usgs.gov, quoting
contract G16PC00022, task order 140G0221F0253, Work Unit ID 222535) should ask for
(a) the magnitude of the vertical adjustment described on p. 15 of the Lidar
Mapping Report, and (b) the project's line-to-line calibration solution. A draft:

> Subject: MN SE Driftless 2021 B21, WU 222535 – vertical adjustment magnitude and
> line-to-line calibration values
>
> I am using MN_SEDriftless_2_2021 together with the 2008 Minnesota DNR lidar to
> measure landscape change near Elba, Winona County, and I am trying to account for
> the absolute vertical datum of both epochs.
>
> Section 2.3 of the Lidar Mapping Report (Woolpert, November 2022, WU 222535,
> page 15) states that "Based on the statistical analysis, the lidar data was then
> adjusted to reduce the vertical bias when compared to the survey ground control of
> higher accuracy." Two requests:
>
> 1. The magnitude and form of that adjustment – a single constant for the work
>    unit, or per lift, or per swath – and its value.
> 2. The line-to-line calibration values from section 2.1.6 (pitch, roll, heading,
>    mirror flex and GPS/IMU drift). The factory certificates in Attachment 1 are
>    published; the project solution is not.
>
> Either would be useful on its own. I am happy to say what I find.

**Not a Data Practices Act matter on either side** – both are ordinary public
records requests, and the 2008 one is a retrieval problem rather than a disclosure
problem.

---

## Appendix A – what was fetched

All requests spaced 1–3 s apart, one at a time. **No lidar was downloaded**; every
point cloud read was already on disk. Total fetched by me: **about 10.7 MB**;
the sub-agent's archive sweep is marked separately.

| what | host | bytes |
|---|---|---|
| MnGeo `documentation/`, `lidar/`, `projects/`, `tile_index/`, `examples/`, `q250k/`, `q250k/q4342/`, `county/winona/{,laz/,geodatabase/}` listings | resources.gisdata.mn.gov | ~40 kB |
| `lidar_semn2008.xml` + `.html` | resources.gisdata.mn.gov | 42,923 |
| `raw_LiDAR_Data_README.rtf`, `about_this_data.rtf`, `readme_first.txt` | resources.gisdata.mn.gov | 101,248 |
| **`flight_lines_and_dates.zip`** | resources.gisdata.mn.gov | **759,353** |
| USGS staged `metadata/` and `LPC/Projects/` indexes | rockyweb.usgs.gov | 229 kB |
| `MN_SE_Driftless_2021_B21` root + 5 work-unit + 10 sub-directory listings | rockyweb.usgs.gov | ~90 kB |
| `USGS_MN_SE_Driftless_2021_B21_Project_Report.pdf` | rockyweb.usgs.gov | 501,262 |
| **`MN_SE_Driftless_2_2021_Lidar_Mapping_Report.pdf`** | rockyweb.usgs.gov | **1,123,464** |
| `Attachment1_Sensor_Calibration_Reports.pdf` (WU2) | rockyweb.usgs.gov | 1,833,652 |
| `Attachment2_Flight_Logs_…_WU222535.pdf` | rockyweb.usgs.gov | 1,388,800 |
| `USGS_MN_SEDriftless_2_2021_Summary_Report.pdf` | rockyweb.usgs.gov | 164,915 |
| `MN_SEDriftless_2_2021_Classified_Point_Cloud_Metadata.xml` | rockyweb.usgs.gov | 12,522 |
| `SwathPolygon.gdb` (WU2) – 25 component files | rockyweb.usgs.gov | 737,088 |
| `SwathPolygon.gdb/a00000009.gdbtable` for WU1, 3, 4, 5 | rockyweb.usgs.gov | 2,383,171 |
| `interswath.{shp,shx,dbf,prj}` | rockyweb.usgs.gov | 3,057 |
| HEAD requests for sizing (≈20) | rockyweb, resources.gisdata | ~0 |
| NOAA S3 bucket listings (2), `chs.coast.noaa.gov` probes (3) | amazonaws, noaa.gov | ~40 kB |
| `metadata/legacy/` index, `legacy/MN_SEMN_2008/`, `legacy/MN_PINECO_2006/` listings | rockyweb.usgs.gov | 115 kB |
| `legacy/MN_SEMN_2008/semn_metadata.{xml,docx}` | rockyweb.usgs.gov | 33,675 |
| `mn_lidar_certification-se_mn_10-22-2012.pdf` | files.dnr.state.mn.us | 81,131 |
| `elev_09feb.pdf` (Digital Elevation Committee minutes) | mngeo.state.mn.us | 13,514 |
| Wayback/archive searches, county and vendor sites, ScienceBase, MN GIS/LIS | various | *(sub-agent)* |
| NOAA InPort 68818, MnGeo `lidar_2008-2012.html` | fisheries.noaa.gov, mngeo.state.mn.us | *(WebFetch)* |
| gen1 tiles at Elba, gen2 subpatch, county validation PDFs | **local disk / prior session** | 0 |

Scratchpad artefacts (session-scoped): `docmine/` holding `fl/` (the flight-line
shapefile), `SwathPolygon.gdb/`, `isw/`, the four gen2 PDFs and their `pdftotext
-layout` extractions, `semn2008.txt`, `vendor_lpc.txt`.

## Appendix B – verification status

**Verified by me, this session:**

| claim | how |
|---|---|
| 204 SE Minnesota flight lines, all N–S, zero within 45° of E–W | shapefile read, heading histogram |
| line spacing 962.84 ± 0.16 m at Elba, and a constant 43.3063" of longitude project-wide | 145 consecutive gaps at N 4 883 678, repeated at two other northings, then converted to geographic |
| `point_source_id` = `OBJECTID` − 2644 | 12 lines, mean −6.9 m, sd 18.9 m against 962.8 m spacing |
| Elba lines 133–138 flown 2008-11-25 | `Date_Flown` field **and** gen1 GPS week-seconds → Tuesday 18:00–20:38 GPS |
| the `point_source_id` ↔ `OBJECTID` mapping survives a date-boundary test | the layer puts line 128 on 11-26 and line 129 on 11-25; the data's GPS day-of-week flips from 2.964 (Tue) to 3.068 (Wed) between exactly those two lines, in `4342-29-61` |
| metadata's "November 18 - 24" window is too narrow | 67 SE MN lines dated 11-25 to 11-27 |
| "Sidelap - 60%" is inconsistent with the same record's 32° FOV and 2 400 m height at 962.8 m spacing | arithmetic, shown in §2.2 |
| gen2 at Elba is `PS_ID` 3040/3041/3042, lift `Day12121_TM515`, bearings 89.8–90.8° | swath polygons read with GDAL; `point_source_id` read from `data/after/3dep2021_subpatch.laz` |
| 181 swaths in WU2, all `TYPE = PROJECT`; the same for all five work units | full attribute dump (WU2) + domain-string scan of each work unit's feature table |
| interswath polygon FID 14 is 1.8 km from Elba at RMSDz 0.009 m | `interswath.shp` geometry + attributes, cross-checked against Table 3.4.1 |
| Elba is in WU2 | `scripts/fetch_3dep_curl.py:19` **and** `data/ept_boundaries.json` containment |
| gen2 had a vertical bias adjustment applied | Lidar Mapping Report p. 15, verbatim |
| the 143 LCPs calibrated gen2; NVA/VVA were held out | vendor FGDC XML, Ground Conditions, verbatim |
| no `MN_SEMN_2008` on rockyweb | 404 on two paths + absent from both directory indexes |
| the county validation reports contain no method narrative | full text search of `winona_val.pdf` |
| `MN_SEMN_2008` **is** on rockyweb, under `legacy/` | directory listed, HTTP 200; `semn_metadata.xml` and `.docx` downloaded and read |
| USGS publishes `Pine_LIDAR_Final_Report.doc` for MN_PINECO_2006 and nothing equivalent for SE MN | both legacy directories listed |
| the FEMA certification letter gives Freeborn 0.144 m / 126 points and transposes Houston relative to MnGeo | PDF downloaded, text extracted, compared county by county |
| the DEC minutes of 2009-02-26 say nothing quantitative | PDF downloaded and searched |

**Explicitly unverified:**

| claim | status |
|---|---|
| the value of either epoch's bias adjustment | **not found. Must be requested.** |
| whether gen1's TerraMatch correction was one solution or per lift | the record says "applied to the entire dataset" and nothing more |
| whether the 2008 cross line 10010 is one of many | only one is on disk; there is no index, so the count is unknown |
| "Surface Z" = TIN interpolation | **inferred** from the metadata's "of the TIN" wording; the validation reports define nothing |
| which 2008 lines are Optech and which are Leica | no public source assigns sensors to lines |
| whether Houston is 0.110 m / 134 pts or 0.134 m / 110 pts | two published sources disagree; `Houston_county_validation_report.pdf` would settle it and has not been parsed |
| the 2009 archived DNR project page and the MN GIS/LIS article | **sub-agent findings, not re-fetched by me** |
| the absence of the 2008 RFP anywhere public | sub-agent's search; a negative I have not independently repeated |
| the `point_source_id` = `OBJECTID` − 2644 mapping outside lines 128–145 | verified only over the Elba corridor; extrapolation elsewhere |
| the p. 9 "May 31, 2021" collection start | contradicted by four other sources; treated as a report error |
