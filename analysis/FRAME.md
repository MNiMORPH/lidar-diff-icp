# FRAME — the current state. Read this first; the dated FRAME_*.md files are superseded.

Written 2026-09-06. Verify every structural claim against git and the files before relying
on it: this is a map of ROLES and NEXT ACTIONS, not a results log.

    ./lidar-icp/bin/python -m pytest -q        expect 432 passed
    git rev-list --count origin/main..HEAD     expect ~115, all UNPUSHED
    lidar-diff-workflow --tile data/derived/elba --check

## The pipeline is six steps

    completeness -> base -> slope -> ridge_mask -> convexity -> curvature

`base` is `pipeline.difference_dem` via `run_all_sites`, and it produces the shipped
product: `dod.npy`, `lod.npy`, `z_after.npy`, `change.npy`, the three GeoTIFFs,
`corrections.json`, `regions.json`. Everything else is downstream.

**It takes no cloud arguments** — `base` gets them from the Site record in `sites.py`.

`lidar-diff-workflow --tile <dir> --run` executes it, stopping at the first failure, and
refuses before running anything if a selected step lacks an argument. `--only`, `--force`,
`--dry-run`.

## Alongside, reached over to, not part of it

    --with vegetation_correction     10 steps: class2_spread, q2_fit, dod_cover, lod_cover,
                                     and their feeders pfs_cover, gen1_angles, beam_table,
                                     nearground, nearground_split, canopy_struct

The correction is MEASURED AND NOT ADOPTED: on open ground the calibrated curve was WORSE
than the median (RMS 52.5 vs 49.1 mm, held out on 227 NVA marks). Its outputs are
`dod_cover_q2.npy` / `lod_cover_q2.npy`, never `dod.npy` / `lod.npy`.

The split was measured, not judged: a step is alongside when its products are read only by
`q2cover.py`, the correction's own library. `slope`/`ridge_mask`/`convexity`/`curvature`
look similar but are read by `refcells.py`, which 20+ scripts use, so they are pipeline.

`analysis/tools/` holds nine hand-run tools, declared in `workflow.TOOLS` with each
script's own docstring line.

## The product / record split — measured 2026-09-06, and the reason the package moved

`analysis/` held 157 scripts. **30 were reachable from the pipeline, tests or tools; 127
were not.** The boundary ran through the middle of every directory, which is why the folder
names stopped telling anyone anything:

    analysis/ridgelines                61 scripts -- 3 reachable
    analysis/slope_bias                14         -- 1
    analysis/investigations            24         -- 0
    analysis (top level)               34         -- 8
    analysis/tools, modules/veg_corr    9 + 9     -- all

The three reachable ridgeline scripts were the producers for `ridge_mask`, `convexity` and
`curvature`: **half the six-step pipeline was inside a folder that reads as a concluded
investigation.** They are now `lidar_diff_icp.steps`, run as `python -m`. Verified
byte-identical on elba across the move, all 11 products -- and the baseline was checked by
re-running the steps BEFORE any move, so the gate is known to mean something.

`script_of` had to learn `-m` at the same time: it found producers by matching `.py` paths
in the command string, so a module command names no source and the step silently loses
CODE-STALE detection. `_MODULE_RE` resolves a `-m` module name to its file under `src/`.

**Why the package and not an `analysis/lib`:** a test reached its subject with a
`sys.path.insert` naming the subject's directory, so a move was a code change and the move
plus the test fix could not be split. In the package the insert is DELETED, not repointed.
That closed the blocker the previous plan was stuck on. 38 such hacks remain under
`analysis/`, and `analysis/` still has no `__init__.py` anywhere.

**One step keeps a file-path invocation:** `pfs_cover` runs under the conda env, which does
not have `lidar_diff_icp` installed (verified: ModuleNotFoundError). Its module imports
nothing from the package, so the path form is correct there, not a workaround.

## The swath alignment — see analysis/SWATH_ALIGNMENT_METHOD.md

Eight steps, and step 7 is the one that is easy to miss: the free network deliberately
leaves the absolute level free, and **the ground-control datum is what fixes it**, exactly
(42.40 mm of zero-line spread -> < 1e-9).

Since 2026-09-05 the network is weighted by **1/variance, not by cell count**. That changed
a shipped product: Battle Creek's swath 1102 moved +13.7 -> +27.1 mm, because its only edge
was a 36-cell sliver with a 43:1 lever arm returning -3.4640 m, which under `w = n` entered
at weight zero and left the swath at an `lstsq` minimum-norm value.

`across_track_tie` (the intercept tie) is OURS: the regression is USGS/ASPRS standard
(Sampath et al. 2016) but they constrain the intercept to zero and report the slope. See
`analysis/ACROSS_TRACK_TIE_LITERATURE.md`. **Open**: our single covariate estimates a
COMMON-MODE roll only; regressing on `tan θ_ref` and `tan θ_src` separately would test
whether `k` carries a differential bias.

## Sites

All six rebuilt on current code. Five are byte-identical across the reorganization
(dod/lod/z_after), which is the evidence that moving 55 files changed no science.

    site         stable_sigma  median LoD   gen2 completeness
    elba            0.055 m      0.113 m       0.978
    whitewater      0.055        0.105         0.971
    carlton         0.041        0.078         0.989
    cook            0.088        0.160         0.936
    mnrv            0.056        0.103         UNMEASURABLE (bbox spans two 3DEP projects)
    battlecreek     0.044        0.090         0.998

## ⚠ BUG: the wrong geoid model is applied at three sites (2026-09-07)

`references.geoid_difference` defaults to `before_geoid="us_noaa_geoid03_conus.tif"` and
`pipeline.py:622` calls it as `geoid_difference(bounds, 26915)` with **no per-site
override**. `Site` has no geoid field. So every site is differenced GEOID03 → GEOID18.

But gen1 is not one survey. The six sites fall in **four acquisitions on two geoids**,
each quoted from that project's MnGeo metadata page:

    site         county    project              gen1 geoid   applied   correct    ERROR
    elba         winona    lidar_semn2008       Geoid03      +67.28    +67.28     +0.00
    whitewater   wabasha   lidar_semn2008       Geoid03      +57.92    +57.92     +0.00
    mnrv         lesueur   lidar_swmn2010       Geoid03      +67.23    +67.23     +0.00
    battlecreek  ramsey    lidar_metro2011      Geoid09      +71.85    +16.98    +54.87
    cook         cook      lidar_arrowhead2011  Geoid09      +51.68    +25.29    +26.39
    carlton      carlton   lidar_duluth2012     Geoid09      +78.14    +50.39    +27.75

`geoid_datum` is ADDED to gen1, so adding too much makes gen1 read HIGH and the DoD
(`gen2 − gen1`) read LOW by the error. **Battle Creek's DoD is ~54.87 mm too low against
an LoD of 90 mm** — 61% of its own detection threshold, and more than that site's entire
zero-line lever (37.50 mm).

FIXED 2026-09-07 (commits dc20fb3, b239e9c). `acquisitions.py` records the survey, its
geoid and the sentence from its metadata page that asserts it; `Site.gen1_project` names
it; `geoid_difference`'s `before_geoid` default is REMOVED so six callers now state their
frame; `apply_datum` refuses as its first statement. Three regression tests, each shown to
bite. **battlecreek, cook and carlton ALL REBUILT** on the correct frame, pipeline steps
included; elba, whitewater and mnrv are unaffected and were not touched.

## ⚠ The along-track drift fit ABSORBS datum error (found on that rebuild)

The Battle Creek rebuild did NOT move the DoD by the 54.87 mm the datum moved. It moved by
**median +11.85 mm**, non-uniformly. The reason is in the products:

    swath   mean drift, GEOID03 -> GEOID09      absorbed
    1012        -52.16 -> -0.71 mm   +51.45       94%
    1013        -55.36 -> -9.55       +45.81      83%
    1014        -61.42 -> -18.60      +42.82      78%
    1101        -64.63 -> -28.65      +35.98      66%
                                mean +44.02       80%

The per-swath along-track drift is a data-driven `f(gps_time)` fitted against gen2, so it
absorbs any constant or slowly-varying gen1 datum error. Verified unchanged in the same
run: the gen1-internal per-swath dz and the Nuth-Kääb horizontal shift, so the absorption
is specifically the DRIFT stage.

Two consequences. **This is why the geoid bug survived** — the pipeline was 80%
self-correcting, so a wrong frame did not look like an offset. And it **contaminates the
project goal**: a reusable statewide per-swath GPS-drift correction. Battle Creek's
published drift curves carried ~50 mm of geoid bookkeeping rather than instrument drift,
and curves fitted at sites on different geoids are not comparable.

**It also absorbs a GRADIENT — the part parallel to the flight lines.** Predicted from
flight heading (x,y regressed on gps_time per swath) against the tilt-error direction,
then observed in the rebuilt DoDs:

    site         datum tilt   along-track   expected surviving   observed
    cook           12.35 mm          56%            10.23 mm      11.66 mm
    carlton         7.44             98%             1.65          0.26
    battlecreek     1.96             70%             1.40      not measured

Carlton's flight lines run 98% parallel to its tilt error, so almost none of it reached the
DoD; Cook's are 56% parallel, so most of it did. **How much datum error reaches the DoD
depends on FLIGHT HEADING.** Two sites with identical datum errors give different DoD
errors purely from the angle between their flight lines and the error's gradient — a method
artifact across sites, and it breaks the rule that compared results be produced the same
way.

**Can control marks resolve it? Andy's question, 2026-09-07 — partly, and the part they
resolve is already in the architecture.** The drift curve splits into a LEVEL (its per-swath
mean) and a SHAPE (variation about it), and only one of those is degenerate:

    site        swaths   level spread   median shape   max shape
                         (between,mm)   (sd within)      (p-p)
    elba             4          14.23           9.46       68.90
    whitewater       4          40.20          15.20      106.80
    mnrv             6          34.43          16.01       75.40
    cook             4          41.16          10.37       62.30
    carlton          5          21.36           5.55       56.40
    battlecreek      4          27.94          13.01       70.40

1. **Swath-to-swath differences** are constrained by gen1 overlaps (`align_swaths`). No
   control needed.
2. **Within-swath shape** is real and large — up to 106.80 mm peak-to-peak — and it IS
   identifiable, because nothing else in the pipeline is a function of time. Control marks
   cannot constrain it: Elba has 8 marks on 5 lines, 1.6 per line, and a spline shape
   cannot be fitted from 1.6 points.
3. **The common level** is constrained by none of the above; it is the free-network gauge.
   **Control fixes exactly this**, and that is the existing ground-control datum step.

So marks resolve the degenerate part and already do. What they cannot do is separate
`drift(t)` from a datum gradient ALONG-track: within one swath those are collinear
functions of the same variable. The resolution is that the datum gradient is known
INDEPENDENTLY from the geoid grid — so there is nothing left to separate, provided the
geoid is right. What was absorbed was an ERROR in it, now fixed. The residual exposure is
any UNMODELLED datum gradient (a vendor-introduced tilt, say), which would still be eaten
along-track and would be invisible.

**NEGATIVE RESULT — the drift is not measurably undoing `align_swaths`.** Both stages carry
a free per-swath constant, so they are collinear by construction, and the per-site
correlations looked suggestive (−0.16 to −0.71 at five of six sites). Pooled over 27 swaths
centred within site: **r = −0.2589, p = 0.192, slope −0.1096 ± 0.0818**. Not significant.
The degeneracy is real but is not causing measurable harm.

**Revised recommendation on #63.** The urgent part was the geoid and is fixed. Zero-meaning
the drift is NOT free — that per-swath constant is doing work, and removing it without
compensation shifts the DoD. What is left is a CLAIM problem rather than an engineering
one: the pipeline says `register_gen1` is gen1-internal and `apply_datum` is the only
cross-epoch stage, but the drift fit re-levels each swath against gen2. Say that plainly.

OPEN, Andy's call (task #63): constrain the drift fit so it cannot absorb a datum term —
zero-mean per swath at minimum, possibly zero-slope, or fit it only after an independent
datum. Any of these changes what the pipeline measures. Figures:
`figures/battlecreek_geoid_fix.png`, `figures/drift_absorbs_alongtrack_tilt.png`.

Found while resolving counties for the control-mark work, not by looking for it.

## The geoid and the control are COMPLEMENTARY, not redundant (measured 2026-09-07)

Andy's test: if the geoid difference is only a constant, control marks carry the load and
the geoid term is redundant; if it has a slope, it does what control cannot. Measured — it
has a slope, and the slope is the load-bearing part.

    site          const    tilt E    tilt N   tile km   tilt p-p   field p-p   RMS resid
                     mm     mm/km     mm/km       ExN         mm          mm          mm
    elba          67.19     0.822    -0.463   2.5x3.5       3.71        4.44       0.404
    whitewater    57.88    -0.150     1.620   2.5x3.5       6.05        8.12       0.936
    mnrv          67.25     1.157    -0.151   2.5x3.5       3.46        5.22       0.615
    cook          25.40    -0.500     6.135   2.5x3.6      23.04       22.91       0.687
    carlton       50.81     2.871    -0.566   2.4x3.5       8.92       10.70       1.108
    battlecreek   16.98    -1.750    -1.254   0.6x0.9       2.17        2.23       0.048

**Three scales, and control can only reach one of them.**

*Within a tile.* The tilt is 2.17–23.04 mm peak-to-peak. Elba's control constant — the only
one we have — is **58.70 ± 25.89 mm from 8 marks on 5 lines**, and Elba's tilt is 3.96 mm,
**0.153×** that SE. Marks cannot see it. (At cook the tilt is 23.04 mm, comparable to the
SE — but cook has no marks at all.)

*Between tiles.* elba and whitewater are the SAME survey on the SAME geoid, 12.8 km apart,
and their geoid differences are **+66.87 and +57.18 mm — a 9.70 mm gap**. A control constant
measured at one site does not transfer to its neighbour.

*Statewide*, which is the project goal. GEOID03→GEOID18 over Minnesota ranges
**272.23 mm** (min −79.72, max +192.51, sd 32.25); GEOID09→GEOID18 ranges **139.80 mm**.
No amount of local control supplies that field.

The planar fit is good locally — RMS residual 0.048–1.108 mm — so a plane per tile is
adequate; it will not be across the state.

**This also reframes the geoid bug: the wrong grid got the SHAPE wrong too, not just the
level**, and no control campaign could have caught that:

    site          const err    tilt err, p-p over tile
    battlecreek     54.87 mm            1.96 mm
    cook            26.39               12.35
    carlton         27.75                7.44

**The division of labour.** The geoid supplies the SHAPE — gradient and spatial variation,
deterministic, everywhere, no marks needed. Control supplies the LEVEL — the one constant
per site the geoid cannot give, because that constant is the survey's own vendor error, not
geodesy. Neither substitutes for the other, which is why both are in the pipeline.
Figure: `figures/geoid_slope_vs_control.png`.

## The ground-control datum — measured 2026-09-07, and it is NOT appliable at four sites

README step 6 and `ground_control/FRAME.md` call the datum a **required** pipeline step,
because the gauge choice it removes is worth 42.40 mm at elbaext against a +2.12 mm
correction. Verified on the real product, `gauge_invariance_residual` over all six
candidate gauges: **uncorrected spread 42.40 mm, corrected 3.553e-15 mm.** The machinery
works.

But `absolute_datum_mm` is **None at all seven tiles**, and it cannot simply be filled in.
The bundled 2008 control is eight SE-MN counties (dodge, fillmore, houston, mower, olmsted,
steele, wabasha, winona; 1004 marks). Nearest mark to each site centre:

    elba          1.5 km    35 within 10 km    CAN measure
    whitewater    0.7 km    23                 CAN measure
    mnrv         50.4 km     0                 NO 2008 control bundled
    battlecreek  70.6 km     0                 NO 2008 control bundled
    carlton     243.0 km     0                 NO 2008 control bundled
    cook        430.8 km     0                 NO 2008 control bundled

**What the gauge is worth at each site** (spread of per-swath dz: how far the WHOLE DoD
moves if a different line is pinned. `dod = dod + (g2 - g1)/1000`, a uniform shift):

    cook           24.90 mm   (4 swaths)
    battlecreek    37.50      (5)
    elba           38.60      (4)
    elbaext        42.40      (6)
    carlton        57.40      (5)
    whitewater     90.40      (4)
    mnrv          287.70      (7)   <- 7.5x Elba

**Elba is near the BEST case, not a typical one.** The +2.12 mm correction that made this
look negligible was measured at the site with almost the smallest lever, and the FRAME
already called line 133 "a lucky pin". The sites with the largest levers -- mnrv, whitewater,
carlton -- include two with no control at all. That is the worst pairing available.

So "required" is achievable at two sites of six. The other four need control transcribed
from their own acquisitions' validation reports — task #11, a data problem, not wiring.
Either do #11 or soften the claim.

And the one constant that exists, `SITE_DATUM_elbaext.json`, is gauged on line **133**,
which **elba does not carry** (its psids are 135–138). `on_zero_line()` is arithmetic
WITHIN one product; it does not carry a constant from elbaext's surface onto elba's
independently-solved one, and ties are known to be extent-dependent. elba and whitewater
each need their own `run_site_datum.py` run.

**Stale-number note.** The dz spread quoted throughout the repo was **44.60 mm**, from an elbaext build older than 2026-09-01. It was corrected to **42.40 mm** on 2026-09-07 in README, this file, `SWATH_ALIGNMENT_METHOD.md`, `apply_datum.py`, `pipeline.py`, `ground_control/FRAME.md` and the handoff, together with the per-swath table each quotes. It will move again: elbaext was last built 2026-09-01, BEFORE the 1/variance weighting of 09-05, and rebuilding it also invalidates `SITE_DATUM_elbaext.json`, measured against those corrections three minutes after they were written.

## ⚠ A CACHE SERVED A SUPERSEDED METHOD (2026-09-10)

The swath-constants cache was keyed `tile|res|tie|exclude` — **nothing about the METHOD** —
so it CAN serve constants from a superseded alignment. Commit 9e78f4a (09-05) reweighted
the swath network by 1/variance instead of cell count.

**CORRECTED 2026-09-10, same day: it did not actually happen.** I first wrote that a run
"mixed two alignments in one answer". It did not. Promoting `our_surface` had moved the
cache path into the package, where no cache existed, so that run computed every tile fresh
and bypassed the stale committed cache entirely — by accident, not design. Proof: a
deliberate clean re-solve, cache deleted and keyed on the code, reproduces the earlier run
on **29 of 29 marks** and gives the same mean, −5.5231.

So the mixing was a HAZARD I found, not an event I observed. What the numbers actually
show is simpler: **−4.04 was computed under the OLD weighting and −5.52 is the same 29
marks under the current one.** Four marks move, one by 45.1 mm, and the median stays +2.79
because the movers sit in the tails — which is how a change like this hides.

FIXED (b7bb21c): the key now carries a digest of the whole `coreg` module, and the cache
moved out of the package's bundled-control-data directory — promoting `our_surface` had
moved `_HERE`, so a runtime file was being written beside the control CSVs. Two regression
tests.

**`bridge_mm = −4.04` is the MEAN of 29 marks and rests on the superseded weighting.** It
feeds elba's adopted +58.70 ± 25.89 and was applied to whitewater on 09-10. Re-measuring
clean is in flight; nothing is re-adopted until it lands.

## ITEM 2 SO FAR — the datum is thinly sampled at these tile sizes

    site        c1 delivered      marks lines   geoid    DoD shift   zero_line
    elbaext     +62.74 ± 23.38      8     5    +67.44      +2.18        133
    elba        +44.41 ± 26.32      4     3    +67.28     +20.35        135
    whitewater   −6.52 ± 24.26      5     2    +57.92     +61.92        143

RECORDED, NOT ADOPTED. The spread tracks LINE COVERAGE, not geography: 5 lines +2.18,
3 lines +20.35, 2 lines +61.92 mm. Every run warned — whitewater's says "only 2 flight
line(s) carry a mark: the SE over lines is barely defined" and names 144 and 146 as
carrying none. Applying whitewater's +61.92 on two lines would be worse than applying
nothing.

`control_reach` (new step, 36bff5f) now reports this before anyone runs: it is the gen1
mirror of `completeness`, and it exists because mnrv passed every other step while reaching
ZERO marks. mnrv's 19 tiles are now fetched (618 MB, 19/19 reachable).

## FLOODPLAIN VEGETATION — settled 2026-09-16, LoD not correction

**The primary signal in this DoD is a VEGETATION-STATE DIFFERENCE BETWEEN EPOCHS**, on
hillslopes and floodplain alike (Andy, 2026-09-16). gen1 flew **November 2008** with
floodplain sedge standing as tall dead biomass; gen2 flew **May 2021** with that flattened
by winter and regrowth short. Not a sensor difference — a season difference, and the dates
are fixed for the whole 2008 MN program against 3DEP, so the sign is predictable statewide.

**THE CONTROL, and the result worth keeping.** On flat floodplain (slope ≤ 2°, 15,961 cells
of 5 m), DoD against each epoch's own near-ground spread `p90 − p10`:

    spread source                 ret/cell  spread p50      m    b (mm)      r
    gen1 CSF ground                     17      220 mm  -0.294     +55.0  -0.353
      ^^ the intercept is extrapolated (spread never reaches 0) and not a measured quantity.
      Do not quote it. Only the SLOPE is used.
    gen2 class-2 GROUND only           158       80 mm  -0.101      +1.7  +0.012
    gen2 ALL near-ground (Hg+Hn)       251      120 mm  +0.006     -11.4  +0.070

gen2 shows **nothing**, even with vegetation included and 15× the returns. Terrain
roughness or gridding would appear in both epochs; this appears only in gen1. Andy caught
the first version of this test, which compared gen1's imperfectly-classified ground against
gen2's successfully-classified ground — not like-for-like. The null survives the fix.

**WHY LoD AND NOT A CORRECTION.** `DoD + 0.294 × spread` is fitted and available, but:

    penetration OK    48,574 cells  spread 207 mm  DoD p50   +1.3 -> +62.0
    penetration FAILS  5,097 cells  spread 259 mm  DoD p50 -216.2 -> -140.2

It works where the bias is small and fails where it is large — 35% recovery on the failure
cells, because when no pulse reaches ground the distribution NARROWS onto the vegetation
top and the covariate stops responding. And on the 90% that do penetrate, applying it turns
"no detectable change" into "60 mm of deposition", resting wholly on an intercept that is
extrapolated past the data (see above).

**SHIPPED:** `vegetation_lod.inflate_lod` adds the *un-applied* correction in quadrature on
flat floodplain. `m` and `slope_max_deg` and the failure floor all REFUSE without a value.
Banks keep their LoD, which is the point.

**★★ FIXED 2026-09-20 -- floodplain support cells.** The surface now keeps the
well-measured floodplain cells in its stable set (`correction_surface_keep`; builder flag
`--floodplain-support`): |laplacian| <= curv_max, gen1 reaches gen2, and NO gen1 return in
the 0.15-2.00 m tussock band. 30,845 cells at elbaext, closing the support gap 578 -> 61 m.

COMPARED STRUCTURALLY ONLY -- every field de-gauged by one constant on stable uplands,
because a vertical offset is GAUGE, not error (Andy 2026-09-20):

    five patches vs the surface-free route   delong 70.0 NMAD / 257.3 RMS
                                        +fp support 24.4 /  31.4     (8.2x)
    whole flat floodplain                    delong 66.2 /  93.7
                                        +fp support 32.7 /  50.1     (2.0x)
    floodplain own spread     delong 127.4 | +fp support 69.5 | independent 87.2
    uplands vs independent    delong  37.2 | +fp support 38.2  (unchanged)

Two-fold structural improvement, uplands untouched, 75.9% of upland cells bit-identical.

**★★ THE FLOODPLAIN "EROSION" AT ELBA IS THE DELONG CORRECTION SURFACE (2026-09-17).**
Andy: "the floodplain upstream of the pointed meander bend is still showing significant
erosion." It is an artifact of the ADOPTED route. The 5 largest patches that survive the
vegetation LoD (1,850 cells, 44.098 N -92.009 W):

    applied DeLong correction surface   +261 mm   (rest of floodplain +84, whole grid +94)
    DoD, delong route                   -181 mm
    DoD, INDEPENDENT route                -0 mm   <-- no erosion there at all

Regressing (DoD_delong - DoD_independent) on the applied surface over 35,753 flat floodplain
cells: slope -1.0143, r -0.9932. The route difference IS the surface.

MECHANISM: the IDW surface is fitted on STABLE cells and the valley cut removes 95,506
floodplain cells from that set, so over the valley floor it EXTRAPOLATES. Fraction of flat
floodplain still flagged, by distance to the nearest stable cell: 0.9% (0-25 m), 2.5%
(25-50), 5.3% (50-100), 5.6% (100-200), 17.4% (200+, median DoD -102 mm).

It was invisible because the step applied the grid and discarded it, and the beam table
carries every dz_* term EXCEPT this one. Fixed in 985a24d -- saved as
applied_correction_surface.npy.

⚠️ DeLong was adopted on STABLE-UPLAND skill, all measured where the surface is constrained.
Nothing in that comparison tested the floodplain. DO NOT read floodplain change off the
DeLong route until this is resolved.

**⚠️ SUPERSEDED 2026-09-16 -- `m` NO LONGER EXISTS.** The vegetation term is gen1's own
p90-p10 span, used directly: `lod_veg = hypot(lod, span)`. The fitted slope was the
CONDITIONAL MEAN of DoD given spread, so half the cells exceeded their own term by
construction and the vegetated patches kept flagging. A fitted mean is a BIAS estimate; an
LoD needs an UNCERTAINTY. Thickest-vegetation decile still flagged: 77.8% bare, 50.2% with
the fitted m, 26.7% half-span, 24.0% p50-min, **4.9% with the full span**. Result: delong
LoD 50 -> 220 mm (153,232 -> 147,283 detected), independent 78 -> 237 mm (122,935 ->
117,538). The population discussion below is kept because it explains how the fit went
wrong, not because the fit is still used.

**⚠️ THE 450 mm FAILURE FLOOR IS ARBITRARY AND FLAGGED FOR REPLACEMENT** (Andy 2026-09-16).
It is CIRCULAR: the value that silences 90% of detections in the DoD it then judges. The
span alone cannot replace it -- on failure cells the span is 1.35x the penetrating cells'
while the error is 3.1x, leaving 38% flagged. Needs an instrument external to the DoD.

**[HISTORICAL] THE POPULATION FOR `m` WAS EVERY MEASURABLE FLAT FLOODPLAIN CELL** -- Andy 2026-09-16, "just
use the floodplain cells and those with vegetation to decide on the LoD". It is fitted over
all cells with >=10 gen1 returns, 53,671 of elbaext's 58,767 flat (<=2 deg) floodplain cells.
The earlier m = -0.294 was the SAME fit restricted to the 16,092 cells (30%) that
`nearground_cells_sn.npz` happens to cover -- an artifact built 2026-08-26 for the
return-structure work. A file footprint, not a population, and it undersized the LoD by ~60%.
The difference is WHICH CELLS, not the estimator: on the 16,050 cells both artifacts cover,
raw-quantile and 2 cm-binned spread correlate +0.942 and give m = -0.297 vs -0.251.

    tile                  m (all)   m (cube 30%)   LoD med.   detections
    elbaext_independent   -0.475    -0.325         78->134 mm 122,935 -> 119,885
    elbaext_delong        -0.494    -0.329         NOT RUN -- no stable.npy

**elbaext_delong CANNOT BE EVALUATED YET.** The adopted route's tile was built 2026-09-11 and
the `stable.npy` save was only added 2026-09-16 (c87c069), so `detect_change_standard` had no
stable set. It needs a route rebuild. Until then every detection count on the DeLong route is
from the independent tile, not the adopted one.
(`figures_with_vegetation_lod.py` now REFUSES a missing stable.npy rather than substituting
an empty mask -- the tell was 3,020 cells GAINED after the LoD was WIDENED.)

**IDENTIFYING THE UNRECOVERABLE CELLS** is a yes/no physical test with no threshold: does
any gen1 return reach gen2 (`min(d_mm) > 0`)? 9,460 floodplain cells fail it, in 55 coherent
clusters ≥20 cells holding 57% of them. Pooling proves they are unrecoverable — 0.0% of
11,091 shots reached ground in the largest cluster.

**NOT A QUANTITY WE HAVE.** At low gen1 spread the DoD sits a few tens of mm above zero
(binned medians +32.2 / +25.6 / +23.4). That residual is UNEXPLAINED -- deposition, residual
level error, and something else are all consistent with it, and nothing here separates them.
It is also INSIDE the flat floodplain, which the inflated LoD masks, so no product reads it.
Do not name it, do not correct for it, do not quote a value. (Andy, 2026-09-16: "we do not
know if it is deposition or datum or something else. And it will be masked in any case.")

**STILL OPEN:** `m` is population-sensitive: -0.294 on cube cells,
-0.378 on beam-table cells, so the SHIPPED inflation is probably too small; the slope cut is a patch over an elevation-cut floodplain that includes
banks; and metre-scale outliers (1.3% of cells, p10 −2010 mm in the top spread bin) are
structures, not vegetation, and hijack any mean-based statistic.

## THE ORDER (Andy, 2026-09-11, DeLong) — work down this list

Supersedes the 2026-09-10 order below, which is kept for the record. What changed: Andy
adopted DeLong's correction surface as the pipeline default (`correction_surface=True`,
`along_track_drift=False`, commit 744c9a1) after measurement, and directed that the
marks/corrections subsystem move out of the main pipeline (#68).

**What DISSOLVED, and why.** Old item 2 — "measure and APPLY the datum at all six sites" —
is gone as a pipeline requirement. The correction surface is fit on the stable residual, so
it pulls gen1 onto gen2's frame whatever datum was applied first: a 54.87 mm geoid error
moves the DoD by −54.870 mm with the surface OFF and by **0.000 mm** with it ON
(`test_the_correction_surface_absorbs_a_wrong_geoid_and_hides_it`). The product is
therefore absolute by inheritance, at gen2's level, with no datum step to apply. The marks
now BOUND that level instead of supplying it: gen2 open-ground at Elba is
`constant_mm` −6.5647, `sd_field_mm` 31.10, over 139 marks. Old item 7 (rebuild elbaext) is
subsumed by new item 1.

**The new standing hazard.** The surface absorbs everything smooth — the geoid error above,
the injected scanner roll, any datum mistake. That is why it wins on skill and why it is
robust, and it is also why every diagnostic that works by SEEING a smooth error is now
blind. The battlecreek bug would not have shown up in a product. Items 2, 3 and 6 exist
because of this.

1. **Rebuild all six sites on the new default.** Every `corrections.json` on disk was
   written with drift ON and the surface OFF, so no shipped product matches the code that
   would rebuild it. Largest single consequence of the switch; everything downstream
   depends on it. One site at a time — shared laptop.
2. **Report the stable set's COVER COMPOSITION as a pipeline output.** The inherited level
   is gen2's, and gen2 floats high under canopy: forest stable ground sits −28.61 mm
   against open at Elba. Elba is safe by terrain luck — the slope ≤ 3° cut leaves 0.1%
   forest (1,035 of 1,072,025 returns), 69.7% open — but that is not construction. Where
   flat ground IS forested (cook, arrowhead) the surface would pull gen1 onto a
   canopy-biased level. Must be reported per site, never assumed. NEW.
3. **Verify the correction surface at a SECOND tile (elbaext).** Every number behind this
   decision is one tile and four lines: skill 0.444/0.384/0.283 vs the drift's
   0.160/0.142/0.087, and the drift's per-line damage 5.07 → 19.91 mm. If those do not
   reproduce, the decision needs revisiting — so this gates any claim that it generalises.
4. **Push.** ~107 commits.
5. **Fix the documentation the switch FALSIFIED.** README step 6 still calls
   `absolute_datum` required; it is now inherited, not applied. #63's complaint (the drift
   re-levels each swath against gen2, contradicting "gen1-internal") is moot for the
   default path but is MORE true of the surface, which re-levels everything against gen2.
   Re-lead the framing; do not append a note.
6. **Make the independent checks ROUTINE (#68).** The marks and the geoid are now the only
   instruments that can see what the pipeline has stopped being sensitive to. They must run
   as a standing check rather than ad hoc, or a battlecreek-class error ships silently.
   This is the real content of "placed aside and out of the main pipeline": they stop being
   stages and become the audit.
7. **#65 — control set from the Site's own SURVEY, not a default.** Was a prerequisite for
   applying the datum; now a prerequisite for the CHECK in item 6 being correct. Same fix
   (stem in `acquisitions`, resolved from `Site.gen1_project`, refusing when absent), lower
   urgency, still a silent-wrong-answer bug at 18 bare call sites.
8. **The cleanup queue, unchanged by DeLong:** #55 co-location, #56 the 127 unreached
   scripts, #50 Battle Creek's unexplained 1-cell instability, and Battle Creek's valley
   top (Andy must state an elevation or a fraction ceiling).

The vegetation thread (#25, #19, #16) stays parked, but its RATIONALE changed: it is no
longer about correcting gen1's ground, it is about bounding what canopy costs the inherited
level. That is item 2's question.

## SUPERSEDED — the order of 2026-09-10 (kept for the record)

1. **Control set from the Site's SURVEY, not a default.** `gen1_datum.load_control()`
   resolves `DEFAULT_CONTROL = "mn_dnr_2008_control_semn"` at **18 bare call sites**, and
   `residual_field.GEN1_CSV` hard-codes the same file. Five control sets sit in
   `groundtruth/data/`, so the marks added 2026-09-09 for mnrv, battlecreek, cook and
   carlton are reachable but never reached. Same shape as the geoid default that cost
   +54.87 mm at Ramsey. Fix with the pattern that worked: the stem goes in
   `acquisitions`, resolved from `Site.gen1_project`, refusing when absent. **Task #65 —
   prerequisite for 2.**
2. **Measure and apply the datum at all six sites (#60, and #11 in practice).**
   `absolute_datum_mm` is None at 6 of 6 while README step 6 calls it required. Needs
   Andy's confirmation of `run_site_datum.py`'s eleven parameters. Answers: are the six
   DoDs mutually comparable once each sits on its own surveyed level?
3. **Push.**
4. **#63 as a DOCUMENTATION fix**, not an engineering change: `apply_datum`'s drift stage
   re-levels each swath against gen2, which contradicts the "gen1-internal" claim.
5. **Co-locate each findings document with its investigation (#55)** — ~10 moves.
6. **The 127 unreached scripts (#56)** — Andy's call. By function and by removal-test: a
   reachability scan pointed at load-bearing code four times in one week.
7. **Rebuild elbaext** — the only product still on pre-variance-weighting code. Low value
   until 2, which may supersede `SITE_DATUM_elbaext.json` anyway.
8. **The vegetation thread (#25, #19, #16)** — parked; the correction measured worse than
   the median on open ground and nothing since has changed that.

## The queue

`analysis/NEXT_SESSION.md` holds the prepared plan: the last 7 moves (ready, needs a
one-line fix in three tests first), the decisions waiting, the open science, and what NOT
to redo.

## Open, and needing Andy

0. **Five structural decisions are queued** (tasks #54, #56, #57, #58, and the "alongside"
   question). They block further promotion, so nothing more should move until they are made:
   * **Two ground-control trees.** `src/lidar_diff_icp/groundtruth/` (7 modules) and
     top-level `ground_control/` (29 scripts) BOTH have a `datum.py` -- 230 vs 236 lines,
     **not identical** -- and both have `data/`. Duplication in load-bearing datum code
     decides which answer you get.
   * **`trust/` is a package at the repo root, not under `src/`,** so it imports only
     because everything runs with the repo root as CWD. 43 files import it. This blocks
     moving `crossline_fit.py` + `swath_across_track_test.py`, a coupled pair.
   * **The 127 unreached scripts.** Andy 2026-09-06: delete if we want; investigations were
     exploratory and can go if not needed. Both preservation conditions are ALREADY met --
     `rivernetworkx` and `catchment-dod-balance` are declared git dependencies, installed
     and importable. The value is the findings documents, not the scripts; git holds the
     scripts.
   * **`ground_control/FRAME.md` is a SECOND unmarked "READ FIRST" frame** (2026-08-27),
     missed by the consolidation that marked six others superseded.
   * **Does "alongside" mean outside the package?** The 9 vegetation-correction modules and
     9 tools are all reachable. "Alongside" was about DEPENDENCY DIRECTION, not file
     location, so they could live in the package as a subpackage the pipeline optionally
     reaches over to -- but that is Andy's architecture call, not a cleanup.

## Also open, and needing Andy

1. **Battle Creek's valley top.** Its histogram cut removes 72.4% of the grid (283.6 m) in
   a built environment where graded lots set the modal elevation. Its `stable_sigma` is
   computed on that set. He must state an elevation or a fraction ceiling.
2. **A 1-cell instability at Battle Creek, cause NOT established.** Two runs gave 4,841 vs
   4,842 stable cells at the same valley top with identical `z_after`, moving the DoD by
   <= 0.883 mm through the drift fit. Re-running twice more is byte-identical, so the
   pipeline IS deterministic given identical inputs — something in the input differed and I
   did not find it. Sub-mm, but unexplained.
3. **The deletion candidates** — `analysis/DELETION_CANDIDATES.md`, 26 scripts, nothing
   deleted. The eleven dated 2026-09-01 or later are FALSE POSITIVES and are named there.
4. **105 unpushed commits.**

## Guardrails

* **PLOT WHAT YOU MAKE**, same turn, unprompted.
* **A claim that would change course must be verified in the same turn it is made.**
* **When Andy points at prior work, LIST THE FILES.** A repeated request is evidence the
  premise is wrong.
* **No invented thresholds.** Say it in one sentence and let him decide.
* **A regression test must be shown to FAIL without the fix.** Two vacuous tests were caught
  this way on 2026-09-05; both breaks passed before the scene was fixed.
* **Numbers in prose are pasted from command output.** Four claims went out unverified on
  2026-09-05-06 — mnrv's 0.67, two correlations, and a `swath_tie` defect that was my own
  code. All four were one command away.
* **Evidence of use hides from static scans.** Three times: scripts that COMPUTE a quantity
  rather than loading it; documents naming a file without its path; work recorded only in a
  commit message. Inventory by function, import AND artifact, and read the result.
* Shared laptop: one heavy job at a time, watch `free -h`.
* Commit granularly; never push/tag/release without explicit current-message authorization.
