STATUS: COMPLETE (2026-11-18; 21 sources read in full or at abstract level, 5 data sources, 48 log entries)

# Plant phenology and remote-sensing identification literature for the Whitewater River floodplain (Elba, MN; 44.12 N, 92.00 W)

Compiled 2026-11-18 by a Claude subagent for Andy Wickert. Purpose: tell reed canary grass
(Phalaris arundinacea) and other floodplain and bluff cover types apart from Sentinel-2 (10 m),
Landsat (30 m) and NAIP (0.6 m) imagery by phenology, because the 2008 (leaf-off, November) and
2021 (1 May, green-up) lidar surveys have vegetation-dependent vertical bias.

Standards: every number is pasted from the source read in this session; "abstract only" is
stated where the full text was not read; UNVERIFIED marks a number taken from a secondary
source. Latitude is given with each phenology observation because dates shift with latitude.

## (a) Summary table

See `phenology_literature_table.csv` in this folder (one row per taxon x region x source).
A condensed reading of it is in section (f) at the end.

## (b) Per-source notes

### S01. Bernthal & Willis 2004 (Wisconsin DNR, Landsat-7 pilot) -- FULL TEXT READ (pdftotext)
Bernthal, T.W., and K.G. Willis. 2004. Using Landsat 7 Imagery to Map Invasive Reed Canary Grass (Phalaris arundinacea): A Landscape Level Wetland Monitoring Methodology. Final Report to U.S. EPA Region V, Wetland grant CD975115-01-0. Wisconsin Department of Natural Resources, Madison, March 2004. https://dnr.wisconsin.gov/sites/default/files/topic/Wetlands/rcg_mapping_report.pdf
One Landsat-7 ETM+ scene (Path 24 Row 30, ~43 N, south-central Wisconsin), 30 m, single date. Executive summary: "Its widespread nature and unique spectral signature in fall due to its comparatively later senescence make it feasible to map reed canary grass using spectral classification of satellite imagery." Date choice (sec. 4.1.1): "From many images taken throughout the previous growing season, spring through fall of 1999, we selected data captured on October 15. Imagery from this date was cloud-free, and displayed the highest potential for our objectives." Classes: heavy dominance >=80 % RCG cover; co-dominance 50-79 %; absent/subdominant <50 %. Method: "a combination of unsupervised and supervised classification procedures"; MMU 0.5 acre. Accuracy (sec. 4.2.2, Table 1, 249 plots): "overall accuracy of 71 percent: (89 + 15 + 74) / 249 = 0.71"; user's accuracies "86%, 41%, and 69%" for heavily dominant, co-dominant, absent/subdominant. Error matrix: of 104 pixels classed heavy-dominant, 89 verified heavy, 8 co-dominant, 7 absent/sub; of 37 classed co-dominant, 18 were actually heavily dominant. "Classification of forested, shrub, and open water wetlands was less than satisfactory, due to the spectral influences of woody vegetation and open water." Caveat (4.2.3): "The procedures here distinguish only a single species, reed canary grass, from all other plants ... purple loosestrife (Lythrum salicaria), giant reed grass (Phragmites australis) and hybrid cattails (Typha x glauca) ... will be considered simply 'other vegetation.'" Results: 737,454 acres of wetland analyzed; 79,490 acres (11 %) heavily dominated; 23,378 acres (3 %) co-dominated.
Meaning for Elba: the method that works at 30 m is ONE mid-October image, when RCG is still green and other herbaceous wetland cover has senesced. No spring/multi-date NDVI was used. It never distinguished RCG from Typha/Phragmites, and pasture or hay fields were not in the wetland mask (upland masked out), so the RCG-vs-mown-grass confusion that matters on a hayed floodplain field was NOT tested here.

### S02. Hatch & Bernthal 2008 (Wisconsin DNR, statewide) -- FULL TEXT READ (pdftotext)
Hatch, B.K., and T.W. Bernthal. 2008. Mapping Wisconsin Wetlands Dominated by Reed Canary Grass, Phalaris arundinacea L.: A Landscape Level Assessment. Final Report to U.S. EPA Region V, Wetland grant CD96544501-0. Wisconsin Department of Natural Resources PUB-WT-900-2008, October 2008. https://dnr.wisconsin.gov/sites/default/files/topic/Wetlands/RCGFinalReport10_08.pdf
Phenology statement (sec. 2): "Reed canary grass is one of the first plants to sprout in spring, and remains green well into the fall ... It is easily classified by satellite imagery because of its late senescence period, growing long past other wetland species." Imagery (sec. 4.1.1): twelve Landsat-7 ETM+ and Landsat-5 TM scenes, bands 1,2,3,4,5,7 (no thermal); "As determined by Bernthal and Willis (2004), mid-October is the best for distinguishing reed canary grass from other vegetation because of its late senescence. They concluded that multiple dates and use of Normalized Difference Vegetation Index (NDVI) are not needed to classify reed canary grass in Wisconsin. The scenes acquired were selected using three criteria: 1) as close to mid-October as possible, 2) 10% or less of cloud cover, and 3) within a small range of years, none being older than 1999. The range of dates for the imagery acquired was Oct 09-Oct 31, 1999-2003." Classes collapsed to two: A <50 % RCG, B >=50 % RCG (because the 2004 co-dominant class had 41 % user's accuracy). Confusion noted at class definition: "One exception to this is some wetter sites where reed canary grass and cattails (often a mix of Typha latifolia, T. angustifolia and Typha X glauca) are both near 50% cover." Accuracy (sec. 4.2.2, Tables 1-7, per scene): overall 83, 83, 66, 82, 78, 66, 61 %; RCG user's 89, 75, 72, 92, 82, 84, 84 %; RCG producer's 38, 75, 53, 79, 83, 44, 50 %. Executive summary: "a range of 61-83% overall accuracy and 72-92% user's accuracy." Errors of commission (sec. 4.2.2): "It appears that heavily grazed pasture fields were classified as reed canary grass. They are noticeable on the grid by their distinct geometric shapes. In some cases, a pasture might be reed canary grass, but it was hard to tell during the site visits because of how short the grass blades were cropped. Another point of confusion was cranberry bogs." Woody wetlands: "some deciduous shrubs that are not as woody could have lost their leaves by the second half of October at the time of the Landsat scene capture" so RCG under shrubs was detected. Results (Tables 8-10): 498,250 acres RCG-dominated of 5,065,419 classified wetland acres (9.84 %); emergent wetlands 305,878 acres, 26.64 % of class; forested wetlands 82,218 acres (3.36 %); lowland shrub 109,976 acres (7.40 %). (The 510,000 acre / 26.3 % figures quoted in the search snippets are from a Hatch presentation, not this report; the report says 498,250 and 26.64 %.)
Meaning for Elba: producer's accuracy for RCG was as low as 38-53 % in several scenes, i.e. the mid-October single-date method misses a lot of RCG; grazed PASTURE is the named false positive, which is exactly the hayed-field case on the Whitewater floodplain. The method relies on an a priori wetland mask (WISCLAND) and does not try to separate RCG from upland cool-season grasses.

### S03. Waggy 2010, FEIS species review of Phalaris arundinacea -- READ (web, in three parts)
Waggy, M.A. 2010. Phalaris arundinacea. In: Fire Effects Information System [Online]. USDA Forest Service, Rocky Mountain Research Station, Fire Sciences Laboratory. https://www.fs.usda.gov/database/feis/plants/graminoid/phaaru/all.html
Secondary (review) source; it compiles many primary observations with region labels. Seasonal development, as summarized by the fetch tool from the page (quotes are the page's wording where marked): growth begins "early spring, around April in North America"; "Growth peaks in mid-June and declines by mid-August in the Pacific Northwest, Illinois, and Wisconsin" (refs 148, 133, 157); flowering, Wisconsin: "May to mid-June" (refs 119 Hoffman & Kearns 1997 WDNR fact sheet; 128 Howe 1999); Iowa: "spikelets emerged May 2" (Holt 1954); North Dakota: "around June 20"; Illinois May-July; Great Plains May-June. Seed: "Wisconsin: seed matures and disperses in late June" (Klopatek & Stearns 1978); north-central states: "seed heads appear about the first week of June, ripen the last week of June, and mature about the first week of July" (Wheeler & Hill 1957). Senescence: "After seed dispersal, heads and stems die back to the upper leaves, but the rest of the plant stays green" (Wheeler & Hill 1957); "In Wisconsin, shoots collapse in mid to late summer into a dense mat" (Hoffman & Kearns 1997); "Growth may resume in fall" (refs 157, 119). Fuels: "it remains green into the fall" (refs 119, 277 Tu 2004); "As mature reed canarygrass plants senesce late in the season, litter accumulates and forms thick, impenetrable mats" (refs 119, 294 Weaver 1960); "litter comprised more than 80% of total biomass in unburned plots in May" (ref 129). Height: "2 to 7 feet (0.6-2 m) tall"; Missouri floodplain (Weaver 1960): "shoots typically reach about one foot in April and up to 3.5 feet (1.1 m) by July". No winter standing-dead description.
Meaning for Elba: the FEIS record supports a spring green-up before other wetland herbs, flowering late May to mid-June at ~43-45 N, mid-summer LODGING of the canopy into a mat (so height is not monotone through the season), fall regrowth, and a litter mat in spring that is >80 % of biomass -- that litter mat is what a May lidar sees.


### S04. Minnesota DNR, reed canary grass species page -- READ (web)
Minnesota Department of Natural Resources. Reed canary grass (Phalaris arundinacea). Invasive terrestrial plants, grasses. https://www.dnr.state.mn.us/invasives/terrestrialplants/grasses/reedcanarygrass.html (no last-updated date shown; footer "(c) 2026 Minnesota DNR").
Verbatim: "Reed canary grass is a coarse, cool season grass that grows 2 - 6 feet tall." "Reed canary grass is one of the first grasses to sprout in the spring." "Reed canary grass blooms from May to mid-June." "Seeds ripen in late June." "Flowers ... typically fade to beige in late summer or fall." "If cut during the growing season, reed canary grass will have a second growth spurt in the fall, so multiple mowings per year are necessary." Recommended "mowing in mid-June and October to reduce seed production." Leaves "up to 10 inches long and 1/4 inch-1/3 inch wide." Look-alikes named: native reedgrasses (Calamagrostis; RCG "greater than 6 feet vs. 4 feet") and orchard grass (Dactylis). No statement about fall die-back date.
Meaning for Elba: statewide-general (no latitude given); confirms the May to mid-June flowering and the fall regrowth after a summer cut, which matters because the Whitewater field is hayed in some years and a hayed RCG stand will re-green in September-October while a hayed brome/bluegrass field also regreens, so autumn greenness alone cannot separate them if both are cut.

### S05. Resasco et al. 2007, late-fall Landsat detects an understory shrub by leaf retention -- FULL TEXT READ (pdftotext)
Resasco, J., A.N. Hale, M.C. Henry, and D.L. Gorchov. 2007. Detecting an invasive shrub in a deciduous forest understory using late-fall Landsat sensor imagery. International Journal of Remote Sensing 28(16): 3739-3745. doi:10.1080/01431160701373721. PDF: https://www.colorado.edu/lab/resasco/sites/default/files/attached-files/resasco_et_al_2007_ijrs.pdf
SW Ohio (Oxford, ~39.5 N), Amur honeysuckle (Lonicera maackii), 8 stands, 15 Landsat 5/7 scenes 1999-2006. "Because L. maackii expands its leaves earlier and retains them later than native trees and shrubs (Trisel 1997), it was predicted that it would be distinguishable in early spring and late autumn images when the tree canopy is bare." Result (Table 2): spring images (March) "indicated no difference in green biomass levels between classes"; April and May images "showed higher SARVI2 values (more green biomass) in low/zero L. maackii cover stands ... attributed to greater cover of spring ephemeral herbs"; "October images resulted in no significant differences (06 Oct 2003) or differences in the 'wrong direction' (24 Oct 2004). This trend suggests that October images are too early in the leaf-off season to detect L. maackii, but that early to mid-November is ideal for this purpose. Many tree leaves have senesced in October, but do not drop for a few weeks, thereby obscuring the satellite's 'view' of the understory." November dates with significantly higher SARVI2 in high-cover stands: 04 Nov 1999 (JM 1230), 09 Nov 2001 (1056), 12 Nov 2005 (1177), 22 Nov 2000 (901); 11 Dec 2001 (828) also high. "the best date to detect L. maackii will also vary each year due to annual differences in temperature and precipitation." No classification accuracy reported (separability only).
Meaning for Elba: a template for the understory-shrub problem on the bluff (honeysuckle, buckthorn): the window is AFTER canopy leaf drop and BEFORE shrub leaf drop, early-to-mid November at 39.5 N; at 44 N it will be earlier (UNVERIFIED by how much). Spring is confounded by ephemerals. Note that the 2008 lidar (November) falls in exactly this window.

### S06. Becker, Zmijewski & Crail 2013, Landsat time-series buckthorn mapping, Oak Openings -- ABSTRACT ONLY (via library catalog snippet; publisher page redirects to login)
Becker, R.H., K.A. Zmijewski, and T. Crail. 2013. Seeing the forest for the invasives: mapping buckthorn in the Oak Openings. Biological Invasions 15(2): 315-326. doi:10.1007/s10530-012-0288-8.
From the catalog abstract (UNVERIFIED wording): 49 Landsat TM/ETM+ images 2001-2011 over the ~1,500 km2 Oak Openings region, NW Ohio (~41.6 N); time-series tasseled-cap greenness used to separate buckthorn thicket, "which leafs out early and senesces late", from other vegetation; supervised classification; thicket presence/absence kappa 0.73, overall accuracy 88 % (60 samples); thicket ~940 ha (0.43 % of pixels) 2007-2011 vs ~690 ha (0.31 %) 2001-2006. Species naming in the snippet (common buckthorn = Frangula alnus) is garbled; the paper covers Rhamnus cathartica and Frangula alnus.
Meaning for Elba: the only published Landsat buckthorn map found; it needed a decade-long image stack, not one date, and classified dense THICKET, not scattered understory shrubs.

### S07. Harrington, Brown & Reich 1989, exotic vs native shrub leaf phenology, southern Wisconsin -- ABSTRACT ONLY (UMN Experts page; Springer full text behind login)
Harrington, R.A., B.J. Brown, and P.B. Reich. 1989. Ecophysiology of exotic and native shrubs in Southern Wisconsin. I. Relationship of leaf characteristics, resource availability, and phenology to seasonal patterns of carbon gain. Oecologia 80(3): 356-367. doi:10.1007/BF00379037.
Southern Wisconsin (~43 N). Rhamnus cathartica and Lonicera x bella (exotic) vs Cornus racemosa and Prunus serotina (native), open and understory. The abstract gives no calendar dates; it reports that longer leaf life let shrubs use understory light before canopy closure and that "from 27% to 35% of the annual carbon gain" of P. serotina, R. cathartica and L. x bella occurred before C. racemosa (the shortest-lived leaves) had leafed out; early emergence mattered more than late senescence for carbon gain. Secondary citations found in the same search (not read): a 2008 ESA abstract on Wisconsin shrubs reporting exotics held leaves "3 to 5 weeks longer" than natives in fall, and a 2024 Int. J. Biometeorology urban-Wisconsin study with non-native shrub senescence start "DOY 262" and end "DOY 300" (both UNVERIFIED, snippet only).
Meaning for Elba: the primary Wisconsin source for buckthorn/honeysuckle extended leaf display exists but its dates are in the body, not the abstract; the 3-5 week figure is the working number until the paper is read.

### S08. Wang et al. 2013, phenology-assisted C3/C4 grass classification, Great Plains, MODIS -- ABSTRACT READ; full text PDF located (open) but not yet read
Wang, C., E.R. Hunt Jr., L. Zhang, and H. Guo. 2013. Phenology-assisted classification of C3 and C4 grasses in the U.S. Great Plains and their climate dependency with MODIS time series. Remote Sensing of Environment 138: 90-101. doi:10.1016/j.rse.2013.07.025. Open copy: https://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=3363&context=usdaarsfacpub
500 m MODIS MOD09A1 NDVI 2000-2009; metrics "Start of Season, End of Season, Season Length, Peak NDVI, Peak Date, and sum NDVI"; decision tree for four grass PFTs. Accuracy and thresholds not in the abstract. Related (snippet only, UNVERIFIED): Guo et al. 2003-type Landsat TM multitemporal work in the central Great Plains discriminating cool- and warm-season grasslands at "91.5 percent" with a three-date spring-summer-fall stack vs single date.
Meaning for Elba: the C3 (brome, bluegrass, RCG) vs C4 (bluestem, Indiangrass, switchgrass) split is the best-established phenological separation in grasslands: C3 peaks in late spring/early summer and again in fall, C4 peaks in mid-to-late summer. RCG is C3, so phenology alone will NOT separate it from brome/bluegrass hayfields by functional type; it has to be separated by its wetness, its late-October greenness, or its height/biomass.

### S09. Gilmore et al. 2008, multi-temporal QuickBird + field spectra + lidar, Phragmites / Typha / Spartina patens, Connecticut -- FULL TEXT READ (pdftotext)
Gilmore, M.S., E.H. Wilson, N. Barrett, D.L. Civco, S. Prisloe, J.D. Hurd, and C. Chadwick. 2008. Integrating multi-temporal spectral and structural information to map wetland vegetation in a lower Connecticut River tidal marsh. Remote Sensing of Environment 112, doi:10.1016/j.rse.2008.05.020. PDF: https://media.clear.uconn.edu/publications/research/tech_papers/Gilmore_et_al_JRSE2008.pdf
Ragged Rock Creek marsh, Old Saybrook CT (41.3 N). Field spectra 2004-2006 reduced to QuickBird band ratios vs growing degree days (GDD50): "All species showed a rise in these indices at approximately GDD50 400 (early June) corresponding with the green-up phase of plant growth and a subsequent decline in the indices corresponding to senescence. S. patens and Phragmites reached peak values at ~ GDD50 >= 900 (mid-July) and Typha spp. at ~ GDD50 500 (mid- to late June). Typha spp. NDVI and NIR/red values were generally higher than the other species near the time of their peak while Phragmites values exceeded the other species in mid-August through early September." Red/green ratio (a senescence index): "Typha spp. values exceeded those of the other species from mid-July onward." Classification (object-based, three dates + lidar height): Phragmites objects = high NIR/red in the 12 September 2004 image and lidar height 1.1-4.0 m; Typha objects = high red/green in the 13 August 2006 image and middle lidar heights; S. patens = high green/blue on 20 July 2004 and low height. Fuzzy accuracy: "Phragmites had the highest accuracy with both MAX and RIGHT values of 97%. Typha spp. was the least accurate with 63% MAX and 76% RIGHT." Cites Laba et al. (2005), Hudson River: purple loosestrife, Phragmites and cattail "best differentiated in late August."
Meaning for Elba: the one quantitative phenology-separation study of Typha vs a tall late-green grass (Phragmites here stands in for RCG): Typha peaks in June and begins browning by mid-July, while the tall grass holds NIR into September. Lidar canopy height (1.1-4.0 m) was a load-bearing input; RCG at 0.6-2 m sits in the same height band as Typha, so for RCG the height cue is weaker than for Phragmites.

### S10. Manske 2017/2018 (NDSU Dickinson), Kentucky bluegrass autecology, 67 seasons of phenology -- FULL TEXT READ (pdftotext)
Manske, L.L. 2017. Autecology of Kentucky Bluegrass on the Northern Mixed Grass Prairie. NDSU Dickinson Research Extension Center, Range Research Report DREC 17-1173 (in the 2018 Annual Report). https://www.ndsu.edu/agriculture/sites/default/files/2022-10/28.%202018%20Annual%20Report%20Drec%201173_Kentucky_Bluegrass.pdf
Dickinson, western North Dakota (46.9 N) and the Dakota Sandhills (Ransom/Richland Co., ~46.3 N), 1946-2012. "Vegetative growth of Kentucky bluegrass occurs very early in the spring. Green carryover leaves can be found during snow melt in mid March. Rapid growth resumes during early April." Table 1 (Whitman et al. 1951): 13 Apr "Early leaf greenup"; 25 Apr "Flower stalk developing"; 4 May "Flower stalk emerging"; 16 May "Flowering (Anthesis)"; 6 Jun "Seed developing"; 19 Jun "Seed maturing"; 26 Jun "Drying"; 2 Jul "Seed mature"; 16 Jul "Seed shedding"; 13 Aug, 5 Nov "Drying". 1976-78: "Flower stalk emergence occurred during mid May. The first flower stalks reached anthesis between 25 May and 5 June ... mature seeds were easily removed from stalks during mid to late June. The stalks dried and were easily removed by wind by mid July." Heights (Table 4, ungrazed): basal leaves "mean maximum mature height of 17.4 cm (6.9 in) during late June on the mixed grass prairie and of 25.1 cm (9.9 in) during mid June on the tall grass prairie"; flower stalks "52.2 cm (20.6 in) during late June." Senescence (Table 5): upland lead-tiller leaves "rapid senescence during late June with 42% of the leaves drying"; "Fall tiller growth began in mid August continuing through late August."
Meaning for Elba: at 46-47 N, bluegrass is green by mid-April, flowers 16 May-5 June, is 42 % dry by late June, and regrows from mid-August. A bluegrass/brome field therefore has a spring NDVI peak in May, a summer dip, and an autumn regreening, which is the same two-peak shape the Wisconsin guide gives for RCG (S11). The separation between RCG and these C3 hay grasses is NOT in the shape of the season; the Wisconsin mapping relied on RCG staying green into mid/late October when wetland natives had browned, and it never tested RCG against brome/bluegrass.

### S11. Wisconsin Reed Canary Grass Management Working Group 2009, management guide PUB-FR-428 -- FULL TEXT READ (pdftotext)
Wisconsin Reed Canary Grass Management Working Group. 2009. Reed Canary Grass (Phalaris arundinacea) Management Guide: Recommendations for Landowners and Restoration Professionals. Wisconsin DNR PUB-FR-428. https://wmswcd.org/wp-content/uploads/2015/06/RCG-management.pdf (also served by WDNR at widnr.widen.net)
Life cycle, verbatim: "Growth and productivity peak twice during the growing season, first in late spring and again in late summer. These growth peaks are under separate genetic control, with leaf and inflorescence growth dominating in the spring and stem and rhizome growth dominating during the late summer peak. RCG is one of the first wetland plants to emerge in the spring, enabling it to shade out native species that emerge later in the growing season. RCG can stay green and actively growing well past the first killing frost in autumn." "Once established, RCG emerges in the spring from rhizome reserves ... Since most spring growth occurs aboveground, the rhizome becomes depleted of starch until flowering. After flowering, rhizomes elongate and tiller." "RCG is biennial with respect to flowering ... The new stems that develop from seed or rhizome buds require two years to develop panicles." "both green and purple panicles ... in mid-June"; rhizomes "over 10 feet per year and form a thick impenetrable mat below the soil surface." Management table: mow/hay "before RCG seed heads appear (boot to late boot stage)"; herbicide "In early spring or late fall" when other plants are dormant.
Meaning for Elba: the clearest agency statement of the two-peak season and of greenness "well past the first killing frost" (first frost at Elba is typically early-to-mid October; UNVERIFIED here). Also the biennial-flowering note: a hayed stand cut at boot stage in June will show few panicles.

### S12. Delaney et al. 2024, RCG habitat suitability in Upper Mississippi floodplain forest understories -- ABSTRACT ONLY (USGS page; Wiley 403)
Delaney, J.T., M. Van Appledorn, N.R. De Jager, K.L. Bouska, and J.J. Rohweder. 2024. Spatial differences in predicted Phalaris arundinacea (reed canarygrass) occurrence in floodplain forest understories. Ecosphere 15(12) [issue from search snippet], doi:10.1002/ecs2.70138. https://www.usgs.gov/publications/spatial-differences-predicted-phalaris-arundinacea-reed-canarygrass-occurrence
Ensemble SDM (BART, boosted trees, random forest) over ~41,000 ha of UMR floodplain forest plots; most important predictors "flood frequency, distance to contiguous floodplain, distance to forest edge, and distance to invaded wet meadow"; probability of presence declined with greater flood depth, distance to contiguous floodplain, distance to invaded wet meadow, forest cover and forest height. Companion data release: habitat suitability Pools 3-13.
Meaning for Elba: RCG under the floodplain-forest canopy is predicted from proximity to invaded wet meadow and low canopy cover; it is not mapped from imagery. The Whitewater delta at Weaver is in Pool 5, inside the modeled reach (UNVERIFIED that the tributary valley itself is covered).

### S13. Minnesota Wildflowers (minnesotawildflowers.info), species pages -- READ (web); a Minnesota-specific calendar with no latitude finer than "MN"
Phalaris arundinacea: height "2 to 5 feet"; "Fruiting season: late spring to early summer"; dormant-season diagnosis: the upper vegetative portion "fades to a pale, bleached tan that is diagnostic in the dormant season identification"; dead leaves "turn a bleached tan color and persist through winter"; habitat "wetlands, fields, woods"; Wetland Indicator FACW. https://minnesotawildflowers.info/grass-sedge-rush/reed-canary-grass
Bromus inermis: height "20 to 50 inches"; "Fruiting season: July - September"; "the only rhizomatous Brome in Minnesota." https://www.minnesotawildflowers.info/grass-sedge-rush/smooth-brome
Poa pratensis: height "2 to 36 inches"; "Fruiting season: June - August." https://www.minnesotawildflowers.info/grass-sedge-rush/kentucky-bluegrass
Andropogon gerardii: "2 to 7 feet"; "Fruiting season: July - August." https://www.minnesotawildflowers.info/grass-sedge-rush/big-bluestem
Sorghastrum nutans: "3 to 7 feet"; "Fruiting season: August - September." https://www.minnesotawildflowers.info/grass-sedge-rush/indian-grass
Solidago canadensis: "2 to 5 feet"; bloom "July - October". https://www.minnesotawildflowers.info/flower/canada-goldenrod
Cirsium arvense: "3 to 5 feet"; bloom "June - October". https://www.minnesotawildflowers.info/flower/canada-thistle
Rhamnus cathartica: "to 20 feet"; bloom "June - July"; "the leaves do not change color in the fall, persisting green into late fall"; "In late fall and early winter the persistent green foliage is the most standout characteristic." https://www.minnesotawildflowers.info/shrub/common-buckthorn
Typha latifolia: "3 to 10 feet"; bloom "May - July"; pistillate spikes "persist mostly intact through winter." https://www.minnesotawildflowers.info/flower/broad-leaf-cat-tail
Acer saccharinum: "60 to 100 feet"; bloom "March - April", flowers "in late March to early April, before leaves emerge"; "Leaves turn yellow in fall." https://www.minnesotawildflowers.info/tree/silver-maple
Acer negundo: "30 to 50 feet"; bloom "April - May"; flowers "appear before or as the leaves emerge in early spring"; "Leaves turn a dull yellow in fall." https://www.minnesotawildflowers.info/tree/boxelder
Populus deltoides: "80 to 120 feet"; bloom "April - May". https://www.minnesotawildflowers.info/tree/plains-cottonwood
Salix interior: "6 to 22 feet"; bloom "May - July"; "first blooms in early May but will continue putting out catkins into July." https://www.minnesotawildflowers.info/shrub/sandbar-willow
Cornus sericea: "3 to 8 feet"; bloom "May - August"; dormant-season branches "deep maroon red." https://www.minnesotawildflowers.info/shrub/red-osier-dogwood
Quercus macrocarpa: "to 100 feet"; bloom "early spring"; leaves "turn a golden brown in fall"; editor's comment: bur oak "breaking bud in late April to early May". https://www.minnesotawildflowers.info/tree/bur-oak
Meaning for Elba: the RCG "bleached tan ... persist through winter" line is the only Minnesota statement found on the standing dead canopy, i.e. what the November 2008 lidar saw. The site records fruiting rather than bloom for grasses and gives no green-up or senescence dates.

### S14. Minnesota Seasons (minnesotaseasons.com), reed canary grass and smooth brome -- READ (web)
Reed canary grass: height "24 to 80 inches"; flowering "May to August"; sightings dated 6/19/2022 (Apple Blossom Overlook Park, which is in Winona County ~44.0 N, UNVERIFIED location), 11/5/2021 (Lake Shore, MN) and June 2026; inflorescence "When immature, it is compact and resembles a spike ... When in full flower, it is open and green or slightly purple or reddish ... Later, it turns tan to whitish in color and closes up, becoming narrow and pointed." Look-alikes Dactylis glomerata and Calamagrostis canadensis. https://minnesotaseasons.com/Plants/reed_canary_grass.html
Smooth brome: height "36 to 48 inches"; flowering "Early to late June"; sightings include Whitewater State Park and Kellogg Weaver Dunes SNA (both within 15 km of Elba; no dates). https://minnesotaseasons.com/Plants/smooth_brome.html
Meaning for Elba: smooth brome is documented at Whitewater State Park, so brome is a confirmed candidate for the hayed field, not a hypothetical.

### S15. UMN Extension, reed canary grass -- READ (web)
Gupta, A., A. Rager, and M.M. Weber. Reed canary grass. University of Minnesota Extension, reviewed 2022. https://extension.umn.edu/invasive-species/reed-canary-grass
"grows two to six feet high"; "Typically bloom from May to mid-June"; seed heads "green to purple when in bloom" turning tan-beige; distinguished from bluejoint (Calamagrostis canadensis) by the "Prominent, highly transparent membrane" ligule. No fall statement.

### S16. Kolecka et al. 2018, mowing detection with Sentinel-2 -- ABSTRACT READ (Crossref)
Kolecka, N., C. Ginzler, R. Pazur, B. Price, and P.H. Verburg. 2018. Regional Scale Mapping of Grassland Mowing Frequency with Sentinel-2 Time Series. Remote Sensing 10(8): 1221. doi:10.3390/rs10081221
Canton of Aargau, Switzerland (47.4 N). Per the abstract (summarized by the fetch tool): the best configuration (extra cloud masking + shrunken parcels) identified 77 % of mowing events; standard cloud masking over-detected; sparse time series missed events; over 40 % of the area was mown before 15 June. Search snippets (UNVERIFIED): Kolecka's rule was an NDVI drop > 0.2 between consecutive clear dates; Lobert et al. 2021 (Germany, RSE, doi:10.1016/j.rse.2021.112751) reached F1 0.84 with NDVI + Sentinel-1, mean date error 3.8 days.
Meaning for Elba: whether the field was hayed, and when, in 2008 and 2021 is detectable in principle from an NDVI drop; the Landsat-5 cadence in 2008 (16 days, cloud-limited) will miss many events, Sentinel-2 in 2021 will catch most. A cut date changes what the lidar saw more than the species does.

### S17. Sahramaa, Ihamaki & Jauhiainen 2003, RCG biomass traits, Finland -- FULL TEXT READ (pdftotext); off-region, used only for the overwinter-standing question
Sahramaa, M., H. Ihamaki, and L. Jauhiainen. 2003. Variation in biomass related variables of reed canary grass. Agricultural and Food Science in Finland 12: 213-225. https://journal.fi/afs/article/download/5749/4947/13929
Jokioinen, Finland (60.8 N), 1994-1998, 53 wild populations. Delayed (spring) harvest system: stands left standing over winter and "harvested to stubble height of 10-12 cm in May"; plant height measured "after seed ripening in the end" of the season; overwintering ability up to 95 %; biomass over 13 t DM/ha in year four. The Nordic bioenergy literature (search snippets, UNVERIFIED: Swedish trials 1991-1995 8.5 t DM/ha in autumn, ~7 t DM/ha recovered the next spring) exists because RCG stands STAY STANDING through winter as a dry crop; that is the agronomic basis for the "standing dead canopy in winter" that a leaf-off lidar sees. Latitude 60 N; stem lodging under snow at 44 N is not addressed.

### S18. Manske 2017 (NDSU Dickinson), big bluestem autecology -- FULL TEXT READ (pdftotext)
Manske, L.L. 2017. Autecology of Big Bluestem on the Northern Mixed Grass Prairie. NDSU Dickinson REC Range Research Report DREC 17-1167. https://ndsu.edu/agriculture/sites/default/files/2022-10/18.%202018%20Annual%20Report%20Drec%201167_Big_Bluestem.pdf
Dickinson ND (46.9 N). "Flower period is from late July to late August." "Growth of new Big bluestem leaves is visible by 22 May. Rapid leaf growth starts 6 June and continues through June." "During mid June, leaf height is 55% of maximum at 11.4 cm (4.5 in) tall." Leaf "100% maximum height during early August." Stalk "boot stage between 8 and 17 July, reaches head emergence by 29 July, produces the first flowers on 24 July, mean first flowers occur on 4 August, with a 3 week flower period from late July to mid August. Seeds are developing from 6 August, seeds are maturing around 26 August, reach the mature stage between 30 August and 12 September, and start being shed on 7 September." Table 5 (Goetz 1963, 1955-1962) percent leaf dryness: leaf tip dry 11 Jul; 0-25 % 20 Jul; 25-50 % 7 Aug; 50-75 % 20 Aug; 75-100 % 11 Sep.
Meaning for Elba: the C4 reference. Big bluestem is barely up on 22 May and 11 cm in mid-June, when bluegrass is already flowering and 17-25 cm; its canopy peaks in August. A 1 May image (the 2021 lidar date) sees bare ground/litter on a C4 prairie and a green C3 sward.

### S19. Manske 2017 (NDSU Dickinson), smooth bromegrass autecology -- FULL TEXT READ (pdftotext)
Manske, L.L. 2017. Autecology of Smooth Bromegrass on the Northern Mixed Grass Prairie. NDSU Dickinson REC Range Research Report DREC 17-1172. https://ndsu.edu/agriculture/sites/default/files/2022-10/24.%202018%20Annual%20Report%20Drec%201172_Smooth_Bromegrass.pdf
Dickinson ND (46.9 N). "Smooth bromegrass starts early leaf greenup in mid April (table 1). The smooth bromegrass tillers have three and a half new leaves around early May which is four or so weeks earlier than native cool season grasses ... Early boot stage occurs in mid May, stalk emergence occurs in early June." Table 1 (Whitman et al. 1951): 13 Apr "Early leaf greenup"; 4 May "Active leaf growth"; 16 May "Flower stalk developing"; 6 Jun "Flower stalk emerging"; 13 Jun "Flowering (Anthesis)"; 19 Jun "Seed developing"; 26 Jun "Seed maturing"; 2 Jul "Drying"; 8 Jul "Seed mature"; 26 Aug, 5 Nov "Drying". Flower stalks "40-120 cm (16-47 in) tall." "Flower period is from early to late June." "Most of the lead tillers reach the flower stage during June."
Meaning for Elba: brome and RCG green up in the same weeks (mid-April at 47 N) and flower within a fortnight of each other (brome 13 June; RCG late May to mid-June); heights overlap (brome to 1.2 m, RCG 0.6-2 m). Spring and early-summer imagery cannot separate them; the candidates are soil wetness, late-October greenness, and the taller/denser standing litter of RCG.

### S20. Kurtz & Hansen 2014, USFS FIA assessment of reed canarygrass in northern forests -- FULL TEXT READ (pdftotext)
Kurtz, C.M., and M.H. Hansen. 2014. An assessment of reed canarygrass in northern U.S. forests. Research Note NRS-200. USDA Forest Service, Northern Research Station. https://www.govinfo.gov/content/pkg/GOVPUB-A13-PURL-gpo51806/pdf/GOVPUB-A13-PURL-gpo51806.pdf
FIA plots 2008-2012, 24 states: reed canarygrass "was found on 611 plots (8 percent) and occurred in 21 of the 24 states ... Iowa had the highest percentage of plots with reed canarygrass (45 percent)." No phenology. Confirms that RCG is a forest-understory species across the region, which matters for the bluff-forest lidar bias as well as the floodplain.

### S21. Greenstein 2021, Reed Canarygrass State of the Science Summary (Kenai Peninsula CISMA) -- FULL TEXT READ (pdftotext); secondary review
Greenstein, C. (compiler), K. Schake (ed.). 2021. Reed Canarygrass (Phalaris arundinacea) State of the Science Summary, July 2021. Homer Soil & Water Conservation District / Kenai Peninsula CISMA. https://kenaiinvasives.org/wp-content/uploads/gravity_forms/3-2cb8e5dbe3e34c933a71aba9ddcc3948/2021/09/KPCISMA_ReedCanarygrass_ScienceSummary_Final_2021.pdf
Biology restates WRCGMWG 2009: "can grow up to six feet tall and is biennial in regard to flowering. Growth peaks twice, in late spring and late summer, with a mid-summer dormancy period." Its "Remote mapping" section (pp. 13-14) summarizes S01 and S02 accurately (737,454 acres; 86 %; mid-October; NDVI not needed; 498,250 acres; 61-83 % / 72-92 %; grazed pasture and cranberry bogs as commission errors) and concludes for Alaska that "remote sensing options captured at higher resolution (i.e. imagery collected via fixed wing airplane) during October could be further investigated." It cites no RCG remote-sensing study other than the two Wisconsin reports, which corroborates that none of comparable scope exists as of 2021.

## (c) Regional phenology data sources

### D01. USGS UMRR LTRM 2020 land cover, reed canarygrass remap, Pools 2-13 -- landing page READ
Hoy, E.E., and J.E. Cogan. 2026. Upper Mississippi River Restoration (UMRR) Program Long Term Resource Monitoring (LTRM) - 2020 LCU Mapping Reed Canarygrass Remap - Upper Mississippi River System (UMRS) Pools 02-13. USGS data release, https://doi.org/10.5066/P1YQT2ZH
What it provides: a remap of the 2010 product in which reed canarygrass is a DOMINANCE ATTRIBUTE (">80% cover") on the Wm (wet meadow) and Wms (wet meadow shrub) land cover classes of the 2020 LTRM land cover/land use (LCU) mapping, Pools 2-13 of the Mississippi (Pool 5/5A at the Whitewater mouth lies inside 2-13 but the page does not list pools individually; UNVERIFIED that the Whitewater delta itself is mapped). Imagery dates, resolution, MMU and accuracy are not on the landing page; the LCU series is photo-interpreted aerial imagery (UNVERIFIED from this page). Page notes high inundation in the years before 2020 may have changed the RCG distribution. Useful as an independent, photo-interpreted RCG map adjacent to (not on) the Elba reach, i.e. a training/validation source for a Sentinel-2 classifier 25 km downstream.

### D02. MODIS MCD12Q2 v061 Land Cover Dynamics (land surface phenology), 500 m, 2001-present -- product page READ
https://doi.org/10.5067/MODIS/MCD12Q2.061 ; https://www.earthdata.nasa.gov/data/catalog/lpcloud-mcd12q2-061
Yearly phenometrics from NBAR EVI2, up to two cycles per pixel per year: Greenup (15 % of amplitude), MidGreenup (50 %), Maturity (90 %), Peak, Senescence (last 90 %), MidGreendown (50 %), Dormancy (15 %), EVI_Amplitude, EVI_Minimum, EVI_Area, NumCycles, QA. At 500 m a pixel is 25 ha, which is larger than the Elba field; useful only as a regional climatology of green-up/dormancy dates to anchor the Sentinel-2 dates, not for the field itself.

### D03. MSLSP30NA v011, Multi-Source (HLS) Land Surface Phenology, 30 m, North America, 2016-2019 -- product page READ
https://doi.org/10.5067/COMMUNITY/MUSLI/MSLSP30NA.011 ; https://www.earthdata.nasa.gov/data/catalog/lpcloud-mslsp30na-011 ; methods: Bolton et al. 2020, Remote Sensing of Environment 240:111685, doi:10.1016/j.rse.2020.111685 ; user guide https://lpdaac.usgs.gov/documents/652/MSLSP_User_Guide_V1.pdf
30 m yearly phenometrics from Harmonized Landsat 8 + Sentinel-2A/B (HLS): OGI/OGD (onset greenness increase/decrease), 50PCGI/50PCGD, OGMn/OGMx, Peak, EVIamp, EVImax, EVIarea, numObs, QA, and a second cycle (_2). Known issues: sparse 2016 and early 2017 (Sentinel-2B not yet launched); disturbance (mowing, haying) not modeled, "sharp spectral changes can cause premature detections of senescence and dormancy"; pixels with more than two cycles (alfalfa) poorly characterized. For Elba this covers 2016-2019 only, so it brackets neither lidar date; but it is a ready 30 m green-up/senescence map to compare the field against its neighbors, with the caveat that a hayed field will show a false early senescence.

### D04. PhenoCam Network sites in Minnesota and Wisconsin -- site table READ (first ~2/3 of table)
https://phenocam.nau.edu/webcam/network/table/ ; Marcell: https://phenocam.nau.edu/webcam/sites/marcell (imagery 2014-05-21 to 2025-09-15)
Minnesota: arsmnswanlake1 (45.68, -95.80), arsmorris1/2/3 (45.62, -96.13), boglakepeatland (47.51, -93.49), boundarywaters (47.95, -91.50), cedarcreek/cedarcreek2 (45.40, -93.20), macordway (44.81, -93.03), marcell (47.51, -93.47), rosemount, rosemountc6, rosemountcons, rosemountconv, rosemountg21, rosemountnprs (44.68-44.73, -93.06 to -93.09). Wisconsin: lostcreek (46.08, -89.98), NEON STEI and TREE (45.5, -89.59), NEON CRAM (46.21, -89.48), NEON LIRO (46.00, -89.70), kempnrs (45.84, -89.68), downerwoods (43.08, -87.88), kegonsasolar (42.96, -89.29). Nearest to Elba (44.12 N, 92.00 W): the Rosemount cluster and macordway, ~105-112 km NW, all agricultural or restored prairie (rosemountnprs = Native Prairie Restoration, i.e. a C4-grass reference camera). No PhenoCam on a floodplain wet meadow or an RCG stand; no camera in the Driftless Area. Vegetation-type codes and date ranges are not in the table view; they are on each site page. Remainder of the table (second read): Minnesota adds the SPRUCE cluster at Marcell (47.50-47.51, -93.45) and snakerivermn (Hay-Snake State WMA, 46.12, -93.24); Wisconsin adds sylvania (46.24, -89.35), willowcreek (45.81, -90.08), usdadfrc (USDA Dairy Forage Research Center, Prairie du Sac, 43.34, -89.71; a managed cool-season forage reference ~200 km from Elba) and uwmfieldsta (Saukville, 43.39, -88.02). No Iowa sites; nothing within 150 km of Elba beyond the Rosemount/Ordway cluster.

### D05. USA National Phenology Network (Nature's Notebook) observations, MN + WI + IA, 2009-2026 -- DOWNLOADED and summarized this session

API: https://services.usanpn.org/npn_portal/observations/getSummarizedData.json (request_src=wickert_lit_review, start_date=2009-01-01, end_date=2026-11-18, state MN/WI/IA). Species list https://services.usanpn.org/npn_portal/species/getSpecies.json confirms Nature's Notebook profiles for Phalaris arundinacea (species_id 960), Bromus inermis (1441), Poa pratensis (971), Elymus repens (2313), Andropogon gerardii (129), Sorghastrum nutans (1695), Panicum virgatum (137), Typha latifolia (1697), Solidago canadensis (1643), Cirsium arvense (156), Rhamnus cathartica (89), Lonicera maackii (1246), L. tatarica (770), Acer negundo (2), A. saccharinum (781), A. saccharum (61), Populus deltoides (976), Salix nigra (1494), Cornus sericea (70), Quercus macrocarpa (101), Q. rubra (102), Q. alba (100), Juniperus virginiana (44), Tilia americana (93). Raw JSON for five species saved in `analysis/rcg/usanpn/sp_<species_id>.json`; the derived table `analysis/rcg/usanpn/usanpn_earliest_first_yes_per_site_year.csv` holds, for each (species, phenophase, site, year), the EARLIEST first-yes day of year (a Nature's Notebook 'first yes' is the first 'yes' after a 'no', so a season can have several; the earliest per site-year is the onset). Record counts: reed canarygrass 280, smooth brome 245, Kentucky bluegrass 10, common buckthorn 666, big bluestem 319.

Caveat on representativeness: nearly all RCG and brome records come from ONE site (site_id 25568, 45.14 N, Minnesota), observed by the same person in the same years, which makes the RCG-vs-brome comparison below unusually clean (same place, same eyes, same days) but a sample of one stand each. Other RCG sites: 20617 (43.13 N, WI, 2018).

Per-phenophase medians of the per-site-year earliest first-yes DOY (n site-years; DOY median [min-max]; latitude range):

- reed canarygrass:
    - Flower heads (grasses/sedges): n=13, median DOY 162 [146-238], lat 41.56-45.14
    - Fruits: n=13, median DOY 173 [53-222], lat 44.74-45.14
    - Initial growth (grasses/sedges): n=13, median DOY 112 [95-238], lat 41.56-45.14
    - Leaves (grasses): n=17, median DOY 134 [111-282], lat 41.56-45.17
    - Open flowers (grasses/sedges): n=13, median DOY 167 [153-238], lat 41.56-45.14
    - Pollen release (flowers): n=10, median DOY 176.0 [155-230], lat 44.74-45.14
    - Recent fruit or seed drop: n=12, median DOY 209.0 [155-291], lat 44.74-45.14
    - Ripe fruits: n=13, median DOY 185 [53-291], lat 44.74-45.14
- smooth brome:
    - Flower heads (grasses/sedges): n=11, median DOY 151 [126-183], lat 45.14-45.62
    - Fruits: n=10, median DOY 170.5 [146-278], lat 45.14-45.14
    - Initial growth (grasses/sedges): n=9, median DOY 110 [92-127], lat 45.14-45.14
    - Leaves (grasses): n=12, median DOY 127.5 [109-225], lat 45.14-45.62
    - Open flowers (grasses/sedges): n=9, median DOY 169 [153-234], lat 45.14-45.14
    - Pollen release (flowers): n=7, median DOY 169 [166-234], lat 45.14-45.14
    - Recent fruit or seed drop: n=8, median DOY 214.5 [185-251], lat 45.14-45.14
    - Ripe fruits: n=8, median DOY 178.0 [153-207], lat 45.14-45.14
- Kentucky bluegrass:
    - Flower heads (grasses/sedges): n=1, median DOY 158 [158-158], lat 45.57-45.57
    - Initial growth (grasses/sedges): n=1, median DOY 121 [121-121], lat 47.49-47.49
    - Leaves (grasses): n=2, median DOY 139.5 [121-158], lat 45.57-47.49
    - Open flowers (grasses/sedges): n=1, median DOY 158 [158-158], lat 45.57-45.57
- big bluestem:
    - All leaves withered (herbaceous and grasses): n=1, median DOY 283 [283-283], lat 44.59-44.59
    - Flower heads (grasses/sedges): n=24, median DOY 206.5 [188-305], lat 41.59-47.0
    - Fruits: n=20, median DOY 225.0 [9-305], lat 41.57-47.0
    - Initial growth (grasses/sedges): n=29, median DOY 133 [77-298], lat 41.59-47.0
    - Leaves (grasses): n=38, median DOY 166.5 [99-305], lat 41.55-47.0
    - Open flowers (grasses/sedges): n=23, median DOY 211 [184-305], lat 41.59-47.0
    - Pollen release (flowers): n=16, median DOY 221.0 [201-305], lat 41.57-47.0
    - Recent fruit or seed drop: n=20, median DOY 271.0 [9-305], lat 41.57-47.0
    - Ripe fruits: n=21, median DOY 258 [9-305], lat 41.57-47.0
    - Ripe seeds (grasses): n=2, median DOY 225.5 [211-240], lat 44.59-45.8
- common buckthorn:
    - >=50% of leaves colored (deciduous): n=1, median DOY 323 [323-323], lat 42.78-42.78
    - >=75% of full leaf size (deciduous): n=1, median DOY 323 [323-323], lat 42.78-42.78
    - All leaves colored (deciduous): n=1, median DOY 323 [323-323], lat 42.78-42.78
    - Breaking leaf buds: n=39, median DOY 107 [55-136], lat 42.62-47.24
    - Colored leaves: n=28, median DOY 288.0 [97-320], lat 42.96-47.0
    - Early season leaf expansion: n=33, median DOY 126 [97-300], lat 41.2-47.0
    - Falling leaves: n=26, median DOY 290.5 [50-320], lat 42.96-47.0
    - Flowers or flower buds: n=23, median DOY 139 [80-180], lat 42.62-47.0
    - Fruits: n=43, median DOY 190 [21-361], lat 41.2-47.0
    - Leaves: n=52, median DOY 131.0 [50-323], lat 41.2-47.0
    - Open flowers: n=15, median DOY 147 [97-184], lat 42.96-47.0
    - Recent fruit or seed drop: n=28, median DOY 268.5 [21-361], lat 42.64-47.0
    - Ripe fruits: n=37, median DOY 236 [21-361], lat 42.64-47.0

RCG vs smooth brome at the same site (25568, 45.14 N), earliest first-yes DOY per year:

- Flower heads (grasses/sedges) 2017: RCG 149 / brome 131
- Flower heads (grasses/sedges) 2018: RCG 154 / brome 126
- Flower heads (grasses/sedges) 2019: RCG 167 / brome 167
- Flower heads (grasses/sedges) 2020: RCG 183 / brome 183
- Flower heads (grasses/sedges) 2021: RCG 162 / brome 145
- Flower heads (grasses/sedges) 2022: RCG 163 / brome 131
- Flower heads (grasses/sedges) 2023: RCG 149 / brome 134
- Flower heads (grasses/sedges) 2024: RCG 146 / brome 151
- Flower heads (grasses/sedges) 2025: RCG 156 / brome 156
- Flower heads (grasses/sedges) 2026: RCG - / brome 171
- Initial growth (grasses/sedges) 2017: RCG 112 / brome 127
- Initial growth (grasses/sedges) 2018: RCG 111 / brome 111
- Initial growth (grasses/sedges) 2019: RCG 116 / brome 110
- Initial growth (grasses/sedges) 2020: RCG 207 / brome -
- Initial growth (grasses/sedges) 2021: RCG 101 / brome 92
- Initial growth (grasses/sedges) 2022: RCG 121 / brome 121
- Initial growth (grasses/sedges) 2023: RCG 111 / brome 111
- Initial growth (grasses/sedges) 2024: RCG 104 / brome 104
- Initial growth (grasses/sedges) 2025: RCG 103 / brome 103
- Initial growth (grasses/sedges) 2026: RCG 95 / brome 95
- Open flowers (grasses/sedges) 2017: RCG 176 / brome 166
- Open flowers (grasses/sedges) 2018: RCG 154 / brome 174
- Open flowers (grasses/sedges) 2019: RCG 167 / brome 179
- Open flowers (grasses/sedges) 2020: RCG 186 / brome 186
- Open flowers (grasses/sedges) 2021: RCG 164 / brome 162
- Open flowers (grasses/sedges) 2022: RCG 163 / brome 169
- Open flowers (grasses/sedges) 2023: RCG 155 / brome 161
- Open flowers (grasses/sedges) 2024: RCG 153 / brome 153
- Open flowers (grasses/sedges) 2025: RCG 234 / brome 234

Meaning for Elba (all numbers are the per-site-year medians above): at 45.1 N reed canarygrass and smooth brome start growth at median DOY 112 (n=13) and 110 (n=9), and open flowers at median DOY 167 (n=13) and 169 (n=9); at the shared site the yearly onsets are identical in 6 of 10 years for initial growth, so the spring season cannot separate them. Big bluestem starts at median DOY 133 (n=29) and flowers at median DOY 211 (n=23, 30 July), the C4 signature. Common buckthorn breaks bud at median DOY 107 (n=39, 17 April); its colored leaves have median onset DOY 288.0 (n=28, 15 October) and falling leaves median DOY 290.5 (n=26, 17 October), so it is green when the surrounding canopy has turned, which confirms the late-fall window for shrub detection on the bluff. These are onsets at 42-47 N; Elba (44.1 N) sits inside that band.

## (d) Sources-checked log (including dead ends)

| # | source | result |
|---|--------|--------|
| 1 | Bernthal & Willis 2004 WDNR rcg_mapping_report.pdf | READ full text via pdftotext; S01 |
| 2 | Hatch & Bernthal 2008 WDNR RCGFinalReport10_08.pdf | READ full text via pdftotext; S02 |
| 3 | Waggy 2010 FEIS Phalaris arundinacea | READ (web, 3 fetches); S03 |
| 4 | USDA NRCS Plant Guide pg_phar3.pdf (plants.sc.egov.usda.gov) | DEAD END: HTTP 504 twice; plants.usda.gov redirects to the same host |
| 5 | MN DNR reed canary grass page | READ; S04 |
| 6 | USGS UMRR LTRM 2020 LCU RCG remap landing page (doi 10.5066/P1YQT2ZH) | READ; D01 |
| 7 | USGS pubs 70263847 "Spatial differences in predicted Phalaris arundinacea occurrence in floodplain forest understories" | DEAD END at pubs.usgs.gov (HTTP 403); retry via journal |
| 8 | UMN conservancy "Role of leafing phenology in the invasion of forest ecosystems by Rhamnus cathartica" | DEAD END (HTTP 403); seek Harrington, Brown & Reich 1989 instead |
| 9 | peeref.com Oak Openings buckthorn Landsat phenology 2001-2011 | DEAD END (302 redirect); find via DOI |
| 10 | Lavergne & Molofsky 2004 Crit. Rev. Plant Sci. 23:415-429 doi 10.1080/07352680490505934 | abstract only via search snippet; no phenology numbers in abstract; paywalled |
| 11 | Souder et al. 2023 Not. Bot. Horti Agrobot. 51(3) UAS RCG vigor, Oregon | abstract only via snippet; glyphosate-vigor classes, no RCG-vs-other accuracy; off-region |
| 12 | Resasco et al. 2007 IJRS (colorado.edu PDF) | READ full text via pdftotext; S05 |
| 13 | Becker et al. 2013 Biol. Invasions (Springer, Crossref) | Springer redirects to login; Crossref has no abstract; abstract taken from library catalog snippet; S06 |
| 14 | Harrington, Brown & Reich 1989 Oecologia (Springer; UMN Experts) | Springer login redirect; UMN Experts abstract READ; S07 |
| 15 | 2024 Int. J. Biometeorol. urban Wisconsin shrub senescence doi 10.1007/s00484-024-02693-7 | Springer login redirect; snippet only (DOY 262/300) |
| 16 | Wang et al. 2013 RSE C3/C4 MODIS (digitalcommons.unl.edu) | abstract READ, open PDF located; S08 |
| 17 | Delaney et al. 2024 Ecosphere doi 10.1002/ecs2.70138 (Wiley) | DEAD END HTTP 403 at Wiley; abstract from search snippet only; retry USGS page |
| 18 | Gilmore et al. 2008 RSE (uconn PDF) | PDF saved; text extraction pending |
| 19 | Ghioca-Robrecht et al. 2008 Wetlands 28(4) (BioOne) | DEAD END: empty page returned; snippet: overall 62 %, 14 % of Phragmites was Typha |
| 20 | USA-NPN: usanpn.org/nn/Phalaris_arundinacea and services.usanpn.org getSpeciesByGenus | DEAD END: 404 both; no Phalaris profile confirmed |
| 21 | minnesotaseasons.com reed canary grass | HTTP 429 (rate limited); retry |
| 22 | seasonwatch.umn.edu | homepage only; no dataset description on the landing page |
| 23 | NDSU Dickinson 2018 Kentucky bluegrass report (PDF) | PDF saved; text extraction pending |
| 24 | MCD12Q2 v061 product page | READ; D02 |
| 25 | MSLSP30NA v011 product page | READ; D03 |
| 26 | PhenoCam network site table | READ (both halves); D04 |
| 27 | Gilmore et al. 2008 RSE (uconn PDF) | READ full text via pdftotext; S09 |
| 28 | Manske 2017 DREC 17-1173 Kentucky bluegrass (NDSU PDF) | READ full text; S10 |
| 29 | WRCGMWG 2009 RCG Management Guide PUB-FR-428 (wmswcd.org PDF) | READ full text; S11 (WDNR fact-sheet URL redirects to this same PDF) |
| 30 | Delaney et al. 2024 (USGS publications page) | abstract READ; S12 |
| 31 | minnesotawildflowers.info 15 species pages | READ; S13 (first-guess URLs for cattail, boxelder, cottonwood 404ed; corrected URLs found by search) |
| 32 | minnesotaseasons.com reed canary grass, smooth brome | READ on retry; S14 |
| 33 | extension.umn.edu reed canary grass | READ; S15 |
| 34 | Kolecka et al. 2018 Remote Sensing 10:1221 (MDPI 403; Crossref abstract) | abstract READ; S16 |
| 35 | Lobert et al. 2021 RSE mowing detection | snippet only (F1 0.84, 3.8 d) |
| 36 | Sahramaa et al. 2003 AFS Finland 12:213-225 (journal.fi PDF) | READ full text; S17 |
| 37 | Manske 2017 DREC 17-1167 big bluestem (NDSU PDF) | READ full text; S18 |
| 38 | Manske 2017 DREC 17-1172 smooth bromegrass (NDSU PDF) | READ full text; S19 |
| 39 | Kurtz & Hansen 2014 USFS NRS-200 (govinfo PDF) | READ full text; S20 |
| 40 | Greenstein 2021 Kenai RCG State of the Science (PDF) | READ full text; S21 |
| 41 | MN DNR NWI update, southern MN technical documentation (files.dnr.state.mn.us PDF) | DEAD END: >10 MB, fetch refused; east-central doc snippet: wetland-class overall accuracy 66.87 %, RCG used as a dominant-vegetation descriptor |
| 42 | Wang et al. 2013 full PDF (digitalcommons.unl.edu viewcontent) | DEAD END HTTP 403 |
| 43 | Season Watch /explore | not fetched (homepage had no dataset description); low priority |
| 44 | Klopatek & Stearns 1978 Am. Midl. Nat. 100:320-332 (RCG Wisconsin marsh productivity) | not fetched; cited via FEIS for WI seed maturation late June and fall regrowth |
| 45 | Hoffman & Kearns 1997 WDNR invasive plant manual RCG fact sheet | not fetched separately; cited via FEIS (WI flowering May-mid June; mid/late-summer collapse into a mat) |
| 46 | Oak marcescence (search) | blog/extension-level only (MN red oaks hold 15-20 % of leaves in late Jan; marcescence common in red and white oaks of southern MN); no primary source or remote-sensing study found |
| 48 | USA-NPN species list + getSummarizedData for species 960, 1441, 971, 89, 129 in MN/WI/IA | DOWNLOADED via curl; D05; raw JSON and derived CSV in analysis/rcg/usanpn/ |
| 47 | Floodplain tree leaf-out dates (search) | no dated primary source found; Big River Magazine (Winona State ecologist C. Jefferson): silver maple flowers late Feb-early March, leaves later; USGS songbird study: leaf-out earlier in warm 2010/2012, later in cold 2011/2013 |

## (e) Open questions the literature does not answer

1. No study separates reed canary grass from smooth brome, Kentucky bluegrass or quackgrass with any sensor. The Wisconsin maps (S01, S02) masked uplands with WISCLAND and compared RCG only with wetland natives; their one named false positive is grazed pasture. A hayed floodplain field of C3 grasses is exactly the case that was never tested.
2. No dated senescence record for RCG in the region. Every source says "late", "well past the first killing frost", "green into mid/late October" (S01, S02, S03, S11); none gives a day of year, a latitude and a year together. Nature's Notebook has no senescence phenophase for grasses, so D05 cannot supply it either.
3. The standing dead canopy in winter: Minnesota Wildflowers (S13) says the bleached tan upper portion "persist[s] through winter"; the Nordic delayed-harvest literature (S17) shows RCG stands through winter at 60 N; nothing measured at 44 N says how tall the November litter/standing-dead layer is or whether it lodges before snow. This is the quantity the November 2008 lidar saw.
4. No Sentinel-2, NAIP, UAV or hyperspectral RCG classification with reported accuracy was found for the Upper Midwest (searches 2015-2026, standard and extended). The Kenai review (S21) in 2021 also cites only the two Wisconsin Landsat reports. RCG remote sensing after 2008 is a gap, not a hidden literature.
5. Height as a separator: Gilmore (S09) needed lidar height to split Phragmites (1.1-4.0 m) from Typha; RCG at 0.6-2 m overlaps Typha, brome (to 1.2 m) and goldenrod (0.6-1.5 m), so canopy height alone will not do it.
6. Latitude shift of the November shrub window (S05 is 39.5 N; D05 gives buckthorn leaf-fall onset median DOY 290.5 at 42-47 N) and of the mid-October RCG window (S01, S02 at 42.5-47 N) to the exact Elba dates in 2008 and 2021 was not found; it has to come from the imagery itself.
7. Oak marcescence and floodplain tree leaf-out dates: only blog/extension-level statements were found (log items 46, 47); no primary dated source for the Driftless Area.
8. Mowing in 2008: Landsat-5 cadence (16 d, clouds) may not resolve a June cut on the Elba field; Sentinel-2 in 2021 can (S16). Whether the field was hayed in each lidar year is a question for the imagery and the landowner, not the literature.

## (f) Condensed reading of the table (written after the search)

What separates reed canary grass from other floodplain cover at 10-30 m, by date, according to the sources read:

- Mid-to-late OCTOBER greenness (9-31 Oct at 42.5-47 N) is the only operational signature in the literature: a single Landsat scene then gave 86 % user's accuracy for >=80 % RCG cover and 72-92 % for >=50 % cover, against senesced wetland natives (S01, S02). It fails for RCG under canopy, confuses RCG with grazed pasture and cranberry, and was never tested against brome/bluegrass hayfields.
- SPRING does not separate RCG from the C3 hay grasses: initial growth median DOY 112 (RCG) vs 110 (brome) and open flowers DOY 167 vs 169 at the same Minnesota site (D05); mid-April greenup and June anthesis for brome and bluegrass at 47 N (S10, S19); RCG flowering May to mid-June (S03, S04, S15). Spring DOES separate C3 from C4 cover: big bluestem initial growth median DOY 133, flowers DOY 211 (D05); 11 cm tall in mid-June at 47 N (S18). On 1 May 2021 a C3 field is green and a C4 prairie is litter.
- TYPHA peaks in June and browns from mid-July (red/green ratio highest from mid-July onward), while a tall late-green grass holds NIR into September (S09); Typha stalks and spikes stand through winter (S13). So RCG vs cattail is an August-September contrast; RCG vs brome is not.
- The TWO-PEAK season (late spring and late summer growth peaks, fall regrowth after a cut) is shared by RCG (S03, S11) and bluegrass/brome (S10, S19), so a second NDVI rise in September is not diagnostic of RCG; a hayed field regreens regardless.
- WOODY INVADERS: buckthorn breaks bud at median DOY 107 and does not begin leaf fall until median DOY 290.5 (D05; S13 "persisting green into late fall"); Amur honeysuckle was separable only on 4-22 November Landsat scenes at 39.5 N, not in October, because canopy leaves had senesced but not dropped (S05). The November 2008 lidar is inside this window on the bluff.
- WINTER STANDING DEAD: RCG's upper canopy fades to bleached tan and persists through winter (S13); its litter is >80 % of May biomass (S03); Nordic growers harvest the standing dry crop in spring (S17). No height or lodging measurement at 44 N was found.
