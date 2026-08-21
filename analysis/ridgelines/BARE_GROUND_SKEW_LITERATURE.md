# Why gen1 bare-ground returns are positively skewed — literature synthesis

Question (Andy): the semilog plot shows gen1 bare/open ground returns are skewed (heavy
upper tail); find the literature explanation. Agent search 2026-08-22; citations verified to
exist (several primary PDFs paywalled, corroborated via secondary sources).

## Our empirical finding (the thing to explain)
gen1 (2008 discrete-return, leaf-off Nov), ground-class returns vs local bare-earth plane
(slope-normal, 1 cm):
- **OPEN/BARE:** positive skew (+1.6); Gaussian core σ~70 mm; **power-law tails**, α~3.7 lower
  / ~2.6 upper (heavier above). Mode ~55-65 mm below the plane — but that plane is gen2 bare
  earth, so the offset is the **geoid datum (+67 mm)**, NOT an instrument bias.
- **FOREST:** near-symmetric, broad (σ~120-170 mm), **exponential tails** (λ~106 mm).

## Mechanisms (best supported first)
1. **Heavy UPPER tail = one-sided low-vegetation / crop-residue contamination of the ground
   class.** Pulses that miss true soil register on standing residue/stubble above it — adds
   height, never subtracts → one-sided upper tail on a clean core.
   - Hopkinson et al. (2005): ground *overestimated* +0.07 m (graminoids) to +0.15 m (scrub);
     "laser returns tend to be biased upwards into the more highly reflective foliage."
   - Su & Bork (2006): shrub/forest cover overestimates ground; sign is cover-dependent.
   - Bater & Coops (2009): vegetation positively biases DEM at fine cells (+4.5 cm).
   *Caveat:* the specific Nov corn/soybean-stubble instance is NOT in the peer-reviewed record;
   the low-veg mechanism is solid, this instance is our inference.
2. **Lower tail + (some) below-plane pull = timewalk / range-walk.** Baltsavias (1999):
   spreading the spot over inclined/rough terrain lengthens rise-time to threshold →
   range *increased* → elevation biased *low* (up to ~50 cm at 45°, 1000 m AGL). Detector
   range-walk (Li 2017/2018; Laconte 2019) similar. Biases weak/oblique returns LONG → heavier
   *lower* tail. NOTE: does NOT explain the upper tail, and our central offset is datum, so
   don't over-credit timewalk for the mode.
3. **FOREST symmetric/broad/exponential = sparse leaf-off ground returns diluted by
   random-height understory hits** (Reutebuch 2003: uncut 0.31±0.29 vs clearcut 0.16±0.23 m;
   Hodgson & Bresnahan 2004: ~26 cm deciduous vs ~17-19 cm open). Broadens the core both ways
   rather than adding a one-sided tail.

## Novelty / honesty
- **Power-law tail FORM is not in the lidar-accuracy literature.** Non-Gaussian/skewed/
  heavy-tailed is well documented (Zandbergen 2011 "strong evidence" of non-normality, outliers
  common; Höhle & Höhle 2009 robust NMAD/quantiles; ASPRS 2015 uses 95th-pct for vegetated
  ground *because* errors are skewed). But no source fits a power law or Student-t. **Present
  our power-law characterization as a novel finding, not attributed.**
- Third/fourth moments (skew/kurtosis) are generally NOT tabulated in ALS-DTM accuracy papers;
  they report mean bias + RMSE/σ. Large σ implies skew but is not the same measurement.

## Key citations
- Baltsavias 1999, *ISPRS J.* 54(2-3):199-214 (timewalk).
- Hopkinson et al. 2005, *IAPRS* XXXVI-8/W2:108-113 (upward veg bias, with numbers).
- Su & Bork 2006, *PE&RS* 72(11):1265-1274 (cover-dependent sign).
- Höhle & Höhle 2009, *ISPRS J.* 64(4):398-406 (robust stats for non-normal error).
- Zandbergen 2011, *IJRS* 32(2):409-430 (non-normality, outliers, no slope dependence).
- Hodgson & Bresnahan 2004, *PE&RS* 70(3):331-339; Reutebuch et al. 2003, *Can. J. Remote Sens.*
  29(5):527-535 (forest vs open ground-return spread).
- ASPRS 2015, *PE&RS* 81(3) (95th-pct vegetated accuracy).
- Bater & Coops 2009, *Computers & Geosciences* 35(2):289-300; Gardner 1992, *IEEE TGRS*
  30(5):1061-1072; Mallet & Bretar 2009, *ISPRS J.* 64(1):1-16; Wagner et al. 2006, *ISPRS J.*
  60(2):100-112; Aguilar & Mills 2008, *Photogramm. Record* 23(122):148-169; Li et al. 2017
  *Sensors* 17(10):2369 & 2018 18(4):1156; Laconte et al. 2019 IEEE ICRA; Spaete et al. 2011,
  *Remote Sensing Letters* 2(4):317-326.
