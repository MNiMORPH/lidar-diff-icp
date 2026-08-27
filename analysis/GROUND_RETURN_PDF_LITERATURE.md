# Probability density functions for near-ground lidar returns and for lidar/DTM vertical error

Literature grounding, searched 2026-08-27. Every citation below was checked against
Crossref (DOI, title, venue, volume, pages, year) or read from the primary PDF; where a
claim comes from a secondary source rather than the primary text, it is labelled
**[secondary]**. Where a full text was paywalled and only metadata could be verified,
that is labelled **[metadata only]**.

Companion notes already in this repository, which this document extends rather than
repeats: [`ridgelines/RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md`](ridgelines/RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md)
and [`ridgelines/BARE_GROUND_SKEW_LITERATURE.md`](ridgelines/BARE_GROUND_SKEW_LITERATURE.md).

Throughout, each finding is tagged:
- **[DERIVED]** – follows from a stated physical model, with its assumptions named.
- **[FITTED]** – an empirical form chosen because it fit the data better.
- **[CONVENTIONAL]** – used because it is what the field uses, with no derivation and
  no fit comparison offered.

---

## Summary: what is actually established

**1. There is no established probability density function for the vertical distribution
of discrete lidar returns about a ground surface.** This is the central finding, and it
is a negative one. The Gaussian is *conditionally* derived for the full-waveform *echo*
(Wagner et al. 2006), but the condition is an assumption about the height distribution of
scatterers inside the footprint, and Chauve et al. (2007) state in print that this
assumption has been demonstrated only for large-footprint systems and that for
small-footprint systems "there is no assuming that the height distribution is Gaussian,
even over vegetated areas." Our fixed-σ Gaussian core is **[CONVENTIONAL]**, not derived.
Nobody is in a position to tell us the right family, because nobody has published one.

**2. The non-Gaussian echo families are the wrong shelf to shop from.** Generalised
Gaussian, lognormal, Burr and Nakagami are real and really used (Chauve et al. 2007;
Mallet et al. 2011), but they model the *time-domain shape of one echo* – an
instrument-response convolution – not the population distribution of return heights over
a site. Furthermore, the empirically fitted generalised-Gaussian shape parameter came out
slightly *flatter* than Gaussian rather than heavier-tailed (Chauve et al. 2007, Fig. 8),
so that family would not repair a tail misfit even if the transfer were legitimate.

**3. An exponential tail for the vegetation component is genuinely [DERIVED]**, from
Beer–Lambert–Bouguer extinction, but only under assumptions we should check rather than
inherit: foliage elements randomly positioned, statistically independent, and – the one
that matters for us – vertically *uniform* in areal density through the layer (He & Lyu
2023). Clumping does not break the exponential family; it rescales the rate through a
clumping index (He & Lyu 2023). A near-ground layer with a characteristic top – grass,
crop stubble, a shrub stratum – violates the vertical-uniformity assumption, and that is
a concrete, physical reason a single exponential would misfit and would misfit *worse
where there is more vegetation*.

**4. The DTM-error literature is emphatically non-parametric, and that is a deliberate
choice, not an oversight.** Non-normality of lidar vertical error is well established
(Höhle & Höhle 2009; Zandbergen 2011), and the ASPRS standard writes it into the
definitions. The field's response was to abandon parametric description entirely and
report median, NMAD and sample quantiles. No standard names a distributional family. So
there is no accepted parametric model against which a KS statistic is the natural test.

**5. Someone has fitted an explicit two-component mixture to lidar return elevations –
and hit our problem.** Kalinicheva et al. (2022) fit a Gamma + Gamma mixture by
Expectation–Conditional–Maximisation, and report that it "requires a meaningful
initialization of mixture parameters, which can be achieved by trial-and-error guided by
the resulting likelihood value." That is the same instability. But their split is
ground-plus-low-vegetation versus medium-and-high vegetation, which is *not* our split.
**The ground-versus-low-vegetation mixture we are fitting appears not to have been
published by anyone.**

**6. Our identifiability failure is a known, solved problem in statistics, not a lidar
problem.** A normal mixture with free component scale parameters has an unbounded
likelihood and no consistent unconstrained MLE (Day 1969; McLachlan & Peel 2000).
Hathaway (1985) gives the standard fix: constrain the ratio of the scale parameters. The
Gaussian absorbing the tail when σ is freed is the textbook symptom, and the textbook
remedy is a constraint, not a different family.

### What the next fit should be

In priority order. Options (a) and (b) address a *misspecification*; (c) and (d) address
the *identifiability*. They are independent and can be combined.

**(a) Replace the mixture with a convolution for the vegetation component.** Every
return, ground or vegetation, carries the ranging and registration error. The vegetation
component should therefore be `N(0,σ) ⊛ Exp(λ)` – the exponentially-modified Gaussian –
not a bare `Exp(λ)`:

    p(h) = (1-f)·N(0,σ) + f·EMG(σ, λ)

The bare exponential has a hard edge and its maximum density at h = 0 exactly where the
ground peak also lives, which forces the fit to trade the two components against each
other at the one place they are least distinguishable. This is a consequence of the
measurement model and I regard it as sound, but **I did not find it in the lidar
literature and it is not attributed** – the nearest published relative is Ingram et al.
(2026), who use a bare one-sided exponential likelihood for lidar in vegetation, without
convolving in the ranging error. Treat this as our reasoning, to be tested, not as an
inherited result.

**(b) Let the vegetation layer have a top.** If the near-ground layer has a
characteristic height, the vertical-uniformity assumption behind the exponential fails,
and a Gamma or Weibull tail is the natural generalisation – Gamma has lidar precedent
(Kalinicheva et al. 2022), Weibull has vertical-profile precedent (Section 2). Both nest
the exponential as a special case, so the exponential remains testable as the shape = 1
null rather than being discarded by assertion. Note that a Lomax/power-law tail failing
is *consistent* with this: Lomax is the wrong direction of departure, being heavier than
exponential, where a layer with a top is *lighter*.

**(c) To let σ vary without the peak eating the tail – constrain the ratio, do not free
the parameter.** Hathaway (1985) is the citable standard: bound σ_ground/σ_veg (or
σ_ground below a multiple of the fixed bare-ground value). This is a bound we would be
choosing, so it is a proposal for Andy, not a detail – and the bound's value should be
reported with the result.

**(d) Better: make σ measured rather than fitted, but per-site.** Glennie (2007) gives a
rigorous per-point error propagation from sensor geometry, and Gardner (1992) gives the
slope/beamwidth range-broadening term. Together they predict how σ *should* vary with
terrain and geometry from quantities we already have. Fixing σ per site from that model
keeps the parameter count unchanged, removes the identifiability failure by construction,
and – unlike a free σ – is falsifiable. This is the option I would try first if the
covariates are available.

**A caution on the diagnostic.** Given point 4 above, a KS test against a parametric
model is a strong test the field itself does not apply to this quantity. Rejection at
almost every site is what a KS test does to a two-parameter model at large n. The
misfit's *correlation with vegetation* is the informative signal here, not the KS
statistic's absolute level; that correlation is what points at the vegetation component
being misspecified, which is what (a) and (b) address.

---

## Q1. Ground-return / terrain-echo vertical distribution

### The Gaussian is conditional, and the condition is unverified for small footprints

The standard model decomposes a return waveform into a sum of Gaussians (Hofton et al.
2000; Wagner et al. 2006). Wagner et al. (2006) reported that more than 98% of observed
RIEGL LMS-Q560 waveforms could be fitted this way **[FITTED]**.

The *derivation* is conditional. Chauve et al. (2007, §3.3) state it precisely:

> "Each laser output pulse shape is assumed to be Gaussian, with a specific and calibrated
> width. The collected pulse is therefore a convolution between this Gaussian distribution
> and a 'surface' function, depending on the hit objects. It has been shown that if the
> vertical height distribution of the elements within the diffraction cone follows a
> Gaussian law, the reflected waveform can be approximated by a sum of Gaussians (Zwally
> et al., 2002)."

So: Gaussian scatterer-height distribution ⟹ Gaussian echo **[DERIVED]**. The antecedent
is the thing we care about, and Chauve et al. immediately disclaim it:

> "Moreover, the Gaussian height distribution of the targets has only be statistically
> shown for large-footprint lidar systems (Carabajal et al., 1999). For small-footprint
> systems, there is no assuming that the height distribution is Gaussian, even over
> vegetated areas. Therefore modelling full-waveform lidar data with a sum of Gaussian
> functions can be inaccurate."

This is the direct answer to the question as asked. **A Gaussian for the vertical spread
of returns about the ground is not justified from physics for small-footprint airborne
lidar; it is convention, and the one paper that checked the underlying assumption says so
in print.**

What *is* derived is how the width should scale. Gardner (1992) gives the
range-broadening of a laser altimeter echo by surface slope and beamwidth (the
`tan(slope)` term), and Glennie (2007) gives rigorous 3-D error propagation for kinematic
scanning lidar, producing a per-point covariance from sensor geometry **[DERIVED]**. These
constrain σ; they do not establish the *family*.

### The non-Gaussian echo families are real – and are about echo shape, not return heights

Verified as real, and used for the reasons stated:

- **Generalised Gaussian and lognormal** – Chauve et al. (2007). Their generalised
  Gaussian is `a·exp(−|x−μ|^α / 2σ²)`, with α = √2 the Gaussian, α = 1 the Laplace, α < √2
  peaked and α > √2 flattened **[FITTED]**.
- **Burr and Nakagami** – Mallet et al. (2011), reportedly generalised Gaussian for
  symmetric echoes, Burr for heavily right-skewed and Nakagami for softly skewed echoes
  **[FITTED]**, **[secondary]** (the Burr/Nakagami assignment is consistently reported by
  papers citing Mallet et al. 2011; the primary is paywalled and I verified only its
  metadata).

Three reasons not to import these:

1. **They describe a different quantity.** These are the shape of a single echo in the
   time domain – the transmitted pulse convolved with the target's range profile. Our
   histogram is a population of discrete return heights pooled over a site. The
   convolution structure and the sample space are not the same.
2. **The empirical result points the wrong way.** Chauve et al. (2007, Fig. 8) fitted α
   over more than 15 000 waveforms each on roofs, asphalt and dense vegetation and found
   all three means near 1.55, against √2 ≈ 1.414 for the Gaussian – that is, "the general
   shape of the backscattered echoes is close to a slightly flattened Gaussian." Flattened,
   not heavy-tailed. This family would not fix a tail misfit.
3. **They were reported to destabilise the fit, for our exact reason.** On the lognormal
   and generalised Gaussian, Chauve et al. (2007) write: "it will also increase the number
   of fits that do not converge, just like the Lognormal. It is due to the increasing
   number of degrees of freedom of the function." The lognormal diverged in ~1% of cases
   against ~0.01% for the Gaussian, and globally fit *worse*. Extra shape freedom bought
   instability there too.

### Verdict on Q1

**Not established.** No published PDF for the vertical distribution of discrete returns
about a ground surface. The Gaussian core is conventional; its width scaling is derived;
its family is not. Our fixed-σ Gaussian is as defensible as anything in the literature,
which is to say defensible but unsupported.

---

## Q2. Vegetation / understorey return distribution above ground

### The exponential is genuinely derived – under three assumptions

Beer–Lambert–Bouguer extinction gives gap probability `k(N) = exp(−γ·λ·A_s/cosθ)`
(He & Lyu 2023, Eq. 1) **[DERIVED]**. He & Lyu (2023) state the assumption explicitly: it
holds when "the plant elements are distributed completely randomly and are statistically
independent," arising from a Poisson process in which each crown independently covers a
proportion of previously uncovered ground.

Two refinements matter for us:

- **Clumping does not change the family.** He & Lyu (2023) show that non-random foliage
  arrangement is absorbed into a variable clumping index γ(i), with the exponential form
  retained. So a clumped near-ground layer still gives an exponential tail, with a
  different λ. This is good news: it means a poor exponential fit is *not* explained away
  by clumping.
- **Vertical uniformity is the assumption that bites.** The exponential in *height*
  requires the foliage area density to be vertically uniform through the layer, so that
  cumulative foliage area is linear in height. This is the assumption a near-ground layer
  with a characteristic top violates.

The canonical inversion from returns to foliage profile is MacArthur & Horn (1969),
extended to lidar waveforms by Ni-Meister et al. (2001), who derive exponential decay of
return density with canopy depth via `P_gap(z) = exp(−G·LAI(z)/cosθ)` **[DERIVED]**.

### Alternatives used for near-ground / stratified vegetation

- **Weibull** for vertical foliage and canopy-fuel-load profiles **[FITTED]** – used
  because it can represent a layer with a mode and a top, which the exponential cannot.
  Widely applied in forestry vertical-structure work; I did not find it derived from
  radiative transfer.
- **Gamma** for the lower stratum of lidar return elevations **[FITTED]** – Kalinicheva
  et al. (2022); see Q4.
- **Bare one-sided exponential likelihood** for lidar in vegetation **[CONVENTIONAL/FITTED]**
  – Ingram et al. (2026) use an exponential likelihood for observations whose error is
  systematically one-sided, reporting reduced vegetated-zone RMSE and mean absolute bias
  relative to a heteroscedastic Gaussian baseline.

### Verdict on Q2

**Exponential is justified [DERIVED]**, with named assumptions, one of which – vertical
uniformity of the near-ground layer – is questionable for our case and is a specific,
testable explanation for a vegetation-correlated misfit. Gamma and Weibull are the
established generalisations, both **[FITTED]** rather than derived, and both nest the
exponential.

---

## Q3. DTM/DEM vertical error distributions

### Non-normality is established

- **Höhle & Höhle (2009)** is the standard citation. They proposed replacing RMSE with
  robust measures – **median, normalised median absolute deviation (NMAD), and sample
  quantiles** – for DEMs from laser scanning and automated photogrammetry **[metadata only:
  DOI/venue verified via Crossref; full text paywalled, content confirmed via multiple
  secondary sources]**.
- **Zandbergen (2011)** found "strong evidence" that lidar vertical error is not normally
  distributed and that both major and minor outliers are common; of five land-cover types,
  only urban approximated a normal distribution **[secondary]**.
- **Aguilar et al. (2010)** modelled vertical error in lidar DEMs **[metadata only]**.
- **Bui & Glennie (2023)** estimate gridded DEM uncertainty as varying with terrain
  roughness and point density **[metadata only]** – relevant precedent for a σ that varies
  with roughness.

### The standards make non-normality official – and decline to name a family

ASPRS (2015) writes it into the definitions, verbatim from the standard:

> "vegetated vertical accuracy (VVA) – An estimate of the vertical accuracy, based on the
> 95th percentile, in vegetated terrain where errors do not necessarily approximate a
> normal distribution."

> "non-vegetated vertical accuracy (NVA) – The vertical accuracy at the 95% confidence
> level in non-vegetated open terrain, where errors should approximate a normal
> distribution."

And the rationale:

> "different methods are used in non-vegetated terrain (where errors typically follow a
> normal distribution suitable for RMSE statistical analyses) and vegetated terrain (where
> errors do not necessarily follow a normal distribution). When errors cannot be
> represented by a normal distribution, the 95th percentile value more fairly estimates
> accuracy at a 95% confidence level."

The standard asks for higher moments to be *reported* but specifies no distribution:

> "ASPRS encourages standard deviation, mean error, skew, kurtosis and RMSE to all be
> computed in error analyses in order to more fully evaluate the magnitude and distribution
> of the estimated error."

Percentiles are taken on absolute errors: "For accuracy testing, percentile calculations
are based on the absolute values of the errors, as it is the magnitude of the errors, not
the sign that is of concern." The USGS Lidar Base Specification (current release LBS 2025
rev. A; historically USGS Techniques and Methods 11-B4) adopts the ASPRS accuracy
standards rather than defining its own error model.

### Recommended estimators, consolidated

Median (bias); NMAD (spread); sample quantiles, in particular the 68.3% and 95%
quantiles of the error and the 95th percentile of |error| for vegetated terrain.

### Named parametric families actually proposed

Very few, and none standard. Per the companion note in this repository, the only fitted
heavy-tail forms found are **Laplace / double-exponential** (Zandbergen 2011; Hejmanowska
& Kay 2011), both on pooled mixed-cover data. No standard proposes a family.

### Verdict on Q3

**Non-normality: established.** **Recommended estimators: established and robust /
non-parametric.** **A distributional family: deliberately not established.** The field
chose distribution-free description over parametric modelling. This is why our KS-based
rejection has no published benchmark to be compared against – and it means "the mixture
fails KS" is a weaker indictment than it looks.

---

## Q4. Mixture models for ground vs non-ground returns

### Yes – one, and it reports our instability

**Kalinicheva et al. (2022)** model the elevation of lidar points as a mixture of two
Gamma distributions **[FITTED]**:

> "By plotting the elevation histograms of all points [...] we observe that this empirical
> distribution follows a mixture model of two Gamma distributions. Moreover, we can easily
> interpret its components: the low elevation density peak corresponds to bare soil and low
> vegetation, while the long-tailed high elevation distribution corresponds to medium and
> high vegetation."

Components: Gamma(α_G, β_G) and Gamma(α_NG, β_NG), plus mixture weights ρ_G, ρ_NG, fitted
by Expectation–Conditional–Maximisation with an inner Newton–Raphson step. In their
Figure 6 the fitted weights are 0.55 and 0.45.

They report the identifiability problem, in the same terms as ours:

> "This procedure is entirely unsupervised but requires a meaningful initialization of
> mixture parameters, which can be achieved by trial-and-error guided by the resulting
> likelihood value."

and in their limitations:

> "the ECM algorithm for elevation modelling requires a manual initialization step.
> However, the parameters of the Gamma distributions are intuitive as they relate to the
> moment of the distribution (mean height, deviation) and can be approximated by a
> knowledgeable operator."

Their mitigation is therefore *informative initialisation from physically interpretable
parameters* – not a constraint and not a different family.

**Two important caveats.** (1) Their split is wrong for us: component 1 lumps ground
*with* low vegetation, and component 2 is medium-and-high vegetation. They are separating
the layer we are trying to decompose. (2) Gamma has support on h > 0 only, so it cannot
represent a peak symmetric about zero. Their z is a height above ground; ours is a signed
residual about a fitted surface. A Gamma tail is transferable to our vegetation component;
a Gamma is not transferable to our ground peak.

### Related but not the same thing

- **Liu et al. (2025)** use a Gaussian mixture model for ground filtering in dense
  vegetation, but in *feature space* – elevation residual, geometric features, intensity
  and green leaf index – with Mahalanobis distance for refinement. This is GMM as a
  clustering tool, not a distributional model of return heights. It does not answer the
  question.
- **Kraus & Pfeifer (1998)** is the classic ground/vegetation separation and is worth
  reading for Q5, but it is *not* a mixture model. It is robust iterative interpolation
  with an asymmetric weight function. Kraus & Rieger (1999) describe it: the algorithm
  "estimates the skewness of the error distribution of the laser scanner data in forests
  and assigns small weights to those points that show large positive errors during the
  interpolation with filtering," and it switches between "a skew error distribution
  function, i.e. forested areas, and regions with a symmetric error distribution, i.e.
  non-forested areas." So the field's canonical treatment of exactly our asymmetry is an
  **M-estimator with a one-sided influence function**, and an explicit cover-dependent
  switch between skewed and symmetric – **[CONVENTIONAL]**, algorithmic rather than a
  fitted PDF.

### Verdict on Q4

**Partly established.** One published two-component mixture on lidar return elevations,
Gamma + Gamma, with the same initialisation fragility. **Nobody has published the
ground-versus-low-vegetation mixture we are fitting.** The field's standard tool for that
separation is a one-sided robust weight function, not a mixture.

---

## Q5. Practical recommendation

### Diagnosis first

Our failure has two separable parts, and they have different literatures.

**Misspecification.** Grounded in Q1–Q2: the ground family is unsupported by anyone
(Q1), while the vegetation family is derived but under an assumption a near-ground layer
plausibly violates (Q2). The misfit correlating with vegetation amount points at the
vegetation component, which is exactly where the derivation has a named, checkable weak
assumption.

**Identifiability.** This is not a lidar problem and should not be researched as one. For
a normal mixture with free component scale parameters the likelihood is unbounded and the
unconstrained MLE is not consistent (Day 1969; McLachlan & Peel 2000). A component whose
variance is free will chase and absorb structure. Our "freeing σ makes f unidentifiable"
is the textbook symptom.

### Recommendations

**(1) Convolve, do not mix, the measurement error into the vegetation component.**
`p(h) = (1-f)·N(0,σ) + f·EMG(σ, λ)`, where EMG is `N(0,σ) ⊛ Exp(λ)`. Justification: every
return carries the ranging error, so the vegetation returns are displaced-*and*-noisy, not
displaced-instead-of-noisy. This also removes the bare exponential's density spike at
h = 0, which currently sits on top of the ground peak where the two components are least
separable. **This is our reasoning, not an inherited result** – I did not find it in the
lidar literature; Ingram et al. (2026) is the nearest published relative and uses a bare
exponential. It costs no extra parameters.

**(2) Give the vegetation layer a top: Gamma or Weibull instead of Exp.** Both nest the
exponential (shape = 1), so this is a likelihood-ratio test against our current model
rather than a replacement by assertion. Precedent: Gamma in Kalinicheva et al. (2022);
Weibull in vertical foliage/fuel-profile work. Note this is consistent with the Lomax
result: a layer with a top is *lighter*-tailed than exponential, and Lomax is heavier, so
Lomax pinning at its bound is evidence pointing the same direction.

**(3) If σ must vary, constrain it – do not free it.** Hathaway (1985) is the citable
standard for constrained ML on normal mixtures: bound the ratio of scale parameters. In
our case, bound σ_ground relative to the measured bare-ground value. *The value of that
bound is a choice, not a detail* – it should be Andy's, and it should be printed with
every result that depends on it.

**(4) Preferred: make σ site-specific but measured, not fitted.** Glennie (2007) gives
per-point 3-D error propagation from sensor geometry; Gardner (1992) gives the
slope/beamwidth broadening term; Bui & Glennie (2023) is precedent for uncertainty varying
with roughness and point density. Predicting σ per site from these and *fixing* it removes
the identifiability failure by construction, keeps the parameter count at two, and is
falsifiable in a way a free σ is not. This is the option I would try first if the
geometry covariates are to hand.

**(5) Reconsider the diagnostic.** Per Q3, the field does not fit parametric families to
this quantity and has no benchmark KS level. At our sample sizes a KS test will reject a
two-parameter model almost regardless. The vegetation-correlation of the misfit is the
real signal and should be the headline diagnostic; a KS statistic reported without the
sample size behind it will mislead.

### What the literature does not answer

- No published PDF for discrete-return heights about a ground surface (Q1).
- No published ground-versus-low-vegetation mixture (Q4).
- No published treatment of the ground-peak-absorbs-tail failure *in a lidar context*
  (Q4/Q5) – the statistics literature has it, the lidar literature does not.
- Nothing establishing how σ should vary with *vegetation* specifically, as opposed to
  with slope, roughness and point density.

If the EMG-plus-Gamma direction works, points 1, 2 and 4 above are, as far as this search
can tell, unpublished – the same standing as the power-law bare-ground tail already noted
in the companion documents.

---

## References

All DOIs verified against Crossref on 2026-08-27 unless noted.

**Full-waveform echo modelling**
- Chauve, A., Mallet, C., Bretar, F., Durrieu, S., Pierrot-Deseilligny, M., & Puech, W.
  (2007). Processing full-waveform lidar data: modelling raw signals. *IAPRS* XXXVI
  (Part 3/W52), 102–107. ISPRS Workshop Laser Scanning 2007 / SilviLaser 2007, Espoo.
  https://www.isprs.org/proceedings/xxxvi/3-w52/final_papers/Chauve_2007.pdf
  **[read in full]**
- Wagner, W., Ullrich, A., Ducic, V., Melzer, T., & Studnicka, N. (2006). Gaussian
  decomposition and calibration of a novel small-footprint full-waveform digitising
  airborne laser scanner. *ISPRS J. Photogramm. Remote Sens.* 60(2), 100–112.
  doi:10.1016/j.isprsjprs.2005.12.001
- Hofton, M. A., Minster, J. B., & Blair, J. B. (2000). Decomposition of laser altimeter
  waveforms. *IEEE Trans. Geosci. Remote Sens.* 38(4), 1989–1996. doi:10.1109/36.851780
- Mallet, C., & Bretar, F. (2009). Full-waveform topographic lidar: state-of-the-art.
  *ISPRS J. Photogramm. Remote Sens.* 64(1), 1–16. doi:10.1016/j.isprsjprs.2008.09.007
- Mallet, C., Bretar, F., Roux, M., Soergel, U., & Heipke, C. (2011). Relevance assessment
  of full-waveform lidar data for urban area classification. *ISPRS J. Photogramm. Remote
  Sens.* 66(6), S71–S84. doi:10.1016/j.isprsjprs.2011.09.008 **[metadata only]**
- Carabajal, C., Harding, D., Luthcke, S., Fong, W., Rowton, S., & Frawley, J. (1999).
  Processing of Shuttle Laser Altimeter range and return pulse data. ISPRS Archives
  workshop paper – **not in Crossref**; cited here as it is cited by Chauve et al. (2007).
- Zwally, H. J., Schutz, B., Abdalati, W., Abshire, J., et al. (2002). ICESat's laser
  measurements of polar ice, atmosphere, ocean, and land. *J. Geodynamics* 34(3–4),
  405–445. doi:10.1016/S0264-3707(02)00042-X. Cited via Chauve et al. (2007).

**Range/geometry error models**
- Gardner, C. S. (1992). Ranging performance of satellite laser altimeters. *IEEE Trans.
  Geosci. Remote Sens.* 30(5), 1061–1072. doi:10.1109/36.175341
- Glennie, C. (2007). Rigorous 3D error analysis of kinematic scanning LIDAR systems.
  *J. Applied Geodesy* 1(3). doi:10.1515/jag.2007.017

**Canopy extinction / vertical profiles**
- MacArthur, R. H., & Horn, H. S. (1969). Foliage profile by vertical measurements.
  *Ecology* 50(5), 802–804. doi:10.2307/1933693
- Ni-Meister, W., Jupp, D. L. B., & Dubayah, R. (2001). Modeling lidar waveforms in
  heterogeneous and discrete canopies. *IEEE Trans. Geosci. Remote Sens.* 39(9),
  1943–1958. doi:10.1109/36.951085
- He, Q., & Lyu, D. (2023). A modified Beer–Lambert–Bouguer law for nonrandom
  distributions and its application in gap probability calculations for heterogeneous
  canopies. *J. Advances in Modeling Earth Systems* 15(11), e2022MS003281.
  doi:10.1029/2022MS003281 **[read via full-text fetch]**

**DTM/DEM vertical error**
- Höhle, J., & Höhle, M. (2009). Accuracy assessment of digital elevation models by means
  of robust statistical methods. *ISPRS J. Photogramm. Remote Sens.* 64(4), 398–406.
  doi:10.1016/j.isprsjprs.2009.02.003 **[metadata only]**
- Zandbergen, P. A. (2011). Characterizing the error distribution of lidar elevation data
  for North Carolina. *Int. J. Remote Sens.* 32(2), 409–430.
  doi:10.1080/01431160903474939 **[metadata only]**
- Aguilar, F. J., Mills, J. P., Delgado, J., & Aguilar, M. A. (2010). Modelling vertical
  error in LiDAR-derived digital elevation models. *ISPRS J. Photogramm. Remote Sens.*
  65(1), 103–110. doi:10.1016/j.isprsjprs.2009.09.003 **[metadata only]**
- Bui, L. K., & Glennie, C. L. (2023). Estimation of lidar-based gridded DEM uncertainty
  with varying terrain roughness and point density. *ISPRS Open J. Photogramm. Remote
  Sens.* 7, 100028. doi:10.1016/j.ophoto.2022.100028 **[metadata only]**
- ASPRS (2015). *Positional Accuracy Standards for Digital Geospatial Data*, Edition 1,
  Version 1.0. *Photogramm. Eng. Remote Sens.* 81(3), A1–A26.
  https://florida.asprs.org/images/documents/ASPRS_Positional_Accuracy_Standards_Edition1_Version100_November2014.pdf
  **[read in full; all quotations verbatim]**. Edition 2 (2023) also exists.
- U.S. Geological Survey. *Lidar Base Specification*, LBS 2025 rev. A (online);
  historically USGS Techniques and Methods 11-B4, doi:10.3133/tm11b4.

**Mixture models and ground filtering**
- Kalinicheva, E., Landrieu, L., Mallet, C., & Chehata, N. (2022). Predicting vegetation
  stratum occupancy from airborne LiDAR data with deep learning. *Int. J. Applied Earth
  Observation and Geoinformation* 112, 102863. doi:10.1016/j.jag.2022.102863
  **[read in full via arXiv:2201.08051; quotations verbatim]**
- Kalinicheva, E., Landrieu, L., Mallet, C., & Chehata, N. (2022). Multi-layer modeling of
  dense vegetation from aerial LiDAR scans. *CVPR Workshops 2022*, 1341–1350.
  doi:10.1109/cvprw56347.2022.00140 **[read in full via arXiv:2204.11620]**
- Liu, C., Wang, H., Feng, B., Wang, C., Lei, X., & Chang, J. (2025). Integrating elevation
  frequency histogram and
  multi-feature Gaussian mixture model for ground filtering of UAV LiDAR point clouds in
  densely vegetated areas. *Remote Sensing* 17(18), 3261. doi:10.3390/rs17183261
  **[abstract only]**
- Kraus, K., & Pfeifer, N. (1998). Determination of terrain models in wooded areas with
  airborne laser scanner data. *ISPRS J. Photogramm. Remote Sens.* 53(4), 193–203.
  doi:10.1016/S0924-2716(98)00009-4 **[metadata only]**
- Kraus, K., & Rieger, W. (1999). Processing of laser scanning data for wooded areas.
  In D. Fritsch & R. Spiller (Eds.), *Photogrammetric Week '99*, 221–231. Wichmann.
  https://phowo.ifp.uni-stuttgart.de/publications/phowo99/kraus.pdf **[read in full]**
- Ingram, B., Paredes, R., Díaz, J., Besoaín, F., & Baettig, R. (2026). Spatial
  multi-sensor fusion with heterogeneous error characteristics. *Applied Sciences* 16(13),
  6294. doi:10.3390/app16136294 **[abstract only]**

**Mixture identifiability (statistics)**
- Day, N. E. (1969). Estimating the components of a mixture of normal distributions.
  *Biometrika* 56(3), 463–474. doi:10.1093/biomet/56.3.463
- Hathaway, R. J. (1985). A constrained formulation of maximum-likelihood estimation for
  normal mixture distributions. *Annals of Statistics* 13(2).
  doi:10.1214/aos/1176349557
- McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley Series in Probability
  and Statistics. doi:10.1002/0471721182
- Young, D. S., Chen, X., Hewage, D. C., & Nilo-Poyanco, R. (2019). Finite
  mixture-of-gamma distributions: estimation, inference, and model-based clustering.
  *Advances in Data Analysis and Classification* 13(4), 1053–1082.
  doi:10.1007/s11634-019-00361-y

---

## Method note

Citations were verified by querying the Crossref REST API for each work
(`api.crossref.org/works?query.bibliographic=...`) and checking author, title, venue,
volume, pages and year; primary PDFs were read where openly available. Items marked
**[metadata only]** had DOI, title, venue and pagination verified but the full text was
behind a paywall, so claims about their contents are attributed to secondary sources and
labelled as such. ScienceDirect and MDPI article pages refused scripted access (HTTP 403),
which is why several ISPRS Journal papers are metadata-only here.
