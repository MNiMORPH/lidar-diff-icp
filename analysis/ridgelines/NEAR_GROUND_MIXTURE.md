# The near-ground return distribution: why Gaussian + exponential, and how to fit it

Adopted 2026-09-03. This records the reasoning so it sits with the code rather than in a
conversation, and so the next person can attack the argument rather than re-derive it.

## What a lidar return actually samples

A pulse records a return at height `h` only if it **reached** `h` without being fully
intercepted above, and then **intercepted** something there. The observed return-height
distribution is therefore the material density multiplied by the gap probability,

    P_gap(z) = exp(-G * LAI(z) / cos(theta))

(`RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md`). Returns are a *shadowed* sample of the vertical
structure, weighted toward whatever is highest. Ground is the extreme case: recorded only for
pulses that got all the way down.

**The bias is in the WEIGHT, not the LOCATION.** Occlusion controls how many ground returns
appear; it does not move where they sit. So the ground component's centre is the part of the
distribution least damaged by preferential sampling, while its mixing weight is essentially a
penetration measurement. Under dense cover the estimator gets noisier rather than drifting.

That is also why a MIXTURE beats a PERCENTILE. A percentile of the observed column moves as
soon as the component weights shift, and the weights are exactly what occlusion changes. The
component location does not move.

## The two components

**Ground: a narrow Gaussian.** Range measurement error plus the true surface's roughness
within the cell. Its width is the INSTRUMENTAL floor, sigma ~ 0.011 m, plus roughness.

**Above it: SHAPE UNKNOWN, and deliberately so.** Superseded 2026-09-03 -- see below.

### Why the exponential was wrong (superseded, kept visible)

The first version of this document specified an exponentially-modified Gaussian: Beer-Lambert
extinction convolved with measurement error. `tau` in that model is a MEAN FREE PATH, and

    1 / tau  =  G * a / cos(theta)

so its reciprocal measures leaf area density `a` -- the MacArthur-Horn gap-fraction
inversion, the same physics PyForestScan uses for plant area density above 2 m. That part
stands and is worth keeping.

Two things were wrong with it.

**Orientation.** Beer-Lambert decays from the TOP of the canopy DOWNWARD, not upward from
the ground. An exponential falling away from the ground is not what the physics predicts.

**Constant `a`.** The exponential holds only if leaf area density is uniform over the fitted
range. It is not: herb layer, shrubs, trunk space, crown. Stacking the control marks'
profiles on their SURVEYED ground shows it directly (mean normalised density, x1000):

       h mm   lowveg<=0.02   0.10-0.25   lowveg>0.25
       -100          17.67       23.22          8.39
          0         130.75       45.40         18.25
        +50         101.26       43.85         30.11
       +200           8.09       28.32         23.02
       +400           0.01       12.98         25.46
       +600           0.00        2.53         15.45

Bare ground falls 130.75 -> 0.01 over 400 mm. Under cover the profile is FLAT from +50 to
+500 mm and rises again near +400 -- a layered slab, not a decay. The EMG fits were already
saying so: `tau` stretched to a p90 of 883 mm trying to make an exponential look flat.

Also visible there: below lowveg 0.25 the mode sits at +10 mm above surveyed ground -- one
bin, i.e. AT the ground -- with ~40% of returns below it, as a distribution centred on the
ground should have. Above 0.25 the mode jumps to +70 mm and the mass below falls to 14%. The
dominant peak stops being the ground and becomes the cover, while a minority ground
component survives underneath: shots still get through.

### The model that replaces it: separate by SCALE, not by family

Leave `a(h)` unknown. The MacArthur-Horn inversion never needed a form for it, and the two
components are distinguishable by scale instead:

    ground      the ONLY sharp feature -- instrumental sigma ~11 mm, plus roughness
    vegetation  smooth in height, structured over hundreds of mm

Fit as EM where the vegetation M-step is a SMOOTH of its own responsibility-weighted counts
rather than a parameter update, with the smoother wide enough that the vegetation component
cannot represent a ground-width peak. Robust to `a(h)` doing anything, which is what a
statewide run over 55,296 tiles requires.

**A short window is not an arbitrary choice -- it linearises the inversion** (Andy). Over a
few tens of centimetres the gap probability is nearly constant, so `n(h) = a(h) P_gap(h)`
reduces to `n(h) ~ a(h)`: the near-ground return density IS the near-ground vegetation
structure, with no extinction correction to make. That is the reason to work near the ground,
and it is independent of any convenience.

### Measured: it cannot be fitted on 0.02 m histograms

Run on the control boxes, sweeping the smoother bandwidth:

    bw=0.10   sigma_g median 50.5 mm   w_g median 0.940   mu_g-surveyed +24.1 mm
    bw=0.15   sigma_g median 58.8 mm   w_g median 0.976   mu_g-surveyed +24.0 mm
    bw=0.25   sigma_g median 67.7 mm   w_g median 0.988   mu_g-surveyed +25.6 mm

`w_g` runs to 0.94-0.99 and `sigma_g` floats to 50-68 mm: the Gaussian SWALLOWS the
vegetation instead of separating from it. Scale separation requires the ground component to
be held at the sharp scale, and 11 mm is BELOW ONE 20 mm BIN. The stored histograms cannot
resolve the feature the method depends on, whatever the parameterisation. Re-acquiring the
marks at 1 cm, or as raw slope-normal distances, is a prerequisite rather than a refinement.

## Why ONE form, when a pooled measurement said otherwise

`ground_mixture_fit.py` measured tails on histograms POOLED per land-cover stratum and found
exponential in forest, power-law on open/bare. A mixture of exponentials with different rates
is heavier-tailed than any of its parts -- with Gamma-distributed rates it is exactly Pareto --
so pooling cells of differing extinction rate can manufacture a power law from exponentials.

`tail_form_pooling_test.py` tests this by estimating each cell's own rate, grouping cells by
it, and pooling only within a group. Elba, class-2 returns, origin at each cell's own modal
bin (no chosen height enters):

    open (cover <= 0.1)
        POOLED all cells      exp R2 0.9263   power R2 0.9430   power wins
        lambda group 1/5      exp R2 0.9521   power R2 0.9139   exp wins
      min_tail = 100
        POOLED all cells      exp R2 0.9512   power R2 0.9206   exp wins
        lambda group 1/5      exp R2 0.9846   power R2 0.8360   exp wins
    forest (cover >= 0.5)
        POOLED all cells      exp R2 0.8255   power R2 0.9420   power wins
        lambda group 5/5      exp R2 0.9506   power R2 0.8421   exp wins

Stratifying by rate reverses the verdict in BOTH covers. The power law is a pooling artifact;
one functional form is right, and no cover-dependence of FORM is needed.

Caveats: R2 on log counts ignores Poisson noise in the bins, so this discriminates without
quantifying; the first version of this test used ALL returns and gave power-law everywhere,
which was measuring canopy above the mode -- a different distribution and not the documented
claim; and the pooled FOREST result here disagrees with ground_mixture_fit's, whose window,
strata and 1 cm binning all differ, so the two are not strictly comparable and have not been
reconciled.

## Fitting it when the Gaussian's data are truncated

They are, in three different senses, and only one of them is a real problem.

**1. The upper flank is MIXED, not truncated.** Above the ground mode, ground returns and
vegetation returns overlap and cannot be separated return by return. This is
non-identifiability, and the geometry solves it: **vegetation cannot lie below the ground**, so
the LOWER flank of the observed near-ground distribution is essentially pure ground. Fit the
Gaussian's centre and width from the left side and let the exponential absorb the excess on the
right. The asymmetry that makes the problem hard is what makes it identifiable.

**2. Classification truncates, so do not fit to class-2.** A class-2 histogram has already been
cut by a classifier, in a way that varies with vendor and terrain, and it discards genuine
ground returns while keeping misclassified low vegetation (measured at the control marks: our
class-2 surface sits ~62 mm high in vegetation). Fitting ALL returns avoids importing a
classifier's decision into a physical model. The exponential component is then larger, which
is not a difficulty -- it is the thing being modelled.

**3. Missing ground returns ARE a real truncation, and fitting cannot fix it.** Where
penetration is poor the ground component is not thinned but absent: no pulse reached the
surface. That is missing-not-at-random and no likelihood repairs it. The correct response is to
report it, not to return a number -- a fitted `mu_g` from a near-zero-weight ground component
should REFUSE rather than be quoted, and the weight itself is the diagnostic that says so.

**Likelihood.** Fit on binned counts with a POISSON likelihood, not least-squares on the
counts or on their logs. Bin counts are Poisson; least-squares assumes constant variance and
log-least-squares is undefined at empty bins and badly biased at small ones -- exactly the bins
in the tail that carry the shape information. 1 cm bins, well below sigma_g ~ 0.1 m, keep the
discretisation far from the parameter being estimated while remaining vectorisable by EM across
every cell of a tile at once, which is what makes 55,296 tiles feasible.

## What is NOT settled

- The gaps through which ground is seen under dense cover are not a random subset of the cell:
  they are where vegetation is absent, which may correlate with microtopography.
- Attenuated ground returns under canopy have lower SNR, and range walk biases weak returns LOW
  (Baltsavias 1999, `BARE_GROUND_SKEW_LITERATURE.md`). That moves `mu_g` down for a reason
  unrelated to vegetation, and is not in this model.
- `mu_g` has not been validated against surveyed control. The 1,497 marks exist for exactly
  that test, and until it is run the mixture is a better-motivated estimator, not a verified one.
