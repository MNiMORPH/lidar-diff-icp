# Is the across-track intercept tie in the literature?

Literature search run 2026-09-05 by a subagent with web access, in response to Andy's
question. **I have not read these papers myself.** What follows is the agent's report with
its verification status preserved; every DOI is given so any claim can be checked. Treat
section C as unverified.

---

## A. The short answer

**The regression is published. Taking its INTERCEPT as the vertical tie is not** — the agent
found no instance, across ISPRS, USGS, Crossref, Semantic Scholar and three vendor manuals.

## B. Verified against a retrieved source

### B1. The closest published thing does the same fit and throws away the half we keep

Sampath, A., Heidemann, H.K., Stensaas, G.L. (2016), "Geometric quality assessment of lidar
data based on swath overlap", *Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci.*
XLI-B1, 93–99, doi:10.5194/isprsarchives-XLI-B1-93-2016 — a joint USGS / ASPRS Airborne
Lidar Committee working group.

They fit exactly our line: discrepancy against distance from the centre of the overlap, and
call its slope the **Calibration Quality Line**. But:

> "In the flat regions, the magnitude of vertical bias increases from the center of the
> overlap. … these errors can be modelled as straight line **passing through the center of
> the overlap (where they are minimal)**."

So they **constrain the intercept to zero** and report the slope as a QA metric. We free the
intercept, report it, and discard the slope. Same machinery, opposite output. They also state
our motivation from the RMSD side — that a small overlap or centre-clustered sampling can
make the reported RMSD small regardless.

**On their own evidence we measure a nonzero value where they assume zero.** That makes ours
the less-constrained specification; it does not make it published.

### B2. The field's answer: the across-track term belongs to STRIPS, in a block adjustment

- **Crombaghs, M.J.E., Brügelmann, R., de Min, E.J. (2000)**, IAPRS 33(B3/1), 224–231 — a
  per-strip three-parameter model of {height offset, along-track tilt, **across-track tilt**}.
  That is `{k, –, c}` per *strip*, in a block least-squares adjustment. Nearest standard
  formulation to ours; quoted via Vosselman & Maas (2001), OEEPE Stockholm.
- **Glira, P., Pfeifer, N., Briese, C., Ressl, C. (2015)**, *ISPRS Annals* II-3/W5, 73–80,
  doi:10.5194/isprsannals-ii-3-w5-73-2015 — rigorous ICP strip adjustment, 12+6n parameters,
  with a per-strip roll `Δφᵢ` and a scanner angular-encoder zero-point `Δα`. Implemented in
  OPALS `ModuleStripAdjust`. Their warning applies directly to us:

  > "depending on the assembly of the sensors, the flight configuration, or the terrain
  > geometry, some of these parameters may be completely correlated and therefore not
  > estimable."

- **Kilian, Haala & Englich (1996)**; **Filin (2003)**, ISPRS XXXIV-3/W13 — scan-angle error
  as a system parameter, with non-linear "bending" at swath ends.

### B3. Acceptance testing — which governs 3DEP — models nothing at all

**USGS Lidar Base Specification 2022 rev. A**, "Interswath (Overlap) Consistency": a signed
difference raster summarised to RMSDz, with limits QL0 ≤0.04 m, QL1/QL2 ≤0.08, QL3 ≤0.16 m,
and the instruction that *"test areas should be located such that the full width of the
overlap is represented."*

That instruction is the spec conceding the extent-dependence we identified **and handling it
by sampling design rather than by modelling it.**

### B4. Operational practice already pairs the two terms

GeoCue (2015), TerraScan/TerraMatch calibration guide:

> "a dZ on an individual line basis. When evaluating this shift we generally also look at
> roll in the same manner … this roll error is generally referred to as a dR correction and
> is **always performed in conjunction with a dZ correction** on an individual line basis."

Qualitative, per line, not an intercept estimator — but it is the operational statement of
our position: the vertical bias is never estimated without the across-track term.

## C. NOT verified — the agent's own derivation, and the substantive criticism

**This is the part worth acting on, and it comes from no source.** For a small roll/encoder
error Δφ, the vertical displacement goes as `H·Δφ·tan θ`. So

    dh = k + H·Δφ_ref·tan θ_ref − H·Δφ_src·tan θ_src

Our model `dh = k + c·dtan`, with `dtan = tan θ_ref − tan θ_src`, is that expression **only
if Δφ_ref = Δφ_src**. It estimates a *common-mode* term. If the two lines genuinely differ,
the omitted regressor is correlated with `dtan` and **`k` inherits some of that bias** — so
the extent-independence we want from `k` holds only under the model we fitted.

Our own measurement is the evidence for the criticism: per-pair `c` is heterogeneous at
p ≈ 6e-64. The agent reads that as the diagnostic that a **per-pair** parameterisation is
misspecified, because pairwise `c` mixes two lines' terms.

**The proposed fix, small and inside our framework:** regress `dh` on `tan θ_ref` and
`tan θ_src` as **two separate covariates** plus an intercept. That is the reduced object-space
form of a per-strip roll, and one step from assembling into the block adjustment
Crombaghs/Glira describe. Caveat: within a single pair the two regressors may be
near-collinear — which is *why* the field solves it as a block over many pairs — so the
conditioning must be checked before trusting the split.

Also unverified: Latypov (2002)'s method description; Kumari et al. (2011); Skaloud & Lichti
(2006); the two Habib et al. (2010) papers (bibliographic records only). The agent named
Habib's IEEE TGRS 48(1) 221–236 paper on parallel-strip internal QC as the most likely place
to find a published statement that the discrepancy varies systematically across the overlap.

Not searched: ICESat / ATM / ICESat-2 crossover literature, where a roll-vs-cross-track term
is standard — the agent flagged it as the most plausible unexplored home for our exact
estimator.

## D. What this means for us

1. Our estimator is **not** unsupported — its regression is USGS/ASPRS-standard practice
   (B1) and its concern is the field's (B2, B4). What is unpublished is using the intercept
   as the tie.
2. It is a **pairwise** solution to what the literature treats as a **per-strip block**
   problem. That is a real structural difference, not a presentational one.
3. The two-covariate test (C) is cheap, uses data we already have, and is the next thing to
   run. It would tell us whether our `k` is carrying a differential-roll bias.
