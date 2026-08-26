# gen1's vertical datum from routed sediment continuity

**Question.** How large a UNIFORM vertical raise of gen1 does the routed hillslope
sediment budget demand, and how does that compare with the +53.6 ± 13.0 mm the surveyed
control marks give? Continuity touches no survey data, so it is a genuinely independent
instrument.

**Answer, one line.** *Scene-integrated, at the edge of 1 σ, on the canopy-corrected DoD
with the correlated error envelope and the floodplain excluded:* **−24.5 ± 15.2 mm**
(S1 below). Two other equally defensible scene aggregates give **−5.8 ± 9.4 mm** (S3) and
**−2.3 ± 23.8 mm** (S2). All three are at or below zero: **the routed budget does not
demand that gen1 be raised at all.**

**And the caveat that decides how to read that.** Sediment continuity is a ONE-SIDED
constraint — `V_acc ≤ +z·σ`, and raising gen1 only ever makes it *more* satisfied. So the
strict physics yields a LOWER BOUND on delta and nothing more. As a bound, S1 says
`delta ≥ −24.5 mm`, and **+53.6 mm satisfies it: there is no conflict.** The numbers
become an *estimate* of delta only under the extra assumption that the true hillslope
budget is exactly balanced (zero net export), and that assumption is what puts them in
tension with the marks. Both readings are given below; the choice between them is a
geomorphic judgement, not a statistical one.

Producer: `analysis/mass_balance/datum_from_mass_balance.py` (this task's script);
comparison arithmetic `analysis/mass_balance/compare_to_control.py`. Raw run logs under
`analysis/mass_balance/out/`. Everything below is pasted from those runs.

---

## What was run, and what was rebuilt first

Both prior conditions in the brief were real and both are fixed here.

1. **The shipped `V_acc_*.npy` products are stale** — written 13:34–13:35 on 2026-08-26,
   before the extent-invariant swath-tie rebuild at 17:16–17:39. They are not used.
   `lod_cover_q2.npy` was equally stale (13:29) and was rebuilt first with the project's
   own producer, `analysis/ridgelines/lod_cover_q2.py`, which reports
   `dod.npy (uncorrected) 341,239 cells / 159,588 stable / sigma 51.2 mm / LoD med 98.3 mm`
   and `dod_cover_q2.npy 341,174 / 159,224 / 59.8 / 116.1`.
2. **The prior runs violated the standing floodplain rule.** `floodplain_mask.npy` covers
   **38.35 %** of the grid and its cells are now removed from the flow graph, so the
   budget is a hillslope budget: hillslope flow into a floodplain cell EXITS the network
   rather than accumulating in it. The evaluated population is **209,451 cells =
   5.236 km² = 58.9 % of the grid** (evaluable 214,562 of 217,609 routed = 98.60 %; the
   difference between evaluable and evaluated is cells with no known DoD).
3. **The prior runs also used the wrong error envelope** (`analysis/mass_balance/elba.py:28`
   passes no `corr_sill`/`corr_range`). Everything here uses the correlated N_eff envelope
   with a variogram fitted on the REBUILT DoD, and reports the independent envelope beside
   it as a sensitivity.

The variogram, fitted on the 80 % of cells with `|dod| ≤ lod` exactly as
`catchment-dod-balance/scripts/validate_site.py` does it:

    nugget sd 0.0445 m, sill sd 0.0379 m, range 502 m (total 0.0584 m)

against the pipeline's independently computed `stable_1sigma_m` of 0.0516 m — 13 % apart
on the canopy-corrected DoD. On the uncorrected `dod.npy` the same fit gives
`nugget sd 0.0368, sill sd 0.0337, range 588 m, total 0.0500 m`, 3 % from 0.0516.
The range is a per-site quantity (Elba 502–588 m here; *UNVERIFIED — Carlton 166 m, quoted from
`catchment-dod-balance/docs/validation.md`, not refitted in this task*) and is not
hard-coded anywhere.

Off-map boundary terms come from `catchment_dod_balance.offmap.offmap_terms`, which
reproduces the package's published Elba run exactly: **4444 boundary crossings, 45 off-map
trunks, hillslope length 188 m**.

## The sweep is exact, not sampled — and that is checked, not assumed

A uniform raise `delta` enters the budget linearly:

    V_acc(delta) = V_acc(0) − delta · area · N_up ,        sigma_Vacc unchanged.

The run asserts all three parts of that rather than trusting the algebra. From
`out/run_q2_corr_z1.txt`:

    CHECK  sigma_Vacc invariant under a 50 mm uniform shift: max |change| 0.000e+00 m^3
    CHECK  V_acc shift equals -delta*area*N_up: max |residual| 3.638e-12 m^3
    CHECK  contaminated / known / N_up unchanged by the shift: True
    CHECK  cached flow graph vs weighted_accumulation: max |diff| 0.000e+00 cells

So the per-cell **continuity surplus per unit upstream area**

    S(c) = ( V_acc(c) − V_in_acc(c) − z·sigma_Vacc(c) ) / ( area · N_up(c) )     [m]

is the whole curve: cell `c` stops being flagged as soon as `delta ≥ S(c)`, and

    delta leaving a fraction f of cells flagged  =  quantile(S, 1−f)
                                                =  −balancing_offset(keep = 1−f)["offset"].

That identity is verified in-run against the shipped API rather than claimed:

    CHECK  balancing_offset(keep=0.99) = +98.666 mm vs quantile(S,0.99) = +98.717 mm;
           coverage 98.6%, binding cell (667, 35)

**Prior art reused, not forked.** `dinf_proportions`, `mass_balance`, `balancing_offset`,
`offmap_terms` and `correlated_variance` are all called from
`~/projects/catchment-dod-balance`; the variogram is this repo's
`lidar_diff_icp.variogram`. The only new code is the delta sweep, the null calibration and
the block bootstrap.

---

## THE HEADLINE — three NAMED scene-integrated aggregates

All at `z = 1.0` (the edge of 1 σ, Andy 2026-08-26 — `mass_balance` ships 1.96), on
`dod_cover_q2.npy` with `lod_cover_q2.npy`, correlated envelope, floodplain removed from
the flow graph, over the 209,451-cell / 5.236 km² evaluated hillslope population.
Sign convention throughout: **positive = a constant to ADD to gen1**, which lowers the DoD.

### S1 — scene net volume brought to its own +1 σ  → **−24.5 ± 15.2 mm**

*The statistic:* the scene net volume, `sum over the evaluated hillslope cells of
DoD·area` — equivalently the sum of `V_acc` over the network's outlets — set equal to
`+1 σ` of that same sum under the fitted correlated error model (nugget over every cell,
sill over `N_eff`, the Rolstad/Hugonnet convention `mass_balance` uses internally).

      scene net volume, sum over evaluated hillslope cells of DoD*area, m^3  -51,014
                                      its correlated 1-sigma envelope, m^3   77,171
                                                        N_eff used (cells)      6.6
                       delta putting the scene net volume at +1 sigma (mm)    -24.5
                               the envelope's own share of that delta (mm)    +14.7
                                     the mean-DoD share of that delta (mm)     -9.7

Note what the decomposition says: the delta is the scene-mean DoD (**−9.74 mm**) minus its
own correlated 1 σ (**14.70 mm**). With `N_eff = 6.6`, the fitted 502 m correlation length
treats the whole 5.2 km² hillslope as ~7 independent patches, and that — not the data — is
what makes this criterion weak.

*Its uncertainty, named:* σ of the scene-mean DoD under the fitted correlated model,
**14.70 mm** (the same 14.7 mm above, since the criterion IS mean − 1 σ), combined in
quadrature with the half-range of the variogram-refit sensitivity, **3.90 mm** → **15.2 mm**.
The refit sensitivity, each fit carried all the way through to the delta:

              fit  nugget_sd_m  sill_sd_m  range_m  d_scene_mm
          as used       0.0445     0.0379      502       -24.5
    max_lag 400 m       0.0404     0.0385      300       -18.7
    max_lag 800 m       0.0455     0.0378      572       -26.5
           seed 1       0.0398     0.0434      414       -23.7
         2x pairs       0.0366     0.0454      365       -22.6

A purely empirical alternative to the model σ is the block bootstrap of the scene-mean
DoD. It has **not plateaued** by the largest block tried, so read it as a lower bound:

      L_m  mean_block_dod_mm  se_mean_mm  n_blocks
       25              -9.08        1.12      9722
       50              -9.75        1.50      2742
      100             -10.89        2.39       816
      200              -9.07        3.69       233
      400             -12.64        5.51        63
      800             -17.22        9.26        20

At `L = 800 m` (above the fitted 502 m range) it is 9.26 mm from 20 blocks, against the
model's 14.70 mm. Effective n here is 20 spatial blocks, or `N_eff = 6.6` correlation
patches — not 209,451 cells.

### S2 — total excess deposition volume matched to the fitted error model → **−2.3 ± 23.8 mm**

*The statistic:* `T(delta) = sum over evaluable cells of max(V_acc(delta) − 1·sigma_Vacc, 0)`,
in m³ — the total volume of deposition the budget cannot source. It is compared not to a
chosen threshold but to what the FITTED ERROR MODEL ALONE produces: 100 realisations of
white noise at each cell's own `perror`, plus a spherical Gaussian random field at the
fitted sill and range, plus one coherent draw on the unaccounted-area term exactly as the
envelope adds it, each routed through the same flow graph with no real deposition present.

                                 T observed at delta = 0 (m^3)  +331,036
    T from the error model alone, mean over realisations (m^3)  +357,805
                                its SD over realisations (m^3)  +404,098
                      delta at which T falls to that mean (mm)      -2.3
                             at mean - 1 SD of the null T (mm)   +3060.2
                             at mean + 1 SD of the null T (mm)     -26.1

*Its uncertainty, named:* the one-sided half-width from the delta at the null mean to the
delta at null-mean + 1 SD, **23.8 mm**. The other side is unbounded (the +3060.2 mm entry
is the blunder cell) because the null SD exceeds the null mean. **This aggregate is too
noisy at z = 1 to constrain delta on its own** and is reported for completeness.

### S3 — flagged-cell fraction matched to the fitted error model → **−5.8 ± 9.4 mm**

*The statistic:* the fraction of evaluable cells flagged as unphysical deposition, matched
to the fraction the same routed null realisations produce.

    flag rate produced by error alone, mean over realisations         0.1538
                                     its SD over realisations         0.0354
                     nominal one-sided rate implied by z=1.00         0.1587
       delta at which the observed flag rate falls to it (mm)           -5.8
                   +/- from the SD of the null flag rate (mm)  -13.1 .. +3.4
            block-bootstrap SE of that quantile, L=400 m (mm)            4.5
                   +/- from the SE of the mean null rate (mm)   -6.6 .. -5.0
                              observed flag rate at delta = 0         0.1308

The measured null rate 0.1538 sits within its own SD of the nominal one-sided 0.1587 that
`z = 1` implies — i.e. **the correlated envelope is well calibrated**, which is worth
having as a result in its own right.

*Its uncertainty, named:* half-width of the one-draw SD interval on the null flag rate,
**8.25 mm** (the observed tile is one realisation of its own error field), in quadrature
with the block-bootstrap SE of the quantile at `L = 400 m`, **4.50 mm** → **9.4 mm**.

**Why S1 is the number I would quote** and S2/S3 the support: S1 needs no null model and
no auxiliary assumption about how many cells "should" flag. It is the budget itself
against its own envelope. S3 is the tightest, but it depends on the realism of the
simulated error field; S2 is honest but uninformative at z = 1.

---

## Envelope sensitivity — this choice moves the answer

Same DoD, same population, `z = 1`, only the error model changes:

| aggregate | correlated (N_eff) | independent |
|---|---|---|
| S1, scene net volume at +1 σ | **−24.5 mm** (envelope 77,171 m³, N_eff 6.6) | −10.0 mm (envelope 1,418 m³, N_eff 209,451) |
| S2, excess volume matched to null | −2.3 mm (T obs +331,036 vs null +357,805) | **+36.3 mm** (T obs +794,112 vs null +234,313) |
| S3, flag rate matched to null | −5.8 mm (null 0.1538 ± 0.0354) | +3.3 mm (null 0.1839 ± 0.0200) |

The independent envelope is a documented LOWER BOUND on the true error, so it manufactures
apparent surplus and therefore a larger apparent raise — S2 moves by 39 mm on this choice
alone. **The correlated envelope is the correct one**; the independent column is here only
so the sensitivity is visible.

## Other sensitivities

| variant | S1 | S3 | scene mean DoD |
|---|---|---|---|
| canopy-corrected `dod_cover_q2`, floodplain out (headline) | −24.5 mm | −5.8 mm | −9.74 mm |
| uncorrected `dod.npy`, floodplain out | −18.3 mm | −2.3 mm | −2.91 mm |
| canopy-corrected, floodplain KEPT (rule-violating) | −23.7 mm | +3.0 mm | −12.13 mm |

The canopy correction moves the scene mean DoD by 6.8 mm and S1 by 6.2 mm; the floodplain
rule moves S1 by 0.8 mm and S3 by 8.8 mm, and changes the population from 209,451 to
338,894 cells.

---

## Supporting: the worst-cell reading and the full curve

Kept because it is instructive about how blunder-sensitive the strict reading is, NOT as
the answer. `keep = 1.0` is set by a single cell.

      keep  delta_mm  delta_api_mm  se_mm  binding_cell  cover
    1.0000   +3060.2       +3060.2  239.5     (73, 203)  0.986
    0.9990    +497.6        +497.6  147.1    (619, 235)  0.986
    0.9900     +98.7         +98.7   11.1     (667, 35)  0.986
    0.9750     +56.7         +56.6    7.3    (522, 333)  0.986
    0.9500     +31.2         +31.1    5.6      (90, 74)  0.986
    0.9000      +8.6          +8.6    4.7     (401, 82)  0.986

`delta_api_mm` is `−balancing_offset(keep)["offset"]` from the shipped API; `cover` is its
`evaluable_fraction`. The 0.986 coverage is the same at every `keep`, so no value here is
resting on a collapsed constraint set.

**`keep` is a caller's parameter and is not chosen here.** *UNVERIFIED, quoted from `catchment-dod-balance/docs/validation.md`:* the package's
earlier Elba run gave −0.104 m at `keep = 0.99` (i.e. raise by 104 mm) at `z = 1.96` on the pre-rebuild
products; the rebuilt, floodplain-masked, correlated-envelope, `z = 1` equivalent is
**+98.7 mm**. That the two are close is a coincidence of two offsetting changes (a wider
envelope pulls it down, `z = 1.96 → 1.0` pushes it up: the same run at `z = 1.96` gives
`keep = 0.99 → +67.0 mm`), and neither is the scene answer.

The full curve, `z = 1`, canopy-corrected, correlated, floodplain out
(`net_vol_m3` = scene net volume after the raise):

    f_flagged    keep  delta_mm  se_mm   net_vol_m3
       0.0000  1.0000   +3060.2  230.2  -16,074,905
       0.0005  0.9995    +901.9  221.7   -4,773,440
       0.0010  0.9990    +497.6  159.0   -2,656,728
       0.0025  0.9975    +220.6   43.6   -1,206,299
       0.0050  0.9950    +140.5   18.8     -786,642
       0.0100  0.9900     +98.7   11.5     -567,921
       0.0250  0.9750     +56.7    8.2     -347,759
       0.0500  0.9500     +31.2    5.6     -214,179
       0.0750  0.9250     +17.4    5.0     -142,143
       0.1000  0.9000      +8.6    4.9      -96,015
       0.1500  0.8500      -4.9    4.6      -25,326
       0.2000  0.8000     -15.2    4.3      +28,461
       0.2500  0.7500     -23.9    4.3      +74,241
       0.3000  0.7000     -31.9    4.3     +115,859
       0.4000  0.6000     -45.9    4.3     +189,528
       0.5000  0.5000     -59.3    4.2     +259,542

SE column: block-bootstrap SD of the quantile over `L = 400 m` blocks, 300 replicates.
The SE itself depends on block size and has not plateaued:

      L_m  se_f0025_mm  se_f005_mm  se_f010_mm  n_blocks
       25          1.9          1.3         0.7      9812
       50          3.6          2.1         1.3      2753
      100          5.1          3.4         2.3       817
      200          7.7          5.0         3.3       233
      400          8.5          5.4         4.5        63
      800          9.7          6.8         4.9        20

---

## Comparison with the control-derived +53.6 ± 13.0 mm

    S1: scene net volume brought to +1 sigma (q2 DoD, correlated envelope, floodplain out)
        sigma component: correlated-1sigma-of-the-scene-mean-DoD = 14.70 mm
        sigma component: variogram-refit-half-range = 3.90 mm
      continuity (mass balance)      -24.5 +/- 15.2 mm (quadrature of the above)
      control marks                  +53.6 +/- 13.0 mm
      difference                     +78.1 +/- 20.0 mm  ->  3.90 sigma

    S3: flagged-cell fraction matched to the fitted error model (q2 DoD, correlated, z=1)
        sigma component: one-draw-SD-of-the-null-flag-rate = 8.25 mm
        sigma component: block-bootstrap-SE-of-the-quantile-L400 = 4.50 mm
      continuity (mass balance)       -5.8 +/-  9.4 mm (quadrature of the above)
      control marks                  +53.6 +/- 13.0 mm
      difference                     +59.4 +/- 16.0 mm  ->  3.70 sigma

    S2: total excess deposition volume matched to the fitted error model
        sigma component: one-sided-half-width-to-null-T-plus-1SD = 23.80 mm
      continuity (mass balance)       -2.3 +/- 23.8 mm (quadrature of the above)
      control marks                  +53.6 +/- 13.0 mm
      difference                     +55.9 +/- 27.1 mm  ->  2.06 sigma

**Read as a bound: consistent.** Continuity requires only `delta ≥ −24.5 mm`; +53.6 mm
satisfies that and the two instruments do not disagree.

**Read as an estimate: 2.1–3.9 σ apart**, depending on which scene aggregate is used —
and that reading rests entirely on assuming the true hillslope net export is zero.

**The equivalence that makes the assumption checkable.** The budget and the raise trade
off exactly one for one, so believing the control value is the same as asserting a real
net hillslope lowering over the epoch:

      epoch 2008-11-25 -> 2021-05-01 = 12.43 yr
      S1: 78.1 mm over the epoch = 6.28 mm/yr
      S3: 59.4 mm over the epoch = 4.78 mm/yr
      S2: 55.9 mm over the epoch = 4.50 mm/yr

averaged over the whole 5.236 km² of evaluated hillslope. Whether a Driftless cultivated
hillslope averages 4.5–6.3 mm/yr of net lowering is a geomorphic question this measurement
cannot answer, and I have not checked it against any literature — flagging it as the
decisive open item rather than asserting an answer.

---

## What a uniform raise CANNOT absorb, and what survives the best delta

A uniform raise is a one-parameter model. It cannot absorb:

- **Per-flight-line offsets.** Elba's marks show line structure at ANOVA F = 8.63,
  p < 0.001, with sd 40.7 mm between six line means *(UNVERIFIED here — quoted from
  `analysis/FRAME_2026-08-26-PM.md`, not recomputed in this task)*. A single constant is
  the wrong shape for that, and 40.7 mm of line-to-line scatter is larger than every
  scene-integrated delta above.
- **The canopy term.** Swapping the uncorrected DoD for the `q2 = 0.5 − 0.19·cover`
  product moves the scene mean by 6.8 mm and S1 by 6.2 mm, so the canopy correction and
  the datum are not separable by this measurement.
- **Real spatial change.** Deposition that genuinely happened is indistinguishable here
  from datum error; that is the whole content of the one-sided-bound caveat.

**Structure does survive.** At the S3 delta of −5.8 mm the residual flag rate is strongly
organised by canopy cover and, more weakly, by drainage area:

                    stratum       n   f_res
      drainage 0.00-0.01 ha  107281  0.1155
      drainage 0.01-0.04 ha   64368  0.1928
      drainage 0.04-0.13 ha   32184  0.2100
      drainage 0.13-0.78 ha    8583  0.1386
    drainage 0.78-494.98 ha    2146  0.1240
            cover 0.00-0.05   86179  0.1265
            cover 0.05-0.15   14547  0.1267
            cover 0.15-0.35   78908  0.1699
            cover 0.35-0.60   33849  0.1910
            cover 0.60-1.01    1079  0.3624

Open ground (cover < 0.05) sits at 0.1265 against the 0.1538 null; the densest canopy sits
at 0.3624 — **2.36× the null rate, 2.86× the open-ground residual**. *Hypothesis, not a measurement:* a residual canopy term
survives the `q2` correction at high cover. It is not tested here.

---

## Limits of this measurement, stated plainly

1. **Continuity is one-sided.** It bounds delta from below and can never demand a raise.
   Turning it into a point estimate requires assuming zero true net hillslope export.
2. **The envelope is still a lower bound** in one respect the package documents:
   flow-reconvergence cross-terms are omitted from `sigma_Vacc`.
3. **`N_eff = 6.6`** over 5.2 km² is what limits S1. A longer correlation length would
   make it weaker still; a shorter one (Carlton fits 166 m) would make it much stronger.
   The range is the single most influential fitted number here.
4. **The floodplain mask is `TPI(800 m) < −2 m`**, described in its own producer
   (`analysis/ridgelines/convexity_dod_landcover.py`) as crude. It removes 38.35 % of the
   grid. *UNVERIFIED — the package's own principled `trunk_floodplain` detector selects
   1.85–2.11 % of the tile and overlaps the crude mask at IoU 0.104, quoted from
   `catchment-dod-balance/docs/validation.md`, not recomputed here* — so these are
   different objects; the standing
   project rule names the crude one and that is what was used.
5. **Parameters chosen by me, flagged as such:** `n_boot` (bootstrap replicates, 300) and
   `null_real` (null realisations, 100) — both convergence knobs, neither changes a
   definition. Block sizes are swept, not chosen. `z = 1.0` is Andy's. `keep` is not
   chosen at all. No minimum counts, no thresholds, no filters were introduced.

## Reproduce

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/lod_cover_q2.py
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/mass_balance/datum_from_mass_balance.py \
        --dod dod_cover_q2.npy --lod lod_cover_q2.npy --envelope correlated \
        --floodplain-mode routing --z 1.0 --n-boot 300 --null-real 100 \
        --save analysis/mass_balance/out/S_q2_corr_z1.npz
