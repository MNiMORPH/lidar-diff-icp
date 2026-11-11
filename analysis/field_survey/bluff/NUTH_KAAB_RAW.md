STATUS: COMPLETE

# Nuth & Kaab test on the raw 2008 cloud, tile 4342-29-64

Question: is the ~0.4 m gen1-below-gen2 offset on the 30-degree east-facing bluff
(RTK traverse of 25 Oct 2026, grid window rows 483:550, cols 440:469) a HORIZONTAL
misregistration of gen1 relative to gen2 (aspect-dependent, Nuth & Kaab 2011
signature dh/tan(slope) = a cos(aspect - b) + c), or a slope-dependent incidence
error (dh = f(slope), no aspect dependence)?

Script: `analysis/field_survey/bluff/nuth_kaab_raw_4342-29-64.py`.
Outputs: `data/derived/4342-29-64/bluff_check/`.
Every number below is pasted from the script's printed output (the verbatim log is
the last section of this file).

## What the pipeline applied (read from the code, 2026-11-11)

- `corrections.json` for this tile records `swath_frame: "star"`, `zero_line: 135`,
  `per_swath_internal_alignment_dxdydz_m` = {135: [0,0,0], 136: [0.3176,-0.0821,-0.0179],
  137: [0.5906,-0.1707,-0.026], 138: [1.0359,0.0964,-0.0364]} and
  `per_swath_star_frame_dxdydz_m` = {135: [-0.1695,0.0798,0.0806], 136: [-0.1281,-0.0733,0.0727],
  137: [-0.1106,-0.1265,0.066], 138: [0.1579,0.084,0.0635]}, plus
  `cross_epoch_datum.horizontal_shift_m` = [-0.0373, -0.0338].
- Under `swath_frame="star"` the internal alignment (`coreg.align_swaths`) is SOLVED and
  recorded but NOT applied: `pipeline.py` passes `apply_swath_corr=(swath_frame == "chain")`
  and only adds the shifts when that is true (pipeline.py lines ~623-633, 1237).
  `registration.read_swath_alignment` says the same and refuses to fall back.
- What IS applied, in this order (`pipeline.py` ~710, `chain.run_chain`):
  1. `chain.StarFrame`: `ctx.x[m] += dx; ctx.y[m] += dy; ctx.z[m] += dz` per swath, with
     the values in `per_swath_star_frame_dxdydz_m`. Sign convention: the recorded numbers
     are ADDED to gen1's x, y, z (`starframe.fit_star_shifts` docstring: "ready to ADD to
     the coordinates").
  2. `chain.LateralShift`: a single Nuth & Kaab order-0 (dx, dy) over the whole tile,
     `ctx.x += a0; ctx.y += b0`, recorded as `horizontal_shift_m`. Also ADDED.
- So the "registered gen1 the pipeline differenced" is raw + star + lateral, NOT internal
  alignment + star. The brief's case (iii) is run both ways below and labeled.
- The gen2 reference. `z_after.npy` is the per-cell HORIZONTAL median of class 2 (Z21,
  `groundest.cellstat`); `z_after_differenced.npy` is the SLOPE-NORMAL median read at the
  cell center (Zref, `groundest.slope_normal_ground`), and `dod = Zref - gen1`. The brief
  names `z_after.npy` but describes the slope-normal median; the two differ by up to
  3.69 m in the tile and up to 0.76 m inside the bluff window (checked before writing the
  script). The slope-normal grid is the one that refers to the cell center, so it is the
  one bilinear sampling is valid on, and it is what the DoD was taken against: it is the
  PRIMARY reference here, with slope and aspect computed from it by `np.gradient`.
  `z_after.npy` is run as a sensitivity on the pooled cases. This choice is mine and is
  flagged for Andy.

## Conventions

- dh = z_gen1 - z_gen2(bilinear at the gen1 return's x, y), meters.
- slope, aspect from `np.gradient` of the reference grid (5 m); aspect = DOWNSLOPE
  direction, degrees clockwise from north (east-facing = 90).
- Model: dh / tan(slope) = a cos(aspect - b) + c, fitted by unweighted least squares to
  the medians of 36 ten-degree aspect bins (points with slope > 15 deg, non-floodplain).
  Derivation: a gen1 return truly at p but recorded at p + s has z(p - s) at its label,
  so dh = -s . grad z = tan(slope) (s_x sin(aspect) + s_y cos(aspect)). Hence
  s_x = a sin(b), s_y = a cos(b) is the displacement of gen1 RELATIVE TO gen2 (gen1
  features sit s east/north of where gen2 has them), and the CORRECTION to add to gen1's
  coordinates is (-s_x, -s_y). The script checks this sign numerically by applying the
  correction, resampling and refitting.
- Uncertainties: the LSQ covariance from the scatter of the 36 bin medians about the
  fit, and a 1000-rep bootstrap over bins.
- NMAD = 1.4826 * median |y - median(y)|.
- The slope > 15 deg band edge is Andy's to change; the 36-bin count and the 0-5-...-35-45
  slope bands are from the brief.


## Step 2. Nuth & Kaab fits (slope > 15 deg, non-floodplain), table built from `bluff_check/nuth_kaab_fits.csv`

s_x, s_y = fitted displacement of gen1 RELATIVE TO gen2 (m); the correction to add to gen1 is (-s_x, -s_y). sigma = bootstrap over the 36 bins (the LSQ-covariance sigma is in the CSV and is within ~20 % of it). NMAD columns are before -> after removing the fit, in m. 'registered' = raw + star + lateral, the shifts the pipeline APPLIED (confirmed by the coordinator against pipeline.py line 1237); 'internal+star' is the brief's original wording, run only as a labeled hypothetical. The CSF cache holds raw positions (same xyz minima as the raw file) and was shifted the same way.

| case | n returns | a (m) | b (deg) | c (m) | s_x = dx (m) | s_y = dy (m) | sig s_x | sig s_y | a boot 95 % | NMAD dh/tan before -> after | NMAD dh before -> after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| i raw pooled | 2717095 | 0.1406 | 49.0 | -0.1312 | +0.1061 | +0.0923 | 0.0114 | 0.0062 | [0.1236, 0.1589] | 0.3811 -> 0.3664 | 0.1714 -> 0.1646 |
| ii raw swath 135 | 131627 | 0.1917 | 113.2 | -0.0815 | +0.1762 | -0.0755 | 0.0196 | 0.0148 | [0.1586, 0.2276] | 0.4148 -> 0.3603 | 0.1801 -> 0.1569 |
| ii raw swath 136 | 958456 | 0.1642 | 59.3 | -0.1464 | +0.1412 | +0.0839 | 0.0126 | 0.0074 | [0.1428, 0.1852] | 0.3719 -> 0.3474 | 0.1679 -> 0.1572 |
| ii raw swath 137 | 1069080 | 0.2153 | 52.7 | -0.1062 | +0.1713 | +0.1304 | 0.0139 | 0.0109 | [0.1924, 0.2414] | 0.3759 -> 0.3410 | 0.1665 -> 0.1508 |
| ii raw swath 138 | 557932 | 0.1900 | 235.7 | -0.0561 | -0.1570 | -0.1070 | 0.0107 | 0.0106 | [0.1676, 0.2136] | 0.3954 -> 0.3736 | 0.1835 -> 0.1728 |
| iii-a registered(star+lateral, AS APPLIED) pooled | 2716570 | 0.0194 | 204.2 | +0.0382 | -0.0079 | -0.0177 | 0.0108 | 0.0068 | [0.0082, 0.0378] | 0.3520 -> 0.3520 | 0.1574 -> 0.1574 |
| iii-a registered(star+lateral) swath 135 | 131569 | 0.0340 | 239.0 | +0.1063 | -0.0291 | -0.0175 | 0.0210 | 0.0157 | [0.0088, 0.0776] | 0.3587 -> 0.3599 | 0.1564 -> 0.1567 |
| iii-a registered(star+lateral) swath 136 | 958341 | 0.0308 | 231.1 | +0.0160 | -0.0239 | -0.0193 | 0.0135 | 0.0079 | [0.0147, 0.0541] | 0.3505 -> 0.3495 | 0.1575 -> 0.1572 |
| iii-a registered(star+lateral) swath 137 | 1068913 | 0.0305 | 144.1 | +0.0454 | +0.0179 | -0.0247 | 0.0148 | 0.0110 | [0.0128, 0.0562] | 0.3402 -> 0.3413 | 0.1502 -> 0.1505 |
| iii-a registered(star+lateral) swath 138 | 557747 | 0.0695 | 213.9 | +0.0873 | -0.0387 | -0.0577 | 0.0117 | 0.0109 | [0.0471, 0.0943] | 0.3740 -> 0.3757 | 0.1724 -> 0.1729 |
| iii-b internal+star (HYPOTHETICAL, not applied) pooled | 2716044 | 0.5882 | 100.9 | +0.0198 | +0.5776 | -0.1111 | 0.0138 | 0.0117 | [0.5619, 0.6156] | 0.6524 -> 0.4222 | 0.2893 -> 0.1881 |
| iii-b internal+star (hypothetical) swath 135 | 131557 | 0.0199 | 27.7 | +0.1060 | +0.0093 | +0.0177 | 0.0206 | 0.0155 | [0.0056, 0.0601] | 0.3600 -> 0.3597 | 0.1570 -> 0.1568 |
| iii-b internal+star (hypothetical) swath 136 | 958504 | 0.3411 | 101.6 | -0.0246 | +0.3341 | -0.0687 | 0.0135 | 0.0077 | [0.3156, 0.3691] | 0.4541 -> 0.3485 | 0.2026 -> 0.1573 |
| iii-b internal+star (hypothetical) swath 137 | 1069337 | 0.6753 | 104.1 | -0.0151 | +0.6551 | -0.1641 | 0.0142 | 0.0107 | [0.6463, 0.7009] | 0.6870 -> 0.3411 | 0.3034 -> 0.1509 |
| iii-b internal+star (hypothetical) swath 138 | 556646 | 1.0577 | 85.2 | -0.0014 | +1.0540 | +0.0892 | 0.0101 | 0.0092 | [1.0368, 1.0773] | 0.9590 -> 0.3745 | 0.4682 -> 0.1731 |
| iv CSF cache registered(star+lateral) pooled | 1625190 | 0.0409 | 255.2 | -0.0031 | -0.0395 | -0.0105 | 0.0124 | 0.0082 | [0.0182, 0.0649] | 0.3285 -> 0.3277 | 0.1368 -> 0.1366 |
| iv CSF cache registered swath 135 | 96631 | 0.0715 | 265.8 | +0.0713 | -0.0713 | -0.0053 | 0.0268 | 0.0204 | [0.0264, 0.1246] | 0.3329 -> 0.3373 | 0.1393 -> 0.1413 |
| iv CSF cache registered swath 136 | 645237 | 0.0399 | 262.5 | -0.0155 | -0.0396 | -0.0052 | 0.0147 | 0.0093 | [0.0158, 0.0689] | 0.3405 -> 0.3412 | 0.1438 -> 0.1441 |
| iv CSF cache registered swath 137 | 676547 | 0.0533 | 243.2 | +0.0008 | -0.0476 | -0.0241 | 0.0142 | 0.0106 | [0.0279, 0.0787] | 0.3062 -> 0.3039 | 0.1273 -> 0.1264 |
| iv CSF cache registered swath 138 | 206775 | 0.0712 | 216.1 | +0.0448 | -0.0420 | -0.0575 | 0.0191 | 0.0181 | [0.0373, 0.1127] | 0.3642 -> 0.3632 | 0.1467 -> 0.1460 |
| sens: i raw pooled vs z_after.npy (Z21) | 2713717 | 0.1314 | 50.4 | -0.1250 | +0.1012 | +0.0839 | 0.0102 | 0.0060 | [0.1159, 0.1478] | 0.4148 -> 0.4046 | 0.1852 -> 0.1804 |
| sens: iii-a registered pooled vs z_after.npy (Z21) | 2713277 | 0.0285 | 205.6 | +0.0450 | -0.0123 | -0.0257 | 0.0106 | 0.0062 | [0.0168, 0.0458] | 0.3917 -> 0.3916 | 0.1744 -> 0.1744 |

Reading the table:

- Sign check passed: adding (-s_x, -s_y) to the raw gen1 coordinates and refitting took a from 0.1406 to 0.0023 m.
- The RAW cloud carries a real horizontal misregistration, and it differs by swath: pooled a = 0.1406 m towards b = 49.0 deg, but per swath s = (+0.1762, -0.0755), (+0.1412, +0.0839), (+0.1713, +0.1304), (-0.1570, -0.1070) m for 135, 136, 137, 138. Swath 138, the one under the bluff, is displaced the OPPOSITE way from the other three.
- The pipeline's APPLIED registration (star + lateral) removes it: pooled a = 0.0194 m (boot 95 % [0.0082, 0.0378]), and the per-swath residuals are 0.0305-0.0695 m, with swath 138 the largest at a = 0.0695 m, s = (-0.0387, -0.0577) m. The NMAD of dh/tan(slope) does not move (0.3520 -> 0.3520): there is no aspect signal left to remove.
- The brief's case (iii) as first written, internal alignment + star, is NOT what the pipeline applied and would be badly wrong: pooled a = 0.5882 m (swath 138 alone 1.0577 m). The star frame was fitted on the UNALIGNED coordinates, so adding align_swaths' horizontal solution on top double-counts; and align_swaths' horizontal half is the part the project already found to be wrong (chain.StarFrame docstring).
- The CSF cache (the pipeline's actual ground population), registered: a = 0.0409 m, s = (-0.0395, -0.0105) m; swath 138: a = 0.0712 m.
- Sensitivity to the reference grid: against `z_after.npy` (horizontal median) the raw pooled fit is a = 0.1314 m at 50.4 deg and the registered one a = 0.0285 m at 205.6 deg, against 0.1406 / 49.0 and 0.0194 / 204.2 on the slope-normal grid. Same answer.
- Independent check of the registration, computed in the addendum: the shift the pipeline applied per swath (star + lateral) agrees with -s from the raw per-swath Nuth & Kaab to (-0.031, -0.030), (-0.024, -0.023), (+0.023, -0.030), (-0.036, -0.057) m for 135-138. Two different estimators (StarFrame on CSF ground over stable cells; 36-bin medians of raw last returns on slopes > 15 deg) give the same per-swath horizontal shifts to 2-6 cm.

## Step 3. Median dh (m) by gen2 slope band, pooled and by aspect quadrant, from `bluff_check/slope_band_medians.csv`

'after fit' moves gen1 by (-s_x, -s_y) and RESAMPLES gen2; the linear-removal check in the log agrees with it to 1 mm in every cell of the table.

**raw, no shift**

| band (deg) | n | pooled | E 45-135 | S 135-225 | W 225-315 | N 315-45 |
|---|---:|---:|---:|---:|---:|---:|
| 0-5 | 1629113 | -0.0606 | -0.0622 | -0.0655 | -0.0558 | -0.0585 |
| 5-10 | 1339500 | -0.0501 | -0.0363 | -0.0560 | -0.0699 | -0.0477 |
| 10-15 | 811857 | -0.0524 | -0.0393 | -0.0700 | -0.0750 | -0.0400 |
| 15-20 | 685097 | -0.0468 | -0.0294 | -0.0779 | -0.0745 | -0.0256 |
| 20-25 | 751712 | -0.0433 | -0.0223 | -0.0773 | -0.0866 | -0.0212 |
| 25-30 | 800403 | -0.0644 | -0.0333 | -0.0909 | -0.1275 | -0.0246 |
| 30-35 | 394196 | -0.0663 | -0.0171 | -0.1124 | -0.1229 | -0.0474 |
| 35-45 | 77835 | -0.0218 | +0.0028 | -0.0679 | -0.1197 | +0.0137 |

**registered (star+lateral), before fit**

| band (deg) | n | pooled | E 45-135 | S 135-225 | W 225-315 | N 315-45 |
|---|---:|---:|---:|---:|---:|---:|
| 0-5 | 1629131 | +0.0094 | +0.0022 | +0.0091 | +0.0216 | +0.0074 |
| 5-10 | 1339680 | +0.0167 | +0.0200 | +0.0245 | +0.0156 | +0.0090 |
| 10-15 | 811893 | +0.0116 | +0.0106 | +0.0173 | +0.0161 | +0.0061 |
| 15-20 | 684927 | +0.0176 | +0.0130 | +0.0215 | +0.0281 | +0.0132 |
| 20-25 | 751369 | +0.0186 | +0.0207 | +0.0271 | +0.0230 | +0.0112 |
| 25-30 | 800419 | +0.0090 | +0.0051 | +0.0196 | -0.0022 | +0.0112 |
| 30-35 | 394123 | +0.0082 | +0.0252 | +0.0110 | -0.0028 | -0.0158 |
| 35-45 | 77904 | +0.0390 | +0.0483 | +0.0554 | +0.0222 | -0.0009 |

**registered + pooled fitted shift (+0.008,+0.018) resampled**

| band (deg) | n | pooled | E 45-135 | S 135-225 | W 225-315 | N 315-45 |
|---|---:|---:|---:|---:|---:|---:|
| 0-5 | 1629066 | +0.0094 | +0.0026 | +0.0083 | +0.0213 | +0.0082 |
| 5-10 | 1339680 | +0.0170 | +0.0208 | +0.0227 | +0.0148 | +0.0111 |
| 10-15 | 811936 | +0.0122 | +0.0119 | +0.0141 | +0.0147 | +0.0096 |
| 15-20 | 684913 | +0.0183 | +0.0149 | +0.0166 | +0.0263 | +0.0181 |
| 20-25 | 751346 | +0.0201 | +0.0232 | +0.0209 | +0.0209 | +0.0175 |
| 25-30 | 800488 | +0.0089 | +0.0076 | +0.0118 | -0.0039 | +0.0183 |
| 30-35 | 394093 | +0.0070 | +0.0263 | +0.0019 | -0.0046 | -0.0071 |
| 35-45 | 77919 | +0.0383 | +0.0486 | +0.0451 | +0.0166 | +0.0131 |

**registered + pooled fit removed LINEARLY (check)**

| band (deg) | n | pooled | E 45-135 | S 135-225 | W 225-315 | N 315-45 |
|---|---:|---:|---:|---:|---:|---:|
| 0-5 | 1629131 | +0.0094 | +0.0026 | +0.0083 | +0.0213 | +0.0082 |
| 5-10 | 1339680 | +0.0170 | +0.0208 | +0.0226 | +0.0147 | +0.0111 |
| 10-15 | 811893 | +0.0122 | +0.0119 | +0.0141 | +0.0147 | +0.0096 |
| 15-20 | 684927 | +0.0183 | +0.0149 | +0.0167 | +0.0264 | +0.0181 |
| 20-25 | 751369 | +0.0201 | +0.0232 | +0.0210 | +0.0210 | +0.0174 |
| 25-30 | 800419 | +0.0089 | +0.0075 | +0.0119 | -0.0039 | +0.0181 |
| 30-35 | 394123 | +0.0070 | +0.0261 | +0.0021 | -0.0045 | -0.0073 |
| 35-45 | 77904 | +0.0384 | +0.0485 | +0.0455 | +0.0177 | +0.0125 |

**registered + PER-SWATH fitted shifts resampled**

| band (deg) | n | pooled | E 45-135 | S 135-225 | W 225-315 | N 315-45 |
|---|---:|---:|---:|---:|---:|---:|
| 0-5 | 1629076 | +0.0095 | +0.0028 | +0.0080 | +0.0211 | +0.0085 |
| 5-10 | 1339713 | +0.0171 | +0.0212 | +0.0219 | +0.0144 | +0.0119 |
| 10-15 | 811989 | +0.0126 | +0.0126 | +0.0126 | +0.0145 | +0.0118 |
| 15-20 | 684780 | +0.0190 | +0.0160 | +0.0144 | +0.0262 | +0.0210 |
| 20-25 | 751417 | +0.0213 | +0.0253 | +0.0178 | +0.0200 | +0.0212 |
| 25-30 | 800518 | +0.0095 | +0.0105 | +0.0083 | -0.0069 | +0.0235 |
| 30-35 | 394020 | +0.0073 | +0.0267 | -0.0031 | -0.0070 | +0.0019 |
| 35-45 | 77932 | +0.0368 | +0.0465 | +0.0367 | +0.0158 | +0.0196 |

Reading the tables (registered + pooled fitted shift, each quadrant relative to its own 0-5 deg band, numbers from the log's 'minus that' lines):

- Pooled: 0-5 +0.0000, 5-10 +0.0076, 10-15 +0.0028, 15-20 +0.0089, 20-25 +0.0107, 25-30 -0.0005, 30-35 -0.0024, 35-45 +0.0289. There is NO gen1-reads-low-on-steep-slopes trend in the raw last returns once registered: the 25-35 deg bands sit within 3 mm of the flat ground.
- By quadrant at 25-30 / 30-35 / 35-45 deg: E +0.0050 / +0.0237 / +0.0461; S +0.0035 / -0.0064 / +0.0368; W -0.0252 / -0.0258 / -0.0046; N +0.0101 / -0.0153 / +0.0050. The quadrants are NOT the same: east-facing steep ground reads 2-5 cm HIGH and west-facing reads 2-3 cm LOW relative to their own flat ground. That E/W antisymmetry is the signature of a residual east-west shift confined to the steepest ground (about 0.025 m / tan(30 deg) = 0.04 m), the same size as the per-swath residuals in step 2, and it has the WRONG sign for the bluff: a residual that makes east-facing gen1 read HIGH cannot be what makes the bluff's gen1 read low.
- Before registration (raw), the quadrant pattern is the mirror image and larger: at 30-35 deg E +0.0452, W -0.0671 relative to flat; registration removes most of it.

## Step 4. The bluff window (rows 483:550, cols 440:469), from `bluff_check/window_slope_band_medians.csv` and the log

- Window: x [579695, 579840), y [4885155, 4885490), 67 x 29 = 1943 cells; gen2 slope median 28.6 deg (p10/p90 13.7/34.7), aspect median 116.6 deg; 28 floodplain cells; 1181 cells steeper than 25 deg, median aspect 116.1 deg.
- Swaths: 30494 gen1 raw last returns in the window, swath 138: 29002, swath 137: 1492 (swaths 135, 136: none). Swath 138 reaches all 1943 cells; swath 137 touches 190 cells at its edge. On cells steeper than 25 deg: 17875 returns from 138, 700 from 137. This agrees with the other agent's finding that every 2008 return at the RTK traverse itself is swath 138. The CSF cache holds only 5708 of those 30494 returns (138: 5519, 137: 189), i.e. the pipeline's ground selection discarded 81 % of the window's last returns.

| case | 0-5 | 5-10 | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-45 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| raw, no shift | -0.1233 (n=1029) | -0.1433 (n=683) | -0.1390 (n=1131) | -0.1488 (n=4142) | -0.1786 (n=4574) | -0.2053 (n=6515) | -0.2165 (n=9486) | -0.2366 (n=2548) |
| registered (star+lateral) | -0.0580 (n=1035) | -0.0709 (n=685) | -0.0540 (n=1130) | -0.0714 (n=4137) | -0.0912 (n=4592) | -0.0990 (n=6484) | -0.1032 (n=9501) | -0.1022 (n=2544) |
| registered + pooled fitted shift (+0.008,+0.018) | -0.0581 (n=1035) | -0.0711 (n=685) | -0.0526 (n=1130) | -0.0694 (n=4136) | -0.0899 (n=4590) | -0.0991 (n=6485) | -0.1045 (n=9500) | -0.1005 (n=2547) |
| registered + per-swath fitted shifts | -0.0580 (n=1034) | -0.0697 (n=684) | -0.0481 (n=1131) | -0.0642 (n=4141) | -0.0865 (n=4583) | -0.0936 (n=6482) | -0.0999 (n=9501) | -0.0894 (n=2550) |
| registered, swath 137 only | nan (n=0) | +0.1404 (n=4) | -0.0033 (n=8) | +0.0063 (n=135) | -0.1122 (n=662) | -0.1556 (n=440) | -0.2249 (n=243) | nan (n=0) |
| registered, swath 138 only | -0.0580 (n=1035) | -0.0733 (n=681) | -0.0541 (n=1122) | -0.0735 (n=4002) | -0.0885 (n=3930) | -0.0974 (n=6044) | -0.1017 (n=9258) | -0.1022 (n=2544) |
| raw + raw-pooled fitted shift (-0.106,-0.092) | -0.1241 (n=1024) | -0.1504 (n=689) | -0.1558 (n=1133) | -0.1646 (n=4139) | -0.2007 (n=4574) | -0.2357 (n=6512) | -0.2490 (n=9491) | -0.2888 (n=2548) |

- Registered, slope > 25 deg: n = 18529, median dh = -0.1015 m, median aspect 116.0 deg, median dh/tan(slope) = -0.1671 m. Registered, slope < 10 deg: n = 1720, median dh = -0.0611 m. So the RAW LAST RETURNS on the bluff sit 0.10 m below gen2, of which 0.06 m is the window's local vertical offset shared by its flat ground; the steep-minus-flat part is 0.04 m.
- The fitted shifts do nothing here: after the pooled fit the steep median is -0.1021 m, after the per-swath (138) fit -0.0967 m. The raw-pooled fit applied to the raw cloud makes the window WORSE (35-45 deg: -0.2366 -> -0.2888) because swath 138's shift has the opposite sign from the pooled one; that row is a demonstration, not a candidate correction.
- Swath 137's few edge returns in the window read lower than 138's on the same cells: 25-30 deg -0.1556 vs -0.0974, 30-35 deg -0.2249 vs -0.1017 m (n = 440 / 243 for 137). This is a per-swath difference on the bluff that the per-swath horizontal fits do NOT predict (137's residual s_x = +0.0179 m would make it read HIGHER on an east-facing slope). Left open; see step 5.
- The DoD in the window. Pipeline `dod.npy` (gen2 - gen1, CSF slope-normal ground): finite on 906 of 1943 cells, median +0.1461 m; on steep (> 25 deg) cells +0.2473 m (n = 290 of 1181 finite); on gentle (< 10 deg) cells +0.0556 m. A per-cell proxy from the raw last returns (median over the cell's returns of gen2_bilinear - gen1), registered: finite on 1916 cells, median +0.0928 m; steep +0.1023 m; gentle +0.0684 m. On the 290 steep cells where BOTH are finite: pipeline +0.2473 m, proxy +0.1239 m; proxy after the pooled fitted shift +0.1247 m, after the per-swath shift +0.1195 m. The proxy and the pipeline grid are different estimators (median of all last returns vs CSF-selected slope-normal median), so compare the proxy to itself across shifts: the fitted shifts move the window DoD by at most 4 mm.
- Under my fitted shift the window DoD would therefore be, to within those 4 mm, what it is now; the fitted horizontal registration is not a lever on this bluff.

## Step 5. Which hypothesis the data favor

**The data favor the non-horizontal explanation, and the one number that decides it is the residual Nuth & Kaab amplitude after the pipeline's registration: a = 0.0194 m (bootstrap 95 % [0.0082, 0.0378]) pooled, 0.0695 m for swath 138, against the 0.693 m of shift along the bluff's downslope direction that a 0.40 m offset on a 30 deg slope requires (36x the fit, 18x its upper bound), or 0.352 m for the 0.1917 m steep-minus-gentle offset that `dod.npy` shows in the window.** Projected onto the bluff (aspect 116 deg, 30 deg), the registered residual shifts predict gen1 - gen2 of +0.0004 +/- 0.0059 m (pooled), -0.0055 +/- 0.0067 m (swath 138) and -0.0178 +/- 0.0068 m (CSF population): at most 2 cm of the bluff offset is horizontal misregistration.

What the raw cloud says instead: (1) the raw gen1 last returns on the bluff are only 0.10 m below gen2 (0.04 m steep-minus-flat), not 0.25-0.40 m, and the raw last returns show no gen1-low-on-steep-slopes trend tile-wide in any quadrant; (2) the pipeline's DoD on the same steep window cells is +0.2473 m where the raw-return proxy is +0.1239 m, and the CSF selection kept 19 % of the window's last returns, so the extra 0.12 m, and the blank cells, come from WHICH gen1 returns are called ground on the bluff (the previous commit's finding), not from where they are placed horizontally; (3) tile-wide the raw cloud DID carry per-swath horizontal shifts of 0.16-0.22 m in four different directions, and the pipeline's star + lateral registration removed them to 2-7 cm, reproducing the raw Nuth & Kaab per-swath corrections to 2-6 cm with an independent estimator.

It is also not a simple incidence error of the form dh = f(slope) independent of aspect, at least not in the raw last returns: the steep-band residuals differ by quadrant (E +0.024, W -0.026 m at 30-35 deg) at the few-cm level, and the pooled slope trend is zero. Whatever makes the CSF ground read 0.2-0.4 m low on this bluff is not visible in the raw last-return population as either signature.

What I could not determine:

- Where the '~0.4 m' comes from: `dod.npy` on the window's steep cells has median +0.2473 m over the 290 cells where it is finite; the RTK comparison at individual points may see more. I did not touch the RTK data.
- Why swath 137's edge returns on the bluff read 0.06-0.12 m lower than swath 138's on the same cells. Not horizontal (its residual shift has the wrong sign); candidate causes are scan-angle / far-edge geometry or 137's own vertical frame at its edge. Needs the beam-angle table, outside this brief.
- The residual E/W antisymmetry at 25-35 deg (about 4 cm of east-west shift confined to steep ground) could be a slope-dependent horizontal error (e.g. a residual roll or timing term) rather than a constant shift; a Nuth & Kaab fit per slope band would separate it. Not run; outside the brief.
- Only one bluff aspect is in the window, so the window alone cannot separate a horizontal shift from a vertical offset; the tile-wide fit is what does.
- Note: `bluff_check/local/` and `bluff_check/translation_fit.csv` were written by another process during this run (16:41-16:43); they are not outputs of this script and were not read.

Thresholds that are mine, not the brief's: the primary-reference choice (`z_after_differenced.npy`), the 1000-rep bin bootstrap, and the '> 25 deg' / '< 10 deg' reporting cuts for 'steep' and 'gentle' cells inside the window. All are flagged as proposals.

## Script output (verbatim, appended as produced)

```
# nuth_kaab_raw_4342-29-64.py  run 2026-11-11 16:43:27
corrections.json: swath_frame=star zero_line=135
  per_swath_star_frame_dxdydz_m (APPLIED, added to x,y,z): {135: [-0.1695, 0.0798, 0.0806], 136: [-0.1281, -0.0733, 0.0727], 137: [-0.1106, -0.1265, 0.066], 138: [0.1579, 0.084, 0.0635]}
  per_swath_internal_alignment_dxdydz_m (SOLVED, NOT applied under star): {135: [0.0, 0.0, 0.0], 136: [0.3176, -0.0821, -0.0179], 137: [0.5906, -0.1707, -0.026], 138: [1.0359, 0.0964, -0.0364]}
  cross_epoch_datum.horizontal_shift_m (APPLIED after star, added to x,y): [-0.0373, -0.0338]
grid 699x507 res 5.0 bounds [577495.0, 4882740.0, 580030.0, 4886235.0]; floodplain cells 38152
raw cloud /home/awickert/projects/lidar-diff-icp/data/before/4342-29-64.laz: 7728747 returns, 7254886 last returns, read in 2.3 s
  last returns per swath: {135: 525918, 136: 2555164, 137: 2703087, 138: 1470717}
csf cache /home/awickert/projects/lidar-diff-icp/data/csf_cache/4342-29-64.las: 5996805 points (all class 2), read in 0.2 s

=== STEP 2: Nuth & Kaab fits (slope > 15 deg, non-floodplain) ===
primary reference: Zref(z_after_differenced, slope-normal median); sensitivity reference: Z21(z_after, horizontal median)
[i raw pooled] n=2717095 bins=36 (min bin n=24675) a=0.1406 m b=49.0 deg c=-0.1312 m | s_x=+0.1061±0.0091(lsq)/0.0114(boot) s_y=+0.0923±0.0091/0.0062 | a boot 95% [0.1236,0.1589] | correction to add to gen1 (dx,dy)=(-0.1061,-0.0923) | NMAD dh/tan: 0.3811 -> 0.3664 m | NMAD dh: 0.1714 -> 0.1646 m | median dh -0.0526 m, median slope 24.5, bin rms resid 0.0384
[i raw pooled, after adding (-s_x,-s_y) to gen1 x,y [SIGN CHECK]] n=2716971 bins=36 (min bin n=24688) a=0.0023 m b=228.9 deg c=-0.1312 m | s_x=-0.0017±0.0091(lsq)/0.0114(boot) s_y=-0.0015±0.0091/0.0063 | a boot 95% [0.0020,0.0267] | correction to add to gen1 (dx,dy)=(+0.0017,+0.0015) | NMAD dh/tan: 0.3663 -> 0.3663 m | NMAD dh: 0.1646 -> 0.1646 m | median dh -0.0591 m, median slope 24.5, bin rms resid 0.0384
  sign check: a went 0.1406 -> 0.0023 m; a residual a near zero confirms the sign convention
[ii raw swath 135] n=131627 bins=36 (min bin n=117) a=0.1917 m b=113.2 deg c=-0.0815 m | s_x=+0.1762±0.0180(lsq)/0.0196(boot) s_y=-0.0755±0.0180/0.0148 | a boot 95% [0.1586,0.2276] | correction to add to gen1 (dx,dy)=(-0.1762,+0.0755) | NMAD dh/tan: 0.4148 -> 0.3603 m | NMAD dh: 0.1801 -> 0.1569 m | median dh -0.0747 m, median slope 23.8, bin rms resid 0.0762
[ii raw swath 136] n=958456 bins=36 (min bin n=10246) a=0.1642 m b=59.3 deg c=-0.1464 m | s_x=+0.1412±0.0104(lsq)/0.0126(boot) s_y=+0.0839±0.0104/0.0074 | a boot 95% [0.1428,0.1852] | correction to add to gen1 (dx,dy)=(-0.1412,-0.0839) | NMAD dh/tan: 0.3719 -> 0.3474 m | NMAD dh: 0.1679 -> 0.1572 m | median dh -0.0636 m, median slope 24.7, bin rms resid 0.0440
[ii raw swath 137] n=1069080 bins=36 (min bin n=9294) a=0.2153 m b=52.7 deg c=-0.1062 m | s_x=+0.1713±0.0123(lsq)/0.0139(boot) s_y=+0.1304±0.0123/0.0109 | a boot 95% [0.1924,0.2414] | correction to add to gen1 (dx,dy)=(-0.1713,-0.1304) | NMAD dh/tan: 0.3759 -> 0.3410 m | NMAD dh: 0.1665 -> 0.1508 m | median dh -0.0402 m, median slope 24.0, bin rms resid 0.0520
[ii raw swath 138] n=557932 bins=36 (min bin n=879) a=0.1900 m b=235.7 deg c=-0.0561 m | s_x=-0.1570±0.0107(lsq)/0.0107(boot) s_y=-0.1070±0.0107/0.0106 | a boot 95% [0.1676,0.2136] | correction to add to gen1 (dx,dy)=(+0.1570,+0.1070) | NMAD dh/tan: 0.3954 -> 0.3736 m | NMAD dh: 0.1835 -> 0.1728 m | median dh -0.0538 m, median slope 25.2, bin rms resid 0.0453
[iii-a registered(star+lateral, AS APPLIED) pooled] n=2716570 bins=36 (min bin n=24690) a=0.0194 m b=204.2 deg c=0.0382 m | s_x=-0.0079±0.0092(lsq)/0.0108(boot) s_y=-0.0177±0.0092/0.0068 | a boot 95% [0.0082,0.0378] | correction to add to gen1 (dx,dy)=(+0.0079,+0.0177) | NMAD dh/tan: 0.3520 -> 0.3520 m | NMAD dh: 0.1574 -> 0.1574 m | median dh +0.0146 m, median slope 24.5, bin rms resid 0.0389
[iii-a registered(star+lateral) swath 135] n=131569 bins=36 (min bin n=116) a=0.0340 m b=239.0 deg c=0.1063 m | s_x=-0.0291±0.0184(lsq)/0.0210(boot) s_y=-0.0175±0.0184/0.0157 | a boot 95% [0.0088,0.0776] | correction to add to gen1 (dx,dy)=(+0.0291,+0.0175) | NMAD dh/tan: 0.3587 -> 0.3599 m | NMAD dh: 0.1564 -> 0.1567 m | median dh +0.0273 m, median slope 23.8, bin rms resid 0.0782
[iii-a registered(star+lateral) swath 136] n=958341 bins=36 (min bin n=10270) a=0.0308 m b=231.1 deg c=0.0160 m | s_x=-0.0239±0.0112(lsq)/0.0135(boot) s_y=-0.0193±0.0112/0.0079 | a boot 95% [0.0147,0.0541] | correction to add to gen1 (dx,dy)=(+0.0239,+0.0193) | NMAD dh/tan: 0.3505 -> 0.3495 m | NMAD dh: 0.1575 -> 0.1572 m | median dh +0.0056 m, median slope 24.7, bin rms resid 0.0476
[iii-a registered(star+lateral) swath 137] n=1068913 bins=36 (min bin n=9307) a=0.0305 m b=144.1 deg c=0.0454 m | s_x=+0.0179±0.0130(lsq)/0.0148(boot) s_y=-0.0247±0.0130/0.0110 | a boot 95% [0.0128,0.0562] | correction to add to gen1 (dx,dy)=(-0.0179,+0.0247) | NMAD dh/tan: 0.3402 -> 0.3413 m | NMAD dh: 0.1502 -> 0.1505 m | median dh +0.0159 m, median slope 24.0, bin rms resid 0.0554
[iii-a registered(star+lateral) swath 138] n=557747 bins=36 (min bin n=881) a=0.0695 m b=213.9 deg c=0.0873 m | s_x=-0.0387±0.0116(lsq)/0.0117(boot) s_y=-0.0577±0.0116/0.0109 | a boot 95% [0.0471,0.0943] | correction to add to gen1 (dx,dy)=(+0.0387,+0.0577) | NMAD dh/tan: 0.3740 -> 0.3757 m | NMAD dh: 0.1724 -> 0.1729 m | median dh +0.0246 m, median slope 25.2, bin rms resid 0.0492
[iii-b internal+star (HYPOTHETICAL, not applied) pooled] n=2716044 bins=36 (min bin n=24646) a=0.5882 m b=100.9 deg c=0.0198 m | s_x=+0.5776±0.0130(lsq)/0.0138(boot) s_y=-0.1111±0.0130/0.0117 | a boot 95% [0.5619,0.6156] | correction to add to gen1 (dx,dy)=(-0.5776,+0.1111) | NMAD dh/tan: 0.6524 -> 0.4222 m | NMAD dh: 0.2893 -> 0.1881 m | median dh +0.0443 m, median slope 24.5, bin rms resid 0.0551
[iii-b internal+star (hypothetical) swath 135] n=131557 bins=36 (min bin n=115) a=0.0199 m b=27.7 deg c=0.1060 m | s_x=+0.0093±0.0183(lsq)/0.0206(boot) s_y=+0.0177±0.0183/0.0155 | a boot 95% [0.0056,0.0601] | correction to add to gen1 (dx,dy)=(-0.0093,-0.0177) | NMAD dh/tan: 0.3600 -> 0.3597 m | NMAD dh: 0.1570 -> 0.1568 m | median dh +0.0255 m, median slope 23.8, bin rms resid 0.0776
[iii-b internal+star (hypothetical) swath 136] n=958504 bins=36 (min bin n=10266) a=0.3411 m b=101.6 deg c=-0.0246 m | s_x=+0.3341±0.0109(lsq)/0.0135(boot) s_y=-0.0687±0.0109/0.0077 | a boot 95% [0.3156,0.3691] | correction to add to gen1 (dx,dy)=(-0.3341,+0.0687) | NMAD dh/tan: 0.4541 -> 0.3485 m | NMAD dh: 0.2026 -> 0.1573 m | median dh +0.0090 m, median slope 24.7, bin rms resid 0.0464
[iii-b internal+star (hypothetical) swath 137] n=1069337 bins=36 (min bin n=9277) a=0.6753 m b=104.1 deg c=-0.0151 m | s_x=+0.6551±0.0127(lsq)/0.0142(boot) s_y=-0.1641±0.0127/0.0107 | a boot 95% [0.6463,0.7009] | correction to add to gen1 (dx,dy)=(-0.6551,+0.1641) | NMAD dh/tan: 0.6870 -> 0.3411 m | NMAD dh: 0.3034 -> 0.1509 m | median dh +0.0500 m, median slope 24.0, bin rms resid 0.0538
[iii-b internal+star (hypothetical) swath 138] n=556646 bins=36 (min bin n=861) a=1.0577 m b=85.2 deg c=-0.0014 m | s_x=+1.0540±0.0100(lsq)/0.0101(boot) s_y=+0.0892±0.0100/0.0092 | a boot 95% [1.0368,1.0773] | correction to add to gen1 (dx,dy)=(-1.0540,-0.0892) | NMAD dh/tan: 0.9590 -> 0.3745 m | NMAD dh: 0.4682 -> 0.1731 m | median dh +0.2211 m, median slope 25.2, bin rms resid 0.0426
[iv CSF cache registered(star+lateral) pooled] n=1625190 bins=36 (min bin n=14680) a=0.0409 m b=255.2 deg c=-0.0031 m | s_x=-0.0395±0.0104(lsq)/0.0124(boot) s_y=-0.0105±0.0104/0.0082 | a boot 95% [0.0182,0.0649] | correction to add to gen1 (dx,dy)=(+0.0395,+0.0105) | NMAD dh/tan: 0.3285 -> 0.3277 m | NMAD dh: 0.1368 -> 0.1366 m | median dh +0.0006 m, median slope 22.2, bin rms resid 0.0439
[iv CSF cache registered swath 135] n=96631 bins=36 (min bin n=103) a=0.0715 m b=265.8 deg c=0.0713 m | s_x=-0.0713±0.0241(lsq)/0.0268(boot) s_y=-0.0053±0.0241/0.0204 | a boot 95% [0.0264,0.1246] | correction to add to gen1 (dx,dy)=(+0.0713,+0.0053) | NMAD dh/tan: 0.3329 -> 0.3373 m | NMAD dh: 0.1393 -> 0.1413 m | median dh +0.0139 m, median slope 22.7, bin rms resid 0.1021
[iv CSF cache registered swath 136] n=645237 bins=36 (min bin n=6302) a=0.0399 m b=262.5 deg c=-0.0155 m | s_x=-0.0396±0.0123(lsq)/0.0147(boot) s_y=-0.0052±0.0123/0.0093 | a boot 95% [0.0158,0.0689] | correction to add to gen1 (dx,dy)=(+0.0396,+0.0052) | NMAD dh/tan: 0.3405 -> 0.3412 m | NMAD dh: 0.1438 -> 0.1441 m | median dh -0.0021 m, median slope 22.7, bin rms resid 0.0521
[iv CSF cache registered swath 137] n=676547 bins=36 (min bin n=5228) a=0.0533 m b=243.2 deg c=0.0008 m | s_x=-0.0476±0.0124(lsq)/0.0142(boot) s_y=-0.0241±0.0124/0.0106 | a boot 95% [0.0279,0.0787] | correction to add to gen1 (dx,dy)=(+0.0476,+0.0241) | NMAD dh/tan: 0.3062 -> 0.3039 m | NMAD dh: 0.1273 -> 0.1264 m | median dh +0.0001 m, median slope 21.9, bin rms resid 0.0524
[iv CSF cache registered swath 138] n=206775 bins=36 (min bin n=400) a=0.0712 m b=216.1 deg c=0.0448 m | s_x=-0.0420±0.0189(lsq)/0.0191(boot) s_y=-0.0575±0.0189/0.0181 | a boot 95% [0.0373,0.1127] | correction to add to gen1 (dx,dy)=(+0.0420,+0.0575) | NMAD dh/tan: 0.3642 -> 0.3632 m | NMAD dh: 0.1467 -> 0.1460 m | median dh +0.0045 m, median slope 21.7, bin rms resid 0.0802
[sens: i raw pooled vs z_after.npy (Z21)] n=2713717 bins=36 (min bin n=24845) a=0.1314 m b=50.4 deg c=-0.1250 m | s_x=+0.1012±0.0085(lsq)/0.0102(boot) s_y=+0.0839±0.0085/0.0060 | a boot 95% [0.1159,0.1478] | correction to add to gen1 (dx,dy)=(-0.1012,-0.0839) | NMAD dh/tan: 0.4148 -> 0.4046 m | NMAD dh: 0.1852 -> 0.1804 m | median dh -0.0501 m, median slope 24.5, bin rms resid 0.0359
[sens: iii-a registered pooled vs z_after.npy (Z21)] n=2713277 bins=36 (min bin n=24813) a=0.0285 m b=205.6 deg c=0.0450 m | s_x=-0.0123±0.0086(lsq)/0.0106(boot) s_y=-0.0257±0.0086/0.0062 | a boot 95% [0.0168,0.0458] | correction to add to gen1 (dx,dy)=(+0.0123,+0.0257) | NMAD dh/tan: 0.3917 -> 0.3916 m | NMAD dh: 0.1744 -> 0.1744 m | median dh +0.0173 m, median slope 24.5, bin rms resid 0.0365

-- FIT TABLE (m, deg). s_x,s_y = displacement of gen1 relative to gen2; correction = -s
   case | n | a | b | c | s_x(dx) | s_y(dy) | sig_sx(lsq/boot) | sig_sy(lsq/boot) | NMAD dh/tan before->after | NMAD dh before->after
   i raw pooled | 2717095 | 0.1406 | 49.0 | -0.1312 | +0.1061 | +0.0923 | 0.0091/0.0114 | 0.0091/0.0062 | 0.3811->0.3664 | 0.1714->0.1646
   ii raw swath 135 | 131627 | 0.1917 | 113.2 | -0.0815 | +0.1762 | -0.0755 | 0.0180/0.0196 | 0.0180/0.0148 | 0.4148->0.3603 | 0.1801->0.1569
   ii raw swath 136 | 958456 | 0.1642 | 59.3 | -0.1464 | +0.1412 | +0.0839 | 0.0104/0.0126 | 0.0104/0.0074 | 0.3719->0.3474 | 0.1679->0.1572
   ii raw swath 137 | 1069080 | 0.2153 | 52.7 | -0.1062 | +0.1713 | +0.1304 | 0.0123/0.0139 | 0.0123/0.0109 | 0.3759->0.3410 | 0.1665->0.1508
   ii raw swath 138 | 557932 | 0.1900 | 235.7 | -0.0561 | -0.1570 | -0.1070 | 0.0107/0.0107 | 0.0107/0.0106 | 0.3954->0.3736 | 0.1835->0.1728
   iii-a registered(star+lateral, AS APPLIED) pooled | 2716570 | 0.0194 | 204.2 | +0.0382 | -0.0079 | -0.0177 | 0.0092/0.0108 | 0.0092/0.0068 | 0.3520->0.3520 | 0.1574->0.1574
   iii-a registered(star+lateral) swath 135 | 131569 | 0.0340 | 239.0 | +0.1063 | -0.0291 | -0.0175 | 0.0184/0.0210 | 0.0184/0.0157 | 0.3587->0.3599 | 0.1564->0.1567
   iii-a registered(star+lateral) swath 136 | 958341 | 0.0308 | 231.1 | +0.0160 | -0.0239 | -0.0193 | 0.0112/0.0135 | 0.0112/0.0079 | 0.3505->0.3495 | 0.1575->0.1572
   iii-a registered(star+lateral) swath 137 | 1068913 | 0.0305 | 144.1 | +0.0454 | +0.0179 | -0.0247 | 0.0130/0.0148 | 0.0130/0.0110 | 0.3402->0.3413 | 0.1502->0.1505
   iii-a registered(star+lateral) swath 138 | 557747 | 0.0695 | 213.9 | +0.0873 | -0.0387 | -0.0577 | 0.0116/0.0117 | 0.0116/0.0109 | 0.3740->0.3757 | 0.1724->0.1729
   iii-b internal+star (HYPOTHETICAL, not applied) pooled | 2716044 | 0.5882 | 100.9 | +0.0198 | +0.5776 | -0.1111 | 0.0130/0.0138 | 0.0130/0.0117 | 0.6524->0.4222 | 0.2893->0.1881
   iii-b internal+star (hypothetical) swath 135 | 131557 | 0.0199 | 27.7 | +0.1060 | +0.0093 | +0.0177 | 0.0183/0.0206 | 0.0183/0.0155 | 0.3600->0.3597 | 0.1570->0.1568
   iii-b internal+star (hypothetical) swath 136 | 958504 | 0.3411 | 101.6 | -0.0246 | +0.3341 | -0.0687 | 0.0109/0.0135 | 0.0109/0.0077 | 0.4541->0.3485 | 0.2026->0.1573
   iii-b internal+star (hypothetical) swath 137 | 1069337 | 0.6753 | 104.1 | -0.0151 | +0.6551 | -0.1641 | 0.0127/0.0142 | 0.0127/0.0107 | 0.6870->0.3411 | 0.3034->0.1509
   iii-b internal+star (hypothetical) swath 138 | 556646 | 1.0577 | 85.2 | -0.0014 | +1.0540 | +0.0892 | 0.0100/0.0101 | 0.0100/0.0092 | 0.9590->0.3745 | 0.4682->0.1731
   iv CSF cache registered(star+lateral) pooled | 1625190 | 0.0409 | 255.2 | -0.0031 | -0.0395 | -0.0105 | 0.0104/0.0124 | 0.0104/0.0082 | 0.3285->0.3277 | 0.1368->0.1366
   iv CSF cache registered swath 135 | 96631 | 0.0715 | 265.8 | +0.0713 | -0.0713 | -0.0053 | 0.0241/0.0268 | 0.0241/0.0204 | 0.3329->0.3373 | 0.1393->0.1413
   iv CSF cache registered swath 136 | 645237 | 0.0399 | 262.5 | -0.0155 | -0.0396 | -0.0052 | 0.0123/0.0147 | 0.0123/0.0093 | 0.3405->0.3412 | 0.1438->0.1441
   iv CSF cache registered swath 137 | 676547 | 0.0533 | 243.2 | +0.0008 | -0.0476 | -0.0241 | 0.0124/0.0142 | 0.0124/0.0106 | 0.3062->0.3039 | 0.1273->0.1264
   iv CSF cache registered swath 138 | 206775 | 0.0712 | 216.1 | +0.0448 | -0.0420 | -0.0575 | 0.0189/0.0191 | 0.0189/0.0181 | 0.3642->0.3632 | 0.1467->0.1460
   sens: i raw pooled vs z_after.npy (Z21) | 2713717 | 0.1314 | 50.4 | -0.1250 | +0.1012 | +0.0839 | 0.0085/0.0102 | 0.0085/0.0060 | 0.4148->0.4046 | 0.1852->0.1804
   sens: iii-a registered pooled vs z_after.npy (Z21) | 2713277 | 0.0285 | 205.6 | +0.0450 | -0.0123 | -0.0257 | 0.0086/0.0106 | 0.0086/0.0062 | 0.3917->0.3916 | 0.1744->0.1744

=== STEP 3: median dh by slope band, pooled and by aspect quadrant (non-floodplain) ===
'after fit' = gen1 x,y moved by the fitted correction (-s_x,-s_y) and gen2 RESAMPLED (not a linear removal)
-- median dh (m) by gen2 slope band, raw, no shift
   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN
   0-5       1629113  -0.0606    -0.0622      -0.0655       -0.0558      -0.0585  478261  387849  313856  449147
   5-10      1339500  -0.0501    -0.0363      -0.0560       -0.0699      -0.0477  367127  350772  192120  429481
   10-15      811857  -0.0524    -0.0393      -0.0700       -0.0750      -0.0400  239988  198659  120687  252523
   15-20      685097  -0.0468    -0.0294      -0.0779       -0.0745      -0.0256  198112  163300  111997  211688
   20-25      751712  -0.0433    -0.0223      -0.0773       -0.0866      -0.0212  213755  143274  136730  257953
   25-30      800403  -0.0644    -0.0333      -0.0909       -0.1275      -0.0246  243065  199049  161446  196843
   30-35      394196  -0.0663    -0.0171      -0.1124       -0.1229      -0.0474  155282  115806   46960   76148
   35-45       77835  -0.0218    +0.0028      -0.0679       -0.1197      +0.0137   38452   20139    4145   15099
-- median dh (m) by gen2 slope band, registered (star+lateral), before fit
   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN
   0-5       1629131  +0.0094    +0.0022      +0.0091       +0.0216      +0.0074  478238  387903  313799  449191
   5-10      1339680  +0.0167    +0.0200      +0.0245       +0.0156      +0.0090  367060  350930  192216  429474
   10-15      811893  +0.0116    +0.0106      +0.0173       +0.0161      +0.0061  240049  198644  120654  252546
   15-20      684927  +0.0176    +0.0130      +0.0215       +0.0281      +0.0132  197956  163271  112040  211660
   20-25      751369  +0.0186    +0.0207      +0.0271       +0.0230      +0.0112  213711  143207  136681  257770
   25-30      800419  +0.0090    +0.0051      +0.0196       -0.0022      +0.0112  243077  198966  161502  196874
   30-35      394123  +0.0082    +0.0252      +0.0110       -0.0028      -0.0158  155264  115829   46913   76117
   35-45       77904  +0.0390    +0.0483      +0.0554       +0.0222      -0.0009   38481   20135    4166   15122
-- median dh (m) by gen2 slope band, registered + pooled fitted shift (+0.008,+0.018) resampled
   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN
   0-5       1629066  +0.0094    +0.0026      +0.0083       +0.0213      +0.0082  478187  387913  313811  449155
   5-10      1339680  +0.0170    +0.0208      +0.0227       +0.0148      +0.0111  367094  350917  192152  429517
   10-15      811936  +0.0122    +0.0119      +0.0141       +0.0147      +0.0096  240084  198639  120631  252582
   15-20      684913  +0.0183    +0.0149      +0.0166       +0.0263      +0.0181  197942  163291  112062  211618
   20-25      751346  +0.0201    +0.0232      +0.0209       +0.0209      +0.0175  213714  143200  136636  257796
   25-30      800488  +0.0089    +0.0076      +0.0118       -0.0039      +0.0183  243083  198986  161533  196886
   30-35      394093  +0.0070    +0.0263      +0.0019       -0.0046      -0.0071  155275  115803   46896   76119
   35-45       77919  +0.0383    +0.0486      +0.0451       +0.0166      +0.0131   38482   20151    4165   15121
-- median dh (m) by gen2 slope band, registered + pooled fit removed LINEARLY (check)
   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN
   0-5       1629131  +0.0094    +0.0026      +0.0083       +0.0213      +0.0082  478238  387903  313799  449191
   5-10      1339680  +0.0170    +0.0208      +0.0226       +0.0147      +0.0111  367060  350930  192216  429474
   10-15      811893  +0.0122    +0.0119      +0.0141       +0.0147      +0.0096  240049  198644  120654  252546
   15-20      684927  +0.0183    +0.0149      +0.0167       +0.0264      +0.0181  197956  163271  112040  211660
   20-25      751369  +0.0201    +0.0232      +0.0210       +0.0210      +0.0174  213711  143207  136681  257770
   25-30      800419  +0.0089    +0.0075      +0.0119       -0.0039      +0.0181  243077  198966  161502  196874
   30-35      394123  +0.0070    +0.0261      +0.0021       -0.0045      -0.0073  155264  115829   46913   76117
   35-45       77904  +0.0384    +0.0485      +0.0455       +0.0177      +0.0125   38481   20135    4166   15122
-- median dh (m) by gen2 slope band, registered + PER-SWATH fitted shifts resampled
   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN
   0-5       1629076  +0.0095    +0.0028      +0.0080       +0.0211      +0.0085  478186  387906  313837  449147
   5-10      1339713  +0.0171    +0.0212      +0.0219       +0.0144      +0.0119  367095  350890  192198  429530
   10-15      811989  +0.0126    +0.0126      +0.0126       +0.0145      +0.0118  240027  198694  120625  252643
   15-20      684780  +0.0190    +0.0160      +0.0144       +0.0262      +0.0210  197965  163226  112045  211544
   20-25      751417  +0.0213    +0.0253      +0.0178       +0.0200      +0.0212  213677  143206  136710  257824
   25-30      800518  +0.0095    +0.0105      +0.0083       -0.0069      +0.0235  243064  199026  161537  196891
   30-35      394020  +0.0073    +0.0267      -0.0031       -0.0070      +0.0019  155256  115779   46879   76106
   35-45       77932  +0.0368    +0.0465      +0.0367       +0.0158      +0.0196   38494   20170    4156   15112
   raw, no shift: pooled median dh 0-5 deg -0.0606; minus that: 0-5 +0.0000, 5-10 +0.0105, 10-15 +0.0083, 15-20 +0.0138, 20-25 +0.0174, 25-30 -0.0038, 30-35 -0.0057, 35-45 +0.0388
      E 45-135: 0-5 +0.0000, 5-10 +0.0259, 10-15 +0.0229, 15-20 +0.0329, 20-25 +0.0399, 25-30 +0.0290, 30-35 +0.0452, 35-45 +0.0651
      S 135-225: 0-5 +0.0000, 5-10 +0.0095, 10-15 -0.0046, 15-20 -0.0124, 20-25 -0.0119, 25-30 -0.0255, 30-35 -0.0469, 35-45 -0.0024
      W 225-315: 0-5 +0.0000, 5-10 -0.0141, 10-15 -0.0192, 15-20 -0.0187, 20-25 -0.0308, 25-30 -0.0717, 30-35 -0.0671, 35-45 -0.0639
      N 315-45: 0-5 +0.0000, 5-10 +0.0108, 10-15 +0.0185, 15-20 +0.0329, 20-25 +0.0373, 25-30 +0.0339, 30-35 +0.0111, 35-45 +0.0722
   registered (star+lateral), before fit: pooled median dh 0-5 deg +0.0094; minus that: 0-5 +0.0000, 5-10 +0.0073, 10-15 +0.0023, 15-20 +0.0083, 20-25 +0.0093, 25-30 -0.0004, 30-35 -0.0012, 35-45 +0.0296
      E 45-135: 0-5 +0.0000, 5-10 +0.0178, 10-15 +0.0084, 15-20 +0.0108, 20-25 +0.0185, 25-30 +0.0029, 30-35 +0.0230, 35-45 +0.0461
      S 135-225: 0-5 +0.0000, 5-10 +0.0155, 10-15 +0.0083, 15-20 +0.0125, 20-25 +0.0180, 25-30 +0.0105, 30-35 +0.0020, 35-45 +0.0464
      W 225-315: 0-5 +0.0000, 5-10 -0.0059, 10-15 -0.0055, 15-20 +0.0065, 20-25 +0.0014, 25-30 -0.0237, 30-35 -0.0244, 35-45 +0.0007
      N 315-45: 0-5 +0.0000, 5-10 +0.0016, 10-15 -0.0013, 15-20 +0.0059, 20-25 +0.0038, 25-30 +0.0038, 30-35 -0.0232, 35-45 -0.0083
   registered + pooled fitted shift (+0.008,+0.018) resampled: pooled median dh 0-5 deg +0.0094; minus that: 0-5 +0.0000, 5-10 +0.0076, 10-15 +0.0028, 15-20 +0.0089, 20-25 +0.0107, 25-30 -0.0005, 30-35 -0.0024, 35-45 +0.0289
      E 45-135: 0-5 +0.0000, 5-10 +0.0183, 10-15 +0.0094, 15-20 +0.0124, 20-25 +0.0207, 25-30 +0.0050, 30-35 +0.0237, 35-45 +0.0461
      S 135-225: 0-5 +0.0000, 5-10 +0.0144, 10-15 +0.0058, 15-20 +0.0083, 20-25 +0.0126, 25-30 +0.0035, 30-35 -0.0064, 35-45 +0.0368
      W 225-315: 0-5 +0.0000, 5-10 -0.0065, 10-15 -0.0066, 15-20 +0.0050, 20-25 -0.0003, 25-30 -0.0252, 30-35 -0.0258, 35-45 -0.0046
      N 315-45: 0-5 +0.0000, 5-10 +0.0030, 10-15 +0.0015, 15-20 +0.0099, 20-25 +0.0093, 25-30 +0.0101, 30-35 -0.0153, 35-45 +0.0050
   registered + pooled fit removed LINEARLY (check): pooled median dh 0-5 deg +0.0094; minus that: 0-5 +0.0000, 5-10 +0.0076, 10-15 +0.0028, 15-20 +0.0089, 20-25 +0.0107, 25-30 -0.0005, 30-35 -0.0024, 35-45 +0.0290
      E 45-135: 0-5 +0.0000, 5-10 +0.0183, 10-15 +0.0094, 15-20 +0.0123, 20-25 +0.0206, 25-30 +0.0050, 30-35 +0.0235, 35-45 +0.0459
      S 135-225: 0-5 +0.0000, 5-10 +0.0143, 10-15 +0.0058, 15-20 +0.0084, 20-25 +0.0127, 25-30 +0.0036, 30-35 -0.0062, 35-45 +0.0372
      W 225-315: 0-5 +0.0000, 5-10 -0.0066, 10-15 -0.0066, 15-20 +0.0051, 20-25 -0.0003, 25-30 -0.0252, 30-35 -0.0257, 35-45 -0.0036
      N 315-45: 0-5 +0.0000, 5-10 +0.0030, 10-15 +0.0015, 15-20 +0.0099, 20-25 +0.0093, 25-30 +0.0099, 30-35 -0.0155, 35-45 +0.0043
   registered + PER-SWATH fitted shifts resampled: pooled median dh 0-5 deg +0.0095; minus that: 0-5 +0.0000, 5-10 +0.0077, 10-15 +0.0031, 15-20 +0.0096, 20-25 +0.0118, 25-30 +0.0001, 30-35 -0.0022, 35-45 +0.0273
      E 45-135: 0-5 +0.0000, 5-10 +0.0184, 10-15 +0.0098, 15-20 +0.0132, 20-25 +0.0226, 25-30 +0.0077, 30-35 +0.0239, 35-45 +0.0437
      S 135-225: 0-5 +0.0000, 5-10 +0.0139, 10-15 +0.0045, 15-20 +0.0064, 20-25 +0.0098, 25-30 +0.0003, 30-35 -0.0112, 35-45 +0.0286
      W 225-315: 0-5 +0.0000, 5-10 -0.0066, 10-15 -0.0065, 15-20 +0.0051, 20-25 -0.0010, 25-30 -0.0280, 30-35 -0.0280, 35-45 -0.0053
      N 315-45: 0-5 +0.0000, 5-10 +0.0034, 10-15 +0.0033, 15-20 +0.0125, 20-25 +0.0127, 25-30 +0.0150, 30-35 -0.0066, 35-45 +0.0111

=== STEP 4: the bluff window rows 483:550, cols 440:469 ===
window x [579695.0, 579840.0) y [4885155.0, 4885490.0) = 67x29 cells
  gen2 slope in window: median 28.6 deg, p10/p90 13.7/34.7; aspect median 116.6 deg; floodplain cells 28 of 1943
  cells with slope > 25 deg: 1181, their median aspect 116.1 deg
  gen1 raw last returns in window: 30494; per swath {137: 1492, 138: 29002}
  CSF-cache points in window: 5708; per swath {137: 189, 138: 5519}
  window cells (of 1943) holding at least one raw last return, per swath: {135: 0, 136: 0, 137: 190, 138: 1943}
    swath 135: raw last returns on window cells steeper than 25 deg: 0
    swath 136: raw last returns on window cells steeper than 25 deg: 0
    swath 137: raw last returns on window cells steeper than 25 deg: 700
    swath 138: raw last returns on window cells steeper than 25 deg: 17875
-- window median dh (m) by gen2 slope band, raw, no shift:  0-5: -0.1233 (n=1029) | 5-10: -0.1433 (n=683) | 10-15: -0.1390 (n=1131) | 15-20: -0.1488 (n=4142) | 20-25: -0.1786 (n=4574) | 25-30: -0.2053 (n=6515) | 30-35: -0.2165 (n=9486) | 35-45: -0.2366 (n=2548)
-- window median dh (m) by gen2 slope band, registered (star+lateral):  0-5: -0.0580 (n=1035) | 5-10: -0.0709 (n=685) | 10-15: -0.0540 (n=1130) | 15-20: -0.0714 (n=4137) | 20-25: -0.0912 (n=4592) | 25-30: -0.0990 (n=6484) | 30-35: -0.1032 (n=9501) | 35-45: -0.1022 (n=2544)
-- window median dh (m) by gen2 slope band, registered + pooled fitted shift (+0.008,+0.018):  0-5: -0.0581 (n=1035) | 5-10: -0.0711 (n=685) | 10-15: -0.0526 (n=1130) | 15-20: -0.0694 (n=4136) | 20-25: -0.0899 (n=4590) | 25-30: -0.0991 (n=6485) | 30-35: -0.1045 (n=9500) | 35-45: -0.1005 (n=2547)
-- window median dh (m) by gen2 slope band, registered + per-swath fitted shifts:  0-5: -0.0580 (n=1034) | 5-10: -0.0697 (n=684) | 10-15: -0.0481 (n=1131) | 15-20: -0.0642 (n=4141) | 20-25: -0.0865 (n=4583) | 25-30: -0.0936 (n=6482) | 30-35: -0.0999 (n=9501) | 35-45: -0.0894 (n=2550)
-- window median dh (m) by gen2 slope band, registered, swath 137 only:  0-5: +nan (n=0) | 5-10: +0.1404 (n=4) | 10-15: -0.0033 (n=8) | 15-20: +0.0063 (n=135) | 20-25: -0.1122 (n=662) | 25-30: -0.1556 (n=440) | 30-35: -0.2249 (n=243) | 35-45: +nan (n=0)
-- window median dh (m) by gen2 slope band, registered, swath 138 only:  0-5: -0.0580 (n=1035) | 5-10: -0.0733 (n=681) | 10-15: -0.0541 (n=1122) | 15-20: -0.0735 (n=4002) | 20-25: -0.0885 (n=3930) | 25-30: -0.0974 (n=6044) | 30-35: -0.1017 (n=9258) | 35-45: -0.1022 (n=2544)
-- window median dh (m) by gen2 slope band, raw + raw-pooled fitted shift (-0.106,-0.092):  0-5: -0.1241 (n=1024) | 5-10: -0.1504 (n=689) | 10-15: -0.1558 (n=1133) | 15-20: -0.1646 (n=4139) | 20-25: -0.2007 (n=4574) | 25-30: -0.2357 (n=6512) | 30-35: -0.2490 (n=9491) | 35-45: -0.2888 (n=2548)
  window, registered, slope > 25 deg: n=18529, median dh -0.1015 m, median aspect 116.0 deg, median dh/tan(slope) -0.1671 m (= implied s_x if purely horizontal and east-facing; correction would be +0.1671 m east)
  window, registered, slope < 10 deg: n=1720, median dh -0.0611 m (the local vertical offset)
-- window DoD (gen2 - gen1, m): pipeline dod.npy vs the per-cell proxy (median over raw last returns of gen2_bilinear - gen1)
   pipeline dod.npy: n finite 906, median +0.1461, steep(>25) cells median +0.2473 (n 290), gentle(<10) median +0.0556
   proxy raw, no shift: n finite 1915, median +0.1945, steep(>25) median +0.2155 (n 1181), gentle(<10) median +0.1356, steep minus gentle +0.0799
   proxy registered (star+lateral): n finite 1916, median +0.0928, steep(>25) median +0.1023 (n 1181), gentle(<10) median +0.0684, steep minus gentle +0.0339
   proxy registered + pooled fitted shift: n finite 1916, median +0.0925, steep(>25) median +0.1024 (n 1181), gentle(<10) median +0.0683, steep minus gentle +0.0341
   proxy registered + per-swath fitted shifts: n finite 1916, median +0.0882, steep(>25) median +0.0969 (n 1181), gentle(<10) median +0.0679, steep minus gentle +0.0290
   on the 290 steep cells where BOTH pipeline dod and registered proxy are finite: pipeline median +0.2473, proxy median +0.1239, proxy after pooled fit +0.1247, after per-swath fit +0.1195

=== STEP 5: the deciding numbers ===
registered pooled fit amplitude a = 0.0194 m (boot 95% [0.0082, 0.0378]), direction b = 204.2 deg, i.e. gen1 displaced (-0.0079, -0.0177) m relative to gen2
a horizontal shift explaining the window's steep offset alone would need |s_x| ~ 0.167 m
window steep median dh, registered: -0.1015; after pooled fitted shift: -0.1021; after per-swath shifts: -0.0967
done 2026-11-11 16:43:50
```

## Addendum: derived deciding numbers (computed from nuth_kaab_fits.csv and the window outputs)

```
registered pooled: s=(-0.0079,-0.0177) m; component along the window downslope (aspect 116 deg) +0.0006 m; predicted gen1-minus-gen2 on a 30-deg slope facing 116 deg = tan(30)*that = +0.0004 m (sigma from boot: 0.0059 m)
registered swath 138 (the bluff swath): s=(-0.0387,-0.0577) m; component along the window downslope (aspect 116 deg) -0.0095 m; predicted gen1-minus-gen2 on a 30-deg slope facing 116 deg = tan(30)*that = -0.0055 m (sigma from boot: 0.0067 m)
CSF registered pooled: s=(-0.0395,-0.0105) m; component along the window downslope (aspect 116 deg) -0.0309 m; predicted gen1-minus-gen2 on a 30-deg slope facing 116 deg = tan(30)*that = -0.0178 m (sigma from boot: 0.0068 m)
horizontal shift along 116 deg that would be needed for dh = -0.40 m on a 30-deg slope: 0.693 m
  ... for the pipeline dod.npy steep-minus-gentle window offset (+0.2473 - +0.0556 = +0.1917 m) at the window median slope 28.6 deg: 0.352 m
  ... for the raw-last-return proxy steep-minus-gentle offset, registered (+0.0339 m) at 28.6 deg: 0.062 m
ratio: needed shift for -0.40 m / fitted residual amplitude a (registered pooled 0.0194 m, boot 97.5% 0.0378): 36x the fit, 18x its upper bound
independent check of the pipeline's APPLIED horizontal registration (star+lateral) against -s from the RAW per-swath Nuth & Kaab (case ii):
  swath 135: pipeline applied (dx,dy)=(-0.2068,+0.0460); raw N-K correction (-s_x,-s_y)=(-0.1762,+0.0755); difference (-0.0306,-0.0295) m
  swath 136: pipeline applied (dx,dy)=(-0.1654,-0.1071); raw N-K correction (-s_x,-s_y)=(-0.1412,-0.0839); difference (-0.0242,-0.0232) m
  swath 137: pipeline applied (dx,dy)=(-0.1479,-0.1603); raw N-K correction (-s_x,-s_y)=(-0.1713,-0.1304); difference (+0.0234,-0.0299) m
  swath 138: pipeline applied (dx,dy)=(+0.1206,+0.0502); raw N-K correction (-s_x,-s_y)=(+0.1570,+0.1070); difference (-0.0364,-0.0568) m
```
