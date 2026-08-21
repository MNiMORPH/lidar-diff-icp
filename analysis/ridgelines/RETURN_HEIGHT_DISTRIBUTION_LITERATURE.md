# Lidar return-height distributions — literature synthesis

Agent literature search 2026-08-22 (Andy-requested), mapping our empirical return-height
structure onto established frameworks. Citations verified to title/venue/DOI (most findings
from abstracts, not full text; several primaries paywalled). Companion to
[bare-ground skew synthesis](BARE_GROUND_SKEW_LITERATURE.md).

## Our findings, restated
1. **Layered all-returns profile:** ground spike (d=0) + low-veg/understory shoulder
   (0.1–1 m) + local minimum (~1.5 m) + canopy body to ~30 m (gen1) / ~36 m (gen2).
   gen2 leaf-on dense canopy; gen1 leaf-off 1–2 orders sparser, thinning upward.
2. **Ground return = Gaussian core + tails, form set by LAND COVER:**
   FOREST → symmetric EXPONENTIAL tails (λ~76–106 mm); OPEN/BARE → POWER-LAW tails
   (α~2.6–3.8, upper-heavier in gen1 = Nov stubble). Form conserved across epochs; epoch
   sets scale (gen2 ~2× tighter) + skew.

## (a) Canopy / vertical-profile frameworks — our layered profile is EXPLAINED (gross shape)
- **MacArthur & Horn (1969, *Ecology* 50:802–804):** foliage-profile transform; Poisson gap
  probability → foliage density ∝ d/dz[−ln(gap)]. Beer–Lambert/Poisson extinction.
- **Ni-Meister, Jupp & Dubayah (2001, *IEEE TGRS* 39(9):1943–1958):** GORT waveform model;
  **return density decays EXPONENTIALLY with canopy depth**, P_gap(z)=exp(−G·LAI(z)/cosθ).
- **Harding, Lefsky, Parker & Blair (2001, *RSE* 76:283–297):** canopy height profiles,
  apparent-vs-true (occlusion) correction. Also Lefsky 2002; Drake 2002; Blair & Hofton 1999.
- **Verdict:** gross layered shape (i) PREDICTED; the explicit ~1.5 m understory/canopy
  minimum (ii) consistent-but-not-directly-documented. No paper models the *whole* histogram
  in closed form — it's assembled piecewise.

## (b) Ground-peak shape/width/tails — Gaussian core EXPLAINED; tail forms NOT
- **Wagner et al. (2006, *ISPRS J* 60:100–112); Hofton et al. (2000, *IEEE TGRS* 38:1989);
  Mallet & Bretar (2009, *ISPRS J* 64:1–16):** ground echo = Gaussian; width = pulse ⊗ range
  spread of footprint scatterers.
- **Gaussian width scaling:** footprint×slope broadening ∝ beamwidth·tanθ (**Gardner 1992,
  *IEEE TGRS* 30:1061**) + roughness + unresolved near-ground veg. Matches our tighter-open /
  wider-forest, tighter-gen2 pattern.
- Non-Gaussian departures documented only as **skew/broadening** → generalized-Gaussian,
  lognormal, Nakagami, Burr echo models (**Chauve 2007/2009; Mallet et al. 2011**). Understory
  merges into ground echo (**Crespo-Peremarch 2018, *RSE* 217:400–413**); discrete-return
  dead-zone ~0.6–1 m (**Ussyshkin & Theriault 2011, *Remote Sensing* 3:416–434**).
- **Verdict:** Gaussian core (i) PREDICTED; specific tail FORMS NOT documented.

## (c) Vertical-error distribution / tail form — the novelty
- **Höhle & Höhle (2009, *ISPRS J* 64:398–406):** bare-earth error is non-normal, skewed,
  heavy-tailed, outlier-prone → use robust stats (median, NMAD, quantiles). Fits NO tail.
- **Only named heavy-tail fits land on LAPLACE (double-exponential):** **Zandbergen (2011,
  *IJRS* 32:409–430)** and **Hejmanowska & Kay (2011, *Arch. Photogramm.* 22:201–213)** —
  both pooled/mixed cover, exponential tail, no slope dependence.
- **POWER-LAW verdict (forward-citation sweep, verified):** NO peer-reviewed source fits a
  power-law/Pareto/algebraic (α≈2–4) tail to bare-ground lidar vertical error. The one
  alpha-stable hit (**Sofia, Pirotti & Tarolli 2013**) is on CURVATURE + outlier simulation,
  not the elevation-error histogram → does not overturn the null. Prior "found NONE" corroborated.

## Bottom line on our two claims
- **Exponential FOREST ground-tail:** mechanism (Beer–Lambert near-ground extinction) is
  established, but the *fitted exponential tail on the forest ground return* is undocumented —
  our physically-motivated form. Closest analog = the pooled Laplace fits (Zandbergen 2011;
  Hejmanowska & Kay 2011), NOT forest-resolved → (ii) consistent-but-stronger-than-documented.
- **POWER-LAW bare-ground tail:** (iii) APPARENTLY NOVEL. No one fits an algebraic tail to
  bare-ground error; the field's fitted heavy-tail form is Laplace, which our bare-ground data
  departs from. Defensible novelty. (Caveat: a search null isn't formal proof — a Scopus/WoS
  citation-graph query on papers citing Zandbergen 2011 + Höhle & Höhle 2009 would close it.)

## Leaf-on/off + discrete-return support
- Leaf-off → denser ground returns, better DTMs, difference concentrated near-ground
  (**Næsset 2005 *RSE* 98:356; Wasser 2013 *PLoS ONE* 8:e54776; Simpson 2017 *Remote Sens.*
  9:1101**; USGS 3DEP prefers leaf-off, QL2 ≥2 pts/m²). Residual near-ground understory/residue
  corruption = forest analog of our Nov crop-stubble upper tail.
- Discrete-return under-samples the waveform, biases to strong echoes, merges near-ground
  returns into the ground peak within a ~0.6–1 m dead zone (**Ussyshkin & Theriault 2011;
  Crespo-Peremarch 2020 *For. Ecol. Manage.* 473:118268**) → plausible why our discrete ground
  peaks carry tails rather than clean Gaussians.
