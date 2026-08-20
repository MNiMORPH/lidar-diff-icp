"""Spatial-coherence Bayesian thresholding of a DEM of Difference.

This is the method of **Wheaton, Brasington, Darby & Sear (2010)** -- "Accounting
for uncertainty in DEMs from repeat topographic surveys: improved sediment
budgets", *Earth Surface Processes and Landforms* 35(2): 136-156,
doi:10.1002/esp.1886 -- as first developed in **Wheaton (2008)**, PhD thesis, and
implemented in their **Geomorphic Change Detection (GCD)** software (Wheaton &
Brasington; https://gcd.riverscapes.net, GPL). The algorithm and its structure are
transcribed faithfully from the GCD MATLAB source (github.com/joewheaton/DoD:
``m_3SpatialCoherence.m`` and ``m_3NeighbourhoodClass.m``). **Credit for the method
is Wheaton et al.'s**; only the default weight-threshold calibration noted below is
ours (GCD leaves those user-defined).

The idea (Wheaton et al., 2010): a per-cell probability of real change, from the
propagated-error t-score, is Bayesian-updated by the SPATIAL COHERENCE of its
neighbourhood. A cell embedded in same-sign change is promoted -- so a broad,
low-amplitude coherent patch can pass even where no single cell exceeds its LoD --
while an isolated cell is demoted, so blunders (and the smoothing 'halos' of naive
amplitude-pooling) do not survive. Crucially the neighbourhood evidence is a
COUNT of same-sign cells in a 5x5 window, NOT a smoothing of amplitude, which is
why it localises rather than haloing.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import convolve


def isotropic_counts(dod, valid, window=5):
    """Wheaton's default neighbourhood evidence: the COUNT of same-sign cells in an
    isotropic ``window`` x ``window`` square (GCD m_3NeighbourhoodClass.m). Returns
    ``(ndepos, neros, n)`` -- deposition/erosion same-sign counts and the window's
    max count ``n`` (used to calibrate the coherence-weight bounds)."""
    k = np.ones((window, window))
    ndepos = convolve((valid & (dod > 0)).astype(float), k, mode="constant")
    neros = convolve((valid & (dod < 0)).astype(float), k, mode="constant")
    return ndepos, neros, window * window


def spatial_coherence_probability(dod, perror, *, counts=None, window=5, low=None,
                                  up=None, dof=1000):
    """Wheaton et al. (2010) posterior probability of real change, signed in
    [-1 (erosion) .. +1 (deposition)]. Threshold ``|posterior|`` at a confidence
    level (e.g. 0.95) for a change mask (see :func:`coherence_change`).

    ``dod``   : DEM of Difference (gen2 - gen1), m.
    ``perror``: propagated per-cell 1-sigma error sqrt(sigma_gen1^2 + sigma_gen2^2)
                (m); equals ``lod / 1.96`` for our heteroscedastic LoD.

    The spatial-coherence NEIGHBOURHOOD is pluggable. ``counts`` = ``(ndepos, neros,
    n)`` supplies the same-sign counts from ANY neighbourhood -- the isotropic square
    (default, :func:`isotropic_counts`) or a flow corridor
    (:func:`lidar_diff_icp.areas.flow_corridor_counts`).
    Everything else -- the amplitude prior, the count->weight mapping, and the
    Bayesian posterior -- is Wheaton's, unchanged; only the neighbourhood geometry
    varies. ``window``/``low``/``up`` set the isotropic default and the weight bounds.
    """
    from scipy.stats import t as tdist
    dod = np.asarray(dod, float); perror = np.asarray(perror, float)
    valid = np.isfinite(dod) & np.isfinite(perror) & (perror > 0)

    # Prior probability of real change from the propagated-error t-score
    # (Wheaton et al. 2010, GCD: tscore = DoD / perror; prior = |2*p - 1|).
    tscore = np.where(valid, dod / np.where(perror > 0, perror, 1.0), 0.0)
    priorp = np.abs(2.0 * tdist.cdf(tscore, dof) - 1.0)

    # Neighbourhood COUNT of same-sign cells (the spatial-coherence evidence),
    # default isotropic square; pass `counts` for a flow-aligned neighbourhood.
    ndepos, neros, n = counts if counts is not None else isotropic_counts(dod, valid, window)

    # Coherence weights = P(neighbourhood | change): linearly rescale counts
    # between low..up (GCD leaves these user-defined). Default calibration (ours):
    # a window of n cells has same-sign count ~Binomial(n, 0.5) ~ n/2 +/- sqrt(n)/2
    # under random signs, so weight 0 up to ~1 sigma above chance and 1 at strong
    # agreement (~7/8 of the window).
    if low is None:
        low = 0.80 * n                            # 20 for 5x5: need >=80% same-sign agreement to boost
    if up is None:
        up = 0.96 * n                             # 24 for 5x5: weight 1 at near-total agreement
    wd = np.clip((ndepos - low) / (up - low), 0.0, 1.0)
    we = np.clip((neros - low) / (up - low), 0.0, 1.0)

    # Bayesian posterior (Wheaton et al. 2010): postp = (prior*w) /
    # (prior*w + (1-prior)*(1-w)), signed by erosion/deposition.
    def bayes(pr, w):
        return (pr * w) / np.maximum(pr * w + (1.0 - pr) * (1.0 - w), 1e-12)

    post = np.full(dod.shape, np.nan)
    dep = valid & (dod > 0); ero = valid & (dod < 0)
    post[dep] = bayes(priorp[dep], wd[dep])
    post[ero] = -bayes(priorp[ero], we[ero])
    return post


def coherence_change(dod, perror, *, conf=0.95, **kw):
    """Boolean change mask: ``|Wheaton posterior| > conf`` (default 95%)."""
    post = spatial_coherence_probability(dod, perror, **kw)
    return np.isfinite(post) & (np.abs(post) > conf)


def flow_coherence_change(dod, perror, flowdown, flowup, *, k=12, width=0,
                          seed_z=1.96, conf=0.95, low=None, up=None, dof=1000):
    """Boolean change mask from Wheaton coherence over a FLOW-CORRIDOR neighbourhood
    -- the flow-aligned counterpart of the isotropic :func:`coherence_change`, for
    gullies/rills/channels that the square window suppresses.

    Two decoupled stages (see :mod:`lidar_diff_icp.areas`): (1) build the corridor
    footprint along flow with :func:`~lidar_diff_icp.areas.flow_corridor_counts`,
    seeded at the amplitude-significant cells (``|DoD/perror| > seed_z``) so flow
    connectivity -- not the isotropic count -- supplies the coherence; (2) score each
    footprint with the unchanged Wheaton count statistic
    (:func:`spatial_coherence_probability`) and threshold ``|posterior| > conf``.
    Only seeds can be flagged: flow coherence promotes a genuinely-significant cell
    embedded in a coherent flow line, it does not manufacture significance from noise.

    For patches use isotropic :func:`coherence_change`; OR the two masks for
    patches-plus-gullies. ``flowdown``/``flowup`` (per-cell downstream / dominant-
    upstream flat indices, -1 = none) come from a DEM routing (RichDEM D8/
    D-infinity), passed in so this stays routing-agnostic. ``k`` = corridor half-
    length (2k+1 long); ``width`` = lateral cells each side (channels wider than 1).
    """
    from . import areas
    dod = np.asarray(dod, float); perror = np.asarray(perror, float)
    valid = np.isfinite(dod) & np.isfinite(perror) & (perror > 0)
    tscore = np.where(valid, dod / np.where(perror > 0, perror, 1.0), 0.0)
    seed = valid & (np.abs(tscore) > seed_z)
    counts = areas.flow_corridor_counts(dod, valid, flowdown, flowup, k=k, width=width,
                                        cells=seed)
    post = spatial_coherence_probability(dod, perror, counts=counts, low=low, up=up, dof=dof)
    return seed & np.isfinite(post) & (np.abs(post) > conf)


def ridge_change(dod, perror, *, sigmas=(1, 2, 3), ridge_thresh=0.35, amp_z=2.0):
    """LINEAR change (gullies, rills, channel incision, levees) via a ridge /
    vesselness filter -- the necessary complement to Wheaton spatial coherence,
    which is isotropic and suppresses narrow linear features (a 1-cell-wide gully
    has too few same-sign neighbours in a 5x5 window, so its posterior collapses
    regardless of amplitude).

    A ridge filter is purpose-built for elongated structures: it responds to a
    line of same-sign change of *any* orientation while staying quiet on isolated
    noise. We apply the **Sato et al. (1998)** multi-scale tubular-structure filter
    (``skimage.filters.sato``) to the per-sign t-score ``DoD / perror``, separately
    for erosion and deposition, then keep cells with a strong ridge response
    (> ``ridge_thresh`` of the per-sign max) that also clear an amplitude gate
    (``|t| > amp_z``) -- the gate is what stops a chance-aligned noise streak from
    lighting up. Returns a boolean line-change mask; OR it with
    :func:`coherence_change` for a detector that keeps both patches and lines.

    Reference: Sato, Nakajima, Shiraga, Atsumi, Yoshida, Koller, Gerig & Kikinis
    (1998), "Three-dimensional multi-scale line filter for segmentation and
    visualization of curvilinear structures in medical images", Medical Image
    Analysis 2(2): 143-168, doi:10.1016/S1361-8415(98)80009-1. The related Frangi
    et al. (1998) vesselness filter (``skimage.filters.frangi``) is a drop-in
    alternative. Ridge/valley detection by directional/Hessian filters is the
    established way to extract linear terrain features (e.g. difference-of-rotating-
    Gaussian ridge detection; Hessian structure-tensor methods).
    """
    from skimage.filters import sato
    valid = np.isfinite(dod) & np.isfinite(perror) & (perror > 0)
    sig = np.nan_to_num(np.where(valid, dod / np.maximum(perror, 1e-9), 0.0))
    ero = sato(np.clip(-sig, 0, None), sigmas=sigmas, black_ridges=False)
    dep = sato(np.clip(sig, 0, None), sigmas=sigmas, black_ridges=False)
    ero = ero / (ero.max() + 1e-9); dep = dep / (dep.max() + 1e-9)
    return valid & (((ero > ridge_thresh) & (sig < -amp_z)) |
                    ((dep > ridge_thresh) & (sig > amp_z)))
