"""Two-axis (amplitude x extent) change detection for DEMs of Difference.

Change is real only if it clears BOTH of two independent axes, because they guard
against two different errors that neither can cover for the other:

* **Extent (N_eff)** -- guards against *random, spatially-correlated* noise. A
  region's area-mean uncertainty from random error is ``sigma / sqrt(N_eff)``,
  ``N_eff = area / (pi * L^2)`` (Rolstad 2009; Hugonnet 2022), so a large coherent
  region can be significant where no single cell exceeds its per-cell LoD.
* **Amplitude (systematic floor tau_sys)** -- guards against *systematic* error
  (residual co-registration warp, datum/geoid, penetration bias). It does NOT
  shrink with area, so N_eff-significance alone would bless an arbitrarily small
  mean over a big enough region. The per-cell LoD / stable-sigma scale sets this
  floor; N_eff cannot buy past it.

Combined (the standard spatially-correlated area error) into one signal-to-noise:

    sigma_mean = sqrt( sigma^2 / N_eff  +  tau_sys^2 ) ;   SNR = |mean| / sigma_mean

For a large region SNR is amplitude-limited (tau_sys); for a small one,
extent-limited (sigma/sqrt(N_eff)). Per-cell LoD is NOT demoted: it IS the
amplitude floor, and it also describes the confident *core* vs the
coherence-included *fringe* within a confirmed feature.

The procedure is **propose -> confirm**, and it does NOT hard-gate cells at the
per-cell LoD first (that truncates the low-amplitude fringe of a real feature and
misses faint broad features entirely). It proposes regions by a multi-scale
standardized statistic -- a sharp scarp seeds at the pixel scale, a broad thin
sheet (and the low-amplitude coherent *margin* of a big feature) seeds at a coarse
scale where smoothing beats the noise down -- and takes same-sign connected
components. No percolating grow (which would merge a feature into surrounding
same-sign noise and dilute it). It then confirms each region by SNR against a
**null-calibrated** threshold: the same propose-and-score procedure is run on
Gaussian-random-field nulls with the data's own error correlation length, so the
selection bias of "score the cluster you found" is present in the null too and
cancels. A region is real change iff its SNR exceeds the null's field-wise 95th
percentile.

The DoD is whitened by the per-cell error (lod/1.96) before seeding, so
heteroscedastic error (larger on steep/rough ground) is handled rather than
compared to one global sigma.

Known limitation -- STEEP TERRAIN. Where the LoD model under-predicts the most
extreme steep-cell error (isolated blunders on bluffs), the coarse smoothing scale
can spread such a cell into a small circular seed and over-detect. On the flat
floodplains this is a non-issue; on dissected bluffs, treat small low-amplitude
detections with suspicion (they are a mix of real steep-slope change and this
artifact). A per-cell blunder filter or a slope-aware LoD refinement would close it.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, label


def _corr_length(z, mask, res, max_lag_m=400.0):
    """1/e correlation length (m) of a masked field via radial autocorrelation."""
    H, W = z.shape
    w = mask.astype(float)
    zm = np.where(mask, z - z[mask].mean(), 0.0)
    num = np.fft.ifft2(np.abs(np.fft.fft2(zm)) ** 2).real
    den = np.maximum(np.fft.ifft2(np.abs(np.fft.fft2(w)) ** 2).real, 1e-6)
    ac = np.fft.fftshift(num / den)
    cy, cx = H // 2, W // 2
    ac /= ac[cy, cx]
    yy, xx = np.mgrid[:H, :W]
    rr = np.hypot(xx - cx, yy - cy) * res
    e = np.arange(0, max_lag_m, res)
    prof = np.array([ac[(rr >= e[i]) & (rr < e[i + 1])].mean() for i in range(len(e) - 1)])
    ctr = 0.5 * (e[1:] + e[:-1]); below = prof < 1 / np.e
    return float(ctr[np.argmax(below)]) if below.any() else res


def _systematic_floor(dod, stable, res, block_m):
    """tau_sys: the irreducible coherent-error floor = NMAD of large-block
    stable-ground DoD means (the part that does NOT average away)."""
    b = max(int(round(block_m / res)), 1); H, W = dod.shape; vals = []
    for i in range(0, H - b + 1, b):
        for j in range(0, W - b + 1, b):
            sub = dod[i:i + b, j:j + b]; m = stable[i:i + b, j:j + b] & np.isfinite(sub)
            if m.sum() > 0.5 * b * b:
                vals.append(np.mean(sub[m]))
    vals = np.array(vals)
    return float(1.4826 * np.median(np.abs(vals - np.median(vals)))) if vals.size >= 5 else 0.0


def _masked_smooth(a, mask, sig):
    if sig <= 0:
        return np.where(mask, a, np.nan)
    num = gaussian_filter(np.where(mask, a, 0.0), sig)
    den = gaussian_filter(mask.astype(float), sig)
    return np.where(mask, num / np.maximum(den, 1e-6), np.nan)


def _extract(dod, sig_cell, mm, L, tau, res, kernels, r_s, t_seed, min_cells):
    """Multi-scale significant seed -> same-sign connected regions -> per-region
    SNR. Heteroscedastic: the DoD is WHITENED by the per-cell error ``sig_cell``
    (= lod/1.96) before smoothing, so a cell is scored against its OWN error, not
    a global sigma -- otherwise noisy steep ground over-seeds. The coarse scales
    make a broad, low-amplitude coherent margin seed on its own (no percolating
    grow, which would merge features into same-sign noise).
    Returns (labels, [(id, cellmask, area, mean, neff, sigma_mean, snr), ...])."""
    H, W = dod.shape; ca = res * res
    white = np.where(mm, dod / np.maximum(sig_cell, 1e-6), np.nan)   # unit-variance noise
    T = np.zeros((H, W)); S = np.zeros((H, W))          # best standardized |Z| and its signed smoothed value
    for k, sk in enumerate(kernels):
        Ws = _masked_smooth(white, mm, sk)
        Zs = np.nan_to_num(Ws) / max(r_s[k], 1e-9)
        upd = np.abs(Zs) > T
        T = np.where(upd, np.abs(Zs), T); S = np.where(upd, np.nan_to_num(Ws), S)
    seed = mm & (T > t_seed)
    labels = np.zeros((H, W), int); out = []; nid = 0
    for sgn in (+1, -1):                                 # same-sign regions only
        lab, n = label(seed & (np.sign(S) == sgn), structure=np.ones((3, 3)))
        for i in range(1, n + 1):
            cm = lab == i; a = int(cm.sum())
            if a < min_cells:                            # isolated-spike / blunder guard
                continue
            mean = float(dod[cm].mean())
            neff = max(a * ca / (np.pi * max(L, res) ** 2), 1.0)
            sig_rand = np.sqrt(np.mean(sig_cell[cm] ** 2) / neff)     # region-local random error
            sm_ = np.sqrt(sig_rand ** 2 + tau ** 2)
            nid += 1; labels[cm] = nid
            out.append((nid, cm, a, mean, neff, sm_, abs(mean) / sm_))
    return labels, out


def detect_change(dod, lod, stable, res, *, mask=None, z_field=95.0,
                  scales_m=(0.0, 30.0, 90.0), t_seed=2.5, min_cells=4,
                  sys_block_m=150.0, nsim=80, seed=0):
    """Two-axis, null-calibrated change detection. Returns dict:
    labels (int, 0=bg), regions (list of per-region dicts, largest first),
    change (bool mask), and calibration (sigma, corr_length_m, tau_sys_m,
    snr_thresh)."""
    rng = np.random.default_rng(seed)
    H, W = dod.shape
    mm = np.isfinite(dod) & np.isfinite(lod) & (lod > 0)
    if mask is not None:
        mm = mm & mask
    st = stable & mm
    sig_cell = np.maximum(lod / 1.96, 1e-6)             # per-cell error (heteroscedastic)
    sigma = float(1.4826 * np.median(np.abs(dod[st] - np.median(dod[st]))))
    L = _corr_length(dod / sig_cell, st, res)          # correlation length of the whitened error
    Lc = max(L / res, 0.5)
    tau = _systematic_floor(dod, stable, res, sys_block_m)
    ca = res * res

    kernels = [max(s / res / np.sqrt(2), 0.0) for s in scales_m]
    # smoothed-noise factor r_s from GRF nulls with the data's correlation length
    r_s = np.zeros(len(kernels))
    grfs = [_grf((H, W), Lc / np.sqrt(2), rng) for _ in range(nsim)]
    for g in grfs:
        for k, sk in enumerate(kernels):
            r_s[k] += (gaussian_filter(g, sk) if sk > 0 else g)[mm].std() / nsim

    lab, regs = _extract(dod, sig_cell, mm, L, tau, res, kernels, r_s, t_seed, min_cells)
    # null-calibrate the SNR threshold: same procedure on whitened GRF nulls scaled
    # back by the per-cell error, so the null carries the same heteroscedasticity.
    null_max = []
    for g in grfs:
        _, nr = _extract(g * sig_cell, sig_cell, mm, L, tau, res, kernels, r_s, t_seed, min_cells)
        null_max.append(max((r[6] for r in nr), default=0.0))
    snr_thresh = float(np.percentile(null_max, z_field))

    labels = np.zeros((H, W), int); regions = []; kid = 0
    for (_, cm, a, mean, neff, sm_, snr) in regs:
        if snr <= snr_thresh:
            continue
        kid += 1; labels[cm] = kid
        core = cm & (np.abs(dod) > lod)
        regions.append(dict(id=kid, n_cells=a, area_ha=a * ca / 1e4,
                            mean_m=round(mean, 4), volume_m3=round(float(dod[cm].sum() * ca), 1),
                            sign="deposition" if mean > 0 else "erosion",
                            n_eff=round(neff, 1), snr=round(snr, 1),
                            core_frac=round(float(core.sum() / a), 3)))
    regions.sort(key=lambda r: -r["area_ha"])
    return dict(labels=labels, regions=regions, change=labels > 0,
                sigma=round(sigma, 4), corr_length_m=round(L, 1),
                tau_sys_m=round(float(tau), 4), snr_thresh=round(snr_thresh, 2))


def _grf(shape, sig_cells, rng):
    w = gaussian_filter(rng.standard_normal(shape), max(sig_cells, 1e-3))
    s = w.std()
    return w / (s if s > 0 else 1.0)
