"""Wheaton et al. (2010) spatial-coherence Bayesian thresholding: the neighbourhood
COUNT promotes a coherent sub-LoD patch and demotes an isolated high spike -- the
behaviour that distinguishes it from amplitude smoothing (which would halo)."""
import numpy as np

from lidar_diff_icp.coherence import coherence_change, spatial_coherence_probability


def test_promotes_coherent_subLoD_patch_demotes_isolated_spike():
    rng = np.random.default_rng(0)
    H = W = 60
    perr = 0.10                                     # per-cell error -> LoD ~ 0.196 m
    dod = rng.normal(0, 0.05, (H, W))               # sub-LoD noise
    dod[20:40, 20:40] += 0.15                       # coherent deposition, BELOW per-cell LoD (0.196)
    dod[50, 50] = 0.8                               # isolated high spike (blunder)
    perror = np.full((H, W), perr)

    assert 0.15 < 1.96 * perr                        # the patch really is sub-LoD per cell
    chg = coherence_change(dod, perror, conf=0.95)
    assert chg[25:35, 25:35].mean() > 0.8           # coherent sub-LoD patch -> promoted
    assert not chg[50, 50]                           # isolated spike -> demoted (no coherent neighbours)
    # posterior is signed (erosion negative, deposition positive)
    post = spatial_coherence_probability(dod, perror)
    assert post[30, 30] > 0                          # deposition patch positive
