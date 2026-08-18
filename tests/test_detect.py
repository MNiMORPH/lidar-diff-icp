"""Regression test for the two-axis change detector.

The property that matters and that per-cell-LoD-first cannot do: a broad, coherent
feature whose every cell is BELOW the per-cell LoD must still be detected (its
area-mean is significant via N_eff), while pure spatially-correlated noise is not.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, label

from lidar_diff_icp.detect import detect_change, detect_change_standard


def _field(rng, H=150, W=150, sig=0.08, Lc=4.0):
    n = gaussian_filter(rng.standard_normal((H, W)), Lc / np.sqrt(2))
    return n * (sig / n.std())


def test_detects_subLoD_broad_feature_and_rejects_noise():
    rng = np.random.default_rng(1)
    H = W = 220; res = 5.0; sig = 0.08
    dod = _field(rng, H, W, sig)
    truth = np.zeros((H, W), bool)
    dod[60:160, 60:160] += 0.12          # +0.12 m, BELOW per-cell LoD (0.157); 25 ha
    truth[60:160, 60:160] = True
    lod = np.full((H, W), 1.96 * sig)
    r = detect_change(dod, lod, ~truth, res, nsim=40)

    # per-cell-LoD-first would recover almost none of it as a coherent region:
    exceed = np.abs(dod) > lod
    lab, n = label(exceed & truth, structure=np.ones((3, 3)))
    biggest = (np.bincount(lab.ravel())[1:].max() if n else 0)
    assert biggest < 0.2 * truth.sum()               # the feature is mostly sub-LoD

    det = r["change"]
    assert det[60:120, 60:120].mean() > 0.5          # ...yet the detector finds it
    assert (det & ~truth).sum() / (~truth).sum() < 0.05   # and stays specific


def test_pure_noise_gives_few_false_positives():
    rng = np.random.default_rng(3)
    H = W = 150; res = 5.0; sig = 0.08
    dod = _field(rng, H, W, sig)
    lod = np.full((H, W), 1.96 * sig)
    r = detect_change(dod, lod, np.ones((H, W), bool), res, nsim=40)
    assert r["change"].mean() < 0.03                 # ~no confirmed change on noise


def test_standard_detector_detects_patch_reports_regions():
    """The standard (Wheaton coherence + tau_sys floor) detector finds a clear
    coherent patch, returns per-region stats, and the tau_sys floor keeps the
    stable background out."""
    rng = np.random.default_rng(0)
    H = W = 80; res = 5.0; perr = 0.08
    dod = rng.normal(0, 0.04, (H, W))
    dod[20:50, 20:50] += 0.25                        # strong coherent deposition
    lod = np.full((H, W), 1.96 * perr)
    stable = np.ones((H, W), bool); stable[20:50, 20:50] = False
    r = detect_change_standard(dod, lod, stable, res)
    assert r["change"][30:40, 30:40].mean() > 0.8    # real patch detected
    assert any(g["sign"] == "deposition" and g["area_ha"] > 2 for g in r["regions"])
    assert r["change"][:15, :15].mean() < 0.05       # stable corner stays clean
    assert {"labels", "regions", "change", "sigma", "tau_sys_m"} <= set(r)
