"""The last three stages: the optional correction, the DoD, and the LoD.

All three were inline in difference_dem. The properties pinned here are the ones a silent
change would turn into wrong elevations rather than a crash: that the correction is OFF
unless a curve is named, that `robust_stable` changes what is REPORTED and never what is
measured, and that the LoD says which of its three routes produced it.
"""
import numpy as np
import pytest

from lidar_diff_icp import pipeline, groundest, terrain


RES = 5.0
NX, NY = 120, 100
X0, Y0 = 577492.8, 4882737.6
GRID = (X0, Y0, RES, NX, NY)


def _scene(n=300_000, seed=5, gen2_above_m=0.030, change=None):
    """A gen1 cloud and an INDEPENDENTLY sampled gen2 surface over the same terrain.

    The two epochs must carry their own noise. An earlier version of this helper built
    both from one cloud, which made the DoD an exact constant -- NMAD 0 -- and xdem's
    heteroscedastic model then returned an all-NaN LoD while still reporting itself as the
    xdem method. Every LoD assertion written against that scene passed vacuously, including
    one that was supposed to catch a reintroduced nan_to_num.
    """
    rng = np.random.default_rng(seed)
    x = rng.uniform(X0, X0 + NX * RES, n)
    y = rng.uniform(Y0, Y0 + NY * RES, n)
    base = 300.0 + 7.0 * np.sin((x - X0) / 60.0) + 5.0 * np.cos((y - Y0) / 50.0)
    z = base + rng.normal(0.0, 0.03, n)                 # gen1
    ground_of = lambda a, b, c: groundest.cellstat(a, b, c, "ground", GRID, 0.50)
    z2 = base + rng.normal(0.0, 0.02, n)                # gen2, its own returns
    Z21 = ground_of(x, y, z2)
    Zref = Z21 + gen2_above_m
    if change is not None:
        Zref = Zref + change
    tm = terrain.terrain_masks(Z21, RES, valley_top_m=float(np.nanmin(Z21)) - 1.0,
                               verbose=False)
    return dict(x=x, y=y, z=z, ground=np.ones(n, bool), Z21=Z21, Zref=Zref,
                ground_of=ground_of, tm=tm,
                spread=groundest.cellstat(x, y, z2, "spread", GRID, 0.50),
                count=groundest.cellstat(x, y, z2, "count", GRID, 0.50),
                rough=pipeline.cell_plane_roughness(x, y, z2, X0, Y0, RES, NX, NY))


def _with_real_change(s, frac=0.05, amount=-0.60, seed=1):
    """Put unmistakable erosion on a patch, so the reporting sigma-clip has something to
    remove. Without it the clip is a no-op and every off-vs-on comparison is vacuous."""
    rng = np.random.default_rng(seed)
    patch = rng.random(s["Zref"].shape) < frac
    t = dict(s); t["Zref"] = np.where(patch, s["Zref"] + amount, s["Zref"])
    t["change_patch"] = patch
    return t


def _difference(s, **kw):
    return pipeline.difference(s["Zref"], s["x"], s["y"], s["z"], s["ground"],
                               s["ground_of"], s["tm"]["stable"], GRID, 0.50, **kw)


# --- correct_reference -------------------------------------------------------------------

def test_no_curve_means_no_correction_and_no_grids():
    """The correction is opt-in BY NAME. With curve=None the reference surface comes back
    untouched -- not approximately, the same object's values -- and there are no correction
    grids to mistake for a correction that happened."""
    s = _scene(n=20_000)
    Zc, grids = pipeline.correct_reference(s["Zref"], s["Z21"], None, None, GRID,
                                           verbose=False)
    assert grids is None
    assert np.array_equal(Zc, s["Zref"], equal_nan=True)


# --- difference --------------------------------------------------------------------------

def test_the_dod_is_gen2_minus_gen1():
    """Sign convention, and it is not a detail: a flipped DoD turns deposition into erosion
    everywhere. gen2 is placed a known distance above gen1 here."""
    s = _scene(gen2_above_m=0.030)
    d = _difference(s)
    assert np.nanmedian(d["dod"]) == pytest.approx(0.030, abs=0.002)
    assert np.allclose(d["dod"], s["Zref"] - d["gen1_ground"], equal_nan=True)


def test_robust_stable_changes_what_is_reported_and_not_what_is_measured():
    """THE GUARD FOR THIS STAGE. The sigma-clip refines the REPORTING mask so the quoted
    error is stable-ground error rather than real change bleeding into it. The DoD surface,
    the gen1 ground and the per-cell statistics must be bit-identical either way -- if the
    clip ever reached the measurement, it would be fitting the answer to the mask."""
    s = _with_real_change(_scene(n=200_000))
    off = _difference(s, robust_stable=False)
    on = _difference(s, robust_stable=True)
    # the clip must actually bite, or the comparison below proves nothing
    assert on["stable"].sum() < off["stable"].sum()
    for k in ("dod", "gen1_ground", "spread", "count", "rough"):
        assert np.array_equal(off[k], on[k], equal_nan=True), k
    # only the reporting population and what is derived FROM it may move
    assert on["stable"].sum() <= off["stable"].sum()
    assert off["stable_clip_fraction"] == 0.0
    assert on["stable_geom_n"] == off["stable_geom_n"]


def test_the_clip_fraction_is_what_the_clip_actually_removed():
    s = _with_real_change(_scene(n=200_000))
    on = _difference(s, robust_stable=True)
    removed = 1.0 - on["stable"].sum() / on["stable_geom_n"]
    assert removed > 0.0
    assert on["stable_clip_fraction"] == pytest.approx(removed)


def test_sigma_is_the_nmad_of_the_dod_on_the_reported_stable_mask():
    """Named statistic, not a bare number: sigma is the NMAD of the DoD over the cells the
    product reports as stable."""
    s = _scene(n=200_000)
    d = _difference(s)
    v = d["dod"][d["stable"]]
    assert d["sigma"] == pytest.approx(1.4826 * np.median(np.abs(v - np.median(v))))


# --- estimate_lod ------------------------------------------------------------------------

def _lod(s, d, **kw):
    kw.setdefault("rough_gen1", d["rough"]); kw.setdefault("count_gen1", d["count"])
    kw.setdefault("rough_gen2", s["rough"]); kw.setdefault("count_gen2", s["count"])
    kw.setdefault("spread_gen1", d["spread"]); kw.setdefault("spread_gen2", s["spread"])
    return pipeline.estimate_lod(d["dod"], s["tm"]["slope_deg"], s["tm"]["abs_curv"],
                                 d["stable"], verbose=False, **kw)


def test_the_lod_reports_which_route_produced_it():
    """The three routes are NOT interchangeable and a product must say which it got --
    lod_method travels into corrections.json for exactly that reason."""
    s = _scene(n=200_000)
    lod, method = _lod(s, _difference(s))
    assert lod is not None and lod.shape == (NY, NX)
    assert "heteroscedastic" in method and "standard-error" in method
    # an all-NaN LoD is NOT None, so the fallback chain does not catch it -- assert the
    # model actually produced numbers, or every LoD test below is vacuous
    assert np.isfinite(lod).sum() > 0.9 * lod.size


def test_an_unmeasurable_standard_error_gives_no_lod_rather_than_an_optimistic_one():
    """An UNMEASURABLE term is not a zero one. nan_to_num used to set a missing epoch's
    variance to zero -- the most optimistic value available -- and it landed precisely on
    the cells that deserve the widest limit, because cell_plane_roughness returns NaN where
    a cell holds too few ground points. Plain arithmetic propagates the NaN instead."""
    s = _scene(n=200_000)
    d = _difference(s)
    r = d["rough"].copy()
    hole = np.zeros(r.shape, bool); hole[10:14, 10:14] = True
    r[hole] = np.nan
    sel = hole & np.isfinite(d["dod"])
    assert sel.sum() == 16                       # the assertion below is not over nothing
    base_lod, _ = _lod(s, d)
    assert np.isfinite(base_lod[sel]).all()      # ... and these cells DO get an LoD normally
    lod, _ = _lod(s, d, rough_gen1=r)
    assert np.all(np.isnan(lod[sel]))
