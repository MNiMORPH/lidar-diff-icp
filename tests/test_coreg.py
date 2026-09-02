"""Regression tests for Nuth & Kaeaeb co-registration.

The key failure these guard against is a sign/convention error in the shift
update: an early version added the offset instead of subtracting it, which
recovered the *negative* of the true shift and then diverged. These tests inject
a known shift into a synthetic surface and require it to be recovered (as its
inverse) to sub-centimetre accuracy. They need no downloaded data.
"""
import numpy as np
import pytest

from lidar_diff_icp import coreg


def _synthetic_surface(res=1.0, n=400):
    """A smooth surface with slopes in both directions (varied aspect)."""
    yy, xx = np.mgrid[0:n, 0:n] * res
    return 8.0 * np.sin(xx / 40.0) + 6.0 * np.cos(yy / 55.0) + 0.02 * xx


def test_recovers_known_shift():
    res = 1.0
    z_ref = _synthetic_surface(res)
    for sx, sy, sz in [(1.5, -0.8, 0.05), (-2.3, 1.1, -0.12), (3.0, -3.0, 0.2)]:
        z_src = coreg._shift_grid(z_ref, sx, sy, res) + sz
        c = coreg.nuth_kaab(z_ref, z_src, res)
        assert c.converged, f"did not converge for shift {(sx, sy, sz)}"
        # N&K returns the correction: the inverse of the injected shift
        assert abs(c.dx + sx) < 0.05, f"dx: {c.dx} vs {-sx}"
        assert abs(c.dy + sy) < 0.05, f"dy: {c.dy} vs {-sy}"
        assert abs(c.dz + sz) < 0.02, f"dz: {c.dz} vs {-sz}"


def test_translation_tilt_recovers_known():
    res = 1.0
    z_ref = _synthetic_surface(res)
    n = z_ref.shape[0]
    gy, gx = np.mgrid[0:n, 0:n]
    Xc = (gx + 0.5) * res; Xc = Xc - Xc.mean()
    Yc = (gy + 0.5) * res; Yc = Yc - Yc.mean()
    sx, sy, a0, a1, a2 = 1.5, -0.8, 0.05, 0.001, -0.0005   # tilt 1.0 / -0.5 mm/m
    z_src = coreg._shift_grid(z_ref, sx, sy, res) + (a0 + a1 * Xc + a2 * Yc)
    r = coreg.tie_translation_tilt(z_ref, z_src, res, 0, 0)
    assert r["converged"]
    assert abs(r["dx"] + sx) < 0.05 and abs(r["dy"] + sy) < 0.05
    assert abs(r["c1"] + a1) < 2e-4 and abs(r["c2"] + a2) < 2e-4


def test_polynomial_tie_recovers_quadratic_warp():
    res = 1.0
    z_ref = _synthetic_surface(res)
    n = z_ref.shape[0]
    gy, gx = np.mgrid[0:n, 0:n]
    X = (gx + 0.5) * res; Y = (gy + 0.5) * res
    xm = 0.5 * (X.max() + X.min()); xhr = 0.5 * (X.max() - X.min())
    ym = 0.5 * (Y.max() + Y.min()); yhr = 0.5 * (Y.max() - Y.min())
    Xn = (X - xm) / xhr; Yn = (Y - ym) / yhr
    sxf = 0.8 + 0.4 * Xn - 0.2 * Yn + 0.15 * Xn * Xn
    syf = -0.5 + 0.3 * Yn
    szf = 0.05 + 0.02 * Xn
    z_src = coreg._warp_grid(z_ref, sxf, syf, res) + szf
    r = coreg.tie_polynomial(z_ref, z_src, res, 0, 0, order=2)
    assert r["converged"]
    assert r["nmad_after"] < 0.02
    assert np.sqrt(np.nanmean((r["dx_field"] + sxf) ** 2)) < 0.03


def test_zero_shift_is_zero():
    res = 1.0
    z = _synthetic_surface(res)
    c = coreg.nuth_kaab(z, z.copy(), res)
    assert np.hypot(c.dx, c.dy) < 0.02
    assert abs(c.dz) < 0.01
    assert c.nmad_after <= c.nmad_before + 1e-9


def test_correction_surface_recovers_warp_and_ignores_real_change():
    """DeLong-style correction surface: it must recover a smooth vertical warp
    from STABLE (flat, small-dz) ground, and must NOT absorb a large localized
    change that the masks exclude -- the surface should interpolate smoothly
    through it, leaving the change in the corrected difference."""
    res = 5.0; n = 80
    yy, xx = np.mgrid[0:n, 0:n] * res
    z_src = np.zeros((n, n))                       # flat ground -> all "stable"
    warp = 0.15 * np.sin(2 * np.pi * xx / 300.0)   # smooth spurious warp, <0.7 m
    z_ref = z_src + warp
    z_ref[40:45, 40:45] += 2.0                     # localized real change
    cs = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0,
                                  slope_thresh_deg=3.0, dz_thresh=0.7, radius=400.0)
    C = cs["C"]
    assert not cs["stable"][42, 42]                # change patch excluded from stable
    far = np.ones((n, n), bool); far[35:50, 35:50] = False
    assert np.sqrt(np.nanmean((C[far] - warp[far]) ** 2)) < 0.03  # recovers warp
    corrected = cs["dz"] - C
    assert corrected[42, 42] > 1.5                 # does not absorb the +2 m change


def test_correction_surface_exclude_drops_sources():
    """A caller-supplied exclude mask (e.g. a TPI floodplain buffer) must remove
    those cells from the stable IDW sources so the correction is not fit to
    them -- otherwise real low-slope change there is absorbed."""
    res = 5.0; n = 80
    z_src = np.zeros((n, n)); z_ref = z_src + 0.10   # flat, all otherwise stable
    ex = np.zeros((n, n), bool); ex[:, :40] = True
    a = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0)
    b = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0, exclude=ex)
    assert b["n_stable"] < a["n_stable"]
    assert not b["stable"][10, 10]                   # excluded region dropped
    assert b["stable"][10, 60]                       # kept region retained


def test_along_track_drift_recovers_known():
    """Inject a smooth per-swath drift as a function of gps_time and recover it.
    The residual on stable ground is drift + noise; the spline fit must track the
    drift, not the noise."""
    rng = np.random.default_rng(0); n = 20000
    gps = np.sort(rng.uniform(0.0, 100.0, n))
    swath = np.ones(n, int)
    true = 0.04 * np.sin(np.pi * gps / 100.0)          # smooth half-sine
    change = true + rng.normal(0, 0.02, n)             # stable residual = drift + noise
    stable = np.ones(n, bool)
    drift, curves = coreg.fit_along_track_drift(gps, change, stable, swath)
    assert 1 in curves
    assert np.corrcoef(drift, true)[0, 1] > 0.9
    assert np.sqrt(np.mean((drift - true) ** 2)) < 0.02


def test_tie_falls_back_on_gentle_terrain():
    """On gentle terrain with no real warp the Nuth & Kaeaeb horizontal shift is
    unconstrained; the tie must fall back to a rigid vertical offset, never a
    runaway (dx ~ 1000 m) that makes the fit worse than the input."""
    rng = np.random.default_rng(0); res = 5.0; n = 60
    yy, xx = np.mgrid[0:n, 0:n] * res; cx = cy = n * res / 2
    dome = 105.0 - 1.5e-4 * ((xx - cx) ** 2 + (yy - cy) ** 2)   # < 3 deg everywhere
    z_ref = dome + rng.normal(0, 0.02, (n, n))
    z_src = dome + rng.normal(0, 0.02, (n, n))                  # no warp, just noise
    t = coreg.tie_polynomial(z_ref, z_src, res, 0.0, 0.0, order=2)
    # safe by property, however reached (rigid fallback or early break to
    # identity): no runaway horizontal shift, and never worse than the input.
    assert np.ptp(t["dx_field"]) < 5.0
    assert t["nmad_after"] <= t["nmad_before"] + 1e-9


def test_slope_aspect_planar():
    """A plane rising toward the east has slope arctan(gradient) and aspect 90
    deg (steepest ascent toward east, clockwise from north)."""
    res = 1.0; n = 20
    _, xx = np.mgrid[0:n, 0:n] * res
    z = 0.1 * xx                                   # rises toward +x (east)
    slope, aspect = coreg.slope_aspect(z, res)
    inner = (slice(2, -2), slice(2, -2))           # avoid edge gradients
    assert abs(np.median(slope[inner]) - np.arctan(0.1)) < 1e-6
    assert abs(np.degrees(np.median(aspect[inner])) - 90.0) < 1.0


def test_align_swaths_recovers_vertical_offset():
    """Two overlapping swaths, swath 2 biased +0.3 m in z; the free-network
    alignment must pin swath 1 and recover ~ -0.3 m for swath 2."""
    from lidar_diff_icp import io
    rng = np.random.default_rng(0); n = 15000
    def ground(x, y):
        return 100.0 + 0.05 * x + 0.03 * y         # ~3.4 deg slope
    xa = rng.uniform(0, 120, n); ya = rng.uniform(0, 150, n)     # swath 1
    xb = rng.uniform(80, 200, n); yb = rng.uniform(0, 150, n)    # swath 2 (overlap 80-120)
    za = ground(xa, ya) + rng.normal(0, 0.02, n)
    zb = ground(xb, yb) + 0.3 + rng.normal(0, 0.02, n)           # +0.3 m bias
    x = np.concatenate([xa, xb]); y = np.concatenate([ya, yb]); z = np.concatenate([za, zb])
    ps = np.concatenate([np.ones(n), np.full(n, 2)]).astype(int)
    pc = io.PointCloud(x, y, z, ps, np.zeros_like(z), np.zeros_like(z),
                       np.zeros_like(ps), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)
    assert abs(corr[1][2]) < 0.05                  # reference pinned
    assert abs(corr[2][2] + 0.3) < 0.08            # +0.3 m bias recovered as -0.3


def test_align_swaths_ignores_nan_edge(monkeypatch):
    """A non-adjacent swath pair sharing gridded extent but no actual cells
    returns n=0 with NaN shifts (real multi-swath tiles do this). A single NaN
    observation must NOT poison the least-squares network -- the connected chain
    1-2-3 still determines every swath. Regression for the Carlton all-NaN DoD."""
    from lidar_diff_icp import io
    def mk(dz, n, n_dz=None):                       # a pairwise Coreg observation
        # n_dz defaults to n here: these fakes stand for ordinary cosine-fit observations,
        # where the cells that determined dz ARE the cells of the fit.
        return coreg.Coreg(0.0, 0.0, dz, 0.0, 0.0, 0.0, n, 0.1, 0.05, 3, True,
                           n if n_dz is None else n_dz)
    fake = {(1, 2): mk(-0.02, 1000), (2, 3): mk(-0.03, 1000),
            (1, 3): mk(np.nan, 0)}                  # the poisoning empty-overlap edge
    monkeypatch.setattr(coreg, "coregister_swaths",
                        lambda pc, a, b, res, exclude, tie="overlap_median": fake[(a, b)])
    ps = np.array([1, 1, 2, 2, 3, 3])
    pc = io.PointCloud(np.zeros(6), np.zeros(6), np.zeros(6), ps,
                       np.zeros(6), np.zeros(6), np.zeros(6, int), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)
    assert all(np.isfinite(v).all() for v in corr.values())   # not NaN-poisoned
    assert (1, 3) not in {(e[0], e[1]) for e in edges}        # empty edge dropped
    assert abs(corr[1][2]) < 1e-9                             # reference pinned
    assert abs(corr[2][2] + 0.02) < 1e-6                      # chain determines 2
    assert abs(corr[3][2] + 0.05) < 1e-6                      # ... and 3


if __name__ == "__main__":
    test_recovers_known_shift()
    test_zero_shift_is_zero()
    test_correction_surface_recovers_warp_and_ignores_real_change()
    test_correction_surface_exclude_drops_sources()
    test_along_track_drift_recovers_known()
    test_tie_falls_back_on_gentle_terrain()
    print("coreg regression tests PASS")


def _cut_edge_swaths(n_sw=4, seed=0):
    """Four overlapping N-S strips over a curved surface, each with a known z offset.

    The outermost strip is deliberately cut to half width -- the geometry of the real
    Elba tiles, where the pinned reference swath is the one the tile clips and is
    therefore sampled far off-nadir.
    """
    from lidar_diff_icp.io import PointCloud
    rng = np.random.default_rng(seed)
    xs, ys, zs, ps = [], [], [], []
    for k in range(n_sw):
        x0 = 100.0 * k
        width = 60.0 if k else 30.0
        x = rng.uniform(x0, x0 + 100.0 + width, 40000)
        y = rng.uniform(0.0, 400.0, x.size)
        z = 8.0 * np.sin(x / 40.0) + 6.0 * np.cos(y / 55.0) + 0.02 * x + (0.010 * k - 0.015)
        xs.append(x); ys.append(y); zs.append(z); ps.append(np.full(x.size, k))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    return PointCloud(x=x, y=y, z=z, point_source_id=np.concatenate(ps).astype(np.int32),
                      classification=np.full(x.size, 2, np.uint8),
                      gps_time=np.zeros(x.size), scan_angle=np.zeros(x.size),
                      crs="EPSG:26915")


def test_align_swaths_reference_is_a_gauge_not_an_observation():
    """Choosing a different reference swath shifts every constant by ONE number and
    leaves every swath-to-swath difference untouched.

    This is what makes "the reference swath is cut off-nadir, so every constant is
    measured against a biased reference" a statement about the absolute level only:
    ``align_swaths`` solves the free network first and applies the reference afterwards
    (``c -= c[idx[ref]]``). Re-solving against a symmetrically-sampled interior swath
    therefore cannot move the relative solution, and cannot change what a ground-control
    tie transports. Verified to bite: restricting the network to the edges incident on
    the reference -- i.e. measuring every swath directly AGAINST it, which is what the
    worry describes -- makes this fail. Note what does NOT break it: adding the reference
    to the design as a soft constraint still passes, because the free network has a
    rank-1 null space and every gauge on it gives the same relative solution. That
    robustness is the result, not a weakness of the test.
    """
    pc = _cut_edge_swaths()
    free, _, _ = coreg.align_swaths(pc, res=2.0)
    sw = sorted(free)
    levels = []
    for ref in sw:
        pinned, _, _ = coreg.align_swaths(pc, res=2.0, ref=ref)
        assert pinned[ref] == (0.0, 0.0, 0.0) or abs(pinned[ref][2]) < 1e-12
        for a in sw:
            for b in sw:
                for k in range(3):
                    assert (pinned[b][k] - pinned[a][k]) == pytest.approx(
                        free[b][k] - free[a][k], abs=1e-12), (
                        f"ref={ref} changed the {a}->{b} difference on axis {k}")
        levels.append(pinned[sw[0]][2] - free[sw[0]][2])
    # and the gauge really is doing something: the absolute level moves between choices
    assert max(levels) - min(levels) > 1e-3


def _across_track_pair(x_lo, x_hi, *, k_true=-0.030, c_pair=0.130, h=2500.0,
                       spacing=900.0, seed=1):
    """Two overlapping N-S flight lines whose between-line difference has a known
    across-track ramp, sampled over the easting window ``[x_lo, x_hi]``.

    Built the way the delivered data are: the lines are flown there-and-back, so the
    body-fixed scan angle runs the opposite way on line B, and a ground point at easting
    ``x`` is seen at ``tan A = (x - x_A)/h`` and ``tan B = (x_B - x)/h``. Their SUM is
    then the fixed ``spacing/h`` and only the DIFFERENCE varies across the sidelap, which
    is the geometry measured at Elba (``analysis/SWATH_ACROSS_TRACK_TEST.md`` section 0).

    Line B carries a true constant offset ``k_true`` plus a per-line across-track error,
    split so that the pair coefficient is ``c_pair`` metres per unit tangent. The tie a
    tile should recover is ``k_true``, at ``tan A = tan B`` -- the middle of the sidelap.
    """
    from lidar_diff_icp.io import PointCloud
    rng = np.random.default_rng(seed)
    x_a, x_b = 0.0, spacing                       # the two nadir tracks
    xs, ys, zs, ps, sa = [], [], [], [], []
    for line, (x_track, sign) in enumerate(((x_a, +1.0), (x_b, -1.0))):
        x = rng.uniform(x_lo, x_hi, 200000)
        y = rng.uniform(0.0, 400.0, x.size)
        tan = sign * (x - x_track) / h            # +1: (x-x_A)/h   -1: (x_B-x)/h
        # the same curved ground under both lines (varied aspect, so Nuth & Kaeaeb runs
        # normally and the horizontal solution is identical for the two tie modes).
        # Each line then carries its OWN across-track error c*tan(own scan angle); with
        # both lines given the same c the pair coefficient (c_A + c_B)/2 is c_pair.
        z = (8.0 * np.sin(x / 40.0) + 6.0 * np.cos(y / 55.0) + 0.02 * x
             + (k_true if line == 1 else 0.0)
             + c_pair * tan)
        xs.append(x); ys.append(y); zs.append(z)
        ps.append(np.full(x.size, line + 1))
        sa.append(np.degrees(np.arctan(tan)))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    return PointCloud(x=x, y=y, z=z, point_source_id=np.concatenate(ps).astype(np.int32),
                      classification=np.full(x.size, 2, np.uint8),
                      gps_time=np.zeros(x.size), scan_angle=np.concatenate(sa),
                      crs="EPSG:26915")


def test_swath_tie_intercept_is_extent_invariant():
    """The pairwise tie must not depend on which part of the sidelap the tile covers.

    ``tie="overlap_median"`` -- the shipped estimator -- averages the between-line
    difference over whatever the extent samples, so where that difference has an
    across-track ramp the tie moves with the extent. That is the named mechanism behind
    the elba/elbaext tiles disagreeing by 8.0/9.8/17.4 mm about the same flight lines.
    ``tie="intercept"`` estimates the tie at across-track position zero instead, which is
    a fixed geometric place, and must return the same number on both extents and the
    right one.

    Bites: with ``tie="intercept"`` replaced by the default the last two assertions fail
    by ~29 mm, the size of the ramp across the two windows.
    """
    k_true, c_pair, h, spacing = -0.030, 0.130, 2500.0, 900.0
    west = _across_track_pair(150.0, 500.0, k_true=k_true, c_pair=c_pair, h=h, spacing=spacing)
    wide = _across_track_pair(150.0, 900.0, k_true=k_true, c_pair=c_pair, h=h, spacing=spacing)

    med_w = coreg.coregister_swaths(west, 1, 2, res=2.0).dz
    med_d = coreg.coregister_swaths(wide, 1, 2, res=2.0).dz
    int_w = coreg.coregister_swaths(west, 1, 2, res=2.0, tie="intercept").dz
    int_d = coreg.coregister_swaths(wide, 1, 2, res=2.0, tie="intercept").dz

    # the shipped tie moves with the extent, by c * (mean dtan_1 - mean dtan_2)
    dtan = lambda lo, hi: ((lo + hi) - (0.0 + spacing)) / h   # mean of tanA - tanB
    predicted = c_pair * (dtan(150.0, 500.0) - dtan(150.0, 900.0))
    assert abs(med_w - med_d) > 0.020, "the extent-dependence this test is about is absent"
    assert med_w - med_d == pytest.approx(predicted, abs=0.004)

    # the intercept tie does not, and it recovers the true constant
    assert int_w == pytest.approx(int_d, abs=0.002)
    assert int_w == pytest.approx(-k_true, abs=0.003)
    assert int_d == pytest.approx(-k_true, abs=0.003)


def test_n_dz_travels_with_every_observation(monkeypatch):
    """The population behind dz is now recorded, though nothing filters on it yet.

    Four candidate rules were tested against the real sites and each misclassifies a real
    case -- see analysis/SWATH_TIE_DEGENERACY.md and the comment in align_swaths. This
    pins the plumbing: n_dz is required on Coreg, so no observation can be constructed
    without answering "how many cells determined this?", and it reaches the caller.
    """
    from lidar_diff_icp import io

    def mk(dz, n, n_dz, converged):
        return coreg.Coreg(0.0, 0.0, dz, 0.0, 0.0, 0.0, n, 0.1, 0.09, 20, converged, n_dz)

    fake = {(1, 2): mk(-0.02, 1000, 1000, True),
            (2, 3): mk(-0.0143, 22293, 22293, False),   # SLOW: real, must survive
            (1, 3): mk(-3.4640, 0, 36, False)}          # the extrapolated sliver tie
    monkeypatch.setattr(coreg, "coregister_swaths",
                        lambda pc, a, b, res, exclude, tie="overlap_median": fake[(a, b)])
    ps = np.array([1, 1, 2, 2, 3, 3])
    pc = io.PointCloud(np.zeros(6), np.zeros(6), np.zeros(6), ps,
                       np.zeros(6), np.zeros(6), np.zeros(6, int), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)

    kept = {(e[0], e[1]) for e in edges}
    assert (2, 3) in kept, "a slow but well-populated fit must survive"
    # NOT filtered today: the sliver tie is still admitted, at weight sqrt(n)=0
    assert (1, 3) in kept
    assert next(e for e in edges if (e[0], e[1]) == (1, 3))[5] == 0.0
    # the chain still determines swath 3 from the real edge
    assert abs((corr[3][2] - corr[2][2]) - (-0.0143)) < 1e-9


def test_a_coreg_cannot_be_built_without_stating_its_population():
    import pytest as _pt
    with _pt.raises(TypeError):
        coreg.Coreg(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10, 0.1, 0.1, 3, True)   # no n_dz


def test_the_rigid_fallback_is_kept_but_its_protection_is_nominal(monkeypatch):
    """The case the older comment protects is kept -- and never constrained anything.

    A rigid-fallback edge can have n=0 (the cosine loop broke before any fit, the guard
    then fired on a non-finite nmad). n_dz > 0, so it is rightly NOT dropped. But
    align_swaths weights the network by sqrt(n), deliberately -- coregister_swaths holds n
    at the horizontal fit's count so the two tie modes differ in the vertical estimator
    ALONE, not in the weighting. So such an edge enters at zero weight and constrains
    nothing, exactly as before this change.

    Recorded rather than fixed: making the weight n_dz would change the network solution at
    every site and break that stated design property. Whether a fallback tie should carry
    weight is a separate question.
    """
    from lidar_diff_icp import io

    def mk(dz, n, n_dz):
        return coreg.Coreg(0.0, 0.0, dz, np.nan, np.nan, 0.01, n, 0.1, 0.1, 1, True, n_dz)

    fake = {(1, 2): mk(-0.02, 0, 4200)}     # rigid fallback: n=0, but 4200 cells behind dz
    monkeypatch.setattr(coreg, "coregister_swaths",
                        lambda pc, a, b, res, exclude, tie="overlap_median": fake[(a, b)])
    ps = np.array([1, 1, 2, 2])
    pc = io.PointCloud(np.zeros(4), np.zeros(4), np.zeros(4), ps,
                       np.zeros(4), np.zeros(4), np.zeros(4, int), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)
    assert (1, 2) in {(e[0], e[1]) for e in edges}, "n_dz > 0, so it must not be dropped"
    e = next(x for x in edges if (x[0], x[1]) == (1, 2))
    assert e[5] == 0.0, "its weight is sqrt(n) = 0 -- it enters, but constrains nothing"
    assert corr[2][2] == 0.0, "so swath 2 is left where the minimum-norm solution puts it"
