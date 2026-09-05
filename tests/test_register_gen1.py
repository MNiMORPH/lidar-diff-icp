"""Registering the BEFORE (gen1) cloud into its own frame.

This stage was 82 lines inline in difference_dem, so none of it could be reached without a
real tile and a CSF run. What is pinned here is the contract the later stages depend on:
it is gen1-INTERNAL (no gen2 data, no absolute level), it returns arrays over the WHOLE
cloud rather than just the ground selection, and the per-swath constants it solves are
actually applied to the points it hands back.
"""
import numpy as np
import laspy
import pytest

from lidar_diff_icp import pipeline


BOUNDS = (0.0, 0.0, 460.0, 400.0)
RES = 5.0


def _strips(tmp_path, *, n_sw=4, per=20_000, seed=0, name="before.las", mixed=False):
    """Four overlapping N-S strips over a CURVED surface, each with a known z offset.

    The curvature matters: align_swaths ties a pair through the relationship between dh
    and the local slope, so a single flat plane gives it nothing to fit and every pair is
    dropped as non-finite. This is the geometry tests/test_coreg.py uses.
    """
    rng = np.random.default_rng(seed)
    xs, ys, zs, ps, sa = [], [], [], [], []
    for k in range(n_sw):
        x0 = 100.0 * k
        width = 60.0 if k else 30.0
        x = rng.uniform(x0, x0 + 100.0 + width, per)
        y = rng.uniform(0.0, 400.0, per)
        z = (8.0 * np.sin(x / 40.0) + 6.0 * np.cos(y / 55.0) + 0.02 * x
             + (0.010 * k - 0.015) + rng.normal(0.0, 0.01, per))
        xs.append(x); ys.append(y); zs.append(z)
        ps.append(np.full(per, 10 + k))
        sa.append(np.clip((x - (x0 + 50.0 + width / 2)) / 3.0, -17.0, 17.0))
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    hdr = laspy.LasHeader(version="1.4", point_format=6)
    hdr.offsets = [0.0, 0.0, 0.0]; hdr.scales = [0.001, 0.001, 0.001]
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.point_source_id = np.concatenate(ps).astype(np.uint16)
    las.gps_time = np.linspace(0.0, 900.0, x.size)
    cls = np.full(x.size, 2, np.uint8)
    rn = np.ones(x.size, np.uint8); nr = np.ones(x.size, np.uint8)
    if mixed:
        # a cloud where the two ground selections genuinely disagree: some points are
        # class 1 (unclassified) and some are FIRST of two returns, so neither
        # "class2" nor "last_return" is everything, and neither is a subset of the other.
        cls[rng.random(x.size) < 0.30] = 1
        first_of_two = rng.random(x.size) < 0.20
        rn[first_of_two] = 1; nr[first_of_two] = 2
    las.classification = cls
    las.return_number = rn
    las.number_of_returns = nr
    las.scan_angle = (np.concatenate(sa) / 0.006).astype(np.int16)
    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / name
    las.write(str(p))
    return str(p)


def _reg(path, **kw):
    kw.setdefault("ground_source", "last_return")
    kw.setdefault("before_crs", "EPSG:26915")
    kw.setdefault("verbose", False)
    return pipeline.register_gen1(path, BOUNDS, RES, **kw)


def test_it_returns_the_whole_cloud_not_just_the_ground(tmp_path):
    """`ground` is a MASK over the returned arrays, not a filter already applied to them.
    difference_dem indexes with xc[be]; if the arrays came back pre-filtered, that index
    would silently take a subset of a subset."""
    r = _reg(_strips(tmp_path))
    n = 4 * 20_000
    for k in ("x", "y", "z", "ground", "source_id", "gps_time"):
        assert r[k].shape == (n,), (k, r[k].shape)
    assert r["ground"].dtype == bool


def test_the_solved_swath_constants_are_applied_to_the_points(tmp_path):
    """A constant that is solved but not applied would leave the swaths unaligned while
    corrections.json reported them aligned."""
    p = _strips(tmp_path)
    r = _reg(p)
    raw = laspy.read(p)
    zr = np.asarray(raw.z); ps = np.asarray(raw.point_source_id)
    for s, (dx, dy, dz) in r["swath_corr"].items():
        m = ps == s
        assert np.allclose(r["z"][m] - zr[m], dz, atol=1e-9)


def test_the_zero_line_is_the_lowest_source_id_and_carries_no_correction(tmp_path):
    """The network is solved free and the reference swath's value then subtracted, so the
    zero line's own correction is exactly zero by construction. Its identity is returned
    because it sets the absolute level the whole mosaic inherits -- measured on elbaext the
    six per-swath dz span 44.60 mm, so a different zero line moves every elevation by up to
    that much."""
    r = _reg(_strips(tmp_path))
    assert r["zero_line"] == 10
    assert r["swath_corr"][10] == (0.0, 0.0, 0.0)


def test_boresight_at_zero_is_a_no_op(tmp_path):
    """The opt-in boresight path must not perturb anything when the constant is zero --
    otherwise 'correction applied' and 'correction is zero' would not be the same state."""
    p = _strips(tmp_path)
    off = _reg(p, correct_boresight=False)
    zero = _reg(p, correct_boresight=True, boresight_roll_mm_per_deg=0.0)
    assert np.array_equal(off["z"], zero["z"])
    assert off["boresight"] is None and zero["boresight"] == 0.0


def test_boresight_is_removed_per_point_as_roll_times_scan_angle(tmp_path):
    """mm per degree of scan angle, removed BEFORE the swath network is solved, so the
    empirical alignment cannot absorb it."""
    p = _strips(tmp_path)
    b = 1.5
    r = _reg(p, correct_boresight=True, boresight_roll_mm_per_deg=b)
    assert r["boresight"] == b
    raw = laspy.read(p)
    sa = np.asarray(raw.scan_angle).astype(float) * 0.006
    ps = np.asarray(raw.point_source_id); zr = np.asarray(raw.z)
    # z = (raw - b*sa/1000) + per-swath dz
    for s, (_, _, dz) in r["swath_corr"].items():
        m = ps == s
        assert np.allclose(r["z"][m], zr[m] - b * sa[m] / 1000.0 + dz, atol=1e-9)


def test_ground_source_selects_different_points(tmp_path):
    """'last_return' takes rn == nr (singles included); 'class2' takes the before survey's
    OWN vendor ground class. They are different questions and must not be conflated."""
    p = _strips(tmp_path, mixed=True)
    lr = _reg(p, ground_source="last_return")["ground"]
    c2 = _reg(p, ground_source="class2")["ground"]
    raw = laspy.read(p)
    rn = np.asarray(raw.return_number); nr = np.asarray(raw.number_of_returns)
    assert np.array_equal(lr, rn == nr)
    assert np.array_equal(c2, np.asarray(raw.classification) == 2)
    # neither is everything, and neither contains the other -- so the test can tell them
    # apart, which an all-singles all-class-2 cloud could not
    assert 0 < lr.sum() < lr.size and 0 < c2.sum() < c2.size
    assert (lr & ~c2).any() and (c2 & ~lr).any()


def test_it_reads_no_gen2_data(tmp_path):
    """Structural: the stage is gen1-internal. If a gen2 argument ever appeared here, the
    gen1 tie could absorb a cross-epoch correction -- the exact failure the split of
    register_gen1 from apply_datum exists to prevent."""
    import inspect
    params = set(inspect.signature(pipeline.register_gen1).parameters)
    assert not {p for p in params if "after" in p or "gen2" in p}
