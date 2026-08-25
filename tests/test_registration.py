"""Tests for the cross-epoch registration terms (lidar_diff_icp.registration).

Built on analytically known surfaces: a tilted plane whose gradient and slope are exact,
so the lateral term's magnitude, its tan(slope) growth, its aspect sign change, and its
vanishing on flat ground can all be checked against closed form rather than a golden file.
"""
import json
import math

import numpy as np
import pytest

from lidar_diff_icp import registration as reg

DATUM = {"method": "geoid_difference", "const_m": 0.067,
         "tilt_b_m_per_km": 0.0, "tilt_c_m_per_km": 0.0, "centroid": [0.0, 0.0],
         "horizontal_shift_m": [-0.75, -0.19]}


def plane(slope_deg, azimuth_deg, res=5.0, n=21):
    """Grid tilted by `slope_deg`, falling toward `azimuth_deg` (deg CW from +y/North)."""
    t = math.tan(math.radians(slope_deg))
    ex, ey = math.sin(math.radians(azimuth_deg)), math.cos(math.radians(azimuth_deg))
    yy, xx = np.mgrid[0:n, 0:n] * res
    return -t * (ex * xx + ey * yy)


def test_gradients_match_an_analytic_plane():
    for slope in (0.0, 10.0, 25.0):
        gx, gy, nn = reg.surface_gradients(plane(slope, 90.0), 5.0)
        assert gx[5, 5] == pytest.approx(-math.tan(math.radians(slope)), abs=1e-9)
        assert gy[5, 5] == pytest.approx(0.0, abs=1e-9)
        assert nn[5, 5] == pytest.approx(math.sqrt(math.tan(math.radians(slope))**2 + 1), abs=1e-9)


def test_gradients_fill_nan_before_differencing():
    z = plane(10.0, 90.0); z[4, 4] = np.nan
    gx, gy, nn = reg.surface_gradients(z, 5.0)
    assert np.all(np.isfinite(gx)) and np.all(np.isfinite(gy)) and np.all(np.isfinite(nn))


def test_lateral_term_vanishes_on_flat_ground():
    gx, gy, nn = reg.surface_gradients(plane(0.0, 0.0), 5.0)
    lat = reg.lateral_term(gx, gy, nn, DATUM)
    assert np.allclose(lat, 0.0, atol=1e-9)


def test_lateral_term_grows_as_tan_slope():
    """Magnitude with the shift aligned down the fall line must be |shift|*tan(slope),
    converted to slope-normal by /|n|."""
    dx, dy = DATUM["horizontal_shift_m"]
    mag = math.hypot(dx, dy)
    az = math.degrees(math.atan2(dx, dy))            # azimuth of the shift vector
    prev = 0.0
    for slope in (5.0, 10.0, 20.0, 30.0):
        gx, gy, nn = reg.surface_gradients(plane(slope, az), 5.0)
        lat = reg.lateral_term(gx[5, 5], gy[5, 5], nn[5, 5], DATUM)
        t = math.tan(math.radians(slope))
        assert abs(lat) == pytest.approx(1000 * mag * t / math.sqrt(t*t + 1), rel=1e-6)
        assert abs(lat) > prev                        # monotone in slope
        prev = abs(lat)


def test_lateral_term_changes_sign_with_aspect():
    dx, dy = DATUM["horizontal_shift_m"]
    az = math.degrees(math.atan2(dx, dy))
    gx1, gy1, nn1 = reg.surface_gradients(plane(20.0, az), 5.0)
    gx2, gy2, nn2 = reg.surface_gradients(plane(20.0, az + 180.0), 5.0)
    a = reg.lateral_term(gx1[5, 5], gy1[5, 5], nn1[5, 5], DATUM)
    b = reg.lateral_term(gx2[5, 5], gy2[5, 5], nn2[5, 5], DATUM)
    assert a == pytest.approx(-b, rel=1e-6)
    assert a != pytest.approx(0.0, abs=1.0)
    # ... and averages away over a full turn of aspects, which is why a biased aspect
    # sample (e.g. one flank of a divide) is what makes it dangerous
    vals = []
    for k in range(36):
        g1, g2, nnk = reg.surface_gradients(plane(20.0, 10.0 * k), 5.0)
        vals.append(reg.lateral_term(g1[5, 5], g2[5, 5], nnk[5, 5], DATUM))
    assert abs(np.mean(vals)) < 0.01 * max(abs(v) for v in vals)


def test_geoid_term_is_constant_on_flat_ground_and_matches_the_datum():
    gx, gy, nn = reg.surface_gradients(plane(0.0, 0.0), 5.0)
    g = reg.geoid_term(np.zeros(5), np.zeros(5), nn[5, 5], DATUM)
    assert np.allclose(g, 67.0, atol=1e-9)
    steep = reg.surface_gradients(plane(30.0, 90.0), 5.0)[2][5, 5]
    assert reg.geoid_term(0.0, 0.0, steep, DATUM) == pytest.approx(67.0 / steep, rel=1e-9)


def test_geoid_tilt_is_applied_about_the_centroid():
    d = dict(DATUM, tilt_b_m_per_km=0.001, centroid=[1000.0, 0.0])
    at_centroid = reg.geoid_term(1000.0, 0.0, 1.0, d)
    a_km_east = reg.geoid_term(2000.0, 0.0, 1.0, d)
    assert at_centroid == pytest.approx(67.0, abs=1e-9)
    assert a_km_east - at_centroid == pytest.approx(1.0, abs=1e-9)      # 0.001 m/km = 1 mm


def test_corrected_offset_sums_the_terms_and_is_reversible():
    gx, gy, nn = reg.surface_gradients(plane(20.0, 45.0), 5.0)
    d = np.array([-84.0, -100.0, 12.0])
    x = y = np.zeros(3)
    corr, g, lat = reg.corrected_offset(d, x, y, gx[5, 5], gy[5, 5], nn[5, 5], DATUM)
    assert np.allclose(corr, d + g + lat)
    assert np.allclose(corr - g - lat, d)               # undoable from the stored terms
    assert not np.allclose(lat, 0.0)


def test_read_cross_epoch_datum_prefers_geoid_file_and_fails_loudly(tmp_path):
    with pytest.raises(FileNotFoundError):
        reg.read_cross_epoch_datum(str(tmp_path))
    (tmp_path / "corrections.json").write_text(json.dumps({"cross_epoch_datum": {"const_m": 1.0}}))
    assert reg.read_cross_epoch_datum(str(tmp_path))["const_m"] == 1.0
    (tmp_path / "corrections_geoid.json").write_text(
        json.dumps({"cross_epoch_datum": {"const_m": 2.0}}))
    assert reg.read_cross_epoch_datum(str(tmp_path))["const_m"] == 2.0     # geoid file wins


def test_read_cross_epoch_datum_ignores_a_file_without_the_block(tmp_path):
    (tmp_path / "corrections_geoid.json").write_text(json.dumps({"bounds": [0, 0, 1, 1]}))
    (tmp_path / "corrections.json").write_text(json.dumps({"cross_epoch_datum": {"const_m": 3.0}}))
    assert reg.read_cross_epoch_datum(str(tmp_path))["const_m"] == 3.0


def test_lateral_term_sign_is_pinned_by_geometry():
    """The one thing the magnitude tests cannot catch: a global sign flip, which would turn
    the correction into a doubling of the error.

    Geometry: a plane FALLING toward +x, so gx < 0. Shift gen1 east (dx > 0): it now sits
    over LOWER gen2 ground, so `d = z_gen1 - plane_gen2` must INCREASE. The lateral term
    is therefore positive here, and negative for a westward shift.
    """
    gx, gy, nn = reg.surface_gradients(plane(20.0, 90.0), 5.0)   # falls toward +x
    assert gx[5, 5] < 0 and abs(gy[5, 5]) < 1e-12
    east = dict(DATUM, horizontal_shift_m=[+1.0, 0.0])
    west = dict(DATUM, horizontal_shift_m=[-1.0, 0.0])
    lat_e = reg.lateral_term(gx[5, 5], gy[5, 5], nn[5, 5], east)
    lat_w = reg.lateral_term(gx[5, 5], gy[5, 5], nn[5, 5], west)
    assert lat_e > 0, "shifting gen1 downhill must RAISE the offset"
    assert lat_w < 0
    assert lat_e == pytest.approx(1000 * abs(gx[5, 5]) / nn[5, 5], rel=1e-9)


def test_correction_removes_a_known_synthetic_misregistration():
    """End-to-end: plant a known lateral error, confirm the term cancels it exactly."""
    z = plane(18.0, 30.0)
    gx, gy, nn = reg.surface_gradients(z, 5.0)
    dx, dy = 0.75, -0.2
    # a gen1 return that is truly ON the surface but recorded at (x-dx, y-dy):
    # its measured slope-normal offset is the plane difference over that displacement
    d_meas = -1000.0 * (gx[5, 5] * dx + gy[5, 5] * dy) / nn[5, 5] * -1
    datum = {"const_m": 0.0, "tilt_b_m_per_km": 0.0, "tilt_c_m_per_km": 0.0,
             "centroid": [0.0, 0.0], "horizontal_shift_m": [dx, dy]}
    corr, g, lat = reg.corrected_offset(d_meas, 0.0, 0.0, gx[5, 5], gy[5, 5], nn[5, 5], datum)
    assert g == pytest.approx(0.0)
    assert corr == pytest.approx(0.0, abs=1e-9), "applying the known shift must zero the error"


ALIGN = {135: (0.0, 0.0, 0.0), 136: (0.32, -0.08, -0.024), 138: (1.05, 0.09, -0.044)}
CURVES = {135: (np.array([100.0, 200.0]), np.array([0.000, 0.020])),
          136: (np.array([100.0, 200.0]), np.array([-0.010, 0.010]))}


def test_swath_alignment_reference_swath_contributes_nothing():
    gx, gy, nn = reg.surface_gradients(plane(20.0, 30.0), 5.0)
    t = reg.swath_alignment_term(np.array([135, 135]), gx[5, 5], gy[5, 5], nn[5, 5], ALIGN)
    assert np.allclose(t, 0.0), "the reference flight line is the gauge and must not move"


def test_swath_alignment_matches_closed_form_per_swath():
    gx, gy, nn = reg.surface_gradients(plane(20.0, 30.0), 5.0)
    psid = np.array([135, 136, 138])
    t = reg.swath_alignment_term(psid, gx[5, 5], gy[5, 5], nn[5, 5], ALIGN)
    for i, s in enumerate(psid):
        ax, ay, az = ALIGN[s]
        assert t[i] == pytest.approx(1000 * (az - (gx[5, 5]*ax + gy[5, 5]*ay)) / nn[5, 5], rel=1e-9)
    assert t[0] == 0.0 and abs(t[2]) > abs(t[1])          # 138 is displaced furthest


def test_swath_alignment_vertical_part_survives_on_flat_ground():
    """On flat ground the lateral part vanishes but the dz part must NOT: that vertical
    swath-to-swath disagreement is exactly the internal inconsistency being removed."""
    gx, gy, nn = reg.surface_gradients(plane(0.0, 0.0), 5.0)
    t = reg.swath_alignment_term(np.array([136]), gx[5, 5], gy[5, 5], nn[5, 5], ALIGN)
    assert t[0] == pytest.approx(-24.0, abs=1e-6)


def test_swath_alignment_raises_on_an_unmapped_swath():
    gx, gy, nn = reg.surface_gradients(plane(10.0, 0.0), 5.0)
    with pytest.raises(KeyError, match="137"):
        reg.swath_alignment_term(np.array([135, 137]), gx[5, 5], gy[5, 5], nn[5, 5], ALIGN)


def test_drift_interpolates_within_each_swath_independently():
    psid = np.array([135, 135, 136, 136])
    gt = np.array([100.0, 150.0, 100.0, 150.0])
    t = reg.along_track_drift_term(psid, gt, 1.0, CURVES)
    assert t[0] == pytest.approx(0.0) and t[1] == pytest.approx(10.0)      # 0 -> 20 mm
    assert t[2] == pytest.approx(-10.0) and t[3] == pytest.approx(0.0)     # -10 -> +10 mm


def test_drift_clamps_outside_the_curve_span():
    t = reg.along_track_drift_term(np.array([135, 135]), np.array([0.0, 1e6]), 1.0, CURVES)
    assert t[0] == pytest.approx(0.0) and t[1] == pytest.approx(20.0)


def test_drift_raises_on_a_swath_without_a_curve():
    with pytest.raises(KeyError, match="999"):
        reg.along_track_drift_term(np.array([999]), np.array([100.0]), 1.0, CURVES)


def test_drift_is_slope_normalised():
    flat = reg.along_track_drift_term(np.array([135]), np.array([200.0]), 1.0, CURVES)
    steep_nn = reg.surface_gradients(plane(30.0, 0.0), 5.0)[2][5, 5]
    steep = reg.along_track_drift_term(np.array([135]), np.array([200.0]), steep_nn, CURVES)
    assert steep[0] == pytest.approx(flat[0] / steep_nn, rel=1e-9)


def test_registration_terms_sum_to_d_corr_and_stay_separable(tmp_path):
    (tmp_path / "corrections.json").write_text(json.dumps({
        "cross_epoch_datum": dict(DATUM),
        "per_swath_internal_alignment_dxdydz_m": {str(k): list(v) for k, v in ALIGN.items()},
        "along_track_drift_gpsTime_to_m": {
            str(k): {"gps_time": list(t), "drift_m": list(d)} for k, (t, d) in CURVES.items()}}))
    gx, gy, nn = reg.surface_gradients(plane(15.0, 40.0), 5.0)
    d = np.array([-84.0, -60.0])
    out = reg.registration_terms(d, np.zeros(2), np.zeros(2), np.array([120.0, 180.0]),
                                 np.array([135, 136]), gx[5, 5], gy[5, 5], nn[5, 5], str(tmp_path))
    assert set(out) == {"geoid", "lateral", "swath", "drift", "d_corr"}
    assert np.allclose(out["d_corr"], d + out["geoid"] + out["lateral"] + out["swath"] + out["drift"])
    assert np.allclose(out["swath"][0], 0.0)                      # reference swath
    assert not np.allclose(out["swath"][1], 0.0)
    assert not np.allclose(out["drift"], 0.0)
