"""The shared terrain masks, and the refusals that keep the valley cut honest."""
import numpy as np
import pytest

from lidar_diff_icp import terrain


def _bowl(n=40, res=5.0):
    """A synthetic valley: a low flat floor and a raised rim."""
    y, x = np.mgrid[0:n, 0:n]
    r = np.hypot(x - n / 2, y - n / 2)
    return 200.0 + 40.0 * np.clip((r - 6) / 10.0, 0, 1)


def test_stable_is_low_curvature_above_the_valley_top():
    z = _bowl()
    m = terrain.terrain_masks(z, 5.0, valley_top_m=210.0, curv_max=0.005, verbose=False)
    assert not m["stable"][m["floodplain"]].any()          # never below the valley top
    assert np.all(np.abs(m["laplacian"][m["stable"]]) <= 0.005)
    assert m["valley_top_m"] == 210.0 and m["valley_top_source"] == "stated"


def test_stable_keeps_slope_because_the_lateral_fit_needs_it():
    """The old mask -- (slope<3 & upland) | (5<slope<35 & upland & lap<0) -- left
    whitewater at a MEDIAN SLOPE OF 3.0 deg, nearly flat, which is close to useless for a
    Nuth-Kaeae fit estimated from how dh varies with slope and aspect."""
    z = _bowl()
    m = terrain.terrain_masks(z, 5.0, valley_top_m=210.0, curv_max=0.005, verbose=False)
    assert m["stable"].sum() > 0
    assert m["slope_deg"][m["stable"]].max() > 5.0        # sloping ground survives


def test_a_valley_top_is_never_chosen_for_you():
    with pytest.raises(ValueError, match="will not be chosen for you"):
        terrain.resolve_valley_top(None)


def test_registry_refuses_a_site_it_does_not_know(tmp_path):
    d = tmp_path / "nowhere_site"
    d.mkdir()
    with pytest.raises(ValueError, match="no established valley top"):
        terrain.resolve_valley_top("registry", str(d))


def test_registry_returns_the_established_value():
    v, src = terrain.resolve_valley_top("registry", "data/derived/elba")
    assert v == 230.0 and src == "registry"


def test_report_records_every_cut_with_its_source():
    z = _bowl()
    m = terrain.terrain_masks(z, 5.0, valley_top_m=210.0, curv_max=0.005, verbose=False)
    keys = " ".join(m["report"])
    assert "valley top 210.0 m (stated)" in keys and "laplacian" in keys
    assert m["report"]["kept"] == int(m["stable"].sum())
