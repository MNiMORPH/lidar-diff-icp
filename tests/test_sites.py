"""One record per site, and the refusals that keep it from being invented."""
import pytest

from lidar_diff_icp.sites import SITES, Site, site


def test_a_site_carries_its_own_identity():
    """The gap this closes: difference_dem knew file paths but not WHICH TILE it was
    running, so when terrain.terrain_masks needed identity to resolve a valley top it
    arrived as a parameter named `tile_dir_for_landscape` -- a patch over a missing
    concept."""
    s = site("elba")
    assert s.tile_dir == "data/derived/elba"
    assert s.csf_cache == "data/csf_cache/elba.las"
    assert s.valley_top == "registry"


def test_an_unknown_site_refuses_and_says_what_it_knows():
    with pytest.raises(KeyError, match="no site"):
        site("nowhere")


def test_every_site_states_how_its_valley_cut_is_decided():
    """Never a default: an elevation in metres, 'registry', or 'histogram'."""
    for name, s in SITES.items():
        assert s.valley_top is not None, name
        assert s.valley_top in ("registry", "histogram") or isinstance(
            s.valley_top, (int, float)), (name, s.valley_top)


def test_the_record_is_immutable():
    """A run must not be able to edit the configuration it is running under."""
    with pytest.raises(Exception):
        site("elba").name = "other"


def test_paths_derive_from_the_name_rather_than_being_repeated():
    s = Site("demo", "a.laz", "b.laz", "histogram")
    assert s.tile_dir.endswith("demo") and s.csf_cache.endswith("demo.las")


def test_every_site_streams_so_the_calibrated_route_is_available_everywhere():
    """ground_q="calibrated" REFUSES on the in-memory path -- it grids with pandas
    groupby.quantile, which takes one quantile for all cells, so a per-cell percentile
    cannot be applied there. A site with stream=False is therefore a site the vegetation
    correction cannot run on, and a six-site run raises on it. battlecreek was that site
    until 2026-09-05.

    If a future site genuinely needs the in-memory path, this test should be changed
    deliberately and the site named as correction-incapable -- not left to surface as a
    refusal in the middle of a run."""
    from lidar_diff_icp.sites import SITES
    assert [n for n, s in SITES.items() if not s.stream] == []
