"""Run the DECLARED graph, not a paraphrase of it.

Every other workflow test checks the graph's SHAPE -- that the order satisfies the
requirements, that no two steps claim an output, that staleness propagates. None of them
executes anything, so a step whose command names a script that has been renamed, or whose
`produces` no longer matches what the producer writes, passes every one of them. The graph
would then be a confident description of a pipeline that cannot run.

This drives one real chain: build a tiny tile with difference_dem, write its base products
exactly as run_all_sites does, then run the graph's OWN command string for the next step
and check the product it declares actually appears and the state machine agrees.

It uses a real directory under data/derived because that is what the declared commands
address -- make_slope.py takes a tile NAME, not a path. Substituting a different command
would test a string this repo does not use.
"""
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from lidar_diff_icp import workflow as W
from lidar_diff_icp.pipeline import difference_dem

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))
from test_pipeline import _make_tiles, BOUNDS          # noqa: E402

TILE = "_graph_e2e_test"
TILE_DIR = os.path.join("data", "derived", TILE)


@pytest.fixture
def built_tile(tmp_path):
    """A real tile directory holding real base products, removed afterwards."""
    if os.path.exists(TILE_DIR):
        shutil.rmtree(TILE_DIR)
    os.makedirs(TILE_DIR)
    try:
        before, after = _make_tiles(tmp_path)
        r = difference_dem(before, after, BOUNDS, res=5.0, ground_q=0.10,
                           ground="low_q", ground_source="last_return",
                           after_ground="last_return", geoid_datum=(0.0, 0.0, 0.0),
                           valley_top_m=-1e9)
        # written the way scripts/run_all_sites.py writes them -- same files, same names
        np.save(f"{TILE_DIR}/dod.npy", r["dod"])
        np.save(f"{TILE_DIR}/lod.npy", r["lod"])
        np.save(f"{TILE_DIR}/z_after.npy", r["z_after"])
        with open(f"{TILE_DIR}/corrections.json", "w") as fh:
            json.dump(r["corrections"], fh, indent=2)
        yield TILE_DIR
    finally:
        shutil.rmtree(TILE_DIR, ignore_errors=True)


def test_the_base_step_declares_what_the_pipeline_actually_writes(built_tile):
    """Step("base").produces is a claim about difference_dem's output. If the pipeline
    stopped writing one of these -- or renamed it -- the graph would still order every
    downstream step against a file that never appears."""
    base = next(s for s in W.STEPS if s.name == "base")
    for f in base.produces:
        assert os.path.exists(os.path.join(built_tile, f)), f"base does not produce {f}"
    st = W.state(built_tile)
    assert st["base"][0] == "OK", st["base"]


def test_the_declared_slope_command_runs_and_produces_what_it_says(built_tile):
    """THE POINT OF THIS FILE. Takes the command string out of the graph, substitutes the
    tile the way plan() does, runs it, and checks the declared product appears."""
    cmd = next(c for s, c, _ in W.plan(built_tile) if s.name == "slope")
    assert TILE in cmd, cmd
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    assert p.returncode == 0, f"the declared command failed:\n{cmd}\n{p.stderr[-2000:]}"

    step = next(s for s in W.STEPS if s.name == "slope")
    for f in step.produces:
        assert os.path.exists(os.path.join(built_tile, f)), \
            f"slope claims to produce {f} and did not"

    slope = np.load(os.path.join(built_tile, "slope.npy"))
    assert slope.shape == np.load(os.path.join(built_tile, "z_after.npy")).shape
    assert np.isfinite(slope).any()
    assert np.nanmin(slope) >= 0.0 and np.nanmax(slope) <= 90.0

    st = W.state(built_tile)
    assert st["slope"][0] == "OK", st["slope"]


def test_a_step_whose_inputs_are_absent_reports_missing_not_ok(built_tile):
    """The state machine must not call a step OK because nothing contradicted it."""
    st = W.state(built_tile)
    assert st["convexity"][0] == "MISSING"
    assert st["ridge_mask"][0] == "MISSING"


def test_running_a_step_clears_the_staleness_its_own_input_caused(built_tile):
    """z_after touched after slope.npy makes slope STALE; re-running the DECLARED command
    clears it. This is the loop a user actually performs, end to end."""
    cmd = next(c for s, c, _ in W.plan(built_tile) if s.name == "slope")
    assert subprocess.run(cmd, shell=True, capture_output=True).returncode == 0
    assert W.state(built_tile)["slope"][0] == "OK"

    z = os.path.join(built_tile, "z_after.npy")
    np.save(z, np.load(z) + 1.0)                       # genuinely new content
    assert W.state(built_tile)["slope"][0] == "STALE"

    assert subprocess.run(cmd, shell=True, capture_output=True).returncode == 0
    assert W.state(built_tile)["slope"][0] == "OK"
