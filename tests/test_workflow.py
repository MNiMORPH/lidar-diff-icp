"""The derived-product graph: well-formed, correctly ordered, and honest about staleness.

These exist because both failures they cover happened for real. The convexity producer
rewrites ridgecrest_pixels.npz from scratch and the curvature producer augments it in
place, so running them in the wrong order silently drops three columns. And a q2 fit read
from a beam_offset_table older than corrections.json carries superseded registration --
which no amount of reading the file tells you, because the file looks fine.
"""
import os
import time

import pytest

from lidar_diff_icp import workflow as W


def test_every_requirement_has_a_producer_or_is_a_base_input():
    """The graph refuses to order if a step needs something nothing makes."""
    made = {f for s in W.STEPS for f in s.produces}
    for s in W.STEPS:
        for r in s.requires:
            assert r in made or r in W.BASE_INPUTS, (
                f"{s.name} requires {r}, which no step produces and which is not a base "
                f"input")


def test_the_order_satisfies_every_dependency():
    have = set(W.BASE_INPUTS)
    for s in W.order():
        unmet = [r for r in s.requires if r not in have]
        assert not unmet, f"{s.name} runs before {unmet} exist"
        have.update(s.produces)
    assert len(W.order()) == len(W.STEPS)


def test_curvature_runs_after_convexity():
    """The in-place augmentation. Wrong order = three columns silently gone."""
    names = [s.name for s in W.order()]
    assert names.index("convexity") < names.index("curvature")


def test_no_two_steps_claim_the_same_output():
    seen = {}
    for s in W.STEPS:
        for f in s.produces:
            assert f not in seen, f"{f} claimed by both {seen[f]} and {s.name}"
            seen[f] = s.name


def test_a_missing_producer_is_rejected_not_silently_ordered():
    bad = W.STEPS + (W.Step("invented", produces=("x.npy",),
                            requires=("nothing_makes_this.npy",), command="true"),)
    with pytest.raises(ValueError, match="neither base inputs nor produced"):
        W.order(bad)


def test_a_cycle_is_rejected():
    cyc = (W.Step("a", produces=("a.npy",), requires=("b.npy",), command="true"),
           W.Step("b", produces=("b.npy",), requires=("a.npy",), command="true"))
    with pytest.raises(ValueError, match="cannot order the graph"):
        W.order(cyc)


def _touch(d, name, when=None):
    p = os.path.join(d, name)
    open(p, "wb").close()
    if when is not None:
        os.utime(p, (when, when))
    return p


def test_state_reports_missing_stale_and_ok(tmp_path):
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="true")
    steps = (step,)

    _touch(d, "in.npy", when=1000)
    assert W.state(d, steps)["t"][0] == "MISSING"

    _touch(d, "out.npy", when=2000)
    assert W.state(d, steps)["t"][0] == "OK"

    os.utime(os.path.join(d, "in.npy"), (3000, 3000))
    kind, detail = W.state(d, steps)["t"]
    assert kind == "STALE" and "in.npy" in detail


def test_stale_is_detected_across_the_real_chain(tmp_path):
    """corrections.json newer than beam_offset_table.parquet -- the case that bit us."""
    d = str(tmp_path)
    now = time.time()
    for f in W.BASE_INPUTS:
        _touch(d, f, when=now - 100)
    for s in W.order():
        for f in s.produces:
            _touch(d, f, when=now - 50)
    assert all(k == "OK" for k, _ in W.state(d).values())

    os.utime(os.path.join(d, "corrections.json"), (now, now))
    st = W.state(d)
    assert st["beam_table"][0] == "STALE"
    assert "corrections.json" in st["beam_table"][1]


def test_plan_substitutes_the_tile_and_flags_unsupplied_clouds():
    cmds = W.plan("data/derived/somewhere")
    assert all("{" not in c for _, c, _ in cmds), "every placeholder is substituted"
    assert any("data/derived/somewhere" in c for _, c, _ in cmds)

    # a step that reads a cloud says so, both in the command and in its missing list
    needs_gen2 = [(s, c, m) for s, c, m in cmds if "gen2" in s.needs]
    assert needs_gen2, "some steps read the gen2 cloud"
    assert all("<--gen2 NOT GIVEN>" in c and "gen2" in m for _, c, m in needs_gen2)

    supplied = W.plan("data/derived/somewhere", gen2="X.laz")
    for s, c, m in supplied:
        if "gen2" in s.needs:
            assert "X.laz" in c and "gen2" not in m
            # a step needing BOTH clouds still flags the one still missing
            assert ("<--gen1 NOT GIVEN>" in c) == ("gen1" in s.needs)


def test_an_in_place_augmentation_does_not_make_its_own_step_look_stale():
    """curvature adds columns to ridgecrest_pixels.npz, which convexity owns.

    Declaring that file as a REQUIREMENT made the step report STALE the moment it ran --
    its own write made the input newer than its outputs. It is `mutates`, not `requires`.
    """
    curv = next(s for s in W.STEPS if s.name == "curvature")
    assert "ridgecrest_pixels.npz" in curv.mutates
    assert "ridgecrest_pixels.npz" not in curv.requires
    assert "ridgecrest_pixels.npz" not in curv.produces


def test_every_mutated_file_is_produced_by_an_earlier_step():
    assert W.mutation_order_ok() == []


def test_a_mutation_of_a_later_step_is_caught():
    steps = (W.Step("first", produces=("a.npy",), requires=(), command="true",
                    mutates=("b.npy",)),
             W.Step("second", produces=("b.npy",), requires=("a.npy",), command="true"))
    assert W.mutation_order_ok(steps) != []
