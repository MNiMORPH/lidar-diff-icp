"""The derived-product graph: well-formed, correctly ordered, and honest about staleness.

These exist because both failures they cover happened for real. The convexity producer
rewrites ridgecrest_pixels.npz from scratch and the curvature producer augments it in
place, so running them in the wrong order silently drops three columns. And a q2 fit read
from a beam_offset_table older than corrections.json carries superseded registration --
which no amount of reading the file tells you, because the file looks fine.
"""
import json
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
    with pytest.raises(ValueError, match="produced by no step"):
        W.order(bad)


def test_a_cycle_is_rejected():
    cyc = (W.Step("a", produces=("a.npy",), requires=("b.npy",), command="true"),
           W.Step("b", produces=("b.npy",), requires=("a.npy",), command="true"))
    with pytest.raises(ValueError, match="cannot order the graph"):
        W.order(cyc)


def _touch(d, name, when=None, content=b""):
    """Write a file with given CONTENT and mtime.

    Content matters now: the staleness check arbitrates with a hash, so bumping an mtime
    without changing bytes is deliberately NOT a change. A test that means "this input
    changed" must write different bytes.
    """
    p = os.path.join(d, name)
    with open(p, "wb") as fh:
        fh.write(content)
    if when is not None:
        os.utime(p, (when, when))
    return p


def test_state_reports_missing_stale_and_ok(tmp_path):
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="true")
    steps = (step,)

    _touch(d, "in.npy", when=1000, content=b"v1")
    assert W.state(d, steps)["t"][0] == "MISSING"

    _touch(d, "out.npy", when=2000, content=b"out")
    assert W.state(d, steps)["t"][0] == "OK"

    # a REAL change: new bytes and a newer mtime
    _touch(d, "in.npy", when=3000, content=b"v2 -- different")
    kind, detail = W.state(d, steps)["t"]
    assert kind == "STALE" and "in.npy" in detail


def test_stale_is_detected_across_the_real_chain(tmp_path):
    """corrections.json newer than beam_offset_table.parquet -- the case that bit us."""
    d = str(tmp_path)
    now = time.time()
    for f in W.BASE_INPUTS:
        _touch(d, f, when=now - 100, content=f.encode())
    for s in W.order():
        for f in s.produces:
            _touch(d, f, when=now - 50, content=f.encode())
    assert all(k == "OK" for k, _ in W.state(d).values())

    # corrections.json genuinely re-solved: new bytes, newer mtime
    _touch(d, "corrections.json", when=now, content=b"re-solved constants")
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


def test_script_of_finds_the_producing_source_and_extras():
    conv = next(s for s in W.STEPS if s.name == "convexity")
    assert W.script_of(conv) == ("analysis/ridgelines/convexity_dod_landcover.py",)

    s = W.Step("x", produces=("o.npy",), requires=(), command="python scripts/a.py --tile t",
               code=("src/lidar_diff_icp/b.py",))
    assert W.script_of(s) == ("scripts/a.py", "src/lidar_diff_icp/b.py")


def test_every_step_resolves_at_least_one_source_file():
    """A step whose command names no tracked script cannot be code-checked at all."""
    for s in W.STEPS:
        assert W.script_of(s), f"{s.name}: no source file parsed from its command"


def test_code_state_flags_a_product_older_than_its_source(tmp_path, monkeypatch):
    d = str(tmp_path)
    src = tmp_path / "scripts" / "made_up.py"
    src.parent.mkdir(parents=True)
    src.write_text("# a producer\n")

    step = W.Step("t", produces=("out.npy",), requires=(),
                  command=f"python scripts/made_up.py")
    steps = (step,)

    _touch(d, "out.npy", when=1000)
    monkeypatch.setattr(W, "code_time", lambda p: 500.0)     # source older than the product
    assert W.code_state(d, steps) == {}

    monkeypatch.setattr(W, "code_time", lambda p: 2000.0)    # source NEWER than the product
    assert W.code_state(d, steps) == {"t": ["scripts/made_up.py"]}


def test_a_missing_product_is_not_reported_as_code_stale(tmp_path, monkeypatch):
    """MISSING is state()'s to report; code_state must not double up on it."""
    step = W.Step("t", produces=("absent.npy",), requires=(), command="python scripts/x.py")
    monkeypatch.setattr(W, "code_time", lambda p: 9e9)
    assert W.code_state(str(tmp_path), (step,)) == {}


def test_base_inputs_are_code_checked_too(tmp_path, monkeypatch):
    """A pipeline.py change invalidates dod/lod at every tile, and no Step produces them.

    Since 2026-09-04 the ground-q calibration curve is checked the same way: difference_dem
    takes each cell's percentile from it, so a NEW curve invalidates every tile's ground --
    and nothing inside a tile directory reveals that its ground came from an older one.
    """
    d = str(tmp_path)
    for f in W.BASE_INPUTS:
        _touch(d, f, when=1000)
    monkeypatch.setattr(W, "code_time", lambda p: 500.0)
    assert W.base_code_state(d) == []
    monkeypatch.setattr(W, "code_time", lambda p: 2000.0)
    assert W.base_code_state(d) == sorted(tuple(W.BASE_CODE) + tuple(W.BASE_GLOBAL_INPUTS))


def test_an_identical_rewrite_does_not_cascade_staleness(tmp_path):
    """A producer that re-runs and writes the SAME bytes must not invalidate the world.

    At Carlton this was not hypothetical: re-running the DoD rewrote a byte-identical
    z_after, and a plain mtime check then demanded a six-minute PyForestScan rebuild that
    would have changed nothing.
    """
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="python x.py")
    steps = (step,)

    inp = tmp_path / "in.npy"
    inp.write_bytes(b"the same content")
    os.utime(inp, (1000, 1000))
    _touch(d, "out.npy", when=2000)

    assert W.state(d, steps)["t"][0] == "OK"          # bootstraps the manifest

    inp.write_bytes(b"the same content")              # rewritten, identical, mtime bumped
    os.utime(inp, (3000, 3000))
    assert W.state(d, steps)["t"][0] == "OK", "identical bytes must not read as a change"


def test_genuinely_new_content_still_reports_stale(tmp_path):
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="python x.py")
    steps = (step,)

    inp = tmp_path / "in.npy"
    inp.write_bytes(b"first")
    os.utime(inp, (1000, 1000))
    _touch(d, "out.npy", when=2000)
    assert W.state(d, steps)["t"][0] == "OK"

    inp.write_bytes(b"DIFFERENT content")
    os.utime(inp, (3000, 3000))
    kind, detail = W.state(d, steps)["t"]
    assert kind == "STALE" and "in.npy" in detail


def test_the_manifest_is_written_beside_the_products(tmp_path):
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="python x.py")
    _touch(d, "in.npy", when=1000)
    _touch(d, "out.npy", when=2000)
    W.state(d, (step,))
    man = json.load(open(os.path.join(d, W.MANIFEST)))
    assert "in.npy" in man and "sha256" in man["in.npy"] and "content_time" in man["in.npy"]


def test_rerunning_a_step_clears_its_staleness_even_if_the_bytes_match(tmp_path):
    """Outputs are judged by when they were PRODUCED, inputs by when their content changed.

    Using content_time on both sides was a real bug, found by using the tool: Battle
    Creek's gen1_csf_angles.npz was rebuilt after corrections.json changed, came out
    byte-identical because the change did not reach it, kept its old content_time, and so
    reported stale forever however many times it was regenerated.
    """
    d = str(tmp_path)
    step = W.Step("t", produces=("out.npy",), requires=("in.npy",), command="python x.py")
    steps = (step,)

    _touch(d, "in.npy", when=1000, content=b"v1")
    _touch(d, "out.npy", when=2000, content=b"result")
    assert W.state(d, steps)["t"][0] == "OK"

    _touch(d, "in.npy", when=3000, content=b"v2 -- changed")     # a real input change
    assert W.state(d, steps)["t"][0] == "STALE"

    # re-run the step: same bytes out, but produced NOW
    _touch(d, "out.npy", when=4000, content=b"result")
    assert W.state(d, steps)["t"][0] == "OK", \
        "a step re-run against the new input is current, whatever bytes it produced"


# --- optional module groups ---------------------------------------------------------------

def test_the_vegetation_correction_is_a_leaf_so_it_can_be_switched_off():
    """A group may only be called optional if NOTHING outside it requires its products.
    The correction's outputs are dod_cover_q2.npy and lod_cover_q2.npy -- not dod.npy and
    lod.npy -- and this is what pins that separation: if a base step ever came to require a
    corrected product, the shipped DoD would silently depend on a correction that measured
    WORSE than doing nothing on open ground."""
    assert W.group_is_a_leaf("vegetation_correction")
    members = {s.name for s in W.STEPS if s.group == "vegetation_correction"}
    assert members == {"class2_spread", "q2_fit", "dod_cover", "lod_cover",
                       "cover_calibration"}


def test_every_step_in_an_optional_group_is_itself_optional():
    """A mandatory step inside a switchable module is a contradiction: --skip-optional and
    --skip-group would disagree about whether the product is required."""
    for s in W.STEPS:
        if s.group in W.GROUPS:
            assert s.optional, f"{s.name} is in group {s.group} but is not optional"


def test_every_declared_group_has_a_description_and_members():
    for g in W.GROUPS:
        assert W.GROUPS[g].strip(), g
        assert [s for s in W.STEPS if s.group == g], f"group {g} has no steps"
    for s in W.STEPS:
        assert s.group == "" or s.group in W.GROUPS, f"{s.name}: undeclared group {s.group!r}"


def test_skipping_a_group_drops_exactly_its_steps(tmp_path):
    full = {s.name for s, _, _ in W.plan(tmp_path)}
    cut = {s.name for s, _, _ in W.plan(tmp_path,
                                        skip_groups=("vegetation_correction",))}
    assert full - cut == {s.name for s in W.STEPS if s.group == "vegetation_correction"}
    assert cut  # the base pipeline survives


def test_an_unknown_group_is_refused_not_ignored(tmp_path):
    with pytest.raises(ValueError, match="unknown group"):
        W.plan(tmp_path, skip_groups=("no_such_module",))


def test_a_group_something_depends_on_cannot_be_skipped(tmp_path):
    """Refusing beats dropping: a plan missing a step another step needs is unbuildable,
    and a plan that cannot be run is worse than no plan."""
    steps = W.STEPS + (W.Step("downstream", produces=("z.npy",),
                              requires=("dod_cover_q2.npy",), command="true"),)
    with pytest.raises(ValueError, match="cannot be skipped"):
        W.plan(tmp_path, steps=steps, skip_groups=("vegetation_correction",))
