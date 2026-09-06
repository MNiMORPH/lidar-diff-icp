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

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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

    # THE BASE PIPELINE TAKES NO CLOUD FLAGS. Its only cloud reader is `base`, which gets
    # them from the Site record (run_all_sites --only <site>), so nothing in the default
    # plan needs --gen1/--gen2. Every step that does is alongside, in a group -- which is
    # why this looks at the plan WITH the group reached over to.
    assert not [s for s, _, _ in cmds if {"gen1", "gen2"} & set(s.needs)], \
        "the pipeline gets its clouds from sites.py, not from flags"
    cmds = W.plan("data/derived/somewhere", with_groups=("vegetation_correction",))
    needs_gen2 = [(s, c, m) for s, c, m in cmds if "gen2" in s.needs]
    assert needs_gen2, "the alongside apparatus does read the gen2 cloud"
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
    assert W.script_of(conv) == ("src/lidar_diff_icp/steps/convexity_dod_landcover.py",)

    s = W.Step("x", produces=("o.npy",), requires=(), command="python scripts/a.py --tile t",
               code=("src/lidar_diff_icp/b.py",))
    assert W.script_of(s) == ("scripts/a.py", "src/lidar_diff_icp/b.py")


def test_a_step_run_as_a_module_still_names_its_source():
    """A `-m pkg.mod` command contains no .py path, so the path regex alone finds NOTHING
    and the step loses CODE-STALE detection -- a failure that looks like a product which
    is simply never out of date. That is not hypothetical: the three ridgeline steps moved
    into lidar_diff_icp.steps on 2026-09-06 and would have gone silently unwatched.

    Every source named must also EXIST, or the resolution is wrong in a way that leaves
    code_time() reading a missing file and reporting fresh.
    """
    import os
    s = W.Step("x", produces=("o.npy",), requires=(),
               command="python -m lidar_diff_icp.steps.trace_ridgelines {tile}")
    assert W.script_of(s) == ("src/lidar_diff_icp/steps/trace_ridgelines.py",)

    for step in W.STEPS:
        for src in W.script_of(step):
            assert os.path.exists(os.path.join(_REPO, src)), (
                f"step {step.name!r} names a source that does not exist: {src}")


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
    # Extended 2026-09-06 to take in the FEEDERS (Andy: "if this group is the feeder for
    # the vegetation correction, then they should also sit alongside"). The membership test
    # is what decides which steps are pipeline and which are apparatus, so the rule is
    # written here rather than inferred: a step belongs to the group when its products are
    # read ONLY by q2cover.py -- the correction's own library module -- and by the
    # correction's steps. slope, ridge_mask, convexity and curvature look similar and do
    # NOT belong: refcells.py reads them for the strict stable population that 20+ scripts
    # use, so they are pipeline.
    assert members == {"class2_spread", "q2_fit", "dod_cover", "lod_cover",
                       "pfs_cover", "gen1_angles", "beam_table", "nearground",
                       "nearground_split", "canopy_struct"}


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


def test_a_group_is_omitted_until_the_pipeline_reaches_over_to_it(tmp_path):
    """Andy, 2026-09-06: the vegetation-correction chain sits ALONGSIDE the pipeline, and
    the pipeline reaches over to use it optionally. So the default plan must not contain it
    -- the default is what the project claims, and a correction that measured WORSE than
    doing nothing on open ground must not look like part of producing the DoD."""
    members = {s.name for s in W.STEPS if s.group == "vegetation_correction"}
    base = {s.name for s, _, _ in W.plan(tmp_path)}
    assert not (base & members), "an opt-in group must be absent by default"
    reached = {s.name for s, _, _ in W.plan(tmp_path,
                                            with_groups=("vegetation_correction",))}
    assert reached - base == members
    assert base < reached


def test_an_unknown_group_is_refused_not_ignored(tmp_path):
    with pytest.raises(ValueError, match="unknown group"):
        W.plan(tmp_path, with_groups=("no_such_module",))


def test_a_group_something_outside_it_depends_on_is_refused(tmp_path):
    """A group whose products a base step requires is not optional AT ALL, whatever it is
    flagged -- omitting it by default would leave the pipeline unbuildable. That is a
    declaration error, so it raises whichever way the group is asked for."""
    steps = W.STEPS + (W.Step("downstream", produces=("z.npy",),
                              requires=("dod_cover_q2.npy",), command="true"),)
    with pytest.raises(ValueError, match="not optional at all"):
        W.plan(tmp_path, steps=steps)


def test_completeness_is_the_first_step():
    """It gates everything: a truncated fetch does not look like an error from inside the
    tile, and two sites shipped products built from one before it was noticed. Ordering it
    first is what makes --plan put the question before the answer."""
    assert [s.name for s in W.order()][0] == "completeness"


def test_completeness_applies_no_threshold():
    """The step RECORDS the ratio; it does not decide what ratio is acceptable. A cut-off
    nobody stated is exactly the kind of invented parameter that turns a judgement into a
    silent one."""
    step = next(s for s in W.STEPS if s.name == "completeness")
    assert "RECORDED, never judged" in step.note
    from lidar_diff_icp import completeness
    import inspect
    src = inspect.getsource(completeness)
    for op in ("ratio >", "ratio <", "ratio >=", "ratio <="):
        assert op not in src, f"completeness.py applies a threshold: {op}"


# --- the runner: it executes the graph, and it stops -------------------------------------

def _toy(tmp_path, n=3):
    """A tiny linear graph of steps that touch files, so a run can be observed."""
    steps = []
    for i in range(n):
        prev = () if i == 0 else (f"f{i-1}.npy",)
        steps.append(W.Step(f"s{i}", produces=(f"f{i}.npy",), requires=prev,
                            command=f"touch {tmp_path}/f{i}.npy"))
    return tuple(steps)


def test_the_runner_executes_in_dependency_order(tmp_path):
    ran, failure = W.run(tmp_path, steps=_toy(tmp_path), verbose=False)
    assert failure is None
    assert [s.name for s in ran] == ["s0", "s1", "s2"]
    for i in range(3):
        assert (tmp_path / f"f{i}.npy").exists()


def test_the_runner_stops_at_the_first_failure_and_says_what_was_not_attempted(tmp_path):
    """Andy's choice, 2026-09-06, and the conservative one: a later step consuming a
    half-written product turns one failure into a corrupted tile that LOOKS built."""
    steps = list(_toy(tmp_path, 4))
    steps[1] = W.Step("s1", produces=("f1.npy",), requires=("f0.npy",), command="exit 3")
    ran, failure = W.run(tmp_path, steps=tuple(steps), verbose=False)
    assert [s.name for s in ran] == ["s0"]
    assert failure is not None and failure[0].name == "s1" and failure[1] == 3
    # and it really stopped -- s2 and s3 never ran
    assert not (tmp_path / "f2.npy").exists()
    assert not (tmp_path / "f3.npy").exists()


def test_the_runner_refuses_before_running_anything_when_an_argument_is_missing(tmp_path):
    """Discovering at step 7 that step 8 wanted --gen2 means seven steps of point-cloud
    work thrown away. The check is free, so it happens first."""
    steps = _toy(tmp_path, 2) + (
        W.Step("hungry", produces=("z.npy",), requires=("f1.npy",),
               command="touch {gen2}", needs=("gen2",)),)
    with pytest.raises(ValueError, match="nothing was run"):
        W.run(tmp_path, steps=steps, verbose=False)
    assert not (tmp_path / "f0.npy").exists(), "it must refuse BEFORE running step 1"


def test_a_dry_run_executes_nothing(tmp_path):
    ran, failure = W.run(tmp_path, steps=_toy(tmp_path), dry_run=True, verbose=False)
    assert ran == [] and failure is None
    assert not (tmp_path / "f0.npy").exists()


def test_the_runner_skips_what_is_already_current(tmp_path):
    """only_stale is the point of tracking staleness: do not redo work that is right."""
    steps = _toy(tmp_path, 3)
    W.run(tmp_path, steps=steps, verbose=False)                 # build everything
    ran, _ = W.run(tmp_path, steps=steps, verbose=False)        # nothing left to do
    assert ran == []
    ran, _ = W.run(tmp_path, steps=steps, only_stale=False, verbose=False)
    assert [s.name for s in ran] == ["s0", "s1", "s2"], "--force runs them anyway"


def test_runnable_is_pure(tmp_path):
    """The selection can be shown before anything runs, so it executes nothing itself."""
    sel = W.runnable(tmp_path, steps=_toy(tmp_path))
    assert [s.name for s, _, _, _ in sel] == ["s0", "s1", "s2"]
    assert not any((tmp_path / f"f{i}.npy").exists() for i in range(3))


def test_only_is_a_filter_over_the_full_ordering_not_a_smaller_graph(tmp_path):
    """FOUND BY VALIDATING ON A REAL TILE, 2026-09-06. Passing a subset via steps= makes
    order() refuse -- a subset's requirements are produced OUTSIDE it, so `slope` alone
    fails on z_after.npy, which base makes. "Rebuild just this step" is a normal thing to
    want, and was impossible."""
    steps = _toy(tmp_path, 3)
    W.run(tmp_path, steps=steps, verbose=False)                  # build everything
    sel = W.runnable(tmp_path, steps=steps, only=("s1",), only_stale=False)
    assert [s.name for s, _, _, _ in sel] == ["s1"], "ordering is preserved, s0 filtered out"


def test_only_checks_a_filtered_steps_inputs_on_disk_rather_than_scheduling_them(tmp_path):
    """The guard that makes --only safe: running one step must never silently run its
    producer. A "rebuild just this step" that quietly rebuilds the DoD is worse than one
    that refuses."""
    steps = _toy(tmp_path, 3)
    with pytest.raises(ValueError, match="are not scheduled, and these are absent"):
        W.runnable(tmp_path, steps=steps, only=("s1",))          # f0.npy does not exist yet
    W.run(tmp_path, steps=steps, verbose=False)                  # now it does
    assert W.runnable(tmp_path, steps=steps, only=("s1",), only_stale=False)


def test_only_rejects_an_unknown_step_name(tmp_path):
    with pytest.raises(ValueError, match="unknown step"):
        W.runnable(tmp_path, steps=_toy(tmp_path), only=("no_such_step",))


# --- the declared hand-run tools ----------------------------------------------------------

def test_the_tools_registry_matches_the_tools_directory():
    """The registry and analysis/tools/ must agree, or "what produced this figure?" goes
    back to being unanswerable -- which is how a static sweep nearly deleted two live
    producers on 2026-09-05."""
    import glob
    on_disk = set(glob.glob("analysis/tools/*.py"))
    assert set(W.TOOLS) == on_disk, (
        f"declared but absent: {sorted(set(W.TOOLS) - on_disk)}; "
        f"present but undeclared: {sorted(on_disk - set(W.TOOLS))}")


def test_every_tool_says_what_it_is_for_in_its_own_words():
    """The purpose is each script's OWN first docstring line, not a description written
    into the registry. If a tool cannot say what it is for, that is the finding, and this
    test makes it visible rather than letting the registry paper over it."""
    import ast
    for path, purpose in W.TOOLS.items():
        assert purpose.strip(), f"{path} has no stated purpose"
        doc = ast.get_docstring(ast.parse(open(path).read())) or ""
        first = doc.strip().splitlines()[0].rstrip() if doc.strip() else ""
        assert purpose == first, (
            f"{path}: the registry has drifted from the script's own docstring.\n"
            f"  registry: {purpose}\n  docstring: {first}")


def test_a_tool_is_not_also_a_declared_step():
    """The two declarations answer different questions -- a Step is scheduled and gated, a
    tool is run by hand -- so a script must not claim to be both."""
    in_steps = {t for s in W.STEPS for t in s.command.split() if t.endswith(".py")}
    assert not (set(W.TOOLS) & in_steps), set(W.TOOLS) & in_steps
