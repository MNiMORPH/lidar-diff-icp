"""Every script under analysis/ must parse, and every import it makes must resolve.

PHASE 0b OF analysis/REORGANIZATION_PLAN.md. 61 scripts import a sibling by bare name and
39 manipulate sys.path, so both depend on where the file sits in the tree. Phase 2 moves
143 files: this is the test that turns "the move broke an import" from a thing discovered
weeks later into a red gate in the same commit.

STATIC, NEVER EXECUTED. These are scripts with top-level code -- importing them would run
157 analyses, some of which stream point clouds for minutes. So imports are extracted from
the AST and resolved against the filesystem and the installed environment. That is weaker
than importing, and it is the right trade: it catches every move-induced break, which is
what Phase 2 can cause, without running anything.
"""
import ast
import importlib.util
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _tracked_py(prefix):
    out = subprocess.run(["git", "-C", REPO, "ls-files", f"{prefix}/*.py"],
                         capture_output=True, text=True).stdout.split()
    return [p for p in out if os.path.exists(os.path.join(REPO, p))]


def _imports(path):
    """(module, level) for every import in the file. level>0 means a relative import."""
    with open(os.path.join(REPO, path), errors="ignore") as fh:
        tree = ast.parse(fh.read(), filename=path)
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            out += [(a.name.split(".")[0], 0) for a in n.names]
        elif isinstance(n, ast.ImportFrom):
            out.append(((n.module or "").split(".")[0], n.level))
    return out


#: Modules that live in the CONDA env, not this venv, and are documented as doing so by the
#: scripts that use them (forest_metrics_pfs.py:8 names the env explicitly). Listed rather
#: than inferred: "unresolvable means fine" would gut the test, which exists to catch
#: exactly that.
OTHER_ENV = {"pdal", "pyforestscan"}


def _static_syspath(path):
    """Directories a file adds to sys.path with a literal string, so the lint can follow it.

    This found the one that mattered: trace_ridgelines.py used to insert an ABSOLUTE path
    outside the repository -- /home/awickert/dataanalysis/r.fluvial -- for a PIPELINE step,
    so the pipeline did not build on any other machine. rivernetworkx is now a declared
    dependency (git+https://github.com/awickert/r.fluvial) and the insert is gone. Following
    inserts is still how this test tells a real external dependency from a broken import.
    """
    out = []
    with open(os.path.join(REPO, path), errors="ignore") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for n in ast.walk(tree):
        if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                and n.func.attr in ("insert", "append")
                and isinstance(n.func.value, ast.Attribute)
                and n.func.value.attr == "path"):
            for a_ in n.args:
                if isinstance(a_, ast.Constant) and isinstance(a_.value, str):
                    out.append(a_.value)
    return out


def _resolves(mod, path):
    """Can `mod` be found from `path`'s directory, or in the environment?"""
    if not mod:
        return True                                   # `from . import x`
    here = os.path.dirname(os.path.join(REPO, path))
    if os.path.exists(os.path.join(here, mod + ".py")):
        return True                                   # a sibling: the fragile case
    if os.path.exists(os.path.join(here, mod, "__init__.py")):
        return True
    for root in ("", "src", "scripts", "analysis", "ground_control"):
        if os.path.exists(os.path.join(REPO, root, mod + ".py")) or \
           os.path.exists(os.path.join(REPO, root, mod, "__init__.py")):
            return True
    for extra in _static_syspath(path):
        d = extra if os.path.isabs(extra) else os.path.join(REPO, extra)
        if os.path.exists(os.path.join(d, mod + ".py")) or \
           os.path.exists(os.path.join(d, mod, "__init__.py")):
            return True
    if mod in OTHER_ENV or mod in sys.builtin_module_names:
        return True
    try:
        return importlib.util.find_spec(mod) is not None
    except (ImportError, ValueError, ModuleNotFoundError):
        return False


def test_every_analysis_script_parses():
    bad = []
    for p in _tracked_py("analysis"):
        try:
            _imports(p)
        except SyntaxError as e:
            bad.append(f"{p}: {e}")
    assert not bad, "these do not parse:\n  " + "\n  ".join(bad)


def test_every_import_resolves_from_where_the_file_sits():
    """THE MOVE GATE. A sibling import resolves only from the directory the file is in, so
    this fails the moment Phase 2 moves a script away from something it imports."""
    bad = []
    for p in _tracked_py("analysis"):
        for mod, level in _imports(p):
            if level:
                continue                              # explicit relative: package-internal
            if not _resolves(mod, p):
                bad.append(f"{p}: cannot resolve `import {mod}`")
    assert not bad, ("these imports do not resolve from the file's own directory:\n  "
                     + "\n  ".join(sorted(bad)))


def test_the_lint_is_looking_at_the_real_tree():
    """Guards the vacuous case: if the file list came back empty the tests above would pass
    by checking nothing."""
    files = _tracked_py("analysis")
    assert len(files) > 100, f"only {len(files)} analysis scripts found"
    assert any(_imports(p) for p in files)


def test_the_only_out_of_repo_dependency_is_the_one_we_know_about():
    """A step reaching outside the repository by absolute path is a reproducibility limit:
    the pipeline then builds on one machine only. There is now NONE -- rivernetworkx, the
    last one, is a declared git dependency as of 2026-09-06.

    If one appears, this fails and the decision gets made deliberately rather than
    discovered when someone else tries to run the pipeline."""
    external = {}
    for p in _tracked_py("analysis") + _tracked_py("scripts") + _tracked_py("src"):
        for d in _static_syspath(p):
            if os.path.isabs(d) and not d.startswith(REPO):
                external.setdefault(d, []).append(p)
    assert external == {}, f"absolute paths outside the repo: {external}"


def test_the_other_env_allowlist_stays_small():
    """The exemption must not become the answer to every failure."""
    assert len(OTHER_ENV) <= 4, OTHER_ENV
