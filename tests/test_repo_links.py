"""Every path this repository names in prose must exist.

PHASE 0 OF analysis/REORGANIZATION_PLAN.md, and the reason it comes first: 173 markdown
citations name a script by path, and the reorganization moves 143 scripts. Without this
test a move silently orphans a citation, and the document keeps claiming a reproduction
path that is no longer there. With it, every break fails loudly in the same commit.

It also catches the drift that accumulated before anyone was watching: twelve paths were
already dangling when this was written, seven of them from the pre-package layout
(src/pipeline.py, src/coreg.py, src/detect.py), which had been wrong long enough that
nobody noticed.
"""
import os
import re
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#: A path in prose only counts if it looks like one of ours: a repo-relative .py under a
#: known source root. Bare filenames are excluded -- too ambiguous to act on.
CITE = re.compile(r"(?<![\w/.])((?:analysis|scripts|src|ground_control|trust)"
                  r"/[A-Za-z0-9_./-]+\.py)")


def _tracked_markdown():
    out = subprocess.run(["git", "-C", REPO, "ls-files", "*.md"],
                         capture_output=True, text=True).stdout.split()
    return [os.path.join(REPO, m) for m in out]


#: A cited path may legitimately not exist, in exactly two cases, and BOTH MUST BE SAID OUT
#: LOUD on the same line as the citation:
#:
#:     (removed ...)   a historical record naming a file that is now only in git history
#:     (planned ...)   a design document naming a file it proposes to create
#:
#: Anything else is a broken link. The marker is deliberately explicit rather than inferred:
#: a silent exemption is how the twelve dangling paths accumulated in the first place, and a
#: reader who meets `foo.py (removed 2026-09-05, see git history)` knows where to look,
#: while a reader who meets a bare dead path does not.
EXEMPT = re.compile(r"\((?:removed|planned)\b[^)]*\)")


def citations():
    """[(markdown, cited_path, exempt)] over every tracked .md."""
    found = []
    for m in _tracked_markdown():
        with open(m, errors="ignore") as fh:
            for line in fh:
                hits = CITE.findall(line)
                if not hits:
                    continue
                ex = bool(EXEMPT.search(line))
                for hit in hits:
                    found.append((os.path.relpath(m, REPO), hit, ex))
    return found


def test_every_script_named_in_a_document_exists():
    """A citation is a promise that the reader can go and look. When the file is gone the
    document is asserting a reproduction path that does not exist."""
    dangling = sorted({(m, p) for m, p, ex in citations()
                       if not ex and not os.path.exists(os.path.join(REPO, p))})
    assert not dangling, (
        "these documents name scripts that do not exist. Fix the path, or mark the "
        "citation '(removed ...)' or '(planned ...)' on the same line if that is the "
        "truth:\n  " + "\n  ".join(f"{m}  ->  {p}" for m, p in dangling))


def test_the_lint_actually_finds_citations():
    """Guard against the regex silently matching nothing, which would make the test above
    pass by looking at an empty set -- the vacuous-test failure mode that bit twice today."""
    c = citations()
    assert len(c) > 50, f"only {len(c)} citations found; the pattern is probably broken"
    assert any(p.startswith("analysis/") for _, p, _ in c)


def test_the_exemption_is_narrow():
    """The escape hatch must not swallow the check. Most citations are live links, and if
    the exempt fraction ever grew large the test would be passing by excusing itself."""
    c = citations()
    ex = [x for x in c if x[2]]
    assert len(ex) / len(c) < 0.10, (
        f"{len(ex)} of {len(c)} citations are exempted; the marker is being overused")


def test_an_exempt_citation_still_has_to_be_marked_on_its_own_line():
    """Pins the mechanism: the marker is per line, not per file, so one '(removed)' cannot
    excuse every dead path in a long document."""
    assert EXEMPT.search("`foo/bar.py` (removed 2026-09-05, see git history)")
    assert EXEMPT.search("`foo/bar.py` (planned; Phase 1 creates it)")
    assert not EXEMPT.search("`foo/bar.py` was removed last week")
