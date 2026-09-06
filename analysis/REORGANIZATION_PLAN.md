# Reorganizing `analysis/`: a plan that runs autonomously

Written 2026-09-06 at Andy's instruction. Hard now for better later.

The first draft of this plan needed a person at four points — naming investigation folders,
mapping clusters to findings, choosing deletions, folding documents. **This version removes
every one of those, or defers it explicitly.** What remains is a program with a dry-run mode.

## The problem, measured

    analysis/*.py                    157      declared in the workflow graph   16
    analysis/*.md                     67      undeclared                      143

      imported by other code           7   de facto libraries in the wrong place
      cited in a .md                 107   the reproduction path for a finding
      touched since 2026-09-01        18   live tools from the current work
      none of the above               11   the only clear delete candidates

They are not junk. **Nothing says what they are**, so they cannot be told from junk — which
is how a static sweep of mine nearly deleted two live producers this morning.

## Why it can be autonomous: every destination is COMPUTED

The rule below is total, deterministic, and uses only facts already on disk. It is applied
in order; first match wins. **No step asks a question, and no step encodes my opinion of
what a script is "about".**

    1  named in a workflow Step command or `code`     -> stays put (declared)
    2  imported by src/, or by >= 2 other scripts     -> analysis/lib/
    3  named in a vegetation_correction Step command  -> analysis/modules/vegetation_correction/
    4  calls savefig() AND touched since 2026-09-01   -> analysis/tools/
    5  cited by exactly one .md                       -> analysis/investigations/<doc-stem>/
    6  cited by two or more .md                       -> analysis/investigations/_shared/
    7  anything else                                  -> analysis/investigations/_unsorted/

Rule 5 is the load-bearing trick: **the investigation's folder is named by the document that
cites it**, not by my reading of the script. If the finding was worth writing down, its name
already exists. `_unsorted/` is where the rule admits it does not know, rather than guessing.

## Why it is safe: clusters, gates, and revert-on-red

**Clusters.** 39 scripts manipulate `sys.path` and 61 import a sibling by bare name, so a
move is never one file. The connected components of the sibling-import graph are the units
that move. A component's destination is the destination of its majority member; ties break
to the lowest-numbered rule. Deterministic, and it never splits a component.

**Gates.** After every cluster, in this order: citation lint, import lint, full pytest.

**Revert-on-red is the autonomy policy.** A cluster whose gate fails is `git revert`-ed,
recorded as BLOCKED with the failure text, and the run CONTINUES to the next cluster. It does
not stop and it does not patch around. Blocked clusters are the report.

**Resumable.** `analysis/.reorg_manifest.json` records per cluster: destination, status
(pending / done / blocked), and the failing gate. Re-running continues where it stopped.

## The program

    reorganize_analysis.py --plan     compute and print every destination; move nothing
    reorganize_analysis.py --apply    execute, one commit per cluster, gated
    reorganize_analysis.py --report   what moved, what is blocked, what is unsorted

`--plan` is committed as an artifact BEFORE `--apply` runs, so the intended end state is
reviewable as a diff rather than discovered afterwards.

## Phases

### Phase 0 — safety nets. Must pass on today's tree before anything moves.

    0a  citation lint  every (analysis|scripts|src)/*.py named in a .md exists
                       (fixes the 12 already broken)
    0b  import lint    every analysis/**/*.py parses and resolves its imports
    0c  baseline       import graph + citation graph committed as an artifact

### Phase 1 — the program, and its plan. Nothing moves.

    1a  write scripts/reorganize_analysis.py (planned; this phase creates it) implementing the rule and the clustering
    1b  commit the output of --plan: all 143 destinations, and the cluster list
    1c  TOOLS registry entries for whatever lands in analysis/tools/ and analysis/lib/

### Phase 2 — `--apply`. Autonomous, gated, resumable.

One commit per cluster: `git mv` -> fix `sys.path` depth -> fix sibling imports -> rewrite
citations naming it -> run the three gates -> commit, or revert and mark BLOCKED.

### Phase 3 — report only. NOT autonomous, and deliberately so.

Produces `analysis/DELETION_CANDIDATES.md`: the 11 uncited/unimported/untouched, plus
whatever landed in `_unsorted/`, each with its evidence. **Deletes nothing.** Deletion of a
tracked record is Andy's call, and a 17%-wrong classifier is not allowed to make it.

### Phase 4 — deferred, and named as deferred.

Folding six FRAME documents into one is prose about what is currently true. It is not
mechanical, so it stays a human task rather than being faked by a rule.

## What I will need from Andy: nothing, until Phase 3

Phases 0-2 run start to finish without a decision. Phase 3 produces a list and stops.
Phase 4 waits.

## Honest cost and failure modes

Phase 0 is a few hours. Phase 1 is the program, half a day, and is where care is repaid.
Phase 2 is machine time plus whatever the blocked list costs afterwards.

The likely failure: a cluster whose `sys.path` depth cannot be mechanically fixed because it
walks upward a fixed number of levels. Those revert, land on the blocked list, and are fixed
by hand later. I expect a handful, not dozens — but I will not know until `--plan` runs, and
the plan is committed before `--apply` precisely so that estimate is checkable.
