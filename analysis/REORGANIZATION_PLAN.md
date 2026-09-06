# Reorganizing `analysis/`: a plan

Written 2026-09-06 at Andy's instruction, after the pipeline itself was reduced to six
declared steps and everything else was found to be undeclared. Hard now for better later.

## The problem, measured

    analysis/*.py                    157      in the workflow graph        16
    analysis/*.md                     67      undeclared                  143

Of the 143 undeclared:

    imported by other code             7   de facto libraries in the wrong place
    cited in a .md                   107   the reproduction path for a written finding
    touched since 2026-09-01          18   live tools from the current work
    none of the above                 11   the only clear delete candidates

It is not that the scripts are junk. **Nothing says what they are**, so they cannot be told
apart from junk — which is how a static sweep of mine nearly deleted two live producers, and
how three separate name-matching scans over-reported by 10-100x in one day.

## Two hazards that decide the sequence

    39 scripts manipulate sys.path      they depend on their own depth in the tree
    61 import a sibling by bare name    they depend on their neighbours' location
    173 .md citations name a script     they depend on the path

So a move is never one file. **The sibling-import graph defines the units that can move**, and
the citation graph defines what must be rewritten with them. Both are computable, which is
what makes this safe rather than a judgement call per file.

## Target layout

    src/lidar_diff_icp/        the library                                  (unchanged)
    scripts/                   drivers the graph declares                   (unchanged)
    ground_control/            a distinct instrument, self-contained        (unchanged)
    trust/                     provenance tooling                           (unchanged)
    analysis/
      tools/                   hand-run producers, audits, figure makers -- DECLARED in a
                               TOOLS registry, so provenance is answerable
      modules/
        vegetation_correction/ the alongside module's own scripts, beside the group that
                               declares them
      investigations/          the record of questions asked, one directory per question,
                               each carrying its own findings .md
      findings/                cross-cutting conclusions and the single FRAME

## Phases, each verifiable before the next

### Phase 0 — safety nets, BEFORE anything moves

  0a  citation lint: a test asserting every `(analysis|scripts|src)/*.py` named in a `.md`
      exists. Fixes the 12 already broken and makes every later move fail loudly.
  0b  import lint: a test that every `analysis/**/*.py` still parses and resolves its
      imports. This is what catches the 61 sibling imports when a move breaks one.
  0c  baseline: record the import graph and the citation graph as a committed artifact, so
      "did the move change anything?" is answerable rather than argued.

  GATE: 0a and 0b pass on the tree as it stands today.

### Phase 1 — make the structure explicit WITHOUT moving anything

  1a  compute the connected components of the sibling-import graph. These are the clusters
      that must move together. Publish the list.
  1b  TOOLS registry: declare the 18 live tools and the 7 imported libraries -- name, what
      it produces, what question it serves. Provenance becomes answerable here, before any
      file changes location.
  1c  map each cluster to its findings document, or record that it has none.

  GATE: every one of the 143 is in exactly one bucket, by evidence, with the rule written
  down. No file has moved.

### Phase 2 — move, one cluster at a time

  For each cluster, in one commit: `git mv` (history preserved) -> fix `sys.path` depth ->
  fix sibling imports -> rewrite the citations that name it -> full test suite green.

  GATE per cluster: 0a and 0b still pass, and the test suite is green. A cluster that cannot
  be made green is REVERTED, not patched around.

  Order: `tools/` first (smallest, most used), then `modules/vegetation_correction/`
  (already declared as a group, so its membership is known), then `investigations/`
  (largest, least coupled to the pipeline).

### Phase 3 — delete the residue

  The 11 that are uncited, unimported and untouched since before 2026-09-01. One at a time,
  with the evidence in the commit message, as the penetration sweep was done. Everything
  else is kept.

### Phase 4 — the documents

  Six FRAME files, each saying "read this FIRST", become one. Findings documents move to
  `findings/`. The citation lint from 0a keeps every path in them honest.

## What this does NOT do, and why

  - **No bulk classification-driven deletion.** My classifier was 17% wrong on a
    twelve-item hand check. Deletion happens in Phase 3, on 11 files, individually.
  - **No invented categories.** The clusters come from the import graph and the citation
    graph, not from my reading of what a script is "about".
  - **No move before the lints exist.** Phase 0 is what makes Phase 2 reversible in
    practice rather than in principle.

## Honest cost

Phase 0 is a few hours and pays for itself immediately. Phase 1 is mostly computation plus
the TOOLS entries. Phase 2 is the long one -- 143 files in maybe 15-25 clusters, each with a
test-green gate -- and it is where the "hard now" is. Phase 3 is an hour. Phase 4 is a day
of prose.

The alternative is what we have: a directory where 91% of the contents cannot be told from
debris, which has already cost real errors today.
