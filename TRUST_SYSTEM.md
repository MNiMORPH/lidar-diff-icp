# A trust system for this repo

**Problem being solved:** on 2026-08-25, roughly 80% of a working day went to auditing and
correcting the assistant's output instead of doing science. The binding constraint is not
that the assistant is wrong often – it is that Andy cannot tell *which* claims are wrong
without checking all of them, so the audit cost is paid on every claim, including the
correct ones.

**Design constraint:** ten memory files of prior corrections were reproduced anyway.
Advice does not bind. Every load-bearing piece here must work without the assistant's
judgement in the moment, because the assistant's judgement in the moment is the thing
that failed.

Everything described here is written and tested. Nothing is enabled. See
[What Andy has to do](#what-andy-has-to-do-to-turn-this-on).

---

## Part A. The system

Three mechanisms, in adoption order. Each names the failure class it closes and whether it
BLOCKS or merely REMINDS.

### Tier 1 – `Stop` hook, `core` level. **BLOCKS.** Adopt first.

`trust/hooks/stop_gate.py`, wired to the `Stop` and `SubagentStop` events. It reads the
session transcript, takes the final assistant message, and exits 2 – which Claude Code
feeds back to the model and refuses the stop – if any of:

| Code | Rule | Failure class closed |
|---|---|---|
| N1 | A number in the reply appears nowhere in this session's Bash output | Numbers typed rather than pasted |
| N2 | A number appears only in a subagent report and is not listed under an `UNVERIFIED:` line | The "SE = 10.0 mm" quote |
| P | A run declared a parameter with `src="MINE"` and the reply does not name it | The silent `--minn 20` |
| S | An input file of a quoted run has changed on disk since the run | Stale-input substitution |

Every block message ends with **GROUND TRUTH**: the actual resolved input paths, their
declared roles, and the true size of every selection, taken from the run ledgers. So a
reply that says "the CSF comparison on stable ground" is confronted, at the moment of
writing, with `z_after.npy IS: ... the VENDOR class-2 points (NOT our CSF)` and
`selection stable: keeps 92,404/355,600 (26.0%)`.

**Measured cost.** Replayed over all 2,798 delivered assistant replies in the failed
session `0804ccd6`: `core` blocks **100 replies (4%)**, of which 90 are N1. Inspection of a
random sample of the N1 hits shows they are overwhelmingly genuine – ratios computed
mentally (`1.65, 1.80, 1.54` derived in-head from `0.017/0.011`), a `~1.75 h` flight-line
separation carried from memory, and a bounding box arithmetically constructed in prose.
The one false-positive class found (significant-figure rounding, `~100,000 m³` against a
printed `100,234`) is fixed: hedged numbers are allowed to round.

**Runtime cost:** 0.93 s on a 113 MB transcript. Per turn, once.

### Tier 2 – same hook, `normal` level. **BLOCKS.** Adopt when tier 1 stops firing.

Adds two checks:

- **I** – a sentence asserts a cause or a consequence in the same breath as a number,
  with no `Interpretation:` label. Closes: interpretation emitted with the measurement,
  then defended.
- **V** – more than one table, a table wider than `TRUST_MAX_COLS` (default 4), or a
  reply longer than `TRUST_MAX_LINES` (default 40). Closes: eight columns and two tables
  when one column was asked for.

**Measured cost:** `normal` blocks **382 replies (14%)**, of which I accounts for 261.
That 9% is the real historical rate of fusing interpretation to measurement. Compliance
costs one line, so this is friction, not obstruction – but it is real friction, which is
why it is tier 2 and not tier 1.

### Tier 3 – `trust/provenance.py` in analysis scripts, plus the `PostToolUse` banner gate.

`trust/provenance.py` is the piece that makes the label come from the code that made the
number. A script declares `Run(question)`, then:

- `R.input(path, role=...)` – prints the resolved absolute path, size, mtime and a content
  digest, with a required sentence saying what the numbers in it *mean*. **Raises** if the
  role is missing.
- `R.param(name, value, src="andy"|"repo"|"MINE")` – **raises** unless `src` is one of the
  three, and **raises** if `src="MINE"` has no `why=`. An invented parameter can still be
  used; it cannot be used *silently*.
- `R.mask(name, m, defn=...)` – prints what the selection actually keeps, and warns above
  90%. `R.cuts(name, report)` records the criterion-by-criterion report that
  `reference_cells()` already returns and that the failing script threw away into `_`.
- `R.column(name, definition)` – `R.table()` **raises** on any column not defined.
- `R.done()` writes a JSON ledger that tier 1's P and S checks read.

`trust/hooks/banner_gate.py` (PostToolUse, matcher `Bash`) exits 2 when a `.py` under
`analysis/` or `scripts/` prints a numeric table with no banner. **This REMINDS rather
than BLOCKS** – PostToolUse cannot undo a tool that already ran. What it buys is that the
correction happens in the same turn, before the number reaches prose. Do not oversell it.

`TRUST_LEVEL=strict` then additionally blocks any unit-bearing number that did not come
from a banner-carrying run. **Do not enable strict until scripts carry banners** – its
historical block rate is meaningless, because no historical script had a banner.

**A worked example is committed:** `trust/example_instrumented.py` is
`analysis/ridgelines/gen2_csf_compare.py` with the science untouched and `Run()` calls
added. Run it and read the banner. It reports, unprompted, that the invented
`slope_max=90.0` removed `slope >= 90 deg=0` cells – the parameter disabled a criterion
entirely – and that the population called "stable" is 26.0% of the tile.

### Tier 4 – `trust/Makefile.example`, `trust/check_ledgers.py`. **REMINDS.** Optional.

Narrow. See rejections below.

### The report template. **REMINDS.** Convention, not enforcement.

`trust/report_template.md`. One file per question under `analysis/answers/<slug>.md`,
appended to, reviewed by `git diff` rather than by re-reading prose. The chat reply quotes
one number and links the file. The template puts *Interpretation (unverified)* last and
*What this does NOT say* before it.

---

## Part B. How to include the human's time in the calculus

### The exchange rate, computed rather than asserted

An unverified claim costs the assistant nothing to emit and costs Andy the full price to
check. From the failed session: one invented threshold (`--minn 20`) cost about an hour.
Asking a one-line question costs about 30 seconds of his attention.

That gives a break-even. Inventing a parameter rather than asking is correct only when the
probability of guessing the value he would have chosen exceeds

    30 s / 60 min = 0.008

so the assistant would need to be **99.2% certain** to justify not asking. But a gap is,
by definition, a place where the assistant is *not* certain – if it were certain, there
would be no gap to fill. **A gap is therefore always below break-even. The rule needs no
judgement: at a gap, ask.** This is a derivation, not a maxim, and it is the form in which
it survives contact with the moment.

The same arithmetic kills the "should I compute this or state it" question. Re-running a
script costs seconds of wall clock. Andy auditing one unverified number costs minutes.
The inequality never reverses. **There is no case in which emitting an uncomputed number
is correct**, which means the rule is an absolute rather than a threshold – and absolutes
do not require judgement to apply.

### Making the cost fall on the assistant, at the moment of writing

Right now the cost of an unverified artifact falls entirely on Andy. The Stop hook moves
it. When N1 fires, the turn does not end; the assistant must go and compute the number
before it can finish. The cost of the shortcut is charged in the only currency the
assistant has – **additional required work inside its own turn, before delivery** – and it
is charged automatically, at exactly the moment of the choice, not in hindsight.

This is the whole answer to *"how does this survive the fact that I do not experience the
passage of his time?"* It does not survive as empathy, and any proposal that relies on the
assistant feeling the cost should be rejected on sight. It survives as three mechanical
properties: (1) the cost is **computed by a hook from a ledger**, not estimated by the
assistant; (2) it is **charged in-turn**, so the shortcut is not cheaper for the assistant
either; (3) the escape hatch is **`TRUST_OFF=1`, Andy's environment variable**. There is
deliberately no in-band override, because an override the assistant can reach is an
override the assistant will reach.

### The unit: declared audit cost

The `Cost declared` block in `trust/report_template.md` carries, per answer: numbers
computed in-session, numbers carried in and not re-derived, parameters chosen unasked, and
estimated audit minutes. It is filled from the run ledger, not from self-report, so it
cannot be understated. Its purpose is triage: Andy reads the cost line and decides whether
to audit, instead of auditing to discover the cost.

### What to decline to produce, on time-cost grounds alone

1. **Any column not asked for.** Each extra column is audit surface with zero requested
   value, and it buries the requested column. Eight columns is not thoroughness; it is a
   transfer of eight audits for one answer.
2. **Any subagent number in prose that was not re-derived here.** Either re-run it or tag
   it `UNVERIFIED:`. Relaying is not reporting.
3. **Any interpretation in the same message as a first measurement.** It costs two audits
   instead of one, and if the number is wrong the reasoning audit was wasted entirely.
4. **Any figure of a number that has not been verified.** A figure is the most expensive
   thing to audit per bit of information it carries.
5. **Unrequested robustness checks, sensitivity analyses, and model comparisons.** These
   were the scope-creep failure. They read as diligence and bill as audit.
6. **Any re-derivation of something already recorded in a memory file.** Cite the file.
7. **Any result stated at more precision than the run printed.** Rounding silently is a
   small lie with a full-price audit.

### The honest limit

The gate checks that a number was *produced by a command in this session*. It does not
check that the command computed the right statistic. It would have caught "SE = 10.0 mm"
(N2, subagent-only), and it would have caught the invented `--minn 20` (P), and it makes
the vendor/CSF substitution visible in the block footer – but it cannot, in principle,
catch a correctly-pasted number from a wrongly-specified computation. Nothing here
substitutes for Andy reading the banner. What it does is make the banner exist, put it
next to the number, and refuse to let the number travel without it.

Second limit: in a long session, the corpus is large enough that a short number may match
by coincidence, which weakens N1 over time. `TRUST_WINDOW=40` restricts the corpus to the
40 most recent tool results and makes the gate stricter. It is off by default because it
also rejects legitimate reuse of an earlier result.

---

## Rejected, with reasons

- **`pint`.** Rejected. Unit-aware arrays would not have caught a single failure on the
  list. None of the six was a unit error. Adopting it means touching every array in the
  repo for a class of bug that is not occurring.
- **`pandera`.** Rejected. Dataframe assertions are written by the same judgement that
  failed. The `--minn 20` cell deletion would have been caught only by an assertion the
  assistant would have had to think to write – and the assistant that writes `--minn 20`
  unasked is not the assistant that writes `assert cover.max() > 0.5`. The narrow useful
  version, "print n before and after every cut", is already `R.mask()`, which raises if you
  skip it.
- **Snakemake, and `make` repo-wide.** Rejected as a requirement. It would mean
  restructuring roughly 60 ad-hoc analysis scripts, and it would not have caught the
  substitution that motivated it: the vendor/CSF swap happened *inside* a script that
  legitimately reads both files, which no dependency graph can see. The part worth keeping
  is staleness detection, delivered instead by `R.input()`'s digest plus check S, in about
  15 lines. `trust/Makefile.example` is provided for the narrow case that earns it – an
  answer re-derived more than twice.
- **Quarto / Jupyter inline execution.** Rejected as the primary mechanism, kept as the
  ideal. It is correct that a document holding `{python} f"{se:.1f}"` cannot hold a typed
  digit. But the failure surface is the **chat reply**, not a document, and Quarto cannot
  reach it. The Stop hook enforces the same invariant where the damage occurs. Use Quarto
  for anything that becomes a paper figure or a manuscript number.
- **A `UserPromptSubmit` hook injecting the standing rules each turn.** Rejected. It is
  the ten memory files again, re-injected more often. Re-reading advice is what already
  failed; it consumes context and binds nothing.
- **Any assistant-invokable override token.** Rejected on principle, stated above.

---

## What Andy has to do to turn this on

Nothing in this repo was modified. `.claude/settings.json` was not touched. New files only:

```
TRUST_SYSTEM.md
trust/__init__.py
trust/provenance.py            # the module scripts import
trust/check_ledgers.py         # sweep: which runs read a file that has since changed
trust/report_template.md       # answer-file template
trust/Makefile.example         # optional, tier 4
trust/settings.hooks.example.json
trust/example_instrumented.py  # worked example on real data
trust/hooks/stop_gate.py       # the Stop hook -- the piece that binds
trust/hooks/banner_gate.py     # the PostToolUse banner reminder
```

To enable:

1. Read `trust/settings.hooks.example.json` and **merge its `hooks` and `env` blocks by
   hand** into `.claude/settings.local.json`. Start with `"TRUST_LEVEL": "core"`.
2. Add `.trust/` to `.gitignore` if the run ledgers should not be tracked. (Tracking them
   is defensible: they are the provenance record.)
3. Verify before trusting it:
   ```
   ./.venv/bin/python trust/hooks/stop_gate.py --selftest      # 12 cases
   ./.venv/bin/python trust/hooks/banner_gate.py --selftest    #  5 cases
   ./lidar-icp/bin/python trust/provenance.py --selftest
   env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python trust/example_instrumented.py
   ```
4. To turn it off at any time: `TRUST_OFF=1`, or `"TRUST_LEVEL": "off"`.

Nothing was committed.
