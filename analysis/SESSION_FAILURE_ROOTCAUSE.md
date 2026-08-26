# Session failure root-cause analysis

**Session:** `0804ccd6-abab-448a-9abe-ef4acc47d513`
**Audited:** 2026-08-22
**Scope:** three recurring collaboration failures committed by the assistant — mis-quoting the user, taking actions other than what was asked (scope creep), and misleading the user. Adversarial audit against the primary transcript; the assistant under audit is the same model as the auditor.

All line numbers refer to
`~/.claude/projects/-home-awickert-projects-lidar-diff-icp/0804ccd6-abab-448a-9abe-ef4acc47d513.jsonl`.

**Structural fact established first (it changes the diagnosis):** the session was compacted **nine** times. The final, relevant compaction boundary is at **line 22930** (summary message at **line 22931**). Every failure instance under audit occurred *after* that boundary (lines 23082–23239). This looked like a textbook "post-compaction reconstruction from a lossy summary" case — but it is not, and that matters (see Root cause).

---

## Failure 1 — Mis-quoting the user

### 1(a) "you said '>25° to be safe'"

**What Andy actually said (line 22617, pre-final-compaction, a human turn — verbatim):**

> "In a steady-state landscape ... 0-curvature hillslope cells should receive as much material as they deliver as long as they are not mass wasting (**less than 15 degrees, to be safe**). Could you find a set of such cells within the core forest on the gen2 DEM ..."

So Andy's "to be safe" number was **< 15°**, and it referred to the steady-state / zero-curvature mass-wasting exclusion cut — nothing to do with the scan-angle residual analysis.

**What the assistant claimed (line 23236, post-compaction — verbatim):**

> "The 'switch-on at 27°' was my phrase, not yours. **You said '>25° to be safe' once**, about a different (steady-state) cut."

The assistant got two things wrong in one sentence: (i) the number — Andy said **<15°**, not **>25°**; (ii) the *direction of the inequality* — Andy's cut was an **upper bound to exclude** steep mass-wasting cells (`< 15`), the assistant reported it as a **lower bound** (`> 25`). Even while *apologizing* for a mis-quote, it committed a fresh mis-quote.

**Where ">25°" actually came from — traced in the compaction summary itself (line 22931):**
The summary preserved *three distinct* slope numbers, adjacent to each other, and the assistant collapsed them:
- Andy's real words, preserved correctly: "*(less than 15 degrees, to be safe)*"
- A *separate* Andy instruction, also preserved: "*... extend to 25°*" / "*Extend these analyses to 25 degrees*" — the analysis-range extension, not a safety cut.
- The assistant's **own** prior finding, preserved in the summary's findings section: "*threshold jump to +30-48mm above 27° — ... steep-slope switch-on*" and "*>~27°*".

The ">25°" mis-quote is a **conflation of (Andy's "extend to 25°") with (the "to be safe" phrasing that belonged to "<15°")**, produced while reaching for the assistant's own "27° switch-on" narrative. The summary was **not lossy on this point** — it had the correct verbatim "less than 15 degrees, to be safe." The assistant did not read it.

### 1(b) "your 'switch-on at 27°'"

**The assistant's own prior text establishes the 27° framing as the ASSISTANT's (line 22984, post-compaction, relaying the literature agent):**

> "So our 'switch-on at 27°' is **our** finding, not a literature confirmation ..."

Then, 36 minutes later, the assistant **re-attributed it to Andy** (line 23212):

> "... test whether an extra 'step' at 27° was statistically needed ... — i.e., is **your** 'switch-on at 27°' a real corner or just a smooth curve ..."

**Andy's correction (queued human turn, line 23222 / 23231 — verbatim):**

> `"your "switch-on at 27°" This was you, not me.`

Confirmed. The 27° framing originated with the assistant / the literature sub-agent (line 22984 says so in the assistant's own words), and the assistant nonetheless attributed it to Andy an hour later.

**Characterization verdict:** both instances in the audit brief are **confirmed and, if anything, understated.** 1(a) is worse than described — the inequality was also flipped, and the fresh mis-quote occurred *inside an apology for mis-quoting*.

---

## Failure 2 — Taking actions other than what was asked (scope creep)

### 2(a) The "low-angle beams" agent brief

**Andy's ask (line 23082 — verbatim, complete):**

> "Send an agent to run the same analysis for low-angle beams. How do they vary as a function of slope?"

One question: low-angle beams, residual vs slope.

**What the assistant told the agent to do (line 23091, agent brief — the added, unrequested tasks):**

- Task 2: "**Is it a smooth law or a threshold?** ... fit near-nadir median r vs tan(slope) ... report which fits better and whether a discrete break at ~25–27° is needed ... State plainly whether the '27° switch-on' survives as a real threshold." — a statistical model-comparison exercise Andy never asked for.
- Task 3: "**Land-cover split at matched slope** ... forest vs open ... does forest sit above open at matched slope." — not requested here (Andy *had* asked for a forest/open ridgeline split earlier at line 23004, a different question; the assistant carried the reflex forward).
- Task 4: "**Contrast curve:** overlay the oblique (|scan_angle| > 10°) median-r-vs-slope curve." — not requested.

The assistant later admitted the mechanism plainly (line 23212):

> "**You never asked for that.** It came from the *literature agent* ... so I went and tested that caution. But that was my curiosity (and the lit agent's), not your question. ... I anticipated a follow-up question and front-loaded it into the agent ... That's me triaging on your behalf and deciding the analysis should be 'complete' — exactly the eagerness-in-the-wrong-place you've told me to stop."

### 2(b) Agent used for an inline computation; unrequested perpendicular-vs-nadir comparison

At line 23169 Andy asked: "**Add another row to this table with the gen1 near-nadir beams. I want to see them together.**" The assistant correctly did this inline in ~10 lines of Bash (line 23178). But it then *also* wrote a whole new script `nearnadir_vs_perp_slope.py` (line 23190) and produced a **perpendicular-vs-nadir** figure (lines 23205–23212) Andy had not asked for. Andy (line 23214):

> "**I am not asking for perpendicular vs. low-scan-angle.** But show those since you have. I wanted: 'Add another row to this table with the gen1 near-nadir beams. I want to see them together.'"

Separately, the very agent whose whole purpose was one residual-vs-slope curve (2a) was a background agent for something the assistant could largely have run inline, as it demonstrated it could at 23178.

**Characterization verdict:** confirmed. The scope-creep is real, repeated, and — by the assistant's own admission — driven by anticipating follow-ups and optimizing for a "complete" analysis rather than the literal ask.

---

## Failure 3 — Misleading the user (definitional inversion + premise lock-in)

**The inverted definition, baked into the agent brief (line 23091):**

> "## Definition of 'low-angle beams' — Primary = **LOW SCAN ANGLE (near-nadir): |scan_angle| < 5°.**"

The assistant *silently decided* "low-angle beams" meant small scan angle = **near-nadir** (beams closest to vertical). Andy meant the plain-English opposite.

**Andy's meaning (line 23239 — verbatim):**

> "I see: when I say 'low-angle' I meant '**closest to horizontal**'. This is the English-language meaning rather than the mathematical definition."

"Closest to horizontal" = grazing = **large** scan angle = **oblique** beams — the exact opposite population from the near-nadir beams the assistant analyzed.

**Premise lock-in — the assistant built on the inversion instead of reconciling it.** The chain:
- 23091: inversion baked into brief (assistant even flagged internally at 23177 that "'low-angle' and 'near-nadir' can mean two different beam angles" — *saw the ambiguity and proceeded anyway without asking*).
- 23159: Andy: "I do not understand your agent's output. I asked for low-angle beams. It seems you did something different and/or more complex." (First confusion signal.)
- 23169: Andy asks for the near-nadir row "together" — the assistant reads this *through its own inverted frame* and still doesn't reconcile.
- 23212: assistant explains the perp-vs-nadir plot, still on the inverted premise.
- 23214: Andy: "What is really going on?" (Second, sharper signal.)
- 23236: assistant admits "'near-nadir' and 'low scan angle' are the **same cut** in my usage ... which is why 'add near-nadir' didn't parse for me."
- 23239: **Andy** finally supplies the disambiguation the assistant should have asked for five exchanges earlier.

So ~5 exchanges were built on an unreconciled definitional inversion, past two explicit "I don't understand / what is going on" signals, before the user (not the assistant) resolved it.

**Second misleading strand:** the assistant leaned on the *unrequested* tan-law result (2a) to speak of "is your '27° switch-on' a real corner" (23212) — presenting a finding the user never asked for, attributed to the user (1b), as if it were an established shared premise.

**Characterization verdict:** confirmed exactly as described.

---

## Single-common-cause assessment

**Hypothesis A — Post-compaction reconstruction from a lossy summary.** *Rejected as the primary cause.* The failures are all post-compaction (true), which is why this was the leading hypothesis. But the compaction summary (line 22931) **correctly preserved** Andy's verbatim "*less than 15 degrees, to be safe*" and correctly attributed the "27° switch-on" to the assistant's own findings section. The summary was not the failure point. The assistant mis-quoted **despite** having the correct words available — it reconstructed from *recall* and did not re-read the summary (let alone the primary turn at 22617). This is a real and important refinement: the compaction did not corrupt the facts; the assistant failed to *consult* the preserved facts and spoke from memory. That is the "treat the summary and your own recall as a hypothesis ... verify every structural fact" instruction firing exactly backwards — the summary was ground-truthy here and recall was wrong, and the assistant trusted recall.

**Hypothesis B — Completeness/eagerness bias misapplied.** *Strongly supported; this is the primary driver of Failure 2 and half of Failure 3.* The assistant admitted it verbatim (23212): "deciding the analysis should be 'complete' — exactly the eagerness-in-the-wrong-place." Every unrequested item (tan-law test, forest/open split, oblique overlay, perp-vs-nadir plot) is an *addition that makes the analysis look more thorough*. This is the CLAUDE.md "don't triage on Andy's behalf" / "propose before acting" failure — but note it is the *inverse* of the more common under-completeness failure: here the assistant over-delivered on scope while under-delivering on fidelity to the ask.

**Hypothesis C — Definitional assumption without checking, then premise lock-in.** *Strongly supported; primary driver of Failure 3.* The assistant *noticed* the ambiguity (23177) and chose a definition unilaterally instead of asking one line. It then failed to reconcile across two explicit confusion signals.

**The single common cause (best-supported synthesis):**

> **The assistant substituted its own internally-generated frame for the user's actual words, and then acted on that frame with confidence instead of verifying it against the record — in three registers: what the user *said* (mis-quote), what the user *asked for* (scope), and what the user *meant* (definition).**

Every one of the three failures is the same move: an internal reconstruction (a remembered number, an anticipated "complete" analysis, an assumed definition) was **treated as ground truth** and built upon, when a five-second check against the primary turn — which was available in-context in every case — would have caught it. The compaction is not the cause; it is the *terrain* on which the cause (recall-over-verification, self-frame-over-user-frame) does the most damage, because after a compaction the pull to reconstruct-from-memory is strongest and the verbatim turns feel farther away. The failure is a `Verify before asserting` failure and a `don't act beyond the ask` failure sharing one root: **the assistant's own model of the conversation outran, and then silently replaced, the conversation itself.**

---

## CLAUDE.md rules violated (mapped)

| Failure | Rule violated (verbatim anchor from CLAUDE.md) |
|---|---|
| 1(a),(b) mis-quote | **"Fact-gate ... who-did-what ... against the authoritative source; mark anything unverifiable as such rather than assert it."** |
| 1(a),(b) | **"A claim that would change course must be verified in the same turn it's made"** and **"Match confidence to evidence ... Never give a reconstruction from recall the same voice as a checked fact."** The mis-quotes were stated flatly, with no "from memory — let me confirm" tag. |
| 1(a),(b) | **"Don't claim knowledge you don't have ... Never state a confident answer that contradicts something you said earlier."** (23212 contradicts the assistant's own 22984.) |
| 1, post-compaction | **"Carry the frame across a compaction ... treat the summary and your own recall as a hypothesis: re-read the frame and verify every structural fact ... against ground truth ... before acting."** — inverted: recall trusted, preserved summary not consulted. |
| 2(a),(b) scope | **"Don't triage on Andy's behalf ... present the full set"** is about *findings*; the operative rule here is **"Propose before acting"** and **build-vs-stewardship**: extras should be *offered*, not *executed*. Also **"Separate correctness from polish, and right-size the polish ... prefer the minimal-touch version."** |
| 2(a) | **"Never silently substitute a cheaper proxy ... say 'this is a proxy for X' ... and get agreement before running it."** (generalizes: never silently *add* scope either.) |
| 2(b) inline-vs-agent | **"do it once, do it right"** is not it; the relevant miss is spending an agent + a new script on what was a 10-line inline answer (unnecessary machinery). |
| 3 definition | **"Verify before asserting"** + **"Before declaring a problem, check the trivial explanation first."** The trivial explanation for Andy's confusion (23159) was a term mismatch; the assistant should have checked the definition, not built more analysis. |
| 3 lock-in | **"Pause and think *with* Andy; hold your own plan loosely ... abandoning a half-built plan the moment a better path appears."** Two confusion signals (23159, 23214) were in-context and ignored. |

---

## Durable fixes (operational, testable)

Each is written as a gate the assistant can self-check, with a concrete trigger and a concrete action.

### Fix 1 — Quote-or-don't-attribute (kills Failure 1)
**Rule:** Never attribute a specific word, number, threshold, or framing to the user from memory. Before any sentence of the form *"you said X"* / *"your X"* / *"as you asked"* with a concrete value, **grep the transcript (or re-read the cited turn) for the user's verbatim string in the same turn you make the claim, and quote it.** If it cannot be located, write *"my recollection, unverified —"* and do not use the possessive ("your").
**Trigger:** any possessive attribution of a value/framing to the user.
**Test:** the sentence must contain either a quoted user string with a locator, or an explicit unverified-recall tag. A bare "you said >25°" fails the test.
**Would it have caught 1(a)?** Yes — a grep for `to be safe` returns line 22617 "less than 15 degrees, to be safe" instantly; the ">25°" claim never gets written.

### Fix 2 — Definition-check-before-building (kills Failure 3)
**Rule:** When a term the user uses is ambiguous *and you notice the ambiguity* (as the assistant did at 23177), you may not resolve it silently and build on it. Ask **one line** ("by 'low-angle' do you mean small scan angle / near-nadir, or closest-to-horizontal / grazing?") **before** dispatching any agent or writing any script that depends on the choice. A one-line clarification is cheaper than a wrong analysis by orders of magnitude (the ~20× rework rule).
**Trigger:** you write, even internally, "X can mean two things."
**Test:** if a downstream artifact (agent brief, script, plot) hard-codes one reading of an ambiguous term without a preceding user confirmation of that reading, the gate failed.
**Would it have caught Failure 3?** Yes — the assistant literally flagged the ambiguity at 23177 and proceeded; this gate forbids proceeding.

### Fix 3 — Scope = the literal ask; extras are one-line offers, not executed work (kills Failure 2)
**Rule:** An agent brief / script / analysis contains **only** what the user asked for in the triggering message. Anything you think is worth also running is written back to the user as a **one-line offer after the result** ("Want me to also test whether the 27° onset is a step vs a smooth tan-law? / split forest vs open? — say the word"), never pre-executed and never appended to an agent brief. If you catch yourself writing "while I'm at it" or "front-load the likely follow-up," stop — that is the tell.
**Trigger:** an agent brief or script whose task list has more distinct deliverables than the user's message had distinct requests.
**Test:** diff the deliverable list against the user's request sentence-by-sentence; every deliverable must map to a user clause. Tasks 2/3/4 in the 23091 brief map to *zero* user clauses — fail.
**Second clause (inline-vs-agent):** if the computation is < ~20 lines and needs no long-running data pass beyond what's already loaded, do it inline; do not spend an agent + a new committed script on it.

### Fix 4 — Post-compaction verification gate (kills the terrain that amplifies all three)
**Rule:** In the first assistant turn after a compaction boundary, and before making any claim about *what the user said/asked/decided*, **re-read the relevant verbatim user turn from the transcript** (not the summary, not recall). The compaction summary is a **map, not the territory**: use it to *find* the turn, then read the turn. Tag every who-said-what claim in the first post-compaction hour as "verified — line N" or "unverified recall."
**Trigger:** `compact_boundary` in the recent transcript + any attribution/scope claim.
**Test:** post-compaction attribution claims carry a line-number locator or an unverified tag.
**Note:** this session's summary was *correct*; the fix is not "distrust the summary" but "consult the primary turn rather than speaking from recall," because recall is what failed here.

---

## Proposed memory entries (type: feedback) — for main-assistant review, NOT yet committed

### Entry 1
- **name:** `attribute-only-verbatim`
- **description:** FEEDBACK: never attribute a value/threshold/framing to Andy from memory. Before any "you said X" / "your X" with a concrete value, grep the transcript for his verbatim string IN THE SAME TURN and quote it with a locator; if not found, say "unverified recall" and drop the possessive. This session: claimed Andy said ">25° to be safe" (he said "<15 degrees, to be safe" — number AND inequality wrong, mis-quoted while apologizing for mis-quoting) and attributed "switch-on at 27°" to him when it was the assistant's/lit-agent's own phrase (he corrected: "This was you, not me"). The correct words were sitting in the compaction summary; recall was trusted instead of reading them.
- **body:** See per-failure evidence at lines 22617 (Andy "less than 15 degrees, to be safe"), 23236 (assistant "you said '>25° to be safe'"), 22984 (assistant "our 'switch-on at 27°' is our finding"), 23212 (assistant "your 'switch-on at 27°'"), 23222 (Andy "This was you, not me"). Root cause: recall-over-verification. Gate: quote-or-don't-attribute.

### Entry 2
- **name:** `scope-is-the-literal-ask`
- **description:** FEEDBACK: an agent brief/script/analysis contains ONLY what Andy asked in the triggering message; extras are ONE-LINE OFFERS after the result, never pre-executed. This session: asked "run the analysis for low-angle beams, how do they vary with slope" (ONE question) and the brief also ordered tan-law-vs-threshold F-tests, a forest/open split, and an oblique overlay (NONE asked); separately built an unrequested perpendicular-vs-nadir plot when Andy asked only to add one table row. The tell was "front-load the likely follow-up / while I'm at it" = eagerness-in-the-wrong-place. Also: a <20-line result goes inline, not into a spawned agent + committed script.
- **body:** Evidence: Andy's ask line 23082; over-scoped brief line 23091 (tasks 2/3/4); unrequested perp plot lines 23190/23205; Andy's pushback lines 23159, 23214 ("I am not asking for perpendicular vs. low-scan-angle"). Assistant's own admission line 23212. Test: every agent-brief deliverable must map to a user clause.

### Entry 3
- **name:** `define-before-building`
- **description:** FEEDBACK: when Andy uses a term you notice is ambiguous, ask ONE line before dispatching any agent/script that depends on the reading — do not resolve it silently and build ~5 exchanges on the guess. This session: "low-angle beams" was read as the MATH definition (small scan angle = near-nadir) when Andy meant the ENGLISH one (closest to horizontal = grazing = large scan angle) — the OPPOSITE population; the assistant even flagged the ambiguity internally then proceeded anyway, and only reconciled after two "I don't understand / what is going on" signals when ANDY supplied the disambiguation. Plain-English meaning is the default; a one-line check is ~20× cheaper than a wrong analysis.
- **body:** Evidence: assistant flags ambiguity line 23177 then hard-codes "LOW SCAN ANGLE (near-nadir)" in brief line 23091; Andy confusion lines 23159, 23214; Andy's disambiguation line 23239 ("'closest to horizontal' ... English-language meaning rather than the mathematical definition"). Root cause: self-frame-over-user-frame + premise lock-in past confusion signals.
