<!-- trust/report_template.md -- copy to analysis/answers/<slug>.md, ONE FILE PER QUESTION.
     Append; never rewrite. Andy reviews this with `git diff`, not by re-reading prose.
     The chat reply quotes ONE number from here and links the file. -->

# Q: <the question, in Andy's words, quoted verbatim from the message that asked it>

Asked: <date>  ·  Status: OPEN | ANSWERED | ABANDONED

## Answer

<One sentence. The number, its units, its statistic. Nothing else.>

## Measurement

<Paste the provenance banner and the table VERBATIM from stdout. Do not retype, do not
reformat, do not trim columns. If it is too wide to paste, it is too wide to have run.>

```
== PROVENANCE v1 ==
...
== END PROVENANCE ==
...
== END RUN ledger=... ==
```

Reproduce:

```
<the exact command, copy-pasteable, including the interpreter>
```

## What this does NOT say

<The specific claims a reader would wrongly infer. Name the population the mask actually
kept, the parameters that were defaults, and anything the statistic is not.>

## Interpretation (unverified)

<Separated on purpose, and last, so it can be skipped or shot down without re-auditing
the measurement above. Nothing in this section may be quoted as a result.>

## Cost declared

- numbers computed in-session: <n>
- numbers carried from a subagent or an earlier session, not re-derived: <n>  <- each must be listed under UNVERIFIED: in the reply
- parameters chosen by the assistant unasked: <name = value, why, what it excluded>  (or: none)
- estimated audit time: <minutes>

## Log

- <date> — <what changed and why. Append. Never delete a superseded entry.>
