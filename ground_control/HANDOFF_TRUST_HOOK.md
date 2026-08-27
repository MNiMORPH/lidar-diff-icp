# Handoff: the trust hooks do not run when the session's working directory is a subdirectory

**Written 2026-08-27** by the `ground_control/` session, which is fenced to this directory and
therefore **cannot apply the fix itself**. Everything below was verified by command in that
session; each claim names the check that produced it.

**STATUS: FIXED AND VERIFIED 2026-08-27**, on Andy's explicit authorization to cross the
write boundary. The fix of §4 was applied to `.claude/settings.local.json` (untracked --
globally gitignored) and `trust/settings.hooks.example.json` (tracked). §7 records the proof
that the gate now bites. The rest of this document is kept as the diagnosis, because the
failure mode it describes is not specific to this repository and will recur wherever a session
is started in a subdirectory.

---

## 1. The symptom

Every `Bash` tool call in the session emitted:

```
PostToolUse:Bash hook blocking error from command: "python3 "$CLAUDE_PROJECT_DIR/trust/hooks/banner_gate.py"":
python3: can't open file '/home/awickert/projects/lidar-diff-icp/ground_control/trust/hooks/banner_gate.py':
[Errno 2] No such file or directory
```

## 2. The cause

`.claude/settings.local.json` invokes all three hooks as
`python3 "$CLAUDE_PROJECT_DIR/trust/hooks/<gate>.py"`.

`CLAUDE_PROJECT_DIR` resolves to the session's **primary working directory**, not to the
repository root. This session was started with its working directory set to
`ground_control/`, so the variable expanded to
`/home/awickert/projects/lidar-diff-icp/ground_control` and the hook path became
`.../ground_control/trust/hooks/banner_gate.py`, which does not exist.

The hooks live at `<repo>/trust/hooks/{banner_gate,stop_gate}.py` — verified present and
executable (`ls -la trust/hooks/`: both `-rwxrwxr-x`, 3913 and 18824 bytes).

**This is a configuration bug, not a code bug.** Both hooks pass their own selftests from the
repository root:

```
$ python3 trust/hooks/banner_gate.py --selftest     # exit 0, 5 of 5 OK
$ python3 trust/hooks/stop_gate.py   --selftest     # exit 0, 12 of 12 OK
```

## 3. What was actually lost — read this before deciding urgency

The two hooks are **not** equally affected, and the noisy one is the harmless one.

**`banner_gate.py` (PostToolUse) was a no-op anyway.** Its `main()` begins:

```python
if os.environ.get("TRUST_OFF") == "1" or os.environ.get("TRUST_LEVEL", "").lower() in ("off", "core"):
    return 0
```

`.claude/settings.local.json` sets `"TRUST_LEVEL": "core"`. So at the configured level this
hook returns 0 immediately **by design**, and the broken path cost nothing but noise on every
Bash call. Fixing the path will not make it start firing; it stays inert until `TRUST_LEVEL`
is raised to `normal` or `strict`.

**`stop_gate.py` (Stop / SubagentStop) is the one that enforces, and it is broken.** Its
`main()` returns early only at `LEVEL == "off"` (`trust/hooks/stop_gate.py:321`), so at `core`
it runs its full check — the `N1 N2 P S` tier, "is this number real", the tier
`trust/settings.hooks.example.json` measures as blocking 4% of replies. That check did not run
in this session.

**A precise statement of the failure mode, and its limit.** The broken invocation exits **2**:

```
$ echo '{"transcript_path":"/nonexistent","cwd":"'$PWD'"}' \
    | python3 "$CLAUDE_PROJECT_DIR/trust/hooks/stop_gate.py"; echo $?
python3: can't open file '.../ground_control/trust/hooks/stop_gate.py': [Errno 2] ...
2
```

Exit 2 is the same code the gate uses to mean *block and feed the reason back to the model*.
So the gate was not merely absent — it was replaced by something block-shaped carrying a
message about a missing file. **What Claude Code did with that downstream was not observed
from inside the session**, and this document does not claim it. What is verified is the exit
code and the message.

## 4. The fix

Make the hook path independent of where the session was started, by resolving the repository
root instead of trusting `CLAUDE_PROJECT_DIR` to be it:

```
python3 "$(git -C "${CLAUDE_PROJECT_DIR:-$PWD}" rev-parse --show-toplevel)/trust/hooks/stop_gate.py"
```

Verified in all three conditions:

```
CLAUDE_PROJECT_DIR=<repo>/ground_control  -> /home/awickert/projects/lidar-diff-icp   hook exists: YES
CLAUDE_PROJECT_DIR=<repo>                 -> /home/awickert/projects/lidar-diff-icp   hook exists: YES
CLAUDE_PROJECT_DIR unset                  -> /home/awickert/projects/lidar-diff-icp   hook exists: YES
```

End to end on a Stop payload, from a `ground_control/` working directory:

```
CURRENT  : python3: can't open file '.../ground_control/trust/hooks/stop_gate.py'   exit=2
PROPOSED : (no output)                                                              exit=0   <- ran, found nothing to block
```

### Files to edit

1. **`.claude/settings.local.json`** — three `command` strings (`Stop`, `SubagentStop`,
   `PostToolUse`). Replace `$CLAUDE_PROJECT_DIR/trust/hooks/` with the resolver above in each.
2. **`trust/settings.hooks.example.json`** — the same three, so the documented example stops
   reproducing the bug for the next person who copies it.

## 5. Three decisions that are NOT mine to make

Flagged rather than chosen, per the project's standing rule on invented defaults.

- **Whether to fix the path or the launch convention.** The alternative fix is "always start
  sessions at the repository root," which needs no edit. Against it: `ground_control/HANDOFF.md`
  explicitly scopes a session to this subdirectory, so the convention and the tooling are in
  direct conflict as they stand. Recommendation: fix the path, because it is the option that
  does not rely on remembering something.
- **Whether the gate should fail closed.** Right now a resolver failure produces a Python
  `can't open file` error, which is legible to a human but is not a trust message. A wrapper
  that emitted `TRUST GATE NOT RUNNING — <reason>` on exit 2 would make the failure
  self-describing. This is a real change in behaviour, so it is a proposal, not a fix.
- **Whether `PostToolUse`/`banner_gate` should be configured at all at `TRUST_LEVEL=core`,**
  given §3 shows it returns 0 immediately at that level. Keeping it is harmless once the path
  is fixed and means nothing has to change when the level is raised; removing it drops a hook
  that cannot fire. Either is defensible.

## 6. Verifying the fix

From a working directory *inside* a subdirectory — that is the case that fails today:

```
cd <repo>/ground_control
echo '{"transcript_path":"/nonexistent","cwd":"'$PWD'"}' \
  | python3 "$(git -C "${CLAUDE_PROJECT_DIR:-$PWD}" rev-parse --show-toplevel)/trust/hooks/stop_gate.py"
echo "exit=$?"      # expect 0 and no 'can't open file'
```

Then start a session with its working directory set to `ground_control/` and run any Bash
command: the `PostToolUse:Bash hook blocking error` line should be gone.

**A stronger check, worth doing once:** confirm the Stop gate actually bites after the fix, the
way a regression test is proven to bite. Set `TRUST_LEVEL=normal` and end a turn on a reply
containing a typed, unpasted number; `stop_gate.py --selftest`'s first case
(`typed number blocked`) is the behaviour to expect. A gate that has never been seen to fire
is indistinguishable from one that is still broken.


---

## 7. The fix, applied and proven to bite (2026-08-27)

Applied by replacing `$CLAUDE_PROJECT_DIR/trust/hooks/` with
`$(git -C "${CLAUDE_PROJECT_DIR:-$PWD}" rev-parse --show-toplevel)/trust/hooks/` in three
`command` strings in each of the two files. Both files re-parse as valid JSON; the diff is
6 changed lines each (3 replacements x 2), line counts unchanged at 49 and 59, no line-ending
churn.

**Proven to bite, the way a regression test must be.** The command string was read back out of
`.claude/settings.local.json` rather than retyped, and run from a `ground_control/` working
directory with `CLAUDE_PROJECT_DIR` set to that subdirectory -- the exact condition that
failed -- at the configured `TRUST_LEVEL=core`, against two synthetic transcripts: one whose
reply contains a typed, unpasted number, one clean.

```
                    reply=typed          reply=clean
OLD (broken)        exit=2               exit=2        <- "can't open file", both alike
NEW (fixed)         exit=2               exit=0        <- "TRUST GATE -- this reply cannot
                                                          be delivered as written."
```

The discrimination is the result. The broken hook returned 2 on **both** replies, so it was
not a gate that had stopped working -- it was a gate that had been replaced by one which
fired on everything with a message about a missing file. The fixed hook fires on the typed
number and passes the clean reply.

**Not done, still the maintainer's call:** the two proposals in §5 -- making the gate fail
closed with a self-describing `TRUST GATE NOT RUNNING` message rather than a Python
traceback, and whether to configure `PostToolUse`/`banner_gate` at all at `TRUST_LEVEL=core`
given §3 shows it cannot fire there. Only the path was changed.
