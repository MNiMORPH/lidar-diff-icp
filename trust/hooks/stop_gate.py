#!/usr/bin/env python3
"""Stop hook: refuse to let the turn end while the reply contains unearned claims.

This is the only piece of the trust system that does not depend on the assistant's
judgement in the moment, because it runs *after* the assistant has decided it is done.

It reads the session transcript, takes the final assistant message, and blocks the stop
(exit code 2, message on stderr, which Claude Code feeds back to the model) if any of:

  N1  A number in the reply does not appear verbatim in any Bash tool output from this
      session.                                    -> "typed, not pasted"
  N2  A number appears only in a subagent's report and was never re-derived here, and the
      reply does not list it under an `UNVERIFIED:` line.
  N3  A number carrying a physical unit came from a command that printed no provenance
      banner (trust/provenance.py).               -> unlabelled quantity
  P   A run in this session declared a parameter with src="MINE" (assistant-invented) and
      the reply does not name it.
  S   An input file of a run quoted in this session has changed on disk since the run.
  I   The reply states an interpretation alongside a number without labelling it.
  V   The reply carries more than one table, a table wider than TRUST_MAX_COLS, or is
      longer than TRUST_MAX_LINES.

Every block message ends with the actual input paths of every run in this session, so a
narration that contradicts the files it read is visible at the moment of writing.

Escape hatch is the HUMAN's, deliberately: `TRUST_OFF=1` in the environment, or
`"TRUST_LEVEL": "off"` in settings env. The assistant has no in-band override; that is
the point.

Install: see trust/settings.hooks.example.json.
Test:    python trust/hooks/stop_gate.py --selftest
"""
from __future__ import annotations

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BANNER_START = "== PROVENANCE v1 =="
LEDGER_RE = re.compile(r"== END RUN ledger=(\S+)")

# core   -- N1 N2 P S      : "is this number real?"        (measured 3.6% of replies)
# normal -- + I V           : + interpretation & volume     (measured 14%)
# strict -- + N3            : + provenance banner required for unit-bearing numbers
# off    -- nothing
LEVEL = os.environ.get("TRUST_LEVEL", "normal").lower()
MAX_COLS = int(os.environ.get("TRUST_MAX_COLS", "4"))
MAX_LINES = int(os.environ.get("TRUST_MAX_LINES", "40"))
# Corpus window: how many of the most recent tool results count as "this session's
# output". 0 = the whole session. A long session makes coincidental matches more likely
# (a 113 MB transcript contains most short numbers somewhere), so a window of ~40 makes
# the gate STRICTER, not looser. Set it if replies start slipping through.
WINDOW = int(os.environ.get("TRUST_WINDOW", "0"))

UNIT = r"(?:mm|cm|km|m|deg|°|%|px|pt|pts|cells|cell|sigma|σ|x|×)"
NUM_RE = re.compile(r"(?<![\w.\-])([+-]?\d[\d,]*(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![\w])")
UNIT_AFTER = re.compile(r"^\s*" + UNIT + r"\b")
VERSIONISH = re.compile(r"^\d+\.\d+\.\d+$")
INTERP_RE = re.compile(
    r"\b(because|means|shows that|demonstrat|driven by|caused by|"
    r"suggest|implies|indicat|explain|therefore|hence|consistent with|confirms|"
    r"refutes|so the|which is why)\b", re.I)
INTERP_LABEL = re.compile(r"^\s*(interpretation|interp|hypothesis)\b", re.I | re.M)
UNVERIFIED_LINE = re.compile(r"^\s*UNVERIFIED\b(.*)$", re.I | re.M)
FENCE_RE = re.compile(r"```.*?```", re.S)


# ------------------------------------------------------------------ transcript
def read_transcript(path):
    """-> (final_reply_text, {tool_name: [result_text, ...]})"""
    recs = []
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    recs.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    id2tool, results, tail = {}, {}, []
    for r in recs:
        msg = r.get("message") or {}
        content = msg.get("content")
        if r.get("type") == "assistant" and isinstance(content, list):
            for b in content:
                if b.get("type") == "tool_use":
                    id2tool[b.get("id")] = b.get("name", "?")
        if r.get("type") == "user" and isinstance(content, list):
            for b in content:
                if b.get("type") == "tool_result":
                    c = b.get("content")
                    if isinstance(c, list):
                        c = "\n".join(x.get("text", "") for x in c if isinstance(x, dict))
                    name = id2tool.get(b.get("tool_use_id"), "?")
                    results.setdefault(name, []).append(c if isinstance(c, str) else str(c))
    if WINDOW > 0:
        results = {k: v[-WINDOW:] for k, v in results.items()}
    # final reply = assistant text blocks after the last tool_result / user turn
    for r in reversed(recs):
        t = r.get("type")
        msg = r.get("message") or {}
        content = msg.get("content")
        if t == "assistant" and isinstance(content, list):
            txt = "\n".join(b.get("text", "") for b in content if b.get("type") == "text")
            if any(b.get("type") == "tool_use" for b in content):
                break
            if txt.strip():
                tail.append(txt)
        elif t == "user":
            break
    return "\n".join(reversed(tail)), results


# --------------------------------------------------------------------- numbers
def _norm(s):
    return s.replace(",", "").lstrip("+")


def corpus_floats(text):
    out = set()
    for m in NUM_RE.finditer(text):
        try:
            out.add(float(_norm(m.group(1))))
        except ValueError:
            pass
    return out


def present(tok, raw_corpus, float_corpus, hedged=False):
    """Is `tok` genuinely in the corpus -- verbatim, or as a legitimate rounding?

    `hedged` (the prose said ~, about, roughly) also permits significant-figure
    rounding, so "~100,000" is accepted against a printed 100,234.
    """
    n = _norm(tok)
    if n in raw_corpus:
        return True
    try:
        v = float(n)
    except ValueError:
        return True
    k = len(n.split(".")[1]) if "." in n else 0
    if any(round(c, k) == v for c in float_corpus):
        return True
    if hedged and v != 0:
        import math
        sig = len(n.replace("-", "").replace(".", "").rstrip("0")) or 1
        for c in float_corpus:
            if c == 0:
                continue
            d = sig - int(math.floor(math.log10(abs(c)))) - 1
            if round(c, d) == v:
                return True
    return False


def candidate_numbers(reply):
    """Numbers in the reply that constitute a quantitative CLAIM."""
    body = FENCE_RE.sub(" ", reply)
    out = []
    for m in NUM_RE.finditer(body):
        tok = m.group(1)
        n = _norm(tok)
        after = body[m.end():m.end() + 8]
        has_unit = bool(UNIT_AFTER.match(after))
        digits = n.replace(".", "").replace("-", "")
        line_start = body.rfind("\n", 0, m.start()) + 1
        if re.match(r"^\s*\d+[.)]\s", body[line_start:m.end() + 2]):
            continue                                   # list marker
        if VERSIONISH.match(n):
            continue
        if re.match(r"^\s*(-\d\d|/|:)", after) or body[max(0, m.start()-1):m.start()] in "/-:":
            continue                                   # date or path fragment
        if not has_unit and re.fullmatch(r"(19|20)\d\d", n):
            continue                                   # bare year
        if has_unit or "." in n or len(digits) >= 3:
            before = body[max(0, m.start() - 12):m.start()]
            hedged = bool(re.search(r"(~|\u2248|approx|about|roughly|order of)\s*$", before, re.I))
            out.append((tok, has_unit, m.start(), hedged))
    return out


def sentences(text):
    """Crude sentence split that also treats list items and table rows as units."""
    out = []
    for chunk in re.split(r"\n\s*\n|\n(?=\s*[-*|\d])", text):
        out += [s for s in re.split(r"(?<=[.;!?])\s+", chunk) if s.strip()]
    return out


# ---------------------------------------------------------------------- checks
def markdown_tables(reply):
    tables, cur = [], []
    for line in FENCE_RE.sub(" ", reply).splitlines():
        if line.strip().startswith("|"):
            cur.append(line)
        elif cur:
            tables.append(cur); cur = []
    if cur:
        tables.append(cur)
    return [(len(t), max(r.count("|") - 1 for r in t)) for t in tables if len(t) >= 2]


def session_ledgers(bash_text, cwd):
    out = []
    for m in LEDGER_RE.finditer(bash_text):
        p = m.group(1)
        p = p if os.path.isabs(p) else os.path.join(cwd, p)
        if os.path.exists(p):
            try:
                out.append((p, json.load(open(p))))
            except Exception:
                pass
    return out


def check(reply, results, cwd):
    v = []
    bash = "\n".join(results.get("Bash", []))
    bash_banner = "\n".join(c for c in results.get("Bash", []) if BANNER_START in c)
    agent = "\n".join(results.get("Agent", []) + results.get("Task", []))
    other = "\n".join(t for k, ts in results.items() if k not in ("Bash", "Agent", "Task") for t in ts)

    b_raw, b_f = _norm(bash), corpus_floats(bash)
    ban_raw, ban_f = _norm(bash_banner), corpus_floats(bash_banner)
    ag_raw, ag_f = _norm(agent), corpus_floats(agent)
    ot_raw, ot_f = _norm(other), corpus_floats(other)

    unver = " ".join(m.group(0) for m in UNVERIFIED_LINE.finditer(reply))
    unver_raw, unver_f = _norm(unver), corpus_floats(unver)

    typed, subagent, unbannered = [], [], []
    for tok, has_unit, _, hedged in candidate_numbers(reply):
        if present(tok, b_raw, b_f, hedged):
            if has_unit and LEVEL == "strict" and not present(tok, ban_raw, ban_f, hedged):
                unbannered.append(tok)
            continue
        if present(tok, ag_raw, ag_f, hedged) or present(tok, ot_raw, ot_f, hedged):
            if not present(tok, unver_raw, unver_f, hedged):
                subagent.append(tok)
            continue
        typed.append(tok)

    if typed:
        v.append("N1 TYPED, NOT PASTED -- these numbers appear nowhere in this session's "
                 f"command output: {', '.join(sorted(set(typed))[:12])}\n"
                 "    Compute them and paste the output, or delete the claim. Do not "
                 "retype them from memory or from a figure.")
    if subagent:
        v.append("N2 SUBAGENT NUMBER, NOT RE-DERIVED -- present only in a subagent report "
                 f"or a file read, never computed here: {', '.join(sorted(set(subagent))[:12])}\n"
                 "    Either re-run the computation in this session, or list each on a line "
                 "beginning `UNVERIFIED:` saying where it came from and what statistic it is.")
    if unbannered:
        v.append("N3 UNLABELLED QUANTITY -- these carry physical units but came from a "
                 f"command that printed no provenance banner: {', '.join(sorted(set(unbannered))[:12])}\n"
                 "    Re-run through trust/provenance.py so the inputs, the selection and "
                 "the column definition are printed with the number.")

    for path, led in session_ledgers(bash, cwd):
        for name in led.get("unasked", []):
            if not re.search(rf"\b{re.escape(name)}\b", reply):
                v.append(f"P  UNDISCLOSED INVENTED PARAMETER -- run {os.path.basename(path)} "
                         f"used `{name}`, which you chose unasked, and the reply does not "
                         f"mention it. State it and what it excluded, or drop the parameter.")
        for i in led.get("inputs", []):
            try:
                sys.path.insert(0, os.path.dirname(HERE))
                from provenance import _digest
                if os.path.exists(i["path"]) and _digest(i["path"]) != i["digest"]:
                    v.append(f"S  STALE INPUT -- {i['path']} has changed since the run that "
                             f"produced these numbers. Re-run before reporting.")
            except Exception:
                pass

    fused = [s for s in sentences(FENCE_RE.sub(" ", reply))
             if INTERP_RE.search(s) and candidate_numbers(s)] if LEVEL != "core" else []
    if fused and not INTERP_LABEL.search(reply):
        v.append("I  INTERPRETATION FUSED TO MEASUREMENT -- the reply asserts a cause or a "
                 "consequence in the same breath as a number, unlabelled. Put the number "
                 "first, then a separate line starting `Interpretation:` -- so it can be "
                 "skipped, or shot down, without re-auditing the measurement.\n"
                 f"    Offending sentence: {fused[0].strip()[:180]}")

    tabs = markdown_tables(reply) if LEVEL != "core" else []
    wide = [c for _, c in tabs if c > MAX_COLS]
    if len(tabs) > 1 or wide:
        v.append(f"V  VOLUME -- {len(tabs)} table(s), widest {max([c for _,c in tabs], default=0)} "
                 f"columns (limit {MAX_COLS}, one table). Extra columns are audit surface "
                 "with no requested value, and they bury the column that was asked for. "
                 "Put the full table in the answer file; quote the one column here.")
    n_lines = len([l for l in reply.splitlines() if l.strip()])
    if n_lines > MAX_LINES and LEVEL != "core":
        v.append(f"V  VOLUME -- reply is {n_lines} non-blank lines (limit {MAX_LINES}). "
                 "Move the detail to the answer file and link it.")
    return v


def ledger_footer(results, cwd):
    leds = session_ledgers("\n".join(results.get("Bash", [])), cwd)
    if not leds:
        return ("\nNo provenance-carrying run was recorded in this session. If you are "
                "reporting a measurement, it did not come from one.")
    out = ["\nGROUND TRUTH -- what the runs in this session actually read:"]
    for path, led in leds:
        out.append(f"  run: {led.get('question','?')[:90]}")
        for i in led.get("inputs", []):
            out.append(f"    {i['path']}")
            out.append(f"        IS: {i['role']}")
        for m in led.get("masks", []):
            out.append(f"    selection {m['name']}: keeps {m['kept']:,}/{m['total']:,} "
                       f"({100*m['frac']:.1f}%) -- {m['defn']}")
    out.append("  If the reply names a dataset or a population that is not on this list, "
               "the reply is wrong about its own inputs.")
    return "\n".join(out)


def main():
    if os.environ.get("TRUST_OFF") == "1" or LEVEL == "off":
        return 0
    try:
        data = json.load(sys.stdin)
    except Exception:
        return 0                       # never wedge the session on a malformed payload
    if data.get("stop_hook_active"):
        return 0                       # loop guard: let the human see it and decide
    tp = data.get("transcript_path")
    cwd = data.get("cwd") or os.getcwd()
    if not tp or not os.path.exists(tp):
        return 0
    reply, results = read_transcript(tp)
    if not reply.strip():
        return 0
    v = check(reply, results, cwd)
    if not v:
        return 0
    sys.stderr.write(
        "TRUST GATE -- this reply cannot be delivered as written.\n"
        "Each item below is a cost you are about to transfer to Andy as audit time.\n"
        "Fix it here, in this turn; do not explain it away.\n\n"
        + "\n\n".join(v) + "\n"
        + ledger_footer(results, cwd) + "\n")
    return 2


# ----------------------------------------------------------------------- tests
def _selftest():
    import tempfile
    ok = True

    def mk(reply, bash="", agent=""):
        recs = []
        if bash:
            recs += [{"type": "assistant", "message": {"content": [{"type": "tool_use", "id": "t1", "name": "Bash"}]}},
                     {"type": "user", "message": {"content": [{"type": "tool_result", "tool_use_id": "t1", "content": bash}]}}]
        if agent:
            recs += [{"type": "assistant", "message": {"content": [{"type": "tool_use", "id": "t2", "name": "Agent"}]}},
                     {"type": "user", "message": {"content": [{"type": "tool_result", "tool_use_id": "t2", "content": agent}]}}]
        recs.append({"type": "assistant", "message": {"content": [{"type": "text", "text": reply}]}})
        f = tempfile.NamedTemporaryFile("w", suffix=".jsonl", delete=False)
        for r in recs:
            f.write(json.dumps(r) + "\n")
        f.close()
        return read_transcript(f.name)

    cases = [
        ("typed number blocked", "The RMS per-block SE is 22.0 mm.", "", "", "N1"),
        ("pasted number allowed", "The RMS per-block SE is 22.0 mm.",
         BANNER_START + "\n  SE_rms_mm  22.0\n", "", None),
        ("subagent number blocked", "The SE is 10.0 mm.", "", "SE = 10.0 mm\n", "N2"),
        ("subagent number tagged ok", "The SE is 10.0 mm.\nUNVERIFIED: 10.0 mm is the "
         "subagent's figure, not re-derived and not the RMS per-block SE.", "", "SE = 10.0 mm\n", None),
        ("unbannered unit blocked (strict only)", "Offset is 22.0 mm.", "22.0\n", "", "N3!strict"),
        ("unbannered unit allowed at normal", "Offset is 22.0 mm.", "22.0\n", "", None),
        ("hedged sig-fig rounding allowed", "Erosion is ~100,000 m3.",
         BANNER_START + "\n  erosion_m3  100234.7\n", "", None),
        ("unhedged sig-fig rounding blocked", "Erosion is 100,000 m3.",
         BANNER_START + "\n  erosion_m3  100234.7\n", "", "N1"),
        ("interp fused blocked", "The offset is 22.0 mm, which means canopy drives it.",
         BANNER_START + "\n22.0\n", "", "I"),
        ("interp labelled ok", "The offset is 22.0 mm.\nInterpretation: this means canopy drives it.",
         BANNER_START + "\n22.0\n", "", None),
        ("wide table blocked", "| a | b | c | d | e | f |\n|---|---|---|---|---|---|\n| 1 | 2 | 3 | 4 | 5 | 6 |",
         BANNER_START + "\n1 2 3 4 5 6\n", "", "V"),
        ("prose no numbers ok", "Done. The script is in analysis/ridgelines/.", "", "", None),
    ]
    global LEVEL
    for name, reply, bash, agent, want in cases:
        lvl = "normal"
        if want and want.endswith("!strict"):
            want, lvl = want.split("!")[0], "strict"
        r, res = mk(reply, bash, agent)
        LEVEL = lvl
        got = check(r, res, os.getcwd())
        LEVEL = "normal"
        codes = {x.split()[0] for x in got}
        if want is None:
            good = not got
        else:
            good = want in codes
        print(("OK  " if good else "FAIL") + f" {name}" + ("" if good else f"  -> {codes or 'no violations'}"))
        ok &= good
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else main())
