"""Classify every analysis/ script by EVIDENCE, not by name.

Classes, in precedence order -- the first that matches wins:
  GRAPH      named in a workflow Step's command or `code`
  IMPORTED   another .py imports it as a module
  PRODUCER   writes a file that some other script READS, or that a live doc names
  CITED      named in a .md that is not itself a superseded FRAME
  BROKEN     imports a module that no longer exists
  ORPHAN     none of the above
"""
import ast, glob, os, re, subprocess, sys, json
sys.path.insert(0, "src")
from lidar_diff_icp import workflow as W

ALL = sorted(glob.glob("analysis/**/*.py", recursive=True))
MDS = sorted(glob.glob("**/*.md", recursive=True))

in_graph = set()
for s in W.STEPS:
    in_graph |= {t for t in s.command.split() if t.endswith(".py")}
    in_graph |= {c for c in s.code if c.endswith(".py")}

src = {p: open(p, errors="ignore").read() for p in ALL}
mdtext = {m: open(m, errors="ignore").read() for m in MDS}

# modules that exist to be imported
existing_mods = {os.path.splitext(os.path.basename(p))[0] for p in ALL}
existing_mods |= {os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob("src/lidar_diff_icp/*.py")}

def imports_missing(p):
    bad = []
    for m in re.finditer(r"^\s*from\s+([\w.]+)\s+import|^\s*import\s+([\w.]+)",
                         src[p], re.M):
        mod = (m.group(1) or m.group(2))
        if mod.startswith("lidar_diff_icp."):
            # walk the dotted path: a PACKAGE directory is a module too. Checking only for
            # "<leaf>.py" called 17 working scripts BROKEN, including calibrate_ground_q,
            # which produces the ground-q curve and imports lidar_diff_icp.groundtruth.tie.
            rel = os.path.join("src", *mod.split("."))
            if not (os.path.exists(rel + ".py") or os.path.isdir(rel)):
                bad.append(mod)
    return bad

def imported_by_others(p):
    mod = os.path.splitext(os.path.basename(p))[0]
    hits = [q for q in ALL if q != p and re.search(rf"^\s*(from|import)\s+{mod}\b",
                                                   src[q], re.M)]
    return hits

WRITE = re.compile(r"""np\.savez?\(\s*[^,)]*?["']([^"']+\.(?:npy|npz))["']"""
                   r"""|to_parquet\(\s*[^,)]*?["']([^"']+\.parquet)["']"""
                   r"""|to_csv\(\s*[^,)]*?["']([^"']+\.csv)["']""")
def writes(p):
    return {os.path.basename(g) for m in WRITE.finditer(src[p]) for g in m.groups() if g}

# every artifact any OTHER script reads
reads_any = set()
for p in ALL:
    for m in re.finditer(r"""["']([\w./-]+\.(?:npy|npz|parquet|csv|json))["']""", src[p]):
        reads_any.add(os.path.basename(m.group(1)))
live_needed = {os.path.basename(f) for s in W.STEPS for f in s.requires} | \
              {os.path.basename(f) for s in W.STEPS for f in s.produces}

SUPERSEDED_MD = {m for m in MDS if re.search(r"FRAME_2026-0[89]-0[2467]", m)}
def cited_in(p):
    base = os.path.basename(p)
    return [m for m, t in mdtext.items() if base in t and m not in SUPERSEDED_MD]

rows = []
for p in ALL:
    if p in in_graph:                       cls, why = "GRAPH", "in a Step"
    elif imported_by_others(p):             cls, why = "IMPORTED", ",".join(
                                                os.path.basename(x) for x in imported_by_others(p)[:2])
    elif imports_missing(p):                cls, why = "BROKEN", ",".join(imports_missing(p))
    else:
        w = writes(p)
        useful = (w & live_needed) or (w & reads_any)
        if useful:                          cls, why = "PRODUCER", ",".join(sorted(useful)[:2])
        else:
            c = cited_in(p)
            cls, why = ("CITED", os.path.basename(c[0])) if c else ("ORPHAN", "")
    rows.append((cls, p, why))

order = {"GRAPH":0,"IMPORTED":1,"PRODUCER":2,"CITED":3,"BROKEN":4,"ORPHAN":5}
rows.sort(key=lambda r: (order[r[0]], r[1]))
from collections import Counter
c = Counter(r[0] for r in rows)
print("CLASS COUNTS:", dict(sorted(c.items(), key=lambda kv: order[kv[0]])), f"total {len(rows)}")
json.dump([{"cls":a,"path":b,"why":d} for a,b,d in rows], open(sys.argv[1], "w"), indent=1)
for a,b,d in rows:
    if a in ("BROKEN","ORPHAN"): print(f"  {a:9s} {b}")
