#!/usr/bin/env python3
"""Does the BRIDGE organise by flight line, or scatter?

From the DeLong 2022 re-read (2026-09-11). DeLong could only INFER that the vendor had
applied a ground-control-based adjustment: their reprocessing from raw PO data removed it,
and they saw what was left as flightline-parallel streaks in a DoD, which they patched with
a stable-terrain correction surface of mean 0.20 m. They never measured the lost adjustment.

`bridge = z_delivered - z_ours` measures that same quantity DIRECTLY, per mark.

HYPOTHESIS UNDER TEST: if the bridge is reading the vendor's PER-LINE control tie competing
with our free-network swath re-solve, the bridge should organise BY FLIGHT LINE -- marks on
the same line agreeing with each other more than marks on different lines. If instead it is
our reconstruction's own local noise, it should scatter within a line just as much as
between lines.

DO NOT OVERCLAIM THE PARALLEL. We re-solve on the DELIVERED cloud, not raw PO, so we do not
undo the vendor's tie the way DeLong did. A free-network solve does redistribute level
BETWEEN lines, which is what could compete with a per-line control adjustment -- that is the
mechanism being tested, not an established one.

WHY THIS IS A PEEK AND NOT THE ANSWER. It runs on the marks that ALREADY carry a per-mark
DOMINANT flight line from an earlier run. The bridge records themselves carry only the line
SET present in each window (77 of 88 marks see two lines), which cannot group. So the
population here is not chosen by me and is not random either: it is whatever an earlier
analysis happened to assign. n is small and the answer is indicative only. Getting all 88
needs `same_line.assign_line_from_returns` over the unassigned marks -- an hour-scale run,
and the SAME one that #66's within-line comparison needs.

Reports one-way ANOVA of bridge_mm by line plus ICC(1). No mark is cut and no threshold is
applied; the per-line means and every group size are printed so a result driven by one
group is visible.

    ./lidar-icp/bin/python ground_control/bridge_by_line.py \
        --bridge <cmp_csf.json> --assign ground_control/products/gen1_datum_by_returns_elba.json
"""
import argparse
import json
from collections import defaultdict

import numpy as np

from trust.provenance import Run


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bridge", required=True, help="a run_bridge_wide.py output JSON")
    ap.add_argument("--assign", required=True,
                    help="JSON with marks[].point_id and marks[].line (a DOMINANT line per "
                         "mark, from assign_line_from_returns)")
    ap.add_argument("--permutations", type=int, default=200000,
                    help="permutations for the tightest-group test; 0 to skip")
    ap.add_argument("--seed", type=int, default=20260911,
                    help="fixed so the reported p-value reproduces exactly")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    br = json.load(open(a.bridge))
    asg = json.load(open(a.assign))

    R = Run("Does the bridge (z_delivered - z_ours) organise by FLIGHT LINE, or scatter "
            "within a line as much as between lines?")
    R.input(a.bridge, role="bridge per mark; supplies bridge_mm, the quantity being grouped")
    R.input(a.assign, role="the DOMINANT flight line per mark, from an earlier "
                           "assign_line_from_returns run; supplies the GROUPING only, never "
                           "a value")
    R.param("population", "marks that already carry a dominant line", src="MINE",
            why="NOT a random sample and NOT chosen by value: it is whichever marks an "
                "earlier analysis assigned. The bridge records carry only the line SET in "
                "the window (most marks see two lines), which cannot group. Stated because "
                "it bounds what this can conclude, and it excludes every unassigned mark")
    R.param("grouping", "dominant line, as assigned", src="repo",
            why="taken as given from the assignment run; this script does no assignment")
    R.param("outlier_rule", "NONE APPLIED", src="andy",
            why="no mark cut. Every group size and per-line mean is printed so a result "
                "driven by one group or one mark is visible without cutting it")
    R.column("line", "gen1 point_source_id, the flight line, unitless")
    R.column("n", "marks on that line")
    R.column("mean_mm", "mean bridge on that line, mm")
    R.column("sd_mm", "sd of bridge WITHIN that line, ddof=1, mm; nan for a singleton")
    R.column("F", "one-way ANOVA MS_between / MS_within of bridge_mm by line. Large F = "
                  "lines differ from each other by more than marks differ within a line")
    R.column("tightest_group_perm_p", "P that ANY line looks as internally consistent as "
                                     "the tightest one does, under permutation of the line "
                                     "labels across the same group sizes. Pays for the fact "
                                     "that the tightest group was picked after looking")
    R.column("ICC1", "(MS_b - MS_w) / (MS_b + (k-1) MS_w), k = mean group size. The "
                     "fraction of bridge variance that is BETWEEN lines. Can go negative, "
                     "which means no between-line structure at all; reported as computed")

    line = {m["point_id"]: m["line"] for m in asg["marks"]}
    bridge = {r["point_id"]: r["bridge_mm"] for r in br["marks"]}
    both = sorted(set(line) & set(bridge))

    R.mask("assigned_and_bridged", np.ones(len(both), bool),
           defn=f"marks carrying BOTH a dominant line and a bridge. "
                f"{len(set(bridge) - set(line))} bridged marks have no assignment and are "
                f"NOT counted here; {len(set(line) - set(bridge))} assigned marks have no "
                f"bridge",
           of=len(bridge))

    g = defaultdict(list)
    for p in both:
        g[line[p]].append(bridge[p])

    print(f"\n{'line':>6}{'n':>4}{'mean_mm':>11}{'sd_mm':>10}   marks")
    for ln in sorted(g):
        v = np.array(g[ln], float)
        sd = v.std(ddof=1) if v.size > 1 else float("nan")
        ids = [p for p in both if line[p] == ln]
        print(f"{ln:>6}{v.size:>4}{v.mean():>11.1f}{sd:>10.1f}   {', '.join(ids)}")

    allv = np.array([bridge[p] for p in both], float)
    k_groups = len(g)
    n_tot = allv.size
    grand = allv.mean()
    ss_b = sum(len(v) * (np.mean(v) - grand) ** 2 for v in g.values())
    ss_w = sum(((np.array(v) - np.mean(v)) ** 2).sum() for v in g.values())
    df_b, df_w = k_groups - 1, n_tot - k_groups
    ms_b, ms_w = ss_b / df_b, (ss_w / df_w if df_w > 0 else float("nan"))
    F = ms_b / ms_w if df_w > 0 and ms_w > 0 else float("nan")
    try:
        from scipy.stats import f as fdist
        p = float(fdist.sf(F, df_b, df_w))
    except Exception:
        p = float("nan")
    kbar = n_tot / k_groups
    icc = (ms_b - ms_w) / (ms_b + (kbar - 1) * ms_w) if ms_w == ms_w else float("nan")

    print(f"\n  bridge over all {n_tot} assigned marks: mean {allv.mean():+.1f}  "
          f"sd {allv.std(ddof=1):.1f} mm")
    print(f"  between-line sd of the {k_groups} line means: "
          f"{np.std([np.mean(v) for v in g.values()], ddof=1):.1f} mm")
    print(f"  pooled WITHIN-line sd: {np.sqrt(ms_w):.1f} mm")
    print(f"  one-way ANOVA  F({df_b},{df_w}) = {F:.3f}   p = {p:.4f}")
    print(f"  ICC(1) = {icc:+.3f}  (fraction of bridge variance that is BETWEEN lines)")

    # The ANOVA answers "do the lines differ in LEVEL". It says nothing about a line whose
    # marks are anomalously CONSISTENT -- which is the other way a per-line tie could show.
    # Tested by permuting the line labels across the same group sizes, so the null keeps the
    # observed bridge values and the observed design exactly. "any group" is deliberate: the
    # tightest group is picked AFTER looking, so the test must pay for that choice.
    tight = [(np.std(v, ddof=1), ln) for ln, v in g.items() if len(v) > 1]
    perm_p = float("nan")
    if tight and a.permutations:
        obs = min(t[0] for t in tight)
        sizes = [len(v) for v in g.values()]
        rng = np.random.default_rng(a.seed)
        hits = 0
        for _ in range(a.permutations):
            q = rng.permutation(allv)
            i, best = 0, np.inf
            for sz in sizes:
                if sz > 1:
                    best = min(best, q[i:i + sz].std(ddof=1))
                i += sz
            hits += best <= obs
        perm_p = hits / a.permutations
        ln_t = min(tight)[1]
        print(f"  tightest within-line sd: {obs:.1f} mm (line {ln_t}, n={len(g[ln_t])}); "
              f"P(any group this tight by chance) = {perm_p:.4f} "
              f"[{a.permutations} permutations, seed {a.seed}]")

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(dict(n=n_tot, k_groups=k_groups, F=F, p=p, icc=icc,
                           tightest_group_perm_p=perm_p,
                           ms_between=ms_b, ms_within=ms_w, df=(df_b, df_w),
                           by_line={str(k): v for k, v in g.items()}), fh, indent=1)
            fh.write("\n")
        print(f"    wrote {a.out}")

    R.done(headline=(f"n={n_tot} marks on {k_groups} lines; between-line sd "
                     f"{np.std([np.mean(v) for v in g.values()], ddof=1):.1f} vs "
                     f"within-line {np.sqrt(ms_w):.1f} mm; F({df_b},{df_w})={F:.3f} "
                     f"p={p:.4f}; ICC(1)={icc:+.3f}"))


if __name__ == "__main__":
    main()
