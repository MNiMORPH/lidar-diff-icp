"""Datum with cover as a COVARIATE and per-flight-line offsets kept, read out at OPEN.

The two previous options each threw something away:

  * open-only  -- respects the datum's meaning but has 2 marks inside 10 km and must
                  extrapolate the site's level from control 14-63 km away;
  * all covers pooled -- reaches 10 marks inside 10 km and all 6 lines, but its scatter
                  (sd 74.0 mm) is mostly COVER, so its mean measures the canopy mix of
                  whoever sited the control, not the datum.

This uses every mark for SPATIAL and LINE coverage while still reporting the OPEN-ground
level, by fitting

    tie_ij = mu + L_i + C_j        L = flight line, C = cover class, C_open == 0

so ``mu + L_i`` IS line i's open-ground level, and the site's datum is the
ground-weighted combination of those over the lines that actually build its tile.

Cover enters as a nuisance term to be estimated and removed, never as a filter and never
pooled into the level. The line offsets are free parameters, so the between-line structure
(ANOVA F = 8.63, p < 0.001) is respected rather than averaged over.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", required=True)
    p.add_argument("--weights", required=True)
    p.add_argument("--reference-cover", required=True)
    p.add_argument("--max-km", type=float, default=None)
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    d = json.loads(Path(a.product).read_text())
    w = {int(k): float(v) for k, v in json.loads(Path(a.weights).read_text()).items()}
    marks = d["marks"]

    R = Run("What is the site's OPEN-ground datum if cover is a covariate to be removed "
            "and every flight line keeps its own offset?")
    R.input(a.product, role="per-mark ties with the line their returns assign and their "
                            "MnDNR cover class")
    R.input(a.weights, role="each line's share of the ground cells the pipeline built "
                            "for this tile")
    R.param("model", "tie ~ line + cover, additive, no interaction", src="MINE",
            why="an interaction would need a mark in every line x cover cell; there are "
                "6 lines and 6 covers against 25 marks, so it is not estimable")
    R.param("reference_cover", a.reference_cover, src="andy",
            why="set to zero, so the intercept plus each line effect IS that line's "
                "OPEN-ground level and the datum is read out there")
    R.param("max_km", a.max_km, src="MINE",
            why="optional restriction to the near field; None uses every mark")
    R.column("term", "model coefficient: the intercept (reference line at reference cover), a line offset, or a cover effect relative to the reference cover; all in mm")
    R.column("estimate_mm", "its fitted value, mm")
    R.column("se_mm", "its standard error, mm")
    R.column("line", "flight line point_source_id")
    R.column("w", "that line's share of the tile's ground cells, 0-1")
    R.column("open_level_mm", "mu + L_i, that line's OPEN-ground datum, mm")
    R.column("stat", "name of the statistic")
    R.column("value", "its value; mm unless stated")
    R.banner()

    if a.max_km is not None:
        marks = [m for m in marks if m.get("km", 0) <= a.max_km]
    y = np.array([m["tie_mm"] for m in marks], float)
    lines = sorted({int(m["line"]) for m in marks})
    covers = sorted({m["cover"] for m in marks})
    if a.reference_cover not in covers:
        raise SystemExit(f"reference cover {a.reference_cover} not present in {covers}")
    lref, cref = lines[0], a.reference_cover
    lcols = [l for l in lines if l != lref]
    ccols = [c for c in covers if c != cref]
    X = np.column_stack(
        [np.ones(len(marks))]
        + [[1.0 if int(m["line"]) == l else 0.0 for m in marks] for l in lcols]
        + [[1.0 if m["cover"] == c else 0.0 for m in marks] for c in ccols])
    names = ["mu(line %d, %s)" % (lref, cref)] + [f"line {l}" for l in lcols] \
        + [f"cover {c}" for c in ccols]
    rank = int(np.linalg.matrix_rank(X))
    print(f"  design: {X.shape[0]} marks, {X.shape[1]} parameters, rank {rank}, "
          f"dof {X.shape[0]-rank}")
    if rank < X.shape[1]:
        print("  ** RANK DEFICIENT: some line x cover combination is not identifiable **")
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = X.shape[0] - rank
    s2 = float((resid ** 2).sum() / dof)
    XtXi = np.linalg.pinv(X.T @ X)
    cov = s2 * XtXi
    se = np.sqrt(np.diag(cov))
    R.table(["term", "estimate_mm", "se_mm"],
            [[n, f"{b:+.2f}", f"{s:.2f}"] for n, b, s in zip(names, beta, se)])

    rows, cvec = [], np.zeros(X.shape[1])
    wt = sum(w[l] for l in lines)
    for l in lines:
        c = np.zeros(X.shape[1]); c[0] = 1.0
        if l != lref:
            c[1 + lcols.index(l)] = 1.0
        lvl = float(c @ beta)
        rows.append([l, f"{w[l]/wt:.4f}", f"{lvl:+.2f}"])
        cvec += (w[l] / wt) * c
    R.table(["line", "w", "open_level_mm"], rows)

    val = float(cvec @ beta)
    sd = float(np.sqrt(cvec @ cov @ cvec))
    out = [["marks / lines / covers", f"{len(marks)} / {len(lines)} / {len(covers)}"],
           ["residual sd of the fit", f"{np.sqrt(s2):.2f}"],
           ["DATUM at open, ground-weighted over lines", f"{val:+.2f} +/- {sd:.2f}"]]
    R.table(["stat", "value"], out)
    if a.out:
        Path(a.out).write_text(json.dumps(dict(
            datum_open_mm=val, se_mm=sd, n_marks=len(marks), n_lines=len(lines),
            resid_sd_mm=float(np.sqrt(s2)),
            terms={n: [float(b), float(s)] for n, b, s in zip(names, beta, se)}), indent=1))
        print(f"\n  wrote {a.out}")
    R.done(headline=f"open-ground datum {val:+.2f} +/- {sd:.2f} mm from "
                    f"{len(marks)} marks, cover removed as a covariate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
