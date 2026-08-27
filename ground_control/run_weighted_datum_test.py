"""Is the datum a WEIGHTED COMBINATION of known line constants, or a mean over lines?

The estimator in use takes the mean of the per-line means and quotes
``sd(line means)/sqrt(n_lines)``. That treats the lines as random draws from a
distribution. They are not: a tile's ground is built from SPECIFIC lines in measured
proportions, so if every line's constant were known exactly the tile's datum would be the
deterministic sum ``sum_i w_i c_i`` with NO contribution from the spread between lines.

If that framing holds, the current SE is conservative and the way to improve it is to
determine each LINE's constant better -- which the statewide per-line correction would do
-- rather than to find more marks at the site.

This tests it. ``w_i`` is measured from the ground cells the pipeline actually built;
``var(c_i)`` from that line's own marks.
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
    p.add_argument("--weights", required=True)
    p.add_argument("--datum-product", required=True)
    p.add_argument("--sigma-mark-mm", type=float, required=True)
    p.add_argument("--se-over-lines-mm", type=float, required=True)
    p.add_argument("--value-over-lines-mm", type=float, required=True)
    a = p.parse_args(argv)

    w = {int(k): float(v) for k, v in json.loads(Path(a.weights).read_text()).items()}
    d = json.loads(Path(a.datum_product).read_text())
    marks = d["marks"]

    R = Run("Is the site's datum a deterministic weighted sum of its lines' constants, "
            "rather than a mean over lines with a between-line SE?")
    R.input(a.weights, role="each flight line's share of the GROUND CELLS the pipeline "
                            "built for this tile -- the weights w_i")
    R.input(a.datum_product, role="open-ground per-line marks and ties at this site")
    R.param("sigma_mark_mm", a.sigma_mark_mm, src="repo",
            why="within-line sd of a mark about its own line mean, measured on 25 open "
                "marks across 18 lines out to 40 km; the per-mark scatter that sets "
                "var(c_i) = sigma_mark^2 / n_i")
    R.param("unknown_line_treatment", "reported BOTH ways", src="MINE",
            why="a line with no mark has an unknown constant; renormalising its weight "
                "away ASSUMES it behaves like the others, so the alternative -- charging "
                "it the full between-line variance -- is reported beside it")
    R.column("line", "flight line point_source_id")
    R.column("w", "that line's share of the tile's ground cells, 0-1")
    R.column("n_marks", "open control marks on that line, count")
    R.column("c_mm", "that line's mean tie, mm; blank if no mark")
    R.column("var_c_mm2", "variance of that line's constant, mm^2 = sigma_mark^2/n")
    R.column("stat", "name of the statistic")
    R.column("value", "its value; mm unless the name says otherwise")
    R.banner()

    by = {}
    for m in marks:
        by.setdefault(int(m["line"]), []).append(m["tie_mm"])
    rows = []
    for ln in sorted(w):
        v = by.get(ln)
        rows.append([ln, f"{w[ln]:.4f}", len(v) if v else 0,
                     f"{np.mean(v):+.2f}" if v else "--",
                     f"{a.sigma_mark_mm**2/len(v):.0f}" if v else "--"])
    R.table(["line", "w", "n_marks", "c_mm", "var_c_mm2"], rows)

    known = {ln: np.mean(by[ln]) for ln in by}
    missing = [ln for ln in w if ln not in known]
    wk = np.array([w[ln] for ln in sorted(known)])
    ck = np.array([known[ln] for ln in sorted(known)])
    nk = np.array([len(by[ln]) for ln in sorted(known)])
    vk = a.sigma_mark_mm ** 2 / nk

    wr = wk / wk.sum()                       # renormalised over known lines
    val_r = float((wr * ck).sum())
    sd_r = float(np.sqrt((wr ** 2 * vk).sum()))

    between_var = float(np.var(ck, ddof=1))  # what an unknown line costs
    wm = sum(w[ln] for ln in missing)
    val_u = float((wk * ck).sum() + wm * ck.mean())
    sd_u = float(np.sqrt((wk ** 2 * vk).sum() + wm ** 2 * (between_var + vk.mean())))

    out = [
        ["lines with a mark (count)", f"{len(known)} of {len(w)}"],
        ["lines with NO mark", f"{missing} (weight {wm:.4f})"],
        ["A. mean over lines (current)", f"{a.value_over_lines_mm:+.2f} +/- {a.se_over_lines_mm:.2f}"],
        ["B. weighted, unknown line renormalised away", f"{val_r:+.2f} +/- {sd_r:.2f}"],
        ["C. weighted, unknown line charged its variance", f"{val_u:+.2f} +/- {sd_u:.2f}"],
        ["between-line sd of the known constants", f"{np.sqrt(between_var):.2f}"],
        ["value shift, B minus A", f"{val_r - a.value_over_lines_mm:+.2f}"],
        ["uncertainty change, B minus A", f"{sd_r - a.se_over_lines_mm:+.2f}"],
    ]
    R.table(["stat", "value"], out)
    print()
    print("  What WOULD tighten it: better per-line constants. sd of the weighted sum")
    print("  as sigma_mark (per-mark scatter about its line mean) is scaled:")
    print("    %-22s %12s"%("sigma_mark_mm","weighted_sd_mm"))
    for fac in (1.0, 0.75, 0.5, 0.25):
        vv = (a.sigma_mark_mm * fac) ** 2 / nk
        print("    %-22.1f %12.2f" % (a.sigma_mark_mm * fac,
                                      float(np.sqrt((wr ** 2 * vv).sum()))))
    print("  (equivalently: more marks per line -- var(c_i) = sigma_mark^2 / n_i)")
    print()
    print(f"  Does the weighted framing tighten it?  "
          f"A {a.se_over_lines_mm:.2f} -> B {sd_r:.2f} -> C {sd_u:.2f} mm")
    R.done(headline=f"weighted {val_r:+.2f} +/- {sd_r:.2f} (renorm) vs mean-over-lines "
                    f"{a.value_over_lines_mm:+.2f} +/- {a.se_over_lines_mm:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
