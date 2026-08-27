"""Does the control-derived datum CLOSE with the measured DoD?

The check that withdrew the 2026-08-27 retraction. It exists because I once got the
relation wrong and concluded a correct result was contradicted.

    DoD = c1_ours - c2 - g          NOT   c1_ours - c2

NAVD88 is the DATUM; GEOID03 and GEOID18 are MODELS for converting GPS ellipsoidal
heights to orthometric. Both control sets publish NAVD88 and are directly comparable --
but each epoch's LIDAR z was converted with a different geoid model, so c1 and c2
reference surfaces in different frames. The geoid does not cancel between them; it is
exactly the term the pipeline adds to gen1.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--answer", required=True)
    p.add_argument("--gen2-product", required=True)
    p.add_argument("--geoid-mm", type=float, required=True)
    p.add_argument("--measured-dod-mm", type=float, required=True)
    a = p.parse_args(argv)
    ans = json.loads(Path(a.answer).read_text())
    g2 = json.loads(Path(a.gen2_product).read_text())

    R = Run("Does the control-derived gen1 datum close with the measured DoD on stable "
            "open ground?")
    R.input(a.answer, role="adopted gen1 surface offsets at this site")
    R.input(a.gen2_product, role="gen2's datum on its DELIVERED surface, open ground")
    R.param("geoid_mm", a.geoid_mm, src="repo",
            why="references.geoid_difference over elbaext; the GEOID03->GEOID18 term the "
                "pipeline ADDS to gen1, which is why it does not cancel between c1 and c2")
    R.param("measured_dod_mm", a.measured_dod_mm, src="repo",
            why="median of dod_geoid.npy on stable AND open cells, 116,507 cells")
    R.column("quantity", "what is being reported")
    R.column("value_mm", "millimetres; positive = the surface reads LOW, so ADD it")
    R.banner()

    c1 = float(ans["our_surface_mm"]); c2 = float(g2["constant_mm"])
    pred = c1 - c2 - a.geoid_mm
    c1_pipe = c1 - a.geoid_mm
    corr = c2 - c1_pipe
    rows = [
        ["c1_ours, gen1 on OUR surface", f"{c1:+.2f}"],
        ["c2, gen2 on its DELIVERED surface", f"{c2:+.2f}"],
        ["g, geoid term the pipeline adds to gen1", f"{a.geoid_mm:+.2f}"],
        ["PREDICTED DoD = c1_ours - c2 - g", f"{pred:+.2f}"],
        ["MEASURED DoD, stable AND open", f"{a.measured_dod_mm:+.2f}"],
        ["MISS", f"{abs(pred - a.measured_dod_mm):.4f}"],
        ["c1 as it sits IN the DoD (c1_ours - g)", f"{c1_pipe:+.2f}"],
        ["DoD absolute correction (c2 - c1_pipe)", f"{corr:+.2f}"],
        ["corrected DoD on stable open ground", f"{a.measured_dod_mm + corr:+.2f}"],
    ]
    R.table(["quantity", "value_mm"], rows)
    R.done(headline=f"closure miss {abs(pred-a.measured_dod_mm):.4f} mm; DoD absolute "
                    f"correction {corr:+.2f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
