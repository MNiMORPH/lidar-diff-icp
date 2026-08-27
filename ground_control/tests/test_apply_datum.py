"""The gauge-invariance property, demonstrated rather than asserted."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ground_control"))
import apply_datum as A  # noqa: E402

# elbaext's measured per-swath dz, mm (corrections_geoid.json)
DZ = {133: 0.00, 134: 22.00, 135: 6.20, 136: -9.80, 137: -18.40, 138: -22.60}


def test_gauge_choice_moves_an_uncorrected_elevation():
    d = A.DatumApplication(constant_mm=58.70, sigma_mm=25.89, gauge_ref=133,
                           source="ANSWER_gen1_elba.json")
    unc, _ = A.gauge_invariance_residual(223305.0, DZ, d)
    assert np.ptp(unc) > 44.0, "the gauge must matter before correction"
    assert abs(np.ptp(unc) - 44.60) < 0.01


def test_correction_makes_the_elevation_gauge_invariant():
    d = A.DatumApplication(constant_mm=58.70, sigma_mm=25.89, gauge_ref=133,
                           source="ANSWER_gen1_elba.json")
    unc, cor = A.gauge_invariance_residual(223305.0, DZ, d)
    assert np.ptp(cor) < 1e-9, f"corrected spread {np.ptp(cor)} must be zero"
    assert np.ptp(unc) / max(np.ptp(cor), 1e-12) > 1e9


def test_regauging_is_reversible():
    d = A.DatumApplication(58.70, 25.89, 133, "x")
    back = d.regauged_to(138, DZ).regauged_to(133, DZ)
    assert abs(back.constant_mm - d.constant_mm) < 1e-9
