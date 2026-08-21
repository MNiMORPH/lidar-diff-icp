"""Forward waveform model of slope-induced broadening, standard settings (tentative).
Gaussian beam (1/e diameter D) on a planar slope S. The footprint's downslope extent
maps to a RANGE/elevation interval, so the geometric return distribution is Gaussian
with sigma_geom = sigma_x * tan(S), sigma_x = (D/2)/sqrt(2). The received WAVEFORM is
that convolved with the transmitted pulse (Gaussian, from FWHM_ns): sigma_tot =
sqrt(sigma_geom^2 + sigma_pulse^2).

SOLID part: the spread sigma_geom (and sigma_tot) vs slope, for the real footprints.
CRUX: spread -> BIAS requires the detector to sample OFF the waveform centroid. A
symmetric waveform + peak/centroid detector gives ZERO bias. Any bias comes from the
detector sampling below centre by some offset k (in units of sigma), e.g. a trailing/
last-return convention (ground = far/low edge). Then bias = -k*sigma_geom per epoch,
DoD = (k/2)(sigma_geom1 - sigma_geom2) ... i.e. bias ∝ (D1-D2) tan(S). We SOLVE for
the k that the data requires and ask whether it is physically reasonable -- we do NOT
guess a detector and manufacture a bias."""
import numpy as np
c=299792458.0
D1,D2=0.60,0.35                       # 2008 Gemini-narrow, 2021 TerrainMapper (1/e diam, m)
sx1,sx2=(D1/2)/np.sqrt(2),(D2/2)/np.sqrt(2)
for FWHM_ns in (3.0,5.0,8.0):
    sp=(c*FWHM_ns*1e-9/2)/2.355       # pulse sigma in range/elevation (m)
    print(f"\n--- transmitted pulse FWHM {FWHM_ns} ns -> sigma_pulse {100*sp:.0f} cm (range) ---")
    print(f"{'slope':>6} {'sig_geom 2008':>13} {'sig_geom 2021':>13} {'d(sig_geom)':>11}  (cm)")
    for Sd in (10,20,30,40):
        S=np.radians(Sd); sg1=sx1*np.tan(S); sg2=sx2*np.tan(S)
        print(f"{Sd:>5}d {100*sg1:>12.1f} {100*sg2:>12.1f} {100*(sg1-sg2):>10.1f}")
# the measured aspect-independent bias: DoD = c_meas * tan(S), c_meas = 35 mm
c_meas=0.035
# DoD_model = k*(sx1 - sx2)*tan(S)  ->  k*(sx1-sx2) = c_meas
dsx=sx1-sx2
print(f"\nsigma_x: 2008={100*sx1:.1f} cm  2021={100*sx2:.1f} cm  difference={100*dsx:.1f} cm")
print(f"REQUIRED detector offset k = c_meas/(sx1-sx2) = {c_meas/dsx:.2f}  sigma")
print(f"  interpretation: the ground return sits ~{c_meas/dsx:.2f}*sigma_geom below the footprint centroid.")
print(f"  reference points: centroid/peak detector -> k=0 (NO bias); a far/low-edge or")
print(f"  last-return convention -> k ~ 0.5-1; range-walk on a broadened pulse -> k small-positive.")
print(f"  k={c_meas/dsx:.2f} is modest and physically plausible IF the detector samples below")
print(f"  centroid; it is NOT determined by geometry alone -- it is set by the DISCRIMINATOR.")
# also: single-return fraction matters -- a single bare return is detected at peak (k~0);
# only multi-return / edge-triggered pulses carry the far-edge bias. gen1 is ~87% single.
print(f"\nCAVEAT: gen1 is ~87% single-return; a single bare return is peak-detected (k~0),")
print(f"  so the effective k is diluted -- pushing the required per-pulse offset HIGHER,")
print(f"  or implying the bias is carried by the detector's range-walk on ALL returns, not")
print(f"  by a last-return edge. Distinguishing these needs the sensor's detection method.")
