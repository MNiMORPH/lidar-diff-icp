"""FORWARD model of slope-induced pulse broadening (not a correlation).
Geometry: a laser footprint of diameter D hits a planar slope S. The downslope extent
of the footprint spans an ELEVATION interval, so the returned pulse is spread in
elevation. For a Gaussian beam (1/e^2 diameter D) the elevation distribution of the
returned energy is Gaussian with
        sigma_z(S) = (D/2) * tan(S)                      [aspect-INDEPENDENT]
The recorded bare-earth elevation is a biased draw from that spread: bias = -k*sigma_z,
where k is the detector/return offset (k=0 centroid/unbiased; k>0 if the ground echo
sits below centre -- e.g. last-return/far-edge or range-walk on the broadened pulse).
So per epoch:  bias(S) = -k * (D/2) * tan(S).
DoD (gen2-gen1) = bias2 - bias1 = k/2 * (D1 - D2) * tan(S)  -> gen1 LOW if D1 > D2.
Prediction is LINEAR IN tan(S) and aspect-independent -- test that against the data,
then see what (D1-D2, k) magnitude is required."""
import numpy as np
from scipy.ndimage import distance_transform_edt as edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# --- measured signal: median DoD vs tan(slope) (the aspect-independent term) ---
res=5.0; dod=np.load("data/derived/elba/dod_refdatum.npy"); Z=np.load("data/derived/elba/z_after.npy")
Zf=Z.copy(); nm=~np.isfinite(Zf)
if nm.any(): Zf=Zf[tuple(edt(nm,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,res); slope=np.arctan(np.hypot(gx,gy))
fin=np.isfinite(dod)&np.isfinite(Z)&(np.degrees(slope)>=3)&(np.degrees(slope)<=45)
tanS=np.tan(slope[fin]); d=dod[fin]
# median DoD in tan(S) bins
edges=np.linspace(np.tan(np.radians(3)),np.tan(np.radians(45)),12); ctr=(edges[:-1]+edges[1:])/2
med=np.array([np.median(d[(tanS>=edges[i])&(tanS<edges[i+1])]) for i in range(len(edges)-1)])
# fit DoD = c * tan(S) through origin (robust-ish: fit to the bin medians)
c=np.sum(ctr*med)/np.sum(ctr*ctr)
print(f"MEASURED: DoD = c * tan(S),  c = {1000*c:+.1f} mm  (aspect-independent term)")
print(f"  linearity R^2 (medians vs tan S through origin): {1 - np.sum((med-c*ctr)**2)/np.sum((med-med.mean())**2):.3f}")

# --- what footprint difference / detector offset reproduces c = k/2 * (D1-D2) ? ---
print(f"\nFORWARD constraint:  c = (k/2)*(D1 - D2)  ->  (D1 - D2) = 2c/k = {2*c:.3f}/k m")
print(f"  detector offset k     required (D1-D2):")
for k in (0.5, 0.8, 1.0, 1.28):   # 1.28 ~ if ground echo sits at ~p10 of the spread
    print(f"    k={k:.2f} (ground echo ~{k:.1f} sigma below centre): D1-D2 = {2*c/k*100:+.1f} cm")
print("  (positive D1-D2 = 2008 footprint LARGER than 2021 -> gen1 lower, matching sign)")

# --- illustrate the broadening magnitude for example footprints ---
print("\nElevation spread sigma_z = (D/2) tan(S) [cm], example footprints:")
for D in (0.4,0.6,0.8):
    print(f"  D={D:.1f} m:  " + "  ".join(f"S={s}d:{100*(D/2)*np.tan(np.radians(s)):.1f}" for s in (10,20,30,40)))

# --- figure: measured median DoD vs tan(S) with the linear broadening model ---
fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].plot(ctr,1000*med,"o-",label="measured median DoD")
ax[0].plot(ctr,1000*c*ctr,"r--",label=f"broadening model  c·tan(S), c={1000*c:.0f} mm")
ax[0].set_xlabel("tan(slope)"); ax[0].set_ylabel("DoD (mm)")
ax[0].set_title("Aspect-independent bias vs tan(slope):\nlinear = the broadening prediction"); ax[0].legend(); ax[0].grid(alpha=.3)
Sd=np.linspace(0,45,100)
for D in (0.4,0.6,0.8): ax[1].plot(Sd,100*(D/2)*np.tan(np.radians(Sd)),label=f"footprint D={D} m")
ax[1].set_xlabel("slope (deg)"); ax[1].set_ylabel("elevation spread σ_z (cm)")
ax[1].set_title("Slope-induced return spread σ_z=(D/2)tan(S)"); ax[1].legend(); ax[1].grid(alpha=.3)
fig.savefig("figures/broadening_forward.png",dpi=130,bbox_inches="tight")
print("\nwrote figures/broadening_forward.png")
