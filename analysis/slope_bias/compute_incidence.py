"""Compute the INCIDENCE ANGLE i (ray-to-surface, = angle between the laser ray and
the local surface normal) for every gen1 return:
    cos i = cos θ cos S + sin θ sin S cos Δψ
θ = |scan angle| (off-nadir). S, downslope-azimuth = from the DEM (z_after). Δψ =
scan azimuth − upslope azimuth. Scan azimuth comes from each flight line's geometry:
per-line PCA gives along/across-track; the point's cross-track SIDE gives the ray's
horizontal direction (no reliance on the scan-angle-sign convention). Sanity: at
θ≈0, i must equal S. Saves incidence per return aligned to gen1_returns.npz."""
import numpy as np
from scipy.ndimage import distance_transform_edt as edt
X0,Y0,X1,Y1=577492.8,4882737.6,580032.8,4886237.6; res=5.0
d=np.load("data/derived/elba/gen1_returns.npz")
x=d["x"]; y=d["y"]; th=np.radians(np.abs(d["scan_angle_deg"].astype(np.float64))); psid=d["point_source_id"]
# DEM slope + downslope azimuth (fill holes for continuity)
Z=np.load("data/derived/elba/z_after.npy"); Zf=Z.copy(); nm=~np.isfinite(Zf)
if nm.any(): Zf=Zf[tuple(edt(nm,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,res)                       # gx=dz/dEast(col), gy=dz/dNorth(row)
S_grid=np.arctan(np.hypot(gx,gy))               # slope (rad)
psi_down=np.arctan2(-gx,-gy)                     # azimuth of steepest DESCENT (from N, cw; x=E,y=N)
ny,nx=Z.shape
inb=(x>=X0)&(x<X1)&(y>=Y0)&(y<Y1)
jj=np.clip(((x-X0)/res).astype(int),0,nx-1); ii=np.clip(((y-Y0)/res).astype(int),0,ny-1)
S=np.where(inb,S_grid[ii,jj],np.nan); psi_up=np.where(inb,psi_down[ii,jj]+np.pi,np.nan)
# scan azimuth per return from flight-line cross-track geometry
psi_scan=np.full(len(x),np.nan)
for p in np.unique(psid):
    m=psid==p; X=np.c_[x[m]-x[m].mean(), y[m]-y[m].mean()]
    # PCA: along-track = top eigenvector; across-track = perpendicular
    w,V=np.linalg.eigh(X.T@X); a=V[:,0]          # smallest-variance axis = across-track
    c=X@a                                        # cross-track coordinate (signed)
    ae=a[0]*np.sign(c); an=a[1]*np.sign(c)       # across-track vector toward the point's side
    psi_scan[m]=np.arctan2(ae,an)                # azimuth (from N, cw)
dpsi=psi_scan-psi_up
cosi=np.cos(th)*np.cos(S)+np.sin(th)*np.sin(S)*np.cos(dpsi)
i_deg=np.degrees(np.arccos(np.clip(cosi,-1,1)))
np.savez_compressed("data/derived/elba/gen1_incidence.npz",
    incidence_deg=i_deg.astype(np.float32), scan_deg=np.degrees(th).astype(np.float32),
    slope_deg=np.degrees(S).astype(np.float32), dpsi_deg=np.degrees(dpsi).astype(np.float32),
    in_bounds=inb)
# sanity + summary
ok=inb&np.isfinite(i_deg)
near=ok&(np.degrees(th)<2)
print(f"incidence computed for {ok.sum():,} in-bounds returns (of {len(x):,})")
print(f"SANITY at |scan|<2deg: median(i - S) = {np.median(i_deg[near]-np.degrees(S[near])):+.2f} deg (expect ~0)")
print(f"incidence: median {np.median(i_deg[ok]):.1f}  p90 {np.percentile(i_deg[ok],90):.1f}  max {i_deg[ok].max():.1f} deg")
sd=np.degrees(S[ok])
for lo,hi in [(0,5),(5,15),(15,25),(25,90)]:
    mm=ok&(np.degrees(S)>=lo)&(np.degrees(S)<hi)
    print(f"  slope {lo:2d}-{hi:2d}: median incidence {np.median(i_deg[mm]):.1f} deg  (scan spreads it +/- from slope)")
print("saved -> data/derived/elba/gen1_incidence.npz")
