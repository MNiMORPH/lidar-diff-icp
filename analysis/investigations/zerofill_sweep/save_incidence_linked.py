"""Save incidence angle LINKED TO THE LASER PULSES. One read of the raw gen1 tile ->
(1) a LAS copy with incidence/slope/dpsi attached as extra per-point dimensions
    (physically travels with each return/pulse; usable in any lidar tool), and
(2) a portable pulse-keyed table (gps_time + point_source_id + return_number identify
    the pulse/return) with every geometry field, so anything can be computed later.
Incidence recomputed in the SAME pass as the read -> no cross-file alignment risk."""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt as edt
X0,Y0,X1,Y1=577492.8,4882737.6,580032.8,4886237.6; res=5.0
f=laspy.read("data/before/4342-29-64.laz")
x=np.asarray(f.x); y=np.asarray(f.y)
scan_deg=np.asarray(f.scan_angle_rank).astype(np.float64)   # this tile: rank, in degrees
th=np.radians(np.abs(scan_deg)); psid=np.asarray(f.point_source_id)
# DEM slope + downslope azimuth
Z=np.load("data/derived/elba/z_after.npy"); Zf=Z.copy(); nm=~np.isfinite(Zf)
if nm.any(): Zf=Zf[tuple(edt(nm,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,res); S_grid=np.arctan(np.hypot(gx,gy)); psi_down=np.arctan2(-gx,-gy)
ny,nx=Z.shape; inb=(x>=X0)&(x<X1)&(y>=Y0)&(y<Y1)
jj=np.clip(((x-X0)/res).astype(int),0,nx-1); ii=np.clip(((y-Y0)/res).astype(int),0,ny-1)
S=np.where(inb,S_grid[ii,jj],np.nan); psi_up=np.where(inb,psi_down[ii,jj]+np.pi,np.nan)
# scan azimuth per return from flight-line cross-track geometry
psi_scan=np.full(len(x),np.nan)
for p in np.unique(psid):
    m=psid==p; XY=np.c_[x[m]-x[m].mean(), y[m]-y[m].mean()]
    _,V=np.linalg.eigh(XY.T@XY); a=V[:,0]; c=XY@a
    psi_scan[m]=np.arctan2(a[0]*np.sign(c), a[1]*np.sign(c))
dpsi=psi_scan-psi_up
i_deg=np.degrees(np.arccos(np.clip(np.cos(th)*np.cos(S)+np.sin(th)*np.sin(S)*np.cos(dpsi),-1,1))).astype(np.float32)
slope_deg=np.degrees(S).astype(np.float32); dpsi_deg=np.degrees(dpsi).astype(np.float32)

# (1) LAS with extra per-point dims -> physically linked to each return/pulse
for nm_,arr in [("incidence_deg",i_deg),("slope_deg",slope_deg),("dpsi_deg",dpsi_deg)]:
    f.add_extra_dim(laspy.ExtraBytesParams(name=nm_,type=np.float32))
    setattr(f,nm_,np.nan_to_num(arr,nan=np.float32(-9999.0)))
f.write("data/las_local/gen1_geom.las")

# (2) portable pulse-keyed table (gps_time+psid+return_number = pulse/return id)
np.savez_compressed("data/derived/elba/gen1_pulse_geometry.npz",
    gps_time=np.asarray(f.gps_time), point_source_id=psid,
    return_number=np.asarray(f.return_number).astype(np.uint8),
    number_of_returns=np.asarray(f.number_of_returns).astype(np.uint8),
    x=x, y=y, z=np.asarray(f.z).astype(np.float32),
    scan_angle_deg=scan_deg.astype(np.float32), incidence_deg=i_deg,
    slope_deg=slope_deg, dpsi_deg=dpsi_deg,
    intensity=np.asarray(f.intensity), classification=np.asarray(f.classification).astype(np.uint8),
    in_bounds=inb,
    key=np.array(["pulse/return id = (gps_time, point_source_id, return_number); "
                  "incidence_deg = angle between laser ray and surface normal; -9999=out-of-DEM in LAS"]))
print(f"wrote data/las_local/gen1_geom.las (+incidence_deg,slope_deg,dpsi_deg extra dims) and "
      f"data/derived/elba/gen1_pulse_geometry.npz for {len(x):,} returns")
print(f"incidence_deg: median {np.nanmedian(i_deg[inb]):.1f}  max {np.nanmax(i_deg[inb]):.1f}")
