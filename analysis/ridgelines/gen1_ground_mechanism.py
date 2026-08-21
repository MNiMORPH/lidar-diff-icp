#!/usr/bin/env python3
"""Why does gen1's ground peak sink ~20 mm in denser (core) forest? — a gen1-INTERNAL test.
Stream gen1 class-2 ground returns in forest, keep per-return: slope-normal d, intensity,
number_of_returns (nr), return_number (rn), scan angle, and the cell's penetration + core flag.
Test physical mechanisms:
  (1) range-walk/timewalk: weaker (low-intensity) returns trigger later -> lower d?
  (2) multi-return/path delay: ground as a late/last return of a multi-return pulse -> lower d?
and whether intensity / nr shift enough between all-forest and core-forest to make the 20 mm.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_ground_mechanism.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zff=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
pen=np.load("data/derived/elba_fulldensity/penetration.npy"); fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool)
cf=np.load("data/derived/elba_fulldensity/core_forest.npy")
is_forest=((pen<0.25)&~fld&np.isfinite(pen)).ravel(); is_core=cf.ravel(); penf=pen.ravel()

D=[];I=[];NR=[];RN=[];PN=[];CO=[];AN=[]
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        cl=np.asarray(pts.classification); g=cl==2
        x=np.asarray(pts.x,np.float64)[g]; y=np.asarray(pts.y,np.float64)[g]; z=np.asarray(pts.z,np.float64)[g]
        inten=np.asarray(pts.intensity)[g]; nr=np.asarray(pts.number_of_returns)[g]; rn=np.asarray(pts.return_number)[g]
        ang=np.asarray(pts.scan_angle_rank)[g]
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)
        ix=ix[keep];iy=iy[keep];cell=iy*NX+ix
        fmask=is_forest[cell]
        if not fmask.any(): continue
        cell=cell[fmask]; x=x[keep][fmask]; y=y[keep][fmask]; z=z[keep][fmask]
        xc=X0+((cell%NX)+0.5)*RES; yc=Y0+((cell//NX)+0.5)*RES
        d=(z-(Zff[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*cosf[cell]
        D.append(d); I.append(inten[keep][fmask]); NR.append(nr[keep][fmask]); RN.append(rn[keep][fmask])
        AN.append(np.abs(ang[keep][fmask])); PN.append(penf[cell]); CO.append(is_core[cell])
d=np.concatenate(D); inten=np.concatenate(I).astype(float); nr=np.concatenate(NR); rn=np.concatenate(RN)
ang=np.concatenate(AN).astype(float); pn=np.concatenate(PN); co=np.concatenate(CO)
d=d*1000  # mm
print(f"gen1 class-2 forest returns: {d.size:,}")
def med(x): return np.median(x)
print(f"\nall forest median d = {med(d):+.1f} mm ; core forest = {med(d[co]):+.1f} mm  (shift {med(d[co])-med(d):+.1f})")

print("\n(1) RANGE-WALK: median d by intensity bin (low intensity = weak return):")
e=np.quantile(inten,np.linspace(0,1,9))
for i in range(len(e)-1):
    b=(inten>=e[i])&(inten<e[i+1] if i<len(e)-2 else inten<=e[i+1])
    if b.sum()<500: continue
    print(f"  intensity {e[i]:5.0f}-{e[i+1]:5.0f} (med {np.median(inten[b]):4.0f}): d {med(d[b]):+6.1f} mm  n={b.sum()}")
print(f"  corr(d, intensity) = {np.corrcoef(d,inten)[0,1]:+.3f}")

print("\n(2) MULTI-RETURN: median d by number_of_returns, and single vs multi:")
for k in [1,2,3,4]:
    b=nr==k
    if b.sum()<500: continue
    print(f"  nr={k}: d {med(d[b]):+6.1f} mm  n={b.sum()} ({b.mean()*100:.1f}%)")
single=nr==1; multi=nr>1
print(f"  single(nr=1) d {med(d[single]):+.1f}  vs multi(nr>1) d {med(d[multi]):+.1f}  -> multi is {med(d[multi])-med(d[single]):+.1f} mm")
print(f"  is-ground-the-LAST-return? rn==nr: d {med(d[rn==nr]):+.1f}  vs not: {med(d[rn!=nr]):+.1f}")

print("\n(A) How do these covariates shift from ALL forest to CORE forest?")
print(f"  intensity : all {np.median(inten):.0f}  core {np.median(inten[co]):.0f}")
print(f"  multi-return frac (nr>1): all {multi.mean()*100:.1f}%  core {multi[co].mean()*100:.1f}%")
print(f"  |scan angle|: all {np.median(ang):.0f}  core {np.median(ang[co]):.0f}")

print("\n(B) does controlling for intensity+nr remove the all->core shift?")
# match core vs all within intensity x nr cells; compare residual d shift
def strat_shift():
    ie=np.quantile(inten,np.linspace(0,1,6)); ib=np.clip(np.digitize(inten,ie[1:-1]),0,4)
    tot=0.0; w=0.0
    for iv in range(5):
        for nv in [1,2,3,4]:
            m_all=(ib==iv)&(nr==nv); m_co=m_all&co
            if m_co.sum()<200 or (m_all&~co).sum()<200: continue
            sh=med(d[m_co])-med(d[m_all&~co]); wt=m_co.sum(); tot+=sh*wt; w+=wt
    return tot/w if w else np.nan
print(f"  raw all->core shift: {med(d[co])-med(d[~co]):+.1f} mm")
print(f"  shift WITHIN matched intensity x nr strata: {strat_shift():+.1f} mm  (if ~0 => explained by intensity+nr)")

print("\n(C) residual: SINGLE-return (nr=1), mid-intensity (10-24) ground d vs canopy density:")
base=(nr==1)&(inten>=10)&(inten<=24)
pb=pn[base]; db=d[base]
e=np.quantile(pb,np.linspace(0,1,7))
for i in range(len(e)-1):
    b=(pb>=e[i])&(pb<e[i+1] if i<len(e)-2 else pb<=e[i+1])
    if b.sum()<500: continue
    print(f"  pen {e[i]:.3f}-{e[i+1]:.3f} (med {np.median(pb[b]):.3f}): d {med(db[b]):+6.1f} mm  n={b.sum()}")
print(f"  -> single-return fixed-intensity d still sinks with density? corr(d,pen | nr=1,mid-I) = {np.corrcoef(db,pb)[0,1]:+.3f}")
print(f"  ground-return DENSITY per cell may also matter: n class-2 forest returns all vs core computed above")
