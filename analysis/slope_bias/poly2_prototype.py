"""Prototype + synthetic validation of a WINDOWED 2nd-order-polynomial ground gridder.
Per 5 m cell, fit z = a + b*u + c*v + d*u^2 + e*v^2 + f*uv to ground points in the 3x3
(15 m) window (u,v = offset from THIS cell's centre, in cell units for conditioning),
read the CONSTANT term a = surface value AT cell centre -> curvature-unbiased.
Windowed moments via 9 shifts (di,dj in -1,0,1): a point in its own cell contributes to
each target cell with offset (u0-dj, v0-di). Validate before porting to the pipeline."""
import numpy as np

def poly2_ground(x, y, z, X0, Y0, res, nx, ny, minpts=18):
    ix=((x-X0)/res).astype(int); iy=((y-Y0)/res).astype(int)
    ok=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny); ix,iy,x,y,z=ix[ok],iy[ok],x[ok],y[ok],z[ok]
    u0=(x-(X0+(ix+0.5)*res))/res; v0=(y-(Y0+(iy+0.5)*res))/res   # cell units
    N=nx*ny; pairs=[(a,b) for a in range(6) for b in range(a,6)]  # 21 upper-tri
    M=[np.zeros(N) for _ in range(21)]; R=[np.zeros(N) for _ in range(6)]
    for di in (-1,0,1):
        for dj in (-1,0,1):
            ti=iy+di; tj=ix+dj; m=(ti>=0)&(ti<ny)&(tj>=0)&(tj<nx)
            t=ti[m]*nx+tj[m]; u=u0[m]-dj; v=v0[m]-di; zz=z[m]
            phi=[np.ones_like(u),u,v,u*u,v*v,u*v]
            for k,(a,b) in enumerate(pairs): M[k]+=np.bincount(t,phi[a]*phi[b],N)
            for k in range(6): R[k]+=np.bincount(t,phi[k]*zz,N)
    cnt=M[0]; a0=np.full(N,np.nan)
    idx=np.where(cnt>=minpts)[0]
    if len(idx):
        Mm=np.zeros((len(idx),6,6))
        for k,(a,b) in enumerate(pairs): Mm[:,a,b]=M[k][idx]; Mm[:,b,a]=M[k][idx]
        rhs=np.stack([R[k][idx] for k in range(6)],1)
        det=np.linalg.det(Mm); good=np.abs(det)>1e-6
        if good.any(): a0[idx[good]]=np.linalg.solve(Mm[good],rhs[good])[:,0]
    return a0.reshape(ny,nx)

# ---- synthetic validation ----
rng=np.random.default_rng(0); X0=Y0=0.0; res=5.0; nx=ny=20
n=40000; x=rng.uniform(0,nx*res,n); y=rng.uniform(0,ny*res,n)
def cell_center_field(fn):
    xc=X0+(np.arange(nx)+0.5)*res; yc=Y0+(np.arange(ny)+0.5)*res
    XX,YY=np.meshgrid(xc,yc); return fn(XX,YY)
print("SYNTHETIC CHECKS (poly2 should recover the value AT cell centre):")
# 1 planar tilt
z=100+0.03*x-0.02*y; g=poly2_ground(x,y,z,X0,Y0,res,nx,ny); true=cell_center_field(lambda X,Y:100+0.03*X-0.02*Y)
print(f"  planar : max|poly2-truth| = {np.nanmax(np.abs(g-true)):.4f} m (expect ~0)")
# 2 paraboloid (convex bowl): median-grid is biased, poly2 must hit centre exactly
cx=cy=nx*res/2; z=0.004*((x-cx)**2+(y-cy)**2); g=poly2_ground(x,y,z,X0,Y0,res,nx,ny)
true=cell_center_field(lambda X,Y:0.004*((X-cx)**2+(Y-cy)**2))
# median grid for comparison
import pandas as pd
ixp=((x-X0)/res).astype(int); iyp=((y-Y0)/res).astype(int)
med=pd.Series(z).groupby(iyp*nx+ixp).median(); mg=np.full(nx*ny,np.nan); mg[med.index]=med.values; mg=mg.reshape(ny,ny)
inb=np.isfinite(g)&np.isfinite(true)
print(f"  parab  : max|poly2-truth| = {np.nanmax(np.abs(g-true)[inb]):.4f} m (expect ~0, poly2 exact on a quadratic)")
print(f"  parab  : median-grid bias vs truth = {np.nanmean((mg-true)[inb]):+.4f} m (should be NONZERO -> the curvature bias poly2 removes)")
# 3 with noise
z=0.004*((x-cx)**2+(y-cy)**2)+rng.normal(0,0.05,n); g=poly2_ground(x,y,z,X0,Y0,res,nx,ny)
print(f"  parab+noise(5cm): median|poly2-truth| = {np.nanmedian(np.abs(g-true)[inb]):.4f} m (small, robust)")
