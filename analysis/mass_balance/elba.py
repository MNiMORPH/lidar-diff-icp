"""Run the routed sediment-continuity budget on Elba; save V_acc + sigma; quantify
the down-network bias and the NW / hilltop-deposition concern.

Paths are options so the same budget can be run on an alternative DoD/LoD pair
(e.g. a cover-corrected DoD) on the same grid; the defaults are unchanged.
"""
import argparse, time, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from catchment_dod_balance import dinf_proportions, mass_balance, weighted_accumulation
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--dod", default="data/derived/elba_fulldensity/dod.npy")
ap.add_argument("--lod", default="data/derived/elba_fulldensity/lod.npy")
ap.add_argument("--dem", default="data/derived/elba_fulldensity/z_after.npy")
ap.add_argument("--out-prefix", default="data/derived/elba_fulldensity/V_acc",
                help="writes <prefix>.npy and <prefix>_sigma.npy")
ap.add_argument("--fig", default="figures/elba_fulldensity_massbalance.png")
args = ap.parse_args()

X0,Y0,X1,Y1=577492.8,4882737.6,580032.8,4886237.6; res=5.0; area=res*res
dod=np.load(args.dod); lod=np.load(args.lod); dem=np.load(args.dem)
print(f"DoD {args.dod} | LoD {args.lod} | DEM {args.dem}",flush=True)
ny,nx=dod.shape; perror=np.maximum(lod/1.96,1e-6)
print(f"grid {ny}x{nx} = {ny*nx} cells; finite DoD {np.isfinite(dod).sum()}",flush=True)

t=time.time(); props,valid=dinf_proportions(dem,breach=True); print(f"routing (breach+Dinf): {time.time()-t:.1f}s, valid {int(valid.sum())}",flush=True)
t=time.time(); out=mass_balance(dod,perror,props,valid,res); print(f"mass_balance: {time.time()-t:.1f}s",flush=True)
V=out["V_acc"]; sig=out["sigma_Vacc"]; contam=out["contaminated"]; surplus=out["surplus"]
np.save(args.out_prefix+".npy",V); np.save(args.out_prefix+"_sigma.npy",sig)
print(f"saved {args.out_prefix}.npy , {args.out_prefix}_sigma.npy",flush=True)

vb=valid & np.isfinite(dod)
evaluable=vb & ~contam
darea,_=weighted_accumulation(vb.astype(float),props,valid)     # drainage area in cells
print(f"\ncontaminated (edge/hole upstream): {100*contam[vb].mean():.0f}% of valid | evaluable: {int(evaluable.sum())}",flush=True)
print(f"scene net volume change  sum(DoD*area) = {np.nansum(np.where(vb,dod*area,0)):+.0f} m^3  ({np.nanmean(dod[vb])*1000:+.1f} mm mean DoD)",flush=True)
print(f"surplus (unphysical deposition, |V|>1.96sigma & V>0): {100*surplus[evaluable].mean():.1f}% of evaluable",flush=True)
z=np.where(sig>0,V/sig,np.nan)
print(f"V_acc/sigma over evaluable: median {np.nanmedian(z[evaluable]):+.2f}, %positive {100*np.mean(z[evaluable]>0):.0f}, %>+1.96 {100*np.mean(z[evaluable]>1.96):.0f}",flush=True)

print("\nDOWN-NETWORK PROFILE  (V_acc & V_acc/sigma vs drainage area):",flush=True)
print(f"{'drainage area cells':>20} {'n':>8} {'medV_acc(m3)':>13} {'medV/sig':>9} {'%V>0':>6}",flush=True)
edges=[1,10,100,1000,10000,100000,1e9]
for i in range(len(edges)-1):
    m=evaluable&(darea>=edges[i])&(darea<edges[i+1])
    if m.sum()<30: continue
    print(f"{edges[i]:9.0f}-{edges[i+1]:9.0f} {m.sum():8d} {np.nanmedian(V[m]):13.0f} {np.nanmedian(z[m]):9.2f} {100*np.mean(V[m]>0):6.0f}",flush=True)

print("\nQUADRANT budgets (Andy's NW concern) -- mean DoD & headwater(ridge) deposition:",flush=True)
I,J=np.mgrid[0:ny,0:nx]
head=evaluable&(darea<20)                                   # near-ridge / headwater cells
quad={"NW":(I<ny/2)&(J<nx/2),"NE":(I<ny/2)&(J>=nx/2),"SW":(I>=ny/2)&(J<nx/2),"SE":(I>=ny/2)&(J>=nx/2)}
print(f"{'quad':>5} {'meanDoD(mm)':>11} {'headwater meanDoD(mm)':>22} {'%head deposition':>17}",flush=True)
for q,mk in quad.items():
    mm=vb&mk; hh=head&mk
    print(f"{q:>5} {1000*np.nanmean(dod[mm]):11.1f} {1000*np.nanmean(dod[hh]):22.1f} {100*np.mean(dod[hh]>0):17.0f}",flush=True)

hs=hillshade(dem,res,X0,Y0); ext=(X0,X1,Y0,Y1)
fig,ax=plt.subplots(1,2,figsize=(20,13))
ax[0].imshow(hs,extent=ext,origin='lower',cmap='gray')
im=ax[0].imshow(np.where(evaluable,z,np.nan),extent=ext,origin='lower',cmap='RdBu',vmin=-4,vmax=4,alpha=0.7)
ax[0].set_title("V_acc / sigma  (blue>0 = deposition surplus / bias; red = erosion-dominated, physical)"); fig.colorbar(im,ax=ax[0],shrink=.5,label="V_acc/sigma")
ax[1].imshow(hs,extent=ext,origin='lower',cmap='gray')
yy,xx=np.where(surplus); ax[1].scatter(X0+(xx+.5)*res,Y0+(yy+.5)*res,s=3,c='blue',marker='s',label=f'surplus ({int(surplus.sum())})')
yy,xx=np.where(contam); ax[1].scatter(X0+(xx+.5)*res,Y0+(yy+.5)*res,s=1,c='0.6',marker='s',alpha=.3,label='contaminated (excluded)')
ax[1].legend(loc='upper left'); ax[1].set_title("unphysical-deposition flag & excluded cells")
for a in ax: a.set_xlabel("Easting")
fig.savefig(args.fig,dpi=110,bbox_inches='tight'); print(f"\nwrote {args.fig}",flush=True)
