#!/usr/bin/env python3
"""Same treatment as gen1 (moments, robust stats, and the log-derivative form-finder) for
gen2 (2021 3DEP, leaf-on/green-up, ~24x denser).  s(d)=d ln f/dd names the form without
imposing one: straight => Gaussian, flat plateau => exponential tail, 1/x hyperbola =>
power-law tail.  Reuses ground_fine_pooled.npz.  Prints a gen1-vs-gen2 form comparison.

    ./lidar-icp/bin/python analysis/ridgelines/gen2_form_logderiv.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]

def wq(c,q): cdf=np.cumsum(c)/c.sum(); return np.interp(q,cdf,fc)
def logderiv(c, h=0.02):                     # tighter bandwidth: gen2 is much narrower
    N=c.sum(); f=c/N/FW; ok=c>0; x=fc[ok]; y=np.log(f[ok]); w0=c[ok]
    s=np.full(fc.shape,np.nan)
    for i,d0 in enumerate(fc):
        k=np.exp(-0.5*((x-d0)/h)**2)*w0
        if k.sum()<50: continue
        X=x-d0; Sw=k.sum(); Sx=(k*X).sum(); Sxx=(k*X*X).sum(); Sy=(k*y).sum(); Sxy=(k*X*y).sum()
        den=Sw*Sxx-Sx*Sx
        if abs(den)>1e-12: s[i]=(Sw*Sxy-Sx*Sy)/den
    return s

def analyze(key,lbl):
    c=Z[key].astype(float); N=c.sum(); p=c/N
    mean=(p*fc).sum(); std=np.sqrt((p*(fc-mean)**2).sum())
    skew=(p*((fc-mean)/std)**3).sum(); kurt=(p*((fc-mean)/std)**4).sum()-3
    med=wq(c,.5); mode=fc[np.argmax(c)]
    order=np.argsort(np.abs(fc-med)); nmad=1.4826*np.interp(.5,np.cumsum(c[order])/N,np.abs(fc-med)[order])
    s=logderiv(c)
    core=(fc>mode-0.05)&(fc<mode+0.05)&np.isfinite(s); A=np.polyfit(fc[core],s[core],1)
    sigG=np.sqrt(-1/A[0]) if A[0]<0 else np.nan
    print(f"\n=== gen2 {lbl} (n={N:,.0f}) ===")
    print(f"  mode {mode:+.3f}  median {med:+.3f}  mean {mean:+.3f} m")
    print(f"  std {std*1000:.0f}  NMAD {nmad*1000:.0f} mm   skew {skew:+.2f}  exkurt {kurt:+.1f}")
    print(f"  pctl: p05 {wq(c,.05):+.3f} p10 {wq(c,.10):+.3f} p25 {wq(c,.25):+.3f} p50 {med:+.3f} "
          f"p75 {wq(c,.75):+.3f} p90 {wq(c,.90):+.3f} p95 {wq(c,.95):+.3f}")
    print(f"  tail balance: low(<mode-0.15) {p[fc<mode-0.15].sum()*100:.1f}%  high(>mode+0.15) {p[fc>mode+0.15].sum()*100:.1f}%")
    print(f"  core |d-mode|<0.05: s slope {A[0]:+.0f}/m -> Gaussian sig {sigG*1000:.0f} mm")
    for side,rng in [("low",(-0.30,-0.12)),("high",(0.12,0.30))]:
        m=(fc>=mode+rng[0])&(fc<=mode+rng[1])&np.isfinite(s)
        if m.sum()<3: print(f"  {side} tail: too few"); continue
        sm=s[m]; dm=fc[m]-mode
        cv_plateau=np.std(sm)/abs(np.mean(sm)); cv_hyper=np.std(sm*dm)/abs(np.mean(sm*dm))
        form="EXPONENTIAL" if cv_plateau<cv_hyper else "POWER-LAW"
        print(f"  {side:4s} tail |d-mode|0.12-0.30: s~{np.median(sm):+.0f}/m (lam {1000/abs(np.median(sm)):.0f}mm,CV {cv_plateau:.2f}) | "
              f"s*(d-m)~{np.median(sm*dm):+.2f}(alpha,CV {cv_hyper:.2f}) -> {form}")
    return dict(c=c,mode=mode,s=s,sigG=sigG,A=A,skew=skew,nmad=nmad,std=std,med=med,N=N)

F=analyze("gen2_1","FOREST"); O=analyze("gen2_2","OPEN")

print("\n=== gen1 vs gen2 (from prior run) ===")
print("  gen1 forest: sigma 173mm, NMAD 134, skew +0.72, SYMMETRIC EXPONENTIAL tails (lam~106)")
print("  gen1 open:   sigma 102mm, NMAD 64,  skew +1.58, ASYMMETRIC POWER-LAW tails (alpha~2.6-3.7)")
print(f"  gen2 forest: sigma {F['std']*1000:.0f}mm, NMAD {F['nmad']*1000:.0f}, skew {F['skew']:+.2f}")
print(f"  gen2 open:   sigma {O['std']*1000:.0f}mm, NMAD {O['nmad']*1000:.0f}, skew {O['skew']:+.2f}")

# ---- figure: s(d) for gen2, with gen1 overlaid for comparison -----------------------
G1=np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
def ld(key,h):
    c=G1[key].astype(float); return logderiv_key(c,h)
def logderiv_key(c,h):
    N=c.sum(); f=c/N/FW; ok=c>0; x=fc[ok]; y=np.log(f[ok]); w0=c[ok]; s=np.full(fc.shape,np.nan)
    for i,d0 in enumerate(fc):
        k=np.exp(-0.5*((x-d0)/h)**2)*w0
        if k.sum()<50: continue
        X=x-d0; Sw=k.sum();Sx=(k*X).sum();Sxx=(k*X*X).sum();Sy=(k*y).sum();Sxy=(k*X*y).sum()
        den=Sw*Sxx-Sx*Sx
        if abs(den)>1e-12: s[i]=(Sw*Sxy-Sx*Sy)/den
    return s
fig,axes=plt.subplots(1,2,figsize=(14,6))
for ax,(R,g1key,lbl) in zip(axes,[(F,"gen1_1","FOREST"),(O,"gen1_2","OPEN")]):
    m=np.isfinite(R["s"]); ax.plot(fc[m]-R["mode"],R["s"][m],"C3.-",ms=3,lw=1,label="gen2  s(d)")
    s1=logderiv_key(G1[g1key].astype(float),0.035); mo1=fc[np.argmax(G1[g1key])]; m1=np.isfinite(s1)
    ax.plot(fc[m1]-mo1,s1[m1],"C0.-",ms=2,lw=.8,alpha=.6,label="gen1  s(d)")
    ax.axhline(0,color="k",lw=.5); ax.axvline(0,color="k",ls=":",lw=.7)
    ax.set_xlim(-0.4,0.4); ax.set_ylim(-60,60)
    ax.set_xlabel("d - mode (m)"); ax.set_ylabel("s(d)=dln f/dd (1/m)")
    ax.set_title(f"{lbl}: gen2 (red) vs gen1 (blue) — steeper core, tail form")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.suptitle("gen2 ground-return form via log-derivative (vs gen1)",y=1.0)
fig.savefig("figures/refdatum/gen2_form_logderiv.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen2_form_logderiv.png")
