#!/usr/bin/env python3
"""Find the functional form of gen1's ground-return distribution WITHOUT imposing one.

The empirical log-derivative  s(d) = d ln f / dd  names the form region-by-region:
    Gaussian        f ~ exp(-(d-mu)^2/2sig^2)  ->  s(d) = -(d-mu)/sig^2   (STRAIGHT LINE)
    exponential     f ~ exp(-|d-mu|/lam)        ->  s(d) = -/+ 1/lam        (FLAT PLATEAU)
    power law       f ~ |d-c|^-alpha            ->  s(d) = -alpha/(d-c)     (1/x HYPERBOLA)
So plotting s(d) vs d shows: straight => Gaussian, plateau => exponential tail, hyperbola
=> power-law tail.  We read the parameters off the plateaus/slopes; we do not fit a shape.

s(d) is estimated by a count-WEIGHTED local linear regression of ln f (so sparse tail bins,
which are noisy, get little weight).  Bandwidth h in metres.

    ./lidar-icp/bin/python analysis/ridgelines/gen1_form_logderiv.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]

def logderiv(c, h=0.035):
    """count-weighted local-linear slope of ln(density) at each bin centre."""
    N=c.sum(); f=c/N/FW
    ok=c>0; x=fc[ok]; y=np.log(f[ok]); w0=c[ok]            # weight by counts
    s=np.full(fc.shape,np.nan); se=np.full(fc.shape,np.nan)
    for i,d0 in enumerate(fc):
        k=np.exp(-0.5*((x-d0)/h)**2)*w0
        if k.sum()<50: continue
        X=x-d0; Sw=k.sum(); Sx=(k*X).sum(); Sxx=(k*X*X).sum()
        Sy=(k*y).sum(); Sxy=(k*X*y).sum()
        den=Sw*Sxx-Sx*Sx
        if abs(den)<1e-12: continue
        b=(Sw*Sxy-Sx*Sy)/den                              # local slope = s(d)
        s[i]=b
    return s

def analyze(key,lbl):
    c=Z[key].astype(float); mode=fc[np.argmax(c)]; s=logderiv(c)
    # read forms: core linearity, tail plateaus, tail hyperbola
    core=(fc>mode-0.08)&(fc<mode+0.08)&np.isfinite(s)
    A=np.polyfit(fc[core],s[core],1)                      # s ~ A0*d + A1 ; Gaussian => sig=sqrt(-1/A0)
    sigG=np.sqrt(-1/A[0]) if A[0]<0 else np.nan
    lo=(fc<mode-0.18)&(fc>mode-0.45)&np.isfinite(s)       # lower tail region
    hi=(fc>mode+0.18)&(fc<mode+0.45)&np.isfinite(s)
    print(f"\n=== gen1 {lbl} (mode {mode:+.3f}) ===")
    print(f"  CORE |d-mode|<0.08: s(d) slope {A[0]:+.1f} /m  -> if Gaussian, sig = {sigG*1000:.0f} mm; "
          f"curvature check below")
    for m,side in [(lo,'low'),(hi,'high')]:
        if m.sum()<3: continue
        sm=s[m]; dm=fc[m]-mode
        plateau=np.median(sm)                             # if flat => exponential 1/lam
        lam=1/abs(plateau) if plateau!=0 else np.inf
        # power-law test: is s ~ -alpha/(d-mode)? then s*(d-mode) = -alpha = const
        alpha_prod=np.median(sm*dm)
        cv_plateau=np.std(sm)/abs(np.mean(sm))            # flat if small
        cv_hyper=np.std(sm*dm)/abs(np.mean(sm*dm))        # 1/x if this is small
        form="EXPONENTIAL (flat s)" if cv_plateau<cv_hyper else "POWER-LAW (s~1/d)"
        print(f"  {side:4s} tail |d-mode| 0.18-0.45: s~{plateau:+.1f}/m (lam {lam*1000:.0f}mm, CV {cv_plateau:.2f}) | "
              f"s*(d-mode)~{alpha_prod:+.2f} (alpha, CV {cv_hyper:.2f})  -> {form}")
    return dict(c=c,mode=mode,s=s,sigG=sigG,A=A)

F=analyze("gen1_1","FOREST"); O=analyze("gen1_2","OPEN")

# ---- figure: s(d) vs d, with reference forms ----------------------------------------
fig,axes=plt.subplots(1,2,figsize=(14,6))
for ax,(R,lbl) in zip(axes,[(F,"FOREST"),(O,"OPEN")]):
    m=np.isfinite(R["s"]); ax.plot(fc[m],R["s"][m],"C0.-",ms=3,lw=1,label="empirical  s(d)=dln f/dd")
    # Gaussian reference through core
    dd=np.linspace(R["mode"]-0.45,R["mode"]+0.45,200)
    ax.plot(dd,R["A"][0]*dd+R["A"][1],"C2--",lw=1,alpha=.7,
            label=f"Gaussian (straight, sig~{R['sigG']*1000:.0f}mm)")
    ax.axvline(R["mode"],color="k",ls=":",lw=.7); ax.axhline(0,color="k",lw=.5)
    ax.set_xlim(R["mode"]-0.5,R["mode"]+0.5); ax.set_ylim(-40,40)
    ax.set_xlabel("slope-normal d (m)"); ax.set_ylabel("s(d) = d ln f / dd  (1/m)")
    ax.set_title(f"gen1 {lbl}: straight=Gaussian, plateau=exp, hyperbola=power-law")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.suptitle("Reading the functional form from the log-derivative (no imposed shape)",y=1.0)
fig.savefig("figures/refdatum/gen1_form_logderiv.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_form_logderiv.png")
