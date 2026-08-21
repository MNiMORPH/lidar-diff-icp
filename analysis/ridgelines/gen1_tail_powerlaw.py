#!/usr/bin/env python3
"""Are gen1's ground-return tails POWER LAWS?  The core+halo (2-Gaussian) fit was a hack for
heavy tails; the excess kurtosis (+5.8 forest, +11.9 open) says the tails are not Gaussian.
Here we (1) discriminate tail FORM and (2) fit a proper heavy-tailed law.

Discriminator on the tail region (|d-mode| in [0.12, 0.55] m), upper and lower separately:
  - power law   => log(density) linear in LOG|d-mode|      (slope = -exponent alpha)
  - exponential => log(density) linear in |d-mode|
Compare R^2 of the linear fit in each space; higher R^2 wins.

Fit (binned Poisson MLE), compare AIC:
  Gaussian(mu,sig)   |   Student-t(mu,sig,nu)  [tails ~ |x|^-(nu+1), a power law]
Student-t nu -> small = heavy power-law tails; nu -> inf = Gaussian.  Its tail exponent
alpha_t = nu+1 should match the measured log-log slope.

    ./lidar-icp/bin/python analysis/ridgelines/gen1_tail_powerlaw.py
"""
import numpy as np
from scipy.optimize import minimize
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]

def loglin(u, y):                       # linear fit y = a + b*u, return b, R^2
    b,a = np.polyfit(u, y, 1); yh=a+b*u
    r2 = 1 - np.sum((y-yh)**2)/np.sum((y-np.mean(y))**2); return b, r2

def tail_form(c, mode, side):
    d = c/c.sum()/FW
    if side=="low": m=(fc<=mode-0.12)&(fc>=mode-0.55)&(c>0); u=mode-fc[m]
    else:           m=(fc>=mode+0.12)&(fc<=mode+0.55)&(c>0); u=fc[m]-mode
    y=np.log(d[m])
    b_exp,r2_exp = loglin(u, y)                 # semilog: exponential
    b_pow,r2_pow = loglin(np.log(u), y)         # log-log: power law
    return dict(side=side, n=m.sum(), lam=-1/b_exp, r2_exp=r2_exp,
                alpha=-b_pow, r2_pow=r2_pow, form=("power-law" if r2_pow>r2_exp else "exponential"))

def fit(key, lbl):
    c=Z[key].astype(float); N=c.sum(); mode=fc[np.argmax(c)]
    def nll(f): ff=np.clip(f,1e-300,None)*FW; return -(c*np.log(ff)).sum()
    rG=minimize(lambda x:nll(stats.norm.pdf(fc,x[0],np.exp(x[1]))),[mode,np.log(0.1)],
                method="Nelder-Mead",options=dict(maxiter=20000,xatol=1e-7,fatol=1e-3))
    aG=4+2*rG.fun
    rT=minimize(lambda x:nll(stats.t.pdf(fc,np.exp(x[2])+1.0,loc=x[0],scale=np.exp(x[1]))),
                [mode,np.log(0.06),np.log(2.0)],method="Nelder-Mead",
                options=dict(maxiter=20000,xatol=1e-7,fatol=1e-3))
    muT,sigT,nuT=rT.x[0],np.exp(rT.x[1]),np.exp(rT.x[2])+1.0; aT=6+2*rT.fun
    tl=tail_form(c,mode,"low"); tu=tail_form(c,mode,"high")
    print(f"\n=== gen1 {lbl} (n={N:,.0f}, mode {mode:+.3f}) ===")
    print(f"  Gaussian  sig {np.exp(rG.x[1])*1000:.0f} mm            AIC {aG:,.0f}")
    print(f"  Student-t mu {muT:+.4f}  sig {sigT*1000:.0f} mm  nu {nuT:.2f}  "
          f"(tail exponent nu+1 = {nuT+1:.2f})  AIC {aT:,.0f}  dAIC {aT-aG:+,.0f}")
    for t in (tl,tu):
        print(f"  {t['side']:4s} tail (n={t['n']}): power-law alpha {t['alpha']:.2f} (R2 {t['r2_pow']:.3f}) | "
              f"exp lambda {t['lam']*1000:.0f}mm (R2 {t['r2_exp']:.3f})  -> {t['form'].upper()}")
    return dict(c=c,N=N,mode=mode,muT=muT,sigT=sigT,nuT=nuT,tl=tl,tu=tu,sigG=np.exp(rG.x[1]))

print("Testing tail form (power law vs exponential vs Gaussian) on gen1 ground returns:")
F=fit("gen1_1","FOREST"); O=fit("gen1_2","OPEN")

# ---- figure: log-log tails (power-law test) + Student-t overlay ---------------------
fig,axes=plt.subplots(2,2,figsize=(13,11))
for col,(R,lbl) in enumerate([(F,"FOREST"),(O,"OPEN")]):
    c=R["c"]; d=c/R["N"]/FW; mode=R["mode"]
    # top row: log-log both tails with power-law slope
    ax=axes[0,col]
    for side,col2,t in [("low","C0",R["tl"]),("high","C3",R["tu"])]:
        if side=="low": m=(fc<mode)&(c>0); u=mode-fc[m]
        else:           m=(fc>mode)&(c>0); u=fc[m]-mode
        ax.loglog(u,d[m],".",ms=3,color=col2,label=f"{side} tail (alpha={t['alpha']:.2f}, R2={t['r2_pow']:.2f})")
        uu=np.array([0.12,0.55]); A=d[m][np.argmin(np.abs(u-0.12))]
        ax.loglog(uu,A*(uu/0.12)**(-t['alpha']),"-",color=col2,lw=1,alpha=.6)
    ax.set_xlabel("|d - mode| (m)"); ax.set_ylabel("density (1/m)")
    ax.set_title(f"gen1 {lbl}: tails in log-log (straight => power law)"); ax.legend(fontsize=8); ax.grid(alpha=.3,which="both")
    # bottom row: linear-d, log-density, data vs Gaussian vs Student-t
    ax=axes[1,col]
    ax.semilogx(np.where(d>0,d,np.nan),fc,"C0",lw=1.6,label="data")
    ax.semilogx(stats.norm.pdf(fc,mode,R["sigG"]),fc,"C2--",lw=1.2,label=f"Gaussian sig {R['sigG']*1000:.0f}mm")
    ax.semilogx(stats.t.pdf(fc,R["nuT"],loc=R["muT"],scale=R["sigT"]),fc,"k--",lw=1.4,
                label=f"Student-t nu={R['nuT']:.1f}")
    ax.axhline(0,color="k",lw=.5); ax.set_ylim(-0.6,0.5); ax.set_xlim(1e-2,2e1)
    ax.set_xlabel("density (1/m, log)"); ax.set_ylabel("slope-normal d (m)")
    ax.set_title(f"gen1 {lbl}: Gaussian misses tails, Student-t captures them"); ax.legend(fontsize=8); ax.grid(alpha=.3,which="both")
fig.suptitle("gen1 ground-return tails: power-law heavy tails, fit by a Student-t (not stacked Gaussians)",y=0.995)
fig.savefig("figures/refdatum/gen1_tail_powerlaw.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_tail_powerlaw.png")
