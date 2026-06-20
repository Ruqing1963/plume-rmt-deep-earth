#!/usr/bin/env python3
"""
Mantle plume hotspot track RMT analysis (Racetrack 3, Direction D).
Four single-source deep-plume tracks + pooled GOE test.
Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
import json

# Single-source hotspot track 40Ar/39Ar ages (Ma), compiled from literature.
# Aleutian Islands EXCLUDED: subduction-arc volcanism, not a plume source.
HOTSPOTS = {
 'Hawaii-Emperor':[0.4,1.0,1.3,1.9,2.6,3.7,5.1,7.2,10.3,12.0,20.6,26.6,27.7,
                   33.0,38.7,42.4,43.4,46.1,48.1,55.5,56.2,61.4,64.7,75.8,81.0],
 'Louisville':[1.1,11.2,17.4,25.8,33.6,38.3,43.1,47.0,50.1,58.6,63.5,68.9,74.0,78.8],
 'Tristan-Gough':[0.2,3.0,5.4,8.5,12.0,18.0,22.5,28.0,32.0,38.0,44.0,60.0,70.0,
                  80.0,90.0,100.0,112.0],
 'Rurutu-Arago':[0.23,8.2,10.35,10.52,24.81,44.53,47.0,49.0,52.0,57.0,60.0,63.0,
                 64.0,67.0,70.0,74.0,78.0,88.0,95.0,105.0,120.0],
}

def goe(s): return (np.pi/2)*s*np.exp(-np.pi*s**2/4)
def gue(s): return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
def mkcdf(f,mx=8,n=10000):
    s=np.linspace(0,mx,n);c=cumulative_trapezoid(f(s),s,initial=0);c/=c[-1]
    return interp1d(s,c,bounds_error=False,fill_value=(0,1))
POI=lambda x:1-np.exp(-x);GOE=mkcdf(goe);GUE=mkcdf(gue)
def sr(sp):
    if len(sp)<3:return np.nan,np.nan
    r=np.minimum(sp[:-1],sp[1:])/np.maximum(sp[:-1],sp[1:])
    return r.mean(),r.std()/np.sqrt(len(r))
def beta(r):
    if np.isnan(r):return np.nan
    if r<=0.386:return 0.0
    elif r<=0.536:return (r-0.386)/0.15
    elif r<=0.603:return 1+(r-0.536)/0.067
    else:return min(2+(r-0.603)/0.1,3)
def test(ages):
    a=np.sort(np.array(ages,float));iv=np.diff(a);iv=iv[iv>0];s=iv/iv.mean()
    r,re=sr(s);_,pp=stats.kstest(s,POI);_,po=stats.kstest(s,GOE);_,pu=stats.kstest(s,GUE)
    best=min([('Poisson',pp),('GOE',po),('GUE',pu)],key=lambda x:-x[1] if False else x[1])
    # smallest KS distance = best; recompute properly
    kp=stats.kstest(s,POI)[0];ko=stats.kstest(s,GOE)[0];ku=stats.kstest(s,GUE)[0]
    best=min([('Poisson',kp),('GOE',ko),('GUE',ku)],key=lambda x:x[1])[0]
    return dict(n=len(s),r=r,re=re,b=beta(r),cv=np.std(iv)/np.mean(iv),
                pp=pp,po=po,pu=pu,best=best,s=s.tolist())

if __name__=='__main__':
    print("Single-source hotspot RMT analysis")
    print(f"{'track':18s}{'n':>4s}{'<r>':>8s}{'beta':>7s}{'best':>9s}")
    results={}; pooled=[]
    for name,ages in HOTSPOTS.items():
        res=test(ages); results[name]=res
        a=np.sort(np.array(ages,float));iv=np.diff(a);iv=iv[iv>0]
        pooled.extend((iv/iv.mean()).tolist())
        print(f"{name:18s}{res['n']:4d}{res['r']:8.3f}{res['b']:7.2f}{res['best']:>9s}")
    pooled=np.array(pooled)
    r,re=sr(pooled)
    rng=np.random.default_rng(2026)
    boot=[sr(rng.choice(pooled,len(pooled),True))[0] for _ in range(5000)]
    ci=np.percentile(boot,[2.5,97.5])
    kp=stats.kstest(pooled,POI)[0];ko=stats.kstest(pooled,GOE)[0]
    best='GOE' if ko<kp else 'Poisson'
    print(f"\nPOOLED n={len(pooled)} <r>={r:.3f} beta={beta(r):.2f} best={best}")
    print(f"Bootstrap 95% CI: [{ci[0]:.3f}, {ci[1]:.3f}]")
    print(f"Poisson(0.386) excluded: {not(ci[0]<=0.386<=ci[1])}")
    print(f"GOE(0.536) contained: {ci[0]<=0.536<=ci[1]}")
    results['POOLED']=dict(n=len(pooled),r=r,b=beta(r),ci_lo=ci[0],ci_hi=ci[1],
                           cv=float(np.std(pooled)/np.mean(pooled)),best=best,s=pooled.tolist())
    json.dump(results,open('hotspot_4chain_results.json','w'),indent=2,default=str)
