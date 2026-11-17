"""Both RS3 rovers under base B1 (13:20-13:30 CDT): the only minutes both teams were on the same ground.
B1's position error is common to both rovers, so rover-minus-rover is valid; the lidar is read at the raw
position and, as a sensitivity, at the position shifted by the class base B's tie (B1's own error is unknown)."""
import pandas as pd, numpy as np, pyproj, glob, os, json, sys
RAW="data/field/2026-10-25/raw/"
b=pd.read_csv(RAW+"rover_RS3_82434D3F6234D057_hillslope_T4a_T4b.csv"); b["team"]="T4a"
m=pd.read_csv(RAW+"rover_RS3_8243D2AFB94BDAE3_hillslope_T3a_T3b_T3c.csv"); m["team"]="T3a"
d=pd.concat([b,m],ignore_index=True)
d=d[np.isclose(d["Base longitude"],-92.000382,atol=1e-4)].copy()
to_utm=pyproj.Transformer.from_crs("EPSG:6318","EPSG:6344",always_xy=True)
d["E"],d["N"]=to_utm.transform(d.Longitude.values,d.Latitude.values)
tg=pyproj.Transformer.from_crs("EPSG:6319","EPSG:6318+5703",always_xy=True)
_,_,d["navd88"]=tg.transform(d.Longitude.values,d.Latitude.values,d["Ellipsoidal height"].values)
chk=tg.transform(-91.99999,44.12,188.669)[2]; print("geoid check (an arbitrary nearby lon/lat, 188.669 ell):",round(chk,3),"(expect ~218.7 if GEOID18 loaded)")
tiles={}
for c in glob.glob("data/derived/4*/corrections.json"):
    t=os.path.basename(os.path.dirname(c)); j=json.load(open(c)); tiles[t]=(j["bounds"],float(j.get("res") or j.get("res_m")))
def sample(name,x,y):
    for t,(bb,res) in tiles.items():
        if bb[0]<=x<bb[2] and bb[1]<=y<bb[3]:
            a=np.load(f"data/derived/{t}/{name}",mmap_mode="r")
            fy,fx=(y-bb[1])/res-0.5,(x-bb[0])/res-0.5; i0,j0=int(np.floor(fy)),int(np.floor(fx)); wy,wx=fy-i0,fx-j0
            i0=min(max(i0,0),a.shape[0]-2); j0=min(max(j0,0),a.shape[1]-2); q=np.asarray(a[i0:i0+2,j0:j0+2],float)
            if not np.isfinite(q).all(): return float(a[int((y-bb[1])/res),int((x-bb[0])/res)])
            return float((1-wy)*((1-wx)*q[0,0]+wx*q[0,1])+wy*((1-wx)*q[1,0]+wx*q[1,1]))
    return np.nan
for tag,(dE,dN) in (("raw",(0,0)),("Bshift",(0.771,-1.112))):
    d[f"z21_{tag}"]=[sample("z_after.npy",x+dE,y+dN) for x,y in zip(d.E,d.N)]
    d[f"l21mr_{tag}"]=d[f"z21_{tag}"]-d.navd88
d["slope"]=[sample("slope.npy",x,y) for x,y in zip(d.E,d.N)]
d["lowveg"]=[sample("lowveg_grid.npy",x,y) for x,y in zip(d.E,d.N)]
pd.set_option("display.width",250)
print(d[["team","Name","Averaging start","Solution status","Samples","E","N","Ellipsoidal height","navd88","z21_raw","l21mr_raw","l21mr_Bshift","slope","lowveg"]].round(3).to_string())
print("\nteam medians of 2021 lidar minus RTK (mm), B1 points, n:")
for t,s in d.groupby("team"):
    print(f"  {t:8s} raw {1000*s.l21mr_raw.median():+.0f}  Bshift {1000*s.l21mr_Bshift.median():+.0f}  n={len(s)}  (excluding the 3 lot points 579982 E: ", end="")
    s2=s[s.E<579950]; print(f"raw {1000*s2.l21mr_raw.median():+.0f} Bshift {1000*s2.l21mr_Bshift.median():+.0f} n={len(s2)})")
print("\nnearest cross-team pairs (B1 points): distance, time gap, RTK dh (T4a-T3a), lidar dz at the two spots, dh-dz = rover difference on the same ground:")
B=d[d.team=="T4a"]; M=d[d.team=="T3a"]
for _,r in M.iterrows():
    dist=np.hypot(B.E-r.E,B.N-r.N); k=dist.idxmin(); q=B.loc[k]
    dt=(pd.to_datetime(q["Averaging start"].split(" UTC")[0])-pd.to_datetime(r["Averaging start"].split(" UTC")[0])).total_seconds()
    dh=q["Ellipsoidal height"]-r["Ellipsoidal height"]; dz=q.z21_raw-r.z21_raw; dzs=q.z21_Bshift-r.z21_Bshift
    print(f"  T3a {int(r.Name)} / T4a {int(q.Name)}: {dist.min():5.2f} m, {dt:+5.0f} s, dh {1000*dh:+6.0f} mm, lidar dz raw {1000*dz:+6.0f} / Bshift {1000*dzs:+6.0f}, dh-dz raw {1000*(dh-dz):+5.0f} / Bshift {1000*(dh-dzs):+5.0f} mm;  l21-RTK: T3a {1000*r.l21mr_raw:+5.0f} T4a {1000*q.l21mr_raw:+5.0f}")
d.to_csv(sys.argv[1] if len(sys.argv)>1 else "/dev/null",index=False)
