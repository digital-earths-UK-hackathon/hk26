import huracanpy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy import stats
import argparse
import iris
import metpy.calc as mpcalc
from metpy.units import units
from scipy import interpolate
import pickle
from pathlib import Path
import classes
import warnings
from multiprocessing import Pool
warnings.filterwarnings("ignore")


def load_file(date,var_name,sim="n1280_10km-CoMA9",relon=True,stream="c",p=slice(None,None)):
    date_str="%04d%02d%02d"%(date.year,date.month,date.day)
    # data in 12 hour chunks, labelled 0 or 12, so get appropriate value.
    tidx=int(date.hour/12)*12

    if sim=="CTC_km4p4_RAL3P3":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_data/10km-GAL9-nest/CTC_km4p4_RAL3P3/field.pp/apver{}.pp/CTC_km4p4_RAL3P3.n1280_GAL9_nest.apver{}_{}T{:02d}.pp".format(
                            stream,stream,date_str,tidx), var_name)[-1]).sel(latitude=lat_range)
    elif sim=="n2560_RAL3p3":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_data/5km-RAL3/glm/field.pp/apver{}.pp/glm.n2560_RAL3p3.apver{}_{}T{:02d}.pp".format(
                            stream,stream,date_str,tidx), var_name)[-1]).sel(latitude=lat_range)
    elif sim=="n2560_RAL3p3_tuned":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_reruns/5km-RAL3p3-tuned/glm/field.pp/apver{}.pp/glm.n2560_RAL3p3_tuned.apver{}_{}T{:02d}00Z.pp".format(
                            stream,stream,date_str,tidx), var_name)[0]).sel(latitude=lat_range)
    elif sim=="n1280_10km-CoMA9":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_reruns/10km-CoMA9/glm/field.pp/apver{}.pp/glm.n1280_CoMA9_v2.apver{}_{}T{:02d}00Z.pp".format(
                            stream,stream,date_str,tidx), var_name)[0]).sel(latitude=lat_range)
    elif sim=="n2560_5km-CoMA9":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_reruns/5km-CoMA9/glm/field.pp/apver{}.pp/glm.n2560_CoMA9_hier_v2.apver{}_{}T{:02d}00Z.pp".format(
                            stream,stream,date_str,tidx), var_name)[0]).sel(latitude=lat_range)
    elif sim=="n1280_GAL9":
        ds = xr.DataArray.from_iris(iris.load("/gws/nopw/j04/kscale/DYAMOND3_reruns/10km-GAL9/glm/field.pp/apver{}.pp/glm.n1280_GAL9_v2.apver{}_{}T{:02d}00Z.pp".format(
                            stream,stream,date_str,tidx), var_name)[-1]).sel(latitude=lat_range)

    if stream=="c" or stream=="d":
        ds=ds.sel(pressure=p)
        
    if relon:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby("longitude").sel(longitude=lon_range)

    return ds


def get_mcs_tracks(dirpath):
    pkl_paths = sorted(dirpath.glob('tracks_*.p'))
    print([p.name for p in pkl_paths])
        
    tracks=[]
    for i in np.arange(5,8): # THIS HARD CODES TO JJA!!
        with pkl_paths[i].open('rb') as f:
            mtracks=pickle.load(f)
            mtracks=[track for track in mtracks if track.is_in_region(dict(lons=lon_bounds,lats=lat_bounds))]
            tracks=tracks+mtracks
    
    df={}
    df["time"]=[storm.time for track in tracks for storm in track.get_storms()]
    df["stormID"]=[storm.storm for track in tracks for storm in track.get_storms()]
    df["clon"]=[storm.centroidlon for track in tracks for storm in track.get_storms()]
    df["clat"]=[storm.centroidlat for track in tracks for storm in track.get_storms()]
    df["area"]=[storm.area for track in tracks for storm in track.get_storms()]
    df["rain_max"]=[storm.maxpr for track in tracks for storm in track.get_storms()]
    df["rain_mean"]=[storm.meanpr for track in tracks for storm in track.get_storms()]
    df["tmin"]=[storm.minTb for track in tracks for storm in track.get_storms()]
    
    df=pd.DataFrame(df)
    df["time"]=pd.to_datetime(df["time"])
    df["area"]=100*df["area"]
    df=df[df["area"]>20000]
    return df


# Get 0.75 degree box centred on lon/lat location
def get_env_var(da,lat,lon,time):
    da=da.sel(latitude=slice(lat-0.375,lat+0.375),longitude=slice(lon-0.375,lon+0.375),time=time)
    return float(da.mean())

def parallelise(idx):
    period_seg=period[idx*fact:(idx+1)*fact]
    tab=[]
    for env_t in period_seg:
        print(env_t)
        u=load_file(env_t.replace(hour=11+offset),"x_wind",p=[600,850,925])
        shear1= u.sel(pressure=850) - u.sel(pressure=600)
        shear2= u.sel(pressure=925) - u.sel(pressure=600)
        q925=1000*load_file(env_t.replace(hour=11+offset),"specific_humidity",p=925)
        tcw=load_file(env_t.replace(hour=11+offset),"m01s30i461",stream="b")
        
        date_storms=env_mcs[env_mcs["time"].dt.dayofyear==env_t.dayofyear]
    
        date_storms["tcw"]=date_storms.apply(lambda x: get_env_var(tcw,x["clat"],x["clon"],env_t),axis=1)
        date_storms["q925"]=date_storms.apply(lambda x: get_env_var(q925,x["clat"],x["clon"],env_t),axis=1)
        date_storms["ushear600_850"]=date_storms.apply(lambda x: get_env_var(shear1,x["clat"],x["clon"],env_t),axis=1)
        date_storms["ushear600_925"]=date_storms.apply(lambda x: get_env_var(shear2,x["clat"],x["clon"],env_t),axis=1)
        tab.append(date_storms)

    tab=pd.concat(tab)
    return tab


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sim", required=True) # 5km-CoMA9 OR 10km-CoMA9
parser.add_argument("-r", "--region", required=False, default="WAf") # "WAf", "US" or "India" - easy to extend
parser.add_argument("-t", "--tracks", default="calum") # "calum" or "ben", referring to different simpleTrack data. ben only for WAf 10km!
args = parser.parse_args()

sim=args.sim
region=args.region
if region=="WAf":
    lat_bounds=(5,25)
    lon_bounds=(-15,25)
    period=pd.date_range("2020-07-01 12:00","2020-09-30 12:00")
    offset=0
elif region=="US":
    lat_bounds=(30,45)
    lon_bounds=(-100,-50)
    offset=-6
    period=pd.date_range(f"2020-07-01 {12+offset}:00",f"2020-09-30 {12+offset}:00")
elif region=="India":
    lat_bounds=(5,30)
    lon_bounds=(70,91)    
    offset=6
    period=pd.date_range(f"2020-07-01 {12+offset}:00",f"2020-09-30 {12+offset}:00")
    
lat_range=slice(lat_bounds[0],lat_bounds[-1])
lon_range=slice(lon_bounds[0],lon_bounds[-1])

um_res=int(2560/(int(sim.split("-")[0][:-2])/5))

if args.tracks=="calum":
    basepath = Path('/gws/nopw/j04/kscale/USERS/cscullio/DYAMOND3/data_reruns/simpleTrack/')
    dirpath_tpl = '{dataset}/'
    track_str=""
    if sim=="5km-CoMA9":
        coma9_mcs_tracks = get_mcs_tracks( basepath / dirpath_tpl.format(dataset=f"n{um_res}_CoMA9_hier_v2"))
    elif sim=="10km-CoMA9":
        coma9_mcs_tracks = get_mcs_tracks( basepath / dirpath_tpl.format(dataset=f"n{um_res}_CoMA9"))
        
    env_mcs = coma9_mcs_tracks[coma9_mcs_tracks.time.dt.hour.isin(np.arange(16+offset,22+offset))]
    print("Loaded MCS tracks")

if args.tracks=="ben":
    mcs=pd.read_csv(f"/gws/nopw/j04/kscale/USERS/bmaybee/simpleTrack_MCS_outputs/dyamond3_n{um_res}_{sim}_ea-waf_MCS_tracks_rain_reformat.csv")
        
    mcs["time"]=pd.to_datetime(mcs.time)
    mcs["rain_mean"]=mcs["rainvol"]/(mcs["area"]/11.1**2)
    env_mcs=mcs[(mcs.area>20000) & (mcs.time.dt.hour.isin(np.arange(16+offset,22+offset)))
                & (mcs.clat>5) & (mcs.clat<25) & (mcs.clon>-15) & (mcs.clon<25)]
    track_str="BMtracks"

psize=16
fact=int(len(period)/(psize-1))
p=Pool(psize)
out=p.map(parallelise,np.arange(psize))
out=pd.concat(out)
out.sort_values(by="time").to_csv(f"~/Huracan/Hackathons/{sim}_{region}_MCS_envs{track_str}.csv",index=False)