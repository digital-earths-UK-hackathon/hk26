import huracanpy
import xarray as xr
import numpy as np

hk26_colours = {
    '10kmCoMorph': 'violet',
    '5kmCoMorph': 'darkviolet',
    '10kmGAL': 'b',
    '5kmRAL' : 'darkorange',
    'IBTrACS': 'k',
}

hk26_linestyles = {
    '10kmCoMorph': '-',
    '5kmCoMorph': '--',
    '10kmGAL': '-',
    '5kmRAL' : '--',
    'IBTrACS' : '-',
}

def compute_intensification_rate(tracks, var, timedelta):
    """This function computes the 24-hr intensification rate for each track in the provided tracks object.
    Works with either wind or pressure variable
    Uses huracanpy's built in get_rate function to compute the rate of change of
    variable between adjacent records (units specified) for the same track_id.
    The result is assigned to a new variable 'wind_ir_24hr' or 'pres_ir_24hr' in the tracks object.
    tracks: huracanpy tracks dataset
    var (str): name of variable to get rate of change for ('sfcwind_max' or 'psl_min')
    timedelta (int): data frequency in hours
    """
    dt = timedelta * 3600  # 'timedelta' hours in seconds
    
    if var == 'sfcwind_max':
        output_var = 'wind_ir_24hr'
        final_units = 'knots/24hr'
    elif var == 'psl_min':
        output_var = 'pres_ir_24hr'
        final_units = 'hPa/24hr'
        
    rate = tracks.hrcn.get_rate(var_name=var, centering='backward')
    # units of rate are s^-1
    intens_24h = (rate.groupby(tracks.track_id).map(
            lambda x: (x.rolling(record=int(24/timedelta), min_periods=int(24/timedelta)).sum() * dt)))
    
    # Assign with proper variable name and add units attribute
    tracks = tracks.assign(**{output_var: intens_24h})
    tracks[output_var].attrs['units'] = final_units
    print(f'24hr intensification rate calculated, units are {final_units}')    
    return tracks

def load_syclops():
    """
    Helper function to load tracked TC data from Met Office Model Hierarchy and JTWC subset of IBTrACS
    """
    track_dir = '/gws/nopw/j04/kscale/USERS/cscullio/DYAMOND3/data_reruns/syclops/'
    dataset_paths = {'10kmCoMorph': 'n1280CoMA9v2',
                     '10kmGAL': 'n1280GAL9v2', 
                     '5kmCoMorph': 'n2560CoMA9', 
                     '5kmRAL': 'n2560RAL3p3regridv2'}
    tracks = {}
    for name, path in dataset_paths.items():
        model_track = huracanpy.load(track_dir + path + '/' + path + '_track_20200201-20210301_6h_tc_psl.csv')
        # drop storms that form before March 2020 so we have exactly 12 months of data and no spin-up artefacts. 
        model_track = huracanpy.trackswhere(model_track, model_track.track_id, lambda x: x.time[0] >= np.datetime64('2020-03-01'))
        # convert units from m/s to knots when loading in for consistency with IBTrACS
        model_track['sfcwind_max'] *= 1.944
        model_track['psl_min'] *= 1/100
        tracks[name] = model_track 
        print(f'{name} TC tracks loaded successfully')
        
    # also add IBTrACS data from huracanpy
    obs_track = huracanpy.load(source="ibtracs", ibtracs_subset="usa")
    # subset storms that overlap with model dates + prune timestamps to match model
    obs_tracks_overlap = huracanpy.trackswhere(obs_track, obs_track.track_id, lambda x:(x.time[0] > np.datetime64('2020-03-01')) &
                                               (x.time[0] < np.datetime64('2021-03-01')))
    model_times = tracks['10kmCoMorph'].time.dt.time
    obs_tracks_clean = obs_tracks_overlap.where(obs_tracks_overlap.time.dt.time.isin([model_times]), drop=True)  
    obs_tracks_clean = obs_tracks_clean.rename_vars({'wind' : 'sfcwind_max' , 'slp' : 'psl_min'})
    tracks['IBTrACS'] = obs_tracks_clean
    # renamed the IBTrACS variables to be consistent with the model output. 
    print("JTWC tracks loaded successfully")
    return tracks