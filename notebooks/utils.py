import math as maths

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import easygems.healpix as egh


def hp_mods(ds):
    """Convert from CF-compliant to be compatible with egh, and attach lat/lon coords"""
    return ds.rename({'healpix_index': 'cell'}).pipe(egh.attach_coords)


def plot_all_fields(ds_plot):
    """Plot all fields for a given dataset. Assumes that each field is 2D - i.e. sel(time=..., [pressure=...]) has been applied"""
    zoom = ds_plot.crs.attrs['refinement_level']
    projection = ccrs.Robinson(central_longitude=0)
    # Do not plot orog, land surf.
    plot_vars = [(name, da) for name, da in ds_plot.data_vars.items() if name not in {'orog', 'sftlf', 'weights'}]
    rows = maths.ceil(len(plot_vars) / 5)
    fig, axes = plt.subplots(rows, 5, figsize=(30, rows * 20 / 6), subplot_kw={'projection': projection}, layout='constrained')
    if 'pressure' in ds_plot.coords:
        plt.suptitle(f'{ds_plot.simulation} z{zoom} @{float(ds_plot.pressure)}hPa')
    else:
        plt.suptitle(f'{ds_plot.simulation} z{zoom}')

    for ax, (name, da) in zip(axes.flatten(), plot_vars):
        if name == 'mrsol':
            da = da.isel(depth=0)
            name = 'mrsol@depth=0'
        time = pd.Timestamp(ds_plot.time.values.item())

        if abs(da.max() + da.min()) / (da.max() - da.min()) < 0.5:
            # data looks like it needs a diverging cmap.
            # figure out some nice bounds.
            pl, pu = np.percentile(da.values[~np.isnan(da.values)], [2, 98])
            vmax = np.abs([pl, pu]).max()
            kwargs = dict(
                cmap='bwr',
                vmin=-vmax,
                vmax=vmax,
            )
        else:
            kwargs = {}
        ax.set_title(f'time: {time} - {name}')
        ax.set_global()
        im = egh.healpix_show(da, ax=ax, **kwargs);
        long_name = da.long_name

        plt.colorbar(im, label=f'{long_name} ({da.attrs.get("units", "-")})')
        ax.coastlines()


def haversine(lat1, lon1, lat2, lon2):
    # Earth radius in kilometers
    R = 6371.0

    # Convert degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Differences
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = np.sin(dlat / 2)**2 + \
        np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2

    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return R * c