"""
Plot MCS track-timestep frequency over the WAM region on a 1x1 degree grid,
masked to land only.

Usage:
    python plot_mcs_track_density.py --model um_glm_n2560_RAL3p3_tuned_hk26
"""
import argparse
from io import BytesIO
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import numpy as np
import requests
import xarray as xr

import models

LON_MIN, LON_MAX = -25, 25
LAT_MIN, LAT_MAX = -2, 20
RES = 1.0


def load_tracks(model):
    print('Loading MCS track statistics...')
    response = requests.get(models.stats_url(model), stream=True)
    return xr.open_dataset(BytesIO(response.content), mask_and_scale=True)


def compute_density(ds):
    lat    = ds.meanlat.values       # (tracks, times)
    lon    = ds.meanlon.values
    lon180 = (lon + 180) % 360 - 180

    # Land filter: same threshold as calc_mcs_entrainment.py (pf_landfrac > 0.8)
    print('Loading pf_landfrac...')
    ds.pf_landfrac.load()
    mean_lf = np.nanmean(ds.pf_landfrac.values, axis=1)  # mean over times
    is_land_track = mean_lf > 0.8   # (tracks,) boolean

    bins_lon = np.arange(LON_MIN, LON_MAX + RES, RES)
    bins_lat = np.arange(LAT_MIN, LAT_MAX + RES, RES)

    # Broadcast land mask to (tracks, times) and combine with spatial bounds
    land_broad = is_land_track[:, np.newaxis]
    valid = (land_broad &
             np.isfinite(lon180) & np.isfinite(lat) &
             (lon180 >= LON_MIN) & (lon180 <= LON_MAX) &
             (lat >= LAT_MIN) & (lat <= LAT_MAX))

    freq, _, _ = np.histogram2d(lon180[valid], lat[valid],
                                bins=[bins_lon, bins_lat])

    lon_centres = (bins_lon[:-1] + bins_lon[1:]) / 2
    lat_centres = (bins_lat[:-1] + bins_lat[1:]) / 2

    return lon_centres, lat_centres, freq.T


def plot_density(lon_centres, lat_centres, freq, model_display, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(9, 5), subplot_kw={'projection': proj},
                           layout='constrained')
    fig.suptitle(f'MCS track-timestep frequency (land MCS, pf_landfrac > 0.8) — {model_display}')

    lon2d, lat2d = np.meshgrid(lon_centres, lat_centres)
    masked = np.ma.masked_where(freq == 0, freq)
    im = ax.pcolormesh(lon2d, lat2d, masked, cmap='YlOrRd', transform=proj)
    plt.colorbar(im, ax=ax, orientation='vertical', label='Track-timesteps per cell')

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.OCEAN, facecolor='lightblue', alpha=0.3)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)

    out = Path(output_dir) / 'mcs_track_density.png'
    fig.savefig(out, dpi=150)
    print(f'Saved {out}')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    models.add_model_arg(parser)
    parser.add_argument('--output', default=None, help='Output directory')
    args = parser.parse_args()

    model_display = models.MODELS[args.model]['display']
    output_dir    = args.output or str(models.figs_dir(args.model))

    ds = load_tracks(args.model)
    lon_centres, lat_centres, freq = compute_density(ds)
    plot_density(lon_centres, lat_centres, freq, model_display, output_dir)


if __name__ == '__main__':
    main()
