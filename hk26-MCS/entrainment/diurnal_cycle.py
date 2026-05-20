"""
Plot the diurnal cycle of CAPE and LNB over a box centred on Burkina Faso,
split by JJA and DJF seasons, alongside a map of the WAM analysis region.

Usage:
    python diurnal_cycle.py
    python diurnal_cycle.py --zarr entrainment_wam.zarr --output diurnal_cycle.png
"""
import argparse

import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import xarray as xr

# WAM region (entrainment calculation domain), in [-180, 180] lon convention
WAM_LAT_MIN, WAM_LAT_MAX = 2, 15
WAM_LON_MIN, WAM_LON_MAX = -20, 20

# Burkina Faso analysis box, in [-180, 180] lon convention
BF_LAT_MIN, BF_LAT_MAX = 9, 15
BF_LON_MIN, BF_LON_MAX = -5, 5


def load_box(zarr_path):
    """Load zarr store and select cells in the Burkina Faso box."""
    ds = xr.open_zarr(zarr_path)
    # Convert [0, 360] lon to [-180, 180] for masking
    lon180 = (ds.lon + 180) % 360 - 180
    mask = (
        (ds.lat > BF_LAT_MIN) & (ds.lat < BF_LAT_MAX) &
        (lon180 > BF_LON_MIN) & (lon180 < BF_LON_MAX)
    )
    mask = mask.compute()
    print(f'Cells in Burkina Faso box: {mask.values.sum()}')
    return ds, ds.isel(cell=mask)


def diurnal_cycle(da):
    return da.groupby('time.hour').mean(), da.groupby('time.hour').std()


def compute_diurnal_cycles(zarr_path):
    _, ds_box = load_box(zarr_path)
    ds_mean = ds_box[['cape', 'lnb']].mean(dim='cell').compute()

    jja = ds_mean.time.dt.month.isin([6, 7, 8])
    djf = ds_mean.time.dt.month.isin([12, 1, 2])

    return {
        'cape_jja': diurnal_cycle(ds_mean.cape.isel(time=jja)),
        'cape_djf': diurnal_cycle(ds_mean.cape.isel(time=djf)),
        'lnb_jja':  diurnal_cycle(ds_mean.lnb.isel(time=jja)),
        'lnb_djf':  diurnal_cycle(ds_mean.lnb.isel(time=djf)),
    }


def plot_region_map(ax):
    """Map showing the WAM region and Burkina Faso box."""
    ax.set_extent([WAM_LON_MIN - 5, WAM_LON_MAX + 5, WAM_LAT_MIN - 3, WAM_LAT_MAX + 3],
                  crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND, facecolor='lightgrey', alpha=0.4)

    crs = ccrs.PlateCarree()

    # WAM region
    wam = mpatches.Rectangle(
        (WAM_LON_MIN, WAM_LAT_MIN),
        WAM_LON_MAX - WAM_LON_MIN, WAM_LAT_MAX - WAM_LAT_MIN,
        linewidth=1.5, edgecolor='steelblue', facecolor='steelblue',
        alpha=0.15, transform=crs, label='WAM region',
    )
    ax.add_patch(wam)

    # Burkina Faso box
    bf = mpatches.Rectangle(
        (BF_LON_MIN, BF_LAT_MIN),
        BF_LON_MAX - BF_LON_MIN, BF_LAT_MAX - BF_LAT_MIN,
        linewidth=2, edgecolor='tab:orange', facecolor='none',
        transform=crs, label='Burkina Faso box',
    )
    ax.add_patch(bf)

    ax.legend(loc='lower right', fontsize=8)
    ax.set_title('Analysis regions')


def plot(cycles, output):
    hours = cycles['cape_jja'][0].hour.values
    stem = output.rsplit('.', 1)
    map_output = f'{stem[0]}_map.{stem[1]}' if len(stem) == 2 else output + '_map'
    dc_output = f'{stem[0]}_dc.{stem[1]}' if len(stem) == 2 else output + '_dc'

    # Figure 1: region map
    fig_map, ax_map = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6, 5))
    plot_region_map(ax_map)
    fig_map.tight_layout()
    fig_map.savefig(map_output, dpi=100)
    print(f'Saved {map_output}')
    plt.close(fig_map)

    # Figure 2: diurnal cycles
    fig_dc, (ax_cape, ax_lnb) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    fig_dc.suptitle(
        f'Diurnal cycle over Burkina Faso ({BF_LAT_MIN}–{BF_LAT_MAX}°N, '
        f'{BF_LON_MIN}–{BF_LON_MAX}°E)'
    )
    for ax, var, ylabel, invert in [
        (ax_cape, 'cape', 'CAPE (J kg$^{-1}$)', False),
        (ax_lnb,  'lnb',  'LNB pressure (hPa)', True),
    ]:
        for season, color in [('jja', 'tab:orange'), ('djf', 'tab:blue')]:
            mean, std = cycles[f'{var}_{season}']
            ax.plot(hours, mean, label=season.upper(), color=color)
            ax.fill_between(hours, mean - std, mean + std, alpha=0.2, color=color)
        ax.set_xlabel('Hour (UTC)')
        ax.set_ylabel(ylabel)
        ax.set_xticks(hours)
        ax.legend()
        ax.set_title(ylabel.split(' (')[0])
        if invert:
            ax.invert_yaxis()
    fig_dc.savefig(dc_output, dpi=100)
    print(f'Saved {dc_output}')
    plt.close(fig_dc)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--zarr', default='entrainment_wam.zarr', help='Path to zarr store')
    parser.add_argument('--output', default='figs/diurnal_cycle.png', help='Output plot path')
    args = parser.parse_args()

    from pathlib import Path
    Path(args.output).parent.mkdir(exist_ok=True)
    cycles = compute_diurnal_cycles(args.zarr)
    plot(cycles, args.output)


if __name__ == '__main__':
    main()
