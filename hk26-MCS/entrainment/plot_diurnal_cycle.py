"""
Plot the diurnal cycle of CAPE, LNB, w_eff and Tb_diff over a sub-box,
split by JJA and DJF seasons, alongside a map of the analysis region.

Usage:
    python plot_diurnal_cycle.py --model um_glm_n2560_RAL3p3_tuned_hk26
    python plot_diurnal_cycle.py --zarr data/um_glm_n2560_RAL3p3_tuned_hk26/entrainment_wam.zarr --output figs/
"""
import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import xarray as xr

import models

# Burkina Faso sub-box (fixed; within the WAM region), in [-180, 180] convention
BF_LAT_MIN, BF_LAT_MAX = 9, 15
BF_LON_MIN, BF_LON_MAX = -5, 5


def load_box(zarr_path):
    """Load zarr store and select cells in the Burkina Faso box."""
    ds     = xr.open_zarr(zarr_path)
    lon180 = (ds.lon + 180) % 360 - 180
    mask   = (
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
    ds_mean   = ds_box[['cape', 'lnb', 'w_eff', 'tb_diff']].mean(dim='cell').compute()

    jja = ds_mean.time.dt.month.isin([6, 7, 8])
    djf = ds_mean.time.dt.month.isin([12, 1, 2])

    cycles = {}
    for var in ['cape', 'lnb', 'w_eff', 'tb_diff']:
        for season, mask in [('jja', jja), ('djf', djf)]:
            cycles[f'{var}_{season}'] = diurnal_cycle(ds_mean[var].isel(time=mask))
    return cycles


def plot_region_map(ax, region_cfg, model_display, region_display):
    """Map showing the analysis region and Burkina Faso sub-box."""
    # Convert WAM lon to [-180, 180] for the map extent
    wam_lon_min = region_cfg['lon_min'] - 360 if region_cfg['lon_min'] > 180 else region_cfg['lon_min']
    wam_lon_max = region_cfg['lon_max']
    wam_lat_min = region_cfg['lat_min']
    wam_lat_max = region_cfg['lat_max']

    ax.set_extent([wam_lon_min - 5, wam_lon_max + 5, wam_lat_min - 3, wam_lat_max + 3],
                  crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND, facecolor='lightgrey', alpha=0.4)

    crs = ccrs.PlateCarree()
    ax.add_patch(mpatches.Rectangle(
        (wam_lon_min, wam_lat_min),
        wam_lon_max - wam_lon_min, wam_lat_max - wam_lat_min,
        linewidth=1.5, edgecolor='steelblue', facecolor='steelblue',
        alpha=0.15, transform=crs, label=f'{region_display} region',
    ))
    ax.add_patch(mpatches.Rectangle(
        (BF_LON_MIN, BF_LAT_MIN),
        BF_LON_MAX - BF_LON_MIN, BF_LAT_MAX - BF_LAT_MIN,
        linewidth=2, edgecolor='tab:orange', facecolor='none',
        transform=crs, label='Burkina Faso box',
    ))
    ax.legend(loc='lower right', fontsize=8)
    ax.set_title(f'Analysis regions — {model_display}')


def plot(cycles, output_dir, stem, model_display, region_display):
    hours = cycles['cape_jja'][0].hour.values
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    def out(suffix):
        return str(Path(output_dir) / f'{stem}{suffix}.png')

    # Figure 1: region map — needs region_cfg, handled by caller passing ax
    # (map is plotted in plot_all)

    # Figure 2: diurnal cycles of CAPE and LNB
    fig_dc, (ax_cape, ax_lnb) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    fig_dc.suptitle(
        f'Diurnal cycle — Burkina Faso ({BF_LAT_MIN}–{BF_LAT_MAX}°N, '
        f'{BF_LON_MIN}–{BF_LON_MAX}°E) — {model_display}'
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
    fig_dc.savefig(out('_dc'), dpi=100)
    print(f'Saved {out("_dc")}')
    plt.close(fig_dc)

    # Figure 3: entrainment proxies
    fig_entr, (ax_weff, ax_tbdiff) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    fig_entr.suptitle(
        f'Entrainment proxies — Burkina Faso ({BF_LAT_MIN}–{BF_LAT_MAX}°N, '
        f'{BF_LON_MIN}–{BF_LON_MAX}°E) — {model_display}'
    )
    for ax, var, ylabel in [
        (ax_weff,   'w_eff',   'w / $\\sqrt{\\mathrm{CAPE}}$ (-)'),
        (ax_tbdiff, 'tb_diff', '$T_b - T_{\\mathrm{LNB}}$ (K)'),
    ]:
        for season, color in [('jja', 'tab:orange'), ('djf', 'tab:blue')]:
            mean, std = cycles[f'{var}_{season}']
            ax.plot(hours, mean, label=season.upper(), color=color)
            ax.fill_between(hours, mean - std, mean + std, alpha=0.2, color=color)
        ax.axhline(0, color='k', linewidth=0.8, linestyle='--')
        ax.set_xlabel('Hour (UTC)')
        ax.set_ylabel(ylabel)
        ax.set_xticks(hours)
        ax.legend()
        ax.set_title(ylabel.split(' (')[0])
    fig_entr.savefig(out('_entr'), dpi=100)
    print(f'Saved {out("_entr")}')
    plt.close(fig_entr)


def plot_map(region_cfg, model_display, region_display, output_dir, stem):
    fig_map, ax_map = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6, 5))
    plot_region_map(ax_map, region_cfg, model_display, region_display)
    fig_map.tight_layout()
    path = str(Path(output_dir) / f'{stem}_map.png')
    fig_map.savefig(path, dpi=100)
    print(f'Saved {path}')
    plt.close(fig_map)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--model',  default=None, choices=list(models.MODELS),
                        help='Model key; sets default --zarr, --output, and plot titles')
    parser.add_argument('--region', default='wam', choices=list(models.REGIONS),
                        help='Analysis region (default: wam)')
    parser.add_argument('--zarr',   default=None, help='Path to entrainment zarr store')
    parser.add_argument('--output', default=None, help='Output directory')
    args = parser.parse_args()

    region_cfg     = models.REGIONS[args.region]
    region_display = region_cfg['display']

    if args.model:
        model_display = models.MODELS[args.model]['display']
        zarr_path     = args.zarr   or str(models.data_dir(args.model) / f'entrainment_{args.region}.zarr')
        output_dir    = args.output or str(models.figs_dir(args.model))
        stem          = f'diurnal_cycle_{args.model}_{args.region}'
    else:
        model_display = ''
        zarr_path     = args.zarr   or f'entrainment_{args.region}.zarr'
        output_dir    = args.output or 'figs'
        stem          = f'diurnal_cycle_{args.region}'

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cycles = compute_diurnal_cycles(zarr_path)
    plot_map(region_cfg, model_display, region_display, output_dir, stem)
    plot(cycles, output_dir, stem, model_display, region_display)


if __name__ == '__main__':
    main()
