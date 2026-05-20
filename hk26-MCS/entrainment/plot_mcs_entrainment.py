"""
Plot entrainment statistics linked to MCS tracks.

Two figures:
  1. Lifecycle composite — each entrainment variable vs normalised lifecycle
     fraction (0 = start, 1 = end), split by JJA / DJF and duration category.
  2. MCS entrainment diurnal cycle — compare MCS-cell entrainment vs the
     all-cell background diurnal cycle, split by JJA / DJF.

Usage:
    python plot_mcs_entrainment.py
    python plot_mcs_entrainment.py --input mcs_entrainment_wam.nc --output figs/
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

ENTR_VARS = ['cape', 'cin', 'lnb', 'w_eff', 'tb_diff']
VAR_LABELS = {
    'cape':    'CAPE (J kg$^{-1}$)',
    'cin':     'CIN (J kg$^{-1}$)',
    'lnb':     'LNB pressure (hPa)',
    'w_eff':   'w / $\\sqrt{\\mathrm{CAPE}}$ (m s$^{-1}$ (J kg$^{-1}$)$^{-0.5}$)',
    'tb_diff': '$T_b - T_{\\mathrm{LNB}}$ (K)',
}
VAR_INVERT = {'lnb': True}   # invert y-axis for pressure

N_LIFECYCLE_BINS = 20        # number of lifecycle fraction bins


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_mcs_entrainment(nc_path):
    return xr.open_dataset(nc_path)


def season_mask(times):
    """Return boolean arrays for JJA and DJF from a 1-D array of datetime64."""
    months = pd.DatetimeIndex(times).month
    return (np.isin(months, [6, 7, 8]),
            np.isin(months, [12, 1, 2]))


# ---------------------------------------------------------------------------
# Figure 1: lifecycle composites
# ---------------------------------------------------------------------------

def compute_lifecycle_composite(ds, var, season_bool=None):
    """
    Composite mean and IQR of `var_mean` as a function of lifecycle fraction.
    Only uses timesteps where n_wam_cells > 0.

    Returns (bin_centres, composite_median, composite_q25, composite_q75).
    """
    mean_vals = ds[f'{var}_mean'].values          # (tracks, times_3h)
    n_cells   = ds['n_wam_cells'].values           # (tracks, times_3h)
    dur       = ds['track_duration_3h'].values     # (tracks,)

    bins  = np.linspace(0, 1, N_LIFECYCLE_BINS + 1)
    mids  = (bins[:-1] + bins[1:]) / 2

    all_vals = [[] for _ in range(N_LIFECYCLE_BINS)]

    for ti in range(ds.sizes['tracks']):
        d = int(dur[ti])
        if d <= 0:
            continue
        for li in range(min(d, ds.sizes['times_3h'])):
            if n_cells[ti, li] == 0:
                continue
            v = mean_vals[ti, li]
            if not np.isfinite(v):
                continue
            # Lifecycle fraction at the MIDPOINT of this step
            frac = (li + 0.5) / d

            if season_bool is not None and not season_bool[ti]:
                continue

            bin_idx = int(frac * N_LIFECYCLE_BINS)
            bin_idx = min(bin_idx, N_LIFECYCLE_BINS - 1)
            all_vals[bin_idx].append(v)

    medians = np.array([np.nanmedian(b) if b else np.nan for b in all_vals])
    q25     = np.array([np.nanpercentile(b, 25) if b else np.nan for b in all_vals])
    q75     = np.array([np.nanpercentile(b, 75) if b else np.nan for b in all_vals])

    return mids, medians, q25, q75


def get_track_season_bool(ds, season='jja'):
    """Return boolean array (n_tracks,) indicating JJA or DJF tracks by start time."""
    starts = pd.DatetimeIndex(ds['start_basetime'].values)
    months = starts.month
    if season == 'jja':
        return np.isin(months, [6, 7, 8])
    return np.isin(months, [12, 1, 2])


def plot_lifecycle(ds, output_dir):
    jja_mask = get_track_season_bool(ds, 'jja')
    djf_mask = get_track_season_bool(ds, 'djf')

    n_vars = len(ENTR_VARS)
    fig, axes = plt.subplots(1, n_vars, figsize=(4 * n_vars, 4), layout='constrained')
    fig.suptitle('MCS entrainment lifecycle composite (WAM region)')

    for ax, var in zip(axes, ENTR_VARS):
        for mask, label, color in [
            (jja_mask, 'JJA', 'tab:orange'),
            (djf_mask, 'DJF', 'tab:blue'),
        ]:
            mids, med, q25, q75 = compute_lifecycle_composite(ds, var, mask)
            ax.plot(mids, med, color=color, label=label)
            ax.fill_between(mids, q25, q75, alpha=0.2, color=color)

        ax.set_xlabel('Lifecycle fraction')
        ax.set_ylabel(VAR_LABELS[var])
        ax.set_title(var)
        ax.legend(fontsize=8)
        if VAR_INVERT.get(var):
            ax.invert_yaxis()
        if var in ('cin', 'tb_diff', 'w_eff'):
            ax.axhline(0, color='k', linewidth=0.8, linestyle='--')

    out = Path(output_dir) / 'mcs_lifecycle_entrainment.png'
    fig.savefig(out, dpi=100)
    print(f'Saved {out}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: MCS entrainment diurnal cycle
# ---------------------------------------------------------------------------

def compute_mcs_diurnal_cycle(ds, var):
    """
    Group per-track, per-step mean values by UTC hour.
    Only includes steps with n_wam_cells > 0.
    Returns (hours, jja_mean, jja_std, djf_mean, djf_std).
    """
    mean_vals   = ds[f'{var}_mean'].values      # (tracks, times_3h)
    n_cells     = ds['n_wam_cells'].values
    base_time   = ds['base_time'].values         # (tracks, times_3h), datetime64

    hours = np.arange(0, 24, 3)
    jja_by_hour = {h: [] for h in hours}
    djf_by_hour = {h: [] for h in hours}

    jja_bool = get_track_season_bool(ds, 'jja')
    djf_bool = get_track_season_bool(ds, 'djf')

    for ti in range(ds.sizes['tracks']):
        is_jja = jja_bool[ti]
        is_djf = djf_bool[ti]
        if not (is_jja or is_djf):
            continue

        for li in range(ds.sizes['times_3h']):
            if n_cells[ti, li] == 0:
                continue
            v = mean_vals[ti, li]
            if not np.isfinite(v):
                continue
            t = base_time[ti, li]
            if t == np.datetime64('NaT'):
                continue
            h = pd.Timestamp(t).hour

            if is_jja:
                jja_by_hour[h].append(v)
            if is_djf:
                djf_by_hour[h].append(v)

    def agg(by_hour):
        means = np.array([np.nanmean(by_hour[h]) if by_hour[h] else np.nan for h in hours])
        stds  = np.array([np.nanstd(by_hour[h])  if by_hour[h] else np.nan for h in hours])
        return means, stds

    return hours, *agg(jja_by_hour), *agg(djf_by_hour)


def plot_mcs_diurnal_cycle(ds, output_dir):
    n_vars = len(ENTR_VARS)
    fig, axes = plt.subplots(1, n_vars, figsize=(4 * n_vars, 4), layout='constrained')
    fig.suptitle('MCS entrainment diurnal cycle (WAM region, MCS cells only)')

    for ax, var in zip(axes, ENTR_VARS):
        hours, jja_m, jja_s, djf_m, djf_s = compute_mcs_diurnal_cycle(ds, var)
        for mean, std, label, color in [
            (jja_m, jja_s, 'JJA', 'tab:orange'),
            (djf_m, djf_s, 'DJF', 'tab:blue'),
        ]:
            ax.plot(hours, mean, color=color, label=label)
            ax.fill_between(hours, mean - std, mean + std, alpha=0.2, color=color)

        ax.set_xlabel('Hour (UTC)')
        ax.set_ylabel(VAR_LABELS[var])
        ax.set_title(var)
        ax.set_xticks(hours)
        ax.legend(fontsize=8)
        if VAR_INVERT.get(var):
            ax.invert_yaxis()
        if var in ('cin', 'tb_diff', 'w_eff'):
            ax.axhline(0, color='k', linewidth=0.8, linestyle='--')

    out = Path(output_dir) / 'mcs_dc_entrainment.png'
    fig.savefig(out, dpi=100)
    print(f'Saved {out}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',  default='mcs_entrainment_wam.nc', help='Input NetCDF')
    parser.add_argument('--output', default='figs',                   help='Output directory')
    args = parser.parse_args()

    Path(args.output).mkdir(exist_ok=True)

    print(f'Loading {args.input}...')
    ds = load_mcs_entrainment(args.input)
    print(ds)

    plot_lifecycle(ds, args.output)
    plot_mcs_diurnal_cycle(ds, args.output)


if __name__ == '__main__':
    main()
