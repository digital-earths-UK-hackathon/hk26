"""
Plot entrainment statistics linked to MCS tracks.

Three figures:
  1. Lifecycle composite — each entrainment variable vs normalised lifecycle
     fraction (0 = start, 1 = end), split by JJA / DJF and duration category.
  2. MCS entrainment diurnal cycle — compare MCS-cell entrainment vs the
     all-cell background diurnal cycle, split by JJA / DJF.
  3. Distributions — histograms of all key variables (JJA vs DJF).

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

import models

MODEL = 'UM N2560 RAL3.3'   # overridden at runtime if --model is given

ENTR_VARS = ['cape', 'cin', 'lnb', 'tb', 'w_eff', 'tb_diff']
VAR_LABELS = {
    'cape':    'CAPE (J kg$^{-1}$)',
    'cin':     'CIN (J kg$^{-1}$)',
    'lnb':     'LNB pressure (hPa)',
    'w_eff':   'w / $\\sqrt{\\mathrm{CAPE}}$ (-)',
    'tb':      '$T_b$ (K)',
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
    Composite median and IQR of `var_mean` as a function of lifecycle fraction,
    giving each track equal weight.

    For each track, valid (fraction, value) pairs are linearly interpolated onto
    the fixed bin-midpoint grid; points outside the track's observed range are
    NaN.  The median and IQR are then taken across tracks at each grid point.

    Returns (bin_centres, composite_median, composite_q25, composite_q75).
    """
    mean_vals = ds[f'{var}_mean'].values       # (tracks, times_3h)
    n_cells   = ds['n_wam_cells'].values       # (tracks, times_3h)
    dur       = ds['track_duration_3h'].values # (tracks,)

    mids = (np.linspace(0, 1, N_LIFECYCLE_BINS + 1)[:-1]
            + np.linspace(0, 1, N_LIFECYCLE_BINS + 1)[1:]) / 2

    track_interps = []

    for ti in range(ds.sizes['tracks']):
        if season_bool is not None and not season_bool[ti]:
            continue
        d = int(dur[ti])
        if d <= 0:
            continue

        fracs, vals = [], []
        for li in range(min(d, ds.sizes['times_3h'])):
            if n_cells[ti, li] == 0:
                continue
            v = mean_vals[ti, li]
            if not np.isfinite(v):
                continue
            fracs.append((li + 0.5) / d)
            vals.append(v)

        if len(fracs) < 2:
            continue

        # Interpolate onto fixed grid; NaN outside the track's observed range.
        interp = np.interp(mids, fracs, vals, left=np.nan, right=np.nan)
        track_interps.append(interp)

    nans = np.full(N_LIFECYCLE_BINS, np.nan)
    if not track_interps:
        return mids, nans, nans, nans

    arr     = np.array(track_interps)           # (n_tracks, N_LIFECYCLE_BINS)
    medians = np.nanmedian(arr, axis=0)
    q25     = np.nanpercentile(arr, 25, axis=0)
    q75     = np.nanpercentile(arr, 75, axis=0)

    return mids, medians, q25, q75


def get_track_season_bool(ds, season='jja'):
    """Return boolean array (n_tracks,) indicating JJA or DJF tracks by start time."""
    starts = pd.DatetimeIndex(ds['start_basetime'].values)
    months = starts.month
    if season == 'jja':
        return np.isin(months, [6, 7, 8])
    return np.isin(months, [12, 1, 2])


def plot_lifecycle(ds, output_dir, stem):
    jja_mask = get_track_season_bool(ds, 'jja')
    djf_mask = get_track_season_bool(ds, 'djf')

    n_vars = len(ENTR_VARS)
    fig, axes = plt.subplots(1, n_vars, figsize=(4 * n_vars, 4), layout='constrained')
    fig.suptitle(f'MCS entrainment lifecycle composite (WAM region) — {MODEL}')

    for ax, var in zip(axes, ENTR_VARS):
        mean_line_vals = []
        for mask, label, color in [
            (jja_mask, f'JJA (N={jja_mask.sum()})', 'tab:orange'),
            (djf_mask, f'DJF (N={djf_mask.sum()})', 'tab:blue'),
        ]:
            mids, med, q25, q75 = compute_lifecycle_composite(ds, var, mask)
            ax.plot(mids, med, color=color, label=label)
            ax.fill_between(mids, q25, q75, alpha=0.2, color=color)
            mean_line_vals.extend(med[np.isfinite(med)])

        ax.set_xlabel('Lifecycle fraction')
        ax.set_ylabel(VAR_LABELS[var])
        ax.set_title(var)
        ax.legend(fontsize=8)
        if VAR_INVERT.get(var):
            ax.invert_yaxis()
        if var in ('cin', 'tb_diff', 'w_eff'):
            ax.axhline(0, color='k', linewidth=0.8, linestyle='--')
        if var == 'w_eff' and mean_line_vals:
            lo, hi = min(mean_line_vals), max(mean_line_vals)
            margin = max((hi - lo) * 0.15, 1e-4)
            ax.set_ylim(lo - margin, hi + margin)

    out = Path(output_dir) / f'{stem}.mcs_lifecycle_entrainment.png'
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


def plot_mcs_diurnal_cycle(ds, output_dir, stem):
    jja_mask = get_track_season_bool(ds, 'jja')
    djf_mask = get_track_season_bool(ds, 'djf')

    n_vars = len(ENTR_VARS)
    fig, axes = plt.subplots(1, n_vars, figsize=(4 * n_vars, 4), layout='constrained')
    fig.suptitle(f'MCS entrainment diurnal cycle (WAM region, MCS cells only) — {MODEL}')

    for ax, var in zip(axes, ENTR_VARS):
        hours, jja_m, jja_s, djf_m, djf_s = compute_mcs_diurnal_cycle(ds, var)
        mean_line_vals = []
        for mean, std, label, color in [
            (jja_m, jja_s, f'JJA (N={jja_mask.sum()})', 'tab:orange'),
            (djf_m, djf_s, f'DJF (N={djf_mask.sum()})', 'tab:blue'),
        ]:
            ax.plot(hours, mean, color=color, label=label)
            ax.fill_between(hours, mean - std, mean + std, alpha=0.2, color=color)
            mean_line_vals.extend(mean[np.isfinite(mean)])

        ax.set_xlabel('Hour (UTC)')
        ax.set_ylabel(VAR_LABELS[var])
        ax.set_title(var)
        ax.set_xticks(hours)
        ax.legend(fontsize=8)
        if VAR_INVERT.get(var):
            ax.invert_yaxis()
        if var in ('cin', 'tb_diff', 'w_eff'):
            ax.axhline(0, color='k', linewidth=0.8, linestyle='--')
        if var == 'w_eff' and mean_line_vals:
            lo, hi = min(mean_line_vals), max(mean_line_vals)
            margin = max((hi - lo) * 0.15, 1e-4)
            ax.set_ylim(lo - margin, hi + margin)

    out = Path(output_dir) / f'{stem}.mcs_dc_entrainment.png'
    fig.savefig(out, dpi=100)
    print(f'Saved {out}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: variable distributions
# ---------------------------------------------------------------------------

DIST_VARS = ['cape_mean', 'cin_mean', 'lnb_mean', 't_lnb_mean', 'w_eff_mean', 'tb_mean', 'tb_diff_mean']
DIST_LABELS = {
    'cape_mean':    'CAPE (J kg$^{-1}$)',
    'cin_mean':     'CIN (J kg$^{-1}$)',
    'lnb_mean':     'LNB pressure (hPa)',
    't_lnb_mean':   '$T_{\\mathrm{LNB}}$ (K)',
    'w_eff_mean':   'w / $\\sqrt{\\mathrm{CAPE}}$ (-)',
    'tb_mean':      '$T_b$ (K)',
    'tb_diff_mean': '$T_b - T_{\\mathrm{LNB}}$ (K)',
}


def _valid_flat(vals, n_valid, season_bool):
    """Return 1-D array of vals where n_valid > 0, finite, and season matches."""
    v = vals[season_bool]       # (n_season_tracks, times_3h)
    n = n_valid[season_bool]
    ok = (n > 0) & np.isfinite(v)
    return v[ok]


def plot_distributions(ds, output_dir, stem):
    jja_mask = get_track_season_bool(ds, 'jja')
    djf_mask = get_track_season_bool(ds, 'djf')
    n_valid = ds['n_wam_cells'].values     # (tracks, times_3h)

    ncols = 4
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 7), layout='constrained')
    axes_flat = axes.ravel()
    fig.suptitle(f'MCS entrainment distributions (WAM region, MCS cells only) — {MODEL}')

    for ax, vname in zip(axes_flat, DIST_VARS):
        vals = ds[vname].values            # (tracks, times_3h)

        # Determine bin range from full distribution (p1–p99 to exclude outliers).
        all_flat = vals[(n_valid > 0) & np.isfinite(vals)]
        if len(all_flat) == 0:
            ax.set_visible(False)
            continue
        lo, hi = np.percentile(all_flat, [1, 99])
        bins = np.linspace(lo, hi, 51)

        for season_bool, label, color in [
            (jja_mask, f'JJA (N={jja_mask.sum()})', 'tab:orange'),
            (djf_mask, f'DJF (N={djf_mask.sum()})', 'tab:blue'),
        ]:
            v_flat = _valid_flat(vals, n_valid, season_bool)
            ax.hist(v_flat, bins=bins, alpha=0.5, density=True, color=color, label=label)

        ax.set_xlabel(DIST_LABELS[vname])
        ax.set_ylabel('Density')

        # Reference lines and annotations.
        if vname == 'tb_mean':
            ax.axvline(241, color='k', linestyle='--', linewidth=1, label='241 K')
            n_above = (all_flat > 241).sum()
            ax.text(0.97, 0.97, f'{100 * n_above / len(all_flat):.1f}% > 241 K',
                    transform=ax.transAxes, ha='right', va='top', fontsize=7)
        elif vname == 'tb_diff_mean':
            ax.axvline(0, color='k', linestyle='--', linewidth=0.8)
            ax.text(0.03, 0.97, '← overshooting', transform=ax.transAxes,
                    ha='left', va='top', fontsize=7)
        elif vname in ('cin_mean', 'w_eff_mean'):
            ax.axvline(0, color='k', linestyle='--', linewidth=0.8)

        # Invert x-axis for pressure (higher pressure = lower altitude).
        if vname == 'lnb_mean':
            ax.invert_xaxis()

        ax.legend(fontsize=8)

    # Hide unused axes (8 panels, 7 variables).
    for ax in axes_flat[len(DIST_VARS):]:
        ax.set_visible(False)

    out = Path(output_dir) / f'{stem}.mcs_distributions.png'
    fig.savefig(out, dpi=100)
    print(f'Saved {out}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figures 4 & 5: x-variable vs entrainment proxies (hexbin)
# ---------------------------------------------------------------------------

ENTR_PROXY_VARS = [
    ('w_eff_mean',   VAR_LABELS['w_eff'],   False, (-0.2, 0.5)),
    ('tb_diff_mean', VAR_LABELS['tb_diff'],  True,  None),
]


def _plot_x_vs_entr(ds, x_vals, x_label, title, output_path, include_line=False):
    jja_mask = get_track_season_bool(ds, 'jja')
    djf_mask = get_track_season_bool(ds, 'djf')
    n_valid  = ds['n_wam_cells'].values

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')
    fig.suptitle(f'{title} — {MODEL}')

    for col, (vname, ylabel, zero_line, ylim) in enumerate(ENTR_PROXY_VARS):
        vals = ds[vname].values
        for row, (season_bool, season_label) in enumerate([
            (jja_mask, 'JJA'), (djf_mask, 'DJF')
        ]):
            ax = axes[row, col]

            ok = season_bool[:, np.newaxis] & (n_valid > 0) & np.isfinite(x_vals) & np.isfinite(vals)
            x  = x_vals[ok]
            y  = vals[ok]

            extent = None
            if ylim is not None:
                extent = [x.min(), x.max(), ylim[0], ylim[1]]
            hb = ax.hexbin(x, y, gridsize=30, cmap='YlOrRd', mincnt=1, bins='log',
                           extent=extent)
            plt.colorbar(hb, ax=ax, label='log10(count)')

            if include_line:
                x_bins = np.linspace(x.min(), x.max(), 31)
                x_mids = (x_bins[:-1] + x_bins[1:]) / 2
                bin_idx = np.clip(np.digitize(x, x_bins) - 1, 0, len(x_mids) - 1)
                means = np.array([y[bin_idx == i].mean() if (bin_idx == i).any() else np.nan
                                  for i in range(len(x_mids))])
                ax.plot(x_mids, means, color='steelblue', linewidth=1.5, label='mean')
                ax.legend(fontsize=8)

            if zero_line:
                ax.axhline(0, color='k', linewidth=0.8, linestyle='--')
            if ylim is not None:
                ax.set_ylim(ylim)

            ax.set_xlabel(x_label)
            ax.set_ylabel(ylabel)
            ax.set_title(f'{season_label} (N={season_bool.sum()})')

    fig.savefig(output_path, dpi=100)
    print(f'Saved {output_path}')
    plt.close(fig)


def plot_rh_scatter(ds, output_dir, stem, include_line=False):
    if 'hur700_mean' not in ds:
        print('hur700_mean not in dataset — skipping RH scatter')
        return
    _plot_x_vs_entr(
        ds,
        x_vals=ds['hur700_mean'].values,
        x_label='RH at 700 hPa (%)',
        title='Entrainment proxies vs RH at 700 hPa',
        output_path=Path(output_dir) / f'{stem}.mcs_rh_scatter.png',
        include_line=include_line,
    )


def plot_shear_scatter(ds, output_dir, stem, include_line=False):
    if 'shear_mean' not in ds:
        print('shear_mean not in dataset — skipping shear scatter')
        return
    _plot_x_vs_entr(
        ds,
        x_vals=ds['shear_mean'].values,
        x_label='Zonal shear u(600) − u(850 hPa) (m s$^{-1}$)',
        title='Entrainment proxies vs zonal wind shear',
        output_path=Path(output_dir) / f'{stem}.mcs_shear_scatter.png',
        include_line=include_line,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    global MODEL

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',  default=None, choices=list(models.MODELS),
                        help='Model key; sets default --input, --output, and plot title')
    parser.add_argument('--region', default='wam', choices=list(models.REGIONS),
                        help='Analysis region (default: wam)')
    parser.add_argument('--surface', default='all', choices=['all', 'land', 'ocean'],
                        help='Surface type suffix for default input filename (default: all)')
    parser.add_argument('--include-line', action='store_true', default=False,
                        help='Overlay mean line on RH/shear scatter plots')
    parser.add_argument('--input',  default=None, help='Input NetCDF (overrides --model default)')
    parser.add_argument('--output', default=None, help='Output directory (overrides --model default)')
    args = parser.parse_args()

    if args.model:
        MODEL = models.MODELS[args.model]['display']
        suffix = f'_{args.surface}' if args.surface != 'all' else ''
        nc_input   = args.input  or str(models.data_dir(args.model) / f'mcs_entrainment_{args.region}{suffix}.nc')
        output_dir = args.output or str(models.figs_dir(args.model))
    else:
        nc_input   = args.input  or f'mcs_entrainment_{args.region}.nc'
        output_dir = args.output or 'figs'

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    stem = Path(nc_input).stem

    print(f'Loading {nc_input}...')
    ds = load_mcs_entrainment(nc_input)
    print(ds)

    plot_lifecycle(ds, output_dir, stem)
    plot_mcs_diurnal_cycle(ds, output_dir, stem)
    plot_distributions(ds, output_dir, stem)
    plot_rh_scatter(ds, output_dir, stem, include_line=args.include_line)
    plot_shear_scatter(ds, output_dir, stem, include_line=args.include_line)


if __name__ == '__main__':
    main()
