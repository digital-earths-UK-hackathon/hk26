"""
Compute per-track entrainment statistics for region-filtered MCS tracks.

Links entrainment_<region>.zarr (3-hourly) with the PyFLEXTRKR MCS pixel mask
(hourly, S3) to produce a NetCDF with dims (tracks, times_3h) following
PyFLEXTRKR output conventions.

Usage:
    python calc_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26
    python calc_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26 --surface land
    python calc_mcs_entrainment.py --model um_glm_n2560_CoMA9_hk26 --region wam
"""
import argparse
import warnings
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

import models

warnings.filterwarnings('ignore', message='.*The return type of `Dataset.dims`.*', category=FutureWarning)

# ceil(650 h / 3 h): max track duration in 3-hourly steps
MAX_TIMES_3H = 217

ENTR_VARS = ['cape', 'cin', 'lnb', 't_lnb', 'w_eff', 'tb', 'tb_diff']

# These are set in main() from --model / --region args before any function uses them.
ZOOM             = None
MASK_URL         = None
STATS_URL        = None
ENTRAINMENT_ZARR = None


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_track_stats():
    """Fetch PyFLEXTRKR track statistics from S3."""
    print('Loading MCS track statistics...')
    response = requests.get(STATS_URL, stream=True)
    dstracks = xr.open_dataset(BytesIO(response.content), mask_and_scale=True)

    # Round times to nearest second (small offset in source data)
    def _round(t):
        return (np.round(t.astype(int) / 1e9) * 1e9).astype('datetime64[ns]')

    for field in ['base_time', 'start_basetime', 'end_basetime']:
        dstracks[field].load()
        tmask = ~np.isnan(dstracks[field].values)
        dstracks[field].values[tmask] = _round(dstracks[field].values[tmask])

    return dstracks


def open_entrainment():
    """Open the local entrainment zarr store."""
    print('Opening entrainment zarr...')
    return xr.open_zarr(ENTRAINMENT_ZARR)


def open_mcs_mask():
    """Open the MCS pixel mask zarr from S3."""
    print('Opening MCS mask zarr...')
    return xr.open_zarr(MASK_URL, chunks={})


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

def compute_wam_positions(entr_ds, mask_ds):
    """
    Compute positional indices of WAM cells within the global mask cell array.
    Returns an int array of shape (n_wam_cells,).
    """
    wam_cells    = entr_ds.cell.values             # HEALPix cell numbers, WAM subset
    global_cells = mask_ds.healpix_index.values    # HEALPix cell numbers, global (0…N-1)

    # global_cells is 0,1,2,...,N-1 so positions == wam_cells, but use searchsorted
    # for correctness in case of non-contiguous ranges.
    # MM: explain searchsorted.
    positions = np.searchsorted(global_cells, wam_cells)
    assert np.all(global_cells[positions] == wam_cells), \
        'WAM cell indices not found in global mask — zoom mismatch?'
    return positions


def filter_region_tracks(dstracks, region_cfg):
    """
    Keep only tracks whose centroid enters the analysis region (with buffer)
    at any point during their lifetime.
    """
    display = region_cfg['display']
    print(f'Filtering tracks to {display} region...')
    dstracks.meanlat.load()
    dstracks.meanlon.load()

    lat   = dstracks.meanlat.values        # (tracks, times)
    lon   = dstracks.meanlon.values        # (tracks, times), in [0, 360]
    lon180 = (lon + 180) % 360 - 180      # convert to [-180, 180]

    in_lat = (lat   > region_cfg['buf_lat_min']) & (lat   < region_cfg['buf_lat_max'])
    in_lon = (lon180 > region_cfg['buf_lon_min']) & (lon180 < region_cfg['buf_lon_max'])
    in_region = (in_lat & in_lon).any(axis=1)

    filtered = dstracks.isel(tracks=in_region)
    print(f'  {int(in_region.sum())} / {dstracks.sizes["tracks"]} tracks pass {display} filter')
    return filtered


LAND_FRAC_THRESHOLD  = 0.8   # mean pf_landfrac above this → land MCS
OCEAN_FRAC_THRESHOLD = 0.2   # mean pf_landfrac below this → ocean MCS


def filter_surface(dstracks, surface):
    """
    Filter tracks by mean land fraction (pf_landfrac, expressed as 0–1).
      'land'  : mean pf_landfrac > 0.8
      'ocean' : mean pf_landfrac < 0.2
      'all'   : no filter (default)
    """
    if surface == 'all':
        return dstracks

    print(f'Filtering tracks by surface type: {surface}...')
    dstracks.pf_landfrac.load()
    mean_lf = np.nanmean(dstracks.pf_landfrac.values, axis=1)  # (tracks,)

    if surface == 'land':
        mask = mean_lf > LAND_FRAC_THRESHOLD
    else:  # ocean
        mask = mean_lf < OCEAN_FRAC_THRESHOLD

    filtered = dstracks.isel(tracks=mask)
    print(f'  {int(mask.sum())} / {dstracks.sizes["tracks"]} tracks pass {surface} filter')
    return filtered


def align_times(entr_ds, mask_ds):
    """
    Find 3-hourly entrainment times that exist in the hourly mask dataset.
    Returns:
        times_3h   : 1-D array of datetime64 values (3-hourly, in overlap)
        mask_indices: corresponding positional indices in mask_ds.time
    """
    entr_times = entr_ds.time.values
    mask_times = mask_ds.time.values

    mask_time_to_idx = {pd.Timestamp(t): i for i, t in enumerate(mask_times)}

    overlap_entr_idx = []
    overlap_mask_idx = []
    for i, t in enumerate(entr_times):
        ts = pd.Timestamp(t)
        if ts in mask_time_to_idx:
            overlap_entr_idx.append(i)
            overlap_mask_idx.append(mask_time_to_idx[ts])

    print(f'Overlap: {len(overlap_entr_idx)} 3-hourly timesteps '
          f'({pd.Timestamp(entr_times[overlap_entr_idx[0]])} – '
          f'{pd.Timestamp(entr_times[overlap_entr_idx[-1]])})')

    return (np.array(overlap_entr_idx),
            np.array(overlap_mask_idx),
            entr_times[overlap_entr_idx])


# ---------------------------------------------------------------------------
# Bincount aggregation
# ---------------------------------------------------------------------------

def bincount_mean_std(mask_int, values, max_label):
    """
    Vectorised per-label mean and std, excluding NaN values in `values`.

    Returns means, stds, counts each of shape (max_label + 1,).
    Index 0 = background (no MCS).
    """
    mask_work = mask_int.copy()
    valid = np.isfinite(values)
    mask_work[~valid] = 0          # exclude NaN cells from aggregation

    vals = np.where(valid, values.astype(np.float64), 0.0)

    minlen = max_label + 1
    # MM: explain bincount.
    counts  = np.bincount(mask_work, minlength=minlen)
    sums    = np.bincount(mask_work, weights=vals,    minlength=minlen)
    sum_sq  = np.bincount(mask_work, weights=vals**2, minlength=minlen)

    with np.errstate(invalid='ignore', divide='ignore'):
        means = np.where(counts > 0, sums / counts, np.nan)
        ex2   = np.where(counts > 0, sum_sq / counts, np.nan)
        stds  = np.sqrt(np.maximum(0.0, ex2 - means**2))

    return means.astype(np.float32), stds.astype(np.float32), counts


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_track_entrainment(entr_ds, mask_ds, dstracks_wam,
                              wam_positions, entr_idxs, mask_idxs, times_3h):
    """
    Main loop: iterate over 3-hourly timesteps, aggregate entrainment per track.

    Returns a dict of output arrays keyed by variable name.
    """
    n_tracks   = dstracks_wam.sizes['tracks']
    max_label  = int(dstracks_wam.tracks.values.max()) + 1  # mask value = track_idx + 1

    # Map mask track number → output row index
    track_nums = dstracks_wam.tracks.values.astype(int)   # original track indices from full dstracks (sparse after WAM filter)
    mask_num_to_out_idx = np.full(max_label + 1, -1, dtype=np.int32)
    for out_i, tn in enumerate(track_nums):
        # MM: Where does the +1 come from?
        mask_num_to_out_idx[tn + 1] = out_i    # mask value = track_idx + 1

    # First 3-hourly step index for each WAM track (index into times_3h array)
    start_times = dstracks_wam.start_basetime.values   # (n_tracks,)
    first_3h_step = np.searchsorted(times_3h, start_times)

    # Allocate output arrays
    nan32 = np.full((n_tracks, MAX_TIMES_3H), np.nan, dtype=np.float32)
    arrays = {}
    for v in ENTR_VARS:
        arrays[f'{v}_mean'] = nan32.copy()
        arrays[f'{v}_std']  = nan32.copy()
    arrays['n_wam_cells'] = np.zeros((n_tracks, MAX_TIMES_3H), dtype=np.int32)

    # base_time: NaT by default, filled during loop
    base_time_out = np.full((n_tracks, MAX_TIMES_3H),
                            np.datetime64('NaT', 'ns'), dtype='datetime64[ns]')

    n_steps = len(entr_idxs)
    for step, (ei, mi, t) in enumerate(zip(entr_idxs, mask_idxs, times_3h)):
        if step % 100 == 0:
            print(f'  {step}/{n_steps}  ({pd.Timestamp(t)})', flush=True)

        # Load mask at this hour → select WAM cells; fill NaN (no MCS) with 0
        mask_global = mask_ds.mcs_mask.isel(time=mi).compute().values
        mask_global = np.nan_to_num(mask_global, nan=0.0)
        mask_wam    = mask_global[wam_positions].astype(np.int32)

        active_nums = np.unique(mask_wam[mask_wam > 0])
        if len(active_nums) == 0:
            continue

        max_active  = int(active_nums.max())

        # Load entrainment at this 3-hourly step
        entr_slice = entr_ds.isel(time=ei)[ENTR_VARS].compute()

        # Aggregate each variable
        var_results = {}
        for v in ENTR_VARS:
            vals = entr_slice[v].values.astype(np.float64)
            means, stds, counts = bincount_mean_std(mask_wam, vals, max_active)
            var_results[v] = (means, stds, counts)

        # Write to output arrays
        for mask_num in active_nums:
            out_i = mask_num_to_out_idx[mask_num] if mask_num <= len(mask_num_to_out_idx) - 1 else -1
            if out_i < 0:
                continue
            li = step - first_3h_step[out_i]
            if li < 0 or li >= MAX_TIMES_3H:
                continue

            for v in ENTR_VARS:
                means, stds, counts = var_results[v]
                if mask_num <= len(means) - 1:
                    arrays[f'{v}_mean'][out_i, li] = means[mask_num]
                    arrays[f'{v}_std'][out_i, li]  = stds[mask_num]
                    arrays['n_wam_cells'][out_i, li] = counts[mask_num]

            base_time_out[out_i, li] = t

    return arrays, base_time_out


# ---------------------------------------------------------------------------
# Build and save output dataset
# ---------------------------------------------------------------------------

def build_output_dataset(arrays, base_time_out, dstracks_wam, times_3h):
    n_tracks = dstracks_wam.sizes['tracks']
    tracks_coord = dstracks_wam.tracks.values
    times_coord  = np.arange(MAX_TIMES_3H)

    data_vars = {}

    # Entrainment mean/std variables
    for v in ENTR_VARS:
        attrs_base = {}
        if v in dstracks_wam.data_vars:
            pass  # could copy attrs from source
        data_vars[f'{v}_mean'] = xr.DataArray(
            arrays[f'{v}_mean'], dims=['tracks', 'times_3h'],
            attrs={'description': f'Mean {v} over MCS cells in WAM region'}
        )
        data_vars[f'{v}_std'] = xr.DataArray(
            arrays[f'{v}_std'], dims=['tracks', 'times_3h'],
            attrs={'description': f'Std {v} over MCS cells in WAM region'}
        )

    data_vars['n_wam_cells'] = xr.DataArray(
        arrays['n_wam_cells'], dims=['tracks', 'times_3h'],
        attrs={'description': 'Number of MCS cells in WAM region at this lifecycle step'}
    )
    data_vars['base_time'] = xr.DataArray(
        base_time_out, dims=['tracks', 'times_3h'],
        attrs={'description': 'UTC time of each 3-hourly lifecycle step (NaT if inactive)'}
    )

    # Track-level metadata from dstracks
    data_vars['track_duration_3h'] = xr.DataArray(
        np.ceil(dstracks_wam.track_duration.values / 3).astype(np.int32),
        dims=['tracks'],
        attrs={'description': 'Track duration in 3-hourly steps', 'units': '3h'}
    )
    data_vars['start_basetime'] = xr.DataArray(
        dstracks_wam.start_basetime.values, dims=['tracks'],
        attrs={'description': 'Track start time'}
    )
    data_vars['end_basetime'] = xr.DataArray(
        dstracks_wam.end_basetime.values, dims=['tracks'],
        attrs={'description': 'Track end time'}
    )
    data_vars['meanlat'] = xr.DataArray(
        dstracks_wam.meanlat.values[:, :MAX_TIMES_3H * 3:3],  # downsample 1h→3h
        dims=['tracks', 'times_3h'],
        attrs={'description': 'MCS centroid latitude (hourly, downsampled to 3-hourly)'}
    )
    data_vars['meanlon'] = xr.DataArray(
        dstracks_wam.meanlon.values[:, :MAX_TIMES_3H * 3:3],
        dims=['tracks', 'times_3h'],
        attrs={'description': 'MCS centroid longitude (hourly, downsampled to 3-hourly)'}
    )

    ds_out = xr.Dataset(
        data_vars,
        coords={'tracks': tracks_coord, 'times_3h': times_coord},
        attrs={
            'description': 'Per-track entrainment statistics for WAM-region MCS',
            'entrainment_source': str(ENTRAINMENT_ZARR),
            'mask_source': MASK_URL,
            'stats_source': STATS_URL,
        }
    )
    return ds_out


def save_output(ds_out, output_path):
    encoding = {v: {'zlib': True, 'complevel': 4}
                for v in ds_out.data_vars if v not in ('base_time', 'start_basetime', 'end_basetime')}
    ds_out.to_netcdf(output_path, encoding=encoding)
    print(f'Saved {output_path}')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    models.add_model_arg(parser)
    models.add_region_arg(parser)
    parser.add_argument('--output', default=None,
                        help='Output NetCDF path (default: data/<model>/mcs_entrainment_<region>[_<surface>].nc)')
    parser.add_argument('--n-timesteps', type=int, default=None, metavar='N',
                        help='Limit to first N timesteps (for testing)')
    parser.add_argument('--surface', choices=['all', 'land', 'ocean'], default='all',
                        help='Filter MCS by mean land fraction: land (>0.8), ocean (<0.2), all (default)')
    args = parser.parse_args()

    region_cfg = models.REGIONS[args.region]

    if args.output is None:
        suffix = f'_{args.surface}' if args.surface != 'all' else ''
        args.output = str(
            models.data_dir(args.model) / f'mcs_entrainment_{args.region}{suffix}.nc'
        )

    # Patch module-level URL/path constants to match the chosen model/region.
    global ENTRAINMENT_ZARR, MASK_URL, STATS_URL, ZOOM
    ZOOM             = models.MODELS[args.model]['zoom']
    MASK_URL         = models.mask_url(args.model)
    STATS_URL        = models.stats_url(args.model)
    ENTRAINMENT_ZARR = models.data_dir(args.model) / f'entrainment_{args.region}.zarr'

    dstracks         = load_track_stats()
    entr_ds          = open_entrainment()
    mask_ds          = open_mcs_mask()

    wam_positions    = compute_wam_positions(entr_ds, mask_ds)
    dstracks_wam     = filter_region_tracks(dstracks, region_cfg)
    dstracks_wam     = filter_surface(dstracks_wam, args.surface)

    entr_idxs, mask_idxs, times_3h = align_times(entr_ds, mask_ds)

    if args.n_timesteps is not None:
        entr_idxs = entr_idxs[:args.n_timesteps]
        mask_idxs = mask_idxs[:args.n_timesteps]
        times_3h  = times_3h[:args.n_timesteps]

    print(f'Processing {len(entr_idxs)} timesteps for '
          f'{dstracks_wam.sizes["tracks"]} WAM tracks...')

    arrays, base_time_out = compute_track_entrainment(
        entr_ds, mask_ds, dstracks_wam,
        wam_positions, entr_idxs, mask_idxs, times_3h
    )

    ds_out = build_output_dataset(arrays, base_time_out, dstracks_wam, times_3h)
    save_output(ds_out, args.output)


if __name__ == '__main__':
    main()
