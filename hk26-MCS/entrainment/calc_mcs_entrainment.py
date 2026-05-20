"""
Compute per-track entrainment statistics for WAM-region MCS tracks.

Links entrainment_wam.zarr (3-hourly) with the PyFLEXTRKR MCS pixel mask
(hourly, S3) to produce a NetCDF with dims (tracks, times_3h) following
PyFLEXTRKR output conventions.

Usage:
    python calc_mcs_entrainment.py
    python calc_mcs_entrainment.py --output mcs_entrainment_wam.nc
"""
import argparse
import warnings
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

warnings.filterwarnings('ignore', message='.*The return type of `Dataset.dims`.*', category=FutureWarning)

# --- Config ---
ENTRAINMENT_ZARR = Path('entrainment_wam.zarr')
OUTPUT_NC = Path('mcs_entrainment_wam.nc')
ZOOM = 9
MAX_TIMES_3H = 217          # ceil(650 h / 3 h): steps 0,3,6,...,648

MASK_URL = (
    f'https://hackathon-o.s3-ext.jc.rl.ac.uk/sim-data/analysis/PyFLEXTRKR/'
    f'um_glm_n2560_RAL3p3_tuned_z9/mcstracking/'
    f'mcs_mask_hp{ZOOM}_20200201.0000_20210301.0000.zarr'
)
STATS_URL = (
    'https://hackathon-o.s3-ext.jc.rl.ac.uk/sim-data/analysis/PyFLEXTRKR/'
    'um_glm_n2560_RAL3p3_tuned_z9/stats/'
    'mcs_tracks_final_20200201.0000_20210301.0000.nc'
)

ENTR_VARS = ['cape', 'cin', 'lnb', 't_lnb', 'w_eff', 'tb', 'tb_diff']

# WAM pre-filter bounds with 5° buffer (lon in [-180, 180] convention)
WAM_LAT_MIN_BUF, WAM_LAT_MAX_BUF = -3, 20
WAM_LON_MIN_BUF, WAM_LON_MAX_BUF = -25, 25


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
    positions = np.searchsorted(global_cells, wam_cells)
    assert np.all(global_cells[positions] == wam_cells), \
        'WAM cell indices not found in global mask — zoom mismatch?'
    return positions


def filter_wam_tracks(dstracks):
    """
    Keep only tracks whose centroid enters the WAM region (with buffer)
    at any point during their lifetime.
    """
    print('Filtering tracks to WAM region...')
    dstracks.meanlat.load()
    dstracks.meanlon.load()

    lat = dstracks.meanlat.values          # (tracks, times)
    lon = dstracks.meanlon.values          # (tracks, times), in [0, 360]

    # Convert lon to [-180, 180]
    lon180 = (lon + 180) % 360 - 180

    in_lat = (lat > WAM_LAT_MIN_BUF) & (lat < WAM_LAT_MAX_BUF)
    in_lon = (lon180 > WAM_LON_MIN_BUF) & (lon180 < WAM_LON_MAX_BUF)
    in_wam = (in_lat & in_lon).any(axis=1)   # True if track ever in WAM bbox

    dstracks_wam = dstracks.isel(tracks=in_wam)
    print(f'  {int(in_wam.sum())} / {dstracks.sizes["tracks"]} tracks pass WAM filter')
    return dstracks_wam


def align_times(entr_ds, mask_ds):
    """
    Find 3-hourly entrainment times that exist in the hourly mask dataset.
    Returns:
        times_3h   : 1-D array of datetime64 values (3-hourly, in overlap)
        mask_indices: corresponding positional indices in mask_ds.time
    """
    entr_times = entr_ds.time.values
    mask_times = mask_ds.time.values

    # Build set of (y, m, d, h) tuples from mask for fast lookup
    def to_tuples(times):
        return {
            (pd.Timestamp(t).year, pd.Timestamp(t).month,
             pd.Timestamp(t).day,  pd.Timestamp(t).hour): i
            for i, t in enumerate(times)
        }

    mask_tuple_to_idx = to_tuples(mask_times)

    overlap_entr_idx = []
    overlap_mask_idx = []
    for i, t in enumerate(entr_times):
        ts = pd.Timestamp(t)
        key = (ts.year, ts.month, ts.day, ts.hour)
        if key in mask_tuple_to_idx:
            overlap_entr_idx.append(i)
            overlap_mask_idx.append(mask_tuple_to_idx[key])

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
    track_nums = dstracks_wam.tracks.values.astype(int)   # 0-indexed track IDs
    mask_num_to_out_idx = np.full(max_label + 1, -1, dtype=np.int32)
    for out_i, tn in enumerate(track_nums):
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(OUTPUT_NC), help='Output NetCDF path')
    parser.add_argument('--n-timesteps', type=int, default=None, metavar='N',
                        help='Limit to first N timesteps (for testing)')
    args = parser.parse_args()

    dstracks      = load_track_stats()
    entr_ds       = open_entrainment()
    mask_ds       = open_mcs_mask()

    wam_positions = compute_wam_positions(entr_ds, mask_ds)
    dstracks_wam  = filter_wam_tracks(dstracks)

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
