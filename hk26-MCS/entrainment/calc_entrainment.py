"""
Compute CAPE, CIN, LNB, w_eff (w/sqrt(CAPE) at 500 hPa), brightness temperature
(Tb from OLR) and Tb_diff (Tb - T_LNB) for an analysis region over the full
N2560 RAL3p3 simulation, writing results into a pre-created zarr store.

Usage:
    # Create zarr store once (also run automatically by submit.py):
    python calc_entrainment.py --init --model um_glm_n2560_RAL3p3_tuned_hk26

    # Process one chunk (called by the SLURM array via submit.py):
    python calc_entrainment.py slurm/tasks/<job>.json $SLURM_ARRAY_TASK_ID
"""
import argparse
import json
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import dask.array as dsa
import easygems.healpix as egh
import intake
import numpy as np
import xarray as xr

import models
from metpy.calc import cape_cin, dewpoint_from_relative_humidity, el, parcel_profile
from metpy.units import units

CHUNK_SIZE = 10     # timesteps per chunk / CPUs per SLURM task
N_WORKERS  = 10

warnings.filterwarnings('ignore', message='.*The return type of `Dataset.dims`.*', category=FutureWarning)
warnings.filterwarnings('ignore', message='.*Relative humidity >120%.*', category=UserWarning)
warnings.filterwarnings('ignore', message='.*divide by zero encountered in log.*', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*invalid value encountered in divide.*', category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def hp_mods(ds):
    return ds.rename({'healpix_index': 'cell'}).pipe(egh.attach_coords)


def _region_mask(ds, region_cfg):
    """Boolean cell mask for the analysis region (handles lon wrap-around)."""
    lon_min, lon_max = region_cfg['lon_min'], region_cfg['lon_max']
    lat_min, lat_max = region_cfg['lat_min'], region_cfg['lat_max']
    if lon_min > lon_max:   # region wraps across 0°/360°
        lon_mask = (ds.lon > lon_min) | (ds.lon < lon_max)
    else:
        lon_mask = (ds.lon > lon_min) & (ds.lon < lon_max)
    return lon_mask & (ds.lat > lat_min) & (ds.lat < lat_max)


def open_region_dataset(model, region_cfg):
    zoom = models.MODELS[model]['zoom']
    cat  = intake.open_catalog(models.CATALOG_URL)['UK']
    ds3h = cat[model](zoom=zoom, time='PT3H').to_dask().pipe(hp_mods)
    return ds3h.isel(cell=_region_mask(ds3h, region_cfg))


def open_region_1h_dataset(model, region_cfg):
    zoom = models.MODELS[model]['zoom']
    cat  = intake.open_catalog(models.CATALOG_URL)['UK']
    ds1h = cat[model](zoom=zoom, time='PT1H').to_dask().pipe(hp_mods)
    return ds1h.isel(cell=_region_mask(ds1h, region_cfg))


# ---------------------------------------------------------------------------
# Zarr store initialisation
# ---------------------------------------------------------------------------

def init_zarr(model, region):
    region_cfg = models.REGIONS[region]
    ds = open_region_dataset(model, region_cfg)

    n_times  = ds.sizes['time']
    n_cells  = ds.sizes['cell']
    n_chunks = (n_times + CHUNK_SIZE - 1) // CHUNK_SIZE

    empty    = dsa.full((n_times, n_cells), np.nan, dtype=np.float32,
                        chunks=(CHUNK_SIZE, n_cells))

    template = xr.Dataset(
        {
            'cape':   xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': 'J kg-1',
                                          'long_name': 'Convective available potential energy'}),
            'cin':    xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': 'J kg-1',
                                          'long_name': 'Convective inhibition'}),
            'lnb':    xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': 'hPa',
                                          'long_name': 'Level of neutral buoyancy pressure'}),
            't_lnb':  xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': 'K',
                                          'long_name': 'Temperature at level of neutral buoyancy'}),
            'w_eff':  xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': '-',
                                          'long_name': 'Vertical velocity at 500 hPa normalised by sqrt(CAPE)'}),
            'tb':     xr.DataArray(empty, dims=['time', 'cell'],
                                   attrs={'units': 'K',
                                          'long_name': 'Brightness temperature derived from OLR'}),
            'tb_diff': xr.DataArray(empty, dims=['time', 'cell'],
                                    attrs={'units': 'K',
                                           'long_name': 'Brightness temperature minus T at LNB'}),
        },
        coords={'time': ds.time, 'cell': ds.cell, 'lat': ds.lat, 'lon': ds.lon},
    )

    zarr_path = models.data_dir(model) / f'entrainment_{region}.zarr'
    zarr_path.parent.mkdir(parents=True, exist_ok=True)
    template.to_zarr(zarr_path, mode='w', zarr_format=2)

    done_dir = models.done_dir(model)
    done_dir.mkdir(parents=True, exist_ok=True)
    models.init_donefile(model, region).touch()

    print(f'Created {zarr_path}  shape=({n_times}, {n_cells})')
    print(f'Submit array 0-{n_chunks - 1}  ({n_chunks} jobs)  via: python submit.py --model {model}')


# ---------------------------------------------------------------------------
# CAPE/CIN/LNB computation
# ---------------------------------------------------------------------------

def calc_profile(ta_col, hur_col, p_hpa):
    """CAPE, CIN (J/kg), LNB pressure (hPa) and T at LNB (K) for one vertical profile.

    Profile must be sorted surface-to-top (descending pressure).
    Returns (nan, nan, nan, nan) when fewer than 15 valid dewpoint levels or metpy raises.
    """
    p  = p_hpa * units.hPa
    t  = ta_col * units.K
    rh = (hur_col / 100.0) * units.dimensionless
    try:
        td   = dewpoint_from_relative_humidity(t, rh)
        mask = ~np.isnan(td.magnitude)
        if mask.sum() < 15:
            return np.nan, np.nan, np.nan, np.nan
        prof       = parcel_profile(p[mask], t[mask][0], td[mask][0])
        c, ci      = cape_cin(p[mask], t[mask], td[mask], prof)
        el_p, lnb_t = el(p[mask], t[mask], td[mask])
        return float(c.magnitude), float(ci.magnitude), float(el_p.magnitude), float(lnb_t.magnitude)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def _worker(args):
    """Pool worker: CAPE/CIN/LNB/T_LNB for all cells at a single timestep."""
    i, ta_t, hur_t, p_hpa = args
    n_cells   = ta_t.shape[1]
    cape_arr  = np.empty(n_cells, dtype=np.float32)
    cin_arr   = np.empty(n_cells, dtype=np.float32)
    lnb_arr   = np.empty(n_cells, dtype=np.float32)
    t_lnb_arr = np.empty(n_cells, dtype=np.float32)
    for i in range(n_cells):
        cape_arr[i], cin_arr[i], lnb_arr[i], t_lnb_arr[i] = calc_profile(
            ta_t[:, i], hur_t[:, i], p_hpa
        )
    return cape_arr, cin_arr, lnb_arr, t_lnb_arr


# ---------------------------------------------------------------------------
# Brightness temperature  (Yang & Slingo 2001 / Minnis & Harrison 1984)
# ---------------------------------------------------------------------------

def olr_to_tb(olr):
    """Convert OLR (W m-2) to IR brightness temperature (K).

    Tf = (OLR/sigma)^0.25,  Tb = (-a + sqrt(a^2 + 4*b*Tf)) / (2*b)
    where a=1.228, b=-1.106e-3 K^-1, sigma=5.67e-8 W m-2 K-4.
    Source: PyFLEXTRKR ftfunctions.py (Yang & Slingo 2001).
    """
    a, b, sigma = 1.228, -1.106e-3, 5.67e-8
    tf = (olr / sigma) ** 0.25
    return (-a + np.sqrt(a**2 + 4 * b * tf)) / (2 * b)


# ---------------------------------------------------------------------------
# Chunk processing
# ---------------------------------------------------------------------------

def compute_chunk(chunk_idx, model, region, n_timesteps=None):
    done_file  = models.chunk_donefile(model, chunk_idx)
    if done_file.exists():
        print(f'Chunk {chunk_idx} already done, skipping.')
        return

    region_cfg = models.REGIONS[region]
    ds         = open_region_dataset(model, region_cfg)
    ds1h       = open_region_1h_dataset(model, region_cfg)
    zarr_path  = models.data_dir(model) / f'entrainment_{region}.zarr'

    t_start  = chunk_idx * CHUNK_SIZE
    t_end    = min(t_start + CHUNK_SIZE, ds.sizes['time'])
    if n_timesteps is not None:
        t_end = min(t_start + n_timesteps, t_end)
    n_chunk  = t_end - t_start

    print(f'Chunk {chunk_idx}: time[{t_start}:{t_end}] ({n_chunk} timesteps)')

    ds_chunk  = ds.isel(time=slice(t_start, t_end))
    ds_desc   = ds_chunk.sortby('pressure', ascending=False)

    ta_np   = ds_desc.ta.compute().values
    hur_np  = ds_desc.hur.compute().values
    wa_np   = ds_chunk.wa.sel(pressure=500).compute().values
    p_hpa   = ds_desc.pressure.values.astype(float)
    rlut_np = ds1h.rlut.sel(time=ds_chunk.time).compute().values

    worker_args = [(i, ta_np[i], hur_np[i], p_hpa) for i in range(n_chunk)]

    results = []
    with Pool(N_WORKERS) as pool:
        for i, result in enumerate(pool.imap(_worker, worker_args)):
            print(f'  {i + 1}/{n_chunk}', flush=True)
            results.append(result)

    cape_out   = np.stack([r[0] for r in results])
    cin_out    = np.stack([r[1] for r in results])
    lnb_out    = np.stack([r[2] for r in results])
    t_lnb_out  = np.stack([r[3] for r in results])

    cape_safe  = np.where(cape_out > 0, cape_out, np.nan)
    w_eff_out  = wa_np / np.sqrt(cape_safe)

    tb_out     = olr_to_tb(rlut_np)
    tb_diff_out = tb_out - t_lnb_out

    ds_out = xr.Dataset({
        'cape':    xr.DataArray(cape_out,                      dims=['time', 'cell']),
        'cin':     xr.DataArray(cin_out,                       dims=['time', 'cell']),
        'lnb':     xr.DataArray(lnb_out,                       dims=['time', 'cell']),
        't_lnb':   xr.DataArray(t_lnb_out,                     dims=['time', 'cell']),
        'w_eff':   xr.DataArray(w_eff_out,                     dims=['time', 'cell']),
        'tb':      xr.DataArray(tb_out.astype(np.float32),     dims=['time', 'cell']),
        'tb_diff': xr.DataArray(tb_diff_out.astype(np.float32), dims=['time', 'cell']),
    })
    ds_out.to_zarr(zarr_path, region={'time': slice(t_start, t_end)})

    done_file.touch()
    print(f'Chunk {chunk_idx} written.')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # Two modes:
    #   --init --model <model> [--region <region>]   — initialise zarr store
    #   <json_file> <task_index>                     — process one SLURM task

    if len(sys.argv) >= 3 and not sys.argv[1].startswith('-'):
        # SLURM task mode: positional json_file task_index
        json_file  = Path(sys.argv[1])
        task_index = int(sys.argv[2])
        task_cfg   = json.loads(json_file.read_text())
        model      = task_cfg['model']
        region     = task_cfg['region']
        chunk      = task_cfg['tasks'][task_index]['chunk']
        compute_chunk(chunk, model, region)
        return

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--init', action='store_true',
                       help='Create empty zarr store (run once, or let submit.py handle this)')
    parser.add_argument('--n-timesteps', type=int, default=None, metavar='N',
                        help='Limit timesteps processed (for testing)')
    models.add_model_arg(parser)
    models.add_region_arg(parser)
    args = parser.parse_args()

    if args.init:
        init_zarr(args.model, args.region)


if __name__ == '__main__':
    main()
