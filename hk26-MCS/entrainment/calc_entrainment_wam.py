"""
Compute CAPE, CIN, LNB, w_eff (w/sqrt(CAPE) at 500 hPa), brightness temperature
(Tb from OLR) and Tb_diff (Tb - T_LNB) for the West African Monsoon (WAM) region
over the full N2560 RAL3p3 simulation, writing results into a pre-created zarr store.

Usage:
    python calc_entrainment_wam.py --init        # create zarr store once before submitting
    python calc_entrainment_wam.py --chunk N     # process chunk N (called by SLURM)
"""
import argparse
import warnings
from multiprocessing import Pool
from pathlib import Path

import dask.array as dsa
import easygems.healpix as egh
import intake
import numpy as np
import xarray as xr

from metpy.calc import cape_cin, dewpoint_from_relative_humidity, el, parcel_profile
from metpy.units import units

# --- Config ---
SIM = 'um_glm_n2560_RAL3p3_tuned_hk26'
ZOOM = 9
CATALOG_URL = 'https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml'
CHUNK_SIZE = 80          # 10 days × 8 timesteps/day at 3-hourly
N_WORKERS = 10
ZARR_PATH = Path('entrainment_wam.zarr')
DONE_DIR = Path('entrainment_done')

# Brightness temperature constants
SIGMA = 5.67e-8      # Stefan-Boltzmann constant (W m-2 K-4)
TB_A = 1.228
TB_B = -1.106e-3

warnings.filterwarnings('ignore', message='.*The return type of `Dataset.dims`.*', category=FutureWarning)
warnings.filterwarnings('ignore', message='.*Relative humidity >120%.*', category=UserWarning)
warnings.filterwarnings('ignore', message='.*divide by zero encountered in log.*', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*invalid value encountered in divide.*', category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def hp_mods(ds):
    return ds.rename({'healpix_index': 'cell'}).pipe(egh.attach_coords)


def _wam_mask(ds):
    return ((ds.lon > 340) | (ds.lon < 20)) & (ds.lat > 2) & (ds.lat < 15)


def open_wam_dataset():
    cat = intake.open_catalog(CATALOG_URL)['UK']
    ds3h = cat[SIM](zoom=ZOOM, time='PT3H').to_dask().pipe(hp_mods)
    return ds3h.isel(cell=_wam_mask(ds3h))


def open_wam_1h_dataset():
    cat = intake.open_catalog(CATALOG_URL)['UK']
    ds1h = cat[SIM](zoom=ZOOM, time='PT1H').to_dask().pipe(hp_mods)
    return ds1h.isel(cell=_wam_mask(ds1h))


# ---------------------------------------------------------------------------
# Zarr store initialisation
# ---------------------------------------------------------------------------

def init_zarr(ds):
    """Create an empty zarr store using xarray/dask with the correct shape and coordinates."""
    n_times = ds.sizes['time']
    n_cells = ds.sizes['cell']
    n_chunks = (n_times + CHUNK_SIZE - 1) // CHUNK_SIZE

    empty = dsa.full((n_times, n_cells), np.nan, dtype=np.float32,
                     chunks=(CHUNK_SIZE, n_cells))

    template = xr.Dataset(
        {
            'cape': xr.DataArray(empty, dims=['time', 'cell'],
                                 attrs={'units': 'J kg-1',
                                        'long_name': 'Convective available potential energy'}),
            'cin': xr.DataArray(empty, dims=['time', 'cell'],
                                attrs={'units': 'J kg-1',
                                       'long_name': 'Convective inhibition'}),
            'lnb': xr.DataArray(empty, dims=['time', 'cell'],
                                attrs={'units': 'hPa',
                                       'long_name': 'Level of neutral buoyancy pressure'}),
            't_lnb': xr.DataArray(empty, dims=['time', 'cell'],
                                  attrs={'units': 'K',
                                         'long_name': 'Temperature at level of neutral buoyancy'}),
            'w_eff': xr.DataArray(empty, dims=['time', 'cell'],
                                  attrs={'units': 'm s-1 (J kg-1)-0.5',
                                         'long_name': 'Vertical velocity at 500 hPa normalised by sqrt(CAPE)'}),
            'tb': xr.DataArray(empty, dims=['time', 'cell'],
                               attrs={'units': 'K',
                                      'long_name': 'Brightness temperature derived from OLR'}),
            'tb_diff': xr.DataArray(empty, dims=['time', 'cell'],
                                    attrs={'units': 'K',
                                           'long_name': 'Brightness temperature minus T at LNB'}),
        },
        coords={
            'time': ds.time,
            'cell': ds.cell,
            'lat': ds.lat,
            'lon': ds.lon,
        },
    )

    template.to_zarr(ZARR_PATH, mode='w', zarr_format=2)

    DONE_DIR.mkdir(exist_ok=True)
    print(f'Created {ZARR_PATH}  shape=({n_times}, {n_cells})')
    print(f'Submit SLURM array 0-{n_chunks - 1}  ({n_chunks} jobs)')


# ---------------------------------------------------------------------------
# CAPE/CIN/LNB computation
# ---------------------------------------------------------------------------

def calc_profile(ta_col, hur_col, p_hpa):
    """CAPE, CIN (J/kg), LNB pressure (hPa) and T at LNB (K) for one vertical profile.

    Profile must be sorted surface-to-top (descending pressure).
    Returns (nan, nan, nan, nan) when fewer than 15 valid dewpoint levels or metpy raises.
    """
    p = p_hpa * units.hPa
    t = ta_col * units.K
    rh = (hur_col / 100.0) * units.dimensionless
    try:
        td = dewpoint_from_relative_humidity(t, rh)
        mask = ~np.isnan(td.magnitude)
        if mask.sum() < 15:
            return np.nan, np.nan, np.nan, np.nan
        prof = parcel_profile(p[mask], t[mask][0], td[mask][0])
        c, ci = cape_cin(p[mask], t[mask], td[mask], prof)
        el_p, lnb_t = el(p[mask], t[mask], td[mask])
        return float(c.magnitude), float(ci.magnitude), float(el_p.magnitude), float(lnb_t.magnitude)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def _worker(args):
    """Pool worker: CAPE/CIN/LNB/T_LNB for all cells at a single timestep.

    args: (ta_t, hur_t, p_hpa)
        ta_t, hur_t  shape (n_levels, n_cells), sorted surface-to-top
        p_hpa        shape (n_levels,),          sorted surface-to-top
    """
    i, ta_t, hur_t, p_hpa = args
    n_cells = ta_t.shape[1]
    cape_arr = np.empty(n_cells, dtype=np.float32)
    cin_arr = np.empty(n_cells, dtype=np.float32)
    lnb_arr = np.empty(n_cells, dtype=np.float32)
    t_lnb_arr = np.empty(n_cells, dtype=np.float32)
    for i in range(n_cells)[:5]:
        cape_arr[i], cin_arr[i], lnb_arr[i], t_lnb_arr[i] = calc_profile(
            ta_t[:, i], hur_t[:, i], p_hpa
        )
    return cape_arr, cin_arr, lnb_arr, t_lnb_arr


# ---------------------------------------------------------------------------
# Chunk processing
# ---------------------------------------------------------------------------

def compute_chunk(chunk_idx, n_timesteps=None):
    done_file = DONE_DIR / f'chunk_{chunk_idx:03d}.done'
    if done_file.exists():
        print(f'Chunk {chunk_idx} already done, skipping.')
        return

    ds = open_wam_dataset()
    ds1h = open_wam_1h_dataset()

    t_start = chunk_idx * CHUNK_SIZE
    t_end = min(t_start + CHUNK_SIZE, ds.sizes['time'])
    if n_timesteps is not None:
        t_end = min(t_start + n_timesteps, t_end)
    n_chunk = t_end - t_start

    print(f'Chunk {chunk_idx}: time[{t_start}:{t_end}] ({n_chunk} timesteps)')

    # ta/hur sorted surface-to-top (descending pressure) for metpy; wa at 500 hPa only.
    ds_chunk = ds.isel(time=slice(t_start, t_end))
    ds_desc = ds_chunk.sortby('pressure', ascending=False)

    ta_np = ds_desc.ta.compute().values              # (n_chunk, n_levels, n_cells)
    hur_np = ds_desc.hur.compute().values
    wa_np = ds_chunk.wa.sel(pressure=500).compute().values   # (n_chunk, n_cells)
    p_hpa = ds_desc.pressure.values.astype(float)   # (n_levels,) descending

    # Select OLR from 1H dataset at the 3H times.
    rlut_np = ds1h.rlut.sel(time=ds_chunk.time).compute().values  # (n_chunk, n_cells)

    worker_args = [(i, ta_np[i], hur_np[i], p_hpa) for i in range(n_chunk)]

    results = []
    with Pool(N_WORKERS) as pool:
        for i, result in enumerate(pool.imap(_worker, worker_args)):
            print(f'  {i + 1}/{n_chunk}', flush=True)
            results.append(result)

    cape_out = np.stack([r[0] for r in results])    # (n_chunk, n_cells)
    cin_out = np.stack([r[1] for r in results])
    lnb_out = np.stack([r[2] for r in results])
    t_lnb_out = np.stack([r[3] for r in results])

    # w_eff = w/sqrt(CAPE) at 500 hPa; mask zero/nan CAPE.
    cape_safe = np.where(cape_out > 0, cape_out, np.nan)
    w_eff_out = wa_np / np.sqrt(cape_safe)          # (n_chunk, n_cells)

    # Brightness temperature from OLR.
    Tf = (rlut_np / SIGMA) ** 0.25
    tb_out = (TB_A + np.sqrt(TB_A**2 - 4 * TB_B * Tf)) / (2 * TB_B)

    tb_diff_out = tb_out - t_lnb_out

    ds_out = xr.Dataset({
        'cape': xr.DataArray(cape_out, dims=['time', 'cell']),
        'cin': xr.DataArray(cin_out, dims=['time', 'cell']),
        'lnb': xr.DataArray(lnb_out, dims=['time', 'cell']),
        't_lnb': xr.DataArray(t_lnb_out, dims=['time', 'cell']),
        'w_eff': xr.DataArray(w_eff_out, dims=['time', 'cell']),
        'tb': xr.DataArray(tb_out.astype(np.float32), dims=['time', 'cell']),
        'tb_diff': xr.DataArray(tb_diff_out.astype(np.float32), dims=['time', 'cell']),
    })
    ds_out.to_zarr(ZARR_PATH, region={'time': slice(t_start, t_end)})

    done_file.touch()
    print(f'Chunk {chunk_idx} written.')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--init', action='store_true',
                       help='Create empty zarr store (run once before submitting)')
    group.add_argument('--chunk', type=int, metavar='N',
                       help='Chunk index to process (set by SLURM_ARRAY_TASK_ID)')
    parser.add_argument('--n-timesteps', type=int, default=None, metavar='N',
                        help='Limit timesteps processed (for testing)')
    args = parser.parse_args()

    if args.init:
        init_zarr(open_wam_dataset())
    else:
        compute_chunk(args.chunk, n_timesteps=args.n_timesteps)


if __name__ == '__main__':
    main()
